#include "test_common.h"
#include "qaws_platform.h"

/* ------------------------------------------------------------------ */
/*  Helper: map normalized distance [0,1] to blue-green-red colormap  */
/* ------------------------------------------------------------------ */
static void distance_color(double t, double* r, double* g, double* b)
{
	if (t < 0.0) t = 0.0;
	if (t > 1.0) t = 1.0;
	if (t < 0.25)
	{
		*r = 0.0;
		*g = t / 0.25;
		*b = 1.0;
	}
	else if (t < 0.5)
	{
		*r = 0.0;
		*g = 1.0;
		*b = 1.0 - (t - 0.25) / 0.25;
	}
	else if (t < 0.75)
	{
		*r = (t - 0.5) / 0.25;
		*g = 1.0;
		*b = 0.0;
	}
	else
	{
		*r = 1.0;
		*g = 1.0 - (t - 0.75) / 0.25;
		*b = 0.0;
	}
}

/* ------------------------------------------------------------------ */
/*  Helper: write heat distance field as OBJ with per-vertex color    */
/* ------------------------------------------------------------------ */
static void obj_heat_distance_mesh(
	obj_writer* w,
	qaws_surface const* surf,
	qaws_tessellation_desc const* tess_desc,
	qaws_scalar const* u_buf,
	qaws_scalar const* v_buf,
	qaws_scalar const* dist,
	unsigned int count)
{
	unsigned int i;
	qaws_scalar max_dist = 0;
	unsigned int first_v, first_n;

	/* Find max distance for normalization */
	for (i = 0; i < count; i++)
	{
		if (dist[i] > max_dist)
			max_dist = dist[i];
	}
	if ((double)max_dist < 1e-10) max_dist = QAWS_LITERAL(1.0);

	/* Write colored vertices + normals */
	first_v = w->vertex_count + 1;
	first_n = w->normal_count + 1;
	for (i = 0; i < count; i++)
	{
		qaws_surface_eval_result r;
		double cr, cg, cb;
		memset(&r, 0, sizeof(r));
		qaws_surface_evaluate(surf, u_buf[i], v_buf[i],
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
		distance_color((double)dist[i] / (double)max_dist, &cr, &cg, &cb);
		obj_vertex_color(w, r.position, cr, cg, cb);
		obj_normal(w, r.normal);
	}

	/* Re-tessellate with the same desc to get triangle indices that
	   match the vertex order from the geodesic distance call. */
	{
		qaws_tessellation_vertex* tv = NULL;
		unsigned int* idx = NULL;
		unsigned int tv_count = 0, idx_count = 0;
		unsigned int tv_cap = 65536, idx_cap = 262144;
		qaws_status s;

		tv = (qaws_tessellation_vertex*)malloc(tv_cap * sizeof(qaws_tessellation_vertex));
		idx = (unsigned int*)malloc(idx_cap * sizeof(unsigned int));
		if (!tv || !idx) { free(tv); free(idx); return; }

		s = qaws_surface_tessellate(surf, tess_desc, tv, tv_cap, &tv_count,
			idx, idx_cap, &idx_count);
		if (s == QAWS_STATUS_OK && tv_count == count)
		{
			unsigned int ti;
			for (ti = 0; ti + 2 < idx_count; ti += 3)
			{
				unsigned int a = first_v + idx[ti];
				unsigned int b = first_v + idx[ti + 1];
				unsigned int c = first_v + idx[ti + 2];
				unsigned int na = first_n + idx[ti];
				unsigned int nb = first_n + idx[ti + 1];
				unsigned int nc = first_n + idx[ti + 2];
				fprintf(w->fp, "f %u//%u %u//%u %u//%u\n",
					a, na, b, nb, c, nc);
			}
		}

		free(tv);
		free(idx);
	}

	(void)first_v; (void)first_n;
}

/* ------------------------------------------------------------------ */
/*  Helper: render surface as regular grid colored by heat distance.   */
/*  Uses inverse-distance-weighted (IDW) interpolation from the K      */
/*  nearest heat vertices in UV space for smooth, tessellation-free    */
/*  coloring.                                                          */
/* ------------------------------------------------------------------ */
static void obj_heat_surface_grid(
	obj_writer* w,
	qaws_surface const* surf,
	qaws_scalar const* heat_u,
	qaws_scalar const* heat_v,
	qaws_scalar const* heat_dist,
	unsigned int heat_count,
	unsigned int grid_n)
{
	unsigned int i, j, vi;
	qaws_scalar max_dist = 0;
	unsigned int first_v = w->vertex_count + 1;
	unsigned int first_n = w->normal_count + 1;
	enum { K = 6 };

	for (vi = 0; vi < heat_count; vi++)
		if (heat_dist[vi] > max_dist)
			max_dist = heat_dist[vi];
	if ((double)max_dist < 1e-10) max_dist = QAWS_LITERAL(1.0);

	/* Emit vertices on a regular grid, colored by IDW-interpolated heat */
	for (i = 0; i <= grid_n; i++)
		for (j = 0; j <= grid_n; j++)
		{
			qaws_scalar u = (qaws_scalar)j / (qaws_scalar)grid_n;
			qaws_scalar v = (qaws_scalar)i / (qaws_scalar)grid_n;
			qaws_surface_eval_result r;
			double cr, cg, cb;
			qaws_scalar d;

			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surf, u, v,
				QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);

			/* IDW interpolation from K nearest heat vertices */
			{
				qaws_scalar best_dsq[K];
				unsigned int best_idx[K];
				int ki;
				qaws_scalar w_sum = QAWS_ZERO, d_sum = QAWS_ZERO;

				for (ki = 0; ki < K; ki++)
				{
					best_dsq[ki] = (qaws_scalar)1e30;
					best_idx[ki] = 0;
				}

				/* Find K nearest in UV space */
				for (vi = 0; vi < heat_count; vi++)
				{
					qaws_scalar du = u - heat_u[vi];
					qaws_scalar dv = v - heat_v[vi];
					qaws_scalar dsq = du * du + dv * dv;
					/* Insert into sorted top-K if closer than worst */
					if (dsq < best_dsq[K - 1])
					{
						int slot = K - 1;
						best_dsq[slot] = dsq;
						best_idx[slot] = vi;
						/* Bubble up */
						while (slot > 0 && best_dsq[slot] < best_dsq[slot - 1])
						{
							qaws_scalar td = best_dsq[slot];
							unsigned int ti2 = best_idx[slot];
							best_dsq[slot] = best_dsq[slot - 1];
							best_idx[slot] = best_idx[slot - 1];
							best_dsq[slot - 1] = td;
							best_idx[slot - 1] = ti2;
							slot--;
						}
					}
				}

				/* Exact hit check */
				if (best_dsq[0] < (qaws_scalar)1e-20)
				{
					d = heat_dist[best_idx[0]];
				}
				else
				{
					/* Inverse-distance weighting: w_i = 1/dist_i^2 */
					for (ki = 0; ki < K && best_dsq[ki] < (qaws_scalar)1e20; ki++)
					{
						qaws_scalar wi = QAWS_ONE / best_dsq[ki];
						w_sum += wi;
						d_sum += wi * heat_dist[best_idx[ki]];
					}
					d = d_sum / w_sum;
				}
			}

			distance_color((double)d / (double)max_dist, &cr, &cg, &cb);
			obj_vertex_color(w, r.position, cr, cg, cb);
			obj_normal(w, r.normal);
		}

	/* Emit quad faces */
	for (i = 0; i < grid_n; i++)
		for (j = 0; j < grid_n; j++)
		{
			unsigned int a = first_v + i * (grid_n + 1) + j;
			unsigned int b = a + 1;
			unsigned int c = a + grid_n + 2;
			unsigned int d_idx = a + grid_n + 1;
			unsigned int na = first_n + (a - first_v);
			unsigned int nb = first_n + (b - first_v);
			unsigned int nc = first_n + (c - first_v);
			unsigned int nd = first_n + (d_idx - first_v);
			fprintf(w->fp, "f %u//%u %u//%u %u//%u %u//%u\n",
				a, na, b, nb, c, nc, d_idx, nd);
		}
}

/* ------------------------------------------------------------------ */
/*  Helper: create multi-bump terrain (11x11 cubic B-spline)           */
/* ------------------------------------------------------------------ */
static qaws_status create_hilly_surface(qaws_surface** out_surf)
{
	/* Cubic B-spline surface, 11x11 control points.
	   Local support (4 spans each) means bumps are distinct.
	   Grid covers [0,10]x[0,10] in XY space.
	   cp[row * 11 + col]  where col -> u, row -> v. */
	qaws_vec3 cp[11 * 11];
	qaws_surface_bspline_desc bd;
	unsigned int i, j;

	/* Base grid at z=0 */
	for (i = 0; i < 11; i++)
		for (j = 0; j < 11; j++)
		{
			cp[i * 11 + j].x = (qaws_scalar)j;
			cp[i * 11 + j].y = (qaws_scalar)i;
			cp[i * 11 + j].z = QAWS_ZERO;
		}

	/* Bump A: center at (col=2, row=2), peak z=2.5 */
	cp[2 * 11 + 2].z = (qaws_scalar)2.5;
	cp[1 * 11 + 2].z = (qaws_scalar)1.0;
	cp[3 * 11 + 2].z = (qaws_scalar)1.0;
	cp[2 * 11 + 1].z = (qaws_scalar)1.0;
	cp[2 * 11 + 3].z = (qaws_scalar)1.0;

	/* Bump B: center at (col=7, row=3), peak z=2.0 */
	cp[3 * 11 + 7].z = (qaws_scalar)2.0;
	cp[2 * 11 + 7].z = (qaws_scalar)0.8;
	cp[4 * 11 + 7].z = (qaws_scalar)0.8;
	cp[3 * 11 + 6].z = (qaws_scalar)0.8;
	cp[3 * 11 + 8].z = (qaws_scalar)0.8;

	/* Bump C: center at (col=3, row=7), peak z=2.2 */
	cp[7 * 11 + 3].z = (qaws_scalar)2.2;
	cp[6 * 11 + 3].z = (qaws_scalar)0.9;
	cp[8 * 11 + 3].z = (qaws_scalar)0.9;
	cp[7 * 11 + 2].z = (qaws_scalar)0.9;
	cp[7 * 11 + 4].z = (qaws_scalar)0.9;

	/* Bump D: center at (col=8, row=8), peak z=3.0 (tallest) */
	cp[8 * 11 + 8].z = (qaws_scalar)3.0;
	cp[7 * 11 + 8].z = (qaws_scalar)1.2;
	cp[9 * 11 + 8].z = (qaws_scalar)1.2;
	cp[8 * 11 + 7].z = (qaws_scalar)1.2;
	cp[8 * 11 + 9].z = (qaws_scalar)1.2;

	memset(&bd, 0, sizeof(bd));
	bd.u_degree = 3;
	bd.v_degree = 3;
	bd.control_points = cp;
	bd.u_point_count = 11;
	bd.v_point_count = 11;
	bd.u_knot_count = 0; /* auto uniform clamped */
	bd.v_knot_count = 0;

	return qaws_surface_create_bspline(&bd, out_surf);
}

/* ------------------------------------------------------------------ */
/*  Helper: tube that sticks to a surface.  Samples the curve, projects*/
/*  each sample back onto the surface via closest-point, then builds   */
/*  the tube mesh from the projected positions.                        */
/* ------------------------------------------------------------------ */
static void obj_tube_on_surface(obj_writer* w, qaws_curve const* curve,
	qaws_surface const* surface,
	unsigned int length_segments, unsigned int circle_segments,
	qaws_scalar radius)
{
	unsigned int si, ci;
	unsigned int first_vi = w->vertex_count + 1;
	unsigned int first_ni = w->normal_count + 1;
	qaws_range range = qaws_curve_get_parameter_range(curve);
	qaws_scalar t_min = range.min_value;
	qaws_scalar t_len = range.max_value - range.min_value;
	qaws_vec3 prev_T = {0,0,0}, prev_N = {0,0,0}, prev_B = {0,0,0};
	int frame_initialized = 0;
	unsigned int rings_emitted = 0;
	qaws_vec3 prev_pos = {0,0,0};

	for (si = 0; si <= length_segments; si++)
	{
		qaws_scalar t = t_min + t_len * (qaws_scalar)si / (qaws_scalar)length_segments;
		qaws_vec3 center, T_vec, N_vec, B_vec;
		qaws_eval_result_3d r;
		qaws_scalar out_u, out_v;
		qaws_vec3 closest;
		qaws_surface_eval_result sr;

		if (qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION, &r) != QAWS_STATUS_OK)
			continue;

		/* Project the curve sample onto the surface */
		{
			int use_projection = 0;
			if (qaws_surface_find_closest_point(surface, r.position,
				&out_u, &out_v, &closest) == QAWS_STATUS_OK)
			{
				/* Sanity check: if the projected point is far from the curve
				   sample, Newton converged to a wrong local minimum. */
				qaws_scalar dx = closest.x - r.position.x;
				qaws_scalar dy = closest.y - r.position.y;
				qaws_scalar dz = closest.z - r.position.z;
				qaws_scalar dist_sq = dx * dx + dy * dy + dz * dz;
				if (dist_sq < radius * radius * (qaws_scalar)100.0)
					use_projection = 1;
			}

			if (use_projection)
			{
				memset(&sr, 0, sizeof(sr));
				qaws_surface_evaluate(surface, out_u, out_v,
					QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &sr);
				center.x = sr.position.x + radius * sr.normal.x;
				center.y = sr.position.y + radius * sr.normal.y;
				center.z = sr.position.z + radius * sr.normal.z;
			}
			else
			{
				center = r.position;
			}
		}

		/* Tangent from consecutive projected positions */
		if (si == 0)
		{
			prev_pos = center;
			continue;
		}

		T_vec.x = center.x - prev_pos.x;
		T_vec.y = center.y - prev_pos.y;
		T_vec.z = center.z - prev_pos.z;
		{
			double len = sqrt((double)(T_vec.x * T_vec.x + T_vec.y * T_vec.y + T_vec.z * T_vec.z));
			if (len < 1e-15)
			{
				prev_pos = center;
				continue;
			}
			T_vec.x = (qaws_scalar)(T_vec.x / len);
			T_vec.y = (qaws_scalar)(T_vec.y / len);
			T_vec.z = (qaws_scalar)(T_vec.z / len);
		}

		if (frame_initialized && v3_dot(prev_T, T_vec) < 0.0)
		{
			T_vec.x = -T_vec.x;
			T_vec.y = -T_vec.y;
			T_vec.z = -T_vec.z;
		}

		if (!frame_initialized)
		{
			build_initial_frame(T_vec, &N_vec, &B_vec);
			frame_initialized = 1;
			/* Emit ring for the previous (skipped) position too */
			{
				qaws_vec3 first_center = prev_pos;
				for (ci = 0; ci < circle_segments; ci++)
				{
					double angle = 2.0 * M_PI * (double)ci / (double)circle_segments;
					double ca = cos(angle), sa = sin(angle);
					qaws_vec3 p, n;
					n.x = (qaws_scalar)(N_vec.x * ca + B_vec.x * sa);
					n.y = (qaws_scalar)(N_vec.y * ca + B_vec.y * sa);
					n.z = (qaws_scalar)(N_vec.z * ca + B_vec.z * sa);
					p.x = first_center.x + radius * n.x;
					p.y = first_center.y + radius * n.y;
					p.z = first_center.z + radius * n.z;
					obj_vertex(w, p);
					obj_normal(w, n);
				}
				rings_emitted++;
			}
		}
		else
		{
			double d = v3_dot(prev_T, T_vec);
			if (d > 1.0) d = 1.0;
			if (d < -1.0) d = -1.0;
			if (d > 0.9999999)
			{
				N_vec = prev_N;
				B_vec = prev_B;
			}
			else
			{
				qaws_vec3 rot_axis = v3_cross(prev_T, T_vec);
				double sin_a, cos_a;
				v3_normalize(&rot_axis);
				cos_a = d;
				sin_a = sqrt(1.0 - d * d);
				N_vec = v3_rotate(prev_N, rot_axis, cos_a, sin_a);
				B_vec = v3_rotate(prev_B, rot_axis, cos_a, sin_a);
				v3_normalize(&N_vec);
				v3_normalize(&B_vec);
			}
		}

		prev_T = T_vec;
		prev_N = N_vec;
		prev_B = B_vec;
		prev_pos = center;

		for (ci = 0; ci < circle_segments; ci++)
		{
			double angle = 2.0 * M_PI * (double)ci / (double)circle_segments;
			double ca = cos(angle), sa = sin(angle);
			qaws_vec3 p, n;
			n.x = (qaws_scalar)(N_vec.x * ca + B_vec.x * sa);
			n.y = (qaws_scalar)(N_vec.y * ca + B_vec.y * sa);
			n.z = (qaws_scalar)(N_vec.z * ca + B_vec.z * sa);
			p.x = center.x + radius * n.x;
			p.y = center.y + radius * n.y;
			p.z = center.z + radius * n.z;
			obj_vertex(w, p);
			obj_normal(w, n);
		}
		rings_emitted++;
	}

	/* Emit quads between consecutive rings */
	{
		unsigned int ri;
		for (ri = 0; ri + 1 < rings_emitted; ri++)
		{
			for (ci = 0; ci < circle_segments; ci++)
			{
				unsigned int c_next = (ci + 1) % circle_segments;
				unsigned int a = first_vi + ri * circle_segments + ci;
				unsigned int b = first_vi + ri * circle_segments + c_next;
				unsigned int c = first_vi + (ri + 1) * circle_segments + c_next;
				unsigned int d_idx = first_vi + (ri + 1) * circle_segments + ci;
				unsigned int na = first_ni + ri * circle_segments + ci;
				unsigned int nb = first_ni + ri * circle_segments + c_next;
				unsigned int nc = first_ni + (ri + 1) * circle_segments + c_next;
				unsigned int nd = first_ni + (ri + 1) * circle_segments + ci;
				fprintf(w->fp, "f %u//%u %u//%u %u//%u %u//%u\n",
					a, na, b, nb, c, nc, d_idx, nd);
			}
		}
	}
}

/* ------------------------------------------------------------------ */
/*  Test: edge-flip geodesic on flat surface (should match straight    */
/*  line distance within mesh discretization error)                    */
/* ------------------------------------------------------------------ */
static void test_geodesic_flip_flat(void)
{
	qaws_surface* surf = NULL;
	qaws_curve* geo = NULL;
	qaws_status s;

	printf("test_geodesic_flip_flat\n");

	/* Flat bilinear surface [0,2]x[0,2] at z=0 */
	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = 0; bd.p00.y = 0; bd.p00.z = 0;
		bd.p10.x = 2; bd.p10.y = 0; bd.p10.z = 0;
		bd.p01.x = 0; bd.p01.y = 2; bd.p01.z = 0;
		bd.p11.x = 2; bd.p11.y = 2; bd.p11.z = 0;
		s = qaws_surface_create_bilinear(&bd, &surf);
		TEST_ASSERT_STATUS(s);
	}

	/* Geodesic from (0.2, 0.2) to (0.8, 0.8) */
	{
		qaws_tessellation_desc td;
		td.max_depth = 5;
		td.curvature_threshold = QAWS_LITERAL(0.01);
		td.max_edge_length = QAWS_LITERAL(0.5); /* Force subdivision on flat surface */

		s = qaws_surface_compute_geodesic_flip(surf,
			(qaws_scalar)0.2, (qaws_scalar)0.2,
			(qaws_scalar)0.8, (qaws_scalar)0.8,
			&td, 0, &geo);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(geo != NULL, "flip geodesic produced a curve");
	}

	/* On a flat surface, the geodesic should be approximately a straight line.
	   Check that the arc length is close to the Euclidean distance. */
	if (geo)
	{
		qaws_scalar arc_len = 0;
		qaws_range range = qaws_curve_get_parameter_range(geo);
		qaws_status as = qaws_curve_compute_arc_length(geo,
			range.min_value, range.max_value, &arc_len);
		TEST_ASSERT_STATUS(as);

		/* Expected straight-line distance on the flat surface:
		   Start at surface(0.2,0.2) = (0.4, 0.4, 0)
		   End at surface(0.8,0.8) = (1.6, 1.6, 0)
		   Distance = sqrt((1.2)^2 + (1.2)^2) = 1.2*sqrt(2) ~ 1.697 */
		{
			qaws_scalar expected = (qaws_scalar)(1.2 * 1.41421356);
			qaws_scalar error = QAWS_FABS(arc_len - expected) / expected;
			printf("  flip flat: arc_len=%.4f expected=%.4f error=%.2f%%\n",
				(double)arc_len, (double)expected, (double)(error * 100));
			TEST_ASSERT(error < (qaws_scalar)0.20,
				"flip geodesic length within 20% of straight line");
		}

		/* Check that points lie on z=0 plane */
		{
			unsigned int i;
			int on_surface = 1;
			for (i = 0; i < 10; i++)
			{
				qaws_scalar t = range.min_value +
					(range.max_value - range.min_value) * (qaws_scalar)i / (qaws_scalar)9;
				qaws_eval_result_3d r;
				memset(&r, 0, sizeof(r));
				qaws_curve_evaluate_3d(geo, t, QAWS_EVAL_FLAG_POSITION, &r);
				if (QAWS_FABS(r.position.z) > (qaws_scalar)0.1)
					on_surface = 0;
			}
			TEST_ASSERT(on_surface, "flip geodesic lies on flat surface (z~0)");
		}
	}

	qaws_curve_destroy(geo);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: edge-flip geodesic on dome surface                           */
/* ------------------------------------------------------------------ */
static void test_geodesic_flip_dome(void)
{
	qaws_surface* surf = NULL;
	qaws_curve* geo = NULL;
	qaws_status s;

	printf("test_geodesic_flip_dome\n");

	/* Biquadratic dome: center raised */
	{
		qaws_surface_biquadratic_desc bd;
		/* Row 0 (v=0): flat bottom */
		bd.control_points[0].x = 0; bd.control_points[0].y = 0; bd.control_points[0].z = 0;
		bd.control_points[1].x = 1; bd.control_points[1].y = 0; bd.control_points[1].z = 0;
		bd.control_points[2].x = 2; bd.control_points[2].y = 0; bd.control_points[2].z = 0;
		/* Row 1 (v=0.5): middle with dome */
		bd.control_points[3].x = 0; bd.control_points[3].y = 1; bd.control_points[3].z = 0;
		bd.control_points[4].x = 1; bd.control_points[4].y = 1; bd.control_points[4].z = 1;
		bd.control_points[5].x = 2; bd.control_points[5].y = 1; bd.control_points[5].z = 0;
		/* Row 2 (v=1): flat top */
		bd.control_points[6].x = 0; bd.control_points[6].y = 2; bd.control_points[6].z = 0;
		bd.control_points[7].x = 1; bd.control_points[7].y = 2; bd.control_points[7].z = 0;
		bd.control_points[8].x = 2; bd.control_points[8].y = 2; bd.control_points[8].z = 0;
		s = qaws_surface_create_biquadratic(&bd, &surf);
		TEST_ASSERT_STATUS(s);
	}

	{
		qaws_tessellation_desc td;
		td.max_depth = 5;
		td.curvature_threshold = QAWS_LITERAL(0.05);
		td.max_edge_length = QAWS_ZERO;

		s = qaws_surface_compute_geodesic_flip(surf,
			(qaws_scalar)0.1, (qaws_scalar)0.1,
			(qaws_scalar)0.9, (qaws_scalar)0.9,
			&td, 0, &geo);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(geo != NULL, "flip geodesic on dome produced a curve");
	}

	if (geo)
	{
		qaws_scalar arc_len = 0;
		qaws_range range = qaws_curve_get_parameter_range(geo);
		qaws_curve_compute_arc_length(geo, range.min_value, range.max_value, &arc_len);
		printf("  dome flip geodesic arc_len=%.4f\n", (double)arc_len);
		TEST_ASSERT(arc_len > QAWS_LITERAL(0.1), "dome geodesic has positive length");

		/* Write OBJ for visual inspection */
		{
			obj_writer obj;
			svg_ensure_output_dir();
			if (obj_open(&obj, OBJ_OUTPUT_DIR "/44_geodesic_flip_dome.obj",
				OBJ_OUTPUT_DIR "/44_geodesic_flip_dome.mtl"))
			{
				obj_material(&obj, "dome", 0.5, 0.5, 0.8);
				obj_group(&obj, "dome_surface");
				obj_use_material(&obj, "dome");
				obj_surface_mesh(&obj, surf, 30, 30);

				obj_material(&obj, "geodesic", 1.0, 0.3, 0.1);
				obj_group(&obj, "geodesic_path");
				obj_use_material(&obj, "geodesic");
				obj_tube_on_surface(&obj, geo, surf, 80, 6, (qaws_scalar)0.02);

				obj_close(&obj);
			}
		}
	}

	qaws_curve_destroy(geo);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: heat method distance field on flat surface                    */
/* ------------------------------------------------------------------ */
static void test_heat_distance_flat(void)
{
	qaws_surface* surf = NULL;
	qaws_status s;
	qaws_scalar* dist = NULL;
	qaws_scalar* u_buf = NULL;
	qaws_scalar* v_buf = NULL;
	unsigned int count = 0;
	unsigned int cap = 65536;

	printf("test_heat_distance_flat\n");

	/* Flat bilinear surface */
	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = 0; bd.p00.y = 0; bd.p00.z = 0;
		bd.p10.x = 2; bd.p10.y = 0; bd.p10.z = 0;
		bd.p01.x = 0; bd.p01.y = 2; bd.p01.z = 0;
		bd.p11.x = 2; bd.p11.y = 2; bd.p11.z = 0;
		s = qaws_surface_create_bilinear(&bd, &surf);
		TEST_ASSERT_STATUS(s);
	}

	dist = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));
	u_buf = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));
	v_buf = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));

	{
		qaws_tessellation_desc td;
		td.max_depth = 5;
		td.curvature_threshold = QAWS_LITERAL(0.01);
		td.max_edge_length = QAWS_LITERAL(0.5); /* Force subdivision on flat surface */

		s = qaws_surface_compute_geodesic_distance(surf,
			(qaws_scalar)0.5, (qaws_scalar)0.5, &td,
			u_buf, v_buf, dist, cap, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count > 0, "heat method returned vertices");
		printf("  heat flat: %u vertices\n", count);
	}

	if (count > 0)
	{
		/* Distance at source should be ~0 */
		{
			unsigned int i;
			qaws_scalar min_dist = dist[0];
			unsigned int min_idx = 0;
			for (i = 1; i < count; i++)
			{
				if (dist[i] < min_dist)
				{
					min_dist = dist[i];
					min_idx = i;
				}
			}
			printf("  min_dist=%.6f at vertex %u (u=%.3f, v=%.3f)\n",
				(double)min_dist, min_idx,
				(double)u_buf[min_idx], (double)v_buf[min_idx]);
			TEST_ASSERT(min_dist < (qaws_scalar)0.01,
				"minimum distance is ~0 (at source)");
		}

		/* All distances should be non-negative */
		{
			unsigned int i;
			int all_nonneg = 1;
			for (i = 0; i < count; i++)
			{
				if (dist[i] < -QAWS_LITERAL(1e-6))
				{
					all_nonneg = 0;
					break;
				}
			}
			TEST_ASSERT(all_nonneg, "all distances non-negative");
		}

		/* Max distance should be reasonable (surface is 2x2, max distance ~ 2*sqrt(2) ~ 2.83) */
		{
			unsigned int i;
			qaws_scalar max_dist = 0;
			for (i = 0; i < count; i++)
			{
				if (dist[i] > max_dist)
					max_dist = dist[i];
			}
			printf("  max_dist=%.4f\n", (double)max_dist);
			TEST_ASSERT(max_dist > QAWS_LITERAL(0.1), "max distance > 0");
			TEST_ASSERT(max_dist < QAWS_LITERAL(10.0), "max distance reasonable");
		}
	}

	free(dist);
	free(u_buf);
	free(v_buf);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: heat method distance field on dome surface                   */
/* ------------------------------------------------------------------ */
static void test_heat_distance_dome(void)
{
	qaws_surface* surf = NULL;
	qaws_status s;
	qaws_scalar* dist = NULL;
	qaws_scalar* u_buf = NULL;
	qaws_scalar* v_buf = NULL;
	unsigned int count = 0;
	unsigned int cap = 65536;

	printf("test_heat_distance_dome\n");

	/* Biquadratic dome */
	{
		qaws_surface_biquadratic_desc bd;
		bd.control_points[0].x = 0; bd.control_points[0].y = 0; bd.control_points[0].z = 0;
		bd.control_points[1].x = 1; bd.control_points[1].y = 0; bd.control_points[1].z = 0;
		bd.control_points[2].x = 2; bd.control_points[2].y = 0; bd.control_points[2].z = 0;
		bd.control_points[3].x = 0; bd.control_points[3].y = 1; bd.control_points[3].z = 0;
		bd.control_points[4].x = 1; bd.control_points[4].y = 1; bd.control_points[4].z = 1;
		bd.control_points[5].x = 2; bd.control_points[5].y = 1; bd.control_points[5].z = 0;
		bd.control_points[6].x = 0; bd.control_points[6].y = 2; bd.control_points[6].z = 0;
		bd.control_points[7].x = 1; bd.control_points[7].y = 2; bd.control_points[7].z = 0;
		bd.control_points[8].x = 2; bd.control_points[8].y = 2; bd.control_points[8].z = 0;
		s = qaws_surface_create_biquadratic(&bd, &surf);
		TEST_ASSERT_STATUS(s);
	}

	dist = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));
	u_buf = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));
	v_buf = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));

	{
		qaws_tessellation_desc td;
		td.max_depth = 5;
		td.curvature_threshold = QAWS_LITERAL(0.05);
		td.max_edge_length = QAWS_ZERO;

		s = qaws_surface_compute_geodesic_distance(surf,
			(qaws_scalar)0.5, (qaws_scalar)0.5, &td,
			u_buf, v_buf, dist, cap, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count > 0, "heat method on dome returned vertices");
		printf("  heat dome: %u vertices\n", count);
	}

	if (count > 0)
	{
		/* Distance at source should be ~0 */
		{
			unsigned int i;
			qaws_scalar min_dist = dist[0];
			for (i = 1; i < count; i++)
			{
				if (dist[i] < min_dist)
					min_dist = dist[i];
			}
			TEST_ASSERT(min_dist < (qaws_scalar)0.05,
				"dome: minimum distance is ~0");
		}

		/* Max distance should be positive */
		{
			unsigned int i;
			qaws_scalar max_dist = 0;
			for (i = 0; i < count; i++)
			{
				if (dist[i] > max_dist)
					max_dist = dist[i];
			}
			printf("  dome max_dist=%.4f\n", (double)max_dist);
			TEST_ASSERT(max_dist > QAWS_LITERAL(0.1), "dome max distance > 0");
		}
	}

	free(dist);
	free(u_buf);
	free(v_buf);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: buffer too small for heat method                             */
/* ------------------------------------------------------------------ */
static void test_heat_distance_buffer_check(void)
{
	qaws_surface* surf = NULL;
	qaws_status s;
	qaws_scalar dist_small[4];
	qaws_scalar u_small[4], v_small[4];
	unsigned int count = 0;

	printf("test_heat_distance_buffer_check\n");

	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = 0; bd.p00.y = 0; bd.p00.z = 0;
		bd.p10.x = 1; bd.p10.y = 0; bd.p10.z = 0;
		bd.p01.x = 0; bd.p01.y = 1; bd.p01.z = 0;
		bd.p11.x = 1; bd.p11.y = 1; bd.p11.z = 0;
		s = qaws_surface_create_bilinear(&bd, &surf);
		TEST_ASSERT_STATUS(s);
	}

	/* Capacity too small should return BUFFER_TOO_SMALL and needed count */
	s = qaws_surface_compute_geodesic_distance(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, NULL,
		u_small, v_small, dist_small, 4, &count);
	TEST_ASSERT(s == QAWS_STATUS_BUFFER_TOO_SMALL,
		"buffer too small returns correct status");
	TEST_ASSERT(count > 4, "needed count > capacity");
	printf("  needed count: %u\n", count);

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: 3-way geodesic comparison (RK4, flip, heat) on bumpy       */
/*  terrain.  One OBJ with all three methods.                          */
/* ------------------------------------------------------------------ */
static void visual_obj_geodesic_comparison(void)
{
	qaws_surface* surf = NULL;
	qaws_curve* geo_rk4 = NULL;
	qaws_curve* geo_flip = NULL;
	qaws_curve* geo_heat = NULL;
	qaws_status s;
	qaws_scalar* dist = NULL;
	qaws_scalar* u_buf = NULL;
	qaws_scalar* v_buf = NULL;
	unsigned int count = 0;
	unsigned int cap = 131072;
	qaws_tessellation_desc td_flip, td_heat;

	/* Endpoints: bottom-left to top-right, across all four bumps */
	qaws_scalar su = QAWS_LITERAL(0.05), sv = QAWS_LITERAL(0.05);
	qaws_scalar eu = QAWS_LITERAL(0.95), ev = QAWS_LITERAL(0.95);

	printf("visual_obj_geodesic_comparison\n");

	s = create_hilly_surface(&surf);
	if (s != QAWS_STATUS_OK) return;

	/* Tessellation for flip geodesic (moderate) */
	td_flip.max_depth = 6;
	td_flip.curvature_threshold = QAWS_LITERAL(0.02);
	td_flip.max_edge_length = QAWS_LITERAL(0.3);

	/* Tessellation for heat method (finer for smooth distance field) */
	td_heat.max_depth = 7;
	td_heat.curvature_threshold = QAWS_LITERAL(0.01);
	td_heat.max_edge_length = QAWS_LITERAL(0.15);

	/* Method 1: RK4 shooting (continuous ODE) — more iterations and steps */
	s = qaws_surface_compute_geodesic(surf, su, sv, eu, ev,
		50, 600, &geo_rk4);
	if (s == QAWS_STATUS_OK)
		printf("  RK4 geodesic: OK\n");
	else
		printf("  RK4 geodesic: failed (status %d)\n", (int)s);

	/* Method 2: Dijkstra + UV smoothing (discrete mesh) */
	s = qaws_surface_compute_geodesic_flip(surf, su, sv, eu, ev,
		&td_flip, 0, &geo_flip);
	if (s == QAWS_STATUS_OK)
		printf("  Flip geodesic: OK\n");
	else
		printf("  Flip geodesic: failed (status %d)\n", (int)s);

	/* Method 3: Heat geodesic path (gradient descent on distance field) */
	s = qaws_surface_compute_geodesic_heat(surf, su, sv, eu, ev,
		&td_heat, &geo_heat);
	if (s == QAWS_STATUS_OK)
		printf("  Heat geodesic: OK\n");
	else
		printf("  Heat geodesic: failed (status %d)\n", (int)s);

	/* Heat distance field from start (fine tessellation) */
	dist = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));
	u_buf = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));
	v_buf = (qaws_scalar*)malloc(cap * sizeof(qaws_scalar));

	s = qaws_surface_compute_geodesic_distance(surf, su, sv, &td_heat,
		u_buf, v_buf, dist, cap, &count);
	if (s == QAWS_STATUS_OK)
		printf("  Heat distance: %u vertices\n", count);
	else
		printf("  Heat distance: failed (status %d)\n", (int)s);

	/* Endpoint verification */
	{
		qaws_surface_eval_result sr_start, sr_end;
		qaws_eval_result_3d tip;
		qaws_range rng;
		double dx, dy, dz;

		memset(&sr_start, 0, sizeof(sr_start));
		memset(&sr_end, 0, sizeof(sr_end));
		qaws_surface_evaluate(surf, su, sv, QAWS_SURFACE_EVAL_POSITION, &sr_start);
		qaws_surface_evaluate(surf, eu, ev, QAWS_SURFACE_EVAL_POSITION, &sr_end);
		printf("  Start pos: (%.4f, %.4f, %.4f)\n",
			(double)sr_start.position.x, (double)sr_start.position.y, (double)sr_start.position.z);
		printf("  End   pos: (%.4f, %.4f, %.4f)\n",
			(double)sr_end.position.x, (double)sr_end.position.y, (double)sr_end.position.z);

		if (geo_rk4)
		{
			rng = qaws_curve_get_parameter_range(geo_rk4);
			memset(&tip, 0, sizeof(tip));
			qaws_curve_evaluate_3d(geo_rk4, rng.min_value, QAWS_EVAL_FLAG_POSITION, &tip);
			dx = (double)(tip.position.x - sr_start.position.x);
			dy = (double)(tip.position.y - sr_start.position.y);
			dz = (double)(tip.position.z - sr_start.position.z);
			printf("  RK4  t=0: err_start=%.6f\n", sqrt(dx*dx+dy*dy+dz*dz));
			qaws_curve_evaluate_3d(geo_rk4, rng.max_value, QAWS_EVAL_FLAG_POSITION, &tip);
			dx = (double)(tip.position.x - sr_end.position.x);
			dy = (double)(tip.position.y - sr_end.position.y);
			dz = (double)(tip.position.z - sr_end.position.z);
			printf("  RK4  t=1: err_end  =%.6f\n", sqrt(dx*dx+dy*dy+dz*dz));
		}
		if (geo_flip)
		{
			rng = qaws_curve_get_parameter_range(geo_flip);
			memset(&tip, 0, sizeof(tip));
			qaws_curve_evaluate_3d(geo_flip, rng.min_value, QAWS_EVAL_FLAG_POSITION, &tip);
			dx = (double)(tip.position.x - sr_start.position.x);
			dy = (double)(tip.position.y - sr_start.position.y);
			dz = (double)(tip.position.z - sr_start.position.z);
			printf("  Flip t=0: err_start=%.6f\n", sqrt(dx*dx+dy*dy+dz*dz));
			qaws_curve_evaluate_3d(geo_flip, rng.max_value, QAWS_EVAL_FLAG_POSITION, &tip);
			dx = (double)(tip.position.x - sr_end.position.x);
			dy = (double)(tip.position.y - sr_end.position.y);
			dz = (double)(tip.position.z - sr_end.position.z);
			printf("  Flip t=1: err_end  =%.6f\n", sqrt(dx*dx+dy*dy+dz*dz));
		}
		if (geo_heat)
		{
			rng = qaws_curve_get_parameter_range(geo_heat);
			memset(&tip, 0, sizeof(tip));
			qaws_curve_evaluate_3d(geo_heat, rng.min_value, QAWS_EVAL_FLAG_POSITION, &tip);
			dx = (double)(tip.position.x - sr_start.position.x);
			dy = (double)(tip.position.y - sr_start.position.y);
			dz = (double)(tip.position.z - sr_start.position.z);
			printf("  Heat t=0: err_start=%.6f\n", sqrt(dx*dx+dy*dy+dz*dz));
			qaws_curve_evaluate_3d(geo_heat, rng.max_value, QAWS_EVAL_FLAG_POSITION, &tip);
			dx = (double)(tip.position.x - sr_end.position.x);
			dy = (double)(tip.position.y - sr_end.position.y);
			dz = (double)(tip.position.z - sr_end.position.z);
			printf("  Heat t=1: err_end  =%.6f\n", sqrt(dx*dx+dy*dy+dz*dz));
		}
	}

	/* Print arc lengths for comparison */
	if (geo_rk4)
	{
		qaws_scalar len = 0;
		qaws_range rng = qaws_curve_get_parameter_range(geo_rk4);
		qaws_curve_compute_arc_length(geo_rk4, rng.min_value, rng.max_value, &len);
		printf("  RK4 arc length:  %.4f\n", (double)len);
	}
	if (geo_flip)
	{
		qaws_scalar len = 0;
		qaws_range rng = qaws_curve_get_parameter_range(geo_flip);
		qaws_curve_compute_arc_length(geo_flip, rng.min_value, rng.max_value, &len);
		printf("  Flip arc length: %.4f\n", (double)len);
	}
	if (geo_heat)
	{
		qaws_scalar len = 0;
		qaws_range rng = qaws_curve_get_parameter_range(geo_heat);
		qaws_curve_compute_arc_length(geo_heat, rng.min_value, rng.max_value, &len);
		printf("  Heat arc length: %.4f\n", (double)len);
	}

	/* Write combined OBJ */
	{
		obj_writer obj;
		svg_ensure_output_dir();
		if (obj_open(&obj, OBJ_OUTPUT_DIR "/44_geodesic_comparison.obj",
			OBJ_OUTPUT_DIR "/44_geodesic_comparison.mtl"))
		{
			obj_comment(&obj, "Geodesic algorithm comparison on multi-bump terrain");
			obj_comment(&obj, "Surface: Heat distance field (blue=near, red=far)");
			obj_comment(&obj, "Red tube:   RK4 shooting method (continuous ODE)");
			obj_comment(&obj, "Green tube: Dijkstra + edge flip (discrete mesh)");
			obj_comment(&obj, "Blue tube:  Heat method + gradient descent");
			obj_comment(&obj, "Green sphere: start   Red sphere: end");

			/* Heat distance-colored surface on a regular grid.
			   This avoids tessellation-dependent faceting artifacts. */
			if (count > 0)
			{
				obj_group(&obj, "heat_distance_surface");
				obj_heat_surface_grid(&obj, surf,
					u_buf, v_buf, dist, count, 80);
			}
			else
			{
				obj_material(&obj, "terrain", 0.6, 0.6, 0.5);
				obj_group(&obj, "terrain");
				obj_use_material(&obj, "terrain");
				obj_surface_mesh(&obj, surf, 80, 80);
			}

			/* RK4 geodesic — red */
			if (geo_rk4)
			{
				obj_material(&obj, "rk4_path", 1.0, 0.1, 0.1);
				obj_group(&obj, "geodesic_rk4");
				obj_use_material(&obj, "rk4_path");
				obj_tube_on_surface(&obj, geo_rk4, surf,
					200, 8, (qaws_scalar)0.06);
			}

			/* Flip geodesic — green */
			if (geo_flip)
			{
				obj_material(&obj, "flip_path", 0.1, 0.9, 0.1);
				obj_group(&obj, "geodesic_flip");
				obj_use_material(&obj, "flip_path");
				obj_tube_on_surface(&obj, geo_flip, surf,
					200, 8, (qaws_scalar)0.06);
			}

			/* Heat geodesic — blue */
			if (geo_heat)
			{
				obj_material(&obj, "heat_path", 0.1, 0.3, 1.0);
				obj_group(&obj, "geodesic_heat");
				obj_use_material(&obj, "heat_path");
				obj_tube_on_surface(&obj, geo_heat, surf,
					200, 8, (qaws_scalar)0.06);
			}

			/* Start marker — green */
			{
				qaws_surface_eval_result sr;
				memset(&sr, 0, sizeof(sr));
				qaws_surface_evaluate(surf, su, sv,
					QAWS_SURFACE_EVAL_POSITION, &sr);
				obj_material(&obj, "start_pt", 0.0, 1.0, 0.0);
				obj_group(&obj, "start_marker");
				obj_use_material(&obj, "start_pt");
				obj_sphere(&obj, sr.position, (qaws_scalar)0.12);
			}

			/* End marker — red */
			{
				qaws_surface_eval_result sr;
				memset(&sr, 0, sizeof(sr));
				qaws_surface_evaluate(surf, eu, ev,
					QAWS_SURFACE_EVAL_POSITION, &sr);
				obj_material(&obj, "end_pt", 1.0, 0.0, 0.0);
				obj_group(&obj, "end_marker");
				obj_use_material(&obj, "end_pt");
				obj_sphere(&obj, sr.position, (qaws_scalar)0.12);
			}

			obj_close(&obj);
			printf("  -> " OBJ_OUTPUT_DIR "/44_geodesic_comparison.obj\n");
		}
	}

	free(dist);
	free(u_buf);
	free(v_buf);
	qaws_curve_destroy(geo_rk4);
	qaws_curve_destroy(geo_flip);
	qaws_curve_destroy(geo_heat);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_44_geodesic_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_geodesic_flip_flat();
	test_geodesic_flip_dome();
	test_heat_distance_flat();
	test_heat_distance_dome();
	test_heat_distance_buffer_check();
	visual_obj_geodesic_comparison();

	printf("PASS: %d, FAIL: %d\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
