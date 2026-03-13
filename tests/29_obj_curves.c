#include "test_common.h"

static void test_obj_helix(void)
{
	/* Hermite helix: spirals upward in 3D */
	obj_writer w;
	unsigned int fi;

	printf("test_obj_helix\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/helix.obj",
		OBJ_OUTPUT_DIR "/helix.mtl"))
		return;

	obj_material(&w, "curve_line", 1.0, 0.9, 0.3);
	obj_material(&w, "tube", 0.2, 0.6, 1.0);
	obj_material(&w, "frame_t", 1.0, 0.2, 0.2);
	obj_material(&w, "frame_n", 0.2, 1.0, 0.2);
	obj_material(&w, "frame_b", 0.2, 0.2, 1.0);
	obj_material(&w, "control_pts", 1.0, 0.4, 0.7);

	{
		/* Build a helix from Hermite segments: 8 points going up while rotating */
		unsigned int num_pts = 9;
		qaws_scalar pts[27];  /* 9 * 3 */
		qaws_scalar der[27];  /* 9 * 3 */
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		unsigned int pi;
		qaws_range range;

		for (pi = 0; pi < num_pts; pi++)
		{
			double angle = 2.0 * M_PI * (double)pi / 4.0;  /* 2.25 full turns */
			double r_val = 2.0;
			pts[pi * 3 + 0] = (qaws_scalar)(r_val * cos(angle));
			pts[pi * 3 + 1] = (qaws_scalar)((double)pi * 0.5);
			pts[pi * 3 + 2] = (qaws_scalar)(r_val * sin(angle));
			/* Derivative: tangent to helix */
			der[pi * 3 + 0] = (qaws_scalar)(-r_val * sin(angle) * 1.5);
			der[pi * 3 + 1] = (qaws_scalar)0.5;
			der[pi * 3 + 2] = (qaws_scalar)(r_val * cos(angle) * 1.5);
		}

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = der;
		desc.point_count = num_pts;
		desc.derivative_count = num_pts;

		if (qaws_curve_create_hermite(&desc, &curve) == QAWS_STATUS_OK)
		{
			/* Tube mesh */
			obj_group(&w, "helix_tube");
			obj_use_material(&w, "tube");
			obj_tube(&w, curve, 128, 12, (qaws_scalar)0.08);

			/* Wireframe polyline */
			obj_group(&w, "helix_wire");
			obj_use_material(&w, "curve_line");
			obj_curve_polyline(&w, curve, 256);

			/* Frenet frames every 45 degrees */
			obj_comment(&w, "Frenet frames along helix");
			range = qaws_curve_get_parameter_range(curve);
			for (fi = 0; fi < 16; fi++)
			{
				qaws_scalar t = range.min_value + (range.max_value - range.min_value)
					* (qaws_scalar)fi / (qaws_scalar)15;
				qaws_vec3 T_vec, N_vec, B_vec;
				qaws_eval_result_3d r;
				if (qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION, &r)
					== QAWS_STATUS_OK &&
					qaws_curve_compute_frenet_frame_3d(curve, t, &T_vec, &N_vec, &B_vec)
					== QAWS_STATUS_OK)
				{
					char grp[64];
					sprintf(grp, "frame_%02u_T", fi);
					obj_group(&w, grp);
					obj_use_material(&w, "frame_t");
					{
						qaws_vec3 tip;
						tip.x = r.position.x + T_vec.x * (qaws_scalar)0.4;
						tip.y = r.position.y + T_vec.y * (qaws_scalar)0.4;
						tip.z = r.position.z + T_vec.z * (qaws_scalar)0.4;
						obj_line_segment(&w, r.position, tip);
					}
					sprintf(grp, "frame_%02u_N", fi);
					obj_group(&w, grp);
					obj_use_material(&w, "frame_n");
					{
						qaws_vec3 tip;
						tip.x = r.position.x + N_vec.x * (qaws_scalar)0.4;
						tip.y = r.position.y + N_vec.y * (qaws_scalar)0.4;
						tip.z = r.position.z + N_vec.z * (qaws_scalar)0.4;
						obj_line_segment(&w, r.position, tip);
					}
					sprintf(grp, "frame_%02u_B", fi);
					obj_group(&w, grp);
					obj_use_material(&w, "frame_b");
					{
						qaws_vec3 tip;
						tip.x = r.position.x + B_vec.x * (qaws_scalar)0.4;
						tip.y = r.position.y + B_vec.y * (qaws_scalar)0.4;
						tip.z = r.position.z + B_vec.z * (qaws_scalar)0.4;
						obj_line_segment(&w, r.position, tip);
					}
				}
			}

			/* Control points as markers */
			obj_group(&w, "control_points");
			obj_use_material(&w, "control_pts");
			for (pi = 0; pi < num_pts; pi++)
			{
				qaws_vec3 cp;
				cp.x = pts[pi * 3 + 0];
				cp.y = pts[pi * 3 + 1];
				cp.z = pts[pi * 3 + 2];
				obj_marker(&w, cp, (qaws_scalar)0.06);
			}

			qaws_curve_destroy(curve);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/helix.obj\n");
}

static void test_obj_trefoil_knot(void)
{
	/* Trefoil knot: (sin t + 2 sin 2t, cos t - 2 cos 2t, -sin 3t) */
	obj_writer w;
	unsigned int pi;

	printf("test_obj_trefoil_knot\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/trefoil_knot.obj",
		OBJ_OUTPUT_DIR "/trefoil_knot.mtl"))
		return;

	obj_material(&w, "knot_tube", 0.9, 0.3, 0.5);
	obj_material(&w, "knot_wire", 1.0, 1.0, 0.3);
	obj_material(&w, "comb", 0.3, 1.0, 0.6);
	obj_material(&w, "frame_t", 1.0, 0.2, 0.2);
	obj_material(&w, "frame_n", 0.2, 1.0, 0.2);
	obj_material(&w, "frame_b", 0.2, 0.2, 1.0);

	{
		/* Build trefoil from Hermite: sample the parametric equation */
		unsigned int num_pts = 25;  /* closed: first == last */
		qaws_scalar pts[75];  /* 25 * 3 */
		qaws_scalar der[75];
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		/* Scale derivatives from analytical param [0,2pi] to Hermite param [0,num_pts-1] */
		double dt_scale = 2.0 * M_PI / (double)(num_pts - 1);

		for (pi = 0; pi < num_pts; pi++)
		{
			double t = 2.0 * M_PI * (double)pi / (double)(num_pts - 1);
			pts[pi * 3 + 0] = (qaws_scalar)(sin(t) + 2.0 * sin(2.0 * t));
			pts[pi * 3 + 1] = (qaws_scalar)(cos(t) - 2.0 * cos(2.0 * t));
			pts[pi * 3 + 2] = (qaws_scalar)(-sin(3.0 * t));
			/* Derivatives scaled to Hermite parameterization */
			der[pi * 3 + 0] = (qaws_scalar)((cos(t) + 4.0 * cos(2.0 * t)) * dt_scale);
			der[pi * 3 + 1] = (qaws_scalar)((-sin(t) + 4.0 * sin(2.0 * t)) * dt_scale);
			der[pi * 3 + 2] = (qaws_scalar)((-3.0 * cos(3.0 * t)) * dt_scale);
		}

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = der;
		desc.point_count = num_pts;
		desc.derivative_count = num_pts;

		if (qaws_curve_create_hermite(&desc, &curve) == QAWS_STATUS_OK)
		{
			qaws_range range = qaws_curve_get_parameter_range(curve);
			unsigned int fi;

			/* Tube mesh */
			obj_group(&w, "knot_tube");
			obj_use_material(&w, "knot_tube");
			obj_tube(&w, curve, 200, 16, (qaws_scalar)0.15);

			/* Wireframe */
			obj_group(&w, "knot_wire");
			obj_use_material(&w, "knot_wire");
			obj_curve_polyline(&w, curve, 400);

			/* Curvature comb */
			obj_group(&w, "curvature_comb");
			obj_use_material(&w, "comb");
			obj_curvature_comb(&w, curve, 80, (qaws_scalar)0.5);

			/* Frenet frames at 8 positions */
			for (fi = 0; fi < 8; fi++)
			{
				qaws_scalar t = range.min_value + (range.max_value - range.min_value)
					* (qaws_scalar)fi / (qaws_scalar)8;
				qaws_vec3 T_vec, N_vec, B_vec;
				qaws_eval_result_3d r;
				if (qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION, &r)
					== QAWS_STATUS_OK &&
					qaws_curve_compute_frenet_frame_3d(curve, t, &T_vec, &N_vec, &B_vec)
					== QAWS_STATUS_OK)
				{
					char grp[64];
					sprintf(grp, "frame_%u_T", fi);
					obj_group(&w, grp); obj_use_material(&w, "frame_t");
					{ qaws_vec3 tip = { r.position.x + T_vec.x * 0.6f,
						r.position.y + T_vec.y * 0.6f,
						r.position.z + T_vec.z * 0.6f };
					obj_line_segment(&w, r.position, tip); }
					sprintf(grp, "frame_%u_N", fi);
					obj_group(&w, grp); obj_use_material(&w, "frame_n");
					{ qaws_vec3 tip = { r.position.x + N_vec.x * 0.6f,
						r.position.y + N_vec.y * 0.6f,
						r.position.z + N_vec.z * 0.6f };
					obj_line_segment(&w, r.position, tip); }
					sprintf(grp, "frame_%u_B", fi);
					obj_group(&w, grp); obj_use_material(&w, "frame_b");
					{ qaws_vec3 tip = { r.position.x + B_vec.x * 0.6f,
						r.position.y + B_vec.y * 0.6f,
						r.position.z + B_vec.z * 0.6f };
					obj_line_segment(&w, r.position, tip); }
				}
			}

			qaws_curve_destroy(curve);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/trefoil_knot.obj\n");
}

static void test_obj_3d_intersection(void)
{
	/* Two 3D curves that intersect, with intersection points marked */
	obj_writer w;

	printf("test_obj_3d_intersection\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/intersection_3d.obj",
		OBJ_OUTPUT_DIR "/intersection_3d.mtl"))
		return;

	obj_material(&w, "tube_a", 0.15, 0.5, 0.8);
	obj_material(&w, "tube_b", 0.8, 0.4, 0.15);
	obj_material(&w, "tube_c", 0.2, 0.7, 0.3);
	obj_material(&w, "tube_d", 0.7, 0.2, 0.7);
	obj_material(&w, "isect_pt", 1.0, 0.1, 0.1);
	obj_material(&w, "self_isect_pt", 1.0, 0.9, 0.1);

	/* --- Example 1: Self-intersecting 3D figure-8 Hermite curve --- */
	{
		/* A figure-8 that lies mostly in XY but dips in Z at the crossing */
		qaws_scalar h_pts[] = {
			0, 0, 0,   2, 2, 1,   0, 0, 0,   -2, 2, -1,   0, 0, 0
		};
		qaws_scalar h_der[] = {
			2, 2, 1,   0, -2, 0,   -2, -2, -1,   0, 2, 0,   2, 2, 1
		};
		qaws_hermite_desc hd;
		qaws_curve* crv = NULL;

		memset(&hd, 0, sizeof(hd));
		hd.dimension = QAWS_DIMENSION_3D;
		hd.degree = 3;
		hd.points = h_pts;
		hd.derivatives = h_der;
		hd.point_count = 5;
		hd.derivative_count = 5;

		if (qaws_curve_create_hermite(&hd, &crv) == QAWS_STATUS_OK)
		{
			qaws_intersection_3d isects[32];
			unsigned int ic = 0, ii;

			obj_comment(&w, "Example 1: Self-intersecting 3D figure-8");
			obj_group(&w, "ex1_tube");
			obj_use_material(&w, "tube_a");
			obj_tube(&w, crv, 100, 10, (qaws_scalar)0.05);

			if (qaws_curve_find_self_intersections_3d(crv, isects, 32, &ic)
				== QAWS_STATUS_OK && ic > 0)
			{
				obj_group(&w, "ex1_self_isect");
				obj_use_material(&w, "self_isect_pt");
				for (ii = 0; ii < ic; ii++)
					obj_sphere(&w, isects[ii].position, (qaws_scalar)0.15);
			}
			{ char buf[64]; sprintf(buf, "Self-intersection: %u point(s)", ic); obj_comment(&w, buf); }
			qaws_curve_destroy(crv);
		}
	}

	/* --- Example 2: Two Bezier curves sharing a point at (0,2,0) --- */
	{
		/* Curve A: from (-3,0,0) through (0,2,0) at t=0.5 to (3,0,0)
		   For a cubic Bezier, C(0.5) = (P0+3P1+3P2+P3)/8.
		   We want C(0.5) = (0,2,0). Choose symmetric CPs:
		   P0=(-3,0,0), P1=(-1,4,0), P2=(1,4,0), P3=(3,0,0)
		   => C(0.5)=(-3+3*(-1)+3*1+3)/8 = 0, (0+12+12+0)/8 = 3 ... not 2.
		   Just pick two Hermite segments that both interpolate (0,2,0). */
		qaws_scalar pa[] = { -3,0,0,  0,2,0,  3,0,0 };
		qaws_scalar da_d[] = { 2,2,0,  2,-2,0,  2,-2,0 };
		qaws_hermite_desc hda;
		qaws_curve* ca = NULL;

		/* Curve B: from (0,0,-3) through (0,2,0) to (0,0,3) — perpendicular path */
		qaws_scalar pb[] = { 0,0,-3,  0,2,0,  0,0,3 };
		qaws_scalar db_d[] = { 0,2,2,  0,-2,2,  0,-2,2 };
		qaws_hermite_desc hdb;
		qaws_curve* cb = NULL;

		memset(&hda, 0, sizeof(hda));
		hda.dimension = QAWS_DIMENSION_3D; hda.degree = 3;
		hda.points = pa; hda.derivatives = da_d;
		hda.point_count = 3; hda.derivative_count = 3;
		memset(&hdb, 0, sizeof(hdb));
		hdb.dimension = QAWS_DIMENSION_3D; hdb.degree = 3;
		hdb.points = pb; hdb.derivatives = db_d;
		hdb.point_count = 3; hdb.derivative_count = 3;

		if (qaws_curve_create_hermite(&hda, &ca) == QAWS_STATUS_OK &&
			qaws_curve_create_hermite(&hdb, &cb) == QAWS_STATUS_OK)
		{
			qaws_intersection_3d isects[32];
			unsigned int ic = 0, ii;

			obj_comment(&w, "Example 2: Two Hermite arches crossing at (0,2,0)");
			obj_group(&w, "ex2_arch_a");
			obj_use_material(&w, "tube_a");
			obj_tube(&w, ca, 80, 10, (qaws_scalar)0.05);
			obj_group(&w, "ex2_arch_b");
			obj_use_material(&w, "tube_b");
			obj_tube(&w, cb, 80, 10, (qaws_scalar)0.05);

			if (qaws_curve_find_intersections_3d(ca, cb, isects, 32, &ic)
				== QAWS_STATUS_OK && ic > 0)
			{
				obj_group(&w, "ex2_isect");
				obj_use_material(&w, "isect_pt");
				for (ii = 0; ii < ic; ii++)
					obj_sphere(&w, isects[ii].position, (qaws_scalar)0.15);
			}
			{ char buf[64]; sprintf(buf, "Arches crossing: %u point(s)", ic); obj_comment(&w, buf); }
		}
		if (ca) qaws_curve_destroy(ca);
		if (cb) qaws_curve_destroy(cb);
	}

	/* --- Example 3: Hermite S-curve vs Catmull-Rom wave in 3D --- */
	{
		/* S-curve along X passing through (0,0,0) */
		qaws_scalar s_pts[] = { -4,0,0, 0,0,0, 4,0,0 };
		qaws_scalar s_der[] = { 2,3,2, 2,-3,-2, 2,3,2 };
		qaws_hermite_desc sd;
		qaws_curve* sc = NULL;

		/* Catmull-Rom wave perpendicular, crossing S-curve */
		qaws_vec3 cr_pts[] = { {0,-4,-2}, {0,-1,1}, {0,1,-1}, {0,4,2} };
		qaws_catmull_rom_desc cd;
		qaws_curve* cc = NULL;

		memset(&sd, 0, sizeof(sd));
		sd.dimension = QAWS_DIMENSION_3D; sd.degree = 3;
		sd.points = s_pts; sd.derivatives = s_der;
		sd.point_count = 3; sd.derivative_count = 3;

		memset(&cd, 0, sizeof(cd));
		cd.dimension = QAWS_DIMENSION_3D;
		cd.control_points = cr_pts; cd.control_point_count = 4;
		cd.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		cd.closed = 0;

		if (qaws_curve_create_hermite(&sd, &sc) == QAWS_STATUS_OK &&
			qaws_curve_create_catmull_rom(&cd, &cc) == QAWS_STATUS_OK)
		{
			qaws_intersection_3d isects[32];
			unsigned int ic = 0, ii;

			obj_comment(&w, "Example 3: Hermite S-curve vs Catmull-Rom wave");
			obj_group(&w, "ex3_s_curve");
			obj_use_material(&w, "tube_c");
			obj_tube(&w, sc, 80, 10, (qaws_scalar)0.05);
			obj_group(&w, "ex3_cr_wave");
			obj_use_material(&w, "tube_d");
			obj_tube(&w, cc, 80, 10, (qaws_scalar)0.05);

			if (qaws_curve_find_intersections_3d(sc, cc, isects, 32, &ic)
				== QAWS_STATUS_OK && ic > 0)
			{
				obj_group(&w, "ex3_isect");
				obj_use_material(&w, "isect_pt");
				for (ii = 0; ii < ic; ii++)
					obj_sphere(&w, isects[ii].position, (qaws_scalar)0.15);
			}
			{ char buf[64]; sprintf(buf, "S-curve vs CR wave: %u point(s)", ic); obj_comment(&w, buf); }
		}
		if (sc) qaws_curve_destroy(sc);
		if (cc) qaws_curve_destroy(cc);
	}

	/* --- Example 4: Two Catmull-Rom splines sharing a point at (1,1,1) --- */
	{
		/* CR A: horizontal wave through (1,1,1) */
		qaws_vec3 ca_pts[] = { {-3,1,1}, {-1,0,0}, {1,1,1}, {3,2,0}, {5,1,1} };
		qaws_catmull_rom_desc cda;
		qaws_curve* ca = NULL;

		/* CR B: vertical wave through (1,1,1) */
		qaws_vec3 cb_pts[] = { {1,-3,1}, {0,-1,0}, {1,1,1}, {2,3,2}, {1,5,1} };
		qaws_catmull_rom_desc cdb;
		qaws_curve* cb = NULL;

		memset(&cda, 0, sizeof(cda));
		cda.dimension = QAWS_DIMENSION_3D;
		cda.control_points = ca_pts; cda.control_point_count = 5;
		cda.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		memset(&cdb, 0, sizeof(cdb));
		cdb.dimension = QAWS_DIMENSION_3D;
		cdb.control_points = cb_pts; cdb.control_point_count = 5;
		cdb.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;

		if (qaws_curve_create_catmull_rom(&cda, &ca) == QAWS_STATUS_OK &&
			qaws_curve_create_catmull_rom(&cdb, &cb) == QAWS_STATUS_OK)
		{
			qaws_intersection_3d isects[32];
			unsigned int ic = 0, ii;

			obj_comment(&w, "Example 4: Two Catmull-Rom splines crossing at (1,1,1)");
			obj_group(&w, "ex4_cr_a");
			obj_use_material(&w, "tube_c");
			obj_tube(&w, ca, 80, 10, (qaws_scalar)0.05);
			obj_group(&w, "ex4_cr_b");
			obj_use_material(&w, "tube_d");
			obj_tube(&w, cb, 80, 10, (qaws_scalar)0.05);

			if (qaws_curve_find_intersections_3d(ca, cb, isects, 32, &ic)
				== QAWS_STATUS_OK && ic > 0)
			{
				obj_group(&w, "ex4_isect");
				obj_use_material(&w, "isect_pt");
				for (ii = 0; ii < ic; ii++)
					obj_sphere(&w, isects[ii].position, (qaws_scalar)0.15);
			}
			{ char buf[64]; sprintf(buf, "Catmull-Roms crossing: %u point(s)", ic); obj_comment(&w, buf); }
		}
		if (ca) qaws_curve_destroy(ca);
		if (cb) qaws_curve_destroy(cb);
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/intersection_3d.obj\n");
}

static void test_obj_3d_families(void)
{
	/* Multiple 3D curve families rendered together */
	obj_writer w;

	printf("test_obj_3d_families\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/families_3d.obj",
		OBJ_OUTPUT_DIR "/families_3d.mtl"))
		return;

	obj_material(&w, "bezier", 0.9, 0.3, 0.4);
	obj_material(&w, "bspline", 0.3, 0.9, 0.4);
	obj_material(&w, "catmull", 0.3, 0.4, 0.9);
	obj_material(&w, "nurbs", 0.9, 0.7, 0.2);
	obj_material(&w, "hermite", 0.7, 0.2, 0.9);
	obj_material(&w, "ctrl_pts", 0.5, 0.5, 0.5);

	/* 1: Bezier - 3D S-curve offset in Z */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,3,1}, {3,-1,2}, {4,0,3} };
		qaws_bezier_desc desc;
		qaws_curve* crv = NULL;
		unsigned int pi;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;
		if (qaws_curve_create_bezier(&desc, &crv) == QAWS_STATUS_OK)
		{
			obj_group(&w, "bezier_tube");
			obj_use_material(&w, "bezier");
			obj_tube(&w, crv, 64, 8, (qaws_scalar)0.04);
			obj_group(&w, "bezier_ctrl");
			obj_use_material(&w, "ctrl_pts");
			for (pi = 0; pi < 4; pi++)
				obj_marker(&w, pts[pi], (qaws_scalar)0.05);
			qaws_curve_destroy(crv);
		}
	}

	/* 2: B-spline - wavy 3D path offset in X */
	{
		qaws_vec3 pts[] = { {5,0,0}, {6,2,1}, {7,-1,2}, {8,1,3}, {9,0,1} };
		qaws_bspline_desc desc;
		qaws_curve* crv = NULL;
		unsigned int pi;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_uniform = 1;
		desc.is_closed = 0;
		if (qaws_curve_create_bspline(&desc, &crv) == QAWS_STATUS_OK)
		{
			obj_group(&w, "bspline_tube");
			obj_use_material(&w, "bspline");
			obj_tube(&w, crv, 64, 8, (qaws_scalar)0.04);
			obj_group(&w, "bspline_ctrl");
			obj_use_material(&w, "ctrl_pts");
			for (pi = 0; pi < 5; pi++)
				obj_marker(&w, pts[pi], (qaws_scalar)0.05);
			qaws_curve_destroy(crv);
		}
	}

	/* 3: Catmull-Rom - interpolating path */
	{
		qaws_vec3 pts[] = { {0,0,4}, {1,2,5}, {3,0,6}, {4,1,5}, {5,-1,4} };
		qaws_catmull_rom_desc desc;
		qaws_curve* crv = NULL;
		unsigned int pi;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.control_points = pts;
		desc.control_point_count = 5;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 0;
		if (qaws_curve_create_catmull_rom(&desc, &crv) == QAWS_STATUS_OK)
		{
			obj_group(&w, "catmull_tube");
			obj_use_material(&w, "catmull");
			obj_tube(&w, crv, 64, 8, (qaws_scalar)0.04);
			obj_group(&w, "catmull_ctrl");
			obj_use_material(&w, "ctrl_pts");
			for (pi = 0; pi < 5; pi++)
				obj_marker(&w, pts[pi], (qaws_scalar)0.05);
			qaws_curve_destroy(crv);
		}
	}

	/* 4: NURBS - weighted 3D curve */
	{
		qaws_vec3 pts[] = { {0,0,8}, {1,3,9}, {3,3,9}, {4,0,8} };
		qaws_scalar wts[] = { 1, 3, 3, 1 };
		qaws_nurbs_desc desc;
		qaws_curve* crv = NULL;
		unsigned int pi;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;
		desc.weights = wts;
		desc.weight_count = 4;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_closed = 0;
		if (qaws_curve_create_nurbs(&desc, &crv) == QAWS_STATUS_OK)
		{
			obj_group(&w, "nurbs_tube");
			obj_use_material(&w, "nurbs");
			obj_tube(&w, crv, 64, 8, (qaws_scalar)0.04);
			obj_group(&w, "nurbs_ctrl");
			obj_use_material(&w, "ctrl_pts");
			for (pi = 0; pi < 4; pi++)
				obj_marker(&w, pts[pi], (qaws_scalar)0.05);
			qaws_curve_destroy(crv);
		}
	}

	/* 5: Hermite - 3D curve with explicit derivatives */
	{
		qaws_scalar pts[] = { 5,0,5, 6,2,7, 8,0,6, 9,1,5 };
		qaws_scalar der[] = { 1,2,1, 1,-1,2, 1,1,-1, 1,-2,0 };
		qaws_hermite_desc desc;
		qaws_curve* crv = NULL;
		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = der;
		desc.point_count = 4;
		desc.derivative_count = 4;
		if (qaws_curve_create_hermite(&desc, &crv) == QAWS_STATUS_OK)
		{
			obj_group(&w, "hermite_tube");
			obj_use_material(&w, "hermite");
			obj_tube(&w, crv, 64, 8, (qaws_scalar)0.04);
			qaws_curve_destroy(crv);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/families_3d.obj\n");
}

static void test_obj_3d_curvature_comb(void)
{
	/* 3D B-spline with curvature comb and control polygon */
	obj_writer w;

	printf("test_obj_3d_curvature_comb\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/curvature_comb_3d.obj",
		OBJ_OUTPUT_DIR "/curvature_comb_3d.mtl"))
		return;

	obj_material(&w, "curve", 0.2, 0.7, 1.0);
	obj_material(&w, "comb", 1.0, 0.4, 0.7);
	obj_material(&w, "ctrl", 0.6, 0.6, 0.6);
	obj_material(&w, "polygon", 0.4, 0.4, 0.4);

	{
		qaws_vec3 pts[] = {
			{0,0,0}, {1,3,1}, {2,-1,2}, {3,2,0}, {4,-2,1}, {5,1,3}, {6,0,0}
		};
		qaws_bspline_desc desc;
		qaws_curve* crv = NULL;
		unsigned int pi;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 7;
		desc.knots = NULL;
		desc.knot_count = 0;
		desc.is_uniform = 1;
		desc.is_closed = 0;

		if (qaws_curve_create_bspline(&desc, &crv) == QAWS_STATUS_OK)
		{
			/* Tube */
			obj_group(&w, "curve_tube");
			obj_use_material(&w, "curve");
			obj_tube(&w, crv, 100, 10, (qaws_scalar)0.03);

			/* Curvature comb: ribbon surface + spine markers */
			obj_group(&w, "curvature_comb");
			obj_use_material(&w, "comb");
			obj_curvature_comb_surface(&w, crv, 100, (qaws_scalar)1.0);

			/* Control polygon */
			obj_group(&w, "control_polygon");
			obj_use_material(&w, "polygon");
			{
				unsigned int first_vi = w.vertex_count + 1;
				for (pi = 0; pi < 7; pi++)
					obj_vertex(&w, pts[pi]);
				fprintf(w.fp, "l");
				for (pi = 0; pi < 7; pi++)
					fprintf(w.fp, " %u", first_vi + pi);
				fprintf(w.fp, "\n");
			}

			/* Control point markers */
			obj_group(&w, "control_points");
			obj_use_material(&w, "ctrl");
			for (pi = 0; pi < 7; pi++)
				obj_marker(&w, pts[pi], (qaws_scalar)0.05);

			qaws_curve_destroy(crv);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/curvature_comb_3d.obj\n");
}

static void test_obj_3d_traversal(void)
{
	/* 3D curve with traversal: markers at uniform time intervals showing
	   the difference between parameter-uniform and arc-length-uniform */
	obj_writer w;

	printf("test_obj_3d_traversal\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/traversal_3d.obj",
		OBJ_OUTPUT_DIR "/traversal_3d.mtl"))
		return;

	obj_material(&w, "curve", 0.3, 0.6, 1.0);
	obj_material(&w, "param_pts", 1.0, 0.3, 0.3);
	obj_material(&w, "arc_pts", 0.3, 1.0, 0.3);
	obj_material(&w, "tube", 0.2, 0.4, 0.7);

	{
		/* 3D curve with varying speed: tight curve then straight */
		qaws_vec3 pts[] = {
			{0,0,0}, {0,3,2}, {2,3,3}, {3,0,1}, {6,0,0}
		};
		qaws_bezier_desc desc;
		qaws_curve* crv = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 4;
		desc.control_points = pts;
		desc.control_point_count = 5;

		if (qaws_curve_create_bezier(&desc, &crv) == QAWS_STATUS_OK)
		{
			qaws_range range = qaws_curve_get_parameter_range(crv);
			unsigned int mi;

			/* Tube */
			obj_group(&w, "curve_tube");
			obj_use_material(&w, "tube");
			obj_tube(&w, crv, 80, 10, (qaws_scalar)0.04);

			/* Parameter-uniform markers */
			obj_group(&w, "param_uniform");
			obj_use_material(&w, "param_pts");
			for (mi = 0; mi <= 20; mi++)
			{
				qaws_scalar t = range.min_value + (range.max_value - range.min_value)
					* (qaws_scalar)mi / (qaws_scalar)20;
				qaws_eval_result_3d r;
				if (qaws_curve_evaluate_3d(crv, t, QAWS_EVAL_FLAG_POSITION, &r)
					== QAWS_STATUS_OK)
					obj_marker(&w, r.position, (qaws_scalar)0.04);
			}

			/* Arc-length-uniform markers via traversal */
			{
				qaws_traversal_desc td;
				qaws_traversal* trav = NULL;
				memset(&td, 0, sizeof(td));
				td.traversal_mode = QAWS_TRAVERSAL_MODE_ARC_LENGTH;
				td.speed = (qaws_scalar)1.0;
				td.clamp_to_domain = 1;

				if (qaws_traversal_create(crv, &td, &trav) == QAWS_STATUS_OK)
				{
					qaws_scalar arc_len = 0;
					qaws_curve_compute_arc_length(crv, range.min_value, range.max_value, &arc_len);

					obj_group(&w, "arc_length_uniform");
					obj_use_material(&w, "arc_pts");
					for (mi = 0; mi <= 20; mi++)
					{
						qaws_scalar s = arc_len * (qaws_scalar)mi / (qaws_scalar)20;
						qaws_eval_result_3d r;
						if (qaws_traversal_evaluate_3d(trav, s, QAWS_EVAL_FLAG_POSITION, &r)
							== QAWS_STATUS_OK)
							obj_marker(&w, r.position, (qaws_scalar)0.04);
					}
					qaws_traversal_destroy(trav);
				}
			}

			qaws_curve_destroy(crv);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/traversal_3d.obj\n");
}

static void test_obj_3d_torus_knot(void)
{
	/* Torus knot (p=2, q=3) - classic 3D knot on a torus surface */
	obj_writer w;
	unsigned int pi;

	printf("test_obj_3d_torus_knot\n");

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/torus_knot.obj",
		OBJ_OUTPUT_DIR "/torus_knot.mtl"))
		return;

	obj_material(&w, "knot", 0.1, 0.8, 0.6);
	obj_material(&w, "torus_wire", 0.3, 0.3, 0.3);

	{
		/* Build (2,3)-torus knot from Hermite segments */
		/* x = (R + r*cos(q*t)) * cos(p*t)
		   y = (R + r*cos(q*t)) * sin(p*t)
		   z = r * sin(q*t) */
		double R = 3.0, r_val = 1.0;
		int p_val = 2, q_val = 3;
		unsigned int num_pts = 49;  /* closed */
		qaws_scalar pts[147];  /* 49 * 3 */
		qaws_scalar der[147];
		qaws_hermite_desc desc;
		qaws_curve* curve = NULL;
		/* Hermite spans are unit-width [i, i+1], so derivatives must be
		   dC/dt_hermite = dC/dt_analytical * dt_analytical/dt_hermite.
		   dt_analytical covers [0, 2pi] over (num_pts-1) spans. */
		double dt_scale = 2.0 * M_PI / (double)(num_pts - 1);

		for (pi = 0; pi < num_pts; pi++)
		{
			double t = 2.0 * M_PI * (double)pi / (double)(num_pts - 1);
			double rr = R + r_val * cos((double)q_val * t);
			double drr = -r_val * (double)q_val * sin((double)q_val * t);
			pts[pi * 3 + 0] = (qaws_scalar)(rr * cos((double)p_val * t));
			pts[pi * 3 + 1] = (qaws_scalar)(rr * sin((double)p_val * t));
			pts[pi * 3 + 2] = (qaws_scalar)(r_val * sin((double)q_val * t));
			/* Derivative scaled to Hermite parameterization */
			der[pi * 3 + 0] = (qaws_scalar)((drr * cos((double)p_val * t)
				- rr * (double)p_val * sin((double)p_val * t)) * dt_scale);
			der[pi * 3 + 1] = (qaws_scalar)((drr * sin((double)p_val * t)
				+ rr * (double)p_val * cos((double)p_val * t)) * dt_scale);
			der[pi * 3 + 2] = (qaws_scalar)((r_val * (double)q_val * cos((double)q_val * t)) * dt_scale);
		}

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = der;
		desc.point_count = num_pts;
		desc.derivative_count = num_pts;

		if (qaws_curve_create_hermite(&desc, &curve) == QAWS_STATUS_OK)
		{
			obj_group(&w, "torus_knot_tube");
			obj_use_material(&w, "knot");
			obj_tube(&w, curve, 300, 12, (qaws_scalar)0.12);

			qaws_curve_destroy(curve);
		}

		/* Draw reference torus wireframe (circles at major/minor radii) */
		obj_group(&w, "torus_reference");
		obj_use_material(&w, "torus_wire");
		{
			/* Major circle */
			unsigned int first_vi = w.vertex_count + 1;
			unsigned int ti;
			for (ti = 0; ti <= 64; ti++)
			{
				double angle = 2.0 * M_PI * (double)ti / 64.0;
				qaws_vec3 v;
				v.x = (qaws_scalar)(R * cos(angle));
				v.y = (qaws_scalar)(R * sin(angle));
				v.z = 0;
				obj_vertex(&w, v);
			}
			fprintf(w.fp, "l");
			for (ti = 0; ti <= 64; ti++)
				fprintf(w.fp, " %u", first_vi + ti);
			fprintf(w.fp, "\n");
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/torus_knot.obj\n");
}

int test_29_obj_curves_main(void)
{
	svg_ensure_output_dir();
	test_obj_helix();
	test_obj_trefoil_knot();
	test_obj_3d_intersection();
	test_obj_3d_families();
	test_obj_3d_curvature_comb();
	test_obj_3d_traversal();
	test_obj_3d_torus_knot();
	return 0;
}
