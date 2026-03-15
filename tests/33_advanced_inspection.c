#include "test_common.h"

/* ------------------------------------------------------------------ */
/*  Helper: create a dome surface                                      */
/* ------------------------------------------------------------------ */
static qaws_surface* make_dome(void)
{
	qaws_surface* surf = NULL;
	qaws_surface_bezier_desc desc;
	qaws_vec3 cp[9];

	cp[0].x = 0; cp[0].y = 0; cp[0].z = 0;
	cp[1].x = 0; cp[1].y = 1; cp[1].z = 0;
	cp[2].x = 0; cp[2].y = 2; cp[2].z = 0;
	cp[3].x = 1; cp[3].y = 0; cp[3].z = 0;
	cp[4].x = 1; cp[4].y = 1; cp[4].z = 2;
	cp[5].x = 1; cp[5].y = 2; cp[5].z = 0;
	cp[6].x = 2; cp[6].y = 0; cp[6].z = 0;
	cp[7].x = 2; cp[7].y = 1; cp[7].z = 0;
	cp[8].x = 2; cp[8].y = 2; cp[8].z = 0;

	desc.u_degree = 2; desc.v_degree = 2;
	desc.control_points = cp;
	desc.u_point_count = 3; desc.v_point_count = 3;
	qaws_surface_create_bezier(&desc, &surf);
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Helper: create a flat surface at z=0                               */
/* ------------------------------------------------------------------ */
static qaws_surface* make_flat_surface(void)
{
	qaws_surface* surf = NULL;
	qaws_surface_bezier_desc desc;
	qaws_vec3 cp[4];

	cp[0].x = 0; cp[0].y = 0; cp[0].z = 0;
	cp[1].x = 0; cp[1].y = 4; cp[1].z = 0;
	cp[2].x = 4; cp[2].y = 0; cp[2].z = 0;
	cp[3].x = 4; cp[3].y = 4; cp[3].z = 0;

	desc.u_degree = 1; desc.v_degree = 1;
	desc.control_points = cp;
	desc.u_point_count = 2; desc.v_point_count = 2;
	qaws_surface_create_bezier(&desc, &surf);
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Test: Closest point on surface                                     */
/* ------------------------------------------------------------------ */
static void test_closest_point_on_surface(void)
{
	qaws_surface* surf = NULL;
	qaws_scalar u, v;
	qaws_vec3 closest;
	qaws_status s;

	printf("test_closest_point_on_surface\n");

	surf = make_flat_surface();
	TEST_ASSERT(surf != NULL, "closest pt surface created");

	/* Point directly above center: closest should be center */
	{
		qaws_vec3 pt = {2, 2, 5};
		s = qaws_surface_find_closest_point(surf, pt, &u, &v, &closest);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(approx_eq_loose(closest.x, (qaws_scalar)2.0), "flat closest x");
		TEST_ASSERT(approx_eq_loose(closest.y, (qaws_scalar)2.0), "flat closest y");
		TEST_ASSERT(fabs(closest.z) < TOLERANCE_LOOSE, "flat closest z");
		TEST_ASSERT(approx_eq_loose(u, (qaws_scalar)0.5), "flat closest u");
		TEST_ASSERT(approx_eq_loose(v, (qaws_scalar)0.5), "flat closest v");
	}

	/* Point near a corner */
	{
		qaws_vec3 pt = {0, 0, 1};
		s = qaws_surface_find_closest_point(surf, pt, &u, &v, &closest);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(fabs(closest.x) < TOLERANCE_LOOSE, "corner closest x");
		TEST_ASSERT(fabs(closest.y) < TOLERANCE_LOOSE, "corner closest y");
	}

	qaws_surface_destroy(surf);

	/* Dome surface: point above the peak */
	surf = make_dome();
	TEST_ASSERT(surf != NULL, "dome closest pt surface created");

	{
		qaws_vec3 pt = {1, 1, 10};
		s = qaws_surface_find_closest_point(surf, pt, &u, &v, &closest);
		TEST_ASSERT_STATUS(s);
		/* Should converge near the dome peak (u~0.5, v~0.5) */
		TEST_ASSERT(u > (qaws_scalar)0.2 && u < (qaws_scalar)0.8, "dome closest u mid");
		TEST_ASSERT(v > (qaws_scalar)0.2 && v < (qaws_scalar)0.8, "dome closest v mid");
		TEST_ASSERT(closest.z > (qaws_scalar)0.1, "dome closest z elevated");
	}

	/* Null args */
	{
		qaws_vec3 pt = {0,0,0};
		s = qaws_surface_find_closest_point(NULL, pt, &u, &v, &closest);
		TEST_ASSERT(s != QAWS_STATUS_OK, "closest null surface rejected");
		s = qaws_surface_find_closest_point(surf, pt, NULL, &v, &closest);
		TEST_ASSERT(s != QAWS_STATUS_OK, "closest null u rejected");
	}

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: Curve-plane intersection                                     */
/* ------------------------------------------------------------------ */
static void test_curve_plane_intersection(void)
{
	qaws_curve* crv = NULL;
	qaws_scalar params[16];
	qaws_vec3 positions[16];
	unsigned int count = 0;
	qaws_status s;

	printf("test_curve_plane_intersection\n");

	/* 3D line crossing z=1 plane */
	{
		qaws_scalar pts[] = {0,0,0, 2,2,2};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 1;
		d.control_points = pts;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &crv);
	}

	{
		qaws_plane plane;
		plane.point.x = 0; plane.point.y = 0; plane.point.z = 1;
		plane.normal.x = 0; plane.normal.y = 0; plane.normal.z = 1;
		s = qaws_curve_find_plane_intersections(crv, &plane,
			params, positions, 16, &count);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count == 1, "line-plane: 1 intersection");
	TEST_ASSERT(approx_eq_loose(params[0], (qaws_scalar)0.5), "line-plane t=0.5");
	TEST_ASSERT(approx_eq_loose(positions[0].z, (qaws_scalar)1.0), "line-plane z=1");

	qaws_curve_destroy(crv);

	/* Cubic curve crossing z=0 plane twice */
	{
		qaws_scalar pts[] = {0,0,-1, 1,0,2, 2,0,-2, 3,0,1};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = pts;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &crv);
	}

	{
		qaws_plane plane;
		plane.point.x = 0; plane.point.y = 0; plane.point.z = 0;
		plane.normal.x = 0; plane.normal.y = 0; plane.normal.z = 1;
		count = 0;
		s = qaws_curve_find_plane_intersections(crv, &plane,
			params, positions, 16, &count);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count >= 2, "cubic-plane: >= 2 intersections");

	/* All intersection points should be at z~0 */
	{
		unsigned int i;
		for (i = 0; i < count; i++)
			TEST_ASSERT(fabs(positions[i].z) < TOLERANCE_LOOSE,
				"intersection z~0");
	}

	qaws_curve_destroy(crv);

	/* Curve parallel to plane: no intersection */
	{
		qaws_scalar pts[] = {0,0,5, 3,0,5};
		qaws_bezier_desc d;
		qaws_plane plane;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 1;
		d.control_points = pts;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &crv);

		plane.point.x = 0; plane.point.y = 0; plane.point.z = 0;
		plane.normal.x = 0; plane.normal.y = 0; plane.normal.z = 1;
		count = 0;
		s = qaws_curve_find_plane_intersections(crv, &plane,
			params, positions, 16, &count);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 0, "parallel: 0 intersections");
		qaws_curve_destroy(crv);
	}

	/* Null args */
	{
		qaws_plane plane;
		plane.point.x = 0; plane.point.y = 0; plane.point.z = 0;
		plane.normal.x = 0; plane.normal.y = 0; plane.normal.z = 1;
		s = qaws_curve_find_plane_intersections(NULL, &plane,
			params, positions, 16, &count);
		TEST_ASSERT(s != QAWS_STATUS_OK, "plane null curve rejected");
	}
}

/* ------------------------------------------------------------------ */
/*  Test: Surface-curve intersection                                   */
/* ------------------------------------------------------------------ */
static void test_surface_curve_intersection(void)
{
	qaws_surface* surf = NULL;
	qaws_curve* crv = NULL;
	qaws_surface_curve_intersection hits[16];
	unsigned int count = 0;
	qaws_status s;

	printf("test_surface_curve_intersection\n");

	/* Flat surface at z=0, vertical line from (2,2,-1) to (2,2,1) */
	surf = make_flat_surface();
	TEST_ASSERT(surf != NULL, "surf-curve surface created");

	{
		qaws_scalar pts[] = {2,2,-1, 2,2,1};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 1;
		d.control_points = pts;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &crv);
	}

	s = qaws_surface_find_curve_intersections(surf, crv, hits, 16, &count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count >= 1, "flat-line: >= 1 intersection");

	if (count >= 1)
	{
		TEST_ASSERT(approx_eq_loose(hits[0].position.x, (qaws_scalar)2.0),
			"flat-line hit x=2");
		TEST_ASSERT(approx_eq_loose(hits[0].position.y, (qaws_scalar)2.0),
			"flat-line hit y=2");
		TEST_ASSERT(fabs(hits[0].position.z) < TOLERANCE_LOOSE,
			"flat-line hit z=0");
		/* t should be ~0.5 (midpoint of [-1,1] line) */
		TEST_ASSERT(approx_eq_loose(hits[0].t, (qaws_scalar)0.5),
			"flat-line hit t=0.5");
	}

	qaws_curve_destroy(crv);
	qaws_surface_destroy(surf);

	/* Dome surface with diagonal line */
	surf = make_dome();
	TEST_ASSERT(surf != NULL, "dome surf-curve surface created");

	{
		qaws_scalar pts[] = {1,1,-1, 1,1,5};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 1;
		d.control_points = pts;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &crv);
	}

	count = 0;
	s = qaws_surface_find_curve_intersections(surf, crv, hits, 16, &count);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count >= 1, "dome-line: >= 1 intersection");

	qaws_curve_destroy(crv);
	qaws_surface_destroy(surf);

	/* Null args */
	{
		s = qaws_surface_find_curve_intersections(NULL, crv, hits, 16, &count);
		TEST_ASSERT(s != QAWS_STATUS_OK, "surf-curve null surface rejected");
	}
}

/* ------------------------------------------------------------------ */
/*  Test: Adaptive tessellation                                        */
/* ------------------------------------------------------------------ */
static void test_adaptive_tessellation(void)
{
	qaws_surface* surf = NULL;
	qaws_tessellation_vertex verts[4096];
	unsigned int indices[24576];
	unsigned int vert_count = 0, idx_count = 0;
	qaws_status s;

	printf("test_adaptive_tessellation\n");

	/* Flat surface: should produce few triangles */
	surf = make_flat_surface();
	TEST_ASSERT(surf != NULL, "tess flat surface created");

	{
		qaws_tessellation_desc desc;
		desc.max_depth = 4;
		desc.curvature_threshold = (qaws_scalar)(0.1);
		desc.max_edge_length = 0;  /* disabled */
		s = qaws_surface_tessellate(surf, &desc,
			verts, 4096, &vert_count,
			indices, 24576, &idx_count);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(vert_count > 0, "tess flat has vertices");
	TEST_ASSERT(idx_count > 0, "tess flat has indices");
	TEST_ASSERT(idx_count % 3 == 0, "tess index count multiple of 3");

	/* Flat surface should not subdivide much */
	TEST_ASSERT(vert_count < 100, "tess flat vertex count reasonable");

	qaws_surface_destroy(surf);

	/* Dome surface: should produce more triangles (higher curvature) */
	surf = make_dome();
	TEST_ASSERT(surf != NULL, "tess dome surface created");

	{
		qaws_tessellation_desc desc;
		desc.max_depth = 5;
		desc.curvature_threshold = (qaws_scalar)(0.05);
		desc.max_edge_length = 0;
		vert_count = 0; idx_count = 0;
		s = qaws_surface_tessellate(surf, &desc,
			verts, 4096, &vert_count,
			indices, 24576, &idx_count);
	}
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(vert_count > 0, "tess dome has vertices");
	TEST_ASSERT(idx_count > 0, "tess dome has indices");

	/* All vertices should have valid positions and normals */
	{
		unsigned int i;
		int all_valid = 1;
		for (i = 0; i < vert_count && i < 100; i++)
		{
			qaws_scalar len = (qaws_scalar)sqrt((double)(
				verts[i].normal.x * verts[i].normal.x +
				verts[i].normal.y * verts[i].normal.y +
				verts[i].normal.z * verts[i].normal.z));
			if (len < (qaws_scalar)0.5) { all_valid = 0; break; }
		}
		TEST_ASSERT(all_valid, "tess normals unit length");
	}

	/* All indices should be in range */
	{
		unsigned int i;
		int all_valid = 1;
		for (i = 0; i < idx_count; i++)
		{
			if (indices[i] >= vert_count) { all_valid = 0; break; }
		}
		TEST_ASSERT(all_valid, "tess indices in range");
	}

	/* Edge length constraint */
	{
		qaws_tessellation_desc desc;
		unsigned int vert2 = 0, idx2 = 0;
		desc.max_depth = 5;
		desc.curvature_threshold = (qaws_scalar)10.0;  /* very loose - won't trigger */
		desc.max_edge_length = (qaws_scalar)0.5;       /* force edge subdivision */
		s = qaws_surface_tessellate(surf, &desc,
			verts, 4096, &vert2,
			indices, 24576, &idx2);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(vert2 > 0, "tess edge-len has vertices");
	}

	/* Null args */
	{
		qaws_tessellation_desc desc;
		desc.max_depth = 3;
		desc.curvature_threshold = (qaws_scalar)(0.1);
		desc.max_edge_length = 0;
		s = qaws_surface_tessellate(NULL, &desc,
			verts, 4096, &vert_count,
			indices, 24576, &idx_count);
		TEST_ASSERT(s != QAWS_STATUS_OK, "tess null surface rejected");
	}

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ adaptive tessellation                                  */
/* ------------------------------------------------------------------ */
static void visual_obj_tessellation(void)
{
	qaws_surface* surf = make_dome();
	qaws_tessellation_vertex verts[4096];
	unsigned int indices[24576];
	unsigned int vert_count = 0, idx_count = 0;
	FILE* fp;

	printf("visual_obj_tessellation\n");
	svg_ensure_output_dir();
	if (!surf) return;

	{
		qaws_tessellation_desc desc;
		desc.max_depth = 5;
		desc.curvature_threshold = (qaws_scalar)(0.05);
		desc.max_edge_length = 0;
		qaws_surface_tessellate(surf, &desc,
			verts, 4096, &vert_count,
			indices, 24576, &idx_count);
	}

	fp = fopen(OBJ_OUTPUT_DIR "/33_adaptive_tessellation.obj", "w");
	if (fp)
	{
		unsigned int i;
		fprintf(fp, "# Adaptive tessellation: %u vertices, %u triangles\n",
			vert_count, idx_count / 3);
		for (i = 0; i < vert_count; i++)
			fprintf(fp, "v %.6f %.6f %.6f\n",
				(double)verts[i].position.x,
				(double)verts[i].position.y,
				(double)verts[i].position.z);
		for (i = 0; i < vert_count; i++)
			fprintf(fp, "vn %.6f %.6f %.6f\n",
				(double)verts[i].normal.x,
				(double)verts[i].normal.y,
				(double)verts[i].normal.z);
		for (i = 0; i < idx_count; i += 3)
			fprintf(fp, "f %u//%u %u//%u %u//%u\n",
				indices[i]+1, indices[i]+1,
				indices[i+1]+1, indices[i+1]+1,
				indices[i+2]+1, indices[i+2]+1);
		fclose(fp);
		printf("  -> " OBJ_OUTPUT_DIR "/33_adaptive_tessellation.obj (%u verts, %u tris)\n",
			vert_count, idx_count / 3);
	}

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: SVG curve-plane intersection                               */
/* ------------------------------------------------------------------ */
static void visual_svg_plane_intersection(void)
{
	svg_writer svg;
	qaws_curve* crv = NULL;
	qaws_scalar params[16];
	qaws_vec3 positions[16];
	unsigned int count = 0;
	unsigned int i;

	printf("visual_svg_plane_intersection\n");
	svg_ensure_output_dir();

	/* 3D cubic curve projected to XZ for visualization */
	{
		qaws_scalar pts[] = {0,0,-2, 1,0,3, 3,0,-3, 4,0,2};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = pts;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &crv);
	}
	if (!crv) return;

	{
		qaws_plane plane;
		plane.point.x = 0; plane.point.y = 0; plane.point.z = 0;
		plane.normal.x = 0; plane.normal.y = 0; plane.normal.z = 1;
		qaws_curve_find_plane_intersections(crv, &plane,
			params, positions, 16, &count);
	}

	/* Draw XZ projection */
	if (!svg_open(&svg, OBJ_OUTPUT_DIR "/33_plane_intersection.svg",
		(qaws_scalar)-1, (qaws_scalar)-4,
		(qaws_scalar)6, (qaws_scalar)8,
		(qaws_scalar)600, (qaws_scalar)400))
	{
		qaws_curve_destroy(crv);
		return;
	}

	/* Sample curve and draw as XZ polyline */
	{
		qaws_vec2 pts_2d[128];
		unsigned int n;
		for (n = 0; n < 128; n++)
		{
			qaws_scalar t = (qaws_scalar)n / (qaws_scalar)127;
			qaws_eval_result_3d er;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(crv, t, QAWS_EVAL_FLAG_POSITION, &er);
			pts_2d[n].x = er.position.x;
			pts_2d[n].y = er.position.z;
		}
		svg_polyline(&svg, pts_2d, 128, "#4488ff", 2);
	}

	/* Draw z=0 plane as horizontal line */
	svg_line(&svg, (qaws_scalar)-1, 0, (qaws_scalar)5, 0, "#666666", 1);

	/* Draw intersection points */
	for (i = 0; i < count; i++)
		svg_dot(&svg, positions[i].x, positions[i].z, "#ff4444", 5);

	{
		char buf[64];
		sprintf(buf, "%u intersections", count);
		svg_label(&svg, (qaws_scalar)0.5, (qaws_scalar)3.5, buf, "#ffffff");
	}

	svg_close(&svg);
	qaws_curve_destroy(crv);
	printf("  -> " OBJ_OUTPUT_DIR "/33_plane_intersection.svg (%u hits)\n", count);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ closest point on surface                               */
/* ------------------------------------------------------------------ */
static void visual_obj_closest_point(void)
{
	obj_writer w;
	qaws_surface* surf = make_dome();
	qaws_vec3 queries[4];
	unsigned int i;

	printf("visual_obj_closest_point\n");
	svg_ensure_output_dir();
	if (!surf) return;

	queries[0].x = 1; queries[0].y = 1; queries[0].z = 5;
	queries[1].x = 0; queries[1].y = 0; queries[1].z = 3;
	queries[2].x = 2; queries[2].y = 2; queries[2].z = 2;
	queries[3].x = 1; queries[3].y = (qaws_scalar)0.5; queries[3].z = 4;

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/33_closest_point.obj",
		OBJ_OUTPUT_DIR "/33_closest_point.mtl"))
	{
		qaws_surface_destroy(surf);
		return;
	}

	/* Surface mesh */
	obj_material(&w, "dome", 0.5, 0.7, 0.9);
	obj_group(&w, "dome_surface");
	obj_use_material(&w, "dome");
	obj_surface_mesh(&w, surf, 24, 24);

	/* Query points (red) and closest points (green) with connecting lines */
	obj_material(&w, "query", 1.0, 0.2, 0.2);
	obj_material(&w, "closest", 0.2, 1.0, 0.3);
	obj_material(&w, "connector", 0.8, 0.8, 0.3);

	obj_material(&w, "normal", 0.3, 0.3, 1.0);

	for (i = 0; i < 4; i++)
	{
		qaws_scalar u, v;
		qaws_vec3 closest;
		if (qaws_surface_find_closest_point(surf, queries[i],
			&u, &v, &closest) == QAWS_STATUS_OK)
		{
			unsigned int j;
			qaws_surface_eval_result nr;

			obj_group(&w, "query_pt");
			obj_use_material(&w, "query");
			obj_sphere(&w, queries[i], (qaws_scalar)0.06);

			obj_group(&w, "closest_pt");
			obj_use_material(&w, "closest");
			obj_sphere(&w, closest, (qaws_scalar)0.06);

			/* Connecting rod: chain of small spheres */
			obj_group(&w, "connector");
			obj_use_material(&w, "connector");
			for (j = 1; j < 10; j++)
			{
				qaws_scalar f = (qaws_scalar)j / (qaws_scalar)10;
				qaws_vec3 p;
				p.x = closest.x + f * (queries[i].x - closest.x);
				p.y = closest.y + f * (queries[i].y - closest.y);
				p.z = closest.z + f * (queries[i].z - closest.z);
				obj_sphere(&w, p, (qaws_scalar)0.02);
			}

			/* Surface normal arrow at closest point (shows orthogonality) */
			memset(&nr, 0, sizeof(nr));
			if (qaws_surface_evaluate(surf, u, v,
				QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DU
				| QAWS_SURFACE_EVAL_DV, &nr) == QAWS_STATUS_OK)
			{
				qaws_vec3 tip;
				tip.x = closest.x + nr.normal.x * (qaws_scalar)0.5;
				tip.y = closest.y + nr.normal.y * (qaws_scalar)0.5;
				tip.z = closest.z + nr.normal.z * (qaws_scalar)0.5;
				obj_group(&w, "normal_arrow");
				obj_use_material(&w, "normal");
				obj_sphere(&w, tip, (qaws_scalar)0.03);
				for (j = 0; j < 5; j++)
				{
					qaws_scalar f = (qaws_scalar)j / (qaws_scalar)5;
					qaws_vec3 p;
					p.x = closest.x + f * nr.normal.x * (qaws_scalar)0.5;
					p.y = closest.y + f * nr.normal.y * (qaws_scalar)0.5;
					p.z = closest.z + f * nr.normal.z * (qaws_scalar)0.5;
					obj_sphere(&w, p, (qaws_scalar)0.015);
				}
			}
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/33_closest_point.obj\n");
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ surface-curve intersection                             */
/* ------------------------------------------------------------------ */
static void visual_obj_surface_curve_intersection(void)
{
	obj_writer w;
	qaws_surface* surf = make_dome();
	qaws_curve* crv = NULL;
	qaws_surface_curve_intersection hits[16];
	unsigned int count = 0;
	unsigned int i;

	printf("visual_obj_surface_curve_intersection\n");
	svg_ensure_output_dir();
	if (!surf) return;

	/* Cubic S-curve piercing the dome */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.0,  (qaws_scalar)2.0,  (qaws_scalar)-1.0,
			(qaws_scalar)0.5,  (qaws_scalar)0.0,  (qaws_scalar) 1.5,
			(qaws_scalar)1.5,  (qaws_scalar)2.0,  (qaws_scalar) 0.5,
			(qaws_scalar)2.0,  (qaws_scalar)0.0,  (qaws_scalar) 3.0
		};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = pts;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &crv);
	}
	if (!crv) { qaws_surface_destroy(surf); return; }

	qaws_surface_find_curve_intersections(surf, crv, hits, 16, &count);

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/33_surface_curve_intersection.obj",
		OBJ_OUTPUT_DIR "/33_surface_curve_intersection.mtl"))
	{
		qaws_curve_destroy(crv); qaws_surface_destroy(surf); return;
	}

	/* Surface */
	obj_material(&w, "dome", 0.5, 0.7, 0.9);
	obj_group(&w, "dome_surface");
	obj_use_material(&w, "dome");
	obj_surface_mesh(&w, surf, 24, 24);

	/* Piercing curve as a visible tube */
	obj_material(&w, "curve", 0.9, 0.6, 0.1);
	obj_group(&w, "piercing_curve");
	obj_use_material(&w, "curve");
	obj_tube(&w, crv, 48, 8, (qaws_scalar)0.03);

	/* Curve endpoint markers */
	{
		qaws_eval_result_3d er;
		memset(&er, 0, sizeof(er));
		qaws_curve_evaluate_3d(crv, 0, QAWS_EVAL_FLAG_POSITION, &er);
		obj_sphere(&w, er.position, (qaws_scalar)0.04);
		memset(&er, 0, sizeof(er));
		qaws_curve_evaluate_3d(crv, 1, QAWS_EVAL_FLAG_POSITION, &er);
		obj_sphere(&w, er.position, (qaws_scalar)0.04);
	}

	/* Intersection markers */
	obj_material(&w, "hit", 1.0, 0.15, 0.15);
	for (i = 0; i < count; i++)
	{
		obj_group(&w, "hit_point");
		obj_use_material(&w, "hit");
		obj_sphere(&w, hits[i].position, (qaws_scalar)0.07);
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/33_surface_curve_intersection.obj (%u hits)\n", count);
	qaws_curve_destroy(crv);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: SVG surface-curve intersection (XZ projection)             */
/* ------------------------------------------------------------------ */
static void visual_svg_surface_curve_intersection(void)
{
	svg_writer svg;
	qaws_surface* surf = make_dome();
	qaws_curve* crv = NULL;
	qaws_surface_curve_intersection hits[16];
	unsigned int count = 0;
	unsigned int i, ui, vi;

	/* Oblique projection constants: screen_x = x + y*CX, screen_y = z + y*CY */
	qaws_scalar proj_cx = (qaws_scalar)0.45;
	qaws_scalar proj_cy = (qaws_scalar)0.3;

	printf("visual_svg_surface_curve_intersection\n");
	svg_ensure_output_dir();
	if (!surf) return;

	/* Diagonal S-curve piercing the dome */
	{
		qaws_scalar pts[] = {
			(qaws_scalar)0.2,  (qaws_scalar)1.8,  (qaws_scalar)-0.5,
			(qaws_scalar)0.6,  (qaws_scalar)0.3,  (qaws_scalar) 2.5,
			(qaws_scalar)1.4,  (qaws_scalar)1.7,  (qaws_scalar) 0.5,
			(qaws_scalar)1.8,  (qaws_scalar)0.2,  (qaws_scalar) 3.0
		};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 3;
		d.control_points = pts;
		d.control_point_count = 4;
		qaws_curve_create_bezier(&d, &crv);
	}
	if (!crv) { qaws_surface_destroy(surf); return; }

	qaws_surface_find_curve_intersections(surf, crv, hits, 16, &count);

	if (!svg_open(&svg, OBJ_OUTPUT_DIR "/33_surface_curve_intersection.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.0,
		(qaws_scalar)4.0, (qaws_scalar)5.0,
		(qaws_scalar)550, (qaws_scalar)550))
	{
		qaws_curve_destroy(crv); qaws_surface_destroy(surf);
		return;
	}

	/* Draw surface iso-u lines (oblique projection) */
	for (ui = 0; ui <= 8; ui++)
	{
		qaws_vec2 pts_2d[33];
		qaws_scalar u = (qaws_scalar)ui / (qaws_scalar)8;
		for (vi = 0; vi <= 32; vi++)
		{
			qaws_scalar v = (qaws_scalar)vi / (qaws_scalar)32;
			qaws_surface_eval_result r;
			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surf, u, v, QAWS_SURFACE_EVAL_POSITION, &r);
			pts_2d[vi].x = r.position.x + r.position.y * proj_cx;
			pts_2d[vi].y = r.position.z + r.position.y * proj_cy;
		}
		svg_polyline(&svg, pts_2d, 33, "#5588bb", 1);
	}

	/* Draw surface iso-v lines */
	for (vi = 0; vi <= 8; vi++)
	{
		qaws_vec2 pts_2d[33];
		qaws_scalar v = (qaws_scalar)vi / (qaws_scalar)8;
		for (ui = 0; ui <= 32; ui++)
		{
			qaws_scalar u = (qaws_scalar)ui / (qaws_scalar)32;
			qaws_surface_eval_result r;
			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surf, u, v, QAWS_SURFACE_EVAL_POSITION, &r);
			pts_2d[ui].x = r.position.x + r.position.y * proj_cx;
			pts_2d[ui].y = r.position.z + r.position.y * proj_cy;
		}
		svg_polyline(&svg, pts_2d, 33, "#5588bb", 1);
	}

	/* Draw curve (oblique projection) */
	{
		qaws_vec2 pts_2d[64];
		unsigned int ci;
		for (ci = 0; ci < 64; ci++)
		{
			qaws_eval_result_3d er;
			qaws_scalar t = (qaws_scalar)ci / (qaws_scalar)63;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(crv, t, QAWS_EVAL_FLAG_POSITION, &er);
			pts_2d[ci].x = er.position.x + er.position.y * proj_cx;
			pts_2d[ci].y = er.position.z + er.position.y * proj_cy;
		}
		svg_polyline(&svg, pts_2d, 64, "#ff8800", 2);
	}

	/* Draw intersection points */
	for (i = 0; i < count; i++)
	{
		qaws_scalar sx = hits[i].position.x + hits[i].position.y * proj_cx;
		qaws_scalar sy = hits[i].position.z + hits[i].position.y * proj_cy;
		svg_dot(&svg, sx, sy, "#ff2222", 6);
	}

	{
		char buf[64];
		sprintf(buf, "%u intersections", count);
		svg_label(&svg, (qaws_scalar)0.2, (qaws_scalar)3.8, buf, "#ffffff");
	}

	svg_close(&svg);
	printf("  -> " OBJ_OUTPUT_DIR "/33_surface_curve_intersection.svg (%u hits)\n", count);
	qaws_curve_destroy(crv);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: SVG closest point on surface (XZ projection)               */
/* ------------------------------------------------------------------ */
static void visual_svg_closest_point(void)
{
	svg_writer svg;
	qaws_surface* surf = make_dome();
	qaws_vec3 queries[3];
	unsigned int i, ui, vi;

	printf("visual_svg_closest_point\n");
	svg_ensure_output_dir();
	if (!surf) return;

	queries[0].x = 1; queries[0].y = 1; queries[0].z = 5;
	queries[1].x = (qaws_scalar)0.2; queries[1].y = 1; queries[1].z = 3;
	queries[2].x = (qaws_scalar)1.8; queries[2].y = 1; queries[2].z = 3;

	if (!svg_open(&svg, OBJ_OUTPUT_DIR "/33_closest_point.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-0.5,
		(qaws_scalar)3.5, (qaws_scalar)6.5,
		(qaws_scalar)500, (qaws_scalar)500))
	{
		qaws_surface_destroy(surf);
		return;
	}

	/* Draw surface wireframe (iso-v lines at y=1, XZ projection) */
	for (vi = 0; vi <= 6; vi++)
	{
		qaws_vec2 pts_2d[33];
		qaws_scalar v = (qaws_scalar)vi / (qaws_scalar)6;
		for (ui = 0; ui <= 32; ui++)
		{
			qaws_scalar u = (qaws_scalar)ui / (qaws_scalar)32;
			qaws_surface_eval_result r;
			memset(&r, 0, sizeof(r));
			qaws_surface_evaluate(surf, u, v, QAWS_SURFACE_EVAL_POSITION, &r);
			pts_2d[ui].x = r.position.x;
			pts_2d[ui].y = r.position.z;
		}
		svg_polyline(&svg, pts_2d, 33, "#5588bb", 1);
	}

	/* Query points, closest points, and connecting lines */
	for (i = 0; i < 3; i++)
	{
		qaws_scalar u, v;
		qaws_vec3 closest;
		if (qaws_surface_find_closest_point(surf, queries[i],
			&u, &v, &closest) == QAWS_STATUS_OK)
		{
			svg_dot(&svg, queries[i].x, queries[i].z, "#ff4444", 5);
			svg_dot(&svg, closest.x, closest.z, "#44ff44", 5);
			svg_line(&svg, queries[i].x, queries[i].z,
				closest.x, closest.z, "#ffcc44", 1);
		}
	}

	svg_label(&svg, (qaws_scalar)0.2, (qaws_scalar)5.8, "red=query green=closest", "#ffffff");
	svg_close(&svg);
	printf("  -> " OBJ_OUTPUT_DIR "/33_closest_point.svg\n");
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ dual tessellation                                      */
/*  Dual mesh: vertices at leaf-quad centers, faces around each        */
/*  shared primal corner.  Complements the primal (fan) tessellation.  */
/* ------------------------------------------------------------------ */

/* Leaf quad from quadtree subdivision */
typedef struct {
	qaws_scalar u0, v0, u1, v1;
} dual_leaf;

/* Shared corner with references to surrounding leaves */
typedef struct {
	qaws_scalar u, v;
	unsigned int leaf_ids[8];
	unsigned int count;
} dual_corner;

static void visual_obj_dual_tessellation(void)
{
	obj_writer w;
	qaws_surface *surf = make_dome();

	/* Quadtree stack */
	qaws_scalar su0, sv0, su1, sv1;
	unsigned int sdepth;

	dual_leaf *leaves;
	unsigned int leaf_count;
	dual_corner *corners;
	unsigned int corner_count;

	/* Stack for iterative subdivision */
	qaws_scalar *stk;   /* 5 scalars per entry: u0,v0,u1,v1,depth */
	unsigned int stk_top;

	qaws_range u_range, v_range;
	unsigned int first_v, first_n;
	unsigned int i, j, k;

	printf("visual_obj_dual_tessellation\n");
	svg_ensure_output_dir();
	if (!surf) return;

	leaves = (dual_leaf *)malloc(4096 * sizeof(dual_leaf));
	corners = (dual_corner *)malloc(8192 * sizeof(dual_corner));
	stk = (qaws_scalar *)malloc(4096 * 5 * sizeof(qaws_scalar));
	if (!leaves || !corners || !stk)
	{
		free(leaves); free(corners); free(stk);
		qaws_surface_destroy(surf);
		return;
	}

	/* ---- Phase 1: Quadtree subdivision ---- */
	u_range = qaws_surface_get_u_range(surf);
	v_range = qaws_surface_get_v_range(surf);
	leaf_count = 0;

	stk[0] = u_range.min_value;
	stk[1] = v_range.min_value;
	stk[2] = u_range.max_value;
	stk[3] = v_range.max_value;
	stk[4] = 0;
	stk_top = 1;

	while (stk_top > 0)
	{
		qaws_scalar u_mid, v_mid;
		qaws_surface_eval_result r00, r10, r01, r11, rc;
		qaws_vec3 bm;
		qaws_scalar dx, dy, dz, flatness;
		int subdivide = 0;
		unsigned int base;

		--stk_top;
		base = stk_top * 5;
		su0 = stk[base]; sv0 = stk[base + 1];
		su1 = stk[base + 2]; sv1 = stk[base + 3];
		sdepth = (unsigned int)stk[base + 4];

		u_mid = (qaws_scalar)0.5 * (su0 + su1);
		v_mid = (qaws_scalar)0.5 * (sv0 + sv1);

		memset(&r00, 0, sizeof(r00));
		memset(&r10, 0, sizeof(r10));
		memset(&r01, 0, sizeof(r01));
		memset(&r11, 0, sizeof(r11));
		memset(&rc, 0, sizeof(rc));
		qaws_surface_evaluate(surf, su0, sv0,
			QAWS_SURFACE_EVAL_POSITION, &r00);
		qaws_surface_evaluate(surf, su1, sv0,
			QAWS_SURFACE_EVAL_POSITION, &r10);
		qaws_surface_evaluate(surf, su0, sv1,
			QAWS_SURFACE_EVAL_POSITION, &r01);
		qaws_surface_evaluate(surf, su1, sv1,
			QAWS_SURFACE_EVAL_POSITION, &r11);
		qaws_surface_evaluate(surf, u_mid, v_mid,
			QAWS_SURFACE_EVAL_POSITION, &rc);

		if (sdepth < 5)
		{
			bm.x = (qaws_scalar)0.25
				* (r00.position.x + r10.position.x
					+ r01.position.x + r11.position.x);
			bm.y = (qaws_scalar)0.25
				* (r00.position.y + r10.position.y
					+ r01.position.y + r11.position.y);
			bm.z = (qaws_scalar)0.25
				* (r00.position.z + r10.position.z
					+ r01.position.z + r11.position.z);
			dx = rc.position.x - bm.x;
			dy = rc.position.y - bm.y;
			dz = rc.position.z - bm.z;
			flatness = (qaws_scalar)sqrt(
				(double)(dx * dx + dy * dy + dz * dz));
			if (flatness > (qaws_scalar)0.05)
				subdivide = 1;
		}

		if (subdivide && stk_top + 4 <= 4096)
		{
			unsigned int b;
			b = stk_top * 5;
			stk[b] = su0; stk[b+1] = sv0;
			stk[b+2] = u_mid; stk[b+3] = v_mid;
			stk[b+4] = (qaws_scalar)(sdepth + 1);
			++stk_top;

			b = stk_top * 5;
			stk[b] = u_mid; stk[b+1] = sv0;
			stk[b+2] = su1; stk[b+3] = v_mid;
			stk[b+4] = (qaws_scalar)(sdepth + 1);
			++stk_top;

			b = stk_top * 5;
			stk[b] = su0; stk[b+1] = v_mid;
			stk[b+2] = u_mid; stk[b+3] = sv1;
			stk[b+4] = (qaws_scalar)(sdepth + 1);
			++stk_top;

			b = stk_top * 5;
			stk[b] = u_mid; stk[b+1] = v_mid;
			stk[b+2] = su1; stk[b+3] = sv1;
			stk[b+4] = (qaws_scalar)(sdepth + 1);
			++stk_top;
		}
		else if (leaf_count < 4096)
		{
			leaves[leaf_count].u0 = su0;
			leaves[leaf_count].v0 = sv0;
			leaves[leaf_count].u1 = su1;
			leaves[leaf_count].v1 = sv1;
			++leaf_count;
		}
	}

	free(stk);

	/* ---- Phase 2: Build corner-to-leaf map ---- */
	/* Step 1: collect all unique corner positions */
	corner_count = 0;
	for (i = 0; i < leaf_count; i++)
	{
		qaws_scalar cvs[4][2];
		cvs[0][0] = leaves[i].u0; cvs[0][1] = leaves[i].v0;
		cvs[1][0] = leaves[i].u1; cvs[1][1] = leaves[i].v0;
		cvs[2][0] = leaves[i].u0; cvs[2][1] = leaves[i].v1;
		cvs[3][0] = leaves[i].u1; cvs[3][1] = leaves[i].v1;

		for (j = 0; j < 4; j++)
		{
			unsigned int ci;
			int found = 0;
			for (ci = 0; ci < corner_count; ci++)
			{
				if (fabs(corners[ci].u - cvs[j][0]) < 1e-10
					&& fabs(corners[ci].v - cvs[j][1]) < 1e-10)
				{
					found = 1;
					break;
				}
			}
			if (!found && corner_count < 8192)
			{
				corners[corner_count].u = cvs[j][0];
				corners[corner_count].v = cvs[j][1];
				corners[corner_count].count = 0;
				++corner_count;
			}
		}
	}

	/* Step 2: for each corner, register ALL leaves whose boundary
	   (any edge, not just corners) passes through it.
	   This correctly handles T-junctions where a finer leaf's corner
	   lies on a coarser leaf's edge interior. */
	{
		double eps = 1e-10;
		unsigned int ci;
		for (ci = 0; ci < corner_count; ci++)
		{
			double cu = (double)corners[ci].u;
			double cv = (double)corners[ci].v;

			for (i = 0; i < leaf_count; i++)
			{
				double lu0 = (double)leaves[i].u0;
				double lv0 = (double)leaves[i].v0;
				double lu1 = (double)leaves[i].u1;
				double lv1 = (double)leaves[i].v1;
				int on_boundary = 0;

				/* On a horizontal edge (top or bottom)? */
				if ((fabs(cv - lv0) < eps || fabs(cv - lv1) < eps)
					&& cu >= lu0 - eps && cu <= lu1 + eps)
					on_boundary = 1;
				/* On a vertical edge (left or right)? */
				if ((fabs(cu - lu0) < eps || fabs(cu - lu1) < eps)
					&& cv >= lv0 - eps && cv <= lv1 + eps)
					on_boundary = 1;

				if (on_boundary && corners[ci].count < 8)
					corners[ci].leaf_ids[corners[ci].count++] = i;
			}
		}
	}

	/* ---- Phase 3: Sort leaves around each corner by angle ---- */
	for (i = 0; i < corner_count; i++)
	{
		if (corners[i].count < 2) continue;
		/* Bubble sort by atan2(center_v - corner_v, center_u - corner_u) */
		for (j = 0; j < corners[i].count; j++)
		{
			for (k = j + 1; k < corners[i].count; k++)
			{
				unsigned int lj = corners[i].leaf_ids[j];
				unsigned int lk = corners[i].leaf_ids[k];
				double aj = atan2(
					(double)((leaves[lj].v0 + leaves[lj].v1) * 0.5
						- corners[i].v),
					(double)((leaves[lj].u0 + leaves[lj].u1) * 0.5
						- corners[i].u));
				double ak = atan2(
					(double)((leaves[lk].v0 + leaves[lk].v1) * 0.5
						- corners[i].v),
					(double)((leaves[lk].u0 + leaves[lk].u1) * 0.5
						- corners[i].u));
				if (ak < aj)
				{
					unsigned int tmp = corners[i].leaf_ids[j];
					corners[i].leaf_ids[j] = corners[i].leaf_ids[k];
					corners[i].leaf_ids[k] = tmp;
				}
			}
		}
	}

	/* ---- Phase 4: Write dual mesh OBJ ---- */
	if (!obj_open(&w, OBJ_OUTPUT_DIR "/33_dual_tessellation.obj",
		OBJ_OUTPUT_DIR "/33_dual_tessellation.mtl"))
	{
		free(leaves); free(corners);
		qaws_surface_destroy(surf);
		return;
	}

	obj_material(&w, "dual", 0.9, 0.5, 0.2);
	obj_group(&w, "dual_mesh");
	obj_use_material(&w, "dual");

	/* Emit dual vertices: one per leaf at quad center */
	first_v = w.vertex_count + 1;
	first_n = w.normal_count + 1;
	for (i = 0; i < leaf_count; i++)
	{
		qaws_scalar u_mid = (qaws_scalar)0.5
			* (leaves[i].u0 + leaves[i].u1);
		qaws_scalar v_mid = (qaws_scalar)0.5
			* (leaves[i].v0 + leaves[i].v1);
		qaws_surface_eval_result r;
		memset(&r, 0, sizeof(r));
		qaws_surface_evaluate(surf, u_mid, v_mid,
			QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
		obj_vertex(&w, r.position);
		obj_normal(&w, r.normal);
	}

	/* Emit dual faces: one per interior corner (>= 3 leaves) */
	for (i = 0; i < corner_count; i++)
	{
		if (corners[i].count < 3)
			continue;
		fprintf(w.fp, "f");
		for (j = 0; j < corners[i].count; j++)
		{
			unsigned int vi = first_v + corners[i].leaf_ids[j];
			unsigned int ni = first_n + corners[i].leaf_ids[j];
			fprintf(w.fp, " %u//%u", vi, ni);
		}
		fprintf(w.fp, "\n");
	}

	obj_close(&w);
	{
		unsigned int face_count = 0;
		for (i = 0; i < corner_count; i++)
			if (corners[i].count >= 3) face_count++;
		printf("  -> " OBJ_OUTPUT_DIR "/33_dual_tessellation.obj (%u leaves, %u dual faces)\n",
			leaf_count, face_count);
	}

	free(leaves);
	free(corners);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_33_advanced_inspection_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_closest_point_on_surface();
	test_curve_plane_intersection();
	test_surface_curve_intersection();
	test_adaptive_tessellation();

	/* Visual output */
	visual_obj_tessellation();
	visual_obj_dual_tessellation();
	visual_svg_plane_intersection();
	visual_obj_closest_point();
	visual_obj_surface_curve_intersection();
	visual_svg_surface_curve_intersection();
	visual_svg_closest_point();

	printf("33_advanced_inspection: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
