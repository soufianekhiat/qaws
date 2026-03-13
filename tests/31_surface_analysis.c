#include "test_common.h"

/* ------------------------------------------------------------------ */
/*  Helper: create a bilinear Bezier surface (flat quad)               */
/* ------------------------------------------------------------------ */
static qaws_surface* make_bilinear_surface(void)
{
	qaws_surface* surf = NULL;
	qaws_surface_bezier_desc desc;
	qaws_vec3 cp[4];

	/* Flat quad: (0,0,0) (2,0,0) (0,3,0) (2,3,0) */
	cp[0].x = 0; cp[0].y = 0; cp[0].z = 0;
	cp[1].x = 0; cp[1].y = 3; cp[1].z = 0;
	cp[2].x = 2; cp[2].y = 0; cp[2].z = 0;
	cp[3].x = 2; cp[3].y = 3; cp[3].z = 0;

	desc.u_degree = 1;
	desc.v_degree = 1;
	desc.control_points = cp;
	desc.u_point_count = 2;
	desc.v_point_count = 2;

	qaws_surface_create_bezier(&desc, &surf);
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Helper: create a curved Bezier surface (dome)                      */
/* ------------------------------------------------------------------ */
static qaws_surface* make_dome_surface(void)
{
	qaws_surface* surf = NULL;
	qaws_surface_bezier_desc desc;
	qaws_vec3 cp[9];

	/* 3x3 biquadratic patch forming a dome */
	cp[0].x = 0; cp[0].y = 0; cp[0].z = 0;
	cp[1].x = 0; cp[1].y = 1; cp[1].z = 0;
	cp[2].x = 0; cp[2].y = 2; cp[2].z = 0;
	cp[3].x = 1; cp[3].y = 0; cp[3].z = 0;
	cp[4].x = 1; cp[4].y = 1; cp[4].z = 2;  /* raised center */
	cp[5].x = 1; cp[5].y = 2; cp[5].z = 0;
	cp[6].x = 2; cp[6].y = 0; cp[6].z = 0;
	cp[7].x = 2; cp[7].y = 1; cp[7].z = 0;
	cp[8].x = 2; cp[8].y = 2; cp[8].z = 0;

	desc.u_degree = 2;
	desc.v_degree = 2;
	desc.control_points = cp;
	desc.u_point_count = 3;
	desc.v_point_count = 3;

	qaws_surface_create_bezier(&desc, &surf);
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Helper: create a wavy revolution surface (sinusoidal profile)      */
/*  Alternating bulge/waist → positive and negative curvature bands.   */
/* ------------------------------------------------------------------ */
static qaws_surface* make_wavy_surface(qaws_curve** out_profile)
{
	qaws_surface* surf = NULL;
	qaws_curve* prof = NULL;
	qaws_scalar pts[10];
	qaws_bezier_desc bd;
	qaws_surface_revolution_desc rd;

	/* Degree-4 Bezier profile approximating a sine wave.
	   x = radius (2.0 ± 0.8), y = height 0..3.
	   Bulge at y~0.75, waist at y~2.25. */
	pts[0] = 2; pts[1] = 0;
	pts[2] = 3; pts[3] = (qaws_scalar)0.75;
	pts[4] = 2; pts[5] = (qaws_scalar)1.5;
	pts[6] = 1; pts[7] = (qaws_scalar)2.25;
	pts[8] = 2; pts[9] = 3;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D;
	bd.degree = 4;
	bd.control_points = pts;
	bd.control_point_count = 5;
	qaws_curve_create_bezier(&bd, &prof);
	if (!prof) return NULL;

	memset(&rd, 0, sizeof(rd));
	rd.profile = prof;
	rd.axis_origin.x = 0; rd.axis_origin.y = 0; rd.axis_origin.z = 0;
	rd.axis_direction.x = 0; rd.axis_direction.y = 0; rd.axis_direction.z = 1;
	rd.angle = 0;
	qaws_surface_create_revolution(&rd, &surf);
	*out_profile = prof;
	return surf;
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ wavy surface with per-vertex curvature color           */
/* ------------------------------------------------------------------ */
static void visual_obj_dome(void)
{
	obj_writer w;
	qaws_curve* prof = NULL;
	qaws_surface* surf = make_wavy_surface(&prof);

	printf("visual_obj_dome\n");
	svg_ensure_output_dir();
	if (!surf) { qaws_curve_destroy(prof); return; }

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/dome_curvature.obj",
		OBJ_OUTPUT_DIR "/dome_curvature.mtl"))
	{
		qaws_surface_destroy(surf);
		qaws_curve_destroy(prof);
		return;
	}

	obj_material(&w, "curvature", 1.0, 1.0, 1.0);
	obj_comment(&w, "Per-vertex color = Gaussian curvature");
	obj_comment(&w, "Blue = negative (saddle), White = zero, Red = positive (dome)");
	obj_group(&w, "wavy_curvature");
	obj_use_material(&w, "curvature");
	obj_surface_mesh_curvature(&w, surf, 48, 48);

	obj_close(&w);
	qaws_surface_destroy(surf);
	qaws_curve_destroy(prof);
	printf("  -> " OBJ_OUTPUT_DIR "/dome_curvature.obj\n");
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ export API roundtrip                                   */
/* ------------------------------------------------------------------ */
static void visual_obj_export_api(void)
{
	qaws_surface* surf = make_dome_surface();
	char* buf;
	unsigned int length = 0;
	FILE* fp;

	printf("visual_obj_export_api\n");
	if (!surf) return;

	buf = (char*)malloc(262144);
	if (!buf) { qaws_surface_destroy(surf); return; }

	if (qaws_surface_export_obj(surf, 32, 32, 1, buf, 262144, &length)
		== QAWS_STATUS_OK)
	{
		fp = fopen(OBJ_OUTPUT_DIR "/dome_curvature_api.obj", "w");
		if (fp)
		{
			fwrite(buf, 1, length, fp);
			fclose(fp);
			printf("  -> " OBJ_OUTPUT_DIR "/dome_curvature_api.obj\n");
		}
	}

	free(buf);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: SVG Hausdorff distance curves                              */
/* ------------------------------------------------------------------ */
static void visual_svg_hausdorff(void)
{
	svg_writer svg;
	qaws_curve* ca = NULL;
	qaws_curve* cb = NULL;
	qaws_scalar dist = 0;
	unsigned int i;
	unsigned int n = 64;
	char label_buf[64];
	qaws_vec2 samples_a[128], samples_b[128];
	unsigned int na, nb;

	printf("visual_svg_hausdorff\n");
	svg_ensure_output_dir();

	/* Quadratic curve vs offset-ish curve */
	{
		qaws_vec2 pts_a[] = { {0, 2}, {3, 5}, {6, 2} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 2;
		d.control_points = pts_a;
		d.control_point_count = 3;
		qaws_curve_create_bezier(&d, &ca);
	}
	{
		qaws_vec2 pts_b[] = { {0, (qaws_scalar)3.5}, {3, (qaws_scalar)6.5}, {6, (qaws_scalar)3.5} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 2;
		d.control_points = pts_b;
		d.control_point_count = 3;
		qaws_curve_create_bezier(&d, &cb);
	}

	qaws_curve_compute_hausdorff_distance_2d(ca, cb, n, &dist);

	/* view: x [-1,7], y [-1,8], SVG 700x500 */
	if (!svg_open(&svg, OBJ_OUTPUT_DIR "/hausdorff_distance.svg",
		(qaws_scalar)-1, (qaws_scalar)-1,
		(qaws_scalar)8, (qaws_scalar)9,
		(qaws_scalar)700, (qaws_scalar)500))
	{
		qaws_curve_destroy(ca);
		qaws_curve_destroy(cb);
		return;
	}

	na = svg_sample_curve(ca, samples_a, 128);
	nb = svg_sample_curve(cb, samples_b, 128);
	svg_polyline(&svg, samples_a, na, "#4488ff", 2);
	svg_polyline(&svg, samples_b, nb, "#ff4444", 2);

	/* Draw distance lines between closest point pairs */
	for (i = 0; i < n; i += 4)
	{
		qaws_scalar t = (qaws_scalar)i / (qaws_scalar)(n - 1);
		qaws_eval_result_2d ra, rb;
		qaws_scalar closest_t;

		qaws_curve_evaluate_2d(ca, t, QAWS_EVAL_FLAG_POSITION, &ra);
		qaws_curve_find_closest_parameter_2d(cb, ra.position, &closest_t);
		qaws_curve_evaluate_2d(cb, closest_t, QAWS_EVAL_FLAG_POSITION, &rb);

		svg_line(&svg, ra.position.x, ra.position.y,
			rb.position.x, rb.position.y, "#888888", 1);
	}

	sprintf(label_buf, "Hausdorff = %.3f", (double)dist);
	svg_label(&svg, (qaws_scalar)0.5, (qaws_scalar)7.5, label_buf, "#ffffff");

	svg_close(&svg);
	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
	printf("  -> " OBJ_OUTPUT_DIR "/hausdorff_distance.svg (dist=%.4f)\n", (double)dist);
}

/* ------------------------------------------------------------------ */
/*  Visual: OBJ bounding box visualization                             */
/* ------------------------------------------------------------------ */
static void visual_obj_bounds(void)
{
	obj_writer w;
	qaws_surface* surf = make_dome_surface();
	qaws_vec3 bmin, bmax;

	printf("visual_obj_bounds\n");
	if (!surf) return;

	qaws_surface_compute_bounds(surf, &bmin, &bmax);

	if (!obj_open(&w,
		OBJ_OUTPUT_DIR "/surface_bounds.obj",
		OBJ_OUTPUT_DIR "/surface_bounds.mtl"))
	{
		qaws_surface_destroy(surf);
		return;
	}

	obj_material(&w, "surf", 0.3, 0.7, 0.3);
	obj_material(&w, "bbox", 1.0, 0.2, 0.2);

	obj_group(&w, "surface");
	obj_use_material(&w, "surf");
	obj_surface_mesh(&w, surf, 24, 24);

	/* Draw bounding box as 12 wireframe edges + 6 translucent quad faces */
	obj_group(&w, "bounding_box");
	obj_use_material(&w, "bbox");
	{
		qaws_vec3 c[8];
		unsigned int first_bv;
		unsigned int ei;
		unsigned int edges[12][2] = {
			{0,1},{1,3},{3,2},{2,0},
			{4,5},{5,7},{7,6},{6,4},
			{0,4},{1,5},{2,6},{3,7}
		};
		c[0].x = bmin.x; c[0].y = bmin.y; c[0].z = bmin.z;
		c[1].x = bmax.x; c[1].y = bmin.y; c[1].z = bmin.z;
		c[2].x = bmin.x; c[2].y = bmax.y; c[2].z = bmin.z;
		c[3].x = bmax.x; c[3].y = bmax.y; c[3].z = bmin.z;
		c[4].x = bmin.x; c[4].y = bmin.y; c[4].z = bmax.z;
		c[5].x = bmax.x; c[5].y = bmin.y; c[5].z = bmax.z;
		c[6].x = bmin.x; c[6].y = bmax.y; c[6].z = bmax.z;
		c[7].x = bmax.x; c[7].y = bmax.y; c[7].z = bmax.z;

		for (ei = 0; ei < 12; ei++)
			obj_line_segment(&w, c[edges[ei][0]], c[edges[ei][1]]);

		/* Box faces as quads (visible in most viewers) */
		first_bv = w.vertex_count + 1;
		for (ei = 0; ei < 8; ei++)
			obj_vertex(&w, c[ei]);
		/* 6 faces: -Z, +Z, -Y, +Y, -X, +X (outward winding) */
		fprintf(w.fp, "f %u %u %u %u\n", first_bv+0, first_bv+2, first_bv+3, first_bv+1);
		fprintf(w.fp, "f %u %u %u %u\n", first_bv+4, first_bv+5, first_bv+7, first_bv+6);
		fprintf(w.fp, "f %u %u %u %u\n", first_bv+0, first_bv+1, first_bv+5, first_bv+4);
		fprintf(w.fp, "f %u %u %u %u\n", first_bv+2, first_bv+6, first_bv+7, first_bv+3);
		fprintf(w.fp, "f %u %u %u %u\n", first_bv+0, first_bv+4, first_bv+6, first_bv+2);
		fprintf(w.fp, "f %u %u %u %u\n", first_bv+1, first_bv+3, first_bv+7, first_bv+5);
	}

	obj_close(&w);
	qaws_surface_destroy(surf);
	printf("  -> " OBJ_OUTPUT_DIR "/surface_bounds.obj\n");
}

/* ------------------------------------------------------------------ */
/*  Test: OBJ surface export                                           */
/* ------------------------------------------------------------------ */
static void test_obj_export(void)
{
	qaws_surface* surf = NULL;
	char buf[8192];
	unsigned int length = 0;
	qaws_status s;

	printf("test_obj_export\n");

	surf = make_bilinear_surface();
	TEST_ASSERT(surf != NULL, "obj export surface created");

	/* Export with normals */
	s = qaws_surface_export_obj(surf, 4, 4, 1, buf, sizeof(buf), &length);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(length > 0, "obj export produced output");

	/* Check that output contains vertex and face data */
	TEST_ASSERT(strstr(buf, "v ") != NULL, "obj has vertices");
	TEST_ASSERT(strstr(buf, "vn ") != NULL, "obj has normals");
	TEST_ASSERT(strstr(buf, "f ") != NULL, "obj has faces");

	/* Export without normals */
	length = 0;
	s = qaws_surface_export_obj(surf, 4, 4, 0, buf, sizeof(buf), &length);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(length > 0, "obj export no-normals produced output");
	TEST_ASSERT(strstr(buf, "vn ") == NULL, "obj no-normals has no vn");
	TEST_ASSERT(strstr(buf, "f ") != NULL, "obj no-normals has faces");

	/* Minimum samples */
	length = 0;
	s = qaws_surface_export_obj(surf, 2, 2, 0, buf, sizeof(buf), &length);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(strstr(buf, "f ") != NULL, "obj 2x2 has faces");

	/* Error: invalid samples */
	s = qaws_surface_export_obj(surf, 1, 4, 0, buf, sizeof(buf), &length);
	TEST_ASSERT(s != QAWS_STATUS_OK, "obj u_samples=1 rejected");
	s = qaws_surface_export_obj(surf, 4, 1, 0, buf, sizeof(buf), &length);
	TEST_ASSERT(s != QAWS_STATUS_OK, "obj v_samples=1 rejected");

	/* Error: null args */
	s = qaws_surface_export_obj(NULL, 4, 4, 0, buf, sizeof(buf), &length);
	TEST_ASSERT(s != QAWS_STATUS_OK, "obj null surface rejected");
	s = qaws_surface_export_obj(surf, 4, 4, 0, NULL, sizeof(buf), &length);
	TEST_ASSERT(s != QAWS_STATUS_OK, "obj null buffer rejected");

	/* Buffer too small */
	{
		char tiny[16];
		unsigned int tiny_len = 0;
		s = qaws_surface_export_obj(surf, 4, 4, 1, tiny, sizeof(tiny), &tiny_len);
		/* Should either truncate or return buffer too small */
		TEST_ASSERT(tiny_len > 0 || s == QAWS_STATUS_BUFFER_TOO_SMALL,
			"obj tiny buffer handled");
	}

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: Surface bounding box                                         */
/* ------------------------------------------------------------------ */
static void test_surface_bounds(void)
{
	qaws_surface* surf = NULL;
	qaws_vec3 bmin, bmax;
	qaws_status s;

	printf("test_surface_bounds\n");

	/* Flat surface: bounds should be exactly the extents */
	surf = make_bilinear_surface();
	TEST_ASSERT(surf != NULL, "bounds surface created");

	s = qaws_surface_compute_bounds(surf, &bmin, &bmax);
	TEST_ASSERT_STATUS(s);

	/* x: [0, 2], y: [0, 3], z: 0 */
	TEST_ASSERT(bmin.x >= (qaws_scalar)-0.01 && bmin.x <= (qaws_scalar)0.01,
		"flat bounds min x ~0");
	TEST_ASSERT(bmax.x >= (qaws_scalar)1.99 && bmax.x <= (qaws_scalar)2.01,
		"flat bounds max x ~2");
	TEST_ASSERT(bmin.y >= (qaws_scalar)-0.01 && bmin.y <= (qaws_scalar)0.01,
		"flat bounds min y ~0");
	TEST_ASSERT(bmax.y >= (qaws_scalar)2.99 && bmax.y <= (qaws_scalar)3.01,
		"flat bounds max y ~3");
	TEST_ASSERT(fabs(bmin.z) < (qaws_scalar)0.01, "flat bounds min z ~0");
	TEST_ASSERT(fabs(bmax.z) < (qaws_scalar)0.01, "flat bounds max z ~0");

	qaws_surface_destroy(surf);

	/* Dome surface: z should go above 0 */
	surf = make_dome_surface();
	TEST_ASSERT(surf != NULL, "dome surface created");

	s = qaws_surface_compute_bounds(surf, &bmin, &bmax);
	TEST_ASSERT_STATUS(s);

	TEST_ASSERT(bmin.z >= (qaws_scalar)-0.01, "dome min z >= 0");
	TEST_ASSERT(bmax.z > (qaws_scalar)0.1, "dome max z > 0 (raised)");

	/* Null args */
	s = qaws_surface_compute_bounds(NULL, &bmin, &bmax);
	TEST_ASSERT(s != QAWS_STATUS_OK, "bounds null surface rejected");
	s = qaws_surface_compute_bounds(surf, NULL, &bmax);
	TEST_ASSERT(s != QAWS_STATUS_OK, "bounds null min rejected");

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: Surface area                                                 */
/* ------------------------------------------------------------------ */
static void test_surface_area(void)
{
	qaws_surface* surf = NULL;
	qaws_scalar area = 0;
	qaws_status s;

	printf("test_surface_area\n");

	/* Flat 2x3 quad: area should be 6.0 */
	surf = make_bilinear_surface();
	TEST_ASSERT(surf != NULL, "area surface created");

	s = qaws_surface_compute_area(surf, &area);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(area, (qaws_scalar)6.0), "flat area ~6.0");

	qaws_surface_destroy(surf);

	/* Dome surface: area should be > 4 (2x2 base, raised center adds area) */
	surf = make_dome_surface();
	TEST_ASSERT(surf != NULL, "dome area surface created");

	s = qaws_surface_compute_area(surf, &area);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(area > (qaws_scalar)4.0, "dome area > 4.0");
	TEST_ASSERT(area < (qaws_scalar)20.0, "dome area reasonable upper bound");

	/* Null args */
	s = qaws_surface_compute_area(NULL, &area);
	TEST_ASSERT(s != QAWS_STATUS_OK, "area null surface rejected");
	s = qaws_surface_compute_area(surf, NULL);
	TEST_ASSERT(s != QAWS_STATUS_OK, "area null out rejected");

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: Surface curvature                                            */
/* ------------------------------------------------------------------ */
static void test_surface_curvature(void)
{
	qaws_surface* surf = NULL;
	qaws_scalar K = 0, H = 0;
	qaws_surface_curvature_result cr;
	qaws_status s;

	printf("test_surface_curvature\n");

	/* Flat surface: both curvatures should be 0 */
	surf = make_bilinear_surface();
	TEST_ASSERT(surf != NULL, "curvature flat surface created");

	s = qaws_surface_compute_gaussian_curvature(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &K);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(K) < TOLERANCE_LOOSE, "flat Gaussian curvature ~0");

	s = qaws_surface_compute_mean_curvature(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &H);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(H) < TOLERANCE_LOOSE, "flat mean curvature ~0");

	s = qaws_surface_compute_principal_curvatures(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &cr);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(fabs(cr.kappa1) < TOLERANCE_LOOSE, "flat kappa1 ~0");
	TEST_ASSERT(fabs(cr.kappa2) < TOLERANCE_LOOSE, "flat kappa2 ~0");
	TEST_ASSERT(fabs(cr.gaussian) < TOLERANCE_LOOSE, "flat cr.gaussian ~0");
	TEST_ASSERT(fabs(cr.mean) < TOLERANCE_LOOSE, "flat cr.mean ~0");

	qaws_surface_destroy(surf);

	/* Dome surface at center: should have nonzero curvature */
	surf = make_dome_surface();
	TEST_ASSERT(surf != NULL, "curvature dome surface created");

	s = qaws_surface_compute_gaussian_curvature(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &K);
	TEST_ASSERT_STATUS(s);
	/* Dome has positive Gaussian curvature at center */
	TEST_ASSERT(K > (qaws_scalar)0.0, "dome Gaussian curvature > 0");

	s = qaws_surface_compute_mean_curvature(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &H);
	TEST_ASSERT_STATUS(s);
	/* H can be positive or negative depending on normal direction */
	TEST_ASSERT(fabs(H) > TOLERANCE, "dome mean curvature nonzero");

	s = qaws_surface_compute_principal_curvatures(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &cr);
	TEST_ASSERT_STATUS(s);
	/* K = kappa1 * kappa2 */
	TEST_ASSERT(approx_eq_loose(cr.gaussian, cr.kappa1 * cr.kappa2),
		"K = kappa1 * kappa2");
	/* H = (kappa1 + kappa2) / 2 */
	TEST_ASSERT(approx_eq_loose(cr.mean,
		(cr.kappa1 + cr.kappa2) / (qaws_scalar)2.0),
		"H = (kappa1 + kappa2) / 2");

	/* Null args */
	s = qaws_surface_compute_gaussian_curvature(NULL,
		(qaws_scalar)0.5, (qaws_scalar)0.5, &K);
	TEST_ASSERT(s != QAWS_STATUS_OK, "gaussian null surface rejected");
	s = qaws_surface_compute_mean_curvature(surf,
		(qaws_scalar)0.5, (qaws_scalar)0.5, NULL);
	TEST_ASSERT(s != QAWS_STATUS_OK, "mean null out rejected");

	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: Hausdorff distance                                           */
/* ------------------------------------------------------------------ */
static void test_hausdorff_2d(void)
{
	qaws_curve* ca = NULL;
	qaws_curve* cb = NULL;
	qaws_scalar dist = 0;
	qaws_status s;

	printf("test_hausdorff_2d\n");

	/* Two parallel lines: y=0 and y=1, both from x=0 to x=10 */
	{
		qaws_vec2 pts_a[] = { {0, 0}, {10, 0} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 1;
		d.control_points = pts_a;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &ca);
	}
	{
		qaws_vec2 pts_b[] = { {0, 1}, {10, 1} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 1;
		d.control_points = pts_b;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &cb);
	}

	s = qaws_curve_compute_hausdorff_distance_2d(ca, cb, 64, &dist);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(dist, (qaws_scalar)1.0), "parallel lines hausdorff ~1.0");

	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);

	/* Identical curves: distance should be ~0 */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {3, 1} };
		qaws_bezier_desc d;
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = 2;
		d.control_points = pts;
		d.control_point_count = 3;
		qaws_curve_create_bezier(&d, &ca);
		qaws_curve_create_bezier(&d, &cb);
	}

	s = qaws_curve_compute_hausdorff_distance_2d(ca, cb, 64, &dist);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(dist < TOLERANCE_LOOSE, "identical curves hausdorff ~0");

	/* Null args */
	s = qaws_curve_compute_hausdorff_distance_2d(NULL, cb, 64, &dist);
	TEST_ASSERT(s != QAWS_STATUS_OK, "hausdorff null curve_a rejected");
	s = qaws_curve_compute_hausdorff_distance_2d(ca, cb, 64, NULL);
	TEST_ASSERT(s != QAWS_STATUS_OK, "hausdorff null out rejected");

	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
}

static void test_hausdorff_3d(void)
{
	qaws_curve* ca = NULL;
	qaws_curve* cb = NULL;
	qaws_scalar dist = 0;
	qaws_status s;

	printf("test_hausdorff_3d\n");

	/* Two parallel 3D lines, 2 units apart in z */
	{
		qaws_scalar pts_a[] = {0,0,0, 5,0,0};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 1;
		d.control_points = pts_a;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &ca);
	}
	{
		qaws_scalar pts_b[] = {0,0,2, 5,0,2};
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_3D;
		d.degree = 1;
		d.control_points = pts_b;
		d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &cb);
	}

	s = qaws_curve_compute_hausdorff_distance_3d(ca, cb, 64, &dist);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(approx_eq_loose(dist, (qaws_scalar)2.0), "3D parallel hausdorff ~2.0");

	qaws_curve_destroy(ca);
	qaws_curve_destroy(cb);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_31_surface_analysis_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_obj_export();
	test_surface_bounds();
	test_surface_area();
	test_surface_curvature();
	test_hausdorff_2d();
	test_hausdorff_3d();

	/* Visual output */
	visual_obj_dome();
	visual_obj_export_api();
	visual_svg_hausdorff();
	visual_obj_bounds();

	printf("31_surface_analysis: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
