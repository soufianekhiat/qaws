#include "test_common.h"
#include "qaws_platform.h"

/* ------------------------------------------------------------------ */
/*  Test: curve projection onto surface                                */
/* ------------------------------------------------------------------ */
static void test_curve_projection(void)
{
	qaws_surface* surf = NULL;
	qaws_curve* curve = NULL;
	qaws_curve* uv_curve = NULL;
	qaws_status s;

	printf("test_curve_projection\n");

	/* Create a flat bilinear surface in XY plane at z=0 */
	{
		qaws_surface_bilinear_desc bd;
		bd.p00.x = 0; bd.p00.y = 0; bd.p00.z = 0;
		bd.p10.x = 2; bd.p10.y = 0; bd.p10.z = 0;
		bd.p01.x = 0; bd.p01.y = 2; bd.p01.z = 0;
		bd.p11.x = 2; bd.p11.y = 2; bd.p11.z = 0;
		s = qaws_surface_create_bilinear(&bd, &surf);
		TEST_ASSERT_STATUS(s);
	}

	/* Create a 3D line curve above the surface: (0.5,0.5,1) to (1.5,1.5,1) */
	{
		qaws_scalar pts[6];
		qaws_bezier_desc bd;
		pts[0] = (qaws_scalar)0.5; pts[1] = (qaws_scalar)0.5; pts[2] = (qaws_scalar)1.0;
		pts[3] = (qaws_scalar)1.5; pts[4] = (qaws_scalar)1.5; pts[5] = (qaws_scalar)1.0;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_3D;
		bd.degree = 1;
		bd.control_points = pts;
		bd.control_point_count = 2;
		s = qaws_curve_create_bezier(&bd, &curve);
		TEST_ASSERT_STATUS(s);
	}

	/* Project the curve onto the surface */
	s = qaws_surface_project_curve(surf, curve, 64, &uv_curve);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(uv_curve != NULL, "projection produced a curve");
	TEST_ASSERT(qaws_curve_get_dimension(uv_curve) == QAWS_DIMENSION_2D,
		"projected curve is 2D");

	/* The projected (u,v) should go roughly from (0.25, 0.25) to (0.75, 0.75)
	   since the surface spans [0,2]x[0,2] and the curve goes from (0.5,0.5) to (1.5,1.5) */
	{
		qaws_range range = qaws_curve_get_parameter_range(uv_curve);
		qaws_eval_result_2d er_start, er_end;

		memset(&er_start, 0, sizeof(er_start));
		memset(&er_end, 0, sizeof(er_end));
		qaws_curve_evaluate_2d(uv_curve, range.min_value, QAWS_EVAL_FLAG_POSITION, &er_start);
		qaws_curve_evaluate_2d(uv_curve, range.max_value, QAWS_EVAL_FLAG_POSITION, &er_end);

		TEST_ASSERT(QAWS_FABS(er_start.position.x - (qaws_scalar)0.25) < (qaws_scalar)0.1,
			"projection start u ~ 0.25");
		TEST_ASSERT(QAWS_FABS(er_start.position.y - (qaws_scalar)0.25) < (qaws_scalar)0.1,
			"projection start v ~ 0.25");
		TEST_ASSERT(QAWS_FABS(er_end.position.x - (qaws_scalar)0.75) < (qaws_scalar)0.1,
			"projection end u ~ 0.75");
		TEST_ASSERT(QAWS_FABS(er_end.position.y - (qaws_scalar)0.75) < (qaws_scalar)0.1,
			"projection end v ~ 0.75");
	}

	qaws_curve_destroy(uv_curve);
	qaws_curve_destroy(curve);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Test: geodesic on a flat surface (should be a straight line)       */
/* ------------------------------------------------------------------ */
static void test_geodesic_flat(void)
{
	qaws_surface* surf = NULL;
	qaws_curve* geo = NULL;
	qaws_status s;
	unsigned int i;

	printf("test_geodesic_flat\n");

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

	/* Geodesic from (0.2, 0.2) to (0.8, 0.8) on a flat surface should be a straight line */
	s = qaws_surface_compute_geodesic(surf,
		(qaws_scalar)0.2, (qaws_scalar)0.2,
		(qaws_scalar)0.8, (qaws_scalar)0.8,
		20, 100, &geo);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(geo != NULL, "geodesic produced a curve");

	/* Check that all points lie approximately on z=0 (flat surface) */
	{
		qaws_range range = qaws_curve_get_parameter_range(geo);
		for (i = 0; i <= 10; i++)
		{
			qaws_eval_result_3d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)10;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(geo, t, QAWS_EVAL_FLAG_POSITION, &er);
			TEST_ASSERT(QAWS_FABS(er.position.z) < (qaws_scalar)0.1,
				"geodesic point z ~ 0");
		}
	}

	qaws_curve_destroy(geo);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: geodesic on a dome                                         */
/* ------------------------------------------------------------------ */
static void visual_obj_geodesic(void)
{
	obj_writer w;
	qaws_surface* surf = NULL;
	qaws_curve* geo = NULL;
	qaws_status s;
	unsigned int i;

	printf("visual_obj_geodesic\n");
	svg_ensure_output_dir();

	/* Biquadratic dome */
	{
		qaws_surface_biquadratic_desc bd;
		bd.control_points[0].x = -1; bd.control_points[0].y = -1; bd.control_points[0].z = 0;
		bd.control_points[1].x =  0; bd.control_points[1].y = -1; bd.control_points[1].z = 0;
		bd.control_points[2].x =  1; bd.control_points[2].y = -1; bd.control_points[2].z = 0;
		bd.control_points[3].x = -1; bd.control_points[3].y =  0; bd.control_points[3].z = 0;
		bd.control_points[4].x =  0; bd.control_points[4].y =  0; bd.control_points[4].z = (qaws_scalar)2.0;
		bd.control_points[5].x =  1; bd.control_points[5].y =  0; bd.control_points[5].z = 0;
		bd.control_points[6].x = -1; bd.control_points[6].y =  1; bd.control_points[6].z = 0;
		bd.control_points[7].x =  0; bd.control_points[7].y =  1; bd.control_points[7].z = 0;
		bd.control_points[8].x =  1; bd.control_points[8].y =  1; bd.control_points[8].z = 0;
		s = qaws_surface_create_biquadratic(&bd, &surf);
		if (s != QAWS_STATUS_OK) return;
	}

	s = qaws_surface_compute_geodesic(surf,
		(qaws_scalar)0.2, (qaws_scalar)0.2,
		(qaws_scalar)0.8, (qaws_scalar)0.8,
		20, 200, &geo);
	if (s != QAWS_STATUS_OK) { qaws_surface_destroy(surf); return; }

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/geodesic.obj",
		OBJ_OUTPUT_DIR "/geodesic.mtl"))
	{
		qaws_curve_destroy(geo); qaws_surface_destroy(surf); return;
	}

	obj_material(&w, "dome", 0.3, 0.6, 0.9);
	obj_material(&w, "geo", 1.0, 0.0, 0.0);

	obj_group(&w, "dome");
	obj_use_material(&w, "dome");
	obj_surface_mesh(&w, surf, 32, 32);

	obj_group(&w, "geodesic");
	obj_use_material(&w, "geo");
	{
		qaws_range range = qaws_curve_get_parameter_range(geo);
		for (i = 0; i <= 50; i++)
		{
			qaws_eval_result_3d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)50;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(geo, t, QAWS_EVAL_FLAG_POSITION, &er);
			obj_sphere(&w, er.position, (qaws_scalar)0.03);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/geodesic.obj\n");

	qaws_curve_destroy(geo);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Visual: curve projection onto a dome surface                       */
/* ------------------------------------------------------------------ */
static void visual_obj_projection(void)
{
	obj_writer w;
	qaws_surface* surf = NULL;
	qaws_curve* curve = NULL;
	qaws_curve* uv_curve = NULL;
	qaws_status s;
	unsigned int i;

	printf("visual_obj_projection\n");
	svg_ensure_output_dir();

	/* Biquadratic dome */
	{
		qaws_surface_biquadratic_desc bd;
		bd.control_points[0].x = -1; bd.control_points[0].y = -1; bd.control_points[0].z = 0;
		bd.control_points[1].x =  0; bd.control_points[1].y = -1; bd.control_points[1].z = 0;
		bd.control_points[2].x =  1; bd.control_points[2].y = -1; bd.control_points[2].z = 0;
		bd.control_points[3].x = -1; bd.control_points[3].y =  0; bd.control_points[3].z = 0;
		bd.control_points[4].x =  0; bd.control_points[4].y =  0; bd.control_points[4].z = (qaws_scalar)2.0;
		bd.control_points[5].x =  1; bd.control_points[5].y =  0; bd.control_points[5].z = 0;
		bd.control_points[6].x = -1; bd.control_points[6].y =  1; bd.control_points[6].z = 0;
		bd.control_points[7].x =  0; bd.control_points[7].y =  1; bd.control_points[7].z = 0;
		bd.control_points[8].x =  1; bd.control_points[8].y =  1; bd.control_points[8].z = 0;
		s = qaws_surface_create_biquadratic(&bd, &surf);
		if (s != QAWS_STATUS_OK) return;
	}

	/* 3D cubic S-curve above the dome */
	{
		qaws_scalar pts[12];
		qaws_bezier_desc bd;
		pts[0]  = (qaws_scalar)-0.7; pts[1]  = (qaws_scalar)-0.5; pts[2]  = (qaws_scalar)2.0;
		pts[3]  = (qaws_scalar)-0.2; pts[4]  = (qaws_scalar) 0.8; pts[5]  = (qaws_scalar)2.5;
		pts[6]  = (qaws_scalar) 0.3; pts[7]  = (qaws_scalar)-0.6; pts[8]  = (qaws_scalar)1.8;
		pts[9]  = (qaws_scalar) 0.8; pts[10] = (qaws_scalar) 0.5; pts[11] = (qaws_scalar)2.2;
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_3D;
		bd.degree = 3;
		bd.control_points = pts;
		bd.control_point_count = 4;
		s = qaws_curve_create_bezier(&bd, &curve);
		if (s != QAWS_STATUS_OK) { qaws_surface_destroy(surf); return; }
	}

	/* Project curve onto surface */
	s = qaws_surface_project_curve(surf, curve, 64, &uv_curve);
	if (s != QAWS_STATUS_OK || !uv_curve)
	{
		qaws_curve_destroy(curve); qaws_surface_destroy(surf); return;
	}

	if (!obj_open(&w, OBJ_OUTPUT_DIR "/curve_projection.obj",
		OBJ_OUTPUT_DIR "/curve_projection.mtl"))
	{
		qaws_curve_destroy(uv_curve); qaws_curve_destroy(curve);
		qaws_surface_destroy(surf); return;
	}

	obj_material(&w, "dome", 0.3, 0.6, 0.9);
	obj_material(&w, "source", 1.0, 1.0, 0.0);
	obj_material(&w, "projected", 1.0, 0.0, 0.0);

	/* Surface mesh */
	obj_group(&w, "dome");
	obj_use_material(&w, "dome");
	obj_surface_mesh(&w, surf, 32, 32);

	/* Source 3D curve */
	obj_group(&w, "source_curve");
	obj_use_material(&w, "source");
	{
		qaws_range range = qaws_curve_get_parameter_range(curve);
		for (i = 0; i <= 50; i++)
		{
			qaws_eval_result_3d er;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)50;
			memset(&er, 0, sizeof(er));
			qaws_curve_evaluate_3d(curve, t, QAWS_EVAL_FLAG_POSITION, &er);
			obj_sphere(&w, er.position, (qaws_scalar)0.03);
		}
	}

	/* Projected curve evaluated on the surface */
	obj_group(&w, "projected_curve");
	obj_use_material(&w, "projected");
	{
		qaws_range range = qaws_curve_get_parameter_range(uv_curve);
		for (i = 0; i <= 50; i++)
		{
			qaws_eval_result_2d er_uv;
			qaws_surface_eval_result ser;
			qaws_scalar t = range.min_value + (range.max_value - range.min_value) *
				(qaws_scalar)i / (qaws_scalar)50;
			memset(&er_uv, 0, sizeof(er_uv));
			qaws_curve_evaluate_2d(uv_curve, t, QAWS_EVAL_FLAG_POSITION, &er_uv);
			memset(&ser, 0, sizeof(ser));
			qaws_surface_evaluate(surf, er_uv.position.x, er_uv.position.y,
				QAWS_SURFACE_EVAL_POSITION, &ser);
			obj_sphere(&w, ser.position, (qaws_scalar)0.03);
		}
	}

	obj_close(&w);
	printf("  -> " OBJ_OUTPUT_DIR "/curve_projection.obj\n");

	qaws_curve_destroy(uv_curve);
	qaws_curve_destroy(curve);
	qaws_surface_destroy(surf);
}

/* ------------------------------------------------------------------ */
/*  Main                                                               */
/* ------------------------------------------------------------------ */
int test_41_curve_projection_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_curve_projection();
	test_geodesic_flat();

	/* Visual output */
	visual_obj_geodesic();
	visual_obj_projection();

	printf("41_curve_projection: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
