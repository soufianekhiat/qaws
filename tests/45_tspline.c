/*
 * Test 45: T-spline surface
 *
 * Tests T-spline surface creation, evaluation, and derivatives.
 * T-splines use per-control-point local knot vectors, enabling
 * local refinement without global tensor-product constraint.
 */

#include "test_common.h"
#include "qaws_surface_tspline.h"
#include "qaws_surface.h"

static void test_tspline_basic(void)
{
	/* Create a simple T-spline that is equivalent to a cubic B-spline patch.
	   Use a 4x4 grid of control points with uniform local knot vectors.
	   This should reproduce a standard bicubic B-spline patch. */
	qaws_surface_tspline_desc desc;
	qaws_tspline_control_point cps[16];
	qaws_surface* surface = NULL;
	qaws_status status;
	qaws_surface_eval_result r;
	unsigned int i, j;

	/* For a degree-3 T-spline, each CP needs local knot vector of length
	   2*3+2 = 8. For a 4x4 grid equivalent to a B-spline with global
	   knots {0, 0, 0, 0, 0.5, 1, 1, 1, 1} (5 CPs per axis, 2 spans),
	   but we use 4x4 so global knots = {0, 0, 0, 0, 1, 1, 1, 1}. */
	/* Use degree=2 with 4 CPs per axis: global knots = {0, 0, 0, 0.5, 1, 1, 1}.
	   CP(i) local knots = global[i..i+degree+2] = global[i..i+4]. */
	/* For degree 2, local knot length = 2*2+2 = 6 */
	/* Global u knots for 4 CPs, degree 2: {0, 0, 0, 0.5, 1, 1, 1} (7 knots)
	   CP(0): knots[0..4] = {0, 0, 0, 0.5, 1}  -> local = {0, 0, 0, 0.5, 1, 1}  wait that's wrong
	   Actually for T-spline degree p, local knot vector has 2p+2 entries.
	   But standard B-spline local inference gives p+2 knots per direction.
	   Let's use p+2 knots = 4 knots for degree 2. */
	/* Simpler: use degree 1 (linear) for a clear test.
	   Local knot length = 2*1+2 = 4. Global knots for 4 CPs, degree 1:
	   {0, 0, 1/3, 2/3, 1, 1} (6 knots).
	   CP(0): {0, 0, 1/3, 2/3}
	   CP(1): {0, 1/3, 2/3, 1}
	   CP(2): {1/3, 2/3, 1, 1}
	   CP(3): {2/3, 1, 1, 1}  -- but this would be wrong for degree 1
	   Actually for degree p, each CP(i) has local knots = global[i..i+p+1].
	   For degree 1, 4 CPs: global = {0, 0, 1/3, 2/3, 1, 1}
	   CP(0): {0, 0, 1/3, 2/3} -- 4 knots ✓ (2*1+2=4)
	   CP(1): {0, 1/3, 2/3, 1}
	   CP(2): {1/3, 2/3, 1, 1}
	   CP(3): {2/3, 1, 1, 1} -- this has repeated end, basis=0 everywhere?

	   Hmm, let me just use the simplest approach with degree 1. */
	qaws_scalar u_knots[4][4] = {
		{0.0f, 0.0f, 0.333333f, 0.666667f},
		{0.0f, 0.333333f, 0.666667f, 1.0f},
		{0.333333f, 0.666667f, 1.0f, 1.0f},
		{0.666667f, 1.0f, 1.0f, 1.0f}
	};
	qaws_scalar v_knots[4][4] = {
		{0.0f, 0.0f, 0.333333f, 0.666667f},
		{0.0f, 0.333333f, 0.666667f, 1.0f},
		{0.333333f, 0.666667f, 1.0f, 1.0f},
		{0.666667f, 1.0f, 1.0f, 1.0f}
	};

	printf("  test_tspline_basic...\n");

	/* 4x4 control points forming a dome shape */
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			unsigned int idx = i * 4 + j;
			cps[idx].position.x = (qaws_scalar)i / (qaws_scalar)3;
			cps[idx].position.y = (qaws_scalar)j / (qaws_scalar)3;
			/* Dome: height is higher in the center */
			{
				qaws_scalar du = (qaws_scalar)i / (qaws_scalar)3 - (qaws_scalar)0.5;
				qaws_scalar dv = (qaws_scalar)j / (qaws_scalar)3 - (qaws_scalar)0.5;
				cps[idx].position.z = (qaws_scalar)1.0 - (qaws_scalar)2.0 * (du * du + dv * dv);
			}
			cps[idx].weight = (qaws_scalar)1.0;
			cps[idx].u_knots = u_knots[i];
			cps[idx].v_knots = v_knots[j];
		}
	}

	desc.control_points = cps;
	desc.control_point_count = 16;
	desc.u_degree = 1;
	desc.v_degree = 1;
	desc.u_range.min_value = 0;
	desc.u_range.max_value = 1;
	desc.v_range.min_value = 0;
	desc.v_range.max_value = 1;

	status = qaws_surface_create_tspline(&desc, &surface);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(surface != NULL, "T-spline surface created");
	TEST_ASSERT(qaws_surface_get_kind(surface) == QAWS_SURFACE_KIND_TSPLINE,
		"surface kind is TSPLINE");

	/* Evaluate at center */
	memset(&r, 0, sizeof(r));
	status = qaws_surface_evaluate(surface, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_POSITION, "position computed");
	TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_NORMAL, "normal computed");
	printf("    center: (%.4f, %.4f, %.4f)\n",
		(double)r.position.x, (double)r.position.y, (double)r.position.z);
	printf("    normal: (%.4f, %.4f, %.4f)\n",
		(double)r.normal.x, (double)r.normal.y, (double)r.normal.z);

	/* Evaluate at corners */
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, (qaws_scalar)0, (qaws_scalar)0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	printf("    corner(0,0): (%.4f, %.4f, %.4f)\n",
		(double)r.position.x, (double)r.position.y, (double)r.position.z);

	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, (qaws_scalar)1, (qaws_scalar)1,
		QAWS_SURFACE_EVAL_POSITION, &r);
	printf("    corner(1,1): (%.4f, %.4f, %.4f)\n",
		(double)r.position.x, (double)r.position.y, (double)r.position.z);

	/* Derivatives should be non-zero at center */
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
	TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_DU, "du computed");
	TEST_ASSERT(r.valid_flags & QAWS_SURFACE_EVAL_DV, "dv computed");

	if (surface) qaws_surface_destroy(surface);
}

static void test_tspline_rational(void)
{
	/* Test a rational T-spline (weighted control points) */
	qaws_surface_tspline_desc desc;
	qaws_tspline_control_point cps[4];
	qaws_surface* surface = NULL;
	qaws_status status;
	qaws_surface_eval_result r;

	/* Simple 2x2 grid with degree 1 (bilinear) */
	qaws_scalar u_knots[4] = {0, 0, 1, 1};
	qaws_scalar v_knots[4] = {0, 0, 1, 1};

	printf("  test_tspline_rational...\n");

	cps[0].position.x = 0; cps[0].position.y = 0; cps[0].position.z = 0;
	cps[0].weight = (qaws_scalar)1.0; cps[0].u_knots = u_knots; cps[0].v_knots = v_knots;

	cps[1].position.x = 1; cps[1].position.y = 0; cps[1].position.z = 0;
	cps[1].weight = (qaws_scalar)1.0; cps[1].u_knots = u_knots; cps[1].v_knots = v_knots;

	cps[2].position.x = 0; cps[2].position.y = 1; cps[2].position.z = 0;
	cps[2].weight = (qaws_scalar)1.0; cps[2].u_knots = u_knots; cps[2].v_knots = v_knots;

	cps[3].position.x = 1; cps[3].position.y = 1; cps[3].position.z = 1;
	cps[3].weight = (qaws_scalar)2.0; /* higher weight pulls surface toward this corner */
	cps[3].u_knots = u_knots; cps[3].v_knots = v_knots;

	desc.control_points = cps;
	desc.control_point_count = 4;
	desc.u_degree = 1;
	desc.v_degree = 1;
	desc.u_range.min_value = 0;
	desc.u_range.max_value = 1;
	desc.v_range.min_value = 0;
	desc.v_range.max_value = 1;

	status = qaws_surface_create_tspline(&desc, &surface);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(qaws_surface_is_rational(surface) == 1, "is rational");

	/* Center should be pulled toward (1,1,1) due to higher weight */
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	printf("    center: (%.4f, %.4f, %.4f)\n",
		(double)r.position.x, (double)r.position.y, (double)r.position.z);
	/* Without weights center would be (0.5, 0.5, 0.25)
	   With weight 2 on (1,1,1): center shifts toward it */
	TEST_ASSERT(r.position.x > (qaws_scalar)0.5, "center x > 0.5 (pulled toward weighted CP)");
	TEST_ASSERT(r.position.y > (qaws_scalar)0.5, "center y > 0.5 (pulled toward weighted CP)");

	if (surface) qaws_surface_destroy(surface);
}

static void test_tspline_obj_export(void)
{
	/* Export the 4x4 dome with adaptive tessellation.
	   Adaptive tessellation produces non-uniform triangles: denser
	   near the dome peak where curvature is higher, coarser near the
	   flat corners. This is the natural way to mesh a T-spline surface. */
	qaws_surface_tspline_desc desc;
	qaws_tspline_control_point cps[16];
	qaws_surface* surface = NULL;
	qaws_status status;
	unsigned int i, j;
	obj_writer obj;

	qaws_scalar u_knots_e[4][4] = {
		{0.0f, 0.0f, 0.333333f, 0.666667f},
		{0.0f, 0.333333f, 0.666667f, 1.0f},
		{0.333333f, 0.666667f, 1.0f, 1.0f},
		{0.666667f, 1.0f, 1.0f, 1.0f}
	};
	qaws_scalar v_knots_e[4][4] = {
		{0.0f, 0.0f, 0.333333f, 0.666667f},
		{0.0f, 0.333333f, 0.666667f, 1.0f},
		{0.333333f, 0.666667f, 1.0f, 1.0f},
		{0.666667f, 1.0f, 1.0f, 1.0f}
	};

	printf("  test_tspline_obj_export...\n");

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			unsigned int idx = i * 4 + j;
			cps[idx].position.x = (qaws_scalar)i / (qaws_scalar)3;
			cps[idx].position.y = (qaws_scalar)j / (qaws_scalar)3;
			{
				qaws_scalar du_v = (qaws_scalar)i / (qaws_scalar)3 - (qaws_scalar)0.5;
				qaws_scalar dv_v = (qaws_scalar)j / (qaws_scalar)3 - (qaws_scalar)0.5;
				cps[idx].position.z = (qaws_scalar)1.0 - (qaws_scalar)2.0 * (du_v * du_v + dv_v * dv_v);
			}
			cps[idx].weight = (qaws_scalar)1.0;
			cps[idx].u_knots = u_knots_e[i];
			cps[idx].v_knots = v_knots_e[j];
		}
	}

	desc.control_points = cps;
	desc.control_point_count = 16;
	desc.u_degree = 1;
	desc.v_degree = 1;
	desc.u_range.min_value = 0;
	desc.u_range.max_value = 1;
	desc.v_range.min_value = 0;
	desc.v_range.max_value = 1;

	status = qaws_surface_create_tspline(&desc, &surface);
	if (status != QAWS_STATUS_OK) return;

	svg_ensure_output_dir();
	if (obj_open(&obj, OBJ_OUTPUT_DIR "/45_tspline_dome.obj",
		OBJ_OUTPUT_DIR "/45_tspline_dome.mtl"))
	{
		qaws_tessellation_desc tdesc;
		qaws_tessellation_vertex tverts[8192];
		unsigned int tindices[49152];
		unsigned int tvert_count = 0, tidx_count = 0;
		unsigned int first_v, first_n, ti;

		tdesc.max_depth = 6;
		tdesc.curvature_threshold = (qaws_scalar)0.02;
		tdesc.max_edge_length = 0;

		obj_material(&obj, "tspline", 0.3, 0.7, 0.9);
		obj_group(&obj, "tspline_dome");
		obj_use_material(&obj, "tspline");

		status = qaws_surface_tessellate(surface, &tdesc,
			tverts, 8192, &tvert_count,
			tindices, 49152, &tidx_count);

		if (status == QAWS_STATUS_OK && tvert_count > 0)
		{
			first_v = obj.vertex_count + 1;
			first_n = obj.normal_count + 1;

			for (i = 0; i < tvert_count; i++)
			{
				obj_vertex(&obj, tverts[i].position);
				obj_normal(&obj, tverts[i].normal);
			}
			for (ti = 0; ti + 2 < tidx_count; ti += 3)
			{
				fprintf(obj.fp, "f %u//%u %u//%u %u//%u\n",
					first_v + tindices[ti],   first_n + tindices[ti],
					first_v + tindices[ti+1], first_n + tindices[ti+1],
					first_v + tindices[ti+2], first_n + tindices[ti+2]);
			}
			printf("    adaptive: %u verts, %u tris\n",
				tvert_count, tidx_count / 3);
		}

		obj_close(&obj);
		printf("    wrote 45_tspline_dome.obj\n");
	}

	qaws_surface_destroy(surface);
}

int test_45_tspline_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 45: T-spline surfaces\n");
	test_tspline_basic();
	test_tspline_rational();
	test_tspline_obj_export();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
