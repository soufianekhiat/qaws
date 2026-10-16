#include "test_common.h"

static int scalar_is_finite(qaws_scalar v)
{
	return (v == v) && (v - v == 0);
}

struct stream_counter { unsigned int count; };

static int stream_count_cb_2d(
	qaws_scalar param, qaws_vec2 const* pos, void* ud)
{
	struct stream_counter* c = (struct stream_counter*)ud;
	(void)param; (void)pos;
	c->count++;
	return 1;
}

static int stream_count_cb_3d(
	qaws_scalar param, qaws_vec3 const* pos, void* ud)
{
	struct stream_counter* c = (struct stream_counter*)ud;
	(void)param; (void)pos;
	c->count++;
	return 1;
}

static int stream_early_stop_cb_2d(
	qaws_scalar param, qaws_vec2 const* pos, void* ud)
{
	struct stream_counter* c = (struct stream_counter*)ud;
	(void)param; (void)pos;
	c->count++;
	return (c->count < 3) ? 1 : 0;
}

static void test_sampling(void)
{
	printf("test_sampling\n");

	qaws_vec2 points[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc;
	desc.dimension = QAWS_DIMENSION_2D;
	desc.degree = 1;
	desc.control_points = points;
	desc.control_point_count = 2;

	qaws_curve* curve = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve, &desc);
	TEST_ASSERT_STATUS(s);

	qaws_sampling_desc sdesc;
	sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
	sdesc.sample_count = 11;
	sdesc.error_tolerance = 0;
	sdesc.use_adaptive_sampling = 0;

	qaws_vec2 samples[11];
	unsigned int count = 0;
	s = qaws_curve_sample_2d(samples, &count, curve, 11, &sdesc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count == 11, "sample count == 11");
	TEST_ASSERT(approx_eq(samples[0].x, (qaws_scalar)0.0), "first sample x");
	TEST_ASSERT(approx_eq(samples[10].x, (qaws_scalar)1.0), "last sample x");
	TEST_ASSERT(approx_eq(samples[5].x, (qaws_scalar)0.5), "mid sample x");

	qaws_curve_destroy(curve);
}

static void test_sampling_extended(void)
{
	printf("test_sampling_extended\n");

	/* 2D: get_sample_count */
	qaws_vec2 points2d[] = { {0, 0}, {1, 1} };
	qaws_bezier_desc desc2;
	desc2.dimension = QAWS_DIMENSION_2D;
	desc2.degree = 1;
	desc2.control_points = points2d;
	desc2.control_point_count = 2;

	qaws_curve* curve2 = NULL;
	qaws_status s = qaws_curve_create_bezier(&curve2, &desc2);
	TEST_ASSERT_STATUS(s);

	qaws_sampling_desc sdesc;
	sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
	sdesc.sample_count = 5;
	sdesc.error_tolerance = 0;
	sdesc.use_adaptive_sampling = 0;

	unsigned int count = 0;
	s = qaws_curve_get_sample_count(&count, curve2, &sdesc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count == 5, "get_sample_count == 5");

	/* sample_eval_2d */
	qaws_eval_result_2d evals[5];
	unsigned int eval_count = 0;
	s = qaws_curve_sample_eval_2d(evals, &eval_count, curve2, 5, &sdesc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(eval_count == 5, "sample_eval count == 5");
	TEST_ASSERT(approx_eq(evals[0].position.x, (qaws_scalar)0.0), "eval sample 0 x");
	TEST_ASSERT(approx_eq(evals[4].position.x, (qaws_scalar)1.0), "eval sample 4 x");

	qaws_curve_destroy(curve2);

	/* 3D: sample_3d */
	qaws_vec3 points3d[] = { {0, 0, 0}, {1, 1, 1} };
	qaws_bezier_desc desc3;
	desc3.dimension = QAWS_DIMENSION_3D;
	desc3.degree = 1;
	desc3.control_points = points3d;
	desc3.control_point_count = 2;

	qaws_curve* curve3 = NULL;
	s = qaws_curve_create_bezier(&curve3, &desc3);
	TEST_ASSERT_STATUS(s);

	qaws_vec3 samples3[5];
	unsigned int count3 = 0;
	s = qaws_curve_sample_3d(samples3, &count3, curve3, 5, &sdesc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(count3 == 5, "sample_3d count == 5");
	TEST_ASSERT(approx_eq(samples3[0].z, (qaws_scalar)0.0), "3d sample start z");
	TEST_ASSERT(approx_eq(samples3[4].z, (qaws_scalar)1.0), "3d sample end z");

	/* sample_eval_3d */
	qaws_eval_result_3d evals3[5];
	unsigned int eval3_count = 0;
	s = qaws_curve_sample_eval_3d(evals3, &eval3_count, curve3, 5, &sdesc);
	TEST_ASSERT_STATUS(s);
	TEST_ASSERT(eval3_count == 5, "sample_eval_3d count == 5");

	qaws_curve_destroy(curve3);
}

static void test_adaptive_eval_sampling(void)
{
	printf("test_adaptive_eval_sampling\n");

	qaws_vec2 pts[] = { {0,0}, {0,4}, {4,4}, {4,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 3;
	bdesc.control_points = pts;
	bdesc.control_point_count = 4;

	status = qaws_curve_create_bezier(&curve, &bdesc);
	TEST_ASSERT_STATUS(status);

	/* Uniform eval sampling baseline */
	{
		qaws_sampling_desc sdesc;
		qaws_eval_result_2d samples[2048];
		unsigned int count = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 0;

		status = qaws_curve_sample_eval_2d(
			samples, &count, curve, 2048, &sdesc);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count == 8, "uniform eval gives 8");
	}

	/* Adaptive eval sampling should produce more points and include derivatives */
	{
		qaws_sampling_desc sdesc;
		qaws_eval_result_2d samples[2048];
		unsigned int count = 0;
		int has_d1 = 1;
		unsigned int i;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.001;

		status = qaws_curve_sample_eval_2d(
			samples, &count, curve, 2048, &sdesc);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count >= 8, "adaptive eval >= initial");

		/* Check that derivatives are populated */
		for (i = 0; i < count && i < 5; ++i)
		{
			if (samples[i].d1.x == (qaws_scalar)0.0 &&
				samples[i].d1.y == (qaws_scalar)0.0)
				has_d1 = 0;
		}
		TEST_ASSERT(has_d1, "adaptive eval has D1 data");
	}

	/* 3D adaptive eval */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {0,4,2}, {4,4,2}, {4,0,0} };
		qaws_curve *curve3d = NULL;
		qaws_eval_result_3d samples3d[2048];
		unsigned int count3d = 0;
		qaws_sampling_desc sdesc;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.control_points = pts3d;

		status = qaws_curve_create_bezier(&curve3d, &bdesc);
		TEST_ASSERT_STATUS(status);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.01;

		status = qaws_curve_sample_eval_3d(
			samples3d, &count3d, curve3d, 2048, &sdesc);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count3d >= 8, "3D adaptive eval >= initial");

		qaws_curve_destroy(curve3d);
	}

	/* Buffer too small: adaptive sampling into a tiny buffer */
	{
		qaws_sampling_desc sdesc;
		qaws_eval_result_2d tiny_buf[2];
		unsigned int count_needed = 0;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 64;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.0001;

		status = qaws_curve_sample_eval_2d(
			tiny_buf, &count_needed, curve, 2, &sdesc);
		TEST_ASSERT(status == QAWS_STATUS_BUFFER_TOO_SMALL,
			"tiny buffer returns BUFFER_TOO_SMALL");
		TEST_ASSERT(count_needed > 2,
			"buffer too small reports needed count > capacity");
	}

	qaws_curve_destroy(curve);
}

static void test_adaptive_sampling(void)
{
	printf("test_adaptive_sampling\n");

	/* Create a curve with high curvature variation: cubic Bezier */
	qaws_vec2 pts[] = { {0,0}, {0,4}, {4,4}, {4,0} };
	qaws_bezier_desc bdesc;
	qaws_curve *curve = NULL;
	qaws_status status;
	qaws_sampling_desc sdesc;
	qaws_vec2 positions[2048];
	unsigned int count = 0;

	bdesc.dimension = QAWS_DIMENSION_2D;
	bdesc.degree = 3;
	bdesc.control_points = pts;
	bdesc.control_point_count = 4;

	status = qaws_curve_create_bezier(&curve, &bdesc);
	TEST_ASSERT_STATUS(status);

	/* Uniform sampling for baseline */
	memset(&sdesc, 0, sizeof(sdesc));
	sdesc.sample_count = 8;
	sdesc.use_adaptive_sampling = 0;

	status = qaws_curve_sample_2d(positions, &count, curve, 2048, &sdesc);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(count == 8, "uniform gives exactly 8 samples");

	/* Adaptive sampling with tight tolerance should produce more points */
	memset(&sdesc, 0, sizeof(sdesc));
	sdesc.sample_count = 8;
	sdesc.use_adaptive_sampling = 1;
	sdesc.error_tolerance = (qaws_scalar)0.001;

	status = qaws_curve_sample_2d(positions, &count, curve, 2048, &sdesc);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(count >= 8, "adaptive produces >= initial samples");

	/* Adaptive with loose tolerance should produce fewer subdivisions */
	{
		unsigned int count_tight = count;
		unsigned int count_loose;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)1.0;

		status = qaws_curve_sample_2d(
			positions, &count_loose, curve, 2048, &sdesc);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count_loose <= count_tight,
			"loose tolerance <= tight tolerance samples");
	}

	/* Verify all points lie on the curve (roughly) */
	{
		unsigned int i;
		int all_valid = 1;
		for (i = 0; i < count; ++i)
		{
			if (positions[i].x < (qaws_scalar)-0.1 ||
				positions[i].x > (qaws_scalar)4.1 ||
				positions[i].y < (qaws_scalar)-0.1 ||
				positions[i].y > (qaws_scalar)4.1)
			{
				all_valid = 0;
				break;
			}
		}
		TEST_ASSERT(all_valid, "adaptive points within curve bounds");
	}

	/* Test 3D adaptive sampling */
	{
		qaws_vec3 pts3d[] = { {0,0,0}, {0,4,2}, {4,4,2}, {4,0,0} };
		qaws_curve *curve3d = NULL;
		qaws_vec3 positions3d[2048];
		unsigned int count3d = 0;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts3d;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&curve3d, &bdesc);
		TEST_ASSERT_STATUS(status);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.sample_count = 8;
		sdesc.use_adaptive_sampling = 1;
		sdesc.error_tolerance = (qaws_scalar)0.01;

		status = qaws_curve_sample_3d(
			positions3d, &count3d, curve3d, 2048, &sdesc);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(count3d >= 8, "3D adaptive produces >= initial samples");

		qaws_curve_destroy(curve3d);
	}

	qaws_curve_destroy(curve);
}

static void test_curvature_weighted_sampling(void)
{
	printf("test_curvature_weighted_sampling\n");

	/* 1. 2D curvature-weighted: S-curve shape with varying curvature */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 3}, {3, -3}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[20];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_2d(samples, &count,
			curve, 20, 20);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 20, "curvature weighted 2d count == 20");

		/* All positions are finite */
		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"curv weighted 2d sample x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"curv weighted 2d sample y finite");
		}

		/* Spacing should be non-uniform: check that not all intervals
		   are the same length. Compute min/max spacing. */
		{
			qaws_scalar min_gap = (qaws_scalar)1e10;
			qaws_scalar max_gap = (qaws_scalar)0.0;
			for (i = 1; i < count; i++) {
				qaws_scalar dx = samples[i].x - samples[i - 1].x;
				qaws_scalar dy = samples[i].y - samples[i - 1].y;
				qaws_scalar gap = dx * dx + dy * dy;
				if (gap < min_gap) min_gap = gap;
				if (gap > max_gap) max_gap = gap;
			}
			TEST_ASSERT(max_gap > min_gap * (qaws_scalar)1.5,
				"curvature weighted spacing is non-uniform");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. 3D curvature-weighted */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,3,1}, {3,-3,2}, {4,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec3 samples[15];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_3d(samples, &count,
			curve, 15, 15);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count == 15, "curvature weighted 3d count == 15");

		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"curv weighted 3d x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"curv weighted 3d y finite");
			TEST_ASSERT(scalar_is_finite(samples[i].z),
				"curv weighted 3d z finite");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Null args */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[10];
		unsigned int count = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_2d(samples, &count,
			NULL, 10, 10);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"curv weighted null curve => INVALID_ARGUMENT");

		s = qaws_curve_sample_curvature_weighted_2d(NULL, &count,
			curve, 10, 10);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"curv weighted null out => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}

	/* 4. sample_count < 2 */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[10];
		unsigned int count = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_curvature_weighted_2d(samples, &count,
			curve, 1, 10);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"curv weighted count<2 => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}
}

static void test_feature_preserving_sampling(void)
{
	printf("test_feature_preserving_sampling\n");

	/* 1. 2D S-curve with inflection: CPs (0,0),(1,3),(2,-3),(3,0) */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 3}, {2, -3}, {3, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[50];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_feature_preserving_2d(samples, &count,
			curve, 5, 50);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 5,
			"feature preserving 2d count >= base count");

		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"feature preserving 2d x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"feature preserving 2d y finite");
		}

		qaws_curve_destroy(curve);
	}

	/* 2. 3D curve */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,3,1}, {2,-3,2}, {3,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec3 samples[50];
		unsigned int count = 0;
		unsigned int i;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_feature_preserving_3d(samples, &count,
			curve, 5, 50);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(count >= 5,
			"feature preserving 3d count >= base count");

		for (i = 0; i < count; i++) {
			TEST_ASSERT(scalar_is_finite(samples[i].x),
				"feature preserving 3d x finite");
			TEST_ASSERT(scalar_is_finite(samples[i].y),
				"feature preserving 3d y finite");
			TEST_ASSERT(scalar_is_finite(samples[i].z),
				"feature preserving 3d z finite");
		}

		qaws_curve_destroy(curve);
	}

	/* 3. Null args */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_vec2 samples[10];
		unsigned int count = 0;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		s = qaws_curve_sample_feature_preserving_2d(samples, &count,
			NULL, 5, 50);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"feature preserving null curve => INVALID_ARGUMENT");

		s = qaws_curve_sample_feature_preserving_2d(NULL, &count,
			curve, 5, 50);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"feature preserving null out => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}
}

static void test_streaming_sampling(void)
{
	printf("test_streaming_sampling\n");

	/* 1. 2D streaming: count invocations */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {3, 1}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;
		struct stream_counter counter;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 10;
		sdesc.error_tolerance = 0;
		sdesc.use_adaptive_sampling = 0;

		counter.count = 0;
		s = qaws_curve_sample_streaming_2d(curve, &sdesc,
			stream_count_cb_2d, &counter);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(counter.count == 10,
			"streaming 2d callback count == 10");

		qaws_curve_destroy(curve);
	}

	/* 2. Early termination: callback returns 0 after 3 calls */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 2}, {3, 1}, {4, 0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;
		struct stream_counter counter;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 20;
		sdesc.error_tolerance = 0;
		sdesc.use_adaptive_sampling = 0;

		counter.count = 0;
		s = qaws_curve_sample_streaming_2d(curve, &sdesc,
			stream_early_stop_cb_2d, &counter);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(counter.count == 3,
			"streaming early stop at 3 calls");

		qaws_curve_destroy(curve);
	}

	/* 3. 3D streaming */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,2,1}, {3,1,2}, {4,0,0} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;
		struct stream_counter counter;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 8;
		sdesc.error_tolerance = 0;
		sdesc.use_adaptive_sampling = 0;

		counter.count = 0;
		s = qaws_curve_sample_streaming_3d(curve, &sdesc,
			stream_count_cb_3d, &counter);
		TEST_ASSERT_STATUS(s);
		TEST_ASSERT(counter.count == 8,
			"streaming 3d callback count == 8");

		qaws_curve_destroy(curve);
	}

	/* 4. Null callback => INVALID_ARGUMENT */
	{
		qaws_vec2 pts[] = { {0, 0}, {1, 1}, {2, 0}, {3, 1} };
		qaws_bezier_desc bdesc;
		qaws_curve* curve = NULL;
		qaws_status s;
		qaws_sampling_desc sdesc;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		s = qaws_curve_create_bezier(&curve, &bdesc);
		TEST_ASSERT_STATUS(s);

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.traversal_mode = QAWS_TRAVERSAL_MODE_PARAMETER;
		sdesc.sample_count = 10;

		s = qaws_curve_sample_streaming_2d(curve, &sdesc, NULL, NULL);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT,
			"streaming null callback => INVALID_ARGUMENT");

		qaws_curve_destroy(curve);
	}
}

int test_15_sampling_main(void)
{
	g_pass = 0; g_fail = 0;

	test_sampling();
	test_sampling_extended();
	test_adaptive_eval_sampling();
	test_adaptive_sampling();
	test_curvature_weighted_sampling();
	test_feature_preserving_sampling();
	test_streaming_sampling();

	printf("15_sampling: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
