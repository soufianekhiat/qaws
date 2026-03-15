#include "test_common.h"

static void test_svg_batch_eval(void)
{
	/*
	 * SVG showing batch evaluation:
	 * - A Bezier curve drawn via single-point sampling (cyan)
	 * - Batch-evaluated points overlaid as dots (magenta)
	 * - Tangent arrows at each batch point (gold)
	 * - A second curve (Catmull-Rom) with batch-sampled positions (green dots)
	 * This visually proves batch eval produces correct results.
	 */
	svg_writer svg;
	printf("test_svg_batch_eval\n");

	if (!svg_open(&svg, SVG_OUTPUT_DIR "/23_batch_eval.svg",
		(qaws_scalar)-0.5, (qaws_scalar)-1.0, (qaws_scalar)6.0, (qaws_scalar)4.0,
		(qaws_scalar)700, (qaws_scalar)500))
		return;

	/* --- Bezier curve with batch eval points + tangents --- */
	{
		qaws_vec2 cp[] = {{0,0}, {1,3}, {3,3}, {5,0}};
		qaws_bezier_desc desc;
		qaws_curve* curve = NULL;
		qaws_vec2 curve_pts[SVG_SAMPLES];
		unsigned int n;
		unsigned int batch_n = 16;
		qaws_scalar params[16];
		qaws_eval_result_2d results[16];
		qaws_vec2 batch_dots[16];
		unsigned int i;
		qaws_scalar arrow_scale = 0.15f;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;
		qaws_curve_create_bezier(&desc, &curve);

		/* Draw curve via standard sampling */
		n = svg_sample_curve(curve, curve_pts, SVG_SAMPLES);
		svg_polyline(&svg, curve_pts, n, "#66d9ef", 2.0f);

		/* Draw control polygon (dashed) */
		svg_polyline_dashed(&svg, cp, 4, "#555577", 1.0f, "4,4");
		svg_dots(&svg, cp, 4, "#888899", 3.0f);

		/* Generate batch parameters (uniform) */
		for (i = 0; i < batch_n; i++)
			params[i] = (qaws_scalar)i / (qaws_scalar)(batch_n - 1);

		/* Batch evaluate with position + D1 */
		qaws_curve_evaluate_batch_2d(curve, params, batch_n,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, results);

		/* Draw batch points as dots */
		for (i = 0; i < batch_n; i++)
		{
			batch_dots[i].x = results[i].position.x;
			batch_dots[i].y = results[i].position.y;
		}
		svg_dots(&svg, batch_dots, batch_n, "#ff79c6", 4.0f);

		/* Draw tangent arrows at each batch point */
		for (i = 0; i < batch_n; i++)
		{
			qaws_scalar tx = results[i].d1.x * arrow_scale;
			qaws_scalar ty = results[i].d1.y * arrow_scale;
			svg_line(&svg,
				results[i].position.x, results[i].position.y,
				results[i].position.x + tx, results[i].position.y + ty,
				"#f1fa8c", 1.5f);
		}

		svg_label(&svg, 0.0f, 3.3f, "Bezier: batch pos (magenta) + tangents (gold)", "#8888aa");
		qaws_curve_destroy(curve);
	}

	/* --- Catmull-Rom with batch eval --- */
	{
		qaws_vec2 cr_cp[] = {{0,-0.5f}, {1,1.5f}, {2.5f,2.0f}, {3.5f,0.5f}, {5,1.0f}};
		qaws_catmull_rom_desc desc;
		qaws_curve* curve = NULL;
		qaws_vec2 curve_pts[SVG_SAMPLES];
		unsigned int n;
		unsigned int batch_n = 20;
		qaws_scalar params[20];
		qaws_eval_result_2d results[20];
		qaws_vec2 batch_dots[20];
		qaws_range range;
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.control_points = cr_cp;
		desc.control_point_count = 5;
		desc.parameterization = QAWS_PARAMETERIZATION_CENTRIPETAL;
		desc.closed = 0;
		qaws_curve_create_catmull_rom(&desc, &curve);

		/* Draw curve */
		n = svg_sample_curve(curve, curve_pts, SVG_SAMPLES);
		svg_polyline(&svg, curve_pts, n, "#50fa7b", 2.0f);

		/* Batch evaluate along the curve */
		range = qaws_curve_get_parameter_range(curve);
		for (i = 0; i < batch_n; i++)
			params[i] = range.min_value +
				(range.max_value - range.min_value) * (qaws_scalar)i / (qaws_scalar)(batch_n - 1);

		qaws_curve_evaluate_batch_2d(curve, params, batch_n,
			QAWS_EVAL_FLAG_POSITION, results);

		for (i = 0; i < batch_n; i++)
		{
			batch_dots[i].x = results[i].position.x;
			batch_dots[i].y = results[i].position.y;
		}
		svg_dots(&svg, batch_dots, batch_n, "#bd93f9", 3.5f);
		svg_dots(&svg, cr_cp, 5, "#ffffff", 3.0f);

		svg_label(&svg, 0.0f, -0.8f, "Catmull-Rom: batch positions (purple)", "#8888aa");
		qaws_curve_destroy(curve);
	}

	svg_close(&svg);
	printf("  -> " SVG_OUTPUT_DIR "/23_batch_eval.svg\n");
}

static void test_batch_eval(void)
{
	/* Test batch evaluation of 2D and 3D curves */
	qaws_curve* curve2d = NULL;
	qaws_curve* curve3d = NULL;
	qaws_status s;

	printf("test_batch_eval\n");

	/* 2D Bezier */
	{
		qaws_vec2 pts[] = {{0, 0}, {1, 2}, {3, 1}, {4, 0}};
		qaws_bezier_desc desc;
		qaws_scalar params[5] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
		qaws_eval_result_2d results[5];
		qaws_eval_result_2d single;
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;
		s = qaws_curve_create_bezier(&desc, &curve2d);
		TEST_ASSERT(s == QAWS_STATUS_OK, "batch: create 2d bezier");

		s = qaws_curve_evaluate_batch_2d(curve2d, params, 5,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, results);
		TEST_ASSERT(s == QAWS_STATUS_OK, "batch: batch eval 2d");

		/* Verify each result matches single eval */
		for (i = 0; i < 5; i++)
		{
			qaws_curve_evaluate_2d(curve2d, params[i],
				QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &single);
			TEST_ASSERT(
				(float)fabs(results[i].position.x - single.position.x) < 1e-6f &&
				(float)fabs(results[i].position.y - single.position.y) < 1e-6f,
				"batch: 2d position matches single");
			TEST_ASSERT(
				(float)fabs(results[i].d1.x - single.d1.x) < 1e-6f &&
				(float)fabs(results[i].d1.y - single.d1.y) < 1e-6f,
				"batch: 2d d1 matches single");
		}
		qaws_curve_destroy(curve2d);
	}

	/* 3D Bezier */
	{
		qaws_vec3 pts[] = {{0, 0, 0}, {1, 2, 1}, {3, 1, 2}, {4, 0, 0}};
		qaws_bezier_desc desc;
		qaws_scalar params[4] = {0.0f, 0.33f, 0.66f, 1.0f};
		qaws_eval_result_3d results[4];
		qaws_eval_result_3d single;
		unsigned int i;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;
		s = qaws_curve_create_bezier(&desc, &curve3d);
		TEST_ASSERT(s == QAWS_STATUS_OK, "batch: create 3d bezier");

		s = qaws_curve_evaluate_batch_3d(curve3d, params, 4,
			QAWS_EVAL_FLAG_POSITION, results);
		TEST_ASSERT(s == QAWS_STATUS_OK, "batch: batch eval 3d");

		for (i = 0; i < 4; i++)
		{
			qaws_curve_evaluate_3d(curve3d, params[i],
				QAWS_EVAL_FLAG_POSITION, &single);
			TEST_ASSERT(
				(float)fabs(results[i].position.x - single.position.x) < 1e-6f &&
				(float)fabs(results[i].position.y - single.position.y) < 1e-6f &&
				(float)fabs(results[i].position.z - single.position.z) < 1e-6f,
				"batch: 3d position matches single");
		}
		qaws_curve_destroy(curve3d);
	}

	/* Batch surface eval */
	{
		qaws_surface* surf = NULL;
		qaws_surface_bezier_desc sdesc;
		qaws_vec3 cps[4] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,1}};
		qaws_scalar u_params[3] = {0.0f, 0.5f, 1.0f};
		qaws_scalar v_params[3] = {0.0f, 0.5f, 1.0f};
		qaws_surface_eval_result sresults[3];
		qaws_surface_eval_result ssingle;
		unsigned int i;

		sdesc.u_degree = 1;
		sdesc.v_degree = 1;
		sdesc.u_point_count = 2;
		sdesc.v_point_count = 2;
		sdesc.control_points = cps;
		s = qaws_surface_create_bezier(&sdesc, &surf);
		TEST_ASSERT(s == QAWS_STATUS_OK, "batch: create surface");

		s = qaws_surface_evaluate_batch(surf, u_params, v_params, 3,
			QAWS_SURFACE_EVAL_POSITION, sresults);
		TEST_ASSERT(s == QAWS_STATUS_OK, "batch: surface batch eval");

		for (i = 0; i < 3; i++)
		{
			qaws_surface_evaluate(surf, u_params[i], v_params[i],
				QAWS_SURFACE_EVAL_POSITION, &ssingle);
			TEST_ASSERT(
				(float)fabs(sresults[i].position.x - ssingle.position.x) < 1e-6f &&
				(float)fabs(sresults[i].position.y - ssingle.position.y) < 1e-6f &&
				(float)fabs(sresults[i].position.z - ssingle.position.z) < 1e-6f,
				"batch: surface position matches single");
		}
		qaws_surface_destroy(surf);
	}

	/* Error handling */
	{
		qaws_scalar params[2] = {0.0f, 1.0f};
		qaws_eval_result_2d res[2];
		s = qaws_curve_evaluate_batch_2d(NULL, params, 2,
			QAWS_EVAL_FLAG_POSITION, res);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "batch: null curve");

		s = qaws_surface_evaluate_batch(NULL, params, params, 2,
			QAWS_SURFACE_EVAL_POSITION, NULL);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "batch: null surface");
	}
}

/* Custom allocator tracking for testing */
typedef struct test_alloc_stats
{
	unsigned int alloc_count;
	unsigned int dealloc_count;
	unsigned long total_bytes;
} test_alloc_stats;

static void* test_alloc_fn(unsigned long size, void* user_data)
{
	test_alloc_stats* stats = (test_alloc_stats*)user_data;
	stats->alloc_count++;
	stats->total_bytes += size;
	return malloc((size_t)size);
}

static void test_dealloc_fn(void* ptr, void* user_data)
{
	test_alloc_stats* stats = (test_alloc_stats*)user_data;
	stats->dealloc_count++;
	free(ptr);
}

static void test_custom_allocator(void)
{
	test_alloc_stats stats;
	qaws_allocator alloc;
	qaws_curve* curve = NULL;
	qaws_status s;

	printf("test_custom_allocator\n");

	memset(&stats, 0, sizeof(stats));
	alloc.alloc = test_alloc_fn;
	alloc.dealloc = test_dealloc_fn;
	alloc.user_data = &stats;

	/* Create bezier with custom allocator */
	{
		qaws_vec3 pts[] = {{0,0,0}, {1,2,1}, {3,1,2}, {4,0,0}};
		qaws_bezier_desc desc;
		qaws_eval_result_3d result;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_create_bezier_ex(&desc, &alloc, &curve);
		TEST_ASSERT(s == QAWS_STATUS_OK, "allocator: create bezier_ex");
		TEST_ASSERT(stats.alloc_count > 0, "allocator: allocs happened");
		TEST_ASSERT(stats.total_bytes > 0, "allocator: bytes tracked");

		/* Evaluate should work normally */
		s = qaws_curve_evaluate_3d(curve, 0.5f, QAWS_EVAL_FLAG_POSITION, &result);
		TEST_ASSERT(s == QAWS_STATUS_OK, "allocator: eval works");
		TEST_ASSERT(result.valid_flags & QAWS_EVAL_FLAG_POSITION,
			"allocator: position valid");

		{
			unsigned int alloc_before = stats.alloc_count;
			unsigned int dealloc_before = stats.dealloc_count;
			qaws_curve_destroy(curve);
			TEST_ASSERT(stats.dealloc_count > dealloc_before,
				"allocator: deallocs on destroy");
			TEST_ASSERT(stats.alloc_count == alloc_before,
				"allocator: no allocs on destroy");
			/* Verify all allocs are matched by deallocs */
			TEST_ASSERT(stats.alloc_count == stats.dealloc_count,
				"allocator: alloc/dealloc balanced");
		}
	}

	/* Test with hermite _ex */
	memset(&stats, 0, sizeof(stats));
	{
		qaws_vec2 pts[] = {{0, 0}, {4, 4}};
		qaws_vec2 der[] = {{2, 0}, {2, 0}};
		qaws_hermite_desc desc;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.point_count = 2;
		desc.derivatives = der;
		desc.derivative_count = 2;

		s = qaws_curve_create_hermite_ex(&desc, &alloc, &curve);
		TEST_ASSERT(s == QAWS_STATUS_OK, "allocator: create hermite_ex");
		TEST_ASSERT(stats.alloc_count > 0, "allocator: hermite allocs");
		qaws_curve_destroy(curve);
		TEST_ASSERT(stats.alloc_count == stats.dealloc_count,
			"allocator: hermite balanced");
	}

	/* Test with NULL allocator (should use default) */
	{
		qaws_vec2 pts[] = {{0, 0}, {1, 1}, {2, 0}};
		qaws_bezier_desc desc;
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		s = qaws_curve_create_bezier_ex(&desc, NULL, &curve);
		TEST_ASSERT(s == QAWS_STATUS_OK, "allocator: null allocator OK");
		qaws_curve_destroy(curve);
	}
}

static void test_inline_curves(void)
{
	qaws_status s;

	printf("test_inline_curves\n");

	/* Inline bezier 2D */
	{
		qaws_curve_inline ic;
		qaws_curve* curve;
		qaws_eval_result_2d result;
		qaws_vec2 pts[] = {{0, 0}, {1, 2}, {3, 1}, {4, 0}};
		qaws_bezier_desc desc;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_init_bezier_inline(&desc, &ic);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: init bezier 2d");

		curve = qaws_curve_inline_get_curve(&ic);
		TEST_ASSERT(curve != NULL, "inline: get_curve non-null");

		s = qaws_curve_evaluate_2d(curve, 0.0f,
			QAWS_EVAL_FLAG_POSITION, &result);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: eval at t=0");
		TEST_ASSERT((float)fabs(result.position.x - 0.0f) < 1e-6f &&
			(float)fabs(result.position.y - 0.0f) < 1e-6f,
			"inline: 2d start point");

		s = qaws_curve_evaluate_2d(curve, 1.0f,
			QAWS_EVAL_FLAG_POSITION, &result);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: eval at t=1");
		TEST_ASSERT((float)fabs(result.position.x - 4.0f) < 1e-6f &&
			(float)fabs(result.position.y - 0.0f) < 1e-6f,
			"inline: 2d end point");

		/* Mid-point evaluation with D1 */
		s = qaws_curve_evaluate_2d(curve, 0.5f,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: eval at t=0.5 with d1");
		TEST_ASSERT(result.valid_flags & QAWS_EVAL_FLAG_D1,
			"inline: d1 valid at 0.5");
		/* No qaws_curve_destroy needed! */
	}

	/* Inline bezier 3D */
	{
		qaws_curve_inline ic;
		qaws_curve* curve;
		qaws_eval_result_3d result;
		qaws_vec3 pts[] = {{0,0,0}, {1,2,1}, {3,1,2}, {4,0,0}};
		qaws_bezier_desc desc;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 3;
		desc.control_points = pts;
		desc.control_point_count = 4;

		s = qaws_curve_init_bezier_inline(&desc, &ic);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: init bezier 3d");

		curve = qaws_curve_inline_get_curve(&ic);
		s = qaws_curve_evaluate_3d(curve, 0.0f,
			QAWS_EVAL_FLAG_POSITION, &result);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: 3d eval");
		TEST_ASSERT((float)fabs(result.position.x - 0.0f) < 1e-6f &&
			(float)fabs(result.position.z - 0.0f) < 1e-6f,
			"inline: 3d start point");

		s = qaws_curve_evaluate_3d(curve, 1.0f,
			QAWS_EVAL_FLAG_POSITION, &result);
		TEST_ASSERT((float)fabs(result.position.x - 4.0f) < 1e-6f,
			"inline: 3d end x");
	}

	/* Inline polynomial 2D */
	{
		qaws_curve_inline ic;
		qaws_curve* curve;
		qaws_eval_result_2d result;
		/* P(t) = (1, 0) + (2, 3)*t + (0, -1)*t^2 */
		qaws_scalar coeffs[] = {1.0f, 0.0f, 2.0f, 3.0f, 0.0f, -1.0f};
		qaws_polynomial_desc desc;

		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.coefficients = coeffs;
		desc.coefficient_count = 3;
		desc.t_min = 0; desc.t_max = 1;

		s = qaws_curve_init_polynomial_inline(&desc, &ic);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: init polynomial 2d");

		curve = qaws_curve_inline_get_curve(&ic);
		s = qaws_curve_evaluate_2d(curve, 0.0f,
			QAWS_EVAL_FLAG_POSITION, &result);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: poly eval t=0");
		TEST_ASSERT((float)fabs(result.position.x - 1.0f) < 1e-6f &&
			(float)fabs(result.position.y - 0.0f) < 1e-6f,
			"inline: poly start point");

		s = qaws_curve_evaluate_2d(curve, 1.0f,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &result);
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: poly eval t=1");
		/* P(1) = (1+2+0, 0+3-1) = (3, 2) */
		TEST_ASSERT((float)fabs(result.position.x - 3.0f) < 1e-6f &&
			(float)fabs(result.position.y - 2.0f) < 1e-6f,
			"inline: poly end point");
		/* P'(t) = (2, 3) + (0, -2)*t => P'(1) = (2, 1) */
		TEST_ASSERT((float)fabs(result.d1.x - 2.0f) < 1e-6f &&
			(float)fabs(result.d1.y - 1.0f) < 1e-6f,
			"inline: poly d1 at t=1");
	}

	/* Verify inline results match heap-allocated */
	{
		qaws_curve_inline ic;
		qaws_curve* heap_curve = NULL;
		qaws_curve* inline_curve;
		qaws_eval_result_3d r_heap, r_inline;
		qaws_vec3 pts[] = {{0,0,0}, {2,4,1}, {5,1,3}};
		qaws_bezier_desc desc;

		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 2;
		desc.control_points = pts;
		desc.control_point_count = 3;

		qaws_curve_create_bezier(&desc, &heap_curve);
		qaws_curve_init_bezier_inline(&desc, &ic);
		inline_curve = qaws_curve_inline_get_curve(&ic);

		qaws_curve_evaluate_3d(heap_curve, 0.3f,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
			&r_heap);
		qaws_curve_evaluate_3d(inline_curve, 0.3f,
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2,
			&r_inline);

		TEST_ASSERT(
			(float)fabs(r_heap.position.x - r_inline.position.x) < 1e-6f &&
			(float)fabs(r_heap.position.y - r_inline.position.y) < 1e-6f &&
			(float)fabs(r_heap.position.z - r_inline.position.z) < 1e-6f,
			"inline: matches heap position");
		TEST_ASSERT(
			(float)fabs(r_heap.d1.x - r_inline.d1.x) < 1e-6f &&
			(float)fabs(r_heap.d1.y - r_inline.d1.y) < 1e-6f &&
			(float)fabs(r_heap.d1.z - r_inline.d1.z) < 1e-6f,
			"inline: matches heap d1");
		TEST_ASSERT(
			(float)fabs(r_heap.d2.x - r_inline.d2.x) < 1e-6f &&
			(float)fabs(r_heap.d2.y - r_inline.d2.y) < 1e-6f &&
			(float)fabs(r_heap.d2.z - r_inline.d2.z) < 1e-6f,
			"inline: matches heap d2");

		qaws_curve_destroy(heap_curve);
		/* inline_curve: no destroy needed */
	}

	/* Error handling */
	{
		qaws_curve_inline ic;
		qaws_bezier_desc desc;
		desc.dimension = QAWS_DIMENSION_3D;
		desc.degree = 8; /* too high */
		desc.control_points = NULL;
		desc.control_point_count = 9;

		s = qaws_curve_init_bezier_inline(&desc, &ic);
		TEST_ASSERT(s != QAWS_STATUS_OK, "inline: reject degree 8");

		s = qaws_curve_init_bezier_inline(NULL, &ic);
		TEST_ASSERT(s == QAWS_STATUS_INVALID_ARGUMENT, "inline: null desc");

		s = qaws_curve_inline_get_curve(NULL) == NULL ? QAWS_STATUS_OK : QAWS_STATUS_INTERNAL_ERROR;
		TEST_ASSERT(s == QAWS_STATUS_OK, "inline: get_curve null");
	}
}

int test_23_performance_main(void)
{
	g_pass = 0;
	g_fail = 0;

	svg_ensure_output_dir();
	test_svg_batch_eval();
	test_batch_eval();
	test_custom_allocator();
	test_inline_curves();

	printf("23_performance: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
