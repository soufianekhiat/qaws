#include "test_common.h"

static void test_family_conversion(void)
{
	printf("test_family_conversion\n");

	/* Hermite to Bezier: endpoints should match */
	{
		qaws_vec2 pts[] = { {0,0}, {2,0} };
		qaws_vec2 ders[] = { {0,2}, {0,-2} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_2d r_h, r_b;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hermite, &hdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_hermite_to_bezier(&bezier, hermite, 0);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bezier) == QAWS_CURVE_KIND_BEZIER,
			"hermite->bezier kind");
		TEST_ASSERT(qaws_curve_get_degree(bezier) == 3, "hermite->bezier degree");

		/* Compare at start, mid, end */
		qaws_curve_evaluate_2d(&r_h, hermite,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b start x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b start y");

		qaws_curve_evaluate_2d(&r_h, hermite,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b mid x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b mid y");

		qaws_curve_evaluate_2d(&r_h, hermite,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b end x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b end y");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(hermite);
	}

	/* Catmull-Rom to Bezier */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {3,1}, {4,0} };
		qaws_catmull_rom_desc cdesc;
		qaws_curve *cr = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_2d r_cr, r_b;
		qaws_status status;

		cdesc.dimension = QAWS_DIMENSION_2D;
		cdesc.control_points = pts;
		cdesc.control_point_count = 4;
		cdesc.parameterization = QAWS_PARAMETERIZATION_UNIFORM;
		cdesc.closed = 0;

		status = qaws_curve_create_catmull_rom(&cr, &cdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_catmull_rom_to_bezier(&bezier, cr, 0);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bezier) == QAWS_CURVE_KIND_BEZIER,
			"cr->bezier kind");

		/* Compare at mid-span */
		{
			qaws_range cr_range = qaws_curve_get_parameter_range(cr);
			qaws_scalar mid = (cr_range.min_value + cr_range.max_value) / (qaws_scalar)2.0;
			qaws_curve_evaluate_2d(&r_cr, cr,
				mid, QAWS_EVAL_FLAG_POSITION);
		}
		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_cr.position.x, r_b.position.x), "cr2b mid x");
		TEST_ASSERT(approx_eq(r_cr.position.y, r_b.position.y), "cr2b mid y");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(cr);
	}

	/* Bezier to BSpline: should evaluate identically */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *bspline = NULL;
		qaws_eval_result_2d r_b, r_s;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&bezier, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_bezier_to_bspline(&bspline, bezier);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bspline) == QAWS_CURVE_KIND_BSPLINE,
			"bezier->bspline kind");

		/* Compare at multiple points */
		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)0.25, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_s, bspline,
			(qaws_scalar)0.25, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_b.position.x, r_s.position.x), "b2s 0.25 x");
		TEST_ASSERT(approx_eq(r_b.position.y, r_s.position.y), "b2s 0.25 y");

		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)0.75, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_s, bspline,
			(qaws_scalar)0.75, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_b.position.x, r_s.position.x), "b2s 0.75 x");
		TEST_ASSERT(approx_eq(r_b.position.y, r_s.position.y), "b2s 0.75 y");

		qaws_curve_destroy(bspline);
		qaws_curve_destroy(bezier);
	}

	/* BSpline to NURBS: should evaluate identically */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0}, {3,1} };
		qaws_bspline_desc sdesc;
		qaws_curve *bspline = NULL;
		qaws_curve *nurbs = NULL;
		qaws_eval_result_2d r_s, r_n;
		qaws_status status;

		memset(&sdesc, 0, sizeof(sdesc));
		sdesc.dimension = QAWS_DIMENSION_2D;
		sdesc.degree = 2;
		sdesc.control_points = pts;
		sdesc.control_point_count = 4;
		sdesc.is_uniform = 1;

		status = qaws_curve_create_bspline(&bspline, &sdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_bspline_to_nurbs(&nurbs, bspline);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(nurbs) == QAWS_CURVE_KIND_NURBS,
			"bspline->nurbs kind");
		TEST_ASSERT(qaws_curve_is_rational(nurbs), "nurbs is rational");

		/* Compare at midpoint */
		{
			qaws_range sr = qaws_curve_get_parameter_range(bspline);
			qaws_scalar mid = (sr.min_value + sr.max_value) / (qaws_scalar)2.0;
			qaws_curve_evaluate_2d(&r_s, bspline,
				mid, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&r_n, nurbs,
				mid, QAWS_EVAL_FLAG_POSITION);
		}
		TEST_ASSERT(approx_eq(r_s.position.x, r_n.position.x), "s2n mid x");
		TEST_ASSERT(approx_eq(r_s.position.y, r_n.position.y), "s2n mid y");

		qaws_curve_destroy(nurbs);
		qaws_curve_destroy(bspline);
	}

	/* Hermite-to-Bezier derivative check */
	{
		qaws_vec2 pts[] = { {0,0}, {3,0} };
		qaws_vec2 ders[] = { {1,3}, {1,-3} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_2d r_h, r_b;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hermite, &hdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_hermite_to_bezier(&bezier, hermite, 0);
		TEST_ASSERT_STATUS(status);

		/* D1 at start should match */
		qaws_curve_evaluate_2d(&r_h, hermite,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_D1);
		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_D1);
		TEST_ASSERT(approx_eq(r_h.d1.x, r_b.d1.x), "h2b d1 start x");
		TEST_ASSERT(approx_eq(r_h.d1.y, r_b.d1.y), "h2b d1 start y");

		/* D1 at end should match */
		qaws_curve_evaluate_2d(&r_h, hermite,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_D1);
		qaws_curve_evaluate_2d(&r_b, bezier,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_D1);
		TEST_ASSERT(approx_eq(r_h.d1.x, r_b.d1.x), "h2b d1 end x");
		TEST_ASSERT(approx_eq(r_h.d1.y, r_b.d1.y), "h2b d1 end y");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(hermite);
	}

	/* 3D Hermite-to-Bezier */
	{
		qaws_vec3 pts[] = { {0,0,0}, {2,0,1} };
		qaws_vec3 ders[] = { {0,2,1}, {0,-2,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *bezier = NULL;
		qaws_eval_result_3d r_h, r_b;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_3D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;

		status = qaws_curve_create_hermite(&hermite, &hdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_hermite_to_bezier(&bezier, hermite, 0);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_dimension(bezier) == QAWS_DIMENSION_3D,
			"h2b 3d dimension");

		/* Start position */
		qaws_curve_evaluate_3d(&r_h, hermite,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_3d(&r_b, bezier,
			(qaws_scalar)0.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b 3d start x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b 3d start y");
		TEST_ASSERT(approx_eq(r_h.position.z, r_b.position.z), "h2b 3d start z");

		/* Mid position */
		qaws_curve_evaluate_3d(&r_h, hermite,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_3d(&r_b, bezier,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b 3d mid x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b 3d mid y");
		TEST_ASSERT(approx_eq(r_h.position.z, r_b.position.z), "h2b 3d mid z");

		/* End position */
		qaws_curve_evaluate_3d(&r_h, hermite,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_3d(&r_b, bezier,
			(qaws_scalar)1.0, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_h.position.x, r_b.position.x), "h2b 3d end x");
		TEST_ASSERT(approx_eq(r_h.position.y, r_b.position.y), "h2b 3d end y");
		TEST_ASSERT(approx_eq(r_h.position.z, r_b.position.z), "h2b 3d end z");

		qaws_curve_destroy(bezier);
		qaws_curve_destroy(hermite);
	}

	/* Bezier-to-BSpline-to-NURBS chain */
	{
		qaws_vec2 pts[] = { {0,0}, {1,3}, {2,1}, {3,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *bspline = NULL;
		qaws_curve *nurbs = NULL;
		qaws_eval_result_2d r_bez, r_nurbs;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&bezier, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_convert_bezier_to_bspline(&bspline, bezier);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(bspline) == QAWS_CURVE_KIND_BSPLINE,
			"chain b2s kind");

		status = qaws_curve_convert_bspline_to_nurbs(&nurbs, bspline);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_kind(nurbs) == QAWS_CURVE_KIND_NURBS,
			"chain s2n kind");

		/* Compare Bezier and NURBS at several points */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.25) {
			qaws_scalar param = t;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(&r_bez, bezier,
				param, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&r_nurbs, nurbs,
				param, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_bez.position.x, r_nurbs.position.x),
				"chain b->s->n x");
			TEST_ASSERT(approx_eq(r_bez.position.y, r_nurbs.position.y),
				"chain b->s->n y");
		}

		qaws_curve_destroy(nurbs);
		qaws_curve_destroy(bspline);
		qaws_curve_destroy(bezier);
	}

	/* Span index out of range */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1}, {2,0} };
		qaws_vec2 ders[] = { {1,0}, {0,1}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *out = NULL;
		qaws_status status;
		unsigned int span_count;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 3;
		hdesc.derivative_count = 3;

		status = qaws_curve_create_hermite(&hermite, &hdesc);
		TEST_ASSERT_STATUS(status);

		span_count = qaws_curve_get_span_count(hermite);

		/* Request span_index == span_count (one past end) */
		status = qaws_curve_convert_hermite_to_bezier(&out, hermite, span_count);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"span index out of range");

		/* Request a large index */
		status = qaws_curve_convert_hermite_to_bezier(&out, hermite, 999);
		TEST_ASSERT(
			status == QAWS_STATUS_OUT_OF_RANGE ||
			status == QAWS_STATUS_INVALID_ARGUMENT,
			"span index 999 out of range");

		qaws_curve_destroy(hermite);
	}

	/* Error: wrong family */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1} };
		qaws_bezier_desc bdesc;
		qaws_curve *bezier = NULL;
		qaws_curve *out = NULL;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;
		qaws_curve_create_bezier(&bezier, &bdesc);

		status = qaws_curve_convert_hermite_to_bezier(&out, bezier, 0);
		TEST_ASSERT(status == QAWS_STATUS_UNSUPPORTED_OPERATION,
			"wrong family conversion");

		qaws_curve_destroy(bezier);
	}
}

static void test_degree_elevation(void)
{
	printf("test_degree_elevation\n");

	/* Elevate quadratic Bezier to cubic: should evaluate identically */
	{
		qaws_vec2 pts[] = { {0,0}, {1,2}, {2,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *original = NULL;
		qaws_curve *elevated = NULL;
		qaws_eval_result_2d r_orig, r_elev;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&original, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(&elevated, original);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(elevated) == 3, "elevated degree");

		/* Compare at multiple parameter values */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.25) {
			qaws_scalar param = t;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_2d(&r_orig, original,
				param, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&r_elev, elevated,
				param, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_elev.position.x),
				"elevated matches original x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_elev.position.y),
				"elevated matches original y");
		}

		qaws_curve_destroy(elevated);
		qaws_curve_destroy(original);
	}

	/* Double elevation: linear -> quadratic -> cubic */
	{
		qaws_vec2 pts[] = { {0,0}, {4,4} };
		qaws_bezier_desc bdesc;
		qaws_curve *linear = NULL;
		qaws_curve *quad = NULL;
		qaws_curve *cubic = NULL;
		qaws_eval_result_2d r_l, r_c;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 1;
		bdesc.control_points = pts;
		bdesc.control_point_count = 2;

		status = qaws_curve_create_bezier(&linear, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(&quad, linear);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(quad) == 2, "quad degree");

		status = qaws_curve_elevate_degree(&cubic, quad);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(cubic) == 3, "cubic degree");

		/* Still a straight line: mid = (2,2) */
		qaws_curve_evaluate_2d(&r_l, linear,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		qaws_curve_evaluate_2d(&r_c, cubic,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION);
		TEST_ASSERT(approx_eq(r_l.position.x, r_c.position.x),
			"double elevation mid x");
		TEST_ASSERT(approx_eq(r_l.position.y, r_c.position.y),
			"double elevation mid y");

		qaws_curve_destroy(cubic);
		qaws_curve_destroy(quad);
		qaws_curve_destroy(linear);
	}

	/* 3D elevation: quadratic to cubic */
	{
		qaws_vec3 pts[] = { {0,0,0}, {1,2,1}, {2,0,2} };
		qaws_bezier_desc bdesc;
		qaws_curve *original = NULL;
		qaws_curve *elevated = NULL;
		qaws_eval_result_3d r_orig, r_elev;
		qaws_status status;
		qaws_scalar t;

		bdesc.dimension = QAWS_DIMENSION_3D;
		bdesc.degree = 2;
		bdesc.control_points = pts;
		bdesc.control_point_count = 3;

		status = qaws_curve_create_bezier(&original, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(&elevated, original);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(elevated) == 3, "3d elevated degree");
		TEST_ASSERT(qaws_curve_get_dimension(elevated) == QAWS_DIMENSION_3D,
			"3d elevated dimension");

		/* Compare at multiple parameter values */
		for (t = (qaws_scalar)0.0; t <= (qaws_scalar)1.001;
			t += (qaws_scalar)0.2) {
			qaws_scalar param = t;
			if (param > (qaws_scalar)1.0) param = (qaws_scalar)1.0;
			qaws_curve_evaluate_3d(&r_orig, original,
				param, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_3d(&r_elev, elevated,
				param, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r_orig.position.x, r_elev.position.x),
				"3d elevated x");
			TEST_ASSERT(approx_eq(r_orig.position.y, r_elev.position.y),
				"3d elevated y");
			TEST_ASSERT(approx_eq(r_orig.position.z, r_elev.position.z),
				"3d elevated z");
		}

		qaws_curve_destroy(elevated);
		qaws_curve_destroy(original);
	}

	/* Derivative preservation: elevate cubic, verify D1 at t=0.5 */
	{
		qaws_vec2 pts[] = { {0,0}, {1,3}, {3,3}, {4,0} };
		qaws_bezier_desc bdesc;
		qaws_curve *original = NULL;
		qaws_curve *elevated = NULL;
		qaws_eval_result_2d r_orig, r_elev;
		qaws_status status;

		bdesc.dimension = QAWS_DIMENSION_2D;
		bdesc.degree = 3;
		bdesc.control_points = pts;
		bdesc.control_point_count = 4;

		status = qaws_curve_create_bezier(&original, &bdesc);
		TEST_ASSERT_STATUS(status);

		status = qaws_curve_elevate_degree(&elevated, original);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(qaws_curve_get_degree(elevated) == 4, "cubic elevated to 4");

		/* D1 at midpoint should match */
		qaws_curve_evaluate_2d(&r_orig, original,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		qaws_curve_evaluate_2d(&r_elev, elevated,
			(qaws_scalar)0.5, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1);
		TEST_ASSERT(approx_eq(r_orig.position.x, r_elev.position.x),
			"elevated d1 pos x");
		TEST_ASSERT(approx_eq(r_orig.position.y, r_elev.position.y),
			"elevated d1 pos y");
		TEST_ASSERT(approx_eq(r_orig.d1.x, r_elev.d1.x),
			"elevated d1 match x");
		TEST_ASSERT(approx_eq(r_orig.d1.y, r_elev.d1.y),
			"elevated d1 match y");

		qaws_curve_destroy(elevated);
		qaws_curve_destroy(original);
	}

	/* Error: non-Bezier */
	{
		qaws_vec2 pts[] = { {0,0}, {1,1} };
		qaws_vec2 ders[] = { {1,0}, {1,0} };
		qaws_hermite_desc hdesc;
		qaws_curve *hermite = NULL;
		qaws_curve *out = NULL;
		qaws_status status;

		hdesc.dimension = QAWS_DIMENSION_2D;
		hdesc.degree = 3;
		hdesc.points = pts;
		hdesc.derivatives = ders;
		hdesc.point_count = 2;
		hdesc.derivative_count = 2;
		qaws_curve_create_hermite(&hermite, &hdesc);

		status = qaws_curve_elevate_degree(&out, hermite);
		TEST_ASSERT(status == QAWS_STATUS_UNSUPPORTED_OPERATION,
			"elevate non-bezier");

		qaws_curve_destroy(hermite);
	}
}

static void test_degree_reduction(void)
{
	printf("test_degree_reduction\n");

	/* Reduce a degree-elevated Bezier: elevate then reduce should give approx same curve */
	{
		qaws_vec2 cp[] = { {0,0}, {1,3}, {3,3}, {4,0} };
		qaws_bezier_desc desc;
		qaws_curve *src = NULL, *elevated = NULL, *reduced = NULL;
		qaws_eval_result_2d r_src, r_red;
		qaws_range rng;
		int close = 1;
		unsigned int i;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.control_points = cp;
		desc.control_point_count = 4;

		TEST_ASSERT_STATUS(qaws_curve_create_bezier(&src, &desc));

		/* Elevate degree 3 -> 4, then reduce 4 -> 3 */
		TEST_ASSERT_STATUS(qaws_curve_elevate_degree(&elevated, src));
		TEST_ASSERT(qaws_curve_get_degree(elevated) == 4, "elevated degree 4");

		TEST_ASSERT_STATUS(qaws_curve_reduce_degree(&reduced, elevated));
		TEST_ASSERT(reduced != NULL, "reduced created");
		TEST_ASSERT(qaws_curve_get_degree(reduced) == 3, "reduced degree 3");

		/* Compare at multiple parameter values */
		rng = qaws_curve_get_parameter_range(src);
		for (i = 0; i <= 20; ++i) {
			qaws_scalar t = rng.min_value + (rng.max_value - rng.min_value) * (qaws_scalar)i / (qaws_scalar)20;
			qaws_scalar t_red;
			qaws_range rng_red = qaws_curve_get_parameter_range(reduced);
			t_red = rng_red.min_value + (rng_red.max_value - rng_red.min_value) * (qaws_scalar)i / (qaws_scalar)20;

			qaws_curve_evaluate_2d(&r_src, src, t, QAWS_EVAL_FLAG_POSITION);
			qaws_curve_evaluate_2d(&r_red, reduced, t_red, QAWS_EVAL_FLAG_POSITION);

			if (!approx_eq(r_src.position.x, r_red.position.x) ||
				!approx_eq(r_src.position.y, r_red.position.y))
				close = 0;
		}
		TEST_ASSERT(close, "elevate-then-reduce recovers original");

		qaws_curve_destroy(reduced);
		qaws_curve_destroy(elevated);
		qaws_curve_destroy(src);
	}

	/* Reduce quadratic to linear */
	{
		qaws_vec2 cp[] = { {0,0}, {1,2}, {2,0} };
		qaws_bezier_desc desc;
		qaws_curve *src = NULL, *reduced = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 2;
		desc.control_points = cp;
		desc.control_point_count = 3;

		TEST_ASSERT_STATUS(qaws_curve_create_bezier(&src, &desc));
		TEST_ASSERT_STATUS(qaws_curve_reduce_degree(&reduced, src));
		TEST_ASSERT(qaws_curve_get_degree(reduced) == 1, "reduced to linear");

		/* Endpoints should match */
		{
			qaws_eval_result_2d r;
			qaws_range rng = qaws_curve_get_parameter_range(reduced);
			qaws_curve_evaluate_2d(&r, reduced, rng.min_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)0) &&
				approx_eq(r.position.y, (qaws_scalar)0), "reduced start");
			qaws_curve_evaluate_2d(&r, reduced, rng.max_value, QAWS_EVAL_FLAG_POSITION);
			TEST_ASSERT(approx_eq(r.position.x, (qaws_scalar)2) &&
				approx_eq(r.position.y, (qaws_scalar)0), "reduced end");
		}

		qaws_curve_destroy(reduced);
		qaws_curve_destroy(src);
	}

	/* Error: degree 1 cannot be reduced */
	{
		qaws_vec2 cp[] = { {0,0}, {1,1} };
		qaws_bezier_desc desc;
		qaws_curve *src = NULL, *reduced = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 1;
		desc.control_points = cp;
		desc.control_point_count = 2;

		TEST_ASSERT_STATUS(qaws_curve_create_bezier(&src, &desc));
		TEST_ASSERT(qaws_curve_reduce_degree(&reduced, src) == QAWS_STATUS_INVALID_DEGREE,
			"degree 1 reduce fails");
		TEST_ASSERT(reduced == NULL, "degree 1 reduce null");

		qaws_curve_destroy(src);
	}

	/* Error: non-Bezier curve */
	{
		qaws_scalar pts[] = { 0,0,  1,1,  2,0 };
		qaws_scalar derivs[] = { 1,0,  1,0,  1,0 };
		qaws_hermite_desc desc;
		qaws_curve *src = NULL, *reduced = NULL;

		memset(&desc, 0, sizeof(desc));
		desc.dimension = QAWS_DIMENSION_2D;
		desc.degree = 3;
		desc.points = pts;
		desc.derivatives = derivs;
		desc.point_count = 3;
		desc.derivative_count = 3;

		TEST_ASSERT_STATUS(qaws_curve_create_hermite(&src, &desc));
		TEST_ASSERT(qaws_curve_reduce_degree(&reduced, src) == QAWS_STATUS_UNSUPPORTED_OPERATION,
			"hermite reduce unsupported");

		qaws_curve_destroy(src);
	}
}

int test_19_conversion_main(void)
{
	g_pass = 0;
	g_fail = 0;

	test_family_conversion();
	test_degree_elevation();
	test_degree_reduction();

	printf("19_conversion: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail > 0 ? 1 : 0;
}
