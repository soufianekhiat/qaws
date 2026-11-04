/*
 * Test 157: guaranteed ranges of curves over parameter intervals (qaws_curve_enclose)
 *
 *   - a NURBS circle (rational, 2D), a cubic B-spline (several spans, 3D), a Hermite spline (2D), a uniform Catmull-Rom
 *     spline (2D), a cubic Bezier (3D): over intervals of a half, a tenth and a hundredth of the range, every sampled C,
 *     C', C'' inside its range;
 *   - the ranges shrink with the interval (a point: the position to rounding);
 *   - the refusals: an interval past the range, an empty one, an order past 2, an unsupported curve.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

static qaws_curve* ce_circle(void)
{
	static double const cx[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, cy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	static qaws_scalar const kn[12] = { 0, 0, 0, 0.25f, 0.25f, 0.5f, 0.5f, 0.75f, 0.75f, 1, 1, 1 };
	qaws_vec2 cp[9];
	qaws_scalar w[9];
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 9; i++)
	{
		cp[i].x = (qaws_scalar)cx[i], cp[i].y = (qaws_scalar)cy[i];
		w[i] = (qaws_scalar)(i & 1 ? sqrt(0.5) : 1);
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D, d.degree = 2;
	d.control_points = cp, d.control_point_count = 9;
	d.weights = w, d.weight_count = 9;
	d.knots = kn, d.knot_count = 12;
	qaws_curve_create_nurbs(&c, &d);
	return c;
}

static qaws_curve* ce_bspline(void)
{
	static qaws_scalar const kn[10] = { 0, 0, 0, 0, 0.3f, 0.6f, 1, 1, 1, 1 };
	qaws_vec3 cp[6];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 6; i++)
		cp[i].x = (qaws_scalar)i, cp[i].y = (qaws_scalar)sin(1.3 * i), cp[i].z = (qaws_scalar)cos(0.7 * i);
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D, d.degree = 3;
	d.control_points = cp, d.control_point_count = 6;
	d.knots = kn, d.knot_count = 10;
	qaws_curve_create_bspline(&c, &d);
	return c;
}

static qaws_curve* ce_hermite(void)
{
	qaws_vec2 pts[4], der[4];
	qaws_hermite_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 4; i++)
	{
		pts[i].x = (qaws_scalar)i, pts[i].y = (qaws_scalar)(0.3 * i * i);
		der[i].x = (qaws_scalar)(1 + 0.1 * i), der[i].y = (qaws_scalar)(0.6 * i - 0.7);
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D, d.degree = 3;
	d.points = pts, d.point_count = 4;
	d.derivatives = der, d.derivative_count = 4;
	qaws_curve_create_hermite(&c, &d);
	return c;
}

static qaws_curve* ce_catmull_rom(void)
{
	qaws_vec2 pts[6];
	qaws_catmull_rom_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 6; i++)
		pts[i].x = (qaws_scalar)i, pts[i].y = (qaws_scalar)(i & 1 ? 1 : -0.5);
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.control_points = pts, d.control_point_count = 6;
	d.parameterization = QAWS_PARAMETERIZATION_UNIFORM;
	qaws_curve_create_catmull_rom(&c, &d);
	return c;
}

static qaws_curve* ce_bezier(void)
{
	qaws_vec3 cp[4];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 4; i++)
		cp[i].x = (qaws_scalar)i, cp[i].y = (qaws_scalar)((i * 5 % 4) * 0.5), cp[i].z = (qaws_scalar)(i == 2 ? 1 : 0);
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D, d.degree = 3;
	d.control_points = cp, d.control_point_count = 4;
	qaws_curve_create_bezier(&c, &d);
	return c;
}

/* every sample of the interval inside the ranges (a slack for the evaluation's own rounding); the ranges' widest width */
static int ce_check(qaws_curve const* c, double t0, double t1, double* widest)
{
	double lo[9], hi[9];
	unsigned int const dim = qaws_curve_get_dimension(c) == QAWS_DIMENSION_3D ? 3u : 2u;
	unsigned int i, q, k;
	int ok = qaws_curve_enclose(lo, hi, c, t0, t1, 2) == QAWS_STATUS_OK;
	double const slack = sizeof(qaws_scalar) == sizeof(double) ? 1e-9 : 2e-3;
	*widest = 0;
	for (q = 0; ok && q < 3 * dim; q++)
		*widest = fmax(*widest, hi[q] - lo[q]);
	for (i = 0; ok && i <= 16; i++)
	{
		double const t = t0 + (t1 - t0) * i / 16;
		double v[9];
		if (dim == 3)
		{
			qaws_eval_result_3d r;
			qaws_curve_evaluate_3d(&r, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
			v[0] = r.position.x, v[1] = r.position.y, v[2] = r.position.z;
			v[3] = r.d1.x, v[4] = r.d1.y, v[5] = r.d1.z;
			v[6] = r.d2.x, v[7] = r.d2.y, v[8] = r.d2.z;
		}
		else
		{
			qaws_eval_result_2d r;
			qaws_curve_evaluate_2d(&r, c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
			v[0] = r.position.x, v[1] = r.position.y;
			v[2] = r.d1.x, v[3] = r.d1.y;
			v[4] = r.d2.x, v[5] = r.d2.y;
		}
		for (q = 0; q < 3; q++)
			for (k = 0; k < dim; k++)
			{
				double const x = v[q * dim + k], scale = 1 + fabs(x);
				if (x < lo[q * dim + k] - slack * scale || x > hi[q * dim + k] + slack * scale)
					ok = 0;
			}
	}
	return ok;
}

int test_157_curve_enclose_main(void)
{
	qaws_curve* curves[5];
	char const* names[5] = { "a NURBS circle", "a cubic B-spline", "a Hermite spline", "a uniform Catmull-Rom spline", "a cubic Bezier" };
	char msg[300];
	unsigned int k;
	g_pass = 0;
	g_fail = 0;
	printf("Test 157: guaranteed ranges of curves over parameter intervals (qaws_curve_enclose)\n");
	curves[0] = ce_circle();
	curves[1] = ce_bspline();
	curves[2] = ce_hermite();
	curves[3] = ce_catmull_rom();
	curves[4] = ce_bezier();
	TEST_ASSERT(curves[0] && curves[1] && curves[2] && curves[3] && curves[4], "a circle, a B-spline, a Hermite, a Catmull-Rom, a Bezier");
	for (k = 0; k < 5; k++)
	{
		static double const size[3] = { 0.5, 0.1, 0.01 };
		double widths[3] = { 0, 0, 0 }, point_lo[3], point_hi[3], t;
		qaws_range r;
		unsigned int z;
		int all = 1;
		if (!curves[k])
			continue;
		r = qaws_curve_get_parameter_range(curves[k]);
		t = r.min_value + 0.43 * (r.max_value - r.min_value);   /* (inside a span at every size: C'' jumps at knots) */
		for (z = 0; z < 3; z++)
		{
			double const w = size[z] * (r.max_value - r.min_value);
			all = all && ce_check(curves[k], t - w / 2, t + w / 2, &widths[z]);
		}
		qaws_curve_enclose(point_lo, point_hi, curves[k], t, t, 0);
		sprintf(msg, "%s: every sample inside; widest range %.3g, %.3g, %.3g over a half, a tenth, a hundredth of the range; a point's %.1e", names[k],
			widths[0], widths[1], widths[2], fmax(point_hi[0] - point_lo[0], point_hi[1] - point_lo[1]));
		printf("    %s\n", msg);
		TEST_ASSERT(all && widths[1] < widths[0] && widths[2] < widths[1] && point_hi[0] - point_lo[0] < 1e-12, msg);
	}
	/* refusals */
	{
		double lo[9], hi[9];
		qaws_curve* arc = NULL;
		qaws_arc_segment seg;
		qaws_arc_desc ad;
		TEST_ASSERT(qaws_curve_enclose(lo, hi, curves[0], -0.1, 0.2, 0) == QAWS_STATUS_INVALID_ARGUMENT, "an interval past the range refused");
		TEST_ASSERT(qaws_curve_enclose(lo, hi, curves[0], 0.3, 0.2, 0) == QAWS_STATUS_INVALID_ARGUMENT, "an empty interval refused");
		TEST_ASSERT(qaws_curve_enclose(lo, hi, curves[0], 0, 1, 3) == QAWS_STATUS_INVALID_ARGUMENT, "an order past 2 refused");
		memset(&seg, 0, sizeof seg);
		seg.radius = 1, seg.angle_end = 1;
		memset(&ad, 0, sizeof ad);
		ad.dimension = QAWS_DIMENSION_2D, ad.segments = &seg, ad.segment_count = 1;
		qaws_curve_create_arc(&arc, &ad);
		{
			qaws_range const ar = arc ? qaws_curve_get_parameter_range(arc) : qaws_curve_get_parameter_range(curves[0]);
			TEST_ASSERT(!arc || qaws_curve_enclose(lo, hi, arc, ar.min_value, ar.max_value, 0) == QAWS_STATUS_UNSUPPORTED_OPERATION, "an arc unsupported");
		}
		qaws_curve_destroy(arc);
	}
	for (k = 0; k < 5; k++)
		qaws_curve_destroy(curves[k]);
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
