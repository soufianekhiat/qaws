/*
 * Test 150: adaptive curve tessellation (qaws_curve_tessellate_2d / 3d)
 *
 *   - a rational circle at a chord tolerance: every segment's true distance
 *     to the curve (64 points each) against the tolerance;
 *   - an S whose middle lies on its chord (a midpoint test stops there);
 *   - an angle tolerance: the tangent's turning over every segment;
 *   - a kink at a double knot is a vertex; a straight span one segment;
 *   - counts first, then the same vertices every call; a short buffer;
 *   - refusals.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define TS_CAP 4096

static double ts_seg_error_2d(qaws_curve const* c, qaws_scalar const* t, qaws_vec2 const* p, unsigned int n)
{
	double worst = 0;
	unsigned int i, k;
	for (i = 0; i + 1 < n; i++)
		for (k = 1; k < 64; k++)
		{
			qaws_eval_result_2d e;
			double ax = p[i + 1].x - p[i].x, ay = p[i + 1].y - p[i].y, l = ax * ax + ay * ay, f, dx, dy;
			qaws_curve_evaluate_2d(c, t[i] + (t[i + 1] - t[i]) * (qaws_scalar)k / 64, QAWS_EVAL_FLAG_POSITION, &e);
			f = l > 0 ? ((e.position.x - p[i].x) * ax + (e.position.y - p[i].y) * ay) / l : 0;
			f = f < 0 ? 0 : f > 1 ? 1 : f;
			dx = e.position.x - p[i].x - f * ax;
			dy = e.position.y - p[i].y - f * ay;
			worst = fmax(worst, sqrt(dx * dx + dy * dy));
		}
	return worst;
}

static qaws_curve* ts_circle(void)
{
	static qaws_vec2 p[9] = { { 2, 0 }, { 2, 2 }, { 0, 2 }, { -2, 2 }, { -2, 0 }, { -2, -2 }, { 0, -2 }, { 2, -2 }, { 2, 0 } };
	static qaws_scalar knots[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4 };
	qaws_scalar w[9];
	qaws_nurbs_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 9; i++)
		w[i] = (qaws_scalar)(i % 2 ? 0.5 * sqrt(2.0) : 1);
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = p; d.control_point_count = 9;
	d.knots = knots; d.knot_count = 12; d.weights = w; d.weight_count = 9;
	qaws_curve_create_nurbs(&d, &c);
	return c;
}

static void test_chord(void)
{
	static qaws_scalar t[TS_CAP], t2[TS_CAP];
	static qaws_vec2 p[TS_CAP], p2[TS_CAP];
	qaws_curve* c = ts_circle();
	qaws_curve_tessellation_desc d;
	unsigned int n = 0, n2 = 0, n3 = 0, i;
	double err;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.chord_tolerance = (qaws_scalar)1e-3;
	TEST_ASSERT(qaws_curve_tessellate_2d(c, &d, NULL, NULL, 0, &n) == QAWS_STATUS_OK && n > 8, "circle: the count first");
	TEST_ASSERT(qaws_curve_tessellate_2d(c, &d, t, p, TS_CAP, &n2) == QAWS_STATUS_OK && n2 == n, "circle: the same count filled");
	err = ts_seg_error_2d(c, t, p, n);
	snprintf(msg, sizeof msg, "circle of radius 2 at chord 1e-3: %u vertices, true chord error %.3e", n, err);
	printf("  %s\n", msg);
	TEST_ASSERT(err <= 1.05e-3, msg);
	/* the fewest equal segments meeting the tolerance: 2 acos(1 - 1e-3 / 2) each */
	TEST_ASSERT(n - 1 <= 4 * (unsigned int)ceil(2 * M_PI / (2 * acos(1 - 1e-3 / 2)) / 4 * 2), "no more than twice the fewest");
	TEST_ASSERT(t[0] == 0 && t[n - 1] == 4 && p[0].x == 2 && p[0].y == 0, "starts and ends at the range's ends");
	for (i = 1; i < n; i++)
		if (!(t[i] > t[i - 1]))
			break;
	TEST_ASSERT(i == n, "parameters increasing");
	qaws_curve_tessellate_2d(c, &d, t2, p2, TS_CAP, &n2);
	TEST_ASSERT(n2 == n && memcmp(t, t2, sizeof(qaws_scalar) * n) == 0 && memcmp(p, p2, sizeof(qaws_vec2) * n) == 0, "the same vertices every call");
	TEST_ASSERT(qaws_curve_tessellate_2d(c, &d, t2, NULL, 10, &n3) == QAWS_STATUS_BUFFER_TOO_SMALL && n3 == n
		&& memcmp(t, t2, sizeof(qaws_scalar) * 10) == 0, "short buffer: the count and the first vertices");
	qaws_curve_destroy(c);
}

static void test_s_and_angle(void)
{
	static qaws_scalar t[TS_CAP];
	static qaws_vec2 p[TS_CAP];
	/* the S: its middle (1.5, 0) is on its chord */
	qaws_vec2 s[4] = { { 0, 0 }, { 1, 1 }, { 2, -1 }, { 3, 0 } };
	qaws_bezier_desc bd;
	qaws_curve *b = NULL, *c = ts_circle();
	qaws_curve_tessellation_desc d;
	unsigned int n = 0, i;
	double err, worst = 0;
	char msg[200];
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 3; bd.control_points = s; bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &b);
	memset(&d, 0, sizeof d);
	d.chord_tolerance = (qaws_scalar)1e-2;
	qaws_curve_tessellate_2d(b, &d, t, p, TS_CAP, &n);
	err = ts_seg_error_2d(b, t, p, n);
	snprintf(msg, sizeof msg, "S curve at chord 1e-2: %u vertices, true chord error %.3e", n, err);
	printf("  %s\n", msg);
	TEST_ASSERT(n > 4 && err <= 1.05e-2, msg);

	memset(&d, 0, sizeof d);
	d.angle_tolerance = (qaws_scalar)0.1;
	qaws_curve_tessellate_2d(c, &d, t, p, TS_CAP, &n);
	for (i = 0; i + 1 < n; i++)
	{
		qaws_eval_result_2d e0, e1;
		qaws_curve_evaluate_2d(c, t[i] + (qaws_scalar)1e-9, QAWS_EVAL_FLAG_D1, &e0);
		qaws_curve_evaluate_2d(c, t[i + 1] - (qaws_scalar)1e-9, QAWS_EVAL_FLAG_D1, &e1);
		worst = fmax(worst, fabs(atan2(e0.d1.x * e1.d1.y - e0.d1.y * e1.d1.x, e0.d1.x * e1.d1.x + e0.d1.y * e1.d1.y)));
	}
	snprintf(msg, sizeof msg, "circle at angle 0.1: %u vertices, largest turning %.4f", n, worst);
	printf("  %s\n", msg);
	TEST_ASSERT(worst <= 0.1 && n - 1 >= (unsigned int)(2 * M_PI / 0.1), msg);
	qaws_curve_destroy(b);
	qaws_curve_destroy(c);
}

static void test_kink_and_line(void)
{
	static qaws_scalar t[TS_CAP];
	/* a cubic with a triple knot at 1: a corner there; its second span straight */
	qaws_vec3 p[7] = { { 0, 0, 0 }, { 1, 2, 0 }, { 2, 2, 1 }, { 3, 0, 0 }, { 4, 0, 0 }, { 5, 0, 0 }, { 6, 0, 0 } };
	qaws_scalar knots[11] = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2 };
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	qaws_curve_tessellation_desc td;
	unsigned int n = 0, i, has_kink = 0, after = 0;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = p; d.control_point_count = 7; d.knots = knots; d.knot_count = 11;
	TEST_ASSERT(qaws_curve_create_bspline(&d, &c) == QAWS_STATUS_OK, "kinked curve created");
	memset(&td, 0, sizeof td);
	td.chord_tolerance = (qaws_scalar)1e-3;
	td.angle_tolerance = (qaws_scalar)0.2;
	qaws_curve_tessellate_3d(c, &td, t, NULL, TS_CAP, &n);
	for (i = 0; i < n; i++)
	{
		has_kink |= t[i] == 1;
		after += t[i] > 1;
	}
	TEST_ASSERT(has_kink, "the kink at the triple knot is a vertex");
	TEST_ASSERT(after == 1, "the straight span is one segment");
	qaws_curve_destroy(c);
}

static void test_refusals(void)
{
	qaws_curve* c = ts_circle();
	qaws_curve_tessellation_desc d;
	unsigned int n;
	memset(&d, 0, sizeof d);
	TEST_ASSERT(qaws_curve_tessellate_2d(c, &d, NULL, NULL, 0, &n) == QAWS_STATUS_INVALID_ARGUMENT, "no tolerance refused");
	d.chord_tolerance = -1;
	TEST_ASSERT(qaws_curve_tessellate_2d(c, &d, NULL, NULL, 0, &n) == QAWS_STATUS_INVALID_ARGUMENT, "a negative tolerance refused");
	d.chord_tolerance = (qaws_scalar)0.1;
	TEST_ASSERT(qaws_curve_tessellate_3d(c, &d, NULL, NULL, 0, &n) == QAWS_STATUS_INVALID_DIMENSION, "a 2D curve as 3D refused");
	TEST_ASSERT(qaws_curve_tessellate_2d(c, &d, NULL, NULL, 0, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "no count refused");
	qaws_curve_destroy(c);
}

int test_150_curve_tessellation_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 150: adaptive curve tessellation\n");
	test_chord();
	test_s_and_angle();
	test_kink_and_line();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
