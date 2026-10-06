/*
 * Test 76: Exact composite curves and frozen Catmull-Rom preparations
 *
 *   - a composite of a Bezier, a 3-span NURBS and a segment: one exact
 *     curve with every span moved to [i, i + 1]; positions and derivatives
 *     agree with the runtime composite; the NURBS bounds 1/3, 2/3 are not
 *     dyadic and are reported as a parameter quantization (the geometry
 *     stays exact: the span joints are the same exact points)
 *   - chordal and centripetal Catmull-Rom: the runtime's preparation frozen
 *     exactly (QAWS_EXACT_FLAG_PREP_QUANTIZED), values equal to its own
 *     evaluator to rounding
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "exact/qaws_exact_curve.h"
#include <math.h>
#include <string.h>

static double ec_gap(qaws_curve* c, qaws_exact_curve* e, double t0, double t1, unsigned int dim)
{
	double worst = 0;
	unsigned int i, k, j;
	for (i = 0; i <= 200; i++)
	{
		double t = ldexp(nearbyint(ldexp(t0 + (t1 - t0) * i / 200.0, 12)), -12), out[12];
		qaws_eval_result_3d r3;
		qaws_eval_result_2d r2;
		double ref[12];
		unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2;
		qaws_exact_curve_evaluate(e, t, 2, out, NULL);
		if (dim == 3)
		{
			qaws_curve_evaluate_3d(c, (qaws_scalar)t, flags, &r3);
			ref[0] = r3.position.x; ref[1] = r3.position.y; ref[2] = r3.position.z;
			ref[3] = r3.d1.x; ref[4] = r3.d1.y; ref[5] = r3.d1.z;
			ref[6] = r3.d2.x; ref[7] = r3.d2.y; ref[8] = r3.d2.z;
		}
		else
		{
			qaws_curve_evaluate_2d(c, (qaws_scalar)t, flags, &r2);
			ref[0] = r2.position.x; ref[1] = r2.position.y;
			ref[2] = r2.d1.x; ref[3] = r2.d1.y;
			ref[4] = r2.d2.x; ref[5] = r2.d2.y;
		}
		for (k = 0; k <= 2; k++)
			for (j = 0; j < dim; j++)
			{
				double g = fabs(out[k * dim + j] - ref[k * dim + j]) / (1 + fabs(ref[k * dim + j]));
				if (g > worst) worst = g;
			}
	}
	return worst;
}

static void test_composite(void)
{
	qaws_scalar bz[8] = { 0, 0, 1, 2, 2, 2, 3, 0 };
	qaws_scalar nb[10] = { 3, 0, 4, -2, 5, 1, 6, -1, 7, 0 }, nw[5] = { 1, 3, 1, 2, 1 }, nk[8] = { 0, 0, 0, 1, 2, 3, 3, 3 };
	qaws_scalar sg[4] = { 7, 0, 8, 1 };
	qaws_bezier_desc bd;
	qaws_nurbs_desc nd;
	qaws_composite_desc cd;
	qaws_curve* segs[3];
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_report rep;
	qaws_status st;
	double gap;
	unsigned int i;
	int joints = 1;
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D;
	bd.degree = 3;
	bd.control_points = bz;
	bd.control_point_count = 4;
	qaws_curve_create_bezier(&bd, &segs[0]);
	memset(&nd, 0, sizeof(nd));
	nd.dimension = QAWS_DIMENSION_2D;
	nd.degree = 2;
	nd.control_points = nb;
	nd.control_point_count = 5;
	nd.weights = nw;
	nd.weight_count = 5;
	nd.knots = nk;
	nd.knot_count = 8;
	qaws_curve_create_nurbs(&nd, &segs[1]);
	bd.degree = 1;
	bd.control_points = sg;
	bd.control_point_count = 2;
	qaws_curve_create_bezier(&bd, &segs[2]);
	memset(&cd, 0, sizeof(cd));
	cd.dimension = QAWS_DIMENSION_2D;
	cd.segments = segs;
	cd.segment_count = 3;
	qaws_curve_create_composite(&cd, &c);
	st = qaws_exact_curve_prepare(NULL, c, &e, &rep);
	TEST_ASSERT(st == QAWS_STATUS_OK && qaws_exact_curve_span_count(e) == 5, "a composite of three curves: one exact curve of 1 + 3 + 1 spans");
	if (st != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(c);
		return;
	}
	TEST_ASSERT((rep.flags & QAWS_EXACT_FLAG_INPUT_QUANTIZED) && rep.parameter_quantization_error > 0 && rep.parameter_quantization_error < 1e-7,
		"the non-dyadic span bounds 1/3, 2/3 are reported as a parameter quantization");
	/* joints: each span's end point is exactly the next span's start point */
	for (i = 0; i + 1 < e->span_count; i++)
	{
		qaws_exact_span const* a = &e->spans[i];
		qaws_exact_span const* b = &e->spans[i + 1];
		qaws_exact_int l, r;
		unsigned int c2;
		for (c2 = 0; c2 < 2; c2++)
		{
			qaws_exact_int_mul(&l, &a->h[a->degree * 3 + c2], &b->h[2]);
			qaws_exact_int_mul(&r, &b->h[c2], &a->h[a->degree * 3 + 2]);
			joints &= qaws_exact_int_cmp(&l, &r) == 0;
		}
	}
	TEST_ASSERT(joints, "every joint is continuous exactly");
	gap = ec_gap(c, e, 0.0, 0.99, 2);   /* t = 1 is the joint: the next segment owns it */
	printf("    composite: %u spans, parameter quantization %.2e; vs runtime composite on the Bezier: %.2e", e->span_count,
		rep.parameter_quantization_error, gap);
	TEST_ASSERT(gap < (QAWS_SCALAR_IS_FLOAT ? 1e-4 : 1e-12), "on the dyadic segments: equal to the runtime composite (C, C', C'')");
	gap = ec_gap(c, e, 2.0, 3.0, 2);
	printf(", on the segment: %.2e\n", gap);
	TEST_ASSERT(gap < (QAWS_SCALAR_IS_FLOAT ? 1e-4 : 1e-12), "the last segment as well");
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
}

static void test_frozen(void)
{
	static qaws_parameterization const kinds[2] = { QAWS_PARAMETERIZATION_CHORDAL, QAWS_PARAMETERIZATION_CENTRIPETAL };
	static char const* const names[2] = { "chordal", "centripetal" };
	qaws_scalar pts[7 * 3] = { 0, 0, 0, 1, 3, 1, 2, -1, 0, 4, 2, 2, 7, 1, -1, 8, 4, 0, 10, 0, 1 };
	unsigned int k, closed;
	char msg[160];
	for (k = 0; k < 2; k++)
		for (closed = 0; closed < 2; closed++)
		{
			qaws_catmull_rom_desc d;
			qaws_curve* c = NULL;
			qaws_exact_curve* e = NULL;
			qaws_exact_report rep;
			qaws_status st;
			double gap;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_3D;
			d.control_points = pts;
			d.control_point_count = 7;
			d.parameterization = kinds[k];
			d.closed = (int)closed;
			qaws_curve_create_catmull_rom(&d, &c);
			st = qaws_exact_curve_prepare(NULL, c, &e, &rep);
			if (st != QAWS_STATUS_OK)
			{
				sprintf(msg, "%s%s Catmull-Rom: frozen preparation", names[k], closed ? " closed" : "");
				TEST_ASSERT(0, msg);
				qaws_curve_destroy(c);
				continue;
			}
			gap = ec_gap(c, e, 0, (double)e->span_count, 3);
			printf("    %s%s: %u spans, storage %u bits, vs its own evaluator %.2e\n", names[k], closed ? " closed" : "", e->span_count, rep.storage_bits,
				gap);
			sprintf(msg, "%s%s Catmull-Rom: the frozen preparation taken exactly (flagged), equal to its evaluator", names[k], closed ? " closed" : "");
			TEST_ASSERT(rep.flags == QAWS_EXACT_FLAG_PREP_QUANTIZED && gap < (QAWS_SCALAR_IS_FLOAT ? 1e-4 : 1e-12), msg);
			qaws_exact_curve_destroy(e);
			qaws_curve_destroy(c);
		}
}

int test_76_exact_composite_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 76: Exact composite curves and frozen Catmull-Rom preparations\n");
	test_composite();
	test_frozen();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
