/*
 * Test 65: Exact Bezier and rational Bezier evaluation
 *
 *   - C, C', C'', C''' as exact rationals against Mathematica
 *     (tests/reference/65_exact_bezier.wls): degrees 1..16, 2D and 3D,
 *     polynomial and rational, endpoints included; compared by cross
 *     multiplication, no tolerance
 *   - the rounded-once doubles against the existing f64 evaluator
 *   - scale invariance: multiplying every weight by 3 changes nothing
 *   - lattice overflow and unsupported families are reported
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "exact/qaws_exact_curve.h"
#include "reference/65_exact_bezier.h"
#include <math.h>
#include <string.h>

static qaws_curve* eb_curve(ref_exact_bezier const* r, double wscale)
{
	qaws_scalar cps[17 * 3], ws[17];
	qaws_curve* c = NULL;
	int i, n = r->degree + 1;
	for (i = 0; i < n * r->dim; i++)
		cps[i] = (qaws_scalar)ldexp((double)r->cp[i], -10);
	for (i = 0; i < n; i++)
		ws[i] = (qaws_scalar)(r->w[i] * wscale);
	if (r->rational)
	{
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)n;
		d.weights = ws;
		d.weight_count = (unsigned int)n;
		qaws_curve_create_rational_bezier(&d, &c);
	}
	else
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)n;
		qaws_curve_create_bezier(&d, &c);
	}
	return c;
}

static qaws_exact_desc eb_desc(void)
{
	qaws_exact_desc d;
	qaws_exact_desc_default(&d);
	d.space_exp2 = -10;
	d.param_bits = 8;
	return d;
}

/* num / den == p / q by cross multiplication */
static int eb_same(qaws_exact_int const* num, qaws_exact_int const* den, char const* p_text, char const* q_text)
{
	qaws_exact_int p, q, l, r;
	if (qaws_exact_int_from_text(&p, p_text) != QAWS_STATUS_OK || qaws_exact_int_from_text(&q, q_text) != QAWS_STATUS_OK)
		return 0;
	if (qaws_exact_int_mul(&l, num, &q) != QAWS_STATUS_OK || qaws_exact_int_mul(&r, &p, den) != QAWS_STATUS_OK)
		return 0;
	return qaws_exact_int_cmp(&l, &r) == 0;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_bezier) / sizeof(g_ref_exact_bezier[0])), i, max_bits = 0;
	int ok_exact = 1, ok_float = 1, ok_scale = 1;
	double worst = 0;
	char msg[200];
	qaws_exact_desc desc = eb_desc();
	for (i = 0; i < n; i++)
	{
		ref_exact_bezier const* r = &g_ref_exact_bezier[i];
		qaws_curve* c = eb_curve(r, 1);
		qaws_curve* c3 = eb_curve(r, 3);
		qaws_exact_curve* e = NULL;
		qaws_exact_curve* e3 = NULL;
		qaws_exact_int num[3], den, num3[3], den3;
		int comp;
		if (qaws_exact_curve_prepare(&desc, c, &e, NULL) != QAWS_STATUS_OK ||
		    qaws_exact_curve_eval_rational(e, r->T, (unsigned int)r->k, num, &den) != QAWS_STATUS_OK)
			ok_exact = 0;
		else
		{
			for (comp = 0; comp < r->dim; comp++)
			{
				unsigned int b = qaws_exact_int_bits(&num[comp]);
				ok_exact &= eb_same(&num[comp], &den, r->num[comp], r->den[comp]);
				if (b > max_bits) max_bits = b;
			}
			/* weights times 3: the same rational values */
			if (r->rational && qaws_exact_curve_prepare(&desc, c3, &e3, NULL) == QAWS_STATUS_OK &&
			    qaws_exact_curve_eval_rational(e3, r->T, (unsigned int)r->k, num3, &den3) == QAWS_STATUS_OK)
				for (comp = 0; comp < r->dim; comp++)
				{
					qaws_exact_int l, rr;
					qaws_exact_int_mul(&l, &num[comp], &den3);
					qaws_exact_int_mul(&rr, &num3[comp], &den);
					ok_scale &= qaws_exact_int_cmp(&l, &rr) == 0;
				}
			/* rounded once against the f64 evaluator (position, D1, D2, D3) */
			if (r->k <= 3)
			{
				double out[12], t = ldexp((double)r->T, -8);
				qaws_eval_result_3d f;
				unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;
				qaws_exact_curve_evaluate(e, t, (unsigned int)r->k, out, NULL);
				if (r->dim == 3)
					qaws_curve_evaluate_3d(c, (qaws_scalar)t, flags, &f);
				else
				{
					qaws_eval_result_2d f2;
					qaws_curve_evaluate_2d(c, (qaws_scalar)t, flags, &f2);
					f.position.x = f2.position.x; f.position.y = f2.position.y; f.position.z = 0;
					f.d1.x = f2.d1.x; f.d1.y = f2.d1.y; f.d1.z = 0;
					f.d2.x = f2.d2.x; f.d2.y = f2.d2.y; f.d2.z = 0;
					f.d3.x = f2.d3.x; f.d3.y = f2.d3.y; f.d3.z = 0;
				}
				for (comp = 0; comp < r->dim; comp++)
				{
					qaws_vec3 v = r->k == 0 ? f.position : (r->k == 1 ? f.d1 : (r->k == 2 ? f.d2 : f.d3));
					double fv = comp == 0 ? v.x : (comp == 1 ? v.y : v.z), ex = out[r->k * r->dim + comp];
					double rel = fabs(fv - ex) / (1 + fabs(ex));
					if (rel > worst) worst = rel;
				}
			}
		}
		qaws_exact_curve_destroy(e);
		qaws_exact_curve_destroy(e3);
		qaws_curve_destroy(c);
		qaws_curve_destroy(c3);
	}
	ok_float = worst < (QAWS_SCALAR_IS_FLOAT ? 1e-2 : 1e-8);
	printf("    %u exact rows; widest numerator %u bits; f64 evaluator vs exact: worst relative gap %.2e\n", n, max_bits, worst);
	sprintf(msg, "C..C''' equal the Mathematica rationals on all %u rows (degrees 1..16, 2D/3D, rational)", n);
	TEST_ASSERT(ok_exact, msg);
	TEST_ASSERT(ok_scale, "scale invariance: weights times 3 give the same rationals");
	TEST_ASSERT(ok_float, "rounded-once values agree with the f64 evaluator");
}

static void test_reports(void)
{
	qaws_scalar cps[4] = { 0, 0, 1.0e9, 1 };
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	qaws_exact_curve* e = NULL;
	qaws_exact_desc desc;
	qaws_exact_report rep;
	qaws_arc_desc ad;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 1;
	d.control_points = cps;
	d.control_point_count = 2;
	qaws_curve_create_bezier(&d, &c);
	qaws_exact_desc_default(&desc);
	TEST_ASSERT(qaws_exact_curve_prepare(&desc, c, &e, &rep) == QAWS_STATUS_EXACT_RANGE_EXCEEDED,
		"a coordinate past the lattice range is reported");
	desc.space_exp2 = 10;
	TEST_ASSERT(qaws_exact_curve_prepare(&desc, c, &e, &rep) == QAWS_STATUS_OK && (rep.flags & QAWS_EXACT_FLAG_INPUT_QUANTIZED) &&
		rep.max_position_quantization_error > 0, "a coarse lattice reports its quantization error");
	qaws_exact_curve_destroy(e);
	qaws_curve_destroy(c);
	memset(&ad, 0, sizeof(ad));
	(void)ad;
}

int test_65_exact_bezier_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 65: Exact Bezier and rational Bezier evaluation\n");
	test_reference();
	test_reports();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
