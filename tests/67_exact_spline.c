/*
 * Test 67: Exact B-spline and NURBS evaluation (local Bezier extraction)
 *
 *   - C, C', C'', C''' as exact rationals against Mathematica
 *     (tests/reference/67_exact_spline.wls): degrees 1..5 and 7, 2D and 3D,
 *     polynomial and rational, clamped vectors with repeated interior knots
 *     and unclamped vectors, at the domain ends, at every interior knot and
 *     inside spans; compared by cross multiplication, no tolerance
 *   - the rounded-once doubles against the existing f64 evaluator
 *   - the storage width of the extracted spans per degree
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "exact/qaws_exact_curve.h"
#include "reference/67_exact_spline.h"
#include <math.h>
#include <string.h>
#include <stdint.h>

static qaws_curve* es_curve(ref_exact_spline const* r)
{
	qaws_scalar cps[16 * 3], ws[16], knots[32];
	qaws_curve* c = NULL;
	int i;
	for (i = 0; i < r->n * r->dim; i++)
		cps[i] = (qaws_scalar)ldexp((double)r->cp[i], -10);
	for (i = 0; i < r->n; i++)
		ws[i] = (qaws_scalar)r->w[i];
	for (i = 0; i < r->nk; i++)
		knots[i] = (qaws_scalar)((double)r->knots[i] / r->knot_den);
	if (r->rational)
	{
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)r->n;
		d.knots = knots;
		d.knot_count = (unsigned int)r->nk;
		d.weights = ws;
		d.weight_count = (unsigned int)r->n;
		qaws_curve_create_nurbs(&d, &c);
	}
	else
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
		d.degree = (unsigned int)r->degree;
		d.control_points = cps;
		d.control_point_count = (unsigned int)r->n;
		d.knots = knots;
		d.knot_count = (unsigned int)r->nk;
		qaws_curve_create_bspline(&d, &c);
	}
	return c;
}

static int es_same(qaws_exact_int const* num, qaws_exact_int const* den, char const* p_text, char const* q_text)
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
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_spline) / sizeof(g_ref_exact_spline[0])), i, storage[8];
	int ok_exact = 1, ok_prep = 1, bad = -1;
	double worst = 0;
	char msg[240];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	desc.param_bits = 12;
	memset(storage, 0, sizeof(storage));
	for (i = 0; i < n; i++)
	{
		ref_exact_spline const* r = &g_ref_exact_spline[i];
		qaws_curve* c = es_curve(r);
		qaws_exact_curve* e = NULL;
		qaws_exact_report rep;
		qaws_exact_int num[3], den;
		int comp;
		if (!c || qaws_exact_curve_prepare(&desc, c, &e, &rep) != QAWS_STATUS_OK)
		{
			ok_prep = 0;
			qaws_curve_destroy(c);
			continue;
		}
		if (rep.storage_bits > storage[r->degree])
			storage[r->degree] = rep.storage_bits;
		if (qaws_exact_curve_eval_rational(e, r->T, (unsigned int)r->k, num, &den) != QAWS_STATUS_OK)
			ok_exact = 0;
		else
			for (comp = 0; comp < r->dim; comp++)
				if (!es_same(&num[comp], &den, r->num[comp], r->den[comp]))
				{
					ok_exact = 0;
					if (bad < 0) bad = (int)i;
				}
		if (r->k <= 3)
		{
			double out[12], t = ldexp((double)r->T, -9);
			unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;
			qaws_eval_result_3d f;
			memset(&f, 0, sizeof(f));
			qaws_exact_curve_evaluate(e, t, (unsigned int)r->k, out, NULL);
			if (r->dim == 3)
				qaws_curve_evaluate_3d(c, (qaws_scalar)t, flags, &f);
			else
			{
				qaws_eval_result_2d f2;
				qaws_curve_evaluate_2d(c, (qaws_scalar)t, flags, &f2);
				f.position.x = f2.position.x; f.position.y = f2.position.y;
				f.d1.x = f2.d1.x; f.d1.y = f2.d1.y;
				f.d2.x = f2.d2.x; f.d2.y = f2.d2.y;
				f.d3.x = f2.d3.x; f.d3.y = f2.d3.y;
			}
			for (comp = 0; comp < r->dim; comp++)
			{
				qaws_vec3 v = r->k == 0 ? f.position : (r->k == 1 ? f.d1 : (r->k == 2 ? f.d2 : f.d3));
				double fv = comp == 0 ? v.x : (comp == 1 ? v.y : v.z), ex = out[r->k * r->dim + comp];
				double rel = fabs(fv - ex) / (1 + fabs(ex));
				if (rel > worst) worst = rel;
			}
		}
		qaws_exact_curve_destroy(e);
		qaws_curve_destroy(c);
	}
	printf("    %u exact rows; extracted span storage bits by degree: 1:%u 2:%u 3:%u 4:%u 5:%u 7:%u; f64 vs exact: worst %.2e\n", n,
		storage[1], storage[2], storage[3], storage[4], storage[5], storage[7], worst);
	if (bad >= 0)
		printf("    first mismatch: row %d (degree %d, rational %d, T %lld, k %d)\n", bad, g_ref_exact_spline[bad].degree,
			g_ref_exact_spline[bad].rational, g_ref_exact_spline[bad].T, g_ref_exact_spline[bad].k);
	TEST_ASSERT(ok_prep, "every B-spline / NURBS of the reference prepares");
	sprintf(msg, "C..C''' equal the Mathematica rationals on all %u rows (degrees 1..5, 7; repeated and unclamped knots)", n);
	TEST_ASSERT(ok_exact, msg);
	TEST_ASSERT(worst < (QAWS_SCALAR_IS_FLOAT ? 5e-2 : 1e-7), "rounded-once values agree with the f64 evaluator");
}

/* Width of the extracted spans for full-precision inputs: random NURBS of
   every degree with 24-bit knots, 26-bit coordinates and 24-bit weights;
   every one must prepare and evaluate its third derivative within the
   2048-bit integers or report EXACT_RANGE_EXCEEDED (never wrap). */
static void test_widths(void)
{
	unsigned int p, i, report_bits[17], eval_ok[17], max_order[17];
	int ok = 1;
	uint64_t state = 0x853C49E6748FEA9Bull;
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -20;
	for (p = 1; p <= 16; p++)
	{
		qaws_scalar cps[24 * 2], ws[24], knots[48];
		unsigned int n = p + 4, nk = n + p + 1;
		qaws_nurbs_desc d;
		qaws_curve* c = NULL;
		qaws_exact_curve* e = NULL;
		qaws_exact_report rep;
		qaws_status st;
		for (i = 0; i < n * 2; i++)
		{
			state ^= state << 13; state ^= state >> 7; state ^= state << 17;
			cps[i] = (qaws_scalar)ldexp((double)(int64_t)(state % (1u << 26)) - (1 << 25), -20);
		}
		for (i = 0; i < n; i++)
		{
			state ^= state << 13; state ^= state >> 7; state ^= state << 17;
			ws[i] = (qaws_scalar)(1 + (double)(state % 16000000u));
		}
		/* clamped, interior knots random in (0, 1) on the 2^-24 grid */
		for (i = 0; i < nk; i++)
		{
			state ^= state << 13; state ^= state >> 7; state ^= state << 17;
			knots[i] = i <= p ? 0 : (i >= n ? 1 : (qaws_scalar)ldexp((double)(1 + state % ((1u << 24) - 2)), -24));
		}
		for (i = p + 1; i < n; i++)   /* sort the interior */
		{
			unsigned int j;
			for (j = i; j > p + 1 && knots[j - 1] > knots[j]; j--)
			{
				qaws_scalar tmp = knots[j];
				knots[j] = knots[j - 1];
				knots[j - 1] = tmp;
			}
		}
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D;
		d.degree = p;
		d.control_points = cps;
		d.control_point_count = n;
		d.knots = knots;
		d.knot_count = nk;
		d.weights = ws;
		d.weight_count = n;
		qaws_curve_create_nurbs(&d, &c);
		st = qaws_exact_curve_prepare(&desc, c, &e, &rep);
		report_bits[p] = st == QAWS_STATUS_OK ? rep.storage_bits : 0;
		eval_ok[p] = 0;
		max_order[p] = 0;
		if (st == QAWS_STATUS_OK)
		{
			double out[8];
			unsigned int k;
			for (k = 0; k <= 3; k++)
			{
				st = qaws_exact_curve_evaluate(e, 0.37, k, out, NULL);
				if (st != QAWS_STATUS_OK)
					break;
				max_order[p] = k + 1;
			}
			eval_ok[p] = max_order[p] == 4;
		}
		ok &= st == QAWS_STATUS_OK || st == QAWS_STATUS_EXACT_RANGE_EXCEEDED;
		qaws_exact_curve_destroy(e);
		qaws_curve_destroy(c);
	}
	printf("    full-precision NURBS (24-bit knots and weights, 26-bit coordinates): storage bits / C''' evaluates\n     ");
	for (p = 1; p <= 16; p++)
		printf(" p%u:%u/C%u", p, report_bits[p], max_order[p] ? max_order[p] - 1 : 99);
	printf("\n");
	TEST_ASSERT(ok, "every degree 1..16 prepares or reports the integer range, never wraps");
}

int test_67_exact_spline_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 67: Exact B-spline and NURBS evaluation\n");
	test_reference();
	test_widths();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
