/*
 * Test 68: Exact Hermite, uniform Catmull-Rom and polynomial curves
 *
 *   - C, C', C'', C''' as exact rationals against Mathematica
 *     (tests/reference/68_exact_families.wls): cubic Hermite, uniform
 *     Catmull-Rom open and closed, polynomials of degree 1..8 on shifted
 *     dyadic domains; cross multiplication, no tolerance
 *   - the rounded-once doubles against the f64 evaluator of each family
 *     (which also pins the Catmull-Rom span convention)
 *   - polynomial coefficients are taken exactly (no quantization reported)
 *   (chordal / centripetal Catmull-Rom: test 76, frozen preparation)
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "exact/qaws_exact_curve.h"
#include "reference/68_exact_families.h"
#include <math.h>
#include <string.h>

static qaws_curve* ef_curve(ref_exact_family const* r)
{
	qaws_scalar a[9 * 3], b[8 * 3];
	qaws_curve* c = NULL;
	qaws_dimension dim = r->dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
	int i;
	if (r->family == 0)
	{
		qaws_hermite_desc d;
		for (i = 0; i < r->n * r->dim; i++)
		{
			a[i] = (qaws_scalar)ldexp((double)r->p[i], -10);
			b[i] = (qaws_scalar)ldexp((double)r->m[i], -10);
		}
		memset(&d, 0, sizeof(d));
		d.dimension = dim;
		d.degree = 3;
		d.points = a;
		d.derivatives = b;
		d.point_count = (unsigned int)r->n;
		d.derivative_count = (unsigned int)r->n;
		qaws_curve_create_hermite(&d, &c);
	}
	else if (r->family <= 2)
	{
		qaws_catmull_rom_desc d;
		for (i = 0; i < r->n * r->dim; i++)
			a[i] = (qaws_scalar)ldexp((double)r->p[i], -10);
		memset(&d, 0, sizeof(d));
		d.dimension = dim;
		d.control_points = a;
		d.control_point_count = (unsigned int)r->n;
		d.parameterization = QAWS_PARAMETERIZATION_UNIFORM;
		d.closed = r->family == 2;
		qaws_curve_create_catmull_rom(&d, &c);
	}
	else
	{
		qaws_polynomial_desc d;
		for (i = 0; i < r->n * r->dim; i++)
			a[i] = (qaws_scalar)ldexp((double)r->cf[i], -12);
		memset(&d, 0, sizeof(d));
		d.dimension = dim;
		d.degree = (unsigned int)r->n - 1;
		d.coefficients = a;
		d.coefficient_count = (unsigned int)r->n;
		d.t_min = (qaws_scalar)(r->t0n / 16.0);
		d.t_max = (qaws_scalar)(r->t1n / 16.0);
		qaws_curve_create_polynomial(&d, &c);
	}
	return c;
}

static qaws_exact_desc ef_desc(void)
{
	qaws_exact_desc d;
	qaws_exact_desc_default(&d);
	d.space_exp2 = -10;
	return d;
}

static int ef_same(qaws_exact_int const* num, qaws_exact_int const* den, char const* p_text, char const* q_text)
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
	static char const* const names[4] = { "Hermite", "Catmull-Rom open", "Catmull-Rom closed", "polynomial" };
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_family) / sizeof(g_ref_exact_family[0])), i, rows[4] = { 0 }, bad[4] = { 0 };
	double worst[4] = { 0 };
	int ok_float = 1, ok_report = 1, fam;
	char msg[200];
	qaws_exact_desc desc = ef_desc();
	for (i = 0; i < n; i++)
	{
		ref_exact_family const* r = &g_ref_exact_family[i];
		qaws_curve* c = ef_curve(r);
		qaws_exact_curve* e = NULL;
		qaws_exact_report rep;
		qaws_exact_int num[3], den;
		int comp, ok = 1;
		rows[r->family]++;
		if (!c || qaws_exact_curve_prepare(&desc, c, &e, &rep) != QAWS_STATUS_OK || e->param_shift < 6 ||
		    qaws_exact_curve_eval_rational(e, r->T << (e->param_shift - 6), (unsigned int)r->k, num, &den) != QAWS_STATUS_OK)
			ok = 0;
		else
		{
			for (comp = 0; comp < r->dim; comp++)
				ok &= ef_same(&num[comp], &den, r->num[comp], r->den[comp]);
			/* lattice inputs and dyadic coefficients: nothing was rounded */
			ok_report &= rep.quality == QAWS_NUMERIC_EXACT_RATIONAL && rep.flags == QAWS_EXACT_FLAG_NONE;
			{
				double out[12], t = ldexp((double)r->T, -6);
				qaws_eval_result_3d f;
				unsigned int flags = QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3;
				memset(&f, 0, sizeof(f));
				qaws_exact_curve_evaluate(e, t, (unsigned int)r->k, out, NULL);
				if (r->dim == 3)
					qaws_curve_evaluate_3d(c, (qaws_scalar)t, flags, &f);
				else
				{
					qaws_eval_result_2d f2;
					memset(&f2, 0, sizeof(f2));
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
					if (rel > worst[r->family]) worst[r->family] = rel;
				}
			}
		}
		if (!ok)
			bad[r->family]++;
		qaws_exact_curve_destroy(e);
		qaws_curve_destroy(c);
	}
	for (fam = 0; fam < 4; fam++)
	{
		printf("    %-19s %4u rows, %u mismatches; scalar evaluator vs exact: worst relative gap %.2e\n", names[fam], rows[fam], bad[fam], worst[fam]);
		sprintf(msg, "%s: C..C''' equal the Mathematica rationals on all %u rows", names[fam], rows[fam]);
		TEST_ASSERT(rows[fam] > 0 && bad[fam] == 0, msg);
		ok_float &= worst[fam] < (QAWS_SCALAR_IS_FLOAT ? 1e-2 : 1e-8);
	}
	TEST_ASSERT(ok_report, "lattice points and dyadic coefficients: exact, nothing quantized");
	TEST_ASSERT(ok_float, "rounded-once values agree with each family's f64 evaluator");
}


int test_68_exact_families_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 68: Exact Hermite, uniform Catmull-Rom and polynomial curves\n");
	test_reference();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
