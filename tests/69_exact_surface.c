/*
 * Test 69: Exact tensor-product surfaces
 *
 *   - S, Su, Sv, Suu, Suv, Svv and Su x Sv as exact rationals against
 *     Mathematica (tests/reference/69_exact_surface.wls): Bezier patches up
 *     to degree (5, 4), B-spline and NURBS surfaces with repeated interior
 *     knots; cross multiplication, no tolerance
 *   - the rounded-once doubles against the f64 surface evaluator
 *   - patch counts equal the non-empty knot rectangles
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "exact/qaws_exact_surface.h"
#include "reference/69_exact_surface.h"
#include <math.h>
#include <string.h>

static qaws_surface* es_surface(ref_exact_surface_case const* r)
{
	qaws_vec3 cps[8 * 8];
	qaws_scalar ws[8 * 8], ku[16], kv[16];
	qaws_surface* s = NULL;
	int i, n = r->nu * r->nv;
	for (i = 0; i < n; i++)
	{
		cps[i].x = (qaws_scalar)ldexp((double)r->cp[3 * i], -10);
		cps[i].y = (qaws_scalar)ldexp((double)r->cp[3 * i + 1], -10);
		cps[i].z = (qaws_scalar)ldexp((double)r->cp[3 * i + 2], -10);
		ws[i] = (qaws_scalar)r->w[i];
	}
	for (i = 0; i < r->nu + r->p + 1; i++)
		ku[i] = (qaws_scalar)r->ku[i];
	for (i = 0; i < r->nv + r->q + 1; i++)
		kv[i] = (qaws_scalar)r->kv[i];
	if (r->kind == 0)
	{
		qaws_surface_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.u_degree = (unsigned int)r->p;
		d.v_degree = (unsigned int)r->q;
		d.control_points = cps;
		d.u_point_count = (unsigned int)r->nu;
		d.v_point_count = (unsigned int)r->nv;
		qaws_surface_create_bezier(&d, &s);
	}
	else if (r->kind == 1)
	{
		qaws_surface_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.u_degree = (unsigned int)r->p;
		d.v_degree = (unsigned int)r->q;
		d.control_points = cps;
		d.u_point_count = (unsigned int)r->nu;
		d.v_point_count = (unsigned int)r->nv;
		d.u_knots = ku;
		d.u_knot_count = (unsigned int)(r->nu + r->p + 1);
		d.v_knots = kv;
		d.v_knot_count = (unsigned int)(r->nv + r->q + 1);
		qaws_surface_create_bspline(&d, &s);
	}
	else
	{
		qaws_surface_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.u_degree = (unsigned int)r->p;
		d.v_degree = (unsigned int)r->q;
		d.control_points = cps;
		d.u_point_count = (unsigned int)r->nu;
		d.v_point_count = (unsigned int)r->nv;
		d.weights = ws;
		d.u_knots = ku;
		d.u_knot_count = (unsigned int)(r->nu + r->p + 1);
		d.v_knots = kv;
		d.v_knot_count = (unsigned int)(r->nv + r->q + 1);
		qaws_surface_create_nurbs(&d, &s);
	}
	return s;
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

/* non-empty knot intervals of a clamped vector */
static unsigned int es_spans(long long const* k, int deg, int n)
{
	unsigned int c = 0;
	int i;
	for (i = deg; i < n; i++)
		if (k[i] < k[i + 1])
			c++;
	return c;
}

static void test_reference(void)
{
	static char const* const names[7] = { "S", "Su", "Sv", "Suu", "Suv", "Svv", "Su x Sv" };
	static unsigned int const da[6] = { 0, 1, 0, 2, 1, 0 }, db[6] = { 0, 0, 1, 0, 1, 2 };
	unsigned int ncase = (unsigned int)(sizeof(g_ref_exact_surface_cases) / sizeof(g_ref_exact_surface_cases[0]));
	unsigned int nrow = (unsigned int)(sizeof(g_ref_exact_surface_rows) / sizeof(g_ref_exact_surface_rows[0]));
	unsigned int i, bad[7] = { 0 }, rows[7] = { 0 }, max_bits = 0;
	int ok_count = 1, ok_float = 1, k;
	double worst = 0;
	char msg[200];
	qaws_surface* surf[32];
	qaws_exact_surface* ex[32];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	for (i = 0; i < ncase; i++)
	{
		ref_exact_surface_case const* r = &g_ref_exact_surface_cases[i];
		qaws_exact_report rep;
		unsigned int pu = 0, pv = 0;
		surf[i] = es_surface(r);
		ex[i] = NULL;
		if (qaws_exact_surface_prepare(&desc, surf[i], &ex[i], &rep) != QAWS_STATUS_OK)
		{
			ok_count = 0;
			continue;
		}
		if (rep.storage_bits > max_bits) max_bits = rep.storage_bits;
		qaws_exact_surface_patch_count(ex[i], &pu, &pv);
		ok_count &= pu == es_spans(r->ku, r->p, r->nu) && pv == es_spans(r->kv, r->q, r->nv);
	}
	for (i = 0; i < nrow; i++)
	{
		ref_exact_surface_row const* r = &g_ref_exact_surface_rows[i];
		qaws_exact_surface const* e = ex[r->c];
		qaws_exact_int num[3], den;
		qaws_status st;
		int c, ok = 1;
		rows[r->which]++;
		if (!e)
		{
			bad[r->which]++;
			continue;
		}
		{
			int64_t U = (int64_t)r->U << (e->u_shift - 6), V = (int64_t)r->V << (e->v_shift - 6);
			st = r->which < 6 ? qaws_exact_surface_eval_rational(e, U, V, da[r->which], db[r->which], num, &den)
			                  : qaws_exact_surface_normal_rational(e, U, V, num, &den);
		}
		if (st != QAWS_STATUS_OK)
			ok = 0;
		for (c = 0; c < 3 && ok; c++)
			ok &= es_same(&num[c], &den, r->num[c], r->den[c]);
		if (!ok)
			bad[r->which]++;
		/* rounded once against the f64 evaluator */
		if (ok && r->which == 0)
		{
			double out[18], nrm[3], u = ldexp((double)r->U, -6), v = ldexp((double)r->V, -6);
			qaws_surface_eval_result f;
			unsigned int flags = QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DUU |
			                     QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV;
			qaws_vec3 fv[6];
			int j;
			qaws_exact_surface_evaluate(e, u, v, 2, out, nrm);
			qaws_surface_evaluate(surf[r->c], (qaws_scalar)u, (qaws_scalar)v, flags, &f);
			fv[0] = f.position; fv[1] = f.du; fv[2] = f.dv; fv[3] = f.duu; fv[4] = f.duv; fv[5] = f.dvv;
			for (j = 0; j < 6; j++)
			{
				double m = 1 + fabs(out[3 * j]) + fabs(out[3 * j + 1]) + fabs(out[3 * j + 2]);
				double g = fabs(fv[j].x - out[3 * j]) + fabs(fv[j].y - out[3 * j + 1]) + fabs(fv[j].z - out[3 * j + 2]);
				if (g / m > worst) worst = g / m;
			}
		}
	}
	for (k = 0; k < 7; k++)
	{
		printf("    %-8s %4u rows, %u mismatches\n", names[k], rows[k], bad[k]);
		sprintf(msg, "%s equals the Mathematica rationals on all %u rows", names[k], rows[k]);
		TEST_ASSERT(rows[k] > 0 && bad[k] == 0, msg);
	}
	ok_float = worst < (QAWS_SCALAR_IS_FLOAT ? 1e-2 : 1e-8);
	printf("    %u surfaces (Bezier to degree (5, 4), B-spline, NURBS); widest stored integer %u bits; scalar evaluator vs exact: worst relative gap %.2e\n",
		ncase, max_bits, worst);
	TEST_ASSERT(ok_count, "one exact patch per non-empty knot rectangle");
	TEST_ASSERT(ok_float, "rounded-once partials agree with the surface evaluator");
	for (i = 0; i < ncase; i++)
	{
		qaws_exact_surface_destroy(ex[i]);
		qaws_surface_destroy(surf[i]);
	}
}

int test_69_exact_surface_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 69: Exact tensor-product surfaces\n");
	test_reference();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
