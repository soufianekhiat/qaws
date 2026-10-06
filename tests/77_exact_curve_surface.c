/*
 * Test 77: Certified curve / surface intersections
 *
 *   - against Mathematica (tests/reference/77_exact_curve_surface.wls):
 *     3D Bezier and rational Bezier curves of degree 1..3 against Bezier
 *     patches up to bicubic (polynomial and rational) and a 2 x 2-patch
 *     NURBS surface; every hit found once and enclosed (t, u, v)
 *   - a curve lying on the surface is refused
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "reference/77_exact_curve_surface.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static qaws_surface* ecs_surface(ref_exact_curve_surface const* r)
{
	qaws_vec3 cps[16];
	qaws_scalar ws[16], ku[8], kv[8];
	qaws_surface* s = NULL;
	qaws_surface_nurbs_desc d;
	int i;
	for (i = 0; i < r->nu * r->nv; i++)
	{
		cps[i].x = (qaws_scalar)ldexp((double)r->cp[3 * i], -10);
		cps[i].y = (qaws_scalar)ldexp((double)r->cp[3 * i + 1], -10);
		cps[i].z = (qaws_scalar)ldexp((double)r->cp[3 * i + 2], -10);
		ws[i] = (qaws_scalar)r->w[i];
	}
	for (i = 0; i < r->nku; i++) ku[i] = (qaws_scalar)r->ku[i];
	for (i = 0; i < r->nkv; i++) kv[i] = (qaws_scalar)r->kv[i];
	memset(&d, 0, sizeof(d));
	d.u_degree = (unsigned int)r->p;
	d.v_degree = (unsigned int)r->q;
	d.control_points = cps;
	d.u_point_count = (unsigned int)r->nu;
	d.v_point_count = (unsigned int)r->nv;
	d.weights = ws;
	d.u_knots = ku;
	d.u_knot_count = (unsigned int)r->nku;
	d.v_knots = kv;
	d.v_knot_count = (unsigned int)r->nkv;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

static qaws_curve* ecs_curve(ref_exact_curve_surface const* r)
{
	qaws_scalar cps[4 * 3], ws[4];
	qaws_rational_bezier_desc d;
	qaws_curve* c = NULL;
	int i;
	for (i = 0; i < (r->m + 1) * 3; i++)
		cps[i] = (qaws_scalar)ldexp((double)r->ccp[i], -10);
	for (i = 0; i <= r->m; i++)
		ws[i] = (qaws_scalar)r->cw[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = (unsigned int)r->m;
	d.control_points = cps;
	d.control_point_count = (unsigned int)r->m + 1;
	d.weights = ws;
	d.weight_count = (unsigned int)r->m + 1;
	qaws_curve_create_rational_bezier(&d, &c);
	return c;
}

static int ecs_in(double lo, double hi, char const* text)
{
	double v = strtod(text, NULL);
	return lo <= v && v <= hi;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_curve_surface) / sizeof(g_ref_exact_curve_surface[0])), i, total = 0, bad = 0;
	double widest = 0;
	char msg[200];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	for (i = 0; i < n; i++)
	{
		ref_exact_curve_surface const* r = &g_ref_exact_curve_surface[i];
		qaws_surface* s = ecs_surface(r);
		qaws_curve* c = ecs_curve(r);
		qaws_exact_surface* es = NULL;
		qaws_exact_curve* ec = NULL;
		qaws_exact_curve_surface_hit hits[16];
		unsigned int count = 0, k, q;
		qaws_status st;
		int ok;
		qaws_exact_surface_prepare(&desc, s, &es, NULL);
		qaws_exact_curve_prepare(&desc, c, &ec, NULL);
		st = qaws_exact_curve_surface_hits(ec, es, hits, 16, &count);
		ok = st == QAWS_STATUS_OK && count == (unsigned int)r->n;
		/* hits are unordered: each reference hit inside exactly one enclosure */
		for (k = 0; ok && k < (unsigned int)r->n; k++)
		{
			unsigned int found = 0;
			for (q = 0; q < count; q++)
				found += ecs_in(hits[q].t_lo, hits[q].t_hi, r->t[k]) && ecs_in(hits[q].u_lo, hits[q].u_hi, r->u[k]) &&
				         ecs_in(hits[q].v_lo, hits[q].v_hi, r->v[k]);
			ok &= found == 1;
		}
		for (q = 0; q < count; q++)
		{
			if (hits[q].t_hi - hits[q].t_lo > widest) widest = hits[q].t_hi - hits[q].t_lo;
			if (hits[q].u_hi - hits[q].u_lo > widest) widest = hits[q].u_hi - hits[q].u_lo;
			if (hits[q].v_hi - hits[q].v_lo > widest) widest = hits[q].v_hi - hits[q].v_lo;
		}
		total += (unsigned int)r->n;
		if (!ok)
		{
			bad++;
			printf("      case %u (%d x %d%s, curve degree %d): status %d, %u hits, expected %d\n", i, r->p, r->q, r->rational ? " rational" : "", r->m,
				(int)st, count, r->n);
			for (q = 0; q < count && q < 6; q++)
				printf("        t [%.17g, %.17g] u [%.17g, %.17g] v [%.17g, %.17g]\n", hits[q].t_lo, hits[q].t_hi, hits[q].u_lo, hits[q].u_hi, hits[q].v_lo,
					hits[q].v_hi);
			for (k = 0; k < (unsigned int)r->n && k < 6; k++)
				printf("        ref t %s u %s v %s\n", r->t[k], r->u[k], r->v[k]);
		}
		qaws_exact_surface_destroy(es);
		qaws_exact_curve_destroy(ec);
		qaws_surface_destroy(s);
		qaws_curve_destroy(c);
	}
	printf("    %u curve / surface pairs: %u hits, widest enclosure %.1e\n", n, total, widest);
	sprintf(msg, "every curve / surface hit found once and enclosed (t, u, v) on all %u Mathematica pairs", n);
	TEST_ASSERT(bad == 0, msg);
	TEST_ASSERT(widest < 1e-9, "t, u and v enclosed tightly");
}

static void test_on_surface(void)
{
	/* a bilinear patch z = 0 and a segment inside it */
	qaws_vec3 cps[4];
	qaws_scalar seg[6] = { 0.25, 0.25, 0, 0.75, 0.5, 0 };
	qaws_surface_bezier_desc d;
	qaws_bezier_desc bd;
	qaws_surface* s = NULL;
	qaws_curve* c = NULL;
	qaws_exact_surface* es = NULL;
	qaws_exact_curve* ec = NULL;
	qaws_exact_curve_surface_hit hits[4];
	unsigned int count = 0;
	cps[0].x = 0; cps[0].y = 0; cps[0].z = 0;
	cps[1].x = 0; cps[1].y = 1; cps[1].z = 0;
	cps[2].x = 1; cps[2].y = 0; cps[2].z = 0;
	cps[3].x = 1; cps[3].y = 1; cps[3].z = 0;
	memset(&d, 0, sizeof(d));
	d.u_degree = 1;
	d.v_degree = 1;
	d.control_points = cps;
	d.u_point_count = 2;
	d.v_point_count = 2;
	qaws_surface_create_bezier(&d, &s);
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_3D;
	bd.degree = 1;
	bd.control_points = seg;
	bd.control_point_count = 2;
	qaws_curve_create_bezier(&bd, &c);
	qaws_exact_surface_prepare(NULL, s, &es, NULL);
	qaws_exact_curve_prepare(NULL, c, &ec, NULL);
	TEST_ASSERT(qaws_exact_curve_surface_hits(ec, es, hits, 4, &count) == QAWS_STATUS_CERTIFICATION_FAILED,
		"a curve lying on the surface is refused (the subdivision budget runs out)");
	qaws_exact_surface_destroy(es);
	qaws_exact_curve_destroy(ec);
	qaws_surface_destroy(s);
	qaws_curve_destroy(c);
}

int test_77_exact_curve_surface_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 77: Certified curve / surface intersections\n");
	test_reference();
	test_on_surface();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
