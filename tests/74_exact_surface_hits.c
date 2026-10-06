/*
 * Test 74: Certified line / surface intersections (exact ray casting)
 *
 *   - against Mathematica (tests/reference/74_exact_surface_hits.wls):
 *     Bezier patches up to bicubic, polynomial and rational, multi-patch
 *     B-spline and NURBS surfaces; lines through rounded surface points,
 *     random lines, lines through an exact surface point. Every hit found
 *     once, u, v and the line parameter t enclosed, dyadic hits exact
 *   - a line lying on a bilinear patch is refused (overlap)
 */

#include "test_common.h"
#include "qaws_exact.h"
#include "reference/74_exact_surface_hits.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double eh_rational(char const* text)
{
	char const* slash = strchr(text, '/');
	double n = strtod(text, NULL);
	return slash ? n / strtod(slash + 1, NULL) : n;
}

static qaws_surface* eh_surface(ref_exact_surface_hits const* r)
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

static int eh_in(double lo, double hi, char const* text)
{
	double v = strtod(text, NULL);
	return lo <= v && v <= hi;
}

static void test_reference(void)
{
	unsigned int n = (unsigned int)(sizeof(g_ref_exact_surface_hits) / sizeof(g_ref_exact_surface_hits[0])), i, total = 0, exact = 0, bad = 0;
	double widest = 0;
	char msg[200];
	qaws_exact_desc desc;
	qaws_exact_desc_default(&desc);
	desc.space_exp2 = -10;
	for (i = 0; i < n; i++)
	{
		ref_exact_surface_hits const* r = &g_ref_exact_surface_hits[i];
		qaws_surface* s = eh_surface(r);
		qaws_exact_surface* e = NULL;
		qaws_exact_surface_hit hits[16];
		double p0[3], p1[3];
		unsigned int count = 0, k;
		qaws_status st;
		int ok;
		for (k = 0; k < 3; k++)
		{
			p0[k] = ldexp(eh_rational(r->p0[k]), -10);
			p1[k] = ldexp(eh_rational(r->p1[k]), -10);
		}
		qaws_exact_surface_prepare(&desc, s, &e, NULL);
		st = qaws_exact_surface_line_hits(e, p0, p1, hits, 16, &count);
		ok = st == QAWS_STATUS_OK && count == (unsigned int)r->n;
		for (k = 0; ok && k < count; k++)
		{
			ok &= eh_in(hits[k].u_lo, hits[k].u_hi, r->u[k]) && eh_in(hits[k].v_lo, hits[k].v_hi, r->v[k]) && eh_in(hits[k].t_lo, hits[k].t_hi, r->t[k]);
			if (r->dy[k])
				ok &= hits[k].kind == QAWS_EXACT_HIT_POINT;
			exact += hits[k].kind == QAWS_EXACT_HIT_POINT;
			if (hits[k].u_hi - hits[k].u_lo > widest) widest = hits[k].u_hi - hits[k].u_lo;
			if (hits[k].v_hi - hits[k].v_lo > widest) widest = hits[k].v_hi - hits[k].v_lo;
			if (hits[k].t_hi - hits[k].t_lo > widest) widest = hits[k].t_hi - hits[k].t_lo;
		}
		total += (unsigned int)r->n;
		if (!ok)
		{
			bad++;
			printf("      case %u (%d x %d%s, %d x %d points): status %d, %u hits, expected %d\n", i, r->p, r->q, r->rational ? " rational" : "", r->nu,
				r->nv, (int)st, count, r->n);
			for (k = 0; k < count && k < 6; k++)
				printf("        u [%.17g, %.17g] v [%.17g, %.17g] t [%.17g, %.17g]  vs  %s %s %s\n", hits[k].u_lo, hits[k].u_hi, hits[k].v_lo, hits[k].v_hi,
					hits[k].t_lo, hits[k].t_hi, k < (unsigned int)r->n ? r->u[k] : "-", k < (unsigned int)r->n ? r->v[k] : "-",
					k < (unsigned int)r->n ? r->t[k] : "-");
		}
		qaws_exact_surface_destroy(e);
		qaws_surface_destroy(s);
	}
	printf("    %u lines: %u hits (%u exact points), widest enclosure %.1e\n", n, total, exact, widest);
	sprintf(msg, "every line / surface hit found once and enclosed (u, v, t) on all %u Mathematica cases", n);
	TEST_ASSERT(bad == 0, msg);
	TEST_ASSERT(widest < 1e-9, "u, v and t enclosed tightly");
}

static void test_overlap(void)
{
	/* a bilinear patch z = 0 and a line in it */
	qaws_vec3 cps[4];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	qaws_exact_surface* e = NULL;
	qaws_exact_surface_hit hits[4];
	unsigned int count = 0;
	double a[3] = { -1, 0.25, 0 }, b[3] = { 2, 0.75, 0 }, c[3] = { 0.5, 0.5, -1 }, dd[3] = { 0.5, 0.5, 1 };
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
	qaws_exact_surface_prepare(NULL, s, &e, NULL);
	TEST_ASSERT(qaws_exact_surface_line_hits(e, a, b, hits, 4, &count) == QAWS_STATUS_CERTIFICATION_FAILED, "a line lying on the surface is refused");
	TEST_ASSERT(qaws_exact_surface_line_hits(e, c, dd, hits, 4, &count) == QAWS_STATUS_OK && count == 1 && hits[0].kind == QAWS_EXACT_HIT_POINT &&
		hits[0].u_lo == 0.5 && hits[0].v_lo == 0.5 && hits[0].t_lo == 0.5, "a vertical ray through the center: the exact hit (1/2, 1/2), t = 1/2");
	qaws_exact_surface_destroy(e);
	qaws_surface_destroy(s);
}

int test_74_exact_surface_hits_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 74: Certified line / surface intersections\n");
	test_reference();
	test_overlap();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
