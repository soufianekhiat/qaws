/*
 * Test 112: tangent-point (repulsive) energies of curves and meshes
 * (qaws_repulsive.h; the models in docs/repulsive_plan.md)
 *
 *   - against Mathematica (reference/112_repulsive.wls: the discrete sums
 *     differentiated symbolically, 30 digits): the energy, its first and
 *     second derivatives along a direction, the gradient and the Hessian
 *     times the direction, for a closed 3D curve, an open 2D curve
 *     (alpha = 2), a network (a junction, a triangle, an open polyline), an
 *     octahedron and three separate patches (p = 5.5)
 *   - invariance under translation and rotation; the scaling degree
 *     (2 + alpha - beta for curves, 4 - p for meshes), gradients one less
 *   - the adjoint identity <grad, v> = dE[v], the Hessian's symmetry
 *     <u, H v> = <v, H u>, and the HVP against differences of gradients
 *   - a straight curve (|T x d| = 0 everywhere): zero energy, gradient and
 *     Hessian for alpha = 3; the HVP against differences for alpha = 2
 *   - refusals: exponents out of range, bad indices, a zero-length edge, a
 *     zero-area face, coinciding points; outputs untouched on failure
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "reference/112_repulsive.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RT_F32 (sizeof(qaws_scalar) == 4)
#define RT_MAXV 64

typedef struct rt_case
{
	int mesh;
	unsigned int dim, nv, ne;
	double const *x, *dir, *param, *value, *grad, *hv;
	unsigned int const* idx;
} rt_case;

/* a prepared object for a case, its two kinds behind one interface */
typedef struct rt_obj
{
	qaws_repulsive_curves* c;
	qaws_repulsive_mesh* m;
} rt_obj;

static qaws_status rt_make(rt_case const* k, double const* param, rt_obj* o)
{
	memset(o, 0, sizeof(*o));
	if (k->mesh)
	{
		qaws_repulsive_mesh_desc d;
		memset(&d, 0, sizeof(d));
		d.vertex_count = k->nv; d.triangles = k->idx; d.triangle_count = k->ne; d.p = (qaws_scalar)param[0];
		return qaws_repulsive_mesh_create(&d, &o->m);
	}
	else
	{
		qaws_repulsive_curves_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = k->dim; d.vertex_count = k->nv; d.edges = k->idx; d.edge_count = k->ne;
		d.alpha = (qaws_scalar)param[0]; d.beta = (qaws_scalar)param[1];
		return qaws_repulsive_curves_create(&d, &o->c);
	}
}

static void rt_free(rt_obj* o)
{
	qaws_repulsive_curves_destroy(o->c);
	qaws_repulsive_mesh_destroy(o->m);
}

static qaws_status rt_eval(rt_obj const* o, qaws_scalar const* x, qaws_scalar const* v, double* e, double* t, double* t2)
{
	qaws_scalar a = 0, b = 0, c = 0;
	qaws_status st = o->m ? qaws_repulsive_mesh_eval(o->m, x, v, &a, &b, &c) : qaws_repulsive_curves_eval(o->c, x, v, &a, &b, &c);
	if (e) *e = (double)a;
	if (t) *t = (double)b;
	if (t2) *t2 = (double)c;
	return st;
}

/* g = gradient (cleared first) */
static qaws_status rt_grad(rt_obj const* o, qaws_scalar const* x, unsigned int n, qaws_scalar* g)
{
	memset(g, 0, sizeof(qaws_scalar) * n);
	return o->m ? qaws_repulsive_mesh_gradient(o->m, x, 1, g, NULL) : qaws_repulsive_curves_gradient(o->c, x, 1, g, NULL);
}

static qaws_status rt_hvp(rt_obj const* o, qaws_scalar const* x, qaws_scalar const* v, unsigned int n, qaws_scalar* hv)
{
	memset(hv, 0, sizeof(qaws_scalar) * n);
	return o->m ? qaws_repulsive_mesh_hvp(o->m, x, v, hv) : qaws_repulsive_curves_hvp(o->c, x, v, hv);
}

static double rt_rel(double a, double b, double scale)
{
	return fabs(a - b) / (scale > 0 ? scale : 1);
}

static double rt_vec_err(qaws_scalar const* a, double const* b, unsigned int n)
{
	double m = 0, e = 0;
	unsigned int i;
	for (i = 0; i < n; i++) m = fmax(m, fabs(b[i]));
	for (i = 0; i < n; i++) e = fmax(e, fabs((double)a[i] - b[i]));
	return e / (m > 0 ? m : 1);
}

static void test_fixture(rt_case const* k, char const* name)
{
	unsigned int n = k->nv * k->dim, i;
	qaws_scalar x[3 * RT_MAXV], v[3 * RT_MAXV], g[3 * RT_MAXV], hv[3 * RT_MAXV];
	double e = 0, t = 0, t2 = 0, tol = RT_F32 ? 1e-5 : 1e-12, eg, eh;
	rt_obj o;
	qaws_status s0, s1, s2, s3;
	char msg[300];
	for (i = 0; i < n; i++) { x[i] = (qaws_scalar)k->x[i]; v[i] = (qaws_scalar)k->dir[i]; }
	s0 = rt_make(k, k->param, &o);
	s1 = s0 ? s0 : rt_eval(&o, x, v, &e, &t, &t2);
	s2 = s0 ? s0 : rt_grad(&o, x, n, g);
	s3 = s0 ? s0 : rt_hvp(&o, x, v, n, hv);
	eg = rt_vec_err(g, k->grad, n);
	eh = rt_vec_err(hv, k->hv, n);
	sprintf(msg, "%s against Mathematica: E %.2g, dE %.2g, d2E %.2g, gradient %.2g, HVP %.2g (statuses %d %d %d %d)", name,
		rt_rel(e, k->value[0], fabs(k->value[0])), rt_rel(t, k->value[1], fabs(k->value[1])), rt_rel(t2, k->value[2], fabs(k->value[2])),
		eg, eh, (int)s0, (int)s1, (int)s2, (int)s3);
	TEST_ASSERT(!s0 && !s1 && !s2 && !s3 && rt_rel(e, k->value[0], fabs(k->value[0])) <= tol && rt_rel(t, k->value[1], fabs(k->value[1])) <= tol &&
		rt_rel(t2, k->value[2], fabs(k->value[2])) <= tol && eg <= tol && eh <= tol, msg);
	rt_free(&o);
}

/* the rigid motion and scaling checks, the adjoint identity, symmetry and
   the HVP against differences of gradients, on one case */
static void test_identities(rt_case const* k, char const* name)
{
	unsigned int n = k->nv * k->dim, i, c;
	qaws_scalar x[3 * RT_MAXV], y[3 * RT_MAXV], v[3 * RT_MAXV], u[3 * RT_MAXV], g[3 * RT_MAXV], g2[3 * RT_MAXV], hv[3 * RT_MAXV], hu[3 * RT_MAXV];
	double e0 = 0, e1 = 0, e2 = 0, e3 = 0, t = 0, gv = 0, uhv = 0, vhu = 0, deg, cs = 1.7, best = 1e300, h, gm = 0, gsm = 0, tol = RT_F32 ? 1e-6 : 1e-13;
	double R[3][3] = { { 0.36, 0.48, -0.8 }, { -0.8, 0.6, 0 }, { 0.48, 0.64, 0.6 } };
	rt_obj o;
	char msg[300];
	if (k->dim == 2) { R[0][0] = 0.6; R[0][1] = -0.8; R[1][0] = 0.8; R[1][1] = 0.6; }
	deg = k->mesh ? 4 - k->param[0] : 2 + k->param[0] - k->param[1];
	for (i = 0; i < n; i++) { x[i] = (qaws_scalar)k->x[i]; v[i] = (qaws_scalar)k->dir[i]; u[i] = (qaws_scalar)sin(1.3 * i + 0.4); }
	if (rt_make(k, k->param, &o) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "identities: the object");
		return;
	}
	rt_eval(&o, x, v, &e0, &t, NULL);
	/* translation, rotation, scaling */
	for (i = 0; i < k->nv; i++) for (c = 0; c < k->dim; c++) y[k->dim * i + c] = (qaws_scalar)(k->x[k->dim * i + c] + 0.3 + 0.2 * c);
	rt_eval(&o, y, NULL, &e1, NULL, NULL);
	for (i = 0; i < k->nv; i++)
		for (c = 0; c < k->dim; c++)
		{
			double s = 0;
			unsigned int q;
			for (q = 0; q < k->dim; q++) s += R[c][q] * k->x[k->dim * i + q];
			y[k->dim * i + c] = (qaws_scalar)s;
		}
	rt_eval(&o, y, NULL, &e2, NULL, NULL);
	for (i = 0; i < n; i++) y[i] = (qaws_scalar)(cs * k->x[i]);
	rt_eval(&o, y, NULL, &e3, NULL, NULL);
	rt_grad(&o, x, n, g);
	rt_grad(&o, y, n, g2);
	for (i = 0; i < n; i++)
	{
		gm = fmax(gm, fabs((double)g[i]));
		gsm = fmax(gsm, fabs((double)g2[i] - pow(cs, deg - 1) * (double)g[i]));
	}
	sprintf(msg, "%s: translated %.2g, rotated %.2g, scaled by 1.7: E %.2g, gradient %.2g off degree %g", name,
		rt_rel(e1, e0, e0), rt_rel(e2, e0, e0), rt_rel(e3, pow(cs, deg) * e0, pow(cs, deg) * e0), gsm / (pow(cs, deg - 1) * gm), deg);
	TEST_ASSERT(rt_rel(e1, e0, e0) <= 10 * tol && rt_rel(e2, e0, e0) <= 10 * tol && rt_rel(e3, pow(cs, deg) * e0, pow(cs, deg) * e0) <= 10 * tol &&
		gsm <= 10 * tol * pow(cs, deg - 1) * gm, msg);
	/* adjoint identity and symmetry */
	rt_hvp(&o, x, v, n, hv);
	rt_hvp(&o, x, u, n, hu);
	for (i = 0; i < n; i++) { gv += (double)g[i] * v[i]; uhv += (double)u[i] * hv[i]; vhu += (double)v[i] * hu[i]; }
	sprintf(msg, "%s: <grad, v> against dE[v] %.2g, <u, Hv> against <v, Hu> %.2g", name, rt_rel(gv, t, fabs(t)), rt_rel(uhv, vhu, fabs(uhv)));
	TEST_ASSERT(rt_rel(gv, t, fabs(t)) <= 10 * tol && rt_rel(uhv, vhu, fabs(uhv)) <= 10 * tol, msg);
	/* differences of gradients along v (central, swept over h; the best) */
	if (!RT_F32)
	{
		for (h = 1e-3; h >= 1e-7; h *= 0.1)
		{
			double err = 0, m = 0;
			for (i = 0; i < n; i++) y[i] = x[i] + h * v[i];
			rt_grad(&o, y, n, g);
			for (i = 0; i < n; i++) y[i] = x[i] - h * v[i];
			rt_grad(&o, y, n, g2);
			for (i = 0; i < n; i++)
			{
				err = fmax(err, fabs(((double)g[i] - g2[i]) / (2 * h) - hv[i]));
				m = fmax(m, fabs((double)hv[i]));
			}
			best = fmin(best, err / m);
		}
		sprintf(msg, "%s: HVP against central differences of gradients, best %.2g", name, best);
		TEST_ASSERT(best <= 1e-7, msg);
	}
	rt_free(&o);
}

/* a straight open curve: every point on every other edge's line */
static void test_straight(void)
{
	static unsigned int const e[5][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 } };
	qaws_scalar x[18], v[18], g[18], hv[18], g1[18], g2[18], y[18];
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves *c3 = NULL, *c2 = NULL;
	qaws_scalar E = 1;
	double gm = 0, hm = 0, err = 0, m = 0, h = 1e-4;
	unsigned int i;
	char msg[300];
	for (i = 0; i < 6; i++)
	{
		double s = 0.3 * i + 0.02 * i * i;
		x[3 * i] = (qaws_scalar)s; x[3 * i + 1] = 0; x[3 * i + 2] = 0;   /* (on an axis: exactly collinear) */
		v[3 * i] = (qaws_scalar)sin(i + 1.0); v[3 * i + 1] = (qaws_scalar)cos(2.0 * i); v[3 * i + 2] = (qaws_scalar)sin(3.0 * i + 1);
	}
	memset(&d, 0, sizeof(d));
	d.dimension = 3; d.vertex_count = 6; d.edges = &e[0][0]; d.edge_count = 5;
	qaws_repulsive_curves_create(&d, &c3);
	d.alpha = 2; d.beta = (qaws_scalar)4.5;
	qaws_repulsive_curves_create(&d, &c2);
	memset(g, 0, sizeof(g)); memset(hv, 0, sizeof(hv));
	qaws_repulsive_curves_eval(c3, x, NULL, &E, NULL, NULL);
	qaws_repulsive_curves_gradient(c3, x, 1, g, NULL);
	qaws_repulsive_curves_hvp(c3, x, v, hv);
	for (i = 0; i < 18; i++) { gm = fmax(gm, fabs((double)g[i])); hm = fmax(hm, fabs((double)hv[i])); }
	sprintf(msg, "straight curve, alpha = 3: E %g, gradient %g, HVP %g (all 0)", (double)E, gm, hm);
	TEST_ASSERT(c3 && E == 0 && gm == 0 && hm == 0, msg);
	/* alpha = 2: |T x d|^2, a nonzero Hessian at the line */
	memset(hv, 0, sizeof(hv));
	qaws_repulsive_curves_hvp(c2, x, v, hv);
	for (i = 0; i < 18; i++) y[i] = x[i] + (qaws_scalar)h * v[i];
	memset(g1, 0, sizeof(g1)); qaws_repulsive_curves_gradient(c2, y, 1, g1, NULL);
	for (i = 0; i < 18; i++) y[i] = x[i] - (qaws_scalar)h * v[i];
	memset(g2, 0, sizeof(g2)); qaws_repulsive_curves_gradient(c2, y, 1, g2, NULL);
	for (i = 0; i < 18; i++) { err = fmax(err, fabs(((double)g1[i] - g2[i]) / (2 * h) - hv[i])); m = fmax(m, fabs((double)hv[i])); }
	sprintf(msg, "straight curve, alpha = 2: HVP %.3g against central differences, %.2g off", m, err / m);
	TEST_ASSERT(c2 && m > 0 && err / m <= (RT_F32 ? 2e-2 : 1e-6), msg);
	qaws_repulsive_curves_destroy(c3); qaws_repulsive_curves_destroy(c2);
}

static void test_refusals(void)
{
	static unsigned int const e[3][2] = { { 0, 1 }, { 1, 2 }, { 3, 4 } }, bad[1][2] = { { 0, 7 } };
	static unsigned int const t[2][3] = { { 0, 1, 2 }, { 3, 4, 5 } };
	qaws_scalar x[15] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 2, 0, 0, 2, 1, 0 };
	qaws_scalar g[18], keep[18];
	qaws_repulsive_curves_desc d;
	qaws_repulsive_mesh_desc md;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_mesh* m = NULL;
	qaws_status s1, s2, s3, s4, s5, s6, s7;
	unsigned int i, same = 1;
	char msg[300];
	memset(&d, 0, sizeof(d));
	d.dimension = 3; d.vertex_count = 5; d.edges = &e[0][0]; d.edge_count = 3;
	d.alpha = (qaws_scalar)1.5; d.beta = 4;
	s1 = qaws_repulsive_curves_create(&d, &c);                     /* alpha < 2 */
	d.alpha = 3; d.beta = 8;
	s2 = qaws_repulsive_curves_create(&d, &c);                     /* beta > 2 alpha + 1 */
	d.alpha = 0; d.beta = 0; d.edges = &bad[0][0]; d.edge_count = 1;
	s3 = qaws_repulsive_curves_create(&d, &c);                     /* an index out of range */
	d.edges = &e[0][0]; d.edge_count = 3;
	qaws_repulsive_curves_create(&d, &c);
	/* edge (3, 4) on top of edge (0, 1): coinciding points */
	x[9] = 0; x[10] = 0; x[11] = 0; x[12] = 1; x[13] = 0; x[14] = 0;
	for (i = 0; i < 18; i++) g[i] = keep[i] = (qaws_scalar)(0.5 * i);
	s4 = qaws_repulsive_curves_gradient(c, x, 1, g, NULL);
	for (i = 0; i < 15; i++) same &= g[i] == keep[i];
	/* a zero-length edge */
	x[12] = 0;
	s5 = qaws_repulsive_curves_eval(c, x, NULL, NULL, NULL, NULL);
	qaws_repulsive_curves_destroy(c);
	memset(&md, 0, sizeof(md));
	md.vertex_count = 6; md.triangles = &t[0][0]; md.triangle_count = 2; md.p = 4;
	s6 = qaws_repulsive_mesh_create(&md, &m);                      /* p <= 4 */
	md.p = 0;
	qaws_repulsive_mesh_create(&md, &m);
	{
		/* a zero-area face */
		qaws_scalar y[18] = { 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1 };
		s7 = qaws_repulsive_mesh_eval(m, y, NULL, NULL, NULL, NULL);
	}
	sprintf(msg, "refusals: alpha 1.5 %d, beta 8 %d, index %d, coinciding %d (gradient untouched %u), zero edge %d, p 4 %d, zero area %d",
		(int)s1, (int)s2, (int)s3, (int)s4, same, (int)s5, (int)s6, (int)s7);
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT &&
		s4 == QAWS_STATUS_NUMERICAL_FAILURE && same && s5 == QAWS_STATUS_DEGENERATE_CURVE && s6 == QAWS_STATUS_INVALID_ARGUMENT &&
		s7 == QAWS_STATUS_DEGENERATE_CURVE, msg);
	qaws_repulsive_mesh_destroy(m);
}

#define RT_CASE(tag, mesh_, dim_, k) { mesh_, dim_, sizeof(g112_##tag##_x) / sizeof(double) / (dim_), sizeof(g112_##tag##_idx) / sizeof(unsigned int) / (k), \
	g112_##tag##_x, g112_##tag##_dir, g112_##tag##_param, g112_##tag##_value, g112_##tag##_grad, g112_##tag##_hv, g112_##tag##_idx }

int test_112_repulsive_main(void)
{
	static rt_case const cases[5] = { RT_CASE(c1, 0, 3, 2), RT_CASE(c2, 0, 2, 2), RT_CASE(c3, 0, 3, 2), RT_CASE(s1, 1, 3, 3), RT_CASE(s2, 1, 3, 3) };
	static char const* const names[5] = { "closed 3D curve", "open 2D curve (alpha 2)", "network", "octahedron", "three patches (p 5.5)" };
	unsigned int i;
	g_pass = 0;
	g_fail = 0;
	printf("Test 112: tangent-point energies of curves and meshes\n");
	for (i = 0; i < 5; i++) test_fixture(&cases[i], names[i]);
	for (i = 0; i < 5; i++) test_identities(&cases[i], names[i]);
	test_straight();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
