/*
 * Test 144: evaluation consistency fixes
 *
 *   - the rational Bezier core (core/qaws_rational_bezier_core.h, compiled
 *     as C): D2 asked without D1, and D2 of weighted degree-1 curves,
 *     against the runtime rational Bezier curve, 2D and 3D;
 *   - the surface core takes the runtime's flag word: the
 *     QAWS_SURFACE_CORE_* values are the QAWS_SURFACE_EVAL_* ones, and the
 *     Bezier surface core asked for second partials or the normal alone
 *     gives the runtime's values;
 *   - inline curves: a copy (assignment, memcpy) evaluates through
 *     qaws_curve_inline_get_curve after the original is overwritten; the
 *     storage is aligned inside other structs;
 *   - the arc-length curve and the length-matched curve return derivatives
 *     with respect to their own parameter (unit D1; D2 and D3 against
 *     central differences).
 */

#include "test_common.h"
#include "qaws_inline.h"
#include "qaws_platform.h"
#include "qaws_core_types.h"
#include "core/qaws_rational_bezier_core.h"
#include "core/qaws_surface_core.h"
#include "core/qaws_bspline_basis_core.h"
#include "core/qaws_bspline_eval_core.h"
#include "core/qaws_nurbs_eval_core.h"
#include "core/qaws_cubic_poly_core.h"
#include <math.h>
#include <string.h>
#include <stddef.h>

#if QAWS_SCALAR_IS_FLOAT
#define EC_TOL 1e-4
#define EC_FD 2e-2
#else
#define EC_TOL 1e-11
#define EC_FD 2e-5
#endif

static double ec_err2(qaws_vec2 a, qaws_vec2 b)
{
	return (fabs((double)a.x - b.x) + fabs((double)a.y - b.y)) / (1 + fabs((double)b.x) + fabs((double)b.y));
}

static double ec_err3(qaws_vec3 a, qaws_vec3 b)
{
	return (fabs((double)a.x - b.x) + fabs((double)a.y - b.y) + fabs((double)a.z - b.z)) / (1 + fabs((double)b.x) + fabs((double)b.y) + fabs((double)b.z));
}

static void test_rational_core(void)
{
	/* degree 1 and degree 3, 2D and 3D, uneven weights */
	static double const xyz[4][3] = { { 0, 0, 0 }, { 1, 2, -1 }, { 3, 1, 0.5 }, { 4, -1, 2 } };
	static double const w[4] = { 1, 2.5, 0.4, 1.3 };
	double worst = 0, e;
	unsigned int deg, dim, i, c, k;
	char msg[200];
	for (dim = 2; dim <= 3; dim++)
		for (deg = 1; deg <= 3; deg += 2)
		{
			qaws_scalar cp[4 * 3], ws[4], hom[QAWS_CORE_MAX_POINTS * 4];
			qaws_rational_bezier_desc d;
			qaws_curve* rc = NULL;
			for (i = 0; i <= deg; i++)
			{
				ws[i] = (qaws_scalar)w[i];
				for (c = 0; c < dim; c++)
				{
					cp[i * dim + c] = (qaws_scalar)xyz[i][c];
					hom[i * (dim + 1) + c] = (qaws_scalar)(xyz[i][c] * w[i]);
				}
				hom[i * (dim + 1) + dim] = ws[i];
			}
			memset(&d, 0, sizeof d);
			d.dimension = dim == 2 ? QAWS_DIMENSION_2D : QAWS_DIMENSION_3D;
			d.degree = deg;
			d.control_points = cp;
			d.control_point_count = deg + 1;
			d.weights = ws;
			d.weight_count = deg + 1;
			if (qaws_curve_create_rational_bezier(&d, &rc) != QAWS_STATUS_OK)
			{
				TEST_ASSERT(0, "rational Bezier created");
				continue;
			}
			for (k = 0; k <= 10; k++)
			{
				qaws_scalar t = (qaws_scalar)(k / 10.0);
				unsigned int f;
				for (f = 0; f < 2; f++)
				{
					/* D2 alone, then D1 and D2 */
					int flags = f ? 0x7 : 0x5;
					if (dim == 2)
					{
						qaws_eval_result_2d r;
						qaws_eval_2d core = qaws_rational_bezier_eval_2d(hom, (int)deg, t, flags);
						qaws_curve_evaluate_2d(rc, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &r);
						e = ec_err2(core.d2, r.d2) + ec_err2(core.position, r.position) + (f ? ec_err2(core.d1, r.d1) : 0);
					}
					else
					{
						qaws_eval_result_3d r;
						qaws_eval_3d core = qaws_rational_bezier_eval_3d(hom, (int)deg, t, flags);
						qaws_curve_evaluate_3d(rc, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &r);
						e = ec_err3(core.d2, r.d2) + ec_err3(core.position, r.position) + (f ? ec_err3(core.d1, r.d1) : 0);
					}
					worst = e > worst ? e : worst;
				}
			}
			qaws_curve_destroy(rc);
		}
	sprintf(msg, "rational Bezier core against the runtime (D2 with and without D1, degrees 1 and 3, 2D and 3D): %.2g", worst);
	TEST_ASSERT(worst <= EC_TOL * 10, msg);
	printf("  %s\n", msg);
}

static void test_surface_core_flags(void)
{
	qaws_vec3 grid[16];
	qaws_scalar flat[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	double worst = 0;
	unsigned int i, j, k;
	char msg[200];
	TEST_ASSERT(QAWS_SURFACE_CORE_POSITION == QAWS_SURFACE_EVAL_POSITION && QAWS_SURFACE_CORE_DU == QAWS_SURFACE_EVAL_DU
		&& QAWS_SURFACE_CORE_DV == QAWS_SURFACE_EVAL_DV && QAWS_SURFACE_CORE_NORMAL == QAWS_SURFACE_EVAL_NORMAL
		&& QAWS_SURFACE_CORE_DUU == QAWS_SURFACE_EVAL_DUU && QAWS_SURFACE_CORE_DUV == QAWS_SURFACE_EVAL_DUV
		&& QAWS_SURFACE_CORE_DVV == QAWS_SURFACE_EVAL_DVV, "surface core flags: the runtime's values");
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			grid[i * 4 + j].x = (qaws_scalar)i;
			grid[i * 4 + j].y = (qaws_scalar)j;
			grid[i * 4 + j].z = (qaws_scalar)(0.5 * sin(1.1 * i + 0.3) * cos(0.8 * j));
			flat[(i * 4 + j) * 3 + 0] = grid[i * 4 + j].x;
			flat[(i * 4 + j) * 3 + 1] = grid[i * 4 + j].y;
			flat[(i * 4 + j) * 3 + 2] = grid[i * 4 + j].z;
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = grid;
	d.u_point_count = 4;
	d.v_point_count = 4;
	if (qaws_surface_create_bezier(&d, &s) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "Bezier surface created");
		return;
	}
	for (k = 0; k < 25; k++)
	{
		qaws_scalar u = (qaws_scalar)((k % 5) / 4.0), v = (qaws_scalar)((k / 5) / 4.0);
		qaws_surface_eval_result r;
		qaws_surface_eval_3d c;
		double e;
		qaws_surface_evaluate(s, u, v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL
			| QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV, &r);
		c = qaws_surface_bezier_eval_3d(flat, 3, 3, 4, u, v, QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DVV);
		e = ec_err3(c.duu, r.duu) + ec_err3(c.duv, r.duv) + ec_err3(c.dvv, r.dvv);
		c = qaws_surface_bezier_eval_3d(flat, 3, 3, 4, u, v, QAWS_SURFACE_EVAL_NORMAL);
		e += ec_err3(c.normal, r.normal);
		worst = e > worst ? e : worst;
	}
	sprintf(msg, "Bezier surface core with the runtime's flags (second partials, the normal alone): %.2g", worst);
	TEST_ASSERT(worst <= EC_TOL * 10, msg);
	qaws_surface_destroy(s);
}

typedef struct ec_holder
{
	char tag;
	qaws_curve_inline curve;
} ec_holder;

static void test_inline_copy(void)
{
	qaws_vec3 pts[4] = { { 0, 0, 0 }, { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
	qaws_scalar coef[3 * 3] = { 1, 2, 3, 0.5, -1, 2, 0.25, 0.1, -0.3 };
	qaws_bezier_desc bd;
	qaws_polynomial_desc pd;
	qaws_curve_inline a, b, p, q;
	ec_holder h[2];
	unsigned int bad = 0, k;
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_3D;
	bd.degree = 3;
	bd.control_points = pts;
	bd.control_point_count = 4;
	memset(&pd, 0, sizeof pd);
	pd.dimension = QAWS_DIMENSION_3D;
	pd.degree = 2;
	pd.coefficients = coef;
	pd.coefficient_count = 3;
	TEST_ASSERT(qaws_curve_init_bezier_inline(&bd, &a) == QAWS_STATUS_OK && qaws_curve_init_polynomial_inline(&pd, &p) == QAWS_STATUS_OK,
		"inline curves initialised");
	b = a;
	memcpy(&q, &p, sizeof q);
	for (k = 0; k <= 8; k++)
	{
		qaws_eval_result_3d ra, rb, rp, rq;
		qaws_scalar t = (qaws_scalar)(k / 8.0);
		qaws_curve_evaluate_3d(qaws_curve_inline_get_curve(&a), t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &ra);
		qaws_curve_evaluate_3d(qaws_curve_inline_get_curve(&p), t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &rp);
		/* the originals overwritten: the copies must not read them */
		memset(&a, 0x5a, sizeof a);
		memset(&p, 0x5a, sizeof p);
		qaws_curve_evaluate_3d(qaws_curve_inline_get_curve(&b), t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &rb);
		qaws_curve_evaluate_3d(qaws_curve_inline_get_curve(&q), t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, &rq);
		bad += memcmp(&ra, &rb, sizeof ra) != 0 || memcmp(&rp, &rq, sizeof rp) != 0;
		a = b;
		p = q;
	}
	TEST_ASSERT(bad == 0, "inline curves copied by assignment and memcpy evaluate as the originals after those are overwritten");
	TEST_ASSERT(offsetof(ec_holder, curve) % sizeof(double) == 0 && ((size_t)&h[1].curve) % sizeof(void*) == 0,
		"inline storage aligned inside other structs");
}

/* wrapper derivatives against central differences of the wrapper itself */
static double ec_wrapper_fd(qaws_curve const* c, double* unit)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double worst = 0, h = (double)(r.max_value - r.min_value) * EC_FD, len_err = 0;
	unsigned int k;
	for (k = 1; k < 10; k++)
	{
		qaws_scalar s = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)k / 10;
		qaws_eval_result_2d m, p, z;
		qaws_vec2 fd1, fd2, fd3;
		qaws_curve_evaluate_2d(c, s, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, &z);
		qaws_curve_evaluate_2d(c, s - (qaws_scalar)h, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &m);
		qaws_curve_evaluate_2d(c, s + (qaws_scalar)h, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &p);
		fd1.x = (qaws_scalar)(((double)p.position.x - m.position.x) / (2 * h));
		fd1.y = (qaws_scalar)(((double)p.position.y - m.position.y) / (2 * h));
		fd2.x = (qaws_scalar)(((double)p.d1.x - m.d1.x) / (2 * h));
		fd2.y = (qaws_scalar)(((double)p.d1.y - m.d1.y) / (2 * h));
		fd3.x = (qaws_scalar)(((double)p.d2.x - m.d2.x) / (2 * h));
		fd3.y = (qaws_scalar)(((double)p.d2.y - m.d2.y) / (2 * h));
		worst = ec_err2(z.d1, fd1) > worst ? ec_err2(z.d1, fd1) : worst;
		worst = ec_err2(z.d2, fd2) > worst ? ec_err2(z.d2, fd2) : worst;
		worst = ec_err2(z.d3, fd3) > worst ? ec_err2(z.d3, fd3) : worst;
		{
			double l = sqrt((double)z.d1.x * z.d1.x + (double)z.d1.y * z.d1.y);
			len_err = fabs(l - (unit ? *unit : l)) > len_err ? fabs(l - (unit ? *unit : l)) : len_err;
		}
	}
	if (unit)
		*unit = len_err;
	return worst;
}

static void test_arc_length_wrappers(void)
{
	/* a cubic Bezier with uneven speed, and a quadratic for the length match */
	qaws_vec2 pa[4] = { { 0, 0 }, { 0.2f, 1.5f }, { 2.5f, 1.8f }, { 3, 0 } };
	qaws_vec2 pb[3] = { { 0, 0 }, { 4, 3 }, { 5, -1 } };
	qaws_bezier_desc d;
	qaws_curve *a = NULL, *b = NULL, *ra = NULL, *mb = NULL;
	qaws_scalar la = 0, lb = 0;
	double ea, eb, ua = 1, ub;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = pa;
	d.control_point_count = 4;
	qaws_curve_create_bezier(&d, &a);
	d.degree = 2;
	d.control_points = pb;
	d.control_point_count = 3;
	qaws_curve_create_bezier(&d, &b);
	if (!a || !b || qaws_curve_reparameterize_arc_length(a, 2048, &ra) != QAWS_STATUS_OK || qaws_curve_match_arc_length(a, b, 2048, &mb) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "arc-length wrappers created");
		return;
	}
	qaws_curve_compute_arc_length(a, 0, 1, &la);
	qaws_curve_compute_arc_length(b, 0, 1, &lb);
	ea = ec_wrapper_fd(ra, &ua);
	ub = (double)lb / la;
	eb = ec_wrapper_fd(mb, &ub);
	sprintf(msg, "arc-length curve: |D1| - 1 within %.2g, D1 to D3 against differences %.2g; length-matched: |D1| - L_B/L_A within %.2g, "
		"derivatives %.2g", ua, ea, ub, eb);
	/* the wrappers interpolate an arc-length table: the differences see its resolution */
	TEST_ASSERT(ua <= 1e-3 && ub <= 1e-3 && ea <= 2e-2 && eb <= 2e-2, msg);
	printf("  %s\n", msg);
	qaws_curve_destroy(ra);
	qaws_curve_destroy(mb);
	qaws_curve_destroy(a);
	qaws_curve_destroy(b);
}

/* the B-spline, NURBS, cubic and B-spline surface cores (shader code, compiled as C) against the runtime objects: the CPU end of
   the shader validation (examples/eval_gpu_parity.c runs the same cores on a GPU against these) */
static void test_spline_cores(void)
{
	static double const kn[11] = { 0, 0, 0, 0, 0.7, 1.5, 1.5, 3, 3, 3, 3 };   /* n = 7, cubic, a double knot */
	qaws_scalar cp[7 * 3], ws[7], k[11], full[QAWS_CORE_MAX_POINTS * 2];
	qaws_bspline_desc bd;
	qaws_nurbs_desc nd;
	qaws_polynomial_desc pd;
	qaws_curve *bs = NULL, *nu = NULL, *po = NULL;
	qaws_scalar coef[4 * 3];
	double eb = 0, en = 0, ec = 0, es = 0;
	unsigned int i, j, s;
	char msg[300];
	memset(full, 0, sizeof full);
	for (i = 0; i < 7; i++)
	{
		cp[i * 3] = (qaws_scalar)cos(0.9 * i);
		cp[i * 3 + 1] = (qaws_scalar)sin(1.3 * i);
		cp[i * 3 + 2] = (qaws_scalar)(0.2 * i);
		ws[i] = (qaws_scalar)(0.5 + 0.3 * (i % 3));
	}
	for (i = 0; i < 11; i++)
		full[i] = k[i] = (qaws_scalar)kn[i];
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_3D;
	bd.degree = 3;
	bd.control_points = cp;
	bd.control_point_count = 7;
	bd.knots = k;
	bd.knot_count = 11;
	memset(&nd, 0, sizeof nd);
	nd.dimension = QAWS_DIMENSION_3D;
	nd.degree = 3;
	nd.control_points = cp;
	nd.control_point_count = 7;
	nd.knots = k;
	nd.knot_count = 11;
	nd.weights = ws;
	nd.weight_count = 7;
	for (i = 0; i < 12; i++)
		coef[i] = (qaws_scalar)(0.3 * sin(1.7 * i + 0.2));
	memset(&pd, 0, sizeof pd);
	pd.dimension = QAWS_DIMENSION_3D;
	pd.degree = 3;
	pd.coefficients = coef;
	pd.coefficient_count = 4;
	{
		qaws_status s1 = qaws_curve_create_bspline(&bd, &bs), s2 = qaws_curve_create_nurbs(&nd, &nu), s3 = qaws_curve_create_polynomial(&pd, &po);
		sprintf(msg, "spline curves created: B-spline %s, NURBS %s, polynomial %s", qaws_status_to_string(s1), qaws_status_to_string(s2),
			qaws_status_to_string(s3));
	}
	if (!bs || !nu || !po)
	{
		TEST_ASSERT(0, msg);
		qaws_curve_destroy(bs);
		qaws_curve_destroy(nu);
		qaws_curve_destroy(po);
		return;
	}
	for (s = 0; s <= 30; s++)
	{
		qaws_scalar t = (qaws_scalar)(3.0 * s / 30);
		qaws_scalar lcp[QAWS_CORE_MAX_POINTS * 3], lk[QAWS_CORE_MAX_POINTS * 2], lw[QAWS_CORE_MAX_POINTS];
		qaws_eval_result_3d r;
		qaws_eval_3d c;
		qaws_vec3 a, b, cc, d;
		int span = qaws_find_span(full, 3, 7, t);
		for (j = 0; j < 8; j++)
			lk[j] = k[span - 3 + (int)j];
		for (j = 0; j <= 3; j++)
		{
			lcp[j * 3] = cp[(span - 3 + (int)j) * 3];
			lcp[j * 3 + 1] = cp[(span - 3 + (int)j) * 3 + 1];
			lcp[j * 3 + 2] = cp[(span - 3 + (int)j) * 3 + 2];
			lw[j] = ws[span - 3 + (int)j];
		}
		c = qaws_bspline_eval_3d(lcp, lk, 3, 3, t, 7);
		qaws_curve_evaluate_3d(bs, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &r);
		eb = fmax(eb, ec_err3(c.position, r.position) + ec_err3(c.d1, r.d1) + ec_err3(c.d2, r.d2));
		c = qaws_nurbs_eval_3d(lcp, lk, lw, 3, 3, t, 7);
		qaws_curve_evaluate_3d(nu, t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &r);
		en = fmax(en, ec_err3(c.position, r.position) + ec_err3(c.d1, r.d1) + ec_err3(c.d2, r.d2));
		/* the cubic core's a t^3 + b t^2 + c t + d against the monomial curve sum coef[i] t^i on [0, 1] */
		if (s <= 10)
		{
			qaws_scalar u = (qaws_scalar)(s / 10.0);
			a.x = coef[9]; a.y = coef[10]; a.z = coef[11];
			b.x = coef[6]; b.y = coef[7]; b.z = coef[8];
			cc.x = coef[3]; cc.y = coef[4]; cc.z = coef[5];
			d.x = coef[0]; d.y = coef[1]; d.z = coef[2];
			c = qaws_cubic_eval_3d(a, b, cc, d, u, 7);
			qaws_curve_evaluate_3d(po, u, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &r);
			ec = fmax(ec, ec_err3(c.position, r.position) + ec_err3(c.d1, r.d1) + ec_err3(c.d2, r.d2));
		}
	}
	qaws_curve_destroy(bs);
	qaws_curve_destroy(nu);
	qaws_curve_destroy(po);
	/* the B-spline surface core: 6 x 5 points, degrees 3 and 2 */
	{
		static double const uk[10] = { 0, 0, 0, 0, 0.4, 0.4, 1, 1, 1, 1 }, vk[8] = { 0, 0, 0, 0.3, 0.75, 1, 1, 1 };
		qaws_vec3 grid[30];
		qaws_scalar flat[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3], fu[QAWS_CORE_MAX_POINTS * 2], fv[QAWS_CORE_MAX_POINTS * 2];
		qaws_scalar uks[10], vks[8];
		qaws_surface_bspline_desc sd;
		qaws_surface* sf = NULL;
		memset(fu, 0, sizeof fu);
		memset(fv, 0, sizeof fv);
		for (i = 0; i < 30; i++)
		{
			grid[i].x = (qaws_scalar)((i / 5) / 5.0);
			grid[i].y = (qaws_scalar)((i % 5) / 4.0);
			grid[i].z = (qaws_scalar)(0.3 * sin(1.3 * (i / 5)) * cos(0.9 * (i % 5)));
			flat[i * 3] = grid[i].x;
			flat[i * 3 + 1] = grid[i].y;
			flat[i * 3 + 2] = grid[i].z;
		}
		for (i = 0; i < 10; i++)
			fu[i] = uks[i] = (qaws_scalar)uk[i];
		for (i = 0; i < 8; i++)
			fv[i] = vks[i] = (qaws_scalar)vk[i];
		memset(&sd, 0, sizeof sd);
		sd.u_degree = 3;
		sd.v_degree = 2;
		sd.control_points = grid;
		sd.u_point_count = 6;
		sd.v_point_count = 5;
		sd.u_knots = uks;
		sd.u_knot_count = 10;
		sd.v_knots = vks;
		sd.v_knot_count = 8;
		if (qaws_surface_create_bspline(&sd, &sf) == QAWS_STATUS_OK)
		{
			for (s = 0; s < 49; s++)
			{
				qaws_scalar u = (qaws_scalar)((s % 7) / 6.0), v = (qaws_scalar)((s / 7) / 6.0);
				qaws_surface_eval_result r;
				qaws_surface_eval_3d c;
				qaws_surface_evaluate(sf, u, v, 127, &r);
				c = qaws_surface_bspline_eval_3d(flat, fu, fv, 3, 2, 6, 5, qaws_find_span(fu, 3, 6, u), qaws_find_span(fv, 2, 5, v), u, v, 127);
				es = fmax(es, ec_err3(c.position, r.position) + ec_err3(c.du, r.du) + ec_err3(c.dv, r.dv) + ec_err3(c.normal, r.normal)
					+ ec_err3(c.duu, r.duu) + ec_err3(c.duv, r.duv) + ec_err3(c.dvv, r.dvv));
			}
			qaws_surface_destroy(sf);
		}
		else
			es = 1e30;
	}
	sprintf(msg, "cores as C against the runtime: B-spline %.2g, NURBS %.2g, cubic %.2g, B-spline surface %.2g", eb, en, ec, es);
	TEST_ASSERT(eb <= EC_TOL * 10 && en <= EC_TOL * 10 && ec <= EC_TOL * 10 && es <= EC_TOL * 10, msg);
	printf("  %s\n", msg);
}

int test_144_eval_consistency_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 144: evaluation consistency fixes\n");
	test_rational_core();
	test_surface_core_flags();
	test_inline_copy();
	test_arc_length_wrappers();
	test_spline_cores();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
