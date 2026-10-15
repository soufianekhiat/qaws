/*
 * Test 117: curves on a parametric surface (qaws_repulsive.h, R5)
 *
 *   - the energy in the vertices' parameters, E(S(uv)), its first and
 *     second directional derivatives, gradient and Hessian-vector product,
 *     against Mathematica (reference/117_repulsive_surface.wls: the curves'
 *     energy composed with a bicubic Bezier patch, differentiated
 *     symbolically): a closed hexagon; an open polyline and a segment
 *   - on a torus (a surface of revolution of a circle, both parameters
 *     periodic): the energy unchanged by a seam crossing (uv shifted by a
 *     period), the uv gradient against central differences across the seam
 *   - solves on the torus: the energy falls, the total length held, every
 *     vertex on the torus; collision-checked, two loops never touch; pinned
 *     vertices keep their parameters
 *   - the domain: on the patch (not periodic) every vertex stays in [0, 1]^2
 *     while the energy pushes outward; a vertex outside the range, in a
 *     trimmed patch's hole, or at a sphere's pole is refused
 */

#include "test_common.h"
#include "qaws_repulsive.h"
#include "qaws_surface.h"
#include "qaws_surface_bezier.h"
#include "qaws_surface_revolution.h"
#include "qaws_surface_trim.h"
#include "qaws_arc.h"
#include "qaws_bezier.h"
#include "qaws_curve.h"
#include "reference/117_repulsive_surface.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define RS_F32 (sizeof(qaws_scalar) == 4)
#define RS_PI 3.14159265358979323846

static qaws_repulsive_curves* rs_curves(unsigned int nv, unsigned int const* e, unsigned int ne)
{
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = 3; d.vertex_count = nv; d.edges = e; d.edge_count = ne;
	qaws_repulsive_curves_create(&d, &c);
	return c;
}

static qaws_surface* rs_patch(void)
{
	qaws_vec3 p[16];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	unsigned int i;
	for (i = 0; i < 16; i++) { p[i].x = (qaws_scalar)g117_patch[3 * i]; p[i].y = (qaws_scalar)g117_patch[3 * i + 1]; p[i].z = (qaws_scalar)g117_patch[3 * i + 2]; }
	memset(&d, 0, sizeof(d));
	d.u_degree = 3; d.v_degree = 3; d.control_points = p; d.u_point_count = 4; d.v_point_count = 4;
	qaws_surface_create_bezier(&d, &s);
	return s;
}

/* a torus about z: the circle of radius 0.4 at distance 1, revolved */
static qaws_surface* rs_torus(qaws_curve** profile)
{
	qaws_arc_segment a;
	qaws_arc_desc ad;
	qaws_surface_revolution_desc rd;
	qaws_surface* s = NULL;
	memset(&a, 0, sizeof(a));
	a.center[0] = 1; a.radius = (qaws_scalar)0.4; a.angle_start = 0; a.angle_end = (qaws_scalar)(2 * RS_PI);
	memset(&ad, 0, sizeof(ad));
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &a; ad.segment_count = 1;
	*profile = NULL;
	qaws_curve_create_arc(&ad, profile);
	memset(&rd, 0, sizeof(rd));
	rd.profile = *profile; rd.axis_direction.z = 1; rd.angle = 0;
	qaws_surface_create_revolution(&rd, &s);
	return s;
}

static double rs_rel(double a, double b)
{
	return fabs(a - b) / fmax(fabs(b), 1e-300);
}

static double rs_relv(qaws_scalar const* a, double const* b, unsigned int n)
{
	double e = 0, m = 0;
	unsigned int i;
	for (i = 0; i < n; i++) { e = fmax(e, fabs((double)a[i] - b[i])); m = fmax(m, fabs(b[i])); }
	return e / m;
}

static void rs_fixture(qaws_surface const* patch, char const* tag, double const* uvd, double const* dird, unsigned int const* idx, unsigned int nv, unsigned int ne,
	double const* value, double const* grad, double const* hv)
{
	qaws_scalar uv[40], dir[40], g[40], h[40], v = 0, t = 0, t2 = 0;
	qaws_repulsive_surface_desc sd;
	qaws_repulsive_curves* c = rs_curves(nv, idx, ne);
	qaws_status s1, s2, s3;
	double ev, eg, eh, tol = RS_F32 ? 1e-4 : 1e-11;
	unsigned int i;
	char msg[300];
	for (i = 0; i < 2 * nv; i++) { uv[i] = (qaws_scalar)uvd[i]; dir[i] = (qaws_scalar)dird[i]; g[i] = h[i] = 0; }
	memset(&sd, 0, sizeof(sd));
	sd.surface = patch;
	s1 = qaws_repulsive_curves_on_surface_eval(c, &sd, uv, dir, &v, &t, &t2);
	s2 = qaws_repulsive_curves_on_surface_gradient(c, &sd, uv, 1, g, NULL);
	s3 = qaws_repulsive_curves_on_surface_hvp(c, &sd, uv, dir, h);
	ev = fmax(rs_rel((double)v, value[0]), fmax(rs_rel((double)t, value[1]), rs_rel((double)t2, value[2])));
	eg = rs_relv(g, grad, 2 * nv);
	eh = rs_relv(h, hv, 2 * nv);
	sprintf(msg, "against Mathematica, %s: E, dE[dir], d2E[dir, dir] %.2g; gradient %.2g; Hv %.2g (status %d %d %d)", tag, ev, eg, eh, (int)s1, (int)s2, (int)s3);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_OK && ev <= tol && eg <= tol && eh <= tol, msg);
	qaws_repulsive_curves_destroy(c);
}

static void test_fixtures(void)
{
	qaws_surface* patch = rs_patch();
	rs_fixture(patch, "a closed hexagon on a bicubic patch", g117_h1_uv, g117_h1_dir, g117_h1_idx, 6, 6, g117_h1_value, g117_h1_grad, g117_h1_hv);
	rs_fixture(patch, "an open polyline and a segment", g117_h2_uv, g117_h2_dir, g117_h2_idx, 7, 5, g117_h2_value, g117_h2_grad, g117_h2_hv);
	qaws_surface_destroy(patch);
}

/* loops on the torus: count loops of n vertices, loop k around v = v0 +
   k / count (the short way), wavy in u, its u shifted across the seam */
static void rs_loops(qaws_surface const* torus, unsigned int count, unsigned int n, qaws_scalar* uv, unsigned int* e)
{
	qaws_range ru = qaws_surface_get_u_range(torus), rv = qaws_surface_get_v_range(torus);
	double U = (double)ru.max_value - ru.min_value, V = (double)rv.max_value - rv.min_value;
	unsigned int k, i;
	for (k = 0; k < count; k++)
		for (i = 0; i < n; i++)
		{
			double t = (double)i / n;
			uv[2 * (k * n + i)] = (qaws_scalar)((double)ru.min_value + U * (0.93 + 0.11 * k / count + 0.05 * sin(2 * RS_PI * t)));   /* (around u = 1: the seam) */
			uv[2 * (k * n + i) + 1] = (qaws_scalar)((double)rv.min_value + V * t);
			e[2 * (k * n + i)] = k * n + i; e[2 * (k * n + i) + 1] = k * n + (i + 1) % n;
		}
}

static void test_seam(void)
{
	enum { N = 24 };
	qaws_curve* profile = NULL;
	qaws_surface* torus = rs_torus(&profile);
	qaws_range ru = qaws_surface_get_u_range(torus);
	qaws_scalar uv[2 * N], uv2[2 * N], g[2 * N], gp[2 * N], E0 = 0, E1 = 0;
	unsigned int e[2 * N], i;
	qaws_repulsive_curves* c;
	qaws_repulsive_surface_desc sd;
	double eh = RS_F32 ? 1e-3 : 1e-6, err = 0, m = 0;
	char msg[300];
	rs_loops(torus, 2, N / 2, uv, e);
	c = rs_curves(N, e, N);
	memset(&sd, 0, sizeof(sd));
	sd.surface = torus; sd.periodic_u = 1; sd.periodic_v = 1;
	for (i = 0; i < 2 * N; i++) { uv2[i] = uv[i]; g[i] = 0; }
	for (i = 0; i < N; i += 3) uv2[2 * i] = (qaws_scalar)((double)uv2[2 * i] - ((double)ru.max_value - ru.min_value));   /* (a period back) */
	qaws_repulsive_curves_on_surface_gradient(c, &sd, uv, 1, g, &E0);
	qaws_repulsive_curves_on_surface_eval(c, &sd, uv2, NULL, &E1, NULL, NULL);
	/* central differences of the gradient, each coordinate */
	for (i = 0; i < 2 * N; i++)
	{
		qaws_scalar ep = 0, em = 0;
		double hstep = RS_F32 ? 1e-3 : 1e-6;
		memcpy(gp, uv, sizeof(uv));
		gp[i] = (qaws_scalar)((double)uv[i] + hstep);
		qaws_repulsive_curves_on_surface_eval(c, &sd, gp, NULL, &ep, NULL, NULL);
		gp[i] = (qaws_scalar)((double)uv[i] - hstep);
		qaws_repulsive_curves_on_surface_eval(c, &sd, gp, NULL, &em, NULL, NULL);
		err = fmax(err, fabs(((double)ep - em) / (2 * hstep) - (double)g[i]));
		m = fmax(m, fabs((double)g[i]));
	}
	sprintf(msg, "torus, two loops across the u seam: E shifted by a period %.2g off; the uv gradient against central differences %.2g", rs_rel((double)E1, (double)E0), err / m);
	TEST_ASSERT(rs_rel((double)E1, (double)E0) <= (RS_F32 ? 1e-5 : 1e-12) && err / m <= eh * 10, msg);
	qaws_repulsive_curves_destroy(c);
	qaws_surface_destroy(torus);
	qaws_curve_destroy(profile);
}

/* the largest distance from the torus (axis z, radii 1 and 0.4) */
static double rs_off_torus(qaws_repulsive_surface_desc const* sd, qaws_scalar const* uv, unsigned int n)
{
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	double m = 0;
	unsigned int i;
	qaws_repulsive_surface_positions(sd, uv, n, x);
	for (i = 0; i < n; i++)
	{
		double rho = sqrt((double)x[3 * i] * x[3 * i] + (double)x[3 * i + 1] * x[3 * i + 1]) - 1;
		m = fmax(m, fabs(sqrt(rho * rho + (double)x[3 * i + 2] * x[3 * i + 2]) - 0.4));
	}
	free(x);
	return m;
}

static void test_solves(void)
{
	enum { N = 36 };
	qaws_curve* profile = NULL;
	qaws_surface* torus = rs_torus(&profile);
	qaws_scalar uv[2 * N], uv0[2 * N];
	unsigned int e[2 * N], pins[2] = { 0, 12 };
	qaws_repulsive_curves* c;
	qaws_repulsive_surface_desc sd;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_report r;
	qaws_status st;
	char msg[300];
	rs_loops(torus, 3, N / 3, uv, e);
	memcpy(uv0, uv, sizeof(uv));
	c = rs_curves(N, e, N);
	memset(&sd, 0, sizeof(sd));
	sd.surface = torus; sd.periodic_u = 1; sd.periodic_v = 1;
	/* energy-based, total length held */
	memset(&d, 0, sizeof(d));
	d.hold_total_length = 1; d.max_iterations = 30;
	memset(&r, 0, sizeof(r));
	st = qaws_repulsive_surface_curves_solve(c, &sd, &d, uv, &r);
	sprintf(msg, "torus, three loops, length held: status %d, reason %d after %u steps, E %.4g -> %.4g, length residual %.2g, off the torus %.2g",
		(int)st, (int)r.reason, r.iterations, (double)r.initial_energy, (double)r.energy, (double)r.constraint_residual, rs_off_torus(&sd, uv, N));
	TEST_ASSERT(st == QAWS_STATUS_OK && r.iterations >= 5 && r.energy < 0.9 * r.initial_energy && (double)r.constraint_residual <= (RS_F32 ? 1e-4 : 1e-7) &&
		rs_off_torus(&sd, uv, N) <= (RS_F32 ? 1e-5 : 1e-12), msg);
	/* collision-checked, two pins */
	{
		qaws_scalar u0[2] = { uv0[0], uv0[1] }, u1[2] = { uv0[24], uv0[25] };
		memcpy(uv, uv0, sizeof(uv));
		memset(&d, 0, sizeof(d));
		d.hold_total_length = 1; d.max_iterations = 20; d.collision_checked = 1; d.pinned = pins; d.pinned_count = 2;
		memset(&r, 0, sizeof(r));
		st = qaws_repulsive_surface_curves_solve(c, &sd, &d, uv, &r);
		sprintf(msg, "torus, collision-checked, two pins: status %d, reason %d after %u steps, E %.4g -> %.4g; pins moved %.2g",
			(int)st, (int)r.reason, r.iterations, (double)r.initial_energy, (double)r.energy,
			fmax(fabs((double)uv[0] - u0[0]) + fabs((double)uv[1] - u0[1]), fabs((double)uv[24] - u1[0]) + fabs((double)uv[25] - u1[1])));
		TEST_ASSERT(st == QAWS_STATUS_OK && r.iterations >= 3 && r.energy < r.initial_energy &&
			fmax(fabs((double)uv[0] - u0[0]) + fabs((double)uv[1] - u0[1]), fabs((double)uv[24] - u1[0]) + fabs((double)uv[25] - u1[1])) <= (RS_F32 ? 1e-6 : 1e-12), msg);   /* (held by the projection: to rounding) */
	}
	qaws_repulsive_curves_destroy(c);
	qaws_surface_destroy(torus);
	qaws_curve_destroy(profile);
}

static void test_domain(void)
{
	enum { N = 10 };
	qaws_surface* patch = rs_patch();
	qaws_scalar uv[2 * N];
	unsigned int e[2 * N], i;
	qaws_repulsive_curves* c;
	qaws_repulsive_surface_desc sd;
	qaws_repulsive_curves_solve_desc d;
	qaws_repulsive_report r;
	qaws_status st;
	double lo = 1, hi = 0;
	char msg[300];
	/* an open polyline near the patch's corner, its length free: the energy pushes it out */
	for (i = 0; i < N; i++)
	{
		uv[2 * i] = (qaws_scalar)(0.08 + 0.05 * i); uv[2 * i + 1] = (qaws_scalar)(0.06 + 0.02 * sin((double)i));
		e[2 * i] = i; e[2 * i + 1] = i + 1;
	}
	c = rs_curves(N, e, N - 1);
	memset(&sd, 0, sizeof(sd));
	sd.surface = patch;
	memset(&d, 0, sizeof(d));
	d.max_iterations = 40;
	memset(&r, 0, sizeof(r));
	st = qaws_repulsive_surface_curves_solve(c, &sd, &d, uv, &r);
	for (i = 0; i < 2 * N; i++) { lo = fmin(lo, (double)uv[i]); hi = fmax(hi, (double)uv[i]); }
	sprintf(msg, "patch, an open polyline growing: status %d, reason %d after %u steps (rejected %u), E %.4g -> %.4g, parameters in [%.4g, %.4g]",
		(int)st, (int)r.reason, r.iterations, r.rejected_energy, (double)r.initial_energy, (double)r.energy, lo, hi);
	TEST_ASSERT(st == QAWS_STATUS_OK && r.iterations >= 5 && r.energy < r.initial_energy && lo >= 0 && hi <= 1, msg);
	/* refusals: outside the range; a sphere's pole */
	{
		qaws_scalar bad[2 * N];
		qaws_repulsive_surface_curves_solver* s = NULL;
		qaws_status s1, s2, s3, s1b;
		qaws_curve* prof = NULL;
		qaws_surface* sphere = NULL;
		qaws_arc_segment a;
		qaws_arc_desc ad;
		qaws_surface_revolution_desc rd;
		memcpy(bad, uv, sizeof(bad));
		bad[4] = (qaws_scalar)1.2;
		s1 = qaws_repulsive_surface_curves_solver_create(c, &sd, NULL, bad, &s);
		qaws_repulsive_surface_curves_solver_destroy(s);
		s = NULL;
		/* a sphere: the half circle from the south pole to the north, revolved */
		memset(&a, 0, sizeof(a));
		a.radius = 1; a.angle_start = (qaws_scalar)(-RS_PI / 2); a.angle_end = (qaws_scalar)(RS_PI / 2);
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = &a; ad.segment_count = 1;
		qaws_curve_create_arc(&ad, &prof);
		memset(&rd, 0, sizeof(rd));
		rd.profile = prof; rd.axis_direction.z = 1;
		qaws_surface_create_revolution(&rd, &sphere);
		{
			qaws_repulsive_surface_desc ss;
			qaws_range rv = qaws_surface_get_v_range(sphere);
			memset(&ss, 0, sizeof(ss));
			ss.surface = sphere; ss.periodic_u = 1;
			for (i = 0; i < N; i++) { bad[2 * i] = (qaws_scalar)(0.1 * i); bad[2 * i + 1] = (qaws_scalar)(0.5 * ((double)rv.min_value + rv.max_value) + 0.01 * i); }
			bad[1] = rv.max_value;   /* (the north pole) */
			s2 = qaws_repulsive_surface_curves_solver_create(c, &ss, NULL, bad, &s);
			qaws_repulsive_surface_curves_solver_destroy(s);
		}
		qaws_surface_destroy(sphere);
		qaws_curve_destroy(prof);
		{
			/* a hole: the circle of radius 0.15 about (0.5, 0.5) in uv */
			qaws_arc_segment h;
			qaws_curve* hc = NULL;
			qaws_curve const* hl[1];
			qaws_trim_loop loop;
			qaws_surface_trim_desc td;
			qaws_surface* trimmed = NULL;
			qaws_repulsive_surface_desc ts;
			memset(&h, 0, sizeof(h));
			h.center[0] = h.center[1] = (qaws_scalar)0.5; h.radius = (qaws_scalar)0.15; h.angle_end = (qaws_scalar)(2 * RS_PI);
			ad.segments = &h;
			qaws_curve_create_arc(&ad, &hc);
			hl[0] = hc;
			loop.curves = hl; loop.curve_count = 1; loop.is_outer = 0;
			memset(&td, 0, sizeof(td));
			td.base = patch; td.loops = &loop; td.loop_count = 1;
			qaws_surface_create_trimmed(&td, &trimmed);
			memset(&ts, 0, sizeof(ts));
			ts.surface = trimmed;
			memcpy(bad, uv, sizeof(bad));
			s1b = qaws_repulsive_surface_curves_solver_create(c, &ts, NULL, bad, &s);
			qaws_repulsive_surface_curves_solver_destroy(s);
			s = NULL;
			bad[6] = bad[7] = (qaws_scalar)0.52;   /* (in the hole) */
			s3 = qaws_repulsive_surface_curves_solver_create(c, &ts, NULL, bad, &s);
			qaws_repulsive_surface_curves_solver_destroy(s);
			qaws_surface_destroy(trimmed);
			qaws_curve_destroy(hc);
		}
		sprintf(msg, "refusals: a vertex outside the range %d, in a trimmed hole %d (outside it %d), at a sphere's pole %d", (int)s1, (int)s3, (int)s1b, (int)s2);
		TEST_ASSERT(s1 == QAWS_STATUS_OUT_OF_RANGE && s3 == QAWS_STATUS_OUT_OF_RANGE && s1b == QAWS_STATUS_OK && s2 == QAWS_STATUS_OUT_OF_RANGE, msg);
	}
	qaws_repulsive_curves_destroy(c);
	qaws_surface_destroy(patch);
}

int test_117_repulsive_surface_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 117: curves on a parametric surface\n");
	test_fixtures();
	test_seam();
	test_solves();
	test_domain();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
