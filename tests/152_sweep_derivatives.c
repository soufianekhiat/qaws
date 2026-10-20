/*
 * Test 152: analytic derivatives of swept and pipe surfaces
 *
 *   - a profile swept along a helix-like cubic path by its Frenet frame: du,
 *     the mixed derivative and duu against Richardson-extrapolated
 *     differences of the positions (and of dv);
 *   - a pipe around a twisted path: the same, the cross-section at its
 *     radius exactly (an orthonormal frame perpendicular to the tangent),
 *     du continuous across the frame table's samples;
 *   - the normal from the analytic du: unit, perpendicular to du and dv.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

static qaws_curve* sd_path(int twisted)
{
	qaws_vec3 p[9];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 9; i++)
	{
		double a = 0.8 * i;
		p[i].x = (qaws_scalar)(2 * cos(a) + (twisted ? 0.3 * sin(2.1 * a) : 0));
		p[i].y = (qaws_scalar)(2 * sin(a));
		p[i].z = (qaws_scalar)(0.6 * i + (twisted ? 0.4 * cos(1.7 * a) : 0));
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = p; d.control_point_count = 9;
	qaws_curve_create_bspline(&c, &d);
	return c;
}

#define SD_F32 (sizeof(qaws_scalar) != sizeof(double))
/* float bounds, relative to 1 + |value|, at the step below: the measured errors (du 1e-5, duv 6e-5, duu 2e-3) with a
   margin of five to ten for other grids of (u, v) */
#define SD_DU 1e-4
#define SD_DUV 5e-4
#define SD_DUU 1e-2

static qaws_vec3 sd_pos(qaws_surface const* s, double u, double v, unsigned int flag)
{
	qaws_surface_eval_result r;
	qaws_surface_evaluate(&r, s, (qaws_scalar)u, (qaws_scalar)v, flag | QAWS_SURFACE_EVAL_POSITION);
	return flag == QAWS_SURFACE_EVAL_DV ? r.dv : r.position;
}

static double sd_diff(qaws_vec3 a, double x, double y, double z)
{
	return fabs(a.x - x) + fabs(a.y - y) + fabs(a.z - z);
}

/* the largest differences of du, duv, duu from differences of positions (dv for duv) over a grid of (u, v), relative to
   1 + |value| */
static void sd_check(qaws_surface const* s, char const* name, double* e_du, double* e_duv, double* e_duu)
{
	unsigned int i, j;
	*e_du = *e_duv = *e_duu = 0;
	for (i = 0; i < 12; i++)
		for (j = 0; j < 7; j++)
		{
			/* mid-span in the path (six spans): its third derivative jumps at the knots, and so does the Frenet frame's rate */
			/* the step: 1e-3 in double. float: Richardson's truncation (h^4) against the rounding of the positions (eps |P| / h
			   for du, eps |P| / h^2 for duu); measured over h = 0.04 .. 0.005, du and duv are smallest at h = 0.01 (1e-5, 6e-5),
			   duu at 2e-3 there (rounding, rising as 1 / h^2). The step stays inside the span (mid-span is 1 / 12 from the knots) */
			double u = (i + 0.5) / 12, v = 0.141 * j + 0.03, h = SD_F32 ? 0.01 : 1e-3, q[3], m[3];
			qaws_surface_eval_result r;
			qaws_vec3 a, b, c, dd, a2, b2, va, vb, va2, vb2;
			unsigned int k;
			qaws_surface_evaluate(&r, s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DUU
				| QAWS_SURFACE_EVAL_DUV | QAWS_SURFACE_EVAL_DV);
			a = sd_pos(s, u + h, v, 0); b = sd_pos(s, u - h, v, 0);
			a2 = sd_pos(s, u + h / 2, v, 0); b2 = sd_pos(s, u - h / 2, v, 0);
			c = r.position;
			/* Richardson on central differences: (4 D(h/2) - D(h)) / 3 */
			q[0] = (4 * (a2.x - b2.x) / h - (a.x - b.x) / (2 * h)) / 3;
			q[1] = (4 * (a2.y - b2.y) / h - (a.y - b.y) / (2 * h)) / 3;
			q[2] = (4 * (a2.z - b2.z) / h - (a.z - b.z) / (2 * h)) / 3;
			*e_du = fmax(*e_du, sd_diff(r.du, q[0], q[1], q[2]) / (1 + fabs(q[0]) + fabs(q[1]) + fabs(q[2])));
			/* Richardson on second differences */
			for (k = 0; k < 3; k++)
			{
				double pa = k == 0 ? a.x : k == 1 ? a.y : a.z, pb = k == 0 ? b.x : k == 1 ? b.y : b.z, pc = k == 0 ? c.x : k == 1 ? c.y : c.z;
				double pa2 = k == 0 ? a2.x : k == 1 ? a2.y : a2.z, pb2 = k == 0 ? b2.x : k == 1 ? b2.y : b2.z;
				m[k] = (4 * (pa2 - 2 * pc + pb2) / (h * h / 4) - (pa - 2 * pc + pb) / (h * h)) / 3;
			}
			*e_duu = fmax(*e_duu, sd_diff(r.duu, m[0], m[1], m[2]) / (1 + fabs(m[0]) + fabs(m[1]) + fabs(m[2])));
			va = sd_pos(s, u + h, v, QAWS_SURFACE_EVAL_DV); vb = sd_pos(s, u - h, v, QAWS_SURFACE_EVAL_DV);
			va2 = sd_pos(s, u + h / 2, v, QAWS_SURFACE_EVAL_DV); vb2 = sd_pos(s, u - h / 2, v, QAWS_SURFACE_EVAL_DV);
			dd.x = (qaws_scalar)((4 * (va2.x - vb2.x) / h - (va.x - vb.x) / (2 * h)) / 3);
			dd.y = (qaws_scalar)((4 * (va2.y - vb2.y) / h - (va.y - vb.y) / (2 * h)) / 3);
			dd.z = (qaws_scalar)((4 * (va2.z - vb2.z) / h - (va.z - vb.z) / (2 * h)) / 3);
			*e_duv = fmax(*e_duv, sd_diff(r.duv, dd.x, dd.y, dd.z) / (1 + fabs(dd.x) + fabs(dd.y) + fabs(dd.z)));
		}
	printf("  %s: du %.1e, duv %.1e, duu %.1e from differences of the positions\n", name, *e_du, *e_duv, *e_duu);
}

static void test_swept(void)
{
	qaws_vec2 pp[5] = { { 0.4f, 0 }, { 0.4f, 0.4f }, { -0.4f, 0.4f }, { -0.4f, -0.4f }, { 0.3f, -0.3f } };
	qaws_bspline_desc pd;
	qaws_surface_swept_desc d;
	qaws_curve *path = sd_path(0), *prof = NULL;
	qaws_surface* s = NULL;
	double e_du, e_duv, e_duu;
	memset(&pd, 0, sizeof pd);
	pd.dimension = QAWS_DIMENSION_2D; pd.degree = 3; pd.control_points = pp; pd.control_point_count = 5;
	qaws_curve_create_bspline(&prof, &pd);
	memset(&d, 0, sizeof d);
	d.path = path; d.profile = prof;
	TEST_ASSERT(qaws_surface_create_swept(&s, &d) == QAWS_STATUS_OK, "swept surface created");
	sd_check(s, "swept (Frenet)", &e_du, &e_duv, &e_duu);
	TEST_ASSERT(e_du <= (SD_F32 ? SD_DU : 1e-8), "swept: du analytic");
	TEST_ASSERT(e_duv <= (SD_F32 ? SD_DUV : 1e-7), "swept: mixed derivative analytic");
	TEST_ASSERT(e_duu <= (SD_F32 ? SD_DUU : 1e-5), "swept: duu from the analytic du");
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
	qaws_curve_destroy(path);
}

static void test_pipe(void)
{
	qaws_curve* path = sd_path(1);
	qaws_surface_pipe_desc d;
	qaws_surface* s = NULL;
	double e_du, e_duv, e_duu, radius = 0, jump = 0;
	unsigned int i;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.path = path; d.radius_x = (qaws_scalar)0.35;
	TEST_ASSERT(qaws_surface_create_pipe(&s, &d) == QAWS_STATUS_OK, "pipe created");
	sd_check(s, "pipe (rotation-minimizing)", &e_du, &e_duv, &e_duu);
	TEST_ASSERT(e_du <= (SD_F32 ? SD_DU : 1e-7), "pipe: du analytic");
	TEST_ASSERT(e_duv <= (SD_F32 ? SD_DUV : 1e-6), "pipe: mixed derivative analytic");
	TEST_ASSERT(e_duu <= (SD_F32 ? SD_DUU : 1e-5), "pipe: duu from the analytic du");
	for (i = 0; i < 400; i++)
	{
		double u = (i + 0.5) / 400, v = 0.37 * i - floor(0.37 * i);
		qaws_eval_result_3d c;
		qaws_vec3 p = sd_pos(s, u, v, 0);
		qaws_range r = qaws_curve_get_parameter_range(path);
		qaws_curve_evaluate_3d(&c, path, r.min_value + (qaws_scalar)u * (r.max_value - r.min_value), QAWS_EVAL_FLAG_POSITION);
		radius = fmax(radius, fabs(sqrt(((double)p.x - c.position.x) * (p.x - c.position.x) + ((double)p.y - c.position.y) * (p.y - c.position.y)
			+ ((double)p.z - c.position.z) * (p.z - c.position.z)) - 0.35));
	}
	snprintf(msg, sizeof msg, "pipe: the section at its radius (%.1e)", radius);
	/* (pure rounding: the frame is orthonormal; a few eps of the path's size, 4) */
	TEST_ASSERT(radius <= fmax(1e-12, 8 * TEST_EPS * 4), msg);
	/* across the samples of the frame table */
	for (i = 1; i < 256; i += 17)
	{
		qaws_surface_eval_result a, b;
		double u = i / 256.0;
		qaws_surface_evaluate(&a, s, (qaws_scalar)(u - 1e-12), (qaws_scalar)0.3, QAWS_SURFACE_EVAL_DU);
		qaws_surface_evaluate(&b, s, (qaws_scalar)(u + 1e-12), (qaws_scalar)0.3, QAWS_SURFACE_EVAL_DU);
		jump = fmax(jump, sd_diff(a.du, b.du.x, b.du.y, b.du.z));
	}
	snprintf(msg, sizeof msg, "pipe: du continuous at the frame samples (jump %.1e)", jump);
	printf("  %s\n", msg);
	TEST_ASSERT(jump <= 1e-8, msg);
	/* the normal */
	{
		qaws_surface_eval_result r;
		double tu, tv;
		qaws_surface_evaluate(&r, s, (qaws_scalar)0.41, (qaws_scalar)0.77, QAWS_SURFACE_EVAL_NORMAL | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV);
		/* rounding: a few eps, relative to |du| and |dv| for the products with them */
		tu = fmax(1e-12, 8 * TEST_EPS * sqrt((double)r.du.x * r.du.x + (double)r.du.y * r.du.y + (double)r.du.z * r.du.z));
		tv = fmax(1e-12, 8 * TEST_EPS * sqrt((double)r.dv.x * r.dv.x + (double)r.dv.y * r.dv.y + (double)r.dv.z * r.dv.z));
		TEST_ASSERT(fabs(sqrt((double)r.normal.x * r.normal.x + (double)r.normal.y * r.normal.y + (double)r.normal.z * r.normal.z) - 1) <= fmax(1e-12, 8 * TEST_EPS)
			&& fabs((double)r.normal.x * r.du.x + (double)r.normal.y * r.du.y + (double)r.normal.z * r.du.z) <= tu
			&& fabs((double)r.normal.x * r.dv.x + (double)r.normal.y * r.dv.y + (double)r.normal.z * r.dv.z) <= tv, "pipe: unit normal across du and dv");
	}
	qaws_surface_destroy(s);
	qaws_curve_destroy(path);
}

int test_152_sweep_derivatives_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 152: analytic derivatives of swept and pipe surfaces\n");
	test_swept();
	test_pipe();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
