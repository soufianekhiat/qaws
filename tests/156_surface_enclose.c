/*
 * Test 156: guaranteed ranges of surfaces over parameter boxes (qaws_surface_enclose)
 *
 *   - a NURBS torus (rational), a cubic B-spline surface (several spans), a bicubic Bezier patch: over boxes of a half, a
 *     tenth and a hundredth of the domain, every sampled S, S_u, S_v, S_uu, S_uv, S_vv inside its range;
 *   - the ranges shrink with the box (a point box: the position to rounding);
 *   - the refusals: a box past the domain, an empty one, an order past 2, an unsupported surface.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

static qaws_surface* se_torus(void)
{
	static double const cx[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, cy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	static qaws_scalar const kn[12] = { 0, 0, 0, 0.25f, 0.25f, 0.5f, 0.5f, 0.75f, 0.75f, 1, 1, 1 };
	qaws_vec3 cp[81];
	qaws_scalar w[81];
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 9; i++)
		for (j = 0; j < 9; j++)
		{
			double const rad = 2 + 0.7 * cx[j];
			cp[i * 9 + j].x = (qaws_scalar)(rad * cx[i]);
			cp[i * 9 + j].y = (qaws_scalar)(rad * cy[i]);
			cp[i * 9 + j].z = (qaws_scalar)(0.7 * cy[j]);
			w[i * 9 + j] = (qaws_scalar)((i & 1 ? sqrt(0.5) : 1) * (j & 1 ? sqrt(0.5) : 1));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 2, d.v_degree = 2;
	d.control_points = cp, d.u_point_count = 9, d.v_point_count = 9;
	d.weights = w;
	d.u_knots = kn, d.u_knot_count = 12, d.v_knots = kn, d.v_knot_count = 12;
	qaws_surface_create_nurbs(&s, &d);
	return s;
}

static qaws_surface* se_bspline(void)
{
	static qaws_scalar const kn[10] = { 0, 0, 0, 0, 0.3f, 0.6f, 1, 1, 1, 1 };
	qaws_vec3 cp[36];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
		{
			cp[i * 6 + j].x = (qaws_scalar)i;
			cp[i * 6 + j].y = (qaws_scalar)j;
			cp[i * 6 + j].z = (qaws_scalar)sin(1.3 * i + 0.7 * j);
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3, d.v_degree = 3;
	d.control_points = cp, d.u_point_count = 6, d.v_point_count = 6;
	d.u_knots = kn, d.u_knot_count = 10, d.v_knots = kn, d.v_knot_count = 10;
	qaws_surface_create_bspline(&s, &d);
	return s;
}

static qaws_surface* se_bezier(void)
{
	qaws_vec3 cp[16];
	qaws_surface_bezier_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			cp[i * 4 + j].x = (qaws_scalar)i;
			cp[i * 4 + j].y = (qaws_scalar)(j + 0.3 * i);
			cp[i * 4 + j].z = (qaws_scalar)(((i * 3 + j * 5) % 7) * 0.25);
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3, d.v_degree = 3;
	d.control_points = cp, d.u_point_count = 4, d.v_point_count = 4;
	qaws_surface_create_bezier(&s, &d);
	return s;
}

/* every sample of the box inside the ranges (a slack for the evaluation's own rounding); the ranges' widest width */
static int se_check(qaws_surface const* s, double u0, double u1, double v0, double v1, double* widest)
{
	double lo[18], hi[18];
	unsigned int i, j, q, k;
	int ok = qaws_surface_enclose(lo, hi, s, u0, u1, v0, v1, 2) == QAWS_STATUS_OK;
	double const slack = sizeof(qaws_scalar) == sizeof(double) ? 1e-9 : 2e-3;
	*widest = 0;
	for (q = 0; ok && q < 18; q++)
		*widest = fmax(*widest, hi[q] - lo[q]);
	for (i = 0; ok && i <= 8; i++)
		for (j = 0; ok && j <= 8; j++)
		{
			qaws_surface_eval_result r;
			qaws_vec3 const* v[6];
			double const u = u0 + (u1 - u0) * i / 8, w = v0 + (v1 - v0) * j / 8;
			qaws_surface_evaluate(&r, s, (qaws_scalar)u, (qaws_scalar)w,
				QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_DUU | QAWS_SURFACE_EVAL_DUV
					| QAWS_SURFACE_EVAL_DVV);
			v[0] = &r.position, v[1] = &r.du, v[2] = &r.dv, v[3] = &r.duu, v[4] = &r.duv, v[5] = &r.dvv;
			for (q = 0; q < 6; q++)
				for (k = 0; k < 3; k++)
				{
					double const x = k == 0 ? (double)v[q]->x : k == 1 ? (double)v[q]->y : (double)v[q]->z;
					double const scale = 1 + fabs(x);
					if (x < lo[q * 3 + k] - slack * scale || x > hi[q * 3 + k] + slack * scale)
						ok = 0;
				}
		}
	return ok;
}

int test_156_surface_enclose_main(void)
{
	qaws_surface* surfaces[3];
	char const* names[3] = { "a NURBS torus", "a cubic B-spline surface", "a bicubic Bezier patch" };
	char msg[300];
	unsigned int k;
	g_pass = 0;
	g_fail = 0;
	printf("Test 156: guaranteed ranges of surfaces over parameter boxes (qaws_surface_enclose)\n");
	surfaces[0] = se_torus();
	surfaces[1] = se_bspline();
	surfaces[2] = se_bezier();
	TEST_ASSERT(surfaces[0] && surfaces[1] && surfaces[2], "a torus, a B-spline surface, a Bezier patch");
	for (k = 0; k < 3; k++)
	{
		static double const size[3] = { 0.5, 0.1, 0.01 };
		double widths[3] = { 0, 0, 0 }, point_lo[3], point_hi[3];
		unsigned int z;
		int all = 1;
		if (!surfaces[k])
			continue;
		for (z = 0; z < 3; z++)
		{
			double const u0 = 0.37 - size[z] / 2, v0 = 0.61 - size[z] / 2;
			all = all && se_check(surfaces[k], u0, u0 + size[z], v0, v0 + size[z], &widths[z]);
		}
		qaws_surface_enclose(point_lo, point_hi, surfaces[k], 0.37, 0.37, 0.61, 0.61, 0);
		sprintf(msg, "%s: every sample inside; widest range %.3g, %.3g, %.3g over boxes of 0.5, 0.1, 0.01; a point's %.1e", names[k], widths[0],
			widths[1], widths[2], fmax(point_hi[0] - point_lo[0], fmax(point_hi[1] - point_lo[1], point_hi[2] - point_lo[2])));
		printf("    %s\n", msg);
		TEST_ASSERT(all && widths[1] < widths[0] && widths[2] < widths[1] && point_hi[0] - point_lo[0] < 1e-12, msg);
	}
	/* refusals */
	{
		double lo[18], hi[18];
		qaws_surface* ex = NULL;
		qaws_surface_extrusion_desc ed;
		qaws_curve* line = NULL;
		qaws_vec3 pts[2];
		pts[0].x = 0, pts[0].y = 0, pts[0].z = 0, pts[1].x = 1, pts[1].y = 0, pts[1].z = 0;
		TEST_ASSERT(qaws_surface_enclose(lo, hi, surfaces[0], -0.1, 0.2, 0, 1, 0) == QAWS_STATUS_INVALID_ARGUMENT, "a box past the domain refused");
		TEST_ASSERT(qaws_surface_enclose(lo, hi, surfaces[0], 0.3, 0.2, 0, 1, 0) == QAWS_STATUS_INVALID_ARGUMENT, "an empty box refused");
		TEST_ASSERT(qaws_surface_enclose(lo, hi, surfaces[0], 0, 1, 0, 1, 3) == QAWS_STATUS_INVALID_ARGUMENT, "an order past 2 refused");
		{
			qaws_bezier_desc bd;
			memset(&bd, 0, sizeof bd);
			bd.dimension = QAWS_DIMENSION_3D, bd.degree = 1, bd.control_points = pts, bd.control_point_count = 2;
			qaws_curve_create_bezier(&line, &bd);
		}
		memset(&ed, 0, sizeof ed);
		ed.profile = line;
		ed.direction.z = 1;
		ed.length = 1;
		qaws_surface_create_extrusion(&ex, &ed);
		TEST_ASSERT(!ex || qaws_surface_enclose(lo, hi, ex, 0, 0.5, 0, 0.5, 0) == QAWS_STATUS_UNSUPPORTED_OPERATION, "an extrusion unsupported");
		qaws_surface_destroy(ex);
		qaws_curve_destroy(line);
	}
	for (k = 0; k < 3; k++)
		qaws_surface_destroy(surfaces[k]);
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
