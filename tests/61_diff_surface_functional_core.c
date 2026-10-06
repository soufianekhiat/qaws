/*
 * Test 61: Backend-neutral surface functional kernels
 *
 * core/qaws_bspline_surface_functional_core.h compiled on the C backend,
 * driven per quadrature node (value, tangent, second tangent) and per
 * control point (gather of the gradient over the nodes it supports), must
 * reproduce qaws_surface_functional_eval / _gradient for the area, thin
 * plate and Willmore functionals of a nonuniform B-spline surface whose
 * cells do not align with its knots.
 */

#include "test_diff.h"
#include "qaws_diff_functionals.h"
#include "core/qaws_bspline_surface_functional_core.h"

#define SF_NU 5
#define SF_NV 6
#define SF_DEG 3
#define SF_CELLS 3

static qaws_scalar const g_sf_uk[SF_NU + SF_DEG + 1] = { 0, 0, 0, 0, 0.4f, 1, 1, 1, 1 };
static qaws_scalar const g_sf_vk[SF_NV + SF_DEG + 1] = { 0, 0, 0, 0, 0.3f, 0.65f, 1, 1, 1, 1 };

typedef struct sf_window
{
	qaws_scalar cp[QAWS_CORE_SURFACE_WINDOW], dot[QAWS_CORE_SURFACE_WINDOW];
	qaws_scalar uk[QAWS_CORE_MAX_POINTS * 2], vk[QAWS_CORE_MAX_POINTS * 2];
	int ufirst, vfirst;
} sf_window;

static void sf_load(qaws_scalar const* cps, qaws_scalar const* dir, qaws_scalar u, qaws_scalar v, sf_window* w)
{
	qaws_scalar uk[QAWS_CORE_MAX_POINTS * 2], vk[QAWS_CORE_MAX_POINTS * 2];
	int us, vs, i, j, c;
	memset(w, 0, sizeof(*w));
	memset(uk, 0, sizeof(uk));
	memset(vk, 0, sizeof(vk));
	memcpy(uk, g_sf_uk, sizeof(g_sf_uk));
	memcpy(vk, g_sf_vk, sizeof(g_sf_vk));
	us = qaws_find_span(uk, SF_DEG, SF_NU, u);
	vs = qaws_find_span(vk, SF_DEG, SF_NV, v);
	for (i = 0; i < 2 * (SF_DEG + 1); i++)
	{
		w->uk[i] = g_sf_uk[us - SF_DEG + i];
		w->vk[i] = g_sf_vk[vs - SF_DEG + i];
	}
	w->ufirst = us - SF_DEG;
	w->vfirst = vs - SF_DEG;
	for (i = 0; i <= SF_DEG; i++)
		for (j = 0; j <= SF_DEG; j++)
			for (c = 0; c < 3; c++)
			{
				int g = ((w->ufirst + i) * SF_NV + w->vfirst + j) * 3 + c, l = (i * (SF_DEG + 1) + j) * 3 + c;
				w->cp[l] = cps[g];
				w->dot[l] = dir ? dir[g] : 0;
			}
}

static void test_functional(int id, qaws_surface_functional functional, char const* name)
{
	qaws_scalar cps[SF_NU * SF_NV * 3], dir[SF_NU * SF_NV * 3], grad[SF_NU * SF_NV * 3], value = 0, t1 = 0, t2 = 0, gv = 0;
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	qaws_field_view fd, fg;
	qaws_diff_views vd, vg;
	double v = 0, d1 = 0, d2 = 0, tol = QAWS_SCALAR_IS_FLOAT ? 2e-3 : 1e-10, gn = 0, kn = 0;
	int i, j, ci, cj, q, ok_g = 1;
	char msg[200];
	diff_seed(6161u);
	for (i = 0; i < SF_NU; i++)
		for (j = 0; j < SF_NV; j++)
		{
			qaws_scalar* p = &cps[(i * SF_NV + j) * 3];
			p[0] = (qaws_scalar)(0.8 * i);
			p[1] = (qaws_scalar)(0.7 * j);
			p[2] = (qaws_scalar)(0.35 * sin(1.1 * i + 0.8 * j) + 0.1 * i * j / 5.0);
		}
	diff_rand_fill(dir, SF_NU * SF_NV * 3);
	memset(&d, 0, sizeof(d));
	d.u_degree = SF_DEG;
	d.v_degree = SF_DEG;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = SF_NU;
	d.v_point_count = SF_NV;
	d.u_knots = g_sf_uk;
	d.u_knot_count = SF_NU + SF_DEG + 1;
	d.v_knots = g_sf_vk;
	d.v_knot_count = SF_NV + SF_DEG + 1;
	TEST_ASSERT_STATUS(qaws_surface_create_bspline(&d, &s));
	fd = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, dir, SF_NU * SF_NV, 3);
	fg = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, grad, SF_NU * SF_NV, 3);
	vd.fields = &fd;
	vg.fields = &fg;
	vd.field_count = vg.field_count = 1;
	vd.children = vg.children = NULL;
	vd.child_count = vg.child_count = 0;
	memset(grad, 0, sizeof(grad));
	TEST_ASSERT_STATUS(qaws_surface_functional_eval(NULL, s, functional, SF_CELLS, &vd, &value, &t1, &t2));
	TEST_ASSERT_STATUS(qaws_surface_functional_gradient(NULL, s, functional, SF_CELLS, &vg, &gv));

	/* one thread per node, then a sum */
	for (ci = 0; ci < SF_CELLS; ci++)
		for (cj = 0; cj < SF_CELLS; cj++)
			for (q = 0; q < 16; q++)
			{
				qaws_scalar node[3];
				sf_window w;
				qaws_dual1 f;
				qaws_surface_cell_node(0, 1, 0, 1, (qaws_scalar)ci, (qaws_scalar)cj, SF_CELLS, q, node);
				sf_load(cps, dir, node[0], node[1], &w);
				f = qaws_bspline_surface_node(id, w.cp, w.dot, w.uk, w.vk, SF_DEG, SF_DEG, node[0], node[1], node[2]);
				v += f.v;
				d1 += f.t;
				d2 += f.tt;
			}
	sprintf(msg, "%s: node kernels give the runtime value and tangents (%.12g vs %.12g)", name, v, (double)value);
	TEST_ASSERT(diff_close(v, value, tol) && diff_close(d1, t1, tol) && diff_close(d2, t2, tol), msg);

	/* one thread per control point: gather over the nodes whose window holds it */
	for (i = 0; i < SF_NU; i++)
		for (j = 0; j < SF_NV; j++)
		{
			qaws_vec3 acc = qaws_v3_zero();
			int k = (i * SF_NV + j) * 3;
			for (ci = 0; ci < SF_CELLS; ci++)
				for (cj = 0; cj < SF_CELLS; cj++)
					for (q = 0; q < 16; q++)
					{
						qaws_scalar node[3];
						sf_window w;
						qaws_surface_cell_node(0, 1, 0, 1, (qaws_scalar)ci, (qaws_scalar)cj, SF_CELLS, q, node);
						sf_load(cps, NULL, node[0], node[1], &w);
						if (i < w.ufirst || i > w.ufirst + SF_DEG || j < w.vfirst || j > w.vfirst + SF_DEG)
							continue;
						acc = qaws_v3_add(acc, qaws_bspline_surface_node_adjoint_cp(id, w.cp, w.uk, w.vk, SF_DEG, SF_DEG, node[0], node[1],
							i - w.ufirst, j - w.vfirst, node[2]));
					}
			ok_g &= diff_close(acc.x, grad[k], tol) && diff_close(acc.y, grad[k + 1], tol) && diff_close(acc.z, grad[k + 2], tol);
			gn += (double)grad[k] * grad[k] + (double)grad[k + 1] * grad[k + 1] + (double)grad[k + 2] * grad[k + 2];
			kn += (double)acc.x * acc.x + (double)acc.y * acc.y + (double)acc.z * acc.z;
		}
	printf("    %s: value %.15g (runtime %.15g), |g| %.15g (runtime %.15g)\n", name, v, (double)value, sqrt(kn), sqrt(gn));
	sprintf(msg, "%s: gather over control points gives the runtime gradient", name);
	TEST_ASSERT(ok_g, msg);
	qaws_surface_destroy(s);
}

int test_61_diff_surface_functional_core_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 61: Backend-neutral surface functional kernels\n");
	test_functional(QAWS_CORE_SURFACE_AREA, QAWS_FUNCTIONAL_AREA, "area");
	test_functional(QAWS_CORE_SURFACE_THIN_PLATE, QAWS_FUNCTIONAL_THIN_PLATE, "thin plate");
	test_functional(QAWS_CORE_SURFACE_WILLMORE, QAWS_FUNCTIONAL_WILLMORE, "willmore");
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
