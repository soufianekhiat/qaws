/*
 * Test 184: implicit sweeps (qaws_mesh_3d_sweep_implicit), small meshes, coarse grids
 *
 *   - a unit cube translated by (1, 1, 0): closed, its volume within 2% of the Minkowski sum's 3, every vertex of either
 *     sweep within a cell of the other (the explicit one exact here);
 *   - a box turned a quarter about z: the implicit and the explicit sweeps within a cell of each other both ways, their
 *     volumes within 2%;
 *   - a vertex motion (a box's top rising from 1 to 2): closed, volume 2 within 3%;
 *   - refusals: no mesh, a vertex motion of another count.
 */

#include "test_common.h"
#include "qaws_mesh_sweep.h"
#include "qaws_mesh_bvh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static double si_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = m ? qaws_mesh_3d_get_triangle_count(m) : 0, i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		v += ((double)a.x * ((double)b.y * e.z - (double)b.z * e.y) + (double)a.y * ((double)b.z * e.x - (double)b.x * e.z) +
			(double)a.z * ((double)b.x * e.y - (double)b.y * e.x)) / 6;
	}
	return v;
}

/* the farthest vertex of a from b's surface */
static double si_far(qaws_mesh_3d const* a, qaws_mesh_3d const* b)
{
	qaws_mesh_bvh* t = NULL;
	unsigned int n = a && b ? qaws_mesh_3d_get_vertex_count(a) : 0, i;
	double mx = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (3 * (size_t)n + 1));
	qaws_mesh_hit* h = (qaws_mesh_hit*)malloc(sizeof(qaws_mesh_hit) * ((size_t)n + 1));
	if (!n || qaws_mesh_bvh_create_mesh_3d(&t, b, NULL) != QAWS_STATUS_OK)
		mx = 1e30;
	else
	{
		qaws_mesh_3d_get_positions(x, a);
		qaws_mesh_bvh_closest(h, t, x, n, 0, NULL);
		for (i = 0; i < n; i++)
			mx = fmax(mx, (double)h[i].distance);
	}
	qaws_mesh_bvh_destroy(t);
	free(x); free(h);
	return mx;
}

static qaws_mesh_3d* si_box(double x0, double x1, double y0, double y1, double z0, double z1)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_scalar x[24];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? x1 : x0);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? y1 : y0);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? z1 : z0);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

/* the explicit and the implicit sweeps of b by m: volumes, each one's farthest vertex from the other, the cell */
static void si_both(qaws_mesh_3d const* b, qaws_motion const* m, double tol, unsigned int res, double* ve, double* vi, double* ie, double* ei, double* cell,
	int* closed)
{
	qaws_mesh_3d *e = NULL, *i = NULL;
	qaws_sweep_desc d;
	qaws_sweep_report r;
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)tol;
	d.resolution = res;
	qaws_mesh_3d_sweep(&e, NULL, b, m, &d);
	qaws_mesh_3d_sweep_implicit(&i, &r, b, m, &d);
	*ve = si_volume(e);
	*vi = si_volume(i);
	*ie = si_far(i, e);
	*ei = si_far(e, i);
	*closed = r.closed;
	*cell = 1;
	if (e)
	{
		/* (the explicit sweep's box: about the grid's) */
		unsigned int n = qaws_mesh_3d_get_vertex_count(e), j, k;
		double lo[3] = { 1e30, 1e30, 1e30 }, hi[3] = { -1e30, -1e30, -1e30 }, ext = 0;
		for (j = 0; j < n; j++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(e, j);
			double q[3] = { p.x, p.y, p.z };
			for (k = 0; k < 3; k++)
			{
				lo[k] = fmin(lo[k], q[k]);
				hi[k] = fmax(hi[k], q[k]);
			}
		}
		for (k = 0; k < 3; k++)
			ext = fmax(ext, hi[k] - lo[k]);
		*cell = ext / res;
	}
	qaws_mesh_3d_destroy(e);
	qaws_mesh_3d_destroy(i);
}

static void test_rigid(void)
{
	qaws_mesh_3d *c = si_box(0, 1, 0, 1, 0, 1), *b = si_box(1, 2, -0.1, 0.1, 0, 0.6);
	qaws_motion *mt = NULL, *mr = NULL;
	qaws_motion_pose s, e;
	double ve, vi, ie, ei, cell, ve2, vi2, ie2, ei2, cell2;
	int cl, cl2;
	char msg[400];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = e.position.y = 1;
	s.velocity = e.velocity = e.position;
	qaws_motion_create_poses(&mt, &s, &e, 0);
	si_both(c, mt, 1e-3, 32, &ve, &vi, &ie, &ei, &cell, &cl);
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
	e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	qaws_motion_create_poses(&mr, &s, &e, 0);
	si_both(b, mr, 1e-2, 32, &ve2, &vi2, &ie2, &ei2, &cell2, &cl2);
	sprintf(msg, "a cube translated: volume %.4f (3), apart %.4f and %.4f (cell %.4f), closed %d; a box turned a quarter: volumes %.4f and %.4f, apart %.4f and %.4f (cell %.4f), closed %d",
		vi, ie, ei, cell, cl, ve2, vi2, ie2, ei2, cell2, cl2);
	printf("  %s\n", msg);
	TEST_ASSERT(cl && fabs(vi - 3) < 0.06 && ie < cell && ei < cell && cl2 && fabs(vi2 - ve2) < 0.02 * ve2 && ie2 < cell2 && ei2 < cell2, msg);
	qaws_motion_destroy(mt); qaws_motion_destroy(mr);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(b);
}

static void test_vertex_motion(void)
{
	qaws_mesh_3d *b = si_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_motion *mv = NULL, *m7 = NULL;
	qaws_scalar p0[24], p1[24], v[24];
	qaws_sweep_desc d;
	qaws_sweep_report r;
	qaws_status s1, s2;
	unsigned int i;
	double vo;
	char msg[300];
	qaws_mesh_3d_get_positions(p0, b);
	memcpy(p1, p0, sizeof p1);
	for (i = 4; i < 8; i++)
		p1[3 * i + 2] = 2;
	for (i = 0; i < 24; i++)
		v[i] = p1[i] - p0[i];
	qaws_motion_create_vertices(&mv, p0, p1, v, v, NULL, NULL, 8, 0);
	qaws_motion_create_vertices(&m7, p0, p1, v, v, NULL, NULL, 7, 0);
	memset(&d, 0, sizeof d);
	d.resolution = 24;
	qaws_mesh_3d_sweep_implicit(&o, &r, b, mv, &d);
	vo = si_volume(o);
	qaws_mesh_3d_destroy(o);
	o = NULL;
	s1 = qaws_mesh_3d_sweep_implicit(&o, NULL, NULL, mv, &d);
	s2 = qaws_mesh_3d_sweep_implicit(&o, NULL, b, m7, &d);
	sprintf(msg, "a box's top rising: volume %.4f (2), %u samples, closed %d; refusals %d %d", vo, r.slabs, r.closed, (int)s1, (int)s2);
	printf("  %s\n", msg);
	TEST_ASSERT(r.closed && fabs(vo - 2) < 0.06 && s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && !o, msg);
	qaws_motion_destroy(mv); qaws_motion_destroy(m7);
	qaws_mesh_3d_destroy(b);
}

int test_184_sweep_implicit_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 184: implicit sweeps\n");
	test_rigid();
	test_vertex_motion();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
