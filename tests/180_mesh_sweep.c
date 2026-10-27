/*
 * Test 180: swept solids (qaws_mesh_sweep.h), small meshes
 *
 *   - a unit cube translated by (1, 1, 0): its Minkowski sum with the
 *     segment, volume 3 exactly, one slab;
 *   - a thin box turned a quarter about the z axis: the annular sector's
 *     volume (pi / 4 (r1^2 - r0^2) times its height, and the boxes at both
 *     ends half outside it: 2 h), within the chords; ten times thicker (0.2)
 *     and 0.6 high at tolerance 0.002 (32 slabs: its hulls' walls within
 *     1e-5 of one another at the rim) closed and one piece; with simplify the same
 *     volume (to round-off) on fewer faces, still closed;
 *   - a small ball along a Hermite curve: a tube, its volume the ball's
 *     cross-section along the length plus the ball (within the polyhedra's
 *     own error);
 *   - carving (qaws_mesh_sweep_boolean): a cube tool through a block's top
 *     along x, the slot removed exactly (0.5); along an L of two moves, the
 *     removed volume by inclusion and exclusion (0.375), and the intersection
 *     (a separate Boolean) the same volume; stock = result + removed;
 *   - a spin: a box beside the z axis turning two turns about it (one
 *     motion, its turns counted from the angular velocities), drifting
 *     under the tolerance: found a spin (2 turns), swept as its solid of
 *     revolution, the annulus between the box's nearest and farthest
 *     distances from the axis (volume within the chords' 2 tol of the
 *     circumference), closed; a cube on the axis turning two turns: a
 *     cylinder of few faces (the profile's straight sides one edge each), its
 *     NORMAL corner values the faces' (a crease tracing the surface has no
 *     normal of it: the face's own); a 64-gon cylinder drifting 100
 *     tolerances a turn, symmetric about the axis, still its revolution along
 *     the drift (the disc swept along the segment); the cube drifting 10
 *     tolerances a turn with revolve_spins: revolved anyway, its turn's drift
 *     reported, its volume between the true sweep's and the cylinder's along
 *     the drift;
 *   - a vertex motion of another vertex count refused.
 */

#include "test_common.h"
#include "qaws_mesh_sweep.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SW_PI 3.14159265358979323846

static double sw_volume(qaws_mesh_3d const* m)
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

/* the box [x0, x1] x [y0, y1] x [z0, z1], outward */
static qaws_mesh_3d* sw_box(double x0, double x1, double y0, double y1, double z0, double z1)
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

/* an icosahedron split once, on the sphere of radius r: 80 faces */
static qaws_mesh_3d* sw_ball(double r)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2, x[3 * 42];
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int t[3 * 80], nv = 12, nt = 0, i, j, mid[30][3], nm = 0;
	qaws_scalar xs[3 * 42];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 36; i++) x[i] = v0[i];
	for (i = 0; i < 20; i++)
	{
		unsigned int c[3], e;
		for (e = 0; e < 3; e++)
		{
			unsigned int a = (unsigned int)ic[3 * i + e], b = (unsigned int)ic[3 * i + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, s;
			for (s = 0; s < nm && !(mid[s][0] == lo && mid[s][1] == hi); s++)
				;
			if (s == nm)
			{
				mid[nm][0] = lo; mid[nm][1] = hi; mid[nm][2] = nv;
				for (j = 0; j < 3; j++) x[3 * nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
				nv++;
				nm++;
			}
			c[e] = mid[s][2];
		}
		t[3 * nt] = (unsigned int)ic[3 * i]; t[3 * nt + 1] = c[0]; t[3 * nt + 2] = c[2]; nt++;
		t[3 * nt] = (unsigned int)ic[3 * i + 1]; t[3 * nt + 1] = c[1]; t[3 * nt + 2] = c[0]; nt++;
		t[3 * nt] = (unsigned int)ic[3 * i + 2]; t[3 * nt + 1] = c[2]; t[3 * nt + 2] = c[1]; nt++;
		t[3 * nt] = c[0]; t[3 * nt + 1] = c[1]; t[3 * nt + 2] = c[2]; nt++;
	}
	for (i = 0; i < nv; i++)
	{
		double l = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		for (j = 0; j < 3; j++) xs[3 * i + j] = (qaws_scalar)(r * x[3 * i + j] / l);
	}
	memset(&d, 0, sizeof d);
	d.positions = xs; d.vertex_count = nv; d.face_vertices = t; d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static qaws_motion* sw_line(double dx, double dy, double dz)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = (qaws_scalar)dx; e.position.y = (qaws_scalar)dy; e.position.z = (qaws_scalar)dz;
	s.velocity = e.velocity = e.position;
	qaws_motion_create_poses(&m, &s, &e, 0);
	return m;
}

static void test_translation(void)
{
	qaws_mesh_3d *c = sw_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_motion* m = sw_line(1, 1, 0);
	qaws_sweep_report r;
	qaws_status st = qaws_mesh_3d_sweep(&o, &r, c, m, NULL);
	double v = sw_volume(o);
	char msg[300];
	sprintf(msg, "a unit cube moved by (1, 1, 0): volume %.15f (3), %u slab, %u front faces, closed %d", v, r.slabs, r.front_faces, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(v - 3) < 1e-12 && r.slabs == 1 && r.closed, msg);
	qaws_mesh_3d_destroy(c);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
}

static void test_turn(void)
{
	double h = 0.01, r0 = 1, r1 = sqrt(4 + h * h), want = 0.25 * SW_PI * (r1 * r1 - r0 * r0) + 2 * h, v;   /* (the sector, and the boxes at both ends half outside it) */
	qaws_mesh_3d *b = sw_box(1, 2, -h, h, 0, 1), *o = NULL;
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_report r;
	qaws_status st;
	char msg[300];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(SW_PI / 4);   /* (a quarter turn: the half angle) */
	e.rotation[3] = (qaws_scalar)sin(SW_PI / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(SW_PI / 2);
	qaws_motion_create_poses(&m, &s, &e, 0);
	st = qaws_mesh_3d_sweep(&o, &r, b, m, NULL);
	v = sw_volume(o);
	sprintf(msg, "a thin box turned a quarter about z: volume %.5f (the annular sector %.5f), %u slabs, chords within %.1e, closed %d", v, want, r.slabs,
		(double)r.deviation, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(v - want) < 0.01 * want && r.slabs > 4 && r.closed, msg);
	{
		qaws_sweep_desc d;
		qaws_mesh_3d* so = NULL;
		qaws_sweep_report rs;
		double vs;
		memset(&d, 0, sizeof d);
		d.simplify = 1;
		st = qaws_mesh_3d_sweep(&so, &rs, b, m, &d);
		vs = sw_volume(so);
		sprintf(msg, "simplified: volume %.12f (%.12f), %u faces (%u), closed %d", vs, v, so ? qaws_mesh_3d_get_face_count(so) : 0, o ? qaws_mesh_3d_get_face_count(o) : 0, rs.closed);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && so && o && fabs(vs - v) < 1e-9 * v && qaws_mesh_3d_get_face_count(so) < qaws_mesh_3d_get_face_count(o) && rs.closed, msg);
		qaws_mesh_3d_destroy(so);
	}
	{
		/* ten times thicker, 32 slabs: sheets nearer than the Boolean's samples (each sample now brought within 0.4 of the first surface its ray meets) */
		double hh = 0.1, rr1 = sqrt(4 + hh * hh), wt = 0.6 * (0.25 * SW_PI * (rr1 * rr1 - 1) + 2 * hh), vt;
		qaws_mesh_3d *bt = sw_box(1, 2, -hh, hh, 0, 0.6), *ot = NULL;
		qaws_sweep_desc dt;
		qaws_sweep_report rt;
		qaws_mesh_topology tp;
		memset(&dt, 0, sizeof dt);
		dt.tolerance = (qaws_scalar)0.002;
		memset(&tp, 0, sizeof tp);
		st = qaws_mesh_3d_sweep(&ot, &rt, bt, m, &dt);
		vt = sw_volume(ot);
		if (ot)
			qaws_mesh_3d_get_topology(&tp, ot);
		sprintf(msg, "a box 0.2 thick turned a quarter: %u slabs, closed %d, %u piece(s), volume %.5f (the sector %.5f)", rt.slabs, tp.closed, tp.components, vt, wt);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && ot && rt.slabs >= 32 && tp.closed && tp.components == 1 && fabs(vt - wt) < 0.01 * wt, msg);
		qaws_mesh_3d_destroy(bt);
		qaws_mesh_3d_destroy(ot);
	}
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
}

static void test_tube(void)
{
	double rad = 0.2, L = 0, vb, v, want;
	qaws_mesh_3d *b = sw_ball(rad), *o = NULL;
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_report r;
	qaws_sweep_desc d;
	qaws_status st;
	unsigned int k;
	char msg[300];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	e.position.x = 2;
	s.velocity.x = 1; s.velocity.y = 2;
	e.velocity.x = 1; e.velocity.y = -2;
	qaws_motion_create_poses(&m, &s, &e, 0);
	/* the path's length, by fine chords */
	{
		qaws_motion_state q0, q1;
		qaws_motion_evaluate(&q0, m, 0);
		for (k = 1; k <= 2000; k++)
		{
			qaws_motion_evaluate(&q1, m, (qaws_scalar)(k / 2000.0));
			L += sqrt(((double)q1.position.x - q0.position.x) * ((double)q1.position.x - q0.position.x) + ((double)q1.position.y - q0.position.y) * ((double)q1.position.y - q0.position.y));
			q0 = q1;
		}
	}
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.004;
	st = qaws_mesh_3d_sweep(&o, &r, b, m, &d);
	v = sw_volume(o);
	vb = sw_volume(b);
	/* (the polyhedral ball's cross-section scaled as its volume: the tube of a ball of its volume) */
	want = SW_PI * rad * rad * L * pow(vb / (4.0 / 3 * SW_PI * rad * rad * rad), 2.0 / 3) + vb;
	sprintf(msg, "a ball (80 faces, r 0.2) along a Hermite curve of length %.4f: volume %.4f (the tube %.4f), %u slabs, chords within %.1e, closed %d",
		L, v, want, r.slabs, (double)r.deviation, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(v - want) < 0.04 * want && r.closed, msg);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
}

static qaws_motion* sw_move(double x0, double y0, double z0, double x1, double y1, double z1)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	s.position.x = (qaws_scalar)x0; s.position.y = (qaws_scalar)y0; s.position.z = (qaws_scalar)z0;
	e.position.x = (qaws_scalar)x1; e.position.y = (qaws_scalar)y1; e.position.z = (qaws_scalar)z1;
	s.velocity.x = e.velocity.x = (qaws_scalar)(x1 - x0);
	s.velocity.y = e.velocity.y = (qaws_scalar)(y1 - y0);
	s.velocity.z = e.velocity.z = (qaws_scalar)(z1 - z0);
	qaws_motion_create_poses(&m, &s, &e, 0);
	return m;
}

static void test_carving(void)
{
	qaws_mesh_3d *block = sw_box(0, 4, 0, 2, 0, 1), *tool = sw_box(-0.25, 0.25, -0.25, 0.25, -0.25, 0.25), *o1 = NULL, *o2 = NULL, *o3 = NULL;
	qaws_motion *slot = sw_move(-1, 1, 1, 5, 1, 1), *la = sw_move(-1, 1, 1, 2, 1, 1), *lb = sw_move(2, 1, 1, 2, 3, 1);
	qaws_motion const* one[1];
	qaws_motion const* two[2];
	qaws_sweep_boolean_report r1, r2, r3;
	qaws_status s1, s2, s3;
	char msg[400];
	one[0] = slot;
	two[0] = la;
	two[1] = lb;
	s1 = qaws_mesh_sweep_boolean(&o1, &r1, tool, one, block, 1, QAWS_BOOLEAN_DIFFERENCE, NULL);
	s2 = qaws_mesh_sweep_boolean(&o2, &r2, tool, two, block, 2, QAWS_BOOLEAN_DIFFERENCE, NULL);
	s3 = qaws_mesh_sweep_boolean(&o3, &r3, tool, two, block, 2, QAWS_BOOLEAN_INTERSECTION, NULL);
	sprintf(msg, "carving a 4 x 2 x 1 block with a 0.5 cube: a slot along its top removing %.15f (0.5); an L of two moves removing %.15f (0.375), its intersection %.15f; stock %.6f = result %.6f + removed; closed %d %d %d",
		(double)r1.changed_volume, (double)r2.changed_volume, (double)r3.changed_volume, (double)r2.stock_volume, (double)r2.result_volume, r1.closed, r2.closed, r3.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_OK && fabs((double)r1.changed_volume - 0.5) < 1e-12 &&
		fabs((double)r2.changed_volume - 0.375) < 1e-12 && fabs((double)r3.changed_volume - 0.375) < 1e-12 && r2.moves == 2 &&
		fabs((double)r2.stock_volume - (double)r2.result_volume - (double)r2.changed_volume) < 1e-12 && r1.closed && r2.closed && r3.closed, msg);
	qaws_mesh_3d_destroy(block); qaws_mesh_3d_destroy(tool);
	qaws_mesh_3d_destroy(o1); qaws_mesh_3d_destroy(o2); qaws_mesh_3d_destroy(o3);
	qaws_motion_destroy(slot); qaws_motion_destroy(la); qaws_motion_destroy(lb);
}

static void test_spin(void)
{
	double h = 0.5, r0 = 0.5, r1 = sqrt(1 + 0.25 * 0.25), want = SW_PI * (r1 * r1 - r0 * r0) * h, v, tol = 0.002;
	qaws_mesh_3d *b = sw_box(0.5, 1, -0.25, 0.25, 0, h), *o = NULL;
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	qaws_sweep_desc d;
	qaws_sweep_report r;
	qaws_status st;
	char msg[300];
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;   /* (two whole turns: the same rotation at both ends) */
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(4 * SW_PI);
	e.position.x = (qaws_scalar)0.001;
	s.velocity.x = e.velocity.x = (qaws_scalar)0.001;
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)tol;
	memset(&r, 0, sizeof r);
	st = qaws_motion_create_poses(&m, &s, &e, 0);
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_3d_sweep(&o, &r, b, m, &d);
	v = sw_volume(o);
	sprintf(msg, "a spin (%.2f turns): volume %.5f (the annulus %.5f), %u faces, closed %d", (double)r.turns, v, want, o ? qaws_mesh_3d_get_face_count(o) : 0, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs((double)r.turns - 2) < 1e-6 && fabs(v - want) < 2 * SW_PI * r1 * h * 2 * tol && r.closed, msg);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
}

/* a 64-gon cylinder about z, radius r, height 1 */
static qaws_mesh_3d* sw_cylinder(double r)
{
	unsigned int n = 64, i, k = 0, f[3 * 256];
	qaws_scalar x[3 * 130];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < n; i++)
	{
		double u = 2 * SW_PI * i / n;
		x[3 * i] = x[3 * (n + i)] = (qaws_scalar)(r * cos(u));
		x[3 * i + 1] = x[3 * (n + i) + 1] = (qaws_scalar)(r * sin(u));
		x[3 * i + 2] = 0;
		x[3 * (n + i) + 2] = 1;
	}
	x[6 * n] = x[6 * n + 1] = x[6 * n + 2] = 0;
	x[6 * n + 3] = x[6 * n + 4] = 0;
	x[6 * n + 5] = 1;
	for (i = 0; i < n; i++)
	{
		unsigned int a = i, b = (i + 1) % n, c = n + (i + 1) % n, e = n + i;
		f[k++] = a; f[k++] = b; f[k++] = c;
		f[k++] = a; f[k++] = c; f[k++] = e;
		f[k++] = 2 * n; f[k++] = b; f[k++] = a;
		f[k++] = 2 * n + 1; f[k++] = e; f[k++] = c;
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 2 * n + 2; d.face_vertices = f; d.face_count = 4 * n;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

/* two turns about z drifting dx along x */
static qaws_motion* sw_spin(double dx)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(4 * SW_PI);
	e.position.x = (qaws_scalar)dx;
	s.velocity.x = e.velocity.x = (qaws_scalar)dx;
	qaws_motion_create_poses(&m, &s, &e, 0);
	return m;
}

static void test_spin_more(void)
{
	double tol = 0.002, r = 0.5, dd = 0.4, area = 32 * r * r * sin(2 * SW_PI / 64), want = area + 2 * r * dd, v, least = 2;
	qaws_mesh_3d *cube = sw_box(-0.5, 0.5, -0.5, 0.5, 0, 1), *cyl = sw_cylinder(r), *o = NULL, *oc = NULL;
	qaws_motion *m0 = sw_spin(0.001), *m1 = sw_spin(dd);
	qaws_sweep_desc d;
	qaws_sweep_report r0, r1, r2;
	qaws_motion* m2 = sw_spin(0.04);
	qaws_mesh_3d* o2 = NULL;
	double v2;
	unsigned int nt = qaws_mesh_3d_get_triangle_count(cube), i, k, c[3], id, nf = 0;
	qaws_scalar* nn = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * nt);
	qaws_attribute_desc ad;
	char msg[400];
	/* the cube's faces' normals per corner */
	for (i = 0; i < nt; i++)
	{
		qaws_vec3 a, b, e;
		double u[3], w[3], n[3], l;
		qaws_mesh_3d_get_triangle(c, cube, i);
		a = qaws_mesh_3d_get_vertex(cube, c[0]); b = qaws_mesh_3d_get_vertex(cube, c[1]); e = qaws_mesh_3d_get_vertex(cube, c[2]);
		u[0] = b.x - a.x; u[1] = b.y - a.y; u[2] = b.z - a.z;
		w[0] = e.x - a.x; w[1] = e.y - a.y; w[2] = e.z - a.z;
		n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (k = 0; k < 9; k++) nn[9 * i + k] = (qaws_scalar)(n[k % 3] / l);
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "n"; ad.element = QAWS_MESH_CORNER; ad.components = 3; ad.rule = QAWS_ATTRIBUTE_NORMAL;
	qaws_mesh_3d_add_attribute(&id, cube, nn, &ad);
	free(nn);
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)tol;
	memset(&r0, 0, sizeof r0);
	memset(&r1, 0, sizeof r1);
	if (qaws_mesh_3d_sweep(&o, &r0, cube, m0, &d) == QAWS_STATUS_OK && o && qaws_mesh_3d_find_attribute(&id, o, "n") == QAWS_STATUS_OK)
	{
		unsigned int no = qaws_mesh_3d_get_triangle_count(o);
		qaws_scalar* nv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * no);
		nf = no;
		qaws_mesh_3d_get_attribute(nv, o, id);
		for (i = 0; i < no; i++)
		{
			qaws_vec3 a, b, e;
			double u[3], w[3], n[3], l;
			qaws_mesh_3d_get_triangle(c, o, i);
			a = qaws_mesh_3d_get_vertex(o, c[0]); b = qaws_mesh_3d_get_vertex(o, c[1]); e = qaws_mesh_3d_get_vertex(o, c[2]);
			u[0] = b.x - a.x; u[1] = b.y - a.y; u[2] = b.z - a.z;
			w[0] = e.x - a.x; w[1] = e.y - a.y; w[2] = e.z - a.z;
			n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
			l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
			for (k = 0; k < 3; k++)
				least = fmin(least, (n[0] * nv[9 * i + 3 * k] + n[1] * nv[9 * i + 3 * k + 1] + n[2] * nv[9 * i + 3 * k + 2]) / l);
		}
		free(nv);
	}
	qaws_mesh_3d_sweep(&oc, &r1, cyl, m1, &d);
	v = sw_volume(oc);
	d.revolve_spins = 1;
	memset(&r2, 0, sizeof r2);
	qaws_mesh_3d_sweep(&o2, &r2, cube, m2, &d);
	v2 = sw_volume(o2);
	{
		/* (the cube's cylinder, radius sqrt(1/2), along 0.04: the bound; the true sweep is within a turn's drift of it) */
		double hi = SW_PI * 0.5 + 2 * sqrt(0.5) * 0.04;
		char m2s[200];
		sprintf(m2s, "revolve_spins: %.2f turns, a turn drifts %.4f, volume %.5f (at most %.5f), closed %d", (double)r2.turns, (double)r2.spin_drift, v2, hi, r2.closed);
		printf("  %s\n", m2s);
		TEST_ASSERT(o2 && fabs((double)r2.turns - 2) < 1e-4 && fabs((double)r2.spin_drift - 0.02) < 1e-4 && v2 <= hi && v2 > hi - 2 * SW_PI * sqrt(0.5) * 2 * tol && r2.closed, m2s);
	}
	qaws_mesh_3d_destroy(o2);
	qaws_motion_destroy(m2);
	sprintf(msg, "a cube's two turns: %u faces, NORMALs against the faces at least %.4f; a cylinder drifting 100 tolerances a turn (%.2f turns): volume %.5f (the disc along the segment %.5f), closed %d",
		nf, least, (double)r1.turns, v, want, r1.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(nf > 0 && nf < 400 && least > 0.999 && fabs((double)r1.turns - 2) < 1e-6 && fabs(v - want) < 2 * SW_PI * r * 2 * tol + 1e-3 && r1.closed, msg);
	qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(cyl); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(oc);
	qaws_motion_destroy(m0); qaws_motion_destroy(m1);
}

static void test_refusals(void)
{
	qaws_mesh_3d *c = sw_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_scalar p[9] = { 0 };
	qaws_motion* m = NULL;
	qaws_status st;
	char msg[200];
	qaws_motion_create_vertices(&m, p, p, p, p, NULL, NULL, 3, 0);
	st = qaws_mesh_3d_sweep(&o, NULL, c, m, NULL);
	sprintf(msg, "refused: a vertex motion of 3 vertices for an 8-vertex mesh (%d)", (int)st);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_INVALID_ARGUMENT && !o, msg);
	qaws_mesh_3d_destroy(c);
	qaws_motion_destroy(m);
}

int test_180_mesh_sweep_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 180: swept solids\n");
	test_translation();
	test_turn();
	test_tube();
	test_carving();
	test_spin();
	test_spin_more();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
