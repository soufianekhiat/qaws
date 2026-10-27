/*
 * Test 183: continuous collision detection (qaws_mesh_3d_first_contact), small meshes, tolerance 1e-4
 *
 *   - a ball of radius 0.5 falling from 2 by 3 onto a slab at 0: it touches at 0.5 (its pole), the time's interval
 *     holding it, the point near the origin;
 *   - a cube of side 0.2 shot 10 along x through a wall 0.001 thick: it touches at 0.19 (not passed through), and
 *     passing beside it: no contact;
 *   - a rod turning a quarter about z toward a wall: its corner on the wall at the interval's end (the motion's own
 *     rotation, rational in t);
 *   - a small cube inside a big one: contact at 0;
 *   - a vertex motion (a box's top rising into a ceiling): it touches at 0.5;
 *   - straight moves exactly: the ball at 0.5, the cube at the wall at 0.19 (the double nearest 19/100), two cubes
 *     flying at each other at 0.7, the ball and a slab rising to meet it at 0.375, all marked exact;
 *   - a rod turning toward a wall moving toward it (not straight): its corner on the wall within the tolerance at the
 *     interval's end;
 *   - a mesh against itself (qaws_mesh_3d_first_self_contact): two boxes of one mesh, one moved 1 toward the other by
 *     a vertex motion, meeting at 0.8; a box turning (rigid): never; a strip whose vertices go straight from flat to
 *     a ring: x(s, t) = (1 - t) s + t r sin(s / r) is a graph over x while dx/ds = 1 - t + t cos(s / r) > 0, so before
 *     0.5, and pinches at s = pi r at 0.5: the first contact just past 0.5, between triangles sharing no vertex.
 */

#include "test_common.h"
#include "qaws_mesh_sweep.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_mesh_3d* ct_box(double x0, double x1, double y0, double y1, double z0, double z1)
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

/* a UV ball of radius r at the origin, poles on z */
static qaws_mesh_3d* ct_ball(double r, unsigned int nu, unsigned int nv)
{
	unsigned int n = nu * (nv - 1) + 2, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = 3.14159265358979323846 * j / nv, ph = 2 * 3.14159265358979323846 * i / nu;
			unsigned int q = (j - 1) * nu + i;
			x[3 * q] = (qaws_scalar)(r * sin(th) * cos(ph));
			x[3 * q + 1] = (qaws_scalar)(r * sin(th) * sin(ph));
			x[3 * q + 2] = (qaws_scalar)(r * cos(th));
		}
	x[3 * (n - 2)] = x[3 * (n - 2) + 1] = 0;
	x[3 * (n - 2) + 2] = (qaws_scalar)r;
	x[3 * (n - 1)] = x[3 * (n - 1) + 1] = 0;
	x[3 * (n - 1) + 2] = (qaws_scalar)-r;
	for (i = 0; i < nu; i++)
	{
		unsigned int i1 = (i + 1) % nu;
		f[k++] = n - 2; f[k++] = i; f[k++] = i1;
		for (j = 1; j + 1 < nv; j++)
		{
			unsigned int a = (j - 1) * nu + i, b = (j - 1) * nu + i1, c = j * nu + i1, e = j * nu + i;
			f[k++] = a; f[k++] = e; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = b;
		}
		f[k++] = n - 1; f[k++] = (nv - 2) * nu + i1; f[k++] = (nv - 2) * nu + i;
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = n; d.face_vertices = f; d.face_count = k / 3;
	qaws_mesh_3d_create(&m, &d);
	free(x); free(f);
	return m;
}

static qaws_motion* ct_line(double x0, double y0, double z0, double dx, double dy, double dz)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	s.position.x = (qaws_scalar)x0; s.position.y = (qaws_scalar)y0; s.position.z = (qaws_scalar)z0;
	e.position.x = (qaws_scalar)(x0 + dx); e.position.y = (qaws_scalar)(y0 + dy); e.position.z = (qaws_scalar)(z0 + dz);
	s.velocity.x = e.velocity.x = (qaws_scalar)dx;
	s.velocity.y = e.velocity.y = (qaws_scalar)dy;
	s.velocity.z = e.velocity.z = (qaws_scalar)dz;
	qaws_motion_create_poses(&m, &s, &e, 0);
	return m;
}

static void test_translations(void)
{
	qaws_contact_report r, r2, r3;
	qaws_sweep_desc d;
	qaws_mesh_3d *b = ct_ball(0.5, 24, 12), *slab = ct_box(-2, 2, -2, 2, -1, 0), *c = ct_box(-0.1, 0.1, -0.1, 0.1, -0.1, 0.1), *wall = ct_box(2, 2.001, -1, 1, -1, 1);
	qaws_motion *mb = ct_line(0, 0, 2, 0, 0, -3), *mc = ct_line(0, 0, 0, 10, 0, 0), *mm = ct_line(0, 2, 0, 10, 0, 0);
	qaws_status s1, s2, s3;
	char msg[400];
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)1e-4;
	s1 = qaws_mesh_3d_first_contact(&r, b, mb, slab, NULL, &d);
	s2 = qaws_mesh_3d_first_contact(&r2, c, mc, wall, NULL, &d);
	s3 = qaws_mesh_3d_first_contact(&r3, c, mm, wall, NULL, &d);
	sprintf(msg, "a ball on a slab: contact in [%.6f, %.6f] (0.5) at (%.4f %.4f %.4f); a cube through a thin wall: [%.6f, %.6f] (0.19); beside it: contact %d",
		(double)r.time, (double)r.time_upper, (double)r.point.x, (double)r.point.y, (double)r.point.z, (double)r2.time, (double)r2.time_upper, r3.contact);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_OK && r.contact && r.exact && (double)r.time == 0.5 && r2.exact && (double)r2.time == 0.19 &&
		(double)r.time <= 0.5 + 1e-9 && (double)r.time_upper >= 0.5 - 1e-9 &&
		(double)r.time_upper - (double)r.time < 1e-4 && fabs((double)r.point.x) + fabs((double)r.point.y) + fabs((double)r.point.z) < 1e-3 && r2.contact &&
		(double)r2.time <= 0.19 + 1e-6 && (double)r2.time_upper >= 0.19 - 1e-6 && !r3.contact, msg);
	qaws_motion_destroy(mb); qaws_motion_destroy(mc); qaws_motion_destroy(mm);
	qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(slab); qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(wall);
}

static void test_more(void)
{
	qaws_contact_report r, r2, r3;
	qaws_sweep_desc d;
	qaws_mesh_3d *rod = ct_box(0, 1, -0.01, 0.01, -0.01, 0.01), *wall = ct_box(-2, 2, 0.5, 0.6, -1, 1), *sm = ct_box(-0.1, 0.1, -0.1, 0.1, -0.1, 0.1),
		*big = ct_box(-1, 1, -1, 1, -1, 1), *bx = ct_box(-0.2, 0.2, -0.2, 0.2, 0, 0.2), *ceil = ct_box(-1, 1, -1, 1, 0.7, 0.8);
	qaws_motion_pose s, e;
	qaws_motion *mr = NULL, *ms = ct_line(0, 0, 0, 0.1, 0, 0), *mv = NULL;
	qaws_motion_state at;
	qaws_scalar p0[24], p1[24], v[24];
	double ylo = 0, yhi = 0;
	unsigned int i;
	char msg[400];
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)1e-4;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
	e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	qaws_motion_create_poses(&mr, &s, &e, 0);
	qaws_mesh_3d_first_contact(&r, rod, mr, wall, NULL, &d);
	/* the rod's corner (1, 0.01) at the interval's ends: under the wall, then on it */
	qaws_motion_evaluate(&at, mr, r.time);
	ylo = (double)at.rotation[3] + 0.01 * (double)at.rotation[4];
	qaws_motion_evaluate(&at, mr, r.time_upper);
	yhi = (double)at.rotation[3] + 0.01 * (double)at.rotation[4];
	qaws_mesh_3d_first_contact(&r2, sm, ms, big, NULL, &d);
	qaws_mesh_3d_get_positions(p0, bx);
	memcpy(p1, p0, sizeof p1);
	for (i = 4; i < 8; i++)
		p1[3 * i + 2] = (qaws_scalar)1.2;
	for (i = 0; i < 24; i++)
		v[i] = p1[i] - p0[i];
	qaws_motion_create_vertices(&mv, p0, p1, v, v, NULL, NULL, 8, 0);
	qaws_mesh_3d_first_contact(&r3, bx, mv, ceil, NULL, &d);
	sprintf(msg, "a rod turning into a wall: [%.6f, %.6f], its corner at y %.6f then %.6f (the wall at 0.5); inside at 0: contact %d at %.1f; a deforming box: [%.6f, %.6f] (0.5)",
		(double)r.time, (double)r.time_upper, ylo, yhi, r2.contact, (double)r2.time, (double)r3.time, (double)r3.time_upper);
	printf("  %s\n", msg);
	TEST_ASSERT(r.contact && ylo <= 0.5 + 1e-9 && yhi >= 0.5 - 1e-4 && r2.contact && r2.time == 0 && r3.contact && (double)r3.time <= 0.5 + 1e-9 &&
		(double)r3.time_upper >= 0.5 - 1e-4, msg);
	qaws_motion_destroy(mr); qaws_motion_destroy(ms); qaws_motion_destroy(mv);
	qaws_mesh_3d_destroy(rod); qaws_mesh_3d_destroy(wall); qaws_mesh_3d_destroy(sm); qaws_mesh_3d_destroy(big); qaws_mesh_3d_destroy(bx); qaws_mesh_3d_destroy(ceil);
}

static void test_both_moving(void)
{
	qaws_contact_report r, r2, r3;
	qaws_sweep_desc d;
	qaws_mesh_3d *a = ct_box(-0.1, 0.1, -0.1, 0.1, -0.1, 0.1), *b = ct_box(2.9, 3.1, -0.1, 0.1, -0.1, 0.1), *ball = ct_ball(0.5, 24, 12), *slab = ct_box(-2, 2, -2, 2, -1, 0),
		*rod = ct_box(0, 1, -0.01, 0.01, -0.01, 0.01), *wall = ct_box(-2, 2, 0.8, 0.9, -1, 1);
	qaws_motion *ma = ct_line(0, 0, 0, 2, 0, 0), *mb = ct_line(0, 0, 0, -2, 0, 0), *mball = ct_line(0, 0, 2, 0, 0, -3), *mslab = ct_line(0, 0, 0, 0, 0, 1),
		*mrod = NULL, *mwall = ct_line(0, 0, 0, 0, -0.4, 0);
	qaws_motion_pose s, e;
	qaws_motion_state at;
	double gap = 1;
	char msg[400];
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)1e-4;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
	e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	qaws_motion_create_poses(&mrod, &s, &e, 0);
	qaws_mesh_3d_first_contact(&r, a, ma, b, mb, &d);
	qaws_mesh_3d_first_contact(&r2, ball, mball, slab, mslab, &d);
	qaws_mesh_3d_first_contact(&r3, rod, mrod, wall, mwall, &d);
	if (r3.contact)
	{
		/* the rod's corner (1, 0.01) against the wall's face y = 0.8 - 0.4 t at the interval's end */
		qaws_motion_evaluate(&at, mrod, r3.time_upper);
		gap = (double)at.rotation[3] + 0.01 * (double)at.rotation[4] - (0.8 - 0.4 * (double)r3.time_upper);
	}
	sprintf(msg, "both moving: two cubes at %.17g (exact %d, 0.7), a ball and a rising slab at %.17g (exact %d, 0.375); a rod and a moving wall in [%.6f, %.6f], its corner past the wall by %.1e",
		(double)r.time, r.exact, (double)r2.time, r2.exact, (double)r3.time, (double)r3.time_upper, gap);
	printf("  %s\n", msg);
	TEST_ASSERT(r.contact && r.exact && (double)r.time == 0.7 && r2.contact && r2.exact && (double)r2.time == 0.375 && r3.contact && !r3.exact && gap >= -1e-9 && gap <= 1e-4, msg);
	qaws_motion_destroy(ma); qaws_motion_destroy(mb); qaws_motion_destroy(mball); qaws_motion_destroy(mslab); qaws_motion_destroy(mrod); qaws_motion_destroy(mwall);
	qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(ball); qaws_mesh_3d_destroy(slab); qaws_mesh_3d_destroy(rod); qaws_mesh_3d_destroy(wall);
}

/* a strip of n squares along x (length L, width 0.2): its vertices in x (2 (n + 1)), its mesh */
static qaws_mesh_3d* ct_strip(unsigned int n, double L, qaws_scalar* x)
{
	unsigned int i, *f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i <= n; i++)
	{
		x[6 * i] = x[6 * i + 3] = (qaws_scalar)(L * i / n);
		x[6 * i + 1] = 0;
		x[6 * i + 4] = (qaws_scalar)0.2;
		x[6 * i + 2] = x[6 * i + 5] = 0;
	}
	for (i = 0; i < n; i++)
	{
		unsigned int a = 2 * i, b = 2 * i + 2;
		f[6 * i] = a; f[6 * i + 1] = b; f[6 * i + 2] = b + 1;
		f[6 * i + 3] = a; f[6 * i + 4] = b + 1; f[6 * i + 5] = a + 1;
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 2 * (n + 1); d.face_vertices = f; d.face_count = 2 * n;
	qaws_mesh_3d_create(&m, &d);
	free(f);
	return m;
}

static void test_self(void)
{
	static unsigned int const bf[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_contact_report r, r2, r3;
	qaws_sweep_desc d;
	qaws_mesh_3d_desc md;
	qaws_mesh_3d *two = NULL, *box = ct_box(0, 1, 0, 0.3, 0, 0.3), *strip;
	qaws_motion *mt = NULL, *mr = NULL, *mc = NULL;
	qaws_motion_pose s, e;
	qaws_scalar p0[48], p1[48], v[48], *c0, *c1, *cv;
	unsigned int f[72], i, k, n = 48, shared = 0;
	double L = 2, rad = 2 / (2 * 3.14159265358979323846 * 1.15);
	char msg[400];
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)1e-4;
	/* two boxes in one mesh: [-0.1, 0.1]^3 and [0.9, 1.1] x [-0.1, 0.1]^2, the second moving -1 along x */
	for (i = 0; i < 16; i++)
	{
		unsigned int j = i % 8;
		p0[3 * i] = (qaws_scalar)(((j == 1 || j == 2 || j == 5 || j == 6) ? 0.1 : -0.1) + (i >= 8 ? 1 : 0));
		p0[3 * i + 1] = (qaws_scalar)((j == 2 || j == 3 || j == 6 || j == 7) ? 0.1 : -0.1);
		p0[3 * i + 2] = (qaws_scalar)(j >= 4 ? 0.1 : -0.1);
	}
	for (i = 0; i < 72; i++)
		f[i] = bf[i % 36] + (i >= 36 ? 8 : 0);
	memcpy(p1, p0, sizeof p1);
	for (i = 8; i < 16; i++)
		p1[3 * i] -= 1;
	for (i = 0; i < 48; i++)
		v[i] = p1[i] - p0[i];
	memset(&md, 0, sizeof md);
	md.positions = p0; md.vertex_count = 16; md.face_vertices = f; md.face_count = 24;
	qaws_mesh_3d_create(&two, &md);
	qaws_motion_create_vertices(&mt, p0, p1, v, v, NULL, NULL, 16, 0);
	qaws_mesh_3d_first_self_contact(&r, two, mt, &d);
	/* a long box turning a quarter: rigid, never meets itself */
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
	e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	qaws_motion_create_poses(&mr, &s, &e, 0);
	qaws_mesh_3d_first_self_contact(&r2, box, mr, &d);
	/* a strip curled onto a circle 1.15 times shorter than it, each vertex straight from flat to its place on the circle:
	   it pinches at 0.5 */
	c0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (n + 1));
	c1 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (n + 1));
	cv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (n + 1));
	strip = ct_strip(n, L, c0);
	for (i = 0; i < 2 * (n + 1); i++)
	{
		double a = (double)c0[3 * i] / rad;
		c1[3 * i] = (qaws_scalar)(rad * sin(a));
		c1[3 * i + 1] = c0[3 * i + 1];
		c1[3 * i + 2] = (qaws_scalar)(rad * (1 - cos(a)));
		for (k = 0; k < 3; k++)
			cv[3 * i + k] = c1[3 * i + k] - c0[3 * i + k];
	}
	qaws_motion_create_vertices(&mc, c0, c1, cv, cv, NULL, NULL, 2 * (n + 1), 0);
	qaws_mesh_3d_first_self_contact(&r3, strip, mc, &d);
	if (r3.contact)
	{
		unsigned int t[6];
		qaws_mesh_3d_get_triangle(t, strip, r3.face);
		qaws_mesh_3d_get_triangle(t + 3, strip, r3.obstacle_face);
		for (i = 0; i < 9; i++)
			shared |= t[i / 3] == t[3 + i % 3];
	}
	sprintf(msg, "itself: two boxes of one mesh in [%.6f, %.6f] (0.8); a turning box: contact %d; a strip pinching (at 0.5): contact %d in [%.6f, %.6f], triangles %u and %u (shared vertex %u) at (%.4f, %.4f) in x, z (the cusp at (0.4348, 0.2768))",
		(double)r.time, (double)r.time_upper, r2.contact, r3.contact, (double)r3.time, (double)r3.time_upper, r3.face, r3.obstacle_face, shared, (double)r3.point.x, (double)r3.point.z);
	printf("  %s\n", msg);
	TEST_ASSERT(r.contact && (double)r.time <= 0.8 + 1e-9 && (double)r.time_upper >= 0.8 - 1e-4 && !r2.contact && r3.contact && (double)r3.time >= 0.5 - 1e-9 && (double)r3.time < 0.505 &&
		!shared && fabs((double)r3.point.x - 1 / 2.3) < 0.01 && fabs((double)r3.point.z - rad) < 0.01, msg);
	free(c0); free(c1); free(cv);
	qaws_motion_destroy(mt); qaws_motion_destroy(mr); qaws_motion_destroy(mc);
	qaws_mesh_3d_destroy(two); qaws_mesh_3d_destroy(box); qaws_mesh_3d_destroy(strip);
}

int test_183_mesh_contact_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 183: continuous collision detection\n");
	test_translations();
	test_more();
	test_both_moving();
	test_self();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
