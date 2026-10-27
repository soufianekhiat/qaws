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
 *   - a vertex motion (a box's top rising into a ceiling): it touches at 0.5.
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
	s1 = qaws_mesh_3d_first_contact(&r, b, mb, slab, &d);
	s2 = qaws_mesh_3d_first_contact(&r2, c, mc, wall, &d);
	s3 = qaws_mesh_3d_first_contact(&r3, c, mm, wall, &d);
	sprintf(msg, "a ball on a slab: contact in [%.6f, %.6f] (0.5) at (%.4f %.4f %.4f); a cube through a thin wall: [%.6f, %.6f] (0.19); beside it: contact %d",
		(double)r.time, (double)r.time_upper, (double)r.point.x, (double)r.point.y, (double)r.point.z, (double)r2.time, (double)r2.time_upper, r3.contact);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_OK && r.contact && (double)r.time <= 0.5 + 1e-9 && (double)r.time_upper >= 0.5 - 1e-9 &&
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
	qaws_mesh_3d_first_contact(&r, rod, mr, wall, &d);
	/* the rod's corner (1, 0.01) at the interval's ends: under the wall, then on it */
	qaws_motion_evaluate(&at, mr, r.time);
	ylo = (double)at.rotation[3] + 0.01 * (double)at.rotation[4];
	qaws_motion_evaluate(&at, mr, r.time_upper);
	yhi = (double)at.rotation[3] + 0.01 * (double)at.rotation[4];
	qaws_mesh_3d_first_contact(&r2, sm, ms, big, &d);
	qaws_mesh_3d_get_positions(p0, bx);
	memcpy(p1, p0, sizeof p1);
	for (i = 4; i < 8; i++)
		p1[3 * i + 2] = (qaws_scalar)1.2;
	for (i = 0; i < 24; i++)
		v[i] = p1[i] - p0[i];
	qaws_motion_create_vertices(&mv, p0, p1, v, v, NULL, NULL, 8, 0);
	qaws_mesh_3d_first_contact(&r3, bx, mv, ceil, &d);
	sprintf(msg, "a rod turning into a wall: [%.6f, %.6f], its corner at y %.6f then %.6f (the wall at 0.5); inside at 0: contact %d at %.1f; a deforming box: [%.6f, %.6f] (0.5)",
		(double)r.time, (double)r.time_upper, ylo, yhi, r2.contact, (double)r2.time, (double)r3.time, (double)r3.time_upper);
	printf("  %s\n", msg);
	TEST_ASSERT(r.contact && ylo <= 0.5 + 1e-9 && yhi >= 0.5 - 1e-4 && r2.contact && r2.time == 0 && r3.contact && (double)r3.time <= 0.5 + 1e-9 &&
		(double)r3.time_upper >= 0.5 - 1e-4, msg);
	qaws_motion_destroy(mr); qaws_motion_destroy(ms); qaws_motion_destroy(mv);
	qaws_mesh_3d_destroy(rod); qaws_mesh_3d_destroy(wall); qaws_mesh_3d_destroy(sm); qaws_mesh_3d_destroy(big); qaws_mesh_3d_destroy(bx); qaws_mesh_3d_destroy(ceil);
}

int test_183_mesh_contact_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 183: continuous collision detection\n");
	test_translations();
	test_more();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
