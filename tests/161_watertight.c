/*
 * Test 161: watertight surfaces (qaws_mesh_make_watertight)
 *
 *   - a sphere soup with three holes: closed, every vertex within a cell of
 *     the sphere, the volume near the ball's; only grid points near the
 *     input evaluated one by one; the holes filled first (tangent to the
 *     sphere) at under half the winding number's own caps' distance;
 *   - two cubes crossing (a self-intersecting soup): their union, closed,
 *     its vertices on the union's surface (sharp edges and corners kept by
 *     the crossings' planes), the volume the union's, no face crossing
 *     another;
 *   - an open square sheet with an offset: a closed shell round it, its
 *     vertices at the offset from the sheet;
 *   - refusals: no faces, an index out of range, a negative offset.
 */

#include "test_common.h"
#include "qaws_mesh_watertight.h"
#include "qaws_mesh_repair.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define WT_PI 3.14159265358979323846

static double wt_volume(qaws_mesh_3d const* m, double* far_out, double (*dist)(double const*))
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), i;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	double v = 0, far = 0;
	qaws_mesh_3d_get_positions(p, m);
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		double a[3], b[3], c[3];
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			a[k] = p[3 * t[3 * i] + k];
			b[k] = p[3 * t[3 * i + 1] + k];
			c[k] = p[3 * t[3 * i + 2] + k];
		}
		v += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
	}
	for (i = 0; i < nv; i++)
	{
		double q[3] = { p[3 * i], p[3 * i + 1], p[3 * i + 2] };
		far = fmax(far, dist(q));
	}
	*far_out = far;
	free(p);
	free(t);
	return v;
}

static double const wt_holes[3][3] = { { 1, 0, 0 }, { 0, 0.7071, 0.7071 }, { -0.6, -0.8, 0 } };

static double wt_sphere_dist_all(double const* q)
{
	return fabs(sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) - 1);
}

/* the distance to the sphere, away from the holes (their caps close them inside it) */
static double wt_sphere_dist(double const* q)
{
	unsigned int h;
	for (h = 0; h < 3; h++)
		if ((q[0] - wt_holes[h][0]) * (q[0] - wt_holes[h][0]) + (q[1] - wt_holes[h][1]) * (q[1] - wt_holes[h][1]) + (q[2] - wt_holes[h][2]) * (q[2] - wt_holes[h][2]) < 0.45 * 0.45)
			return 0;
	return fabs(sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) - 1);
}

/* the distance to the surface of the union of [0,1]^3 and [0.4,1.4] x [0.3,1.3] x [0.2,1.2] */
static double wt_box_sd(double const* q, double const* lo, double const* hi)
{
	double out = 0, in = HUGE_VAL;
	unsigned int k;
	int inside = 1;
	for (k = 0; k < 3; k++)
	{
		double e = q[k] < lo[k] ? lo[k] - q[k] : q[k] > hi[k] ? q[k] - hi[k] : 0;
		out += e * e;
		if (e > 0)
			inside = 0;
		in = fmin(in, fmin(q[k] - lo[k], hi[k] - q[k]));
	}
	return inside ? -in : sqrt(out);
}

static double wt_union_dist(double const* q)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 0.4, 0.3, 0.2 }, hi1[3] = { 1.4, 1.3, 1.2 };
	return fabs(fmin(wt_box_sd(q, lo0, hi0), wt_box_sd(q, lo1, hi1)));
}

static double wt_sheet_dist(double const* q)
{
	/* the distance to [0,1]^2 x {0} less the offset 0.05 */
	double dx = q[0] < 0 ? -q[0] : q[0] > 1 ? q[0] - 1 : 0, dy = q[1] < 0 ? -q[1] : q[1] > 1 ? q[1] - 1 : 0;
	return fabs(sqrt(dx * dx + dy * dy + q[2] * q[2]) - 0.05);
}

/* a box [lo, hi] as 12 triangles (a soup: its own 8 points), outward */
static void wt_box(double const* lo, double const* hi, qaws_scalar* x, unsigned int* nv, unsigned int* t, unsigned int* nt)
{
	static unsigned int const bt[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5 };
	static int const bx[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
	unsigned int i, k;
	for (i = 0; i < 8; i++)
		for (k = 0; k < 3; k++)
			x[3 * (*nv + i) + k] = (qaws_scalar)(bx[3 * i + k] ? hi[k] : lo[k]);
	for (i = 0; i < 36; i++)
		t[3 * *nt + i] = *nv + bt[i];
	*nv += 8;
	*nt += 12;
}

static void test_sphere(void)
{
	unsigned int na = 40, nb = 20, nv = 2 + na * (nb - 1), i, j, t = 0, kept = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nb);
	unsigned int* k2 = (unsigned int*)malloc(sizeof(unsigned int) * 6 * na * nb);
	qaws_watertight_report r;
	qaws_mesh_3d* m = NULL;
	qaws_status s;
	double far = 0, vol;
	char msg[256];
	x[0] = 0; x[1] = 0; x[2] = 1;
	x[3] = 0; x[4] = 0; x[5] = -1;
	for (j = 1; j < nb; j++)
		for (i = 0; i < na; i++)
		{
			double th = WT_PI * j / nb, ph = 2 * WT_PI * i / na;
			unsigned int v = 2 + (j - 1) * na + i;
			x[3 * v] = (qaws_scalar)(sin(th) * cos(ph));
			x[3 * v + 1] = (qaws_scalar)(sin(th) * sin(ph));
			x[3 * v + 2] = (qaws_scalar)cos(th);
		}
	for (i = 0; i < na; i++)
	{
		f[t++] = 0; f[t++] = 2 + i; f[t++] = 2 + (i + 1) % na;
		f[t++] = 1; f[t++] = 2 + (nb - 2) * na + (i + 1) % na; f[t++] = 2 + (nb - 2) * na + i;
	}
	for (j = 1; j + 1 < nb; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = 2 + (j - 1) * na + i, b = 2 + (j - 1) * na + (i + 1) % na, c = b + na, d = a + na;
			f[t++] = a; f[t++] = d; f[t++] = c;
			f[t++] = a; f[t++] = c; f[t++] = b;
		}
	for (i = 0; i < t / 3; i++)
	{
		double g[3] = { 0, 0, 0 };
		unsigned int k, h, drop = 0;
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				g[j] += x[3 * f[3 * i + k] + j] / 3;
		for (h = 0; h < 3; h++)
			drop |= (g[0] - wt_holes[h][0]) * (g[0] - wt_holes[h][0]) + (g[1] - wt_holes[h][1]) * (g[1] - wt_holes[h][1]) + (g[2] - wt_holes[h][2]) * (g[2] - wt_holes[h][2])
				< 0.3 * 0.3;
		if (!drop)
		{
			memcpy(k2 + 3 * kept, f + 3 * i, sizeof(unsigned int) * 3);
			kept++;
		}
	}
	s = qaws_mesh_make_watertight(&m, &r, x, k2, NULL, nv, kept, NULL);
	vol = m ? wt_volume(m, &far, wt_sphere_dist) : 0;
	sprintf(msg, "sphere with holes: closed %d, %u faces, cell %.4f, farthest %.4f from the sphere (off the holes), volume %.4f (ball %.4f), %u of %u points evaluated",
		r.closed, m ? qaws_mesh_3d_get_face_count(m) : 0, (double)r.cell, far, vol, 4 * WT_PI / 3, r.evaluated, r.grid[0] * r.grid[1] * r.grid[2]);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && far < r.cell && fabs(vol - 4 * WT_PI / 3) < 0.05 * 4 * WT_PI / 3
		&& r.evaluated < r.grid[0] * r.grid[1] * r.grid[2] / 2 && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, msg);
	/* the caps: the fills (tangent to the sphere round the holes) against the winding number's own caps */
	{
		double far_fill = 0, far_cap = 0;
		qaws_watertight_desc d;
		qaws_mesh_3d* mc = NULL;
		memset(&d, 0, sizeof d);
		d.flags = QAWS_WATERTIGHT_WINDING_CAPS;
		qaws_mesh_make_watertight(&mc, NULL, x, k2, NULL, nv, kept, &d);
		if (m)
			wt_volume(m, &far_fill, wt_sphere_dist_all);
		if (mc)
			wt_volume(mc, &far_cap, wt_sphere_dist_all);
		sprintf(msg, "sphere with holes, holes included: %u filled first, farthest %.4f from the sphere (winding caps %.4f)", r.holes_filled, far_fill, far_cap);
		printf("  %s\n", msg);
		TEST_ASSERT(r.holes_filled == 3 && far_fill < 0.5 * far_cap && far_fill < 0.02, msg);
		qaws_mesh_3d_destroy(mc);
	}
	qaws_mesh_3d_destroy(m);
	free(x);
	free(f);
	free(k2);
}

static void test_union(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 0.4, 0.3, 0.2 }, hi1[3] = { 1.4, 1.3, 1.2 };
	qaws_scalar x[48];
	unsigned int t[72], nv = 0, nt = 0;
	qaws_watertight_report r;
	qaws_mesh_3d* m = NULL;
	qaws_status s;
	double far = 0, vol, want = 2 - 0.6 * 0.7 * 0.8;
	char msg[256];
	wt_box(lo0, hi0, x, &nv, t, &nt);
	wt_box(lo1, hi1, x, &nv, t, &nt);
	s = qaws_mesh_make_watertight(&m, &r, x, t, NULL, nv, nt, NULL);
	vol = m ? wt_volume(m, &far, wt_union_dist) : 0;
	sprintf(msg, "two cubes crossing: closed %d, %u faces, %u sharp cells, farthest %.2e from the union, volume %.6f (union %.6f)", r.closed,
		m ? qaws_mesh_3d_get_face_count(m) : 0, r.sharp, far, vol, want);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && far < 1e-3 && fabs(vol - want) < 1e-3 && r.sharp > 0, msg);
	/* no face of the result crosses another */
	if (m)
	{
		unsigned int n = qaws_mesh_3d_get_vertex_count(m), ntr = qaws_mesh_3d_get_triangle_count(m);
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
		unsigned int* tr = (unsigned int*)malloc(sizeof(unsigned int) * 3 * ntr);
		qaws_mesh_repair_desc d;
		qaws_mesh_repair_report rr;
		qaws_mesh_3d* chk = NULL;
		qaws_mesh_3d_get_positions(p, m);
		qaws_mesh_3d_get_triangles(tr, m);
		memset(&d, 0, sizeof d);
		d.flags = QAWS_MESH_REPAIR_FIND_INTERSECTIONS;
		d.weld_tolerance = -1;
		qaws_mesh_repair(&chk, &rr, p, tr, NULL, n, ntr, &d);
		sprintf(msg, "two cubes crossing: the union's faces cross %u times", rr.self_intersections);
		printf("  %s\n", msg);
		TEST_ASSERT(chk && rr.self_intersections == 0, msg);
		qaws_mesh_3d_destroy(chk);
		free(p);
		free(tr);
	}
	qaws_mesh_3d_destroy(m);
}

static void test_sheet(void)
{
	unsigned int n = 8, nv = (n + 1) * (n + 1), i, j, t = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_watertight_desc d;
	qaws_watertight_report r;
	qaws_mesh_3d* m = NULL;
	qaws_status s;
	double far = 0, vol;
	char msg[256];
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * (n + 1) + i)] = (qaws_scalar)((double)i / n);
			x[3 * (j * (n + 1) + i) + 1] = (qaws_scalar)((double)j / n);
			x[3 * (j * (n + 1) + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * (n + 1) + i;
			f[t++] = a; f[t++] = a + 1; f[t++] = a + n + 2;
			f[t++] = a; f[t++] = a + n + 2; f[t++] = a + n + 1;
		}
	memset(&d, 0, sizeof d);
	d.offset = (qaws_scalar)0.05;
	d.resolution = 48;
	s = qaws_mesh_make_watertight(&m, &r, x, f, NULL, nv, t / 3, &d);
	vol = m ? wt_volume(m, &far, wt_sheet_dist) : 0;
	sprintf(msg, "open sheet, offset 0.05: closed %d, %u faces, cell %.4f, farthest %.4f off the offset, volume %.4f (slab 0.1 + rim)", r.closed,
		m ? qaws_mesh_3d_get_face_count(m) : 0, (double)r.cell, far, vol);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && far < 0.5 * r.cell && vol > 0.1 && vol < 0.1 + 4 * 0.05 * 0.05 * WT_PI / 2 + 0.02, msg);
	qaws_mesh_3d_destroy(m);
	free(x);
	free(f);
}

static void test_refusals(void)
{
	qaws_scalar x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	unsigned int t[3] = { 0, 1, 2 }, bad[3] = { 0, 1, 5 };
	qaws_watertight_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	TEST_ASSERT(qaws_mesh_make_watertight(&m, NULL, x, t, NULL, 3, 0, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m, "refusals: no faces");
	TEST_ASSERT(qaws_mesh_make_watertight(&m, NULL, x, bad, NULL, 3, 1, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m, "refusals: an index out of range");
	d.offset = -1;
	TEST_ASSERT(qaws_mesh_make_watertight(&m, NULL, x, t, NULL, 3, 1, &d) == QAWS_STATUS_INVALID_ARGUMENT && !m, "refusals: a negative offset");
}

/* the enclosed volume of a closed mesh */
static double wt_points_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		qaws_vec3 a, b, c;
		qaws_mesh_3d_get_triangle(t, m, i);
		a = qaws_mesh_3d_get_vertex(m, t[0]);
		b = qaws_mesh_3d_get_vertex(m, t[1]);
		c = qaws_mesh_3d_get_vertex(m, t[2]);
		v += ((double)a.x * ((double)b.y * c.z - (double)b.z * c.y) - (double)a.y * ((double)b.x * c.z - (double)b.z * c.x) +
			(double)a.z * ((double)b.x * c.y - (double)b.y * c.x)) / 6;
	}
	return v;
}

/* oriented point clouds: a sphere (Fibonacci points), a cube sampled on its faces; normals turned inward refused */
static void test_points(void)
{
	unsigned int n = 2000, i, k, nc = 0;
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n), *nr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	qaws_scalar *cx = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 6 * 400), *cn = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 6 * 400);
	qaws_mesh_3d *s = NULL, *c = NULL, *z = NULL;
	qaws_watertight_report rs, rc;
	qaws_watertight_desc d;
	qaws_status st, sc, sz, s3;
	double vs = 0, vc = 0, rmin = 1e30, rmax = 0;
	char msg[400];
	for (i = 0; i < n; i++)
	{
		double zz = 1 - 2 * (i + 0.5) / n, r = sqrt(1 - zz * zz), ph = 2.399963229728653 * i;
		x[3 * i] = nr[3 * i] = (qaws_scalar)(r * cos(ph));
		x[3 * i + 1] = nr[3 * i + 1] = (qaws_scalar)(r * sin(ph));
		x[3 * i + 2] = nr[3 * i + 2] = (qaws_scalar)zz;
	}
	for (k = 0; k < 6; k++)
	{
		unsigned int a = k / 2, u = (a + 1) % 3, w = (a + 2) % 3, j;
		double sg = k % 2 ? 0.5 : -0.5;
		for (i = 0; i < 20; i++)
			for (j = 0; j < 20; j++)
			{
				cx[3 * nc + a] = (qaws_scalar)sg;
				cx[3 * nc + u] = (qaws_scalar)((i + 0.5) / 20 - 0.5);
				cx[3 * nc + w] = (qaws_scalar)((j + 0.5) / 20 - 0.5);
				cn[3 * nc] = cn[3 * nc + 1] = cn[3 * nc + 2] = 0;
				cn[3 * nc + a] = (qaws_scalar)(sg > 0 ? 1 : -1);
				nc++;
			}
	}
	memset(&d, 0, sizeof d);
	d.resolution = 32;
	memset(&rs, 0, sizeof rs);
	memset(&rc, 0, sizeof rc);
	st = qaws_mesh_make_watertight_points(&s, &rs, x, nr, n, &d);
	sc = qaws_mesh_make_watertight_points(&c, &rc, cx, cn, nc, &d);
	if (s)
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(s);
		vs = wt_points_volume(s);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 q = qaws_mesh_3d_get_vertex(s, i);
			double r = sqrt((double)q.x * q.x + (double)q.y * q.y + (double)q.z * q.z);
			rmin = fmin(rmin, r);
			rmax = fmax(rmax, r);
		}
	}
	if (c)
		vc = wt_points_volume(c);
	for (i = 0; i < 3 * n; i++) nr[i] = -nr[i];
	sz = qaws_mesh_make_watertight_points(&z, NULL, x, nr, n, &d);
	s3 = qaws_mesh_make_watertight_points(&z, NULL, x, nr, 3, &d);
	sprintf(msg, "point clouds: a sphere's %u points: closed %d, volume %.4f (4.1888), radii %.4f..%.4f (cell %.3f); a cube's %u: closed %d, volume %.4f (1), %u sharp vertices; inward normals %d, 3 points %d",
		n, rs.closed, vs, rmin, rmax, (double)rs.cell, nc, rc.closed, vc, rc.sharp, (int)sz, (int)s3);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && s && rs.closed && fabs(vs - 4.18879) < 0.03 * 4.18879 && rmin > 1 - (double)rs.cell && rmax < 1 + (double)rs.cell &&
		sc == QAWS_STATUS_OK && c && rc.closed && fabs(vc - 1) < 0.03 && rc.sharp > 0 && sz == QAWS_STATUS_INVALID_TOPOLOGY && s3 == QAWS_STATUS_INVALID_ARGUMENT && !z, msg);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(c);
	free(x);
	free(nr);
	free(cx);
	free(cn);
}

int test_161_watertight_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 161: watertight surfaces\n");
	test_sphere();
	test_union();
	test_sheet();
	test_points();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
