/*
 * Test 164: geodesics on meshes (qaws_mesh_geodesic, the heat method)
 *
 *   - a unit icosphere, the source at a vertex: the distances against the
 *     great-circle arcs; two sources (the poles): the nearer one's;
 *   - a flat square grid with its border, the source at its centre or inside
 *     a triangle: the Euclidean distances, for each boundary condition;
 *   - two separate squares: the one without a source is unreached;
 *   - paths from the grid's corner and round the sphere: they end at the
 *     source, go down the distance all the way, and their lengths match it;
 *   - a second query on the same object gives the same distances;
 *   - refusals: no mesh, a source outside the triangles, a short buffer.
 */

#include "test_common.h"
#include "qaws_mesh_geodesic.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define GD_PI 3.14159265358979323846

/* the midpoint of edge a b, made once (open addressing over the edge's ends) */
static unsigned int gd_mid(unsigned long long* key, unsigned int* val, unsigned int cap, double* x, unsigned int* nv, unsigned int a, unsigned int b)
{
	unsigned long long k = a < b ? ((unsigned long long)a << 32 | b) + 1 : ((unsigned long long)b << 32 | a) + 1;
	unsigned int i = (unsigned int)((k * 0x9E3779B97F4A7C15ull) >> 40) & (cap - 1), j;
	double l;
	while (key[i] && key[i] != k)
		i = (i + 1) & (cap - 1);
	if (key[i])
		return val[i];
	key[i] = k;
	val[i] = *nv;
	for (j = 0; j < 3; j++)
		x[3 * *nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
	l = sqrt(x[3 * *nv] * x[3 * *nv] + x[3 * *nv + 1] * x[3 * *nv + 1] + x[3 * *nv + 2] * x[3 * *nv + 2]);
	for (j = 0; j < 3; j++)
		x[3 * *nv + j] /= l;
	return (*nv)++;
}

/* a unit icosphere: the icosahedron's faces split in four, levels times; vertex 0 at the north pole, 3 at the south */
static qaws_mesh_3d* gd_icosphere(unsigned int levels)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2, *x, ax[3], c, s;
	unsigned int nt = 20, nv = 12, cap, i, j, lv, *t, *t2, *val;
	unsigned long long* key;
	qaws_scalar* xs;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	for (lv = 0; lv < levels; lv++)
		nt *= 4;
	x = (double*)malloc(sizeof(double) * 3 * (nt / 2 + 2));
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	for (cap = 16; cap < 4 * nt; cap *= 2)
		;
	key = (unsigned long long*)malloc(sizeof(unsigned long long) * cap);
	val = (unsigned int*)malloc(sizeof(unsigned int) * cap);
	/* turned so vertex 0 (-1, p, 0) is at (0, 0, 1): about the axis (0, 0, 1) x v0 */
	{
		double l = sqrt(1 + p * p), n[3] = { -1 / l, p / l, 0 };
		ax[0] = n[1];
		ax[1] = -n[0];
		ax[2] = 0;
		c = n[2];
		s = sqrt(ax[0] * ax[0] + ax[1] * ax[1]);
		ax[0] /= s;
		ax[1] /= s;
	}
	for (i = 0; i < 12; i++)
	{
		double q[3], l = sqrt(v0[3 * i] * v0[3 * i] + v0[3 * i + 1] * v0[3 * i + 1] + v0[3 * i + 2] * v0[3 * i + 2]), k[3], dt;
		for (j = 0; j < 3; j++)
			q[j] = v0[3 * i + j] / l;
		/* Rodrigues: q cos + (ax x q) sin + ax (ax . q)(1 - cos) */
		k[0] = ax[1] * q[2] - ax[2] * q[1];
		k[1] = ax[2] * q[0] - ax[0] * q[2];
		k[2] = ax[0] * q[1] - ax[1] * q[0];
		dt = ax[0] * q[0] + ax[1] * q[1] + ax[2] * q[2];
		for (j = 0; j < 3; j++)
			x[3 * i + j] = q[j] * c + k[j] * s + ax[j] * dt * (1 - c);
	}
	for (i = 0; i < 60; i++)
		t[i] = (unsigned int)ic[i];
	for (nt = 20, lv = 0; lv < levels; lv++)
	{
		memset(key, 0, sizeof(unsigned long long) * cap);
		for (i = 0; i < nt; i++)
		{
			unsigned int a = t[3 * i], b = t[3 * i + 1], cc = t[3 * i + 2];
			unsigned int ab = gd_mid(key, val, cap, x, &nv, a, b), bc = gd_mid(key, val, cap, x, &nv, b, cc), ca = gd_mid(key, val, cap, x, &nv, cc, a);
			unsigned int const q[12] = { a, ab, ca, ab, b, bc, ca, bc, cc, ab, bc, ca };
			memcpy(t2 + 12 * i, q, sizeof q);
		}
		nt *= 4;
		memcpy(t, t2, sizeof(unsigned int) * 3 * nt);
	}
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	for (i = 0; i < 3 * nv; i++)
		xs[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = xs;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	free(x); free(t); free(t2); free(key); free(val); free(xs);
	return m;
}

/* squares of n x n cells side 1 at x offsets 0, 2, ... (count of them), each cell two triangles */
static qaws_mesh_3d* gd_grids(unsigned int n, unsigned int count)
{
	unsigned int nv = (n + 1) * (n + 1), i, j, q, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv * count);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n * count);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (q = 0; q < count; q++)
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				unsigned int v = q * nv + j * (n + 1) + i;
				x[3 * v] = (qaws_scalar)(2.0 * q + (double)i / n);
				x[3 * v + 1] = (qaws_scalar)((double)j / n);
				x[3 * v + 2] = 0;
			}
	for (q = 0; q < count; q++)
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int v = q * nv + j * (n + 1) + i;
				t[k++] = v;
				t[k++] = v + 1;
				t[k++] = v + n + 2;
				t[k++] = v;
				t[k++] = v + n + 2;
				t[k++] = v + n + 1;
			}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv * count;
	d.face_vertices = t;
	d.face_count = 2 * n * n * count;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

/* the location of vertex v in one of its triangles */
static qaws_mesh_location gd_at_vertex(qaws_mesh_3d const* m, unsigned int v)
{
	qaws_mesh_location l;
	unsigned int f, k, c[3];
	l.patch = QAWS_MESH_NONE;
	l.a = l.b = 0;
	for (f = 0; f < qaws_mesh_3d_get_triangle_count(m); f++)
	{
		qaws_mesh_3d_get_triangle(c, m, f);
		for (k = 0; k < 3; k++)
			if (c[k] == v)
			{
				l.patch = f;
				l.a = (qaws_scalar)(k == 1);
				l.b = (qaws_scalar)(k == 2);
				return l;
			}
	}
	return l;
}

static double gd_arc(qaws_vec3 a, qaws_vec3 b)
{
	double c = (double)a.x * (double)b.x + (double)a.y * (double)b.y + (double)a.z * (double)b.z;
	return acos(c > 1 ? 1 : c < -1 ? -1 : c);
}

static void test_sphere(void)
{
	qaws_mesh_3d* m = gd_icosphere(4);
	qaws_mesh_geodesic* g = NULL;
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), v;
	qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	qaws_scalar* d2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	qaws_mesh_location src[2];
	qaws_vec3 n = qaws_mesh_3d_get_vertex(m, 0), s = qaws_mesh_3d_get_vertex(m, 3);
	double emax = 0, esum = 0, e2max = 0;
	char msg[256];
	qaws_status st = qaws_mesh_geodesic_create(&g, m, NULL);
	TEST_ASSERT(st == QAWS_STATUS_OK && g, "icosphere: created");
	if (!g)
		goto done;
	src[0] = gd_at_vertex(m, 0);
	src[1] = gd_at_vertex(m, 3);
	st = qaws_mesh_geodesic_distance(d, g, src, 1);
	for (v = 0; v < nv; v++)
	{
		double e = fabs((double)d[v] - gd_arc(n, qaws_mesh_3d_get_vertex(m, v)));
		emax = e > emax ? e : emax;
		esum += e;
	}
	sprintf(msg, "icosphere (%u vertices), source at the pole: distance error max %.4f, mean %.4f, at the source %.2e", nv, emax, esum / nv,
		(double)d[0]);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && emax < 0.03 * GD_PI && esum / nv < 0.01 * GD_PI && d[0] == 0, msg);
	qaws_mesh_geodesic_distance(d2, g, src, 1);
	TEST_ASSERT(memcmp(d, d2, sizeof(qaws_scalar) * nv) == 0, "icosphere: a second query gives the same distances");
	st = qaws_mesh_geodesic_distance(d2, g, src, 2);
	for (v = 0; v < nv; v++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, v);
		double a = gd_arc(n, p), b = gd_arc(s, p), e = fabs((double)d2[v] - (a < b ? a : b));
		e2max = e > e2max ? e : e2max;
	}
	sprintf(msg, "icosphere, sources at both poles: the nearer one's distance, error max %.4f", e2max);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && e2max < 0.03 * GD_PI, msg);
	{
		/* the path from the south pole's neighbourhood back to the north pole, along a meridian */
		qaws_mesh_location from = gd_at_vertex(m, 3);
		qaws_scalar pts[3 * 512];
		unsigned int count = 0, i, down = 1;
		double len = 0;
		from.a = (qaws_scalar)0.3;
		from.b = (qaws_scalar)0.3;
		st = qaws_mesh_geodesic_path(pts, NULL, &count, g, d, src, 1, from, 512);
		for (i = 1; i < count; i++)
			len += sqrt((pts[3 * i] - pts[3 * i - 3]) * (pts[3 * i] - pts[3 * i - 3]) + (pts[3 * i + 1] - pts[3 * i - 2]) * (pts[3 * i + 1] - pts[3 * i - 2])
				+ (pts[3 * i + 2] - pts[3 * i - 1]) * (pts[3 * i + 2] - pts[3 * i - 1]));
		for (i = 1; i < count && down; i++)
			down = pts[3 * i + 2] >= pts[3 * i - 1] - 1e-9;   /* (toward the north pole: z grows) */
		{
			qaws_vec3 a;
			a.x = pts[0];
			a.y = pts[1];
			a.z = pts[2];
			{
				double r = sqrt((double)a.x * a.x + (double)a.y * a.y + (double)a.z * a.z), ex;
				a.x = (qaws_scalar)(a.x / r);
				a.y = (qaws_scalar)(a.y / r);
				a.z = (qaws_scalar)(a.z / r);
				ex = gd_arc(n, a);
				sprintf(msg, "icosphere path: %u points, length %.4f against the arc %.4f, ends at the pole %d, z rising all along %d", count, len, ex,
					count && fabs(pts[3 * count - 1] - n.z) < 1e-12, down);
				printf("  %s\n", msg);
				TEST_ASSERT(st == QAWS_STATUS_OK && count > 2 && fabs(pts[3 * count - 1] - n.z) < 1e-12 && down && fabs(len - ex) < 0.02 * ex, msg);
			}
		}
	}
done:
	qaws_mesh_geodesic_destroy(g);
	qaws_mesh_3d_destroy(m);
	free(d);
	free(d2);
}

static void test_plane(void)
{
	static char const* const bname[3] = { "mean", "Neumann", "Dirichlet" };
	unsigned int n = 40, nv = (n + 1) * (n + 1), v, b;
	qaws_mesh_3d* m = gd_grids(n, 1);
	qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	char msg[256];
	for (b = 0; b < 3; b++)
	{
		qaws_mesh_geodesic* g = NULL;
		qaws_mesh_geodesic_desc desc;
		qaws_mesh_location src = gd_at_vertex(m, (n / 2) * (n + 1) + n / 2);
		double emax = 0, esum = 0;
		qaws_status st;
		memset(&desc, 0, sizeof desc);
		desc.boundary = (qaws_geodesic_boundary)b;
		st = qaws_mesh_geodesic_create(&g, m, &desc);
		if (st == QAWS_STATUS_OK)
			st = qaws_mesh_geodesic_distance(d, g, &src, 1);
		for (v = 0; v < nv && st == QAWS_STATUS_OK; v++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, v);
			double e = fabs((double)d[v] - sqrt((p.x - 0.5) * (p.x - 0.5) + (p.y - 0.5) * (p.y - 0.5)));
			emax = e > emax ? e : emax;
			esum += e;
		}
		sprintf(msg, "square grid %u x %u, %s border, source at the centre: error max %.4f, mean %.4f", n, n, bname[b], emax, esum / nv);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && emax < (b == 0 ? 0.03 : 0.07) && esum / nv < (b == 0 ? 0.01 : 0.03), msg);
		if (b == 0)
		{
			/* a source inside a triangle; then the path from the corner (0, 0) */
			qaws_mesh_location in;
			qaws_scalar pts[3 * 512];
			unsigned int count = 0, i, down = 1;
			double len = 0, sx, sy;
			unsigned int c[3];
			in.patch = 2 * (13 * n + 7);
			in.a = (qaws_scalar)0.25;
			in.b = (qaws_scalar)0.4;
			qaws_mesh_3d_get_triangle(c, m, in.patch);
			{
				qaws_vec3 p0 = qaws_mesh_3d_get_vertex(m, c[0]), p1 = qaws_mesh_3d_get_vertex(m, c[1]), p2 = qaws_mesh_3d_get_vertex(m, c[2]);
				sx = 0.35 * p0.x + 0.25 * p1.x + 0.4 * p2.x;
				sy = 0.35 * p0.y + 0.25 * p1.y + 0.4 * p2.y;
			}
			st = qaws_mesh_geodesic_distance(d, g, &in, 1);
			emax = 0;
			for (v = 0; v < nv && st == QAWS_STATUS_OK; v++)
			{
				qaws_vec3 p = qaws_mesh_3d_get_vertex(m, v);
				double e = fabs((double)d[v] - sqrt((p.x - sx) * (p.x - sx) + (p.y - sy) * (p.y - sy)));
				emax = e > emax ? e : emax;
			}
			sprintf(msg, "square grid, source inside a triangle at (%.4f, %.4f): error max %.4f", sx, sy, emax);
			printf("  %s\n", msg);
			TEST_ASSERT(st == QAWS_STATUS_OK && emax < 0.05, msg);
			st = qaws_mesh_geodesic_path(pts, NULL, &count, g, d, &in, 1, gd_at_vertex(m, 0), 512);
			for (i = 1; i < count; i++)
				len += sqrt((pts[3 * i] - pts[3 * i - 3]) * (pts[3 * i] - pts[3 * i - 3]) + (pts[3 * i + 1] - pts[3 * i - 2]) * (pts[3 * i + 1] - pts[3 * i - 2]));
			for (i = 1; i < count && down; i++)
			{
				/* each point nearer the source than the one before */
				double a = (pts[3 * i] - sx) * (pts[3 * i] - sx) + (pts[3 * i + 1] - sy) * (pts[3 * i + 1] - sy);
				double p = (pts[3 * i - 3] - sx) * (pts[3 * i - 3] - sx) + (pts[3 * i - 2] - sy) * (pts[3 * i - 2] - sy);
				down = a <= p + 1e-12;
			}
			sprintf(msg, "square grid path from the corner: %u points, length %.4f against %.4f, ends at the source %d, nearer all along %d", count, len,
				sqrt(sx * sx + sy * sy), count && fabs(pts[3 * count - 3] - sx) < 1e-12 && fabs(pts[3 * count - 2] - sy) < 1e-12, down);
			printf("  %s\n", msg);
			TEST_ASSERT(st == QAWS_STATUS_OK && count > 2 && fabs(pts[3 * count - 3] - sx) < 1e-12 && fabs(pts[3 * count - 2] - sy) < 1e-12 && down
				&& fabs(len - sqrt(sx * sx + sy * sy)) < 0.005 * sqrt(sx * sx + sy * sy), msg);
			{
				unsigned int need = 0;
				st = qaws_mesh_geodesic_path(pts, NULL, &need, g, d, &in, 1, gd_at_vertex(m, 0), 2);
				TEST_ASSERT(st == QAWS_STATUS_BUFFER_TOO_SMALL && need == count, "a short buffer: BUFFER_TOO_SMALL with the count");
			}
		}
		qaws_mesh_geodesic_destroy(g);
	}
	qaws_mesh_3d_destroy(m);
	free(d);
}

/* a grid sheared along x (most of its triangles obtuse, negative cotangent weights): through the intrinsic Delaunay
   triangulation the distances stay those of the plane */
static void test_sheared(void)
{
	unsigned int n = 40, nv = (n + 1) * (n + 1), v;
	qaws_mesh_3d* m = gd_grids(n, 1);
	qaws_mesh_geodesic* g = NULL;
	qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_location src;
	double emax = 0, sx, sy;
	char msg[200];
	qaws_status st;
	qaws_mesh_3d_get_positions(x, m);
	for (v = 0; v < nv; v++)
		x[3 * v] = (qaws_scalar)(x[3 * v] - 2.5 * x[3 * v + 1]);
	qaws_mesh_3d_set_positions(m, x);
	src = gd_at_vertex(m, (n / 2) * (n + 1) + n / 2);
	sx = x[3 * ((n / 2) * (n + 1) + n / 2)];
	sy = x[3 * ((n / 2) * (n + 1) + n / 2) + 1];
	st = qaws_mesh_geodesic_create(&g, m, NULL);
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_geodesic_distance(d, g, &src, 1);
	for (v = 0; v < nv && st == QAWS_STATUS_OK; v++)
	{
		double e = fabs((double)d[v] - sqrt((x[3 * v] - sx) * (x[3 * v] - sx) + (x[3 * v + 1] - sy) * (x[3 * v + 1] - sy)));
		emax = e > emax ? e : emax;
	}
	sprintf(msg, "a grid sheared 2.5 to 1 (obtuse triangles): error max %.4f over a width of 3.5", emax);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && emax < 0.06, msg);
	qaws_mesh_geodesic_destroy(g);
	qaws_mesh_3d_destroy(m);
	free(d);
	free(x);
}

static void test_parts(void)
{
	unsigned int n = 8, nv = (n + 1) * (n + 1), v, far = 1, near = 1;
	qaws_mesh_3d* m = gd_grids(n, 2);
	qaws_mesh_geodesic* g = NULL;
	qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * nv);
	qaws_mesh_location src = gd_at_vertex(m, 0);
	qaws_status st = qaws_mesh_geodesic_create(&g, m, NULL);
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_geodesic_distance(d, g, &src, 1);
	for (v = 0; v < nv; v++)
	{
		near &= d[v] < 2;
		far &= d[nv + v] == (qaws_scalar)HUGE_VAL;
	}
	TEST_ASSERT(st == QAWS_STATUS_OK && near && far, "two squares, the source on the first: the second unreached (HUGE_VAL)");
	qaws_mesh_geodesic_destroy(g);
	qaws_mesh_3d_destroy(m);
	free(d);
}

static void test_refusals(void)
{
	qaws_mesh_3d* m = gd_grids(2, 1);
	qaws_mesh_geodesic* g = NULL;
	qaws_mesh_location bad;
	qaws_scalar d[9];
	unsigned int count;
	bad.patch = 99;
	bad.a = bad.b = 0;
	TEST_ASSERT(qaws_mesh_geodesic_create(&g, NULL, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !g, "no mesh refused");
	qaws_mesh_geodesic_create(&g, m, NULL);
	TEST_ASSERT(qaws_mesh_geodesic_distance(d, g, &bad, 1) == QAWS_STATUS_INVALID_ARGUMENT, "a source outside the triangles refused");
	TEST_ASSERT(qaws_mesh_geodesic_distance(d, g, &bad, 0) == QAWS_STATUS_INVALID_ARGUMENT, "no source refused");
	TEST_ASSERT(qaws_mesh_geodesic_path(NULL, NULL, &count, g, d, &bad, 1, bad, 0) == QAWS_STATUS_INVALID_ARGUMENT, "a path from outside refused");
	qaws_mesh_geodesic_destroy(g);
	qaws_mesh_3d_destroy(m);
}

/* the unit cube's surface, n x n quads a side (each cut in two), vertices shared by their integer coordinates */
static qaws_mesh_3d* gd_cube(unsigned int n)
{
	unsigned int N = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * N * N * N), nv = 0, nf = 0, a, s, u, v, k;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * N * N * N);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (k = 0; k < N * N * N; k++)
	{
		unsigned int q[3];
		q[0] = k % N;
		q[1] = (k / N) % N;
		q[2] = k / (N * N);
		id[k] = 0xFFFFFFFFu;
		if (q[0] % n && q[1] % n && q[2] % n)
			continue;
		id[k] = nv;
		for (a = 0; a < 3; a++) x[3 * nv + a] = (qaws_scalar)((double)q[a] / n);
		nv++;
	}
	for (a = 0; a < 3; a++)
		for (s = 0; s < 2; s++)
			for (u = 0; u < n; u++)
				for (v = 0; v < n; v++)
				{
					unsigned int q4[4], c, b = (a + 1) % 3, e = (a + 2) % 3;
					for (c = 0; c < 4; c++)
					{
						unsigned int q[3];
						q[a] = s * n;
						q[b] = u + (c == 1 || c == 2);
						q[e] = v + (c >= 2);
						q4[c] = id[q[0] + N * (q[1] + N * q[2])];
					}
					if (!s)
					{
						c = q4[1];
						q4[1] = q4[3];
						q4[3] = c;
					}
					f[3 * nf] = q4[0]; f[3 * nf + 1] = q4[1]; f[3 * nf + 2] = q4[2]; nf++;
					f[3 * nf] = q4[0]; f[3 * nf + 1] = q4[2]; f[3 * nf + 2] = q4[3]; nf++;
				}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(id);
	free(x);
	free(f);
	return m;
}

/* the geodesic distance on the unit cube's surface from the corner (0, 0, 0): straight on the three faces round it, else
   the shorter of the two unfoldings across the faces between */
static double gd_cube_distance(double const* p)
{
	unsigned int a;
	for (a = 0; a < 3; a++)
		if (p[a] < 1e-12)
			return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
	for (a = 0; a < 3; a++)
		if (p[a] > 1 - 1e-12)
		{
			double y = p[(a + 1) % 3], z = p[(a + 2) % 3];
			return fmin(sqrt((1 + y) * (1 + y) + z * z), sqrt((1 + z) * (1 + z) + y * y));
		}
	return -1;
}

static void test_exact(void)
{
	qaws_mesh_3d *g = gd_grids(12, 1), *cube = gd_cube(3), *two = gd_grids(4, 2);
	qaws_mesh_location src;
	qaws_scalar *dg, *dc, *d2;
	unsigned int i, ng = qaws_mesh_3d_get_vertex_count(g), nc = qaws_mesh_3d_get_vertex_count(cube), n2 = qaws_mesh_3d_get_vertex_count(two), far = 0;
	double eg = 0, ec = 0, corner = 0, sp[3], plen = 0, perr = 1;
	qaws_status s1, s2, s3;
	char msg[400];
	dg = (qaws_scalar*)malloc(sizeof(qaws_scalar) * ng);
	dc = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc);
	d2 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n2);
	/* a point inside a triangle of the flat grid: the distances Euclidean */
	src.patch = 77;
	src.a = (qaws_scalar)0.3;
	src.b = (qaws_scalar)0.2;
	{
		unsigned int c[3];
		qaws_vec3 p0, p1, p2;
		qaws_mesh_3d_get_triangle(c, g, src.patch);
		p0 = qaws_mesh_3d_get_vertex(g, c[0]);
		p1 = qaws_mesh_3d_get_vertex(g, c[1]);
		p2 = qaws_mesh_3d_get_vertex(g, c[2]);
		sp[0] = 0.5 * (double)p0.x + 0.3 * (double)p1.x + 0.2 * (double)p2.x;
		sp[1] = 0.5 * (double)p0.y + 0.3 * (double)p1.y + 0.2 * (double)p2.y;
		sp[2] = 0;
	}
	s1 = qaws_mesh_3d_geodesic_exact(dg, g, &src, 1, (qaws_scalar)HUGE_VAL);
	for (i = 0; i < ng; i++)
	{
		qaws_vec3 q = qaws_mesh_3d_get_vertex(g, i);
		eg = fmax(eg, fabs((double)dg[i] - sqrt(((double)q.x - sp[0]) * ((double)q.x - sp[0]) + ((double)q.y - sp[1]) * ((double)q.y - sp[1]))));
	}
	/* the cube from its corner (0, 0, 0) */
	for (i = 0; i < nc; i++)
	{
		qaws_vec3 q = qaws_mesh_3d_get_vertex(cube, i);
		if (q.x == 0 && q.y == 0 && q.z == 0)
			break;
	}
	src = gd_at_vertex(cube, i);
	s2 = qaws_mesh_3d_geodesic_exact(dc, cube, &src, 1, (qaws_scalar)HUGE_VAL);
	for (i = 0; i < nc; i++)
	{
		qaws_vec3 q = qaws_mesh_3d_get_vertex(cube, i);
		double pp[3] = { q.x, q.y, q.z }, e = gd_cube_distance(pp);
		ec = fmax(ec, fabs((double)dc[i] - e));
		if (q.x == 1 && q.y == 1 && q.z == 1)
			corner = (double)dc[i];
	}
	/* the path down the exact distances from (1, 2/3, 1/3) on the far face: as long as the unfolding's segment */
	{
		qaws_mesh_geodesic* geo = NULL;
		qaws_scalar pts[3 * 256];
		unsigned int count = 0, k;
		double want = 0;
		qaws_mesh_location from;
		for (i = 0; i < nc; i++)
		{
			qaws_vec3 q = qaws_mesh_3d_get_vertex(cube, i);
			if (q.x == 1 && fabs((double)q.y - 2.0 / 3) < 1e-6 && fabs((double)q.z - 1.0 / 3) < 1e-6)
			{
				double pp[3] = { q.x, q.y, q.z };
				want = gd_cube_distance(pp);
				break;
			}
		}
		from = gd_at_vertex(cube, i);
		qaws_mesh_geodesic_create(&geo, cube, NULL);
		if (geo && qaws_mesh_geodesic_path(pts, NULL, &count, geo, dc, &src, 1, from, 256) == QAWS_STATUS_OK)
			for (k = 1; k < count; k++)
				plen += sqrt(((double)pts[3 * k] - (double)pts[3 * k - 3]) * ((double)pts[3 * k] - (double)pts[3 * k - 3]) +
					((double)pts[3 * k + 1] - (double)pts[3 * k - 2]) * ((double)pts[3 * k + 1] - (double)pts[3 * k - 2]) +
					((double)pts[3 * k + 2] - (double)pts[3 * k - 1]) * ((double)pts[3 * k + 2] - (double)pts[3 * k - 1]));
		perr = fabs(plen - want);
		qaws_mesh_geodesic_destroy(geo);
	}
	/* two grids apart: the second not reached */
	src = gd_at_vertex(two, 0);
	s3 = qaws_mesh_3d_geodesic_exact(d2, two, &src, 1, (qaws_scalar)HUGE_VAL);
	for (i = 0; i < n2; i++) far += (double)d2[i] == HUGE_VAL;
	sprintf(msg, "exact geodesics: a flat grid from inside a triangle, off Euclid by %.1e; the cube from a corner, off its unfoldings by %.1e, the far corner at %.12f (sqrt 5 = 2.236067977500), a path down them %.12f long, off the unfolding by %.1e; two grids apart, %u of %u unreached",
		eg, ec, corner, plen, perr, far, n2);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && s2 == QAWS_STATUS_OK && s3 == QAWS_STATUS_OK && eg < 1e-12 && ec < 1e-12 && fabs(corner - sqrt(5.0)) < 1e-12 && perr < 1e-9 &&
		far == n2 / 2, msg);
	free(dg);
	free(dc);
	free(d2);
	qaws_mesh_3d_destroy(g);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(two);
}

/* saddles everywhere: an icosphere's vertices pushed out and in by turns; the exact distances above the chords, within an
   edge's length of their neighbours' and below the edge graph's */
static void test_exact_saddles(void)
{
	qaws_mesh_3d *s = gd_icosphere(2), *m = NULL;
	unsigned int nv = qaws_mesh_3d_get_vertex_count(s), nt = qaws_mesh_3d_get_triangle_count(s), i, j, k, bad = 0, done_count = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	qaws_scalar* d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	double* g = (double*)malloc(sizeof(double) * nv);
	unsigned char* done = (unsigned char*)calloc(nv, 1);
	qaws_mesh_3d_desc desc;
	qaws_mesh_location src;
	qaws_status st;
	double below = 0;
	char msg[256];
	qaws_mesh_3d_get_positions(x, s);
	qaws_mesh_3d_get_triangles(t, s);
	for (i = 0; i < nv; i++)
		for (k = 0; k < 3; k++) x[3 * i + k] *= (qaws_scalar)(i % 2 ? 1.2 : 0.85);
	memset(&desc, 0, sizeof desc);
	desc.positions = x;
	desc.vertex_count = nv;
	desc.face_vertices = t;
	desc.face_count = nt;
	qaws_mesh_3d_create(&m, &desc);
	src = gd_at_vertex(m, 0);
	st = qaws_mesh_3d_geodesic_exact(d, m, &src, 1, (qaws_scalar)HUGE_VAL);
	/* the edge graph's distances (Dijkstra, quadratic: a small mesh) */
	for (i = 0; i < nv; i++) g[i] = HUGE_VAL;
	g[0] = 0;
	while (done_count < nv)
	{
		unsigned int u = nv;
		for (i = 0; i < nv; i++)
			if (!done[i] && (u == nv || g[i] < g[u]))
				u = i;
		done[u] = 1;
		done_count++;
		for (i = 0; i < 3 * nt; i++)
			if (t[i] == u)
				for (j = 1; j <= 2; j++)
				{
					unsigned int w = t[3 * (i / 3) + (i % 3 + j) % 3];
					double l = sqrt(((double)x[3 * u] - (double)x[3 * w]) * ((double)x[3 * u] - (double)x[3 * w]) + ((double)x[3 * u + 1] - (double)x[3 * w + 1]) * ((double)x[3 * u + 1] - (double)x[3 * w + 1]) + ((double)x[3 * u + 2] - (double)x[3 * w + 2]) * ((double)x[3 * u + 2] - (double)x[3 * w + 2]));
					if (g[u] + l < g[w])
						g[w] = g[u] + l;
					/* (neighbours within an edge of each other, once Dijkstra has passed) */
					if (fabs((double)d[u] - (double)d[w]) > l + 1e-12)
						bad++;
				}
	}
	for (i = 0; i < nv; i++)
	{
		double c = sqrt(((double)x[3 * i] - (double)x[0]) * ((double)x[3 * i] - (double)x[0]) + ((double)x[3 * i + 1] - (double)x[1]) * ((double)x[3 * i + 1] - (double)x[1]) + ((double)x[3 * i + 2] - (double)x[2]) * ((double)x[3 * i + 2] - (double)x[2]));
		bad += (double)d[i] < c - 1e-12 || (double)d[i] > g[i] + 1e-12;
		below += g[i] - (double)d[i];
	}
	sprintf(msg, "exact geodesics on saddles (%u vertices): %u violations, %.3f below the edge graph on average", nv, bad, below / nv);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && bad == 0 && below > 0, msg);
	{
		/* up to a radius: the same distances within it, HUGE_VAL beyond */
		qaws_scalar* dr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		unsigned int in = 0, out = 0, wrong = 0;
		qaws_status sr = qaws_mesh_3d_geodesic_exact(dr, m, &src, 1, (qaws_scalar)1.2), sn = qaws_mesh_3d_geodesic_exact(dr, m, &src, 1, (qaws_scalar)-1);
		for (i = 0; i < nv; i++)
		{
			if ((double)d[i] <= 1.2)
			{
				in++;
				wrong += dr[i] != d[i];
			}
			else
			{
				out++;
				wrong += dr[i] != (qaws_scalar)HUGE_VAL;
			}
		}
		sprintf(msg, "exact geodesics up to 1.2: %u vertices within, %u beyond, %u differing from the full run; a negative radius refused (%d)", in, out, wrong, (int)sn);
		printf("  %s\n", msg);
		TEST_ASSERT(sr == QAWS_STATUS_OK && sn == QAWS_STATUS_INVALID_ARGUMENT && in > 10 && out > 10 && wrong == 0, msg);
		free(dr);
	}
	free(x);
	free(t);
	free(d);
	free(g);
	free(done);
	qaws_mesh_3d_destroy(s);
	qaws_mesh_3d_destroy(m);
}

int test_164_mesh_geodesic_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 164: geodesics on meshes\n");
	test_sphere();
	test_plane();
	test_sheared();
	test_exact();
	test_exact_saddles();
	test_parts();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
