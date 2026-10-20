/*
 * Test 139: repair for farq (qaws_mesh_repair)
 *
 *   - a cube soup (each triangle its own 3 points) welded to the closed cube
 *     of 8 vertices; a cube of quads and an L-shaped prism (a non-convex
 *     cap) cut into triangles;
 *   - a sphere soup with its duplicate points off by 1e-7 welded;
 *   - a sphere with three holes of different sizes: closed, the fills' area
 *     near the area removed, their edges near the surrounding length, no
 *     sliver; the perimeter limit leaves the large ones;
 *   - a box in two halves meshed at different spacings, shifted apart by
 *     less than the gap tolerance: bridged to closed without fills;
 *   - a sphere with half its faces turned: consistent and outward;
 *   - a box with a cavity, and an island inside the cavity: the cavity
 *     turned inward, the island outward; or removed;
 *   - a cap (a corner on its long edge), needles and a zero-area face
 *     removed without opening holes; duplicate and opposite faces;
 *   - a Moebius strip reported non-orientable and left open; two cubes on
 *     one edge: closed, the non-manifold edge left, the directed edges
 *     cancelling through the vertex sources;
 *   - refusals.
 * Every closed result: no directed edge left uncancelled, and the store's
 * own validation.
 */

#include "test_common.h"
#include "qaws_mesh_repair.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_status mr_repair(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf,
	qaws_mesh_repair_desc const* d, qaws_mesh_3d** out, qaws_mesh_repair_report* rep)
{
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv ? nv : 1));
	unsigned int i;
	qaws_status st;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	st = qaws_mesh_repair(out, rep, p, f, fs, nv, nf, d);
	free(p);
	return st;
}

static int mr_ull_cmp(void const* a, void const* b)
{
	unsigned long long x = *(unsigned long long const*)a, y = *(unsigned long long const*)b;
	return x < y ? -1 : x > y;
}

/* the edges whose directed uses do not cancel, the indices mapped through the vertex sources (as farq's point keys merge them) */
static unsigned int mr_uncancelled(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, k, n = 0;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	unsigned long long* e = (unsigned long long*)malloc(sizeof(unsigned long long) * 3 * (nt + 1));
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
		for (k = 0; k < 3; k++)
		{
			unsigned long long a = qaws_mesh_3d_get_vertex_source(m, t[3 * i + k]), b = qaws_mesh_3d_get_vertex_source(m, t[3 * i + (k + 1) % 3]);
			e[3 * i + k] = a < b ? (a << 33 | b << 1) : (b << 33 | a << 1 | 1);
		}
	qsort(e, 3 * nt, sizeof *e, mr_ull_cmp);
	for (i = 0; i < 3 * nt; i = k)
	{
		int net = 0;
		for (k = i; k < 3 * nt && e[k] >> 1 == e[i] >> 1; k++)
			net += e[k] & 1 ? -1 : 1;
		n += net != 0;
	}
	free(t);
	free(e);
	return n;
}

static double mr_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	double v = 0;
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		qaws_vec3 a = qaws_mesh_3d_get_vertex(m, t[3 * i]), b = qaws_mesh_3d_get_vertex(m, t[3 * i + 1]), c = qaws_mesh_3d_get_vertex(m, t[3 * i + 2]);
		v += ((double)a.x * ((double)b.y * c.z - (double)b.z * c.y) - (double)a.y * ((double)b.x * c.z - (double)b.z * c.x)
			+ (double)a.z * ((double)b.x * c.y - (double)b.y * c.x)) / 6;
	}
	free(t);
	return v;
}

/* closed, every directed edge cancelling, the store valid */
static int mr_closed(qaws_mesh_3d const* m, qaws_mesh_repair_report const* r)
{
	return m && r->closed && !r->open_edges && mr_uncancelled(m) == 0 && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK;
}

typedef struct mr_stats
{
	unsigned int fill_faces;
	double fill_area, fill_edge, other_edge, fill_min_angle, min_height;
} mr_stats;

/* the fills (source_face < 0): area, mean edge, smallest angle; the other faces' mean edge; the lowest face of all */
static void mr_measure(qaws_mesh_3d const* m, mr_stats* s)
{
	unsigned int nt = qaws_mesh_3d_get_face_count(m), i, k, a = 0, ne[2] = { 0, 0 };
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	qaws_scalar* src = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nt + 1));
	double le[2] = { 0, 0 };
	memset(s, 0, sizeof *s);
	s->fill_min_angle = 180;
	s->min_height = HUGE_VAL;
	qaws_mesh_3d_get_triangles(t, m);
	qaws_mesh_3d_find_attribute(&a, m, "source_face");
	qaws_mesh_3d_get_attribute(src, m, a);
	for (i = 0; i < nt; i++)
	{
		double p[3][3], l[3], lmax = 0, cr[3], area;
		int fill = (double)src[i] < 0;
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 v = qaws_mesh_3d_get_vertex(m, t[3 * i + k]);
			p[k][0] = (double)v.x;
			p[k][1] = (double)v.y;
			p[k][2] = (double)v.z;
		}
		for (k = 0; k < 3; k++)
		{
			double const *u = p[k], *w = p[(k + 1) % 3];
			l[k] = sqrt((u[0] - w[0]) * (u[0] - w[0]) + (u[1] - w[1]) * (u[1] - w[1]) + (u[2] - w[2]) * (u[2] - w[2]));
			lmax = fmax(lmax, l[k]);
			le[fill] += l[k];
			ne[fill]++;
		}
		cr[0] = (p[1][1] - p[0][1]) * (p[2][2] - p[0][2]) - (p[1][2] - p[0][2]) * (p[2][1] - p[0][1]);
		cr[1] = (p[1][2] - p[0][2]) * (p[2][0] - p[0][0]) - (p[1][0] - p[0][0]) * (p[2][2] - p[0][2]);
		cr[2] = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]);
		area = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
		s->min_height = fmin(s->min_height, lmax > 0 ? 2 * area / lmax : 0);
		if (!fill)
			continue;
		s->fill_faces++;
		s->fill_area += area;
		for (k = 0; k < 3; k++)
		{
			/* the angle opposite edge k + 1 (law of cosines) */
			double x = l[k], y = l[(k + 2) % 3], z = l[(k + 1) % 3];
			double c = (x * x + y * y - z * z) / (2 * x * y);
			s->fill_min_angle = fmin(s->fill_min_angle, acos(c < -1 ? -1 : c > 1 ? 1 : c) * 180 / M_PI);
		}
	}
	s->fill_edge = ne[1] ? le[1] / ne[1] : 0;
	s->other_edge = ne[0] ? le[0] / ne[0] : 0;
	free(t);
	free(src);
}

/* ---- shapes ------------------------------------------------------------- */

static double const mr_cube_x[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
static unsigned int const mr_cube_t[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6,
	1, 6, 5 };

/* each triangle its own three points, moved by up to eps (a fixed sequence) */
static void mr_soup(double const* x, unsigned int const* t, unsigned int nt, double eps, double* sx, unsigned int* st)
{
	unsigned int i, k, s = 12345;
	for (i = 0; i < 3 * nt; i++)
	{
		st[i] = i;
		for (k = 0; k < 3; k++)
		{
			s = s * 1664525u + 1013904223u;
			sx[3 * i + k] = x[3 * t[i] + k] + eps * (((s >> 8) & 0xFFFF) / 32768.0 - 1);
		}
	}
}

/* boxes [lo, hi] appended (12 triangles each, outward) */
static void mr_box(double const* lo, double const* hi, double* x, unsigned int* nv, unsigned int* t, unsigned int* nt)
{
	unsigned int i, k;
	for (i = 0; i < 8; i++)
		for (k = 0; k < 3; k++)
			x[3 * (*nv + i) + k] = mr_cube_x[3 * i + k] ? hi[k] : lo[k];
	for (i = 0; i < 36; i++)
		t[3 * *nt + i] = *nv + mr_cube_t[i];
	*nv += 8;
	*nt += 12;
}

/* the sphere of radius 1: an icosahedron split level times, on the sphere; returns the triangles (x: 3 * (10 * 4^level + 2),
   t: 3 * 20 * 4^level) */
static unsigned int mr_icosphere(unsigned int level, double* x, unsigned int* t, unsigned int* out_nv)
{
	static int const f0[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2;
	double const v0[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	unsigned int nv = 12, nt = 20, l, i, k, cap = 1u << 20;
	unsigned long long* key = (unsigned long long*)malloc(sizeof(unsigned long long) * cap);
	unsigned int* val = (unsigned int*)malloc(sizeof(unsigned int) * cap);
	unsigned int* nt3 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 20 * ((size_t)1 << (2 * level)));
	for (i = 0; i < 36; i++)
		x[i] = v0[i] / sqrt(1 + g * g);
	for (i = 0; i < 60; i++)
		t[i] = (unsigned int)f0[i];
	for (l = 0; l < level; l++)
	{
		unsigned int n = 0;
		memset(key, 0xFF, sizeof(unsigned long long) * cap);
		for (i = 0; i < nt; i++)
		{
			unsigned int mid[3];
			for (k = 0; k < 3; k++)
			{
				unsigned int a = t[3 * i + k], b = t[3 * i + (k + 1) % 3];
				unsigned long long kk = a < b ? (unsigned long long)a << 32 | b : (unsigned long long)b << 32 | a;
				unsigned int h = (unsigned int)((kk * 0x9E3779B97F4A7C15ull) >> 44) & (cap - 1);
				while (key[h] != kk && key[h] != ~0ull)
					h = (h + 1) & (cap - 1);
				if (key[h] != kk)
				{
					double* p = x + 3 * nv;
					double r;
					key[h] = kk;
					val[h] = nv;
					p[0] = x[3 * a] + x[3 * b];
					p[1] = x[3 * a + 1] + x[3 * b + 1];
					p[2] = x[3 * a + 2] + x[3 * b + 2];
					r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
					p[0] /= r;
					p[1] /= r;
					p[2] /= r;
					nv++;
				}
				mid[k] = val[h];
			}
			{
				unsigned int const q[12] = { t[3 * i], mid[0], mid[2], mid[0], t[3 * i + 1], mid[1], mid[2], mid[1], t[3 * i + 2], mid[0], mid[1], mid[2] };
				memcpy(nt3 + 3 * n, q, sizeof q);
				n += 4;
			}
		}
		memcpy(t, nt3, sizeof(unsigned int) * 3 * n);
		nt = n;
	}
	free(key);
	free(val);
	free(nt3);
	*out_nv = nv;
	return nt;
}

/* a grid on the parallelogram o + s U + t V, s, t in [0, 1], about n cells per unit, turned by U x V; appended as its own points */
static void mr_patch(double const* o, double const* U, double const* V, double n, double const* shift, double* x, unsigned int* nv,
	unsigned int* t, unsigned int* nt)
{
	unsigned int nu = (unsigned int)fmax(1, floor(sqrt(U[0] * U[0] + U[1] * U[1] + U[2] * U[2]) * n + 0.5));
	unsigned int nw = (unsigned int)fmax(1, floor(sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]) * n + 0.5));
	unsigned int i, j, k, b = *nv;
	for (i = 0; i <= nu; i++)
		for (j = 0; j <= nw; j++)
		{
			for (k = 0; k < 3; k++)
				x[3 * *nv + k] = o[k] + (double)i / nu * U[k] + (double)j / nw * V[k] + shift[k];
			(*nv)++;
		}
	for (i = 0; i < nu; i++)
		for (j = 0; j < nw; j++)
		{
			unsigned int p00 = b + i * (nw + 1) + j, p10 = p00 + nw + 1, p11 = p10 + 1, p01 = p00 + 1;
			unsigned int q[6] = { p00, p10, p11, p00, p11, p01 };
			memcpy(t + 3 * *nt, q, sizeof q);
			*nt += 2;
		}
}

/* ---- tests -------------------------------------------------------------- */

static void test_soups(void)
{
	double sx[108];
	unsigned int st[36], i;
	/* a cube of quads, and an L-shaped prism: its caps hexagons with a reflex corner */
	static double const qx[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
	static unsigned int const qf[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
	static unsigned int const qs[6] = { 4, 4, 4, 4, 4, 4 };
	static double const L[12] = { 0, 0, 2, 0, 2, 1, 1, 1, 1, 2, 0, 2 };
	double lx[36];
	unsigned int lf[12 + 24], ls[8], n = 0;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_status s;
	char msg[256];
	mr_soup(mr_cube_x, mr_cube_t, 12, 0, sx, st);
	s = mr_repair(sx, 36, st, NULL, 12, NULL, &m, &r);
	sprintf(msg, "a cube soup: %u vertices, %u faces, %u welded, volume %.6f", m ? qaws_mesh_3d_get_vertex_count(m) : 0,
		m ? qaws_mesh_3d_get_face_count(m) : 0, r.welded_vertices, m ? mr_volume(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && qaws_mesh_3d_get_vertex_count(m) == 8 && qaws_mesh_3d_get_face_count(m) == 12
		&& r.welded_vertices == 28 && fabs(mr_volume(m) - 1) < 1e-12 && r.flipped_faces == 0 && r.holes_filled == 0, msg);
	qaws_mesh_3d_destroy(m);
	s = mr_repair(qx, 8, qf, qs, 6, NULL, &m, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && qaws_mesh_3d_get_face_count(m) == 12 && fabs(mr_volume(m) - 1) < 1e-12,
		"a cube of quads: 12 triangles, closed");
	qaws_mesh_3d_destroy(m);
	for (i = 0; i < 12; i++)
	{
		lx[3 * i] = L[2 * (i % 6)];
		lx[3 * i + 1] = L[2 * (i % 6) + 1];
		lx[3 * i + 2] = i < 6 ? 0 : 1;
	}
	for (i = 0; i < 6; i++)
		lf[n++] = 5 - i;
	ls[0] = 6;
	for (i = 0; i < 6; i++)
		lf[n++] = 6 + i;
	ls[1] = 6;
	for (i = 0; i < 6; i++)
	{
		lf[n++] = i;
		lf[n++] = (i + 1) % 6;
		lf[n++] = 6 + (i + 1) % 6;
		lf[n++] = 6 + i;
		ls[2 + i] = 4;
	}
	s = mr_repair(lx, 12, lf, ls, 8, NULL, &m, &r);
	sprintf(msg, "an L-shaped prism: %u triangles, volume %.6f (3)", m ? qaws_mesh_3d_get_face_count(m) : 0, m ? mr_volume(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && qaws_mesh_3d_get_face_count(m) == 20 && fabs(mr_volume(m) - 3) < 1e-12, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_near_duplicates(void)
{
	double* x = (double*)malloc(sizeof(double) * 3 * 700);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 1280);
	double* sx = (double*)malloc(sizeof(double) * 9 * 1280);
	unsigned int* st = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 1280);
	unsigned int nv, nt = mr_icosphere(3, x, t, &nv);
	qaws_mesh_3d *m = NULL, *ref = NULL;
	qaws_mesh_repair_report r;
	qaws_status s;
	char msg[256];
	double v0;
	s = mr_repair(x, nv, t, NULL, nt, NULL, &ref, &r);
	v0 = ref ? mr_volume(ref) : 0;
	mr_soup(x, t, nt, 1e-7, sx, st);
	s = mr_repair(sx, 3 * nt, st, NULL, nt, NULL, &m, &r);
	sprintf(msg, "a sphere soup off by 1e-7: %u of %u vertices, %u faces, volume %.8f (%.8f)", m ? qaws_mesh_3d_get_vertex_count(m) : 0, nv,
		m ? qaws_mesh_3d_get_face_count(m) : 0, m ? mr_volume(m) : 0, v0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && qaws_mesh_3d_get_vertex_count(m) == nv && qaws_mesh_3d_get_face_count(m) == nt
		&& fabs(mr_volume(m) - v0) < 1e-5 && v0 > 4.0 && r.flipped_faces == 0, msg);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(ref);
	free(x); free(t); free(sx); free(st);
}

static void test_holes(void)
{
	static double const dir[3][3] = { { 0, 0, 1 }, { 1, 0, 0 }, { 0, -0.6, -0.8 } };
	static double const ang[3] = { 0.15, 0.35, 0.7 };
	double* x = (double*)malloc(sizeof(double) * 3 * 2600);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 5120);
	unsigned int* kept = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 5120);
	unsigned int nv, nt = mr_icosphere(4, x, t, &nv), nk = 0, i, k, h;
	double removed = 0;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_mesh_repair_desc d;
	mr_stats ms;
	qaws_status s;
	char msg[300];
	for (i = 0; i < nt; i++)
	{
		double c[3] = { 0, 0, 0 }, cr[3], *a = x + 3 * t[3 * i], *b = x + 3 * t[3 * i + 1], *e = x + 3 * t[3 * i + 2], l;
		int out = 0;
		for (k = 0; k < 3; k++)
			c[k] = (a[k] + b[k] + e[k]) / 3;
		l = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
		for (h = 0; h < 3; h++)
			out |= acos((c[0] * dir[h][0] + c[1] * dir[h][1] + c[2] * dir[h][2]) / l) < ang[h];
		if (!out)
		{
			memcpy(kept + 3 * nk++, t + 3 * i, sizeof(unsigned int) * 3);
			continue;
		}
		cr[0] = (b[1] - a[1]) * (e[2] - a[2]) - (b[2] - a[2]) * (e[1] - a[1]);
		cr[1] = (b[2] - a[2]) * (e[0] - a[0]) - (b[0] - a[0]) * (e[2] - a[2]);
		cr[2] = (b[0] - a[0]) * (e[1] - a[1]) - (b[1] - a[1]) * (e[0] - a[0]);
		removed += 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
	}
	s = mr_repair(x, nv, kept, NULL, nk, NULL, &m, &r);
	if (m)
		mr_measure(m, &ms);
	else
		memset(&ms, 0, sizeof ms);
	sprintf(msg, "a sphere with 3 holes: %u filled, fill area %.4f of %.4f removed, fill edges %.4f around %.4f, smallest angle %.1f, %u fill faces",
		r.holes_filled, ms.fill_area, removed, ms.fill_edge, ms.other_edge, ms.fill_min_angle, ms.fill_faces);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && r.holes_filled == 3 && r.holes_left == 0 && fabs(ms.fill_area / removed - 1) < 0.03
		&& fabs(ms.fill_edge / ms.other_edge - 1) < 0.25 && ms.fill_min_angle > 15 && mr_volume(m) > 4.0 && ms.fill_faces == r.fill_faces, msg);
	qaws_mesh_3d_destroy(m);
	/* fills flat, or as the minimal triangulation */
	memset(&d, 0, sizeof d);
	d.flags = QAWS_MESH_REPAIR_FLAT_FILLS;
	s = mr_repair(x, nv, kept, NULL, nk, &d, &m, &r);
	if (m)
		mr_measure(m, &ms);
	sprintf(msg, "flat fills: area %.4f, below the %.4f removed (a sphere's), closed", ms.fill_area, removed);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && ms.fill_area < 0.97 * removed, msg);
	qaws_mesh_3d_destroy(m);
	{
		/* coarse fills add no vertex */
		unsigned char* use = (unsigned char*)calloc(nv, 1);
		unsigned int nu = 0;
		for (i = 0; i < 3 * nk; i++)
			use[kept[i]] = 1;
		for (i = 0; i < nv; i++)
			nu += use[i];
		free(use);
		d.flags = QAWS_MESH_REPAIR_COARSE_FILLS;
		s = mr_repair(x, nv, kept, NULL, nk, &d, &m, &r);
		TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && qaws_mesh_3d_get_vertex_count(m) == nu && r.fill_faces > 0,
			"coarse fills: closed, no vertex added");
		qaws_mesh_3d_destroy(m);
	}
	/* the perimeter limit (the loops' zigzag about 1.0 to 1.35 times the circles'): the two smaller holes filled */
	d.flags = 0;
	d.hole_perimeter = (qaws_scalar)3.5;
	s = mr_repair(x, nv, kept, NULL, nk, &d, &m, &r);
	sprintf(msg, "perimeter up to %.3f: %u filled, %u left, %u boundary edges, closed %d", (double)d.hole_perimeter, r.holes_filled, r.holes_left,
		r.boundary_edges, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.holes_filled == 2 && r.holes_left == 1 && !r.closed && r.boundary_edges > 0
		&& r.boundary_edges == r.open_edges && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, msg);
	qaws_mesh_3d_destroy(m);
	free(x); free(t); free(kept);
}

static void test_gap(void)
{
	/* the box [0, 1]^3 in two halves: below z = 0.5 at 6 cells per unit, above at 5, shifted by delta */
	static double const zero[3] = { 0, 0, 0 }, delta[3] = { 2e-4, 1e-4, 2e-4 };
	double* x = (double*)malloc(sizeof(double) * 3 * 400);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 600);
	unsigned int nv = 0, nt = 0, h;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_mesh_repair_desc d;
	qaws_status s;
	char msg[256];
	for (h = 0; h < 2; h++)
	{
		double const* sh = h ? delta : zero;
		double n = h ? 5 : 6, z0 = h ? 0.5 : 0;
		double o[3], U[3] = { 0, 0, 0.5 }, V[3];
		/* the floor (outward -z) or the ceiling (+z) */
		double fo[3] = { 0, 0, h ? 1 : 0 }, fu[3] = { h ? 1 : 0, h ? 0 : 1, 0 }, fv[3] = { h ? 0 : 1, h ? 1 : 0, 0 };
		mr_patch(fo, fu, fv, n, sh, x, &nv, t, &nt);
		/* x = 0 (u = z, v = y), x = 1 (u = y, v = z), y = 0 (u = x, v = z), y = 1 (u = z, v = x) */
		o[0] = 0; o[1] = 0; o[2] = z0; V[0] = 0; V[1] = 1; V[2] = 0;
		mr_patch(o, U, V, n, sh, x, &nv, t, &nt);
		o[0] = 1; mr_patch(o, V, U, n, sh, x, &nv, t, &nt);
		o[0] = 0; V[0] = 1; V[1] = 0;
		mr_patch(o, V, U, n, sh, x, &nv, t, &nt);
		o[1] = 1; mr_patch(o, U, V, n, sh, x, &nv, t, &nt);
	}
	memset(&d, 0, sizeof d);
	d.gap_tolerance = (qaws_scalar)1e-3;
	s = mr_repair(x, nv, t, NULL, nt, &d, &m, &r);
	sprintf(msg, "a box in halves at 6 and 5 cells, 3e-4 apart: %u stitched, %u filled, %u faces, volume %.6f", r.stitched_vertices, r.holes_filled,
		m ? qaws_mesh_3d_get_face_count(m) : 0, m ? mr_volume(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && r.holes_filled == 0 && r.stitched_vertices >= 40 && fabs(mr_volume(m) - 1) < 1e-3
		&& r.flipped_faces == 0, msg);
	qaws_mesh_3d_destroy(m);
	/* without the gap policy the seam is a hole between two loops: filled as two holes, not stitched */
	d.gap_tolerance = -1;
	d.hole_perimeter = -1;
	s = mr_repair(x, nv, t, NULL, nt, &d, &m, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && !r.closed && r.stitched_vertices == 0 && r.boundary_edges >= 44, "gaps off: left open");
	qaws_mesh_3d_destroy(m);
	free(x); free(t);
}

static void test_orientation(void)
{
	double* x = (double*)malloc(sizeof(double) * 3 * 700);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 1280);
	unsigned int nv, nt = mr_icosphere(3, x, t, &nv), i, turned = 0, sd = 7;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_status s;
	char msg[256];
	for (i = 0; i < nt; i++)
	{
		sd = sd * 1664525u + 1013904223u;
		if ((sd >> 16) & 1)
		{
			unsigned int k = t[3 * i + 1];
			t[3 * i + 1] = t[3 * i + 2];
			t[3 * i + 2] = k;
			turned++;
		}
	}
	s = mr_repair(x, nv, t, NULL, nt, NULL, &m, &r);
	sprintf(msg, "a sphere with %u of %u faces turned: %u turned back, volume %.4f", turned, nt, r.flipped_faces, m ? mr_volume(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && r.flipped_faces == turned && mr_volume(m) > 4.0 && r.non_orientable == 0, msg);
	qaws_mesh_3d_destroy(m);
	/* all turned: turned back outward */
	for (i = 0; i < nt; i++)
	{
		unsigned int k = t[3 * i + 1];
		t[3 * i + 1] = t[3 * i + 2];
		t[3 * i + 2] = k;
	}
	s = mr_repair(x, nv, t, NULL, nt, NULL, &m, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && r.flipped_faces == nt - turned && mr_volume(m) > 4.0, "the other half turned");
	qaws_mesh_3d_destroy(m);
	free(x); free(t);
}

static void test_nested(void)
{
	static double const lo[3][3] = { { 0, 0, 0 }, { 0.3, 0.3, 0.3 }, { 0.4, 0.4, 0.4 } }, hi[3][3] = { { 1, 1, 1 }, { 0.7, 0.7, 0.7 }, { 0.6, 0.6, 0.6 } };
	double x[72];
	unsigned int t[108], nv = 0, nt = 0, b;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_mesh_repair_desc d;
	qaws_status s;
	char msg[256];
	for (b = 0; b < 3; b++)
		mr_box(lo[b], hi[b], x, &nv, t, &nt);
	/* the cavity given outward, the island inward */
	for (b = 72; b < 108; b += 3)
	{
		unsigned int k = t[b + 1];
		t[b + 1] = t[b + 2];
		t[b + 2] = k;
	}
	s = mr_repair(x, 16, t, NULL, 24, NULL, &m, &r);
	sprintf(msg, "a box with a cavity: volume %.6f (0.936), %u nested", m ? mr_volume(m) : 0, r.nested_shells);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && fabs(mr_volume(m) - 0.936) < TEST_BOUND(1e-9, 32) && r.nested_shells == 1, msg);
	qaws_mesh_3d_destroy(m);
	s = mr_repair(x, nv, t, NULL, nt, NULL, &m, &r);
	sprintf(msg, "and an island in the cavity: volume %.6f (0.944), %u nested", m ? mr_volume(m) : 0, r.nested_shells);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && fabs(mr_volume(m) - 0.944) < TEST_BOUND(1e-9, 32) && r.nested_shells == 2, msg);
	qaws_mesh_3d_destroy(m);
	memset(&d, 0, sizeof d);
	d.nested_shells = QAWS_NESTED_SHELLS_REMOVE;
	s = mr_repair(x, nv, t, NULL, nt, &d, &m, &r);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && fabs(mr_volume(m) - 1) < 1e-9 && qaws_mesh_3d_get_face_count(m) == 12
		&& r.nested_shells == 2, "nested shells removed");
	qaws_mesh_3d_destroy(m);
}

static void test_degenerate(void)
{
	/* the cube; a cap: face (0, 2, 1) split at the bottom's centre 8, the cap (0, 2, 8) across its diagonal; needles: the top's
	   diagonal (4, 6) split 1e-4 from 4 (vertex 9); a face of three points 1e-9 apart; a face with a vertex twice */
	double x[39];
	unsigned int t[3 * 20], nt = 0, i;
	static unsigned int const extra[] = { 0, 8, 1, 8, 2, 1, 0, 2, 8, 6, 9, 5, 9, 4, 5, 4, 9, 7, 9, 6, 7, 10, 11, 12, 5, 5, 6 };
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_mesh_repair_desc d;
	mr_stats ms;
	qaws_status s;
	char msg[300];
	double e = 1e-4 / sqrt(2.0);
	memcpy(x, mr_cube_x, sizeof mr_cube_x);
	x[24] = 0.5; x[25] = 0.5; x[26] = 0;
	x[27] = e; x[28] = e; x[29] = 1;
	x[30] = 2; x[31] = 2; x[32] = 2;
	x[33] = 2 + 1e-9; x[34] = 2; x[35] = 2;
	x[36] = 2; x[37] = 2 + 1e-9; x[38] = 2;
	/* the cube without faces 0 (0 2 1) and 2, 3 (the top) */
	for (i = 0; i < 12; i++)
		if (i != 0 && i != 2 && i != 3)
		{
			memcpy(t + 3 * nt, mr_cube_t + 3 * i, sizeof(unsigned int) * 3);
			nt++;
		}
	memcpy(t + 3 * nt, extra, sizeof extra);
	nt += sizeof extra / sizeof *extra / 3;
	memset(&d, 0, sizeof d);
	d.degenerate_tolerance = (qaws_scalar)1e-3;
	s = mr_repair(x, 13, t, NULL, nt, &d, &m, &r);
	if (m)
		mr_measure(m, &ms);
	else
		memset(&ms, 0, sizeof ms);
	sprintf(msg, "a cap, needles, a zero-area face, a vertex twice: %u caps, %u needles, %u dropped, %u faces, %u vertices, lowest %.4f, volume %.6f",
		r.caps, r.needles, r.dropped_faces, m ? qaws_mesh_3d_get_face_count(m) : 0, m ? qaws_mesh_3d_get_vertex_count(m) : 0, ms.min_height,
		m ? mr_volume(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && r.caps == 1 && r.needles >= 1 && r.dropped_faces >= 2 && r.holes_filled == 0
		&& qaws_mesh_3d_get_face_count(m) == 14 && qaws_mesh_3d_get_vertex_count(m) == 9 && ms.min_height > 1e-3 && fabs(mr_volume(m) - 1) < 1e-4, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_duplicates(void)
{
	/* the cube; face 0 twice; face 4 also turned the other way (the pair cancels: a hole, filled); a membrane apart, both ways */
	double x[33];
	unsigned int t[3 * 17], i;
	static unsigned int const extra[15] = { 0, 2, 1, 0, 5, 1, 8, 9, 10, 8, 10, 9, 0, 0, 0 };
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_status s;
	char msg[256];
	memcpy(x, mr_cube_x, sizeof mr_cube_x);
	x[24] = 0.5; x[25] = 0.5; x[26] = 0.2;
	x[27] = 0.6; x[28] = 0.5; x[29] = 0.2;
	x[30] = 0.5; x[31] = 0.6; x[32] = 0.2;
	memcpy(t, mr_cube_t, sizeof mr_cube_t);
	for (i = 0; i < 12; i++)
		t[36 + i] = extra[i];
	s = mr_repair(x, 11, t, NULL, 16, NULL, &m, &r);
	sprintf(msg, "duplicates: %u removed, %u filled, %u faces, volume %.6f", r.duplicate_faces, r.holes_filled, m ? qaws_mesh_3d_get_face_count(m) : 0,
		m ? mr_volume(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && r.duplicate_faces == 5 && r.holes_filled == 1 && r.fill_faces == 1
		&& qaws_mesh_3d_get_face_count(m) == 12 && fabs(mr_volume(m) - 1) < 1e-12, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_moebius_and_non_manifold(void)
{
	/* a Moebius strip: 24 steps round, 3 points across */
	double x[3 * 72 + 3 * 16];
	unsigned int t[3 * 2 * 2 * 24], nt = 0, i, j;
	static double const lo[2][3] = { { 0, 0, 0 }, { 1, 1, 0 } }, hi[2][3] = { { 1, 1, 1 }, { 2, 2, 1 } };
	unsigned int nv = 0;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_status s;
	char msg[256];
	for (i = 0; i < 24; i++)
		for (j = 0; j < 3; j++)
		{
			double u = 2 * M_PI * i / 24, v = 0.3 * ((double)j - 1);
			x[3 * (3 * i + j)] = (1 + v * cos(u / 2)) * cos(u);
			x[3 * (3 * i + j) + 1] = (1 + v * cos(u / 2)) * sin(u);
			x[3 * (3 * i + j) + 2] = v * sin(u / 2);
		}
	for (i = 0; i < 24; i++)
		for (j = 0; j < 2; j++)
		{
			/* the last step joins the first turned over: row j meets row 2 - j */
			unsigned int a = 3 * i + j, b = 3 * i + j + 1, c = i < 23 ? 3 * (i + 1) + j : 2 - j, d = i < 23 ? 3 * (i + 1) + j + 1 : 1 - j;
			unsigned int q[6] = { a, c, d, a, d, b };
			memcpy(t + 3 * nt, q, sizeof q);
			nt += 2;
		}
	s = mr_repair(x, 72, t, NULL, nt, NULL, &m, &r);
	sprintf(msg, "a Moebius strip: %u non-orientable, %u filled, %u left, closed %d", r.non_orientable, r.holes_filled, r.holes_left, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.non_orientable == 1 && r.holes_filled == 0 && r.holes_left == 1 && !r.closed
		&& qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, msg);
	qaws_mesh_3d_destroy(m);
	/* two cubes on one edge, as a soup */
	{
		double sx[3 * 72];
		unsigned int st[72];
		nt = 0;
		mr_box(lo[0], hi[0], x, &nv, t, &nt);
		mr_box(lo[1], hi[1], x, &nv, t, &nt);
		mr_soup(x, t, nt, 0, sx, st);
		s = mr_repair(sx, 72, st, NULL, nt, NULL, &m, &r);
		sprintf(msg, "two cubes on an edge: %u non-manifold, closed %d, %u vertices (%u split), volume %.6f", r.non_manifold_edges, r.closed,
			m ? qaws_mesh_3d_get_vertex_count(m) : 0, m ? qaws_mesh_3d_get_vertex_count(m) - 14 : 0, m ? mr_volume(m) : 0);
		printf("  %s\n", msg);
		TEST_ASSERT(s == QAWS_STATUS_OK && mr_closed(m, &r) && r.non_manifold_edges == 1 && qaws_mesh_3d_get_face_count(m) == 24
			&& fabs(mr_volume(m) - 2) < 1e-12, msg);
		qaws_mesh_3d_destroy(m);
	}
}

static void test_refusals(void)
{
	static double const x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	static unsigned int const bad[3] = { 0, 1, 3 }, ok[3] = { 0, 1, 2 }, two[1] = { 2 };
	double nanx[9];
	qaws_scalar px[9];
	qaws_mesh_3d* m = NULL;
	qaws_mesh_repair_report r;
	qaws_mesh_repair_desc d;
	qaws_status a, b, c, e, f, g;
	unsigned int i;
	for (i = 0; i < 9; i++)
		px[i] = (qaws_scalar)x[i];
	memcpy(nanx, x, sizeof x);
	nanx[4] = sqrt(-1.0 + 0 * nanx[0]);
	a = qaws_mesh_repair(NULL, &r, px, ok, NULL, 3, 1, NULL);
	b = mr_repair(x, 3, bad, NULL, 1, NULL, &m, &r);
	TEST_ASSERT(a == QAWS_STATUS_INVALID_ARGUMENT && b == QAWS_STATUS_INVALID_ARGUMENT && !m && r.invalid_face == 0, "refused: no output, an index out of range");
	c = mr_repair(nanx, 3, ok, NULL, 1, NULL, &m, &r);
	TEST_ASSERT(c == QAWS_STATUS_INVALID_ARGUMENT && r.invalid_vertex == 1 && r.invalid_face == QAWS_MESH_NONE, "refused: a position not finite");
	e = mr_repair(x, 3, ok, two, 1, NULL, &m, &r);
	memset(&d, 0, sizeof d);
	d.nested_shells = (qaws_nested_shells)7;
	f = mr_repair(x, 3, ok, NULL, 1, &d, &m, &r);
	g = mr_repair(x, 3, ok, NULL, 1, NULL, &m, NULL);
	TEST_ASSERT(e == QAWS_STATUS_INVALID_ARGUMENT && f == QAWS_STATUS_INVALID_ARGUMENT && g == QAWS_STATUS_OK && m
		&& qaws_mesh_3d_get_face_count(m) == 2, "refused: a face of 2 vertices, a policy out of range; one triangle filled to a closed pair");
	qaws_mesh_3d_destroy(m);
}

int test_139_mesh_repair_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 139: repair for farq\n");
	test_soups();
	test_near_duplicates();
	test_holes();
	test_gap();
	test_orientation();
	test_nested();
	test_degenerate();
	test_duplicates();
	test_moebius_and_non_manifold();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
