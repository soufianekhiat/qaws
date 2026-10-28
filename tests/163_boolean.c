/*
 * Test 163: mesh Booleans (qaws_mesh_boolean)
 *
 *   - two axis-aligned cubes crossing: union, intersection and difference
 *     closed, their volumes exact (2 - 0.336, 0.336, 0.664);
 *   - a turned cube, two cubes sharing two planes (coplanar faces), two
 *     spheres: closed, V(A u B) + V(A n B) = V(A) + V(B) and
 *     V(A - B) = V(A) - V(A n B) to rounding;
 *   - the self-union of the two cubes given as one soup: the union's volume;
 *     the union resolved again: nothing cut;
 *   - a cube inside another: the union is the outer cube, the difference a
 *     hollow cube (its inner faces turned); disjoint cubes: an empty
 *     intersection; cubes touching face to face: one solid;
 *   - cubes touching along an edge: the union closed (two shells meeting at
 *     the edge, each face paired with its own lobe's), its volume 2;
 *   - four thin solids crossing at shallow angles (a turning box's sweep):
 *     their self-union closed and one piece (patches sampled on the wrong
 *     side of a surface nearer than the sample mended from their
 *     neighbours round the edges where they meet);
 *   - coplanar overlaps flagged: B turned about z sharing A's floor and ceiling,
 *     the union's and intersection's "overlap" faces (all A's) cover 2 V(A n B);
 *     none for the difference or two cubes crossing;
 *   - the two spheres united on 4 threads: the serial result (faces, volume);
 *   - refusals: no faces, an unknown operation.
 */

#include "test_common.h"
#include "qaws_mesh_boolean.h"
#include "test_threads.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static void bo_box(double const* lo, double const* hi, double const* axis, double angle, qaws_scalar* x, unsigned int* nv, unsigned int* t,
	unsigned int* nt)
{
	static unsigned int const bt[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5 };
	static int const bx[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
	double c[3], co = cos(angle), si = sin(angle);
	unsigned int i, k;
	for (k = 0; k < 3; k++)
		c[k] = 0.5 * (lo[k] + hi[k]);
	for (i = 0; i < 8; i++)
	{
		double p[3], q[3], d;
		for (k = 0; k < 3; k++)
			p[k] = (bx[3 * i + k] ? hi[k] : lo[k]) - c[k];
		d = axis[0] * p[0] + axis[1] * p[1] + axis[2] * p[2];
		q[0] = p[0] * co + (axis[1] * p[2] - axis[2] * p[1]) * si + axis[0] * d * (1 - co);
		q[1] = p[1] * co + (axis[2] * p[0] - axis[0] * p[2]) * si + axis[1] * d * (1 - co);
		q[2] = p[2] * co + (axis[0] * p[1] - axis[1] * p[0]) * si + axis[2] * d * (1 - co);
		for (k = 0; k < 3; k++)
			x[3 * (*nv + i) + k] = (qaws_scalar)(q[k] + c[k]);
	}
	for (i = 0; i < 36; i++)
		t[3 * *nt + i] = *nv + bt[i];
	*nv += 8;
	*nt += 12;
}

static void bo_sphere(double const* c, double rad, unsigned int na, unsigned int nb, qaws_scalar* x, unsigned int* nv, unsigned int* t, unsigned int* nt)
{
	unsigned int i, j, v0 = *nv;
	for (i = 0; i < 2; i++)
	{
		x[3 * (v0 + i)] = (qaws_scalar)c[0];
		x[3 * (v0 + i) + 1] = (qaws_scalar)c[1];
		x[3 * (v0 + i) + 2] = (qaws_scalar)(c[2] + (i ? -rad : rad));
	}
	for (j = 1; j < nb; j++)
		for (i = 0; i < na; i++)
		{
			double th = M_PI * j / nb, ph = 2 * M_PI * i / na;
			unsigned int v = v0 + 2 + (j - 1) * na + i;
			x[3 * v] = (qaws_scalar)(c[0] + rad * sin(th) * cos(ph));
			x[3 * v + 1] = (qaws_scalar)(c[1] + rad * sin(th) * sin(ph));
			x[3 * v + 2] = (qaws_scalar)(c[2] + rad * cos(th));
		}
	for (i = 0; i < na; i++)
	{
		unsigned int* f = t + 3 * *nt;
		f[0] = v0; f[1] = v0 + 2 + i; f[2] = v0 + 2 + (i + 1) % na;
		f[3] = v0 + 1; f[4] = v0 + 2 + (nb - 2) * na + (i + 1) % na; f[5] = v0 + 2 + (nb - 2) * na + i;
		*nt += 2;
	}
	for (j = 0; j + 2 < nb; j++)
		for (i = 0; i < na; i++)
		{
			unsigned int a = v0 + 2 + j * na + i, b = v0 + 2 + j * na + (i + 1) % na, cc = a + na, d = b + na;
			unsigned int* f = t + 3 * *nt;
			f[0] = a; f[1] = cc; f[2] = d;
			f[3] = a; f[4] = d; f[5] = b;
			*nt += 2;
		}
	*nv += 2 + na * (nb - 1);
}

static double bo_volume_arrays(qaws_scalar const* p, unsigned int const* t, unsigned int nt)
{
	double v = 0;
	unsigned int i;
	for (i = 0; i < nt; i++)
	{
		double const a[3] = { p[3 * t[3 * i]], p[3 * t[3 * i] + 1], p[3 * t[3 * i] + 2] };
		double const b[3] = { p[3 * t[3 * i + 1]], p[3 * t[3 * i + 1] + 1], p[3 * t[3 * i + 1] + 2] };
		double const c[3] = { p[3 * t[3 * i + 2]], p[3 * t[3 * i + 2] + 1], p[3 * t[3 * i + 2] + 2] };
		v += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
	}
	return v;
}

static double bo_volume(qaws_mesh_3d const* m)
{
	unsigned int nv, nt;
	qaws_scalar* p;
	unsigned int* t;
	double v;
	if (!m)
		return 0;
	nv = qaws_mesh_3d_get_vertex_count(m);
	nt = qaws_mesh_3d_get_triangle_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	qaws_mesh_3d_get_positions(p, m);
	qaws_mesh_3d_get_triangles(t, m);
	v = bo_volume_arrays(p, t, nt);
	free(p);
	free(t);
	return v;
}

static qaws_mesh_3d_desc bo_desc(qaws_scalar const* x, unsigned int nv, unsigned int const* t, unsigned int nt)
{
	qaws_mesh_3d_desc d;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = t;
	d.face_count = nt;
	return d;
}

/* the three operations on A and B: closed, nothing left crossing; their volumes out */
static int bo_ops(char const* name, qaws_mesh_3d_desc const* A, qaws_mesh_3d_desc const* B, double* vol)
{
	static char const* const op_name[3] = { "union", "intersection", "difference" };
	int ok = 1, k;
	for (k = 0; k < 3; k++)
	{
		qaws_boolean_report r;
		qaws_mesh_3d* m = NULL;
		qaws_status s = qaws_mesh_boolean(&m, &r, A, B, (qaws_boolean_op)k, NULL);
		char msg[300];
		/* (one sphere each: an edge made twice by the turning in flat regions would split vertices) */
		int chi = m ? (int)qaws_mesh_3d_get_vertex_count(m) - (int)qaws_mesh_3d_get_edge_count(m) + (int)qaws_mesh_3d_get_face_count(m) : 0;
		vol[k] = bo_volume(m);
		sprintf(msg, "%s, %s: %u pieces kept (%u turned), %u dropped, %u edges turned flat, closed %d, V - E + F = %d, %u crossing left, volume %.15f",
			name, op_name[k], r.kept, r.flipped, r.dropped, r.flat_flips, r.closed, chi, r.intersect.remaining, vol[k]);
		printf("  %s\n", msg);
		TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && chi == 2 && r.intersect.remaining == 0 && r.intersect.unresolved == 0, msg);
		ok &= s == QAWS_STATUS_OK && m != NULL;
		qaws_mesh_3d_destroy(m);
	}
	return ok;
}

static double const bo_z[3] = { 0, 0, 1 };

static void test_cubes(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 0.4, 0.3, 0.2 }, hi1[3] = { 1.4, 1.3, 1.2 };
	qaws_scalar x[48];
	unsigned int t[72], nv = 0, nt = 0;
	qaws_mesh_3d_desc A, B;
	double v[3], w = 0.6 * 0.7 * 0.8;
	char msg[200];
	bo_box(lo0, hi0, bo_z, 0, x, &nv, t, &nt);
	bo_box(lo1, hi1, bo_z, 0, x, &nv, t, &nt);
	A = bo_desc(x, 8, t, 12);
	B = bo_desc(x + 24, 8, t, 12);   /* (B's indices are 0 to 7 too: its own vertices) */
	if (bo_ops("two cubes", &A, &B, v))
	{
		sprintf(msg, "two cubes: volumes %.2e, %.2e, %.2e off", v[0] - (2 - w), v[1] - w, v[2] - (1 - w));
		printf("  %s\n", msg);
		TEST_ASSERT(fabs(v[0] - (2 - w)) < 1e-12 && fabs(v[1] - w) < 1e-12 && fabs(v[2] - (1 - w)) < 1e-12, msg);
	}
	/* the self-union of both as one soup */
	{
		qaws_mesh_3d_desc S = bo_desc(x, 16, t, 24);
		qaws_boolean_report r;
		qaws_mesh_3d* m = NULL;
		qaws_status s = qaws_mesh_boolean(&m, &r, &S, NULL, QAWS_BOOLEAN_UNION, NULL);
		double vs = bo_volume(m);
		sprintf(msg, "self-union of the soup: closed %d, volume %.15f (union %.15f)", r.closed, vs, 2 - w);
		printf("  %s\n", msg);
		TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && fabs(vs - (2 - w)) < 1e-12, msg);
		/* the result resolved again: nothing to cut */
		if (m)
		{
			unsigned int rn = qaws_mesh_3d_get_vertex_count(m), rt = qaws_mesh_3d_get_triangle_count(m);
			qaws_scalar* rp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (rn + 1));
			unsigned int* rtri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (rt + 1));
			qaws_intersect_report ir;
			qaws_mesh_3d* again = NULL;
			qaws_mesh_3d_get_positions(rp, m);
			qaws_mesh_3d_get_triangles(rtri, m);
			s = qaws_mesh_resolve_intersections(&again, &ir, rp, rtri, NULL, rn, rt, NULL);
			sprintf(msg, "the union resolved again: %u pairs, %u points added, %u faces cut", ir.pairs, ir.points, ir.split_faces);
			printf("  %s\n", msg);
			TEST_ASSERT(s == QAWS_STATUS_OK && ir.pairs == 0 && ir.points == 0 && ir.split_faces == 0, msg);
			qaws_mesh_3d_destroy(again);
			free(rp);
			free(rtri);
		}
		qaws_mesh_3d_destroy(m);
	}
}

static void bo_identities(char const* name, qaws_mesh_3d_desc const* A, qaws_mesh_3d_desc const* B)
{
	double v[3], va = bo_volume_arrays(A->positions, A->face_vertices, A->face_count), vb = bo_volume_arrays(B->positions, B->face_vertices, B->face_count);
	char msg[300];
	if (!bo_ops(name, A, B, v))
		return;
	sprintf(msg, "%s: V(A u B) + V(A n B) - V(A) - V(B) = %.2e, V(A - B) - V(A) + V(A n B) = %.2e, V(A n B) = %.6f", name, v[0] + v[1] - va - vb,
		v[2] - va + v[1], v[1]);
	printf("  %s\n", msg);
	TEST_ASSERT(fabs(v[0] + v[1] - va - vb) < 1e-12 * (va + vb) && fabs(v[2] - va + v[1]) < 1e-12 * va && v[1] > 0, msg);
}

static void test_turned(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 0.45, 0.35, 0.25 }, hi1[3] = { 1.45, 1.35, 1.25 };
	static double const lo2[3] = { 0.5, 0.25, 0 }, hi2[3] = { 1.5, 1.25, 1 };
	double ax[3] = { 1 / sqrt(3.0), 1 / sqrt(3.0), 1 / sqrt(3.0) };
	qaws_scalar x[48];
	unsigned int t[72], nv = 0, nt = 0;
	qaws_mesh_3d_desc A, B;
	bo_box(lo0, hi0, bo_z, 0, x, &nv, t, &nt);
	bo_box(lo1, hi1, ax, 0.3, x, &nv, t, &nt);
	A = bo_desc(x, 8, t, 12);
	B = bo_desc(x + 24, 8, t, 12);
	bo_identities("a turned cube", &A, &B);
	nv = 8;
	nt = 12;
	bo_box(lo2, hi2, bo_z, 0.4, x, &nv, t, &nt);
	bo_identities("cubes sharing two planes", &A, &B);
}

/* the area of m's faces whose "overlap" is 1, their count, all from A (source below na); -1 without the attribute */
static double bo_overlap_area(qaws_mesh_3d const* m, unsigned int na, unsigned int* count, int* from_a)
{
	unsigned int id, sid, nf = m ? qaws_mesh_3d_get_triangle_count(m) : 0, nv = m ? qaws_mesh_3d_get_vertex_count(m) : 0, f, k;
	qaws_scalar *ov, *src, *p;
	unsigned int* t;
	double area = 0;
	*count = 0;
	*from_a = 1;
	if (!m || qaws_mesh_3d_find_attribute(&id, m, "overlap") != QAWS_STATUS_OK || qaws_mesh_3d_find_attribute(&sid, m, "source") != QAWS_STATUS_OK)
		return -1;
	ov = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nf + 1));
	src = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nf + 1));
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nf + 1));
	qaws_mesh_3d_get_attribute(ov, m, id);
	qaws_mesh_3d_get_attribute(src, m, sid);
	qaws_mesh_3d_get_positions(p, m);
	qaws_mesh_3d_get_triangles(t, m);
	for (f = 0; f < nf; f++)
	{
		double u[3], w[3];
		if (ov[f] != 1)
			continue;
		for (k = 0; k < 3; k++)
		{
			u[k] = (double)p[3 * t[3 * f + 1] + k] - (double)p[3 * t[3 * f] + k];
			w[k] = (double)p[3 * t[3 * f + 2] + k] - (double)p[3 * t[3 * f] + k];
		}
		area += 0.5 * sqrt(pow(u[1] * w[2] - u[2] * w[1], 2) + pow(u[2] * w[0] - u[0] * w[2], 2) + pow(u[0] * w[1] - u[1] * w[0], 2));
		(*count)++;
		*from_a &= (unsigned int)src[f] < na;
	}
	free(ov); free(src); free(p); free(t);
	return area;
}

/* coplanar overlaps flagged: B turned about z sharing A's floor and ceiling; the union's and the intersection's faces
   standing for both ("overlap" 1, all A's) cover 2 V(A n B) (each plane's overlap, the solids being prisms of height 1);
   the difference has none (the attribute absent), nor two cubes crossing without a shared plane */
static void test_overlap(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo2[3] = { 0.5, 0.25, 0 }, hi2[3] = { 1.5, 1.25, 1 };
	static double const lo1[3] = { 0.4, 0.3, 0.2 }, hi1[3] = { 1.4, 1.3, 1.2 };
	static qaws_boolean_op const ops[3] = { QAWS_BOOLEAN_UNION, QAWS_BOOLEAN_INTERSECTION, QAWS_BOOLEAN_DIFFERENCE };
	qaws_scalar x[48];
	unsigned int t[72], nv = 0, nt = 0, i, cnt[3] = { 0, 0, 0 }, rep[3] = { 0, 0, 0 }, crossing_rep = 1;
	qaws_mesh_3d_desc A, B;
	double ar[3], vi = 0, crossing = 0;
	int fa[3] = { 1, 1, 1 }, fx = 1;
	char msg[400];
	bo_box(lo0, hi0, bo_z, 0, x, &nv, t, &nt);
	bo_box(lo2, hi2, bo_z, 0.4, x, &nv, t, &nt);
	A = bo_desc(x, 8, t, 12);
	B = bo_desc(x + 24, 8, t, 12);
	for (i = 0; i < 3; i++)
	{
		qaws_boolean_report r;
		qaws_mesh_3d* m = NULL;
		qaws_mesh_boolean(&m, &r, &A, &B, ops[i], NULL);
		if (i == 1)
			vi = bo_volume(m);
		ar[i] = bo_overlap_area(m, 12, &cnt[i], &fa[i]);
		rep[i] = r.overlapped;
		qaws_mesh_3d_destroy(m);
	}
	{
		/* (crossing cubes: no shared plane) */
		qaws_boolean_report r;
		qaws_mesh_3d* m = NULL;
		unsigned int c2;
		nv = nt = 0;
		bo_box(lo0, hi0, bo_z, 0, x, &nv, t, &nt);
		bo_box(lo1, hi1, bo_z, 0, x, &nv, t, &nt);
		A = bo_desc(x, 8, t, 12);
		B = bo_desc(x + 24, 8, t, 12);
		qaws_mesh_boolean(&m, &r, &A, &B, QAWS_BOOLEAN_UNION, NULL);
		crossing = bo_overlap_area(m, 12, &c2, &fx);
		crossing_rep = r.overlapped;
		qaws_mesh_3d_destroy(m);
	}
	sprintf(msg, "overlap flagged: union %u faces (%u reported) of area %.15f, intersection %u (%u) of %.15f, 2 V(A n B) %.15f, all A's %d %d; "
		"difference %.0f (%u), crossing cubes %.0f (%u)", cnt[0], rep[0], ar[0], cnt[1], rep[1], ar[1], 2 * vi, fa[0], fa[1], ar[2], rep[2], crossing, crossing_rep);
	printf("  %s\n", msg);
	TEST_ASSERT(vi > 0 && cnt[0] && cnt[0] == rep[0] && cnt[1] == rep[1] && fabs(ar[0] - 2 * vi) < 1e-12 && fabs(ar[1] - 2 * vi) < 1e-12 && fa[0] && fa[1]
		&& ar[2] == -1 && rep[2] == 0 && crossing == -1 && crossing_rep == 0, msg);
}

/* an executor on host threads: count split in n ranges, one thread each */
typedef struct bo_part
{
	qaws_batch_task_fn task;
	void* ctx;
	unsigned int b, e;
} bo_part;

static void bo_part_run(void* p)
{
	bo_part* q = (bo_part*)p;
	if (q->e > q->b)
		q->task(q->ctx, q->b, q->e);
}

static void bo_parallel_for(void* user, unsigned int count, qaws_batch_task_fn task, void* ctx)
{
	unsigned int n = (unsigned int)(size_t)user, k;
	bo_part parts[TEST_MAX_THREADS];
	void* args[TEST_MAX_THREADS];
	for (k = 0; k < n; k++)
	{
		parts[k].task = task;
		parts[k].ctx = ctx;
		parts[k].b = (unsigned int)((unsigned long long)count * k / n);
		parts[k].e = (unsigned int)((unsigned long long)count * (k + 1) / n);
		args[k] = &parts[k];
	}
	test_run_threads(n, bo_part_run, args);
}

static void test_spheres(void)
{
	static double const c0[3] = { 0, 0, 0 }, c1[3] = { 0.7, 0.25, 0.15 };
	unsigned int na = 24, nb = 12, nv = 0, nt = 0, nv0, nt0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (2 + na * nb));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 12 * na * nb);
	qaws_mesh_3d_desc A, B;
	bo_sphere(c0, 1, na, nb, x, &nv, t, &nt);
	nv0 = nv;
	nt0 = nt;
	bo_sphere(c0, 0.8, na, nb, x, &nv, t, &nt);
	{
		unsigned int i;
		for (i = nv0; i < nv; i++)   /* moved: B's vertices from 0 */
		{
			x[3 * i] += (qaws_scalar)c1[0];
			x[3 * i + 1] += (qaws_scalar)c1[1];
			x[3 * i + 2] += (qaws_scalar)c1[2];
		}
		for (i = 3 * nt0; i < 3 * nt; i++)
			t[i] -= nv0;
	}
	A = bo_desc(x, nv0, t, nt0);
	B = bo_desc(x + 3 * nv0, nv - nv0, t + 3 * nt0, nt - nt0);
	bo_identities("two spheres", &A, &B);
	{
		/* the same union on four threads (the operands' BVHs built side by side): the same result */
		qaws_batch_executor ex;
		qaws_intersect_desc id;
		qaws_mesh_3d *s = NULL, *p = NULL;
		double vs, vp;
		char msg[200];
		ex.parallel_for = bo_parallel_for;
		ex.user = (void*)(size_t)4;
		memset(&id, 0, sizeof id);
		id.executor = &ex;
		qaws_mesh_boolean(&s, NULL, &A, &B, QAWS_BOOLEAN_UNION, NULL);
		qaws_mesh_boolean(&p, NULL, &A, &B, QAWS_BOOLEAN_UNION, &id);
		vs = bo_volume(s);
		vp = bo_volume(p);
		sprintf(msg, "two spheres united serially and on 4 threads: %u and %u faces, volumes %.15f and %.15f", s ? qaws_mesh_3d_get_face_count(s) : 0,
			p ? qaws_mesh_3d_get_face_count(p) : 0, vs, vp);
		printf("  %s\n", msg);
		TEST_ASSERT(s && p && qaws_mesh_3d_get_face_count(s) == qaws_mesh_3d_get_face_count(p) && vs == vp, msg);
		qaws_mesh_3d_destroy(s);
		qaws_mesh_3d_destroy(p);
	}
	free(x);
	free(t);
}

static void test_nested(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 0.25, 0.25, 0.25 }, hi1[3] = { 0.75, 0.75, 0.75 };
	static double const lo2[3] = { 2, 0, 0 }, hi2[3] = { 3, 1, 1 };
	qaws_scalar x[48];
	unsigned int t[72], nv = 0, nt = 0;
	qaws_mesh_3d_desc A, B;
	qaws_boolean_report r;
	qaws_mesh_3d* m = NULL;
	qaws_status s;
	char msg[200];
	double v;
	bo_box(lo0, hi0, bo_z, 0, x, &nv, t, &nt);
	bo_box(lo1, hi1, bo_z, 0, x, &nv, t, &nt);
	A = bo_desc(x, 8, t, 12);
	B = bo_desc(x + 24, 8, t, 12);
	s = qaws_mesh_boolean(&m, &r, &A, &B, QAWS_BOOLEAN_UNION, NULL);
	v = bo_volume(m);
	sprintf(msg, "a cube inside another, union: %u faces, volume %.15f", m ? qaws_mesh_3d_get_triangle_count(m) : 0, v);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && qaws_mesh_3d_get_triangle_count(m) == 12 && fabs(v - 1) < 1e-14, msg);
	qaws_mesh_3d_destroy(m);
	m = NULL;
	s = qaws_mesh_boolean(&m, &r, &A, &B, QAWS_BOOLEAN_DIFFERENCE, NULL);
	v = bo_volume(m);
	sprintf(msg, "a cube inside another, difference: %u faces (%u turned), volume %.15f", m ? qaws_mesh_3d_get_triangle_count(m) : 0, r.flipped, v);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && r.flipped == 12 && fabs(v - 0.875) < 1e-14, msg);
	qaws_mesh_3d_destroy(m);
	m = NULL;
	nv = 8;
	nt = 12;
	bo_box(lo2, hi2, bo_z, 0, x, &nv, t, &nt);
	s = qaws_mesh_boolean(&m, &r, &A, &B, QAWS_BOOLEAN_INTERSECTION, NULL);
	sprintf(msg, "disjoint cubes, intersection: status %d, %s", (int)s, m ? "a mesh" : "empty");
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && !m && r.kept == 0, msg);
	qaws_mesh_3d_destroy(m);
	/* touching face to face: the union one solid, the faces between gone */
	{
		static double const lo3[3] = { 1, 0.2, 0.3 }, hi3[3] = { 2, 0.8, 0.7 };
		m = NULL;
		nv = 8;
		nt = 12;
		bo_box(lo3, hi3, bo_z, 0, x, &nv, t, &nt);
		s = qaws_mesh_boolean(&m, &r, &A, &B, QAWS_BOOLEAN_UNION, NULL);
		v = bo_volume(m);
		sprintf(msg, "cubes touching face to face, union: closed %d, volume %.15f (1.24)", r.closed, v);
		printf("  %s\n", msg);
		TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && fabs(v - 1.24) < 1e-14, msg);
		qaws_mesh_3d_destroy(m);
	}
}

/* four convex hulls of a turning box's faces over three slabs (the sweep of figure 60), crossing at shallow angles */
static double const bo_sliver_x[48] = {
	1.6873215416209051, 1.0783997473952096, 0,
	1.6334372557452375, 1.1583966209996763, 0,
	1.6334372557452375, 1.1583966209996763, 0.59999999999999998,
	1.6873215416209051, 1.0783997473952096, 0.59999999999999998,
	1.509739651377463, 1.3155554663558116, 0.59999999999999998,
	1.5713349257130493, 1.2413325707619072, 0.59999999999999998,
	1.5753593163095529, 1.2362212684291991, 0,
	1.5753593163095529, 1.2362212684291991, 0.59999999999999998,
	1.4441883203749217, 1.387198650263423, 0.59999999999999998,
	1.4441883203749217, 1.387198650263423, 0,
	0.68930141120380317, 0.73134367059026739, 0.59999999999999998,
	0.65277249701745454, 0.76412568804981107, 0.59999999999999998,
	1.3747379087228364, 1.4560548349290838, 0,
	1.3747379087228364, 1.4560548349290838, 0.59999999999999998,
	1.5131237380986911, 1.3116617525880077, 0,
	1.5131237380986911, 1.3116617525880077, 0.59999999999999998 };
static unsigned int const bo_sliver_t[96] = {
	0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1, 3, 2, 5, 3, 5, 4, 5, 1, 4, 5, 2, 1,
	1, 6, 7, 1, 7, 2, 1, 2, 8, 1, 8, 6, 2, 7, 4, 2, 4, 8, 4, 6, 8, 4, 7, 6,
	9, 8, 10, 9, 10, 11, 9, 12, 13, 9, 13, 8, 9, 11, 12, 8, 13, 11, 8, 11, 10, 12, 11, 13,
	6, 14, 15, 6, 15, 7, 6, 7, 13, 6, 13, 14, 7, 15, 8, 7, 8, 13, 8, 14, 13, 8, 15, 14
};

static void test_touching(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 1, 1, 0 }, hi1[3] = { 2, 2, 1 };
	qaws_scalar x[48], xs[48];
	unsigned int t[72], nv = 0, nt = 0, i;
	qaws_mesh_3d_desc A, B, S;
	qaws_mesh_3d* m = NULL;
	qaws_boolean_report r;
	qaws_mesh_topology tp;
	qaws_status s;
	double v;
	char msg[300];
	bo_box(lo0, hi0, bo_z, 0, x, &nv, t, &nt);
	bo_box(lo1, hi1, bo_z, 0, x, &nv, t, &nt);
	A = bo_desc(x, 8, t, 12);
	B = bo_desc(x + 24, 8, t, 12);
	memset(&tp, 0, sizeof tp);
	s = qaws_mesh_boolean(&m, &r, &A, &B, QAWS_BOOLEAN_UNION, NULL);
	v = bo_volume(m);
	if (m)
		qaws_mesh_3d_get_topology(&tp, m);
	sprintf(msg, "cubes touching along an edge, union: closed %d (%u boundary loops, %u shells), volume %.15f (2)", r.closed, tp.boundary_loops, tp.components, v);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && tp.closed && tp.components == 2 && fabs(v - 2) < 1e-14, msg);
	qaws_mesh_3d_destroy(m);
	m = NULL;
	for (i = 0; i < 48; i++)
		xs[i] = (qaws_scalar)bo_sliver_x[i];
	S = bo_desc(xs, 16, bo_sliver_t, 32);
	memset(&tp, 0, sizeof tp);
	s = qaws_mesh_boolean(&m, &r, &S, NULL, QAWS_BOOLEAN_UNION, NULL);
	if (m)
		qaws_mesh_3d_get_topology(&tp, m);
	sprintf(msg, "four thin solids crossing at shallow angles, self-union: closed %d (%u boundary loops), %u piece(s), volume %.6f", r.closed, tp.boundary_loops,
		tp.components, bo_volume(m));
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.closed && tp.closed && tp.components == 1, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_refusals(void)
{
	qaws_scalar x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	unsigned int t[3] = { 0, 1, 2 };
	qaws_mesh_3d_desc A = bo_desc(x, 3, t, 1), E = bo_desc(x, 3, t, 0);
	qaws_mesh_3d* m = NULL;
	TEST_ASSERT(qaws_mesh_boolean(&m, NULL, &E, &A, QAWS_BOOLEAN_UNION, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m, "no faces refused");
	TEST_ASSERT(qaws_mesh_boolean(&m, NULL, &A, &A, (qaws_boolean_op)7, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m, "unknown operation refused");
}

int test_163_boolean_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 163: mesh Booleans\n");
	test_cubes();
	test_turned();
	test_spheres();
	test_nested();
	test_touching();
	test_overlap();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
