/*
 * Test 162: self-intersections resolved numerically (qaws_mesh_resolve_intersections)
 *
 *   - two cubes crossing at an angle, three cubes crossing (their cuts
 *     meeting at triple points), two cubes sharing two planes (coplanar
 *     overlaps), two spheres crossing: no pair of result faces crosses
 *     another, every cut kept, the cuts conforming on both sides (the
 *     result's directed edges cancel, vertices matched by position), the
 *     area unchanged, every face's source an input face;
 *     (and nearly sharing them, 1e-12 off: snapped together);
 *   - a cube alone: unchanged;
 *   - refusals: no faces, an index out of range, a negative tolerance.
 */

#include "test_common.h"
#include "qaws_mesh_intersect.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static void ix_box(double const* lo, double const* hi, double const* axis, double angle, qaws_scalar* x, unsigned int* nv, unsigned int* t,
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
		/* Rodrigues about the unit axis */
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

static double ix_area(qaws_scalar const* x, unsigned int const* t, unsigned int nt)
{
	double a = 0;
	unsigned int i, k;
	for (i = 0; i < nt; i++)
	{
		double u[3], v[3], n[3];
		for (k = 0; k < 3; k++)
		{
			u[k] = (double)x[3 * t[3 * i + 1] + k] - (double)x[3 * t[3 * i] + k];
			v[k] = (double)x[3 * t[3 * i + 2] + k] - (double)x[3 * t[3 * i] + k];
		}
		n[0] = u[1] * v[2] - u[2] * v[1];
		n[1] = u[2] * v[0] - u[0] * v[2];
		n[2] = u[0] * v[1] - u[1] * v[0];
		a += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	}
	return a;
}

typedef struct ix_edge
{
	unsigned int a, b;
	int s;
} ix_edge;

static int ix_cmp(void const* p, void const* q)
{
	ix_edge const *a = (ix_edge const*)p, *b = (ix_edge const*)q;
	if (a->a != b->a)
		return a->a < b->a ? -1 : 1;
	return a->b < b->b ? -1 : a->b > b->b;
}

/* directed edges not cancelled, vertices matched by position */
static unsigned int ix_open_edges(qaws_scalar const* x, unsigned int nv, unsigned int const* t, unsigned int nt)
{
	unsigned int* id = (unsigned int*)malloc(sizeof(unsigned int) * (nv + 1));
	ix_edge* e = (ix_edge*)malloc(sizeof(ix_edge) * 3 * (nt + 1));
	unsigned int i, j, open = 0;
	for (i = 0; i < nv; i++)
	{
		id[i] = i;
		for (j = 0; j < i; j++)
			if (x[3 * i] == x[3 * j] && x[3 * i + 1] == x[3 * j + 1] && x[3 * i + 2] == x[3 * j + 2])
			{
				id[i] = id[j];
				break;
			}
	}
	for (i = 0; i < nt; i++)
		for (j = 0; j < 3; j++)
		{
			unsigned int a = id[t[3 * i + j]], b = id[t[3 * i + (j + 1) % 3]];
			e[3 * i + j].a = a < b ? a : b;
			e[3 * i + j].b = a < b ? b : a;
			e[3 * i + j].s = a < b ? 1 : -1;
		}
	qsort(e, 3 * nt, sizeof(ix_edge), ix_cmp);
	for (i = 0; i < 3 * nt; i = j)
	{
		int s = 0;
		for (j = i; j < 3 * nt && e[j].a == e[i].a && e[j].b == e[i].b; j++)
			s += e[j].s;
		open += s != 0;
	}
	free(id);
	free(e);
	return open;
}

static void ix_check(char const* name, qaws_scalar const* x, unsigned int nv, unsigned int const* t, unsigned int nt, int coplanar)
{
	qaws_intersect_report r;
	qaws_mesh_3d* m = NULL;
	qaws_status s = qaws_mesh_resolve_intersections(&m, &r, x, t, NULL, nv, nt, NULL);
	char msg[320];
	if (s == QAWS_STATUS_OK && m)
	{
		unsigned int n = qaws_mesh_3d_get_vertex_count(m), ntr = qaws_mesh_3d_get_triangle_count(m), attr = 0, i, bad = 0, open;
		qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (n + 1));
		qaws_scalar* src = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ntr + 1));
		unsigned int* tr = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (ntr + 1));
		double a0 = ix_area(x, t, nt), a1;
		qaws_mesh_3d_get_positions(p, m);
		qaws_mesh_3d_get_triangles(tr, m);
		a1 = ix_area(p, tr, ntr);
		open = ix_open_edges(p, n, tr, ntr);
		if (qaws_mesh_3d_find_attribute(&attr, m, "source") == QAWS_STATUS_OK && qaws_mesh_3d_get_attribute(src, m, attr) == QAWS_STATUS_OK)
		{
			for (i = 0; i < ntr; i++)
				bad += !(src[i] >= 0 && src[i] < (qaws_scalar)nt);
		}
		else
			bad = 1;
		sprintf(msg, "%s: %u pairs (%u coplanar), %u points added, %u faces cut, %u faces, %u crossing left, %u unresolved, %u open edges, area %.15f (input %.15f)",
			name, r.pairs, r.coplanar_pairs, r.points, r.split_faces, ntr, r.remaining, r.unresolved, open, a1, a0);
		printf("  %s\n", msg);
		TEST_ASSERT(r.pairs > 0 && r.points > 0 && r.split_faces > 0 && (!coplanar || r.coplanar_pairs > 0), msg);
		TEST_ASSERT(r.remaining == 0 && r.unresolved == 0 && open == 0, msg);
		TEST_ASSERT(fabs(a1 - a0) < 1e-12 * a0 && !bad, msg);
		free(p);
		free(src);
		free(tr);
	}
	else
	{
		sprintf(msg, "%s: status %d", name, (int)s);
		TEST_ASSERT(0, msg);
	}
	qaws_mesh_3d_destroy(m);
}

static double const ix_z[3] = { 0, 0, 1 };

static void test_cubes(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 0.45, 0.35, 0.25 }, hi1[3] = { 1.45, 1.35, 1.25 };
	static double const lo2[3] = { -0.3, 0.4, 0.5 }, hi2[3] = { 0.7, 1.4, 1.5 };
	double ax[3] = { 1 / sqrt(3.0), 1 / sqrt(3.0), 1 / sqrt(3.0) }, ay[3] = { 0.6, 0, 0.8 };
	qaws_scalar x[72];
	unsigned int t[108], nv = 0, nt = 0;
	ix_box(lo0, hi0, ix_z, 0, x, &nv, t, &nt);
	ix_box(lo1, hi1, ax, 0.3, x, &nv, t, &nt);
	ix_check("two cubes crossing", x, nv, t, nt, 0);
	ix_box(lo2, hi2, ay, 0.45, x, &nv, t, &nt);
	ix_check("three cubes crossing", x, nv, t, nt, 0);
}

static void test_coplanar(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 }, lo1[3] = { 0.5, 0.25, 0 }, hi1[3] = { 1.5, 1.25, 1 };
	qaws_scalar x[48];
	unsigned int t[72], nv = 0, nt = 0;
	ix_box(lo0, hi0, ix_z, 0, x, &nv, t, &nt);
	ix_box(lo1, hi1, ix_z, 0.4, x, &nv, t, &nt);   /* turned about z: its top and bottom in the first's planes */
	ix_check("two cubes sharing two planes", x, nv, t, nt, 1);
	/* nearly: the second's planes 1e-12 off the first's, tilted by 1e-12 */
	{
		static double const lo2[3] = { 0.5, 0.25, 1e-12 }, hi2[3] = { 1.5, 1.25, 1 + 1e-12 }, ax[3] = { 1, 0, 0 };
		nv = 8;
		nt = 12;
		ix_box(lo2, hi2, ax, 1e-12, x, &nv, t, &nt);
		ix_check("two cubes nearly sharing two planes", x, nv, t, nt, 1);
	}
}

static void ix_sphere(double const* c, double rad, unsigned int na, unsigned int nb, qaws_scalar* x, unsigned int* nv, unsigned int* t, unsigned int* nt)
{
	unsigned int i, j, v0 = *nv;
	x[3 * v0] = (qaws_scalar)c[0];
	x[3 * v0 + 1] = (qaws_scalar)c[1];
	x[3 * v0 + 2] = (qaws_scalar)(c[2] + rad);
	x[3 * v0 + 3] = (qaws_scalar)c[0];
	x[3 * v0 + 4] = (qaws_scalar)c[1];
	x[3 * v0 + 5] = (qaws_scalar)(c[2] - rad);
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

static void test_spheres(void)
{
	static double const c0[3] = { 0, 0, 0 }, c1[3] = { 0.7, 0.25, 0.15 };
	unsigned int na = 24, nb = 12, nv = 0, nt = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * (2 + na * nb));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 12 * na * nb);
	ix_sphere(c0, 1, na, nb, x, &nv, t, &nt);
	ix_sphere(c1, 0.8, na, nb, x, &nv, t, &nt);
	ix_check("two spheres crossing", x, nv, t, nt, 0);
	free(x);
	free(t);
}

static void test_alone(void)
{
	static double const lo0[3] = { 0, 0, 0 }, hi0[3] = { 1, 1, 1 };
	qaws_scalar x[24];
	unsigned int t[36], nv = 0, nt = 0;
	qaws_intersect_report r;
	qaws_mesh_3d* m = NULL;
	qaws_status s;
	char msg[200];
	ix_box(lo0, hi0, ix_z, 0, x, &nv, t, &nt);
	s = qaws_mesh_resolve_intersections(&m, &r, x, t, NULL, nv, nt, NULL);
	sprintf(msg, "a cube alone: %u pairs, %u faces", r.pairs, m ? qaws_mesh_3d_get_triangle_count(m) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(s == QAWS_STATUS_OK && m && r.pairs == 0 && r.points == 0 && qaws_mesh_3d_get_triangle_count(m) == 12, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_refusals(void)
{
	qaws_scalar x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	unsigned int t[3] = { 0, 1, 2 }, bad[3] = { 0, 1, 3 };
	qaws_intersect_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.tolerance = -1;
	TEST_ASSERT(qaws_mesh_resolve_intersections(&m, NULL, x, t, NULL, 3, 0, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m, "no faces refused");
	TEST_ASSERT(qaws_mesh_resolve_intersections(&m, NULL, x, bad, NULL, 3, 1, NULL) == QAWS_STATUS_INVALID_ARGUMENT && !m, "index out of range refused");
	TEST_ASSERT(qaws_mesh_resolve_intersections(&m, NULL, x, t, NULL, 3, 1, &d) == QAWS_STATUS_INVALID_ARGUMENT && !m, "negative tolerance refused");
}

int test_162_intersect_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 162: self-intersections resolved\n");
	test_cubes();
	test_coplanar();
	test_spheres();
	test_alone();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
