/*
 * Test 173: improvement of tetrahedral meshes (qaws_volume_mesh_improve),
 * small meshes
 *
 *   - three flat wedges round a long edge: removed (the 3-2 flip),
 *     two left, the border's six triangles kept;
 *   - a cap (an apex just above a triangle's middle) on a tetrahedron below:
 *     the triangle removed (the 2-3 flip), three round the edge between the
 *     apexes;
 *   - a cube's Kuhn grid (4^3 cubes, 384 tetrahedra), its inside vertices
 *     shaken into slivers: the worst dihedral angle raised past 20 degrees,
 *     the border's triangles and vertices the same, the volume exactly 1,
 *     the result valid; smoothing alone and flips alone each raise it;
 *   - refusals: a mesh without tetrahedra.
 */

#include "test_common.h"
#include "qaws_volume_improve.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_volume_mesh* vt_make(double const* x, unsigned int nv, unsigned int const* t, unsigned int nt)
{
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int i;
	for (i = 0; i < 3 * nv; i++) p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.cell_vertices = t;
	d.cell_count = nt;
	qaws_volume_mesh_create(&m, &d);
	free(p);
	return m;
}

/* the border's triangles, each sorted, the list sorted: a fingerprint */
static int vt_cmp3(void const* a, void const* b)
{
	unsigned int const *x = (unsigned int const*)a, *y = (unsigned int const*)b;
	int i;
	for (i = 0; i < 3; i++)
		if (x[i] != y[i])
			return x[i] < y[i] ? -1 : 1;
	return 0;
}

static unsigned int vt_border(unsigned int* out, qaws_volume_mesh const* m)
{
	static unsigned char const f[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	unsigned int nc = qaws_volume_mesh_get_cell_count(m), *cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc),
		*adj = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc), n = 0, c, k;
	qaws_volume_mesh_get_cells(cv, NULL, m);
	qaws_volume_mesh_get_cell_adjacency(adj, m);
	for (c = 0; c < nc; c++)
		for (k = 0; k < 4; k++)
			if (adj[4 * c + k] == QAWS_MESH_NONE)
			{
				unsigned int* s = out + 3 * n++, t;
				s[0] = cv[4 * c + f[k][0]];
				s[1] = cv[4 * c + f[k][1]];
				s[2] = cv[4 * c + f[k][2]];
				if (s[0] > s[1]) { t = s[0]; s[0] = s[1]; s[1] = t; }
				if (s[1] > s[2]) { t = s[1]; s[1] = s[2]; s[2] = t; }
				if (s[0] > s[1]) { t = s[0]; s[0] = s[1]; s[1] = t; }
			}
	qsort(out, n, sizeof(unsigned int) * 3, vt_cmp3);
	free(cv);
	free(adj);
	return n;
}

static void test_flips(void)
{
	/* a long edge (0, 0, -1) to (0, 0, 1) in a ring of three vertices close to it: flat wedges */
	double x[15] = { 0, 0, -1, 0, 0, 1, 0.2, 0, 0, -0.1, 0.1732, 0, -0.1, -0.1732, 0 };
	unsigned int t[12] = { 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 4, 2 }, i, b0[36], b1[36], n0, n1;
	qaws_volume_mesh *m, *o = NULL;
	qaws_volume_improve_report r;
	qaws_status st;
	char msg[300];
	for (i = 0; i < 3; i++)
	{
		/* each ring cell positive */
		double const *a = x + 3 * t[4 * i], *b = x + 3 * t[4 * i + 1], *c = x + 3 * t[4 * i + 2], *d = x + 3 * t[4 * i + 3];
		double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, v[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] }, w[3] = { d[0] - a[0], d[1] - a[1], d[2] - a[2] };
		if (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]) < 0)
		{
			unsigned int s = t[4 * i + 2];
			t[4 * i + 2] = t[4 * i + 3];
			t[4 * i + 3] = s;
		}
	}
	m = vt_make(x, 5, t, 3);
	n0 = m ? vt_border(b0, m) : 0;
	st = qaws_volume_mesh_improve(&o, &r, m, NULL);
	n1 = o ? vt_border(b1, o) : 0;
	sprintf(msg, "three flat wedges round a long edge: %u left, %u edge removal(s), quality %.3f -> %.3f, border %u triangles the same %d",
		o ? qaws_volume_mesh_get_cell_count(o) : 0, r.edge_removals, (double)r.min_quality_before, (double)r.min_quality_after, n1,
		n0 == n1 && !memcmp(b0, b1, sizeof(unsigned int) * 3 * n0));
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && m && o && qaws_volume_mesh_get_cell_count(o) == 2 && r.edge_removals == 1 && n0 == 6 && n1 == 6 &&
		!memcmp(b0, b1, sizeof(unsigned int) * 18) && r.min_quality_after > r.min_quality_before && qaws_volume_mesh_validate(o) == QAWS_STATUS_OK, msg);
	qaws_volume_mesh_destroy(m);
	qaws_volume_mesh_destroy(o);
	o = NULL;
	{
		/* a cap: apexes (0, 0, 0.1) and (0, 0, -1) on a unit triangle */
		double y[15] = { 0, 0, -1, 0, 0, 0.1, 1, 0, 0, -0.5, 0.866, 0, -0.5, -0.866, 0 };
		unsigned int tt[8] = { 2, 3, 4, 1, 2, 4, 3, 0 };
		m = vt_make(y, 5, tt, 2);
		st = qaws_volume_mesh_improve(&o, &r, m, NULL);
		sprintf(msg, "a cap on a tetrahedron: %u after, %u face removal(s), quality %.3f -> %.3f", o ? qaws_volume_mesh_get_cell_count(o) : 0,
			r.face_removals, (double)r.min_quality_before, (double)r.min_quality_after);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && m && o && qaws_volume_mesh_get_cell_count(o) == 3 && r.face_removals == 1 &&
			r.min_quality_after > r.min_quality_before && qaws_volume_mesh_validate(o) == QAWS_STATUS_OK, msg);
		qaws_volume_mesh_destroy(m);
		qaws_volume_mesh_destroy(o);
	}
}

/* the unit cube's Kuhn grid of n^3 cubes (6 tetrahedra each), the inside vertices shaken by up to amp cells */
static qaws_volume_mesh* vt_kuhn(unsigned int n, double amp)
{
	static unsigned int const perm[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
	unsigned int m1 = n + 1, nv = m1 * m1 * m1, nt = 6 * n * n * n, i, j, k, p, q = 0, seed = 12345u;
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nt);
	qaws_volume_mesh* m;
	for (k = 0; k < m1; k++)
		for (j = 0; j < m1; j++)
			for (i = 0; i < m1; i++)
			{
				double* p3 = x + 3 * (i + m1 * (j + m1 * k));
				int inside = i > 0 && j > 0 && k > 0 && i < n && j < n && k < n;
				unsigned int d;
				p3[0] = (double)i / n;
				p3[1] = (double)j / n;
				p3[2] = (double)k / n;
				for (d = 0; inside && d < 3; d++)
				{
					seed = seed * 1664525u + 1013904223u;
					p3[d] += amp / n * (((seed >> 8) & 0xFFFF) / 32768.0 - 1);
				}
			}
	for (k = 0; k < n; k++)
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
				for (p = 0; p < 6; p++)
				{
					unsigned int c[3] = { i, j, k }, v[4], s, d;
					double u[3], w[3], z[3];
					v[0] = i + m1 * (j + m1 * k);
					for (s = 0; s < 3; s++)
					{
						c[perm[p][s]]++;
						v[s + 1] = c[0] + m1 * (c[1] + m1 * c[2]);
					}
					/* (the orientation of the grid before the shaking: Kuhn cells by their permutation) */
					for (d = 0; d < 3; d++)
					{
						u[d] = (double)(d == perm[p][0]);
						w[d] = u[d] + (double)(d == perm[p][1]);
						z[d] = 1;
					}
					if (u[0] * (w[1] * z[2] - w[2] * z[1]) - u[1] * (w[0] * z[2] - w[2] * z[0]) + u[2] * (w[0] * z[1] - w[1] * z[0]) < 0)
					{
						s = v[1]; v[1] = v[2]; v[2] = s;
					}
					memcpy(t + 4 * q++, v, sizeof v);
				}
	m = vt_make(x, nv, t, nt);
	free(x);
	free(t);
	return m;
}

static void test_grid(void)
{
	qaws_volume_mesh *m = vt_kuhn(4, 0.3), *o = NULL, *s = NULL, *f = NULL;
	qaws_volume_improve_report r, rs, rf;
	qaws_volume_improve_desc d;
	qaws_volume_quality q0, q1;
	unsigned int *b0 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 400), *b1 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * 400), n0 = 0, n1 = 0;
	qaws_status st;
	double moved = 0;
	char msg[400];
	memset(&q0, 0, sizeof q0);
	memset(&rs, 0, sizeof rs);
	memset(&rf, 0, sizeof rf);
	memset(&r, 0, sizeof r);
	memset(&q1, 0, sizeof q1);
	if (m)
	{
		qaws_volume_mesh_measure_quality(&q0, m);
		n0 = vt_border(b0, m);
	}
	st = qaws_volume_mesh_improve(&o, &r, m, NULL);
	if (o)
	{
		unsigned int nv = qaws_volume_mesh_get_vertex_count(o), i, k;
		qaws_scalar *p0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *p1 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_volume_mesh_measure_quality(&q1, o);
		n1 = vt_border(b1, o);
		qaws_volume_mesh_get_positions(p0, m);
		qaws_volume_mesh_get_positions(p1, o);
		for (i = 0; i < nv; i++)
		{
			int border = 0;
			for (k = 0; k < 3; k++) border |= p0[3 * i + k] == 0 || p0[3 * i + k] == 1;
			for (k = 0; border && k < 3; k++) moved = fmax(moved, fabs((double)p1[3 * i + k] - (double)p0[3 * i + k]));
		}
		free(p0);
		free(p1);
	}
	sprintf(msg, "a shaken Kuhn grid (%u tetrahedra): worst dihedral %.2f -> %.2f degrees (largest %.1f -> %.1f), %u edge and %u face removals, %u moves in %u passes; volume %.15f; border %u triangles the same %d, its vertices moved %.1e",
		q0.tetrahedra, (double)q0.min_dihedral, (double)q1.min_dihedral, (double)q0.max_dihedral, (double)q1.max_dihedral, r.edge_removals,
		r.face_removals, r.moves, r.iterations, (double)q1.volume, n1, n0 == n1 && !memcmp(b0, b1, sizeof(unsigned int) * 3 * n0), moved);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && r.min_quality_before > 0 && q0.min_dihedral < 10 && q1.min_dihedral > 20 && q1.max_dihedral < q0.max_dihedral && fabs((double)q1.volume - 1) < 1e-12 &&
		n0 == 192 && n1 == 192 && !memcmp(b0, b1, sizeof(unsigned int) * 3 * n0) && moved == 0 && qaws_volume_mesh_validate(o) == QAWS_STATUS_OK, msg);
	memset(&d, 0, sizeof d);
	d.no_flips = 1;
	qaws_volume_mesh_improve(&s, &rs, m, &d);
	d.no_flips = 0;
	d.no_smoothing = 1;
	qaws_volume_mesh_improve(&f, &rf, m, &d);
	sprintf(msg, "smoothing alone: quality %.3f -> %.3f (%u moves, no flips %d); flips alone: -> %.3f (%u edge, %u face removals, no moves %d); both: %.3f",
		(double)rs.min_quality_before, (double)rs.min_quality_after, rs.moves, rs.edge_removals + rs.face_removals == 0, (double)rf.min_quality_after,
		rf.edge_removals, rf.face_removals, rf.moves == 0, (double)r.min_quality_after);
	printf("  %s\n", msg);
	TEST_ASSERT(s && f && rs.min_quality_after > rs.min_quality_before && rs.edge_removals + rs.face_removals == 0 && rf.min_quality_after > rf.min_quality_before &&
		rf.moves == 0 && r.min_quality_after >= rs.min_quality_after && r.min_quality_after >= rf.min_quality_after, msg);
	qaws_volume_mesh_destroy(s);
	qaws_volume_mesh_destroy(f);
	qaws_volume_mesh_destroy(o);
	qaws_volume_mesh_destroy(m);
	free(b0);
	free(b1);
}

static void test_refusals(void)
{
	static double const x[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
	static unsigned int const h[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	qaws_cell_type ty = QAWS_CELL_HEXAHEDRON;
	qaws_volume_mesh_desc d;
	qaws_volume_mesh *m = NULL, *o = NULL;
	qaws_scalar p[24];
	qaws_status st;
	char msg[200];
	unsigned int i;
	for (i = 0; i < 24; i++) p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = 8;
	d.cell_vertices = h;
	d.cell_types = &ty;
	d.cell_count = 1;
	qaws_volume_mesh_create(&m, &d);
	st = qaws_volume_mesh_improve(&o, NULL, m, NULL);
	sprintf(msg, "refused: a mesh without tetrahedra (%d)", (int)st);
	printf("  %s\n", msg);
	TEST_ASSERT(m && st == QAWS_STATUS_INVALID_ARGUMENT && !o, msg);
	qaws_volume_mesh_destroy(m);
}

int test_173_volume_improve_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 173: improvement of tetrahedral meshes\n");
	test_flips();
	test_grid();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
