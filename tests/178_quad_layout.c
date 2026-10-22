/*
 * Test 178: quad layouts (QAWS_QUADS_LAYOUT), small meshes
 *
 *   - a box of triangles, its edges creases: GLOBAL's quads parted by their
 *     motorcycle graph into the box's 6 faces, each a grid;
 *   - figure 2's blob: patches between the separatrices of its irregular
 *     vertices, each a rectangle of quads (4 corners, none concave), every
 *     face in one, the attribute "patch" numbering them;
 *   - a torus without irregular vertices: one periodic patch round it.
 */

#include "test_common.h"
#include "qaws_quadrangulate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_mesh_3d* ly_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int i;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(p);
	return m;
}

/* the showcase's blob: an icosahedron split once (edges split, then new-old edges flipped), pushed in and out smoothly */
static qaws_mesh_3d* ly_blob(void)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[36];
	double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	qaws_mesh_3d* m;
	unsigned int i, n, ne, nv0, *e, *src, w;
	qaws_scalar* p;
	for (i = 0; i < 36; i++)
		x[i] = v[i] / sqrt(1 + g * g);
	m = ly_make(x, 12, f, 20);
	ne = qaws_mesh_3d_get_edge_count(m);
	nv0 = qaws_mesh_3d_get_vertex_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	src = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
		if (qaws_mesh_3d_split_edge(&w, m, NULL, NULL, e[2 * i], e[2 * i + 1]) == QAWS_STATUS_OK)
		{
			src[2 * (w - nv0)] = e[2 * i];
			src[2 * (w - nv0) + 1] = e[2 * i + 1];
		}
	free(e);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		unsigned int a = e[2 * i], b = e[2 * i + 1], nw, old;
		if ((a >= nv0) == (b >= nv0))
			continue;
		nw = a >= nv0 ? a : b;
		old = a >= nv0 ? b : a;
		if (old != src[2 * (nw - nv0)] && old != src[2 * (nw - nv0) + 1] && qaws_mesh_3d_can_flip_edge(m, a, b))
			qaws_mesh_3d_flip_edge(m, a, b);
	}
	free(e);
	free(src);
	n = qaws_mesh_3d_get_vertex_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < n; i++)
	{
		double q[3], l, s;
		q[0] = (double)p[3 * i]; q[1] = (double)p[3 * i + 1]; q[2] = (double)p[3 * i + 2];
		l = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		q[0] /= l; q[1] /= l; q[2] /= l;
		s = 1 + 0.32 * sin(2.2 * q[0] + 0.4) * cos(1.7 * q[1]) + 0.2 * q[2] * q[2];
		p[3 * i] = (qaws_scalar)(q[0] * 1.25 * s); p[3 * i + 1] = (qaws_scalar)(q[1] * s); p[3 * i + 2] = (qaws_scalar)(q[2] * 0.9 * s);
	}
	qaws_mesh_3d_set_positions(m, p);
	free(p);
	return m;
}
/* a box of sides 1, 0.8, 0.6, each face an n x n grid of triangle pairs, its edges creases */
static qaws_mesh_3d* ly_box(unsigned int n)
{
	unsigned int N = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * N * N * N), nv = 0, nf = 0, a, s, u, v, k;
	double* x = (double*)malloc(sizeof(double) * 3 * N * N * N);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
	static double const size[3] = { 1, 0.8, 0.6 };
	qaws_mesh_3d* m;
	for (k = 0; k < N * N * N; k++)
	{
		unsigned int p[3];
		p[0] = k % N;
		p[1] = (k / N) % N;
		p[2] = k / (N * N);
		id[k] = 0xFFFFFFFFu;
		if (p[0] % n && p[1] % n && p[2] % n)
			continue;
		id[k] = nv;
		for (a = 0; a < 3; a++)
			x[3 * nv + a] = size[a] * ((double)p[a] / n - 0.5);
		nv++;
	}
	for (a = 0; a < 3; a++)
		for (s = 0; s < 2; s++)
			for (u = 0; u < n; u++)
				for (v = 0; v < n; v++)
				{
					unsigned int q[4], c, b = (a + 1) % 3, d = (a + 2) % 3;
					for (c = 0; c < 4; c++)
					{
						unsigned int p[3];
						p[a] = s * n;
						p[b] = u + (c == 1 || c == 2);
						p[d] = v + (c >= 2);
						q[c] = id[p[0] + N * (p[1] + N * p[2])];
					}
					if (!s)
					{
						c = q[1];
						q[1] = q[3];
						q[3] = c;
					}
					f[3 * nf] = q[0]; f[3 * nf + 1] = q[1]; f[3 * nf + 2] = q[2]; nf++;
					f[3 * nf] = q[0]; f[3 * nf + 1] = q[2]; f[3 * nf + 2] = q[3]; nf++;
				}
	m = ly_make(x, nv, f, nf);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(50 * M_PI / 180), 0, 0);
	free(id);
	free(x);
	free(f);
	return m;
}

/* test 125's torus: nu x nv quads of radii R and r split in triangles, every other column shifted half a step */
static qaws_mesh_3d* ly_torus(unsigned int nu, unsigned int nv, double R, double r)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * M_PI * i / nu, v = 2 * M_PI * (j + 0.5 * (i & 1)) / nv;
			x[3 * (j * nu + i)] = (R + r * cos(v)) * cos(u);
			x[3 * (j * nu + i) + 1] = (R + r * cos(v)) * sin(u);
			x[3 * (j * nu + i) + 2] = r * sin(v);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, d = ((j + 1) % nv) * nu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	m = ly_make(x, nu * nv, f, 2 * nu * nv);
	free(x);
	free(f);
	return m;
}

typedef struct ly_res
{
	unsigned int patches, faces, rect, nonquad;
	int st;
} ly_res;

/* LAYOUT at L: its patches, and how many are rectangles of quads (4 corners: vertices of the patch's border in one of its
   faces only; none in three) */
static ly_res ly_run(qaws_mesh_3d const* m, double L)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* q = NULL;
	qaws_mesh_report rep;
	ly_res o;
	memset(&o, 0, sizeof o);
	memset(&d, 0, sizeof d);
	d.method = QAWS_QUADS_LAYOUT;
	d.target_length = (qaws_scalar)L;
	o.st = (int)qaws_mesh_quadrangulate(&q, &rep, m, &d);
	if (q)
	{
		unsigned int id, nf = qaws_mesh_3d_get_face_count(q), nv = qaws_mesh_3d_get_vertex_count(q), nc = qaws_mesh_3d_get_corner_count(q), f, c, p;
		qaws_scalar* pa = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
		unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
		unsigned int* cnt = (unsigned int*)malloc(sizeof(unsigned int) * nv);
		o.faces = nf;
		o.patches = rep.splits;
		qaws_mesh_3d_find_attribute(&id, q, "patch");
		qaws_mesh_3d_get_attribute(pa, q, id);
		qaws_mesh_3d_get_faces(fv, fs, q);
		for (f = 0; f < nf; f++) o.nonquad += fs[f] != 4;
		for (p = 0; p < o.patches; p++)
		{
			unsigned int corners = 0, concave = 0, v, total = 0;
			memset(cnt, 0, sizeof(unsigned int) * nv);
			for (f = 0, c = 0; f < nf; c += fs[f++])
				if ((unsigned int)pa[f] == p)
				{
					unsigned int k;
					for (k = 0; k < fs[f]; k++) cnt[fv[c + k]]++;
					total++;
				}
			/* a vertex inside the patch has 4 of its faces, on its side 2, at a corner 1 (3: a concave corner) */
			for (v = 0; v < nv; v++)
			{
				corners += cnt[v] == 1;
				concave += cnt[v] == 3;
			}
			o.rect += total > 0 && corners == 4 && concave == 0;
		}
		free(pa); free(fv); free(fs); free(cnt);
		qaws_mesh_3d_destroy(q);
	}
	return o;
}

static void test_layouts(void)
{
	qaws_mesh_3d *box = ly_box(4), *blob = ly_blob(), *torus = ly_torus(24, 12, 1, 0.4);
	ly_res a = ly_run(box, 0.2), b = ly_run(blob, 0.35), t = ly_run(torus, 0.25);
	char msg[300];
	sprintf(msg, "a box: %u faces in %u patches, %u rectangles, %u not quads", a.faces, a.patches, a.rect, a.nonquad);
	printf("  %s\n", msg);
	TEST_ASSERT(a.st == QAWS_STATUS_OK && a.patches == 6 && a.rect == 6 && a.nonquad == 0, msg);
	sprintf(msg, "the blob: %u faces in %u patches, %u rectangles, %u not quads", b.faces, b.patches, b.rect, b.nonquad);
	printf("  %s\n", msg);
	TEST_ASSERT(b.st == QAWS_STATUS_OK && b.patches > 6 && b.rect == b.patches && b.nonquad == 0, msg);
	sprintf(msg, "a torus without irregular vertices: %u faces in %u patch, all quads %d", t.faces, t.patches, t.nonquad == 0);
	printf("  %s\n", msg);
	TEST_ASSERT(t.st == QAWS_STATUS_OK && t.patches == 1 && t.nonquad == 0, msg);
	qaws_mesh_3d_destroy(box);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(torus);
}

int test_178_quad_layout_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 178: quad layouts\n");
	test_layouts();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
