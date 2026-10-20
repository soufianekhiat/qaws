/*
 * Test 156: boundary layers (qaws_boundary_layers.h).
 *
 * - A cube of quads, inward and outward: hexahedra only, every vertex moved by
 *   exactly the layers' height across the faces it touches (corners along
 *   their diagonal stretched by sqrt 3), so the layers fill the cube less
 *   the cube shrunk by their thickness (or the grown cube less the cube),
 *   to rounding; the border two closed surfaces.
 * - An open flat grid grown along its normals: area times thickness.
 * - An icosphere: prisms, nothing thinned, the last layer at 1 - thickness.
 * - A slab thinner than the layers asked: the stacks thinned, the two sides'
 *   layers apart, every cell valid.
 * - An L-shaped block: inward directions converging at its convex edges,
 *   the stacks there thinned, every cell valid.
 * - Refusals: no layers, no height, a shrinking ratio, a pentagon.
 */

#include "test_common.h"
#include "qaws_boundary_layers.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* a closed surface of quads on the box lo..hi, n quads along each side, turned outward (vertices welded) */
static qaws_mesh_3d* bl_box(double const* lo, double const* hi, unsigned int n)
{
	unsigned int N = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * N * N * N), nv = 0, nf = 0, a, i, j, k;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * N * N * N);
	unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * 6 * n * n);
	unsigned int* fs = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < N * N * N; i++)
		id[i] = 0xFFFFFFFFu;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			for (k = 0; k < N; k++)
				if (i == 0 || j == 0 || k == 0 || i == n || j == n || k == n)
				{
					id[(i * N + j) * N + k] = nv;
					x[3 * nv] = (qaws_scalar)(lo[0] + (hi[0] - lo[0]) * i / n);
					x[3 * nv + 1] = (qaws_scalar)(lo[1] + (hi[1] - lo[1]) * j / n);
					x[3 * nv + 2] = (qaws_scalar)(lo[2] + (hi[2] - lo[2]) * k / n);
					nv++;
				}
	/* each side: axis a fixed at 0 or n, the quad's corners counter-clockwise seen from outside */
	for (a = 0; a < 3; a++)
	{
		unsigned int b = (a + 1) % 3, c = (a + 2) % 3, side;
		for (side = 0; side < 2; side++)
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++)
				{
					unsigned int q[4][3], t;
					static unsigned int const du[4] = { 0, 1, 1, 0 }, dw[4] = { 0, 0, 1, 1 };
					for (t = 0; t < 4; t++)
					{
						unsigned int u = side ? du[t] : dw[t], w = side ? dw[t] : du[t];   /* (the order turned on the low side) */
						q[t][a] = side * n;
						q[t][b] = i + u;
						q[t][c] = j + w;
						fv[4 * nf + t] = id[(q[t][0] * N + q[t][1]) * N + q[t][2]];
					}
					fs[nf++] = 4;
				}
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = nv; d.face_vertices = fv; d.face_sizes = fs; d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(id); free(x); free(fv); free(fs);
	return m;
}

/* the volume a closed surface encloses (its faces fanned) */
static double bl_enclosed(qaws_mesh_3d const* m)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), i;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt);
	double v = 0;
	qaws_mesh_3d_get_positions(x, m);
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		qaws_scalar const *a = x + 3 * t[3 * i], *b = x + 3 * t[3 * i + 1], *c = x + 3 * t[3 * i + 2];
		v += ((double)a[0] * ((double)b[1] * c[2] - (double)b[2] * c[1]) - (double)a[1] * ((double)b[0] * c[2] - (double)b[2] * c[0])
			+ (double)a[2] * ((double)b[0] * c[1] - (double)b[1] * c[0])) / 6;
	}
	free(x); free(t);
	return v;
}

/* checks shared by the cases: valid, the cells as asked, the volume */
static double bl_check(qaws_volume_mesh const* v, unsigned int prisms, unsigned int hexes, qaws_mesh_topology* top)
{
	qaws_volume_quality q;
	qaws_mesh_3d* b = NULL;
	memset(&q, 0, sizeof q);
	memset(top, 0, sizeof *top);
	if (!v || qaws_volume_mesh_validate(v) != QAWS_STATUS_OK)
		return -1;
	qaws_volume_mesh_measure_quality(&q, v);
	if (q.prisms != prisms || q.hexahedra != hexes || q.tetrahedra || q.pyramids)
		return -2;
	if (qaws_volume_mesh_get_boundary(&b, v) == QAWS_STATUS_OK)
		qaws_mesh_3d_get_topology(top, b);
	qaws_mesh_3d_destroy(b);
	return (double)q.volume;
}

static void test_cube(void)
{
	static double const lo[3] = { -1, -1, -1 }, hi[3] = { 1, 1, 1 };
	qaws_mesh_3d* m = bl_box(lo, hi, 4);
	unsigned int nf = qaws_mesh_3d_get_face_count(m), dir;
	char msg[300];
	for (dir = 0; dir < 2; dir++)
	{
		qaws_boundary_layers_desc d;
		qaws_mesh_report rep;
		qaws_volume_mesh* v = NULL;
		qaws_mesh_topology top;
		qaws_status st;
		double T = 0.02 * (1 + 1.5 + 2.25 + 3.375), vol, want;
		memset(&d, 0, sizeof d);
		memset(&rep, 0, sizeof rep);
		d.layer_count = 4; d.first_height = (qaws_scalar)0.02; d.growth = (qaws_scalar)1.5; d.outward = (int)dir;
		st = qaws_volume_mesh_grow_layers(&v, &rep, m, &d);
		vol = bl_check(v, 0, 4 * nf, &top);
		want = dir ? pow(2 + 2 * T, 3) - 8 : 8 - pow(2 - 2 * T, 3);
		sprintf(msg, "cube %s, 4 layers from 0.02 by 1.5: %s, %u hexahedra, volume %.12f (want %.12f), border %u components closed %d, %u thinned",
			dir ? "outward" : "inward", qaws_status_to_string(st), v ? qaws_volume_mesh_get_cell_count(v) : 0, vol, want, top.components, top.closed, rep.moves);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && fabs(vol - want) < TEST_BOUND(1e-12, 64) * want && top.components == 2 && top.closed && rep.moves == 0, msg);
		qaws_volume_mesh_destroy(v);
	}
	qaws_mesh_3d_destroy(m);
}

static void test_open_grid(void)
{
	enum { N = 6 };
	qaws_scalar x[3 * (N + 1) * (N + 1)];
	unsigned int fv[4 * N * N], fs[N * N], i, j;
	qaws_mesh_3d_desc md;
	qaws_mesh_3d* m = NULL;
	qaws_boundary_layers_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_mesh_topology top;
	qaws_status st;
	double vol, want = 4 * 0.05 * 3;
	char msg[300];
	for (i = 0; i <= N; i++)
		for (j = 0; j <= N; j++)
		{
			x[3 * (i * (N + 1) + j)] = (qaws_scalar)(2.0 * i / N);
			x[3 * (i * (N + 1) + j) + 1] = (qaws_scalar)(2.0 * j / N);
			x[3 * (i * (N + 1) + j) + 2] = 0;
		}
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
		{
			unsigned int f = i * N + j;
			fv[4 * f] = i * (N + 1) + j; fv[4 * f + 1] = (i + 1) * (N + 1) + j; fv[4 * f + 2] = (i + 1) * (N + 1) + j + 1; fv[4 * f + 3] = i * (N + 1) + j + 1;
			fs[f] = 4;
		}
	memset(&md, 0, sizeof md);
	md.positions = x; md.vertex_count = (N + 1) * (N + 1); md.face_vertices = fv; md.face_sizes = fs; md.face_count = N * N;
	qaws_mesh_3d_create(&m, &md);
	memset(&d, 0, sizeof d);
	d.layer_count = 3; d.first_height = (qaws_scalar)0.05; d.growth = 1; d.outward = 1;
	st = qaws_volume_mesh_grow_layers(&v, NULL, m, &d);
	vol = bl_check(v, 0, 3 * N * N, &top);
	sprintf(msg, "open 2 x 2 grid, 3 layers of 0.05 along its normal: %s, volume %.12f (want %.12f), border closed %d", qaws_status_to_string(st), vol, want, top.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs(vol - want) < TEST_BOUND(1e-12, 64) * want && top.closed, msg);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(m);
}

/* an octahedron subdivided level times, on the unit sphere, turned outward */
static qaws_mesh_3d* bl_sphere(unsigned int level)
{
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[8][3] = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 }, { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
	unsigned int cap = 8u << (2 * level), nf = 8, nv = 6, l, i, k;
	double* x = (double*)malloc(sizeof(double) * 3 * (cap + 6));
	unsigned int *t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * cap), *t2 = (unsigned int*)malloc(sizeof(unsigned int) * 3 * cap);
	unsigned int* key = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (cap + 6));
	qaws_scalar* xs;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	memcpy(x, ov, sizeof ov);
	memcpy(t, of, sizeof of);
	for (l = 0; l < level; l++)
	{
		unsigned int n2 = 0, nk = 0;
		for (i = 0; i < nf; i++)
		{
			unsigned int mid[3];
			for (k = 0; k < 3; k++)
			{
				unsigned int a = t[3 * i + k], b = t[3 * i + (k + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, j;
				for (j = 0; j < nk && !(key[3 * j] == lo && key[3 * j + 1] == hi); j++) {}
				if (j == nk)
				{
					double* p = x + 3 * nv;
					double r;
					p[0] = x[3 * a] + x[3 * b]; p[1] = x[3 * a + 1] + x[3 * b + 1]; p[2] = x[3 * a + 2] + x[3 * b + 2];
					r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
					p[0] /= r; p[1] /= r; p[2] /= r;
					key[3 * nk] = lo; key[3 * nk + 1] = hi; key[3 * nk + 2] = nv++;
					nk++;
				}
				mid[k] = key[3 * j + 2];
			}
			{
				unsigned int a = t[3 * i], b = t[3 * i + 1], c = t[3 * i + 2];
				unsigned int const q[4][3] = { { a, mid[0], mid[2] }, { mid[0], b, mid[1] }, { mid[2], mid[1], c }, { mid[0], mid[1], mid[2] } };
				memcpy(t2 + 3 * n2, q, sizeof q);
				n2 += 4;
			}
		}
		memcpy(t, t2, sizeof(unsigned int) * 3 * n2);
		nf = n2;
	}
	xs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	for (i = 0; i < 3 * nv; i++)
		xs[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = xs; d.vertex_count = nv; d.face_vertices = t; d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(x); free(t); free(t2); free(key); free(xs);
	return m;
}

static void test_sphere(void)
{
	qaws_mesh_3d* m = bl_sphere(3);
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nv = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_boundary_layers_desc d;
	qaws_mesh_report rep;
	qaws_volume_mesh* v = NULL;
	qaws_mesh_topology top;
	qaws_status st;
	double T = 0.01 * (1 + 1.2 + 1.44 + 1.728 + 2.0736), rmin = 1e9, rmax = 0;
	char msg[300];
	memset(&d, 0, sizeof d);
	memset(&rep, 0, sizeof rep);
	d.layer_count = 5; d.first_height = (qaws_scalar)0.01;   /* (growth 0: 1.2) */
	st = qaws_volume_mesh_grow_layers(&v, &rep, m, &d);
	bl_check(v, 5 * nf, 0, &top);
	if (v)
	{
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * qaws_volume_mesh_get_vertex_count(v));
		qaws_volume_mesh_get_positions(x, v);
		for (i = 5 * nv; i < 6 * nv; i++)
		{
			double r = sqrt((double)x[3 * i] * x[3 * i] + (double)x[3 * i + 1] * x[3 * i + 1] + (double)x[3 * i + 2] * x[3 * i + 2]);
			rmin = fmin(rmin, r);
			rmax = fmax(rmax, r);
		}
		free(x);
	}
	sprintf(msg, "icosphere (%u triangles), 5 layers from 0.01 by 1.2: %s, %u prisms, the last layer at radius %.5f..%.5f (1 - %.5f), %u thinned, border closed %d in %u",
		nf, qaws_status_to_string(st), v ? qaws_volume_mesh_get_cell_count(v) : 0, rmin, rmax, T, rep.moves, top.closed, top.components);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && fabs(rmin - (1 - T)) < 0.05 * T && fabs(rmax - (1 - T)) < 0.05 * T && rep.moves == 0 && top.closed && top.components == 2, msg);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(m);
}

static void test_slab(void)
{
	static double const lo[3] = { -1, -1, -0.05 }, hi[3] = { 1, 1, 0.05 };
	qaws_mesh_3d* m = bl_box(lo, hi, 8);
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nv = qaws_mesh_3d_get_vertex_count(m), i, crossed = 0;
	qaws_boundary_layers_desc d;
	qaws_mesh_report rep;
	qaws_volume_mesh* v = NULL;
	qaws_mesh_topology top;
	qaws_status st;
	double zt = 1, zb = -1;
	char msg[300];
	memset(&d, 0, sizeof d);
	memset(&rep, 0, sizeof rep);
	d.layer_count = 3; d.first_height = (qaws_scalar)0.04; d.growth = (qaws_scalar)1.5;   /* (0.19 deep in a slab 0.1 thick) */
	st = qaws_volume_mesh_grow_layers(&v, &rep, m, &d);
	bl_check(v, 0, 3 * nf, &top);
	if (v)
	{
		/* the last layer: from the top face it stays above the bottom face's, which stays below */
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * qaws_volume_mesh_get_vertex_count(v));
		qaws_scalar* x0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		qaws_volume_mesh_get_positions(x, v);
		qaws_mesh_3d_get_positions(x0, m);
		for (i = 0; i < nv; i++)
		{
			double z = x[3 * (3 * nv + i) + 2];
			if (fabs((double)x0[3 * i] ) < 0.9 && fabs((double)x0[3 * i + 1]) < 0.9)
			{
				if (x0[3 * i + 2] > 0) zt = fmin(zt, z);
				else zb = fmax(zb, z);
			}
			crossed += fabs(z) >= 0.05;
		}
		free(x); free(x0);
	}
	sprintf(msg, "slab 0.1 thick, layers 0.19 deep asked: %s, %u hexahedra, %u stacks thinned, the top side's last layer down to z %.4f, the bottom's up to %.4f, %u outside, border closed %d",
		qaws_status_to_string(st), v ? qaws_volume_mesh_get_cell_count(v) : 0, rep.moves, zt, zb, crossed, top.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && rep.moves > 0 && zt > zb && zb > -0.05 && zt < 0.05 && crossed == 0 && top.closed, msg);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(m);
}

static void test_l_block(void)
{
	/* an L in x, y (the unit square at (0, 1) left out of the 2 x 2 square), 0.5 high, of unit quads (12 vertices a level) */
	static double const p[8][2] = { { 0, 0 }, { 2, 0 }, { 2, 1 }, { 1, 1 }, { 1, 2 }, { 0, 2 }, { 1, 0 }, { 0, 1 } };
	qaws_scalar x[3 * 16];
	/* the L's cells: (0, 0)-(1, 1): 0 6 3 7, (1, 0)-(2, 1): 6 1 2 3, (0, 1)-(1, 2): 7 3 4 5 */
	static unsigned int const cell[3][4] = { { 0, 6, 3, 7 }, { 6, 1, 2, 3 }, { 7, 3, 4, 5 } };
	static unsigned int const ring[8] = { 0, 6, 1, 2, 3, 4, 5, 7 };   /* the outline counter-clockwise */
	unsigned int fv[4 * 14], fs[14], nf = 0, i, k;
	qaws_mesh_3d_desc md;
	qaws_mesh_3d* m = NULL;
	qaws_boundary_layers_desc d;
	qaws_mesh_report rep;
	qaws_volume_mesh* v = NULL;
	qaws_mesh_topology top;
	qaws_status st;
	char msg[300];
	for (i = 0; i < 8; i++)
		for (k = 0; k < 2; k++)
		{
			x[3 * (k * 8 + i)] = (qaws_scalar)p[i][0];
			x[3 * (k * 8 + i) + 1] = (qaws_scalar)p[i][1];
			x[3 * (k * 8 + i) + 2] = (qaws_scalar)(0.5 * k);
		}
	for (i = 0; i < 3; i++)
	{
		/* bottom (seen from below: turned), top */
		fv[4 * nf] = cell[i][0]; fv[4 * nf + 1] = cell[i][3]; fv[4 * nf + 2] = cell[i][2]; fv[4 * nf + 3] = cell[i][1]; fs[nf++] = 4;
		for (k = 0; k < 4; k++) fv[4 * nf + k] = 8 + cell[i][k];
		fs[nf++] = 4;
	}
	for (i = 0; i < 8; i++)
	{
		unsigned int a = ring[i], b = ring[(i + 1) % 8];
		fv[4 * nf] = a; fv[4 * nf + 1] = b; fv[4 * nf + 2] = 8 + b; fv[4 * nf + 3] = 8 + a; fs[nf++] = 4;
	}
	memset(&md, 0, sizeof md);
	md.positions = x; md.vertex_count = 16; md.face_vertices = fv; md.face_sizes = fs; md.face_count = nf;
	qaws_mesh_3d_create(&m, &md);
	memset(&d, 0, sizeof d);
	memset(&rep, 0, sizeof rep);
	d.layer_count = 4; d.first_height = (qaws_scalar)0.04; d.growth = (qaws_scalar)1.3;   /* (0.23 deep in a block 0.5 high) */
	st = qaws_volume_mesh_grow_layers(&v, &rep, m, &d);
	bl_check(v, 0, 4 * nf, &top);
	sprintf(msg, "L-shaped block 0.5 high, layers 0.23 deep: %s, %u hexahedra, valid %d, %u stacks thinned, border closed %d in %u",
		qaws_status_to_string(st), v ? qaws_volume_mesh_get_cell_count(v) : 0, v && qaws_volume_mesh_validate(v) == QAWS_STATUS_OK, rep.moves, top.closed, top.components);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && v && qaws_volume_mesh_validate(v) == QAWS_STATUS_OK && rep.moves > 0 && top.closed && top.components == 2, msg);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(m);
}

static void test_refusals(void)
{
	static double const lo[3] = { 0, 0, 0 }, hi[3] = { 1, 1, 1 };
	qaws_mesh_3d *m = bl_box(lo, hi, 2), *pent = NULL;
	qaws_boundary_layers_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_status s0, s1, s2, s3, s4;
	qaws_scalar px[15] = { 0, 0, 0, 1, 0, 0, 1.3f, 0.8f, 0, 0.5f, 1.3f, 0, -0.3f, 0.8f, 0 };
	unsigned int pv[5] = { 0, 1, 2, 3, 4 }, ps = 5;
	qaws_mesh_3d_desc md;
	char msg[300];
	memset(&md, 0, sizeof md);
	md.positions = px; md.vertex_count = 5; md.face_vertices = pv; md.face_sizes = &ps; md.face_count = 1;
	qaws_mesh_3d_create(&pent, &md);
	memset(&d, 0, sizeof d);
	d.layer_count = 0; d.first_height = (qaws_scalar)0.1;
	s0 = qaws_volume_mesh_grow_layers(&v, NULL, m, &d);
	d.layer_count = 2; d.first_height = 0;
	s1 = qaws_volume_mesh_grow_layers(&v, NULL, m, &d);
	d.first_height = (qaws_scalar)0.1; d.growth = (qaws_scalar)0.5;
	s2 = qaws_volume_mesh_grow_layers(&v, NULL, m, &d);
	d.growth = 0;
	s3 = qaws_volume_mesh_grow_layers(&v, NULL, pent, &d);
	s4 = qaws_volume_mesh_grow_layers(NULL, NULL, m, &d);
	sprintf(msg, "refused: no layers %s, no height %s, ratio 0.5 %s, a pentagon %s, no output %s", qaws_status_to_string(s0), qaws_status_to_string(s1),
		qaws_status_to_string(s2), qaws_status_to_string(s3), qaws_status_to_string(s4));
	printf("  %s\n", msg);
	TEST_ASSERT(s0 == QAWS_STATUS_INVALID_ARGUMENT && s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT
		&& s4 == QAWS_STATUS_INVALID_ARGUMENT && v == NULL, msg);
	qaws_mesh_3d_destroy(m);
	qaws_mesh_3d_destroy(pent);
}

int test_156_boundary_layers_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 156: boundary layers\n");
	test_cube();
	test_open_grid();
	test_sphere();
	test_slab();
	test_l_block();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
