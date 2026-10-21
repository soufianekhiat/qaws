/*
 * Test 166: deformation and parameterization (qaws_mesh_3d_deform,
 * qaws_mesh_3d_parameterize), on small meshes
 *
 *   - a 12 x 12 plate, its left side held and its right side raised: every
 *     method puts the handles exactly and leaves the held side; ARAP keeps
 *     the edge lengths better than the biharmonic bend it starts from;
 *   - both sides turned by 60 degrees as one: ARAP turns the whole plate
 *     rigidly, Laplacian editing does not; a vertex outside the region stays;
 *   - a flat grid flattened: no distortion at all (LSCM);
 *   - a cap of an icosphere: no triangle turned over; ARAP keeps the areas
 *     better than LSCM;
 *   - a closed icosphere: refused without seams, one chart when cut from pole
 *     to pole; two separate grids: two charts side by side;
 *   - refusals: a seam that is not an edge, BFF (reserved).
 */

#include "test_common.h"
#include "qaws_mesh_deform.h"
#include "qaws_mesh_param.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* count grids of n x n cells, side 1, 2 apart along x */
static qaws_mesh_3d* dt_grid(unsigned int n, unsigned int count)
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
				t[k++] = v; t[k++] = v + 1; t[k++] = v + n + 2;
				t[k++] = v; t[k++] = v + n + 2; t[k++] = v + n + 1;
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

/* an icosahedron split in four twice (162 vertices), unit radius; keep: NULL, or the triangles kept (their centroid's z > 0) */
static qaws_mesh_3d* dt_icosphere(int cap)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2, x[3 * 162];
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int t[3 * 320], t2[3 * 320], nt = 20, nv = 12, i, j, lv, kept = 0, used[162], nu = 0;
	unsigned int mid[162 * 6][3], nm = 0;
	qaws_scalar xs[3 * 162];
	unsigned int ft[3 * 320];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 12; i++)
	{
		double l = sqrt(v0[3 * i] * v0[3 * i] + v0[3 * i + 1] * v0[3 * i + 1] + v0[3 * i + 2] * v0[3 * i + 2]);
		for (j = 0; j < 3; j++)
			x[3 * i + j] = v0[3 * i + j] / l;
	}
	for (i = 0; i < 60; i++)
		t[i] = (unsigned int)ic[i];
	for (lv = 0; lv < 2; lv++)
	{
		nm = 0;
		for (i = 0; i < nt; i++)
		{
			unsigned int c[3], e;
			for (e = 0; e < 3; e++)
			{
				unsigned int a = t[3 * i + e], b = t[3 * i + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, s;
				double l;
				for (s = 0; s < nm && !(mid[s][0] == lo && mid[s][1] == hi); s++)
					;
				if (s == nm)
				{
					mid[nm][0] = lo;
					mid[nm][1] = hi;
					mid[nm++][2] = nv;
					for (j = 0; j < 3; j++)
						x[3 * nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
					l = sqrt(x[3 * nv] * x[3 * nv] + x[3 * nv + 1] * x[3 * nv + 1] + x[3 * nv + 2] * x[3 * nv + 2]);
					for (j = 0; j < 3; j++)
						x[3 * nv + j] /= l;
					nv++;
				}
				c[e] = mid[s][2];
			}
			{
				unsigned int const q[12] = { t[3 * i], c[0], c[2], c[0], t[3 * i + 1], c[1], c[2], c[1], t[3 * i + 2], c[0], c[1], c[2] };
				memcpy(t2 + 12 * i, q, sizeof q);
			}
		}
		nt *= 4;
		memcpy(t, t2, sizeof(unsigned int) * 3 * nt);
	}
	for (i = 0; i < nv; i++)
		used[i] = cap ? 0xFFFFFFFFu : i;
	for (i = 0; i < nt; i++)
	{
		double cz = (x[3 * t[3 * i] + 2] + x[3 * t[3 * i + 1] + 2] + x[3 * t[3 * i + 2] + 2]) / 3;
		if (cap && cz <= 0.1)
			continue;
		for (j = 0; j < 3; j++)
		{
			if (used[t[3 * i + j]] == 0xFFFFFFFFu)
				used[t[3 * i + j]] = nu++;
			ft[3 * kept + j] = used[t[3 * i + j]];
		}
		kept++;
	}
	if (!cap)
		nu = nv;
	for (i = 0; i < nv; i++)
		if (used[i] != 0xFFFFFFFFu)
			for (j = 0; j < 3; j++)
				xs[3 * used[i] + j] = (qaws_scalar)x[3 * i + j];
	memset(&d, 0, sizeof d);
	d.positions = xs;
	d.vertex_count = nu;
	d.face_vertices = ft;
	d.face_count = kept;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static double dt_len_error(qaws_mesh_3d const* m, qaws_scalar const* x0)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), i, k;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	double s = 0;
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		qaws_vec3 a = qaws_mesh_3d_get_vertex(m, e[2 * i]), b = qaws_mesh_3d_get_vertex(m, e[2 * i + 1]);
		double l = sqrt((double)(a.x - b.x) * (a.x - b.x) + (double)(a.y - b.y) * (a.y - b.y) + (double)(a.z - b.z) * (a.z - b.z)), l0 = 0;
		for (k = 0; k < 3; k++)
			l0 += (x0[3 * e[2 * i] + k] - x0[3 * e[2 * i + 1] + k]) * (x0[3 * e[2 * i] + k] - x0[3 * e[2 * i + 1] + k]);
		s += fabs(l / sqrt(l0) - 1);
	}
	free(e);
	return s / ne;
}

static void test_bend(void)
{
	static char const* const name[3] = { "ARAP", "Laplacian", "biharmonic" };
	unsigned int n = 12, nv = (n + 1) * (n + 1), j, k, h[26];
	qaws_scalar tg[3 * 26], x0[3 * 169];
	double err[3];
	char msg[200];
	for (k = 0; k < 3; k++)
	{
		qaws_mesh_3d* m = dt_grid(n, 1);
		qaws_deform_desc d;
		qaws_status st;
		double off = 0;
		qaws_mesh_3d_get_positions(x0, m);
		for (j = 0; j <= n; j++)
		{
			h[j] = j * (n + 1);
			h[n + 1 + j] = j * (n + 1) + n;
			memcpy(tg + 3 * j, x0 + 3 * h[j], 3 * sizeof(qaws_scalar));
			memcpy(tg + 3 * (n + 1 + j), x0 + 3 * h[n + 1 + j], 3 * sizeof(qaws_scalar));
			tg[3 * (n + 1 + j) + 2] += (qaws_scalar)0.4;
		}
		memset(&d, 0, sizeof d);
		d.method = (qaws_deform_method)k;
		st = qaws_mesh_3d_deform(m, h, tg, 2 * (n + 1), &d);
		for (j = 0; j < 2 * (n + 1); j++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, h[j]);
			off = fmax(off, fabs(p.x - tg[3 * j]) + fabs(p.y - tg[3 * j + 1]) + fabs(p.z - tg[3 * j + 2]));
		}
		err[k] = dt_len_error(m, x0);
		sprintf(msg, "plate raised on one side, %s: handles off by %.1e, edge lengths changed by %.4f on average", name[k], off, err[k]);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && off == 0, msg);
		qaws_mesh_3d_destroy(m);
		(void)nv;
	}
	TEST_ASSERT(err[0] < err[2], "ARAP keeps the lengths better than the biharmonic bend");
}

static void test_turn(void)
{
	unsigned int n = 12, nv = (n + 1) * (n + 1), i, j, k, h[26];
	double c = cos(1.0471975511965976), s = sin(1.0471975511965976), err[2];
	qaws_scalar tg[3 * 26], x0[3 * 169];
	char msg[200];
	for (k = 0; k < 2; k++)
	{
		qaws_mesh_3d* m = dt_grid(n, 1);
		qaws_deform_desc d;
		qaws_status st;
		double e = 0;
		qaws_mesh_3d_get_positions(x0, m);
		for (j = 0; j <= n; j++)
		{
			h[j] = j * (n + 1);
			h[n + 1 + j] = j * (n + 1) + n;
		}
		for (j = 0; j < 2 * (n + 1); j++)
		{
			double px = x0[3 * h[j]] - 0.5, py = x0[3 * h[j] + 1] - 0.5;
			tg[3 * j] = (qaws_scalar)(0.5 + c * px - s * py);
			tg[3 * j + 1] = (qaws_scalar)(0.5 + s * px + c * py);
			tg[3 * j + 2] = 0;
		}
		memset(&d, 0, sizeof d);
		d.method = k ? QAWS_DEFORM_LAPLACIAN : QAWS_DEFORM_ARAP;
		st = qaws_mesh_3d_deform(m, h, tg, 2 * (n + 1), &d);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
			double px = x0[3 * i] - 0.5, py = x0[3 * i + 1] - 0.5;
			e = fmax(e, fabs(p.x - (0.5 + c * px - s * py)) + fabs(p.y - (0.5 + s * px + c * py)) + fabs((double)p.z));
		}
		err[k] = e;
		TEST_ASSERT(st == QAWS_STATUS_OK, "plate turned");
		qaws_mesh_3d_destroy(m);
	}
	sprintf(msg, "plate turned 60 degrees by its two sides: off the turned plate by %.2e (ARAP), %.2e (Laplacian)", err[0], err[1]);
	printf("  %s\n", msg);
	TEST_ASSERT(err[0] < 1e-4 && err[1] > 0.1, msg);
	{
		/* a region: the vertices outside it stay */
		qaws_mesh_3d* m = dt_grid(n, 1);
		unsigned char* region = (unsigned char*)calloc(nv, 1);
		unsigned int hv = (n / 2) * (n + 1) + n / 2, fixed = 0;
		qaws_scalar t1[3];
		qaws_deform_desc d;
		qaws_vec3 p;
		for (i = 0; i < nv; i++)
			region[i] = (i % (n + 1)) >= 3 && (i % (n + 1)) <= n - 3 && (i / (n + 1)) >= 3 && (i / (n + 1)) <= n - 3;
		p = qaws_mesh_3d_get_vertex(m, hv);
		t1[0] = p.x;
		t1[1] = p.y;
		t1[2] = (qaws_scalar)0.2;
		memset(&d, 0, sizeof d);
		d.region = region;
		qaws_mesh_3d_deform(m, &hv, t1, 1, &d);
		p = qaws_mesh_3d_get_vertex(m, 0);
		fixed = p.x == 0 && p.y == 0 && p.z == 0;
		p = qaws_mesh_3d_get_vertex(m, (n / 2) * (n + 1) + n / 2 + 1);
		TEST_ASSERT(fixed && p.z > 0.05, "a bump in a region: outside it nothing moves, inside it follows");
		qaws_mesh_3d_destroy(m);
		free(region);
	}
}

static void test_param(void)
{
	char msg[240];
	{
		qaws_mesh_3d* m = dt_grid(8, 1);
		qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * qaws_mesh_3d_get_triangle_count(m));
		qaws_param_report r;
		qaws_status st = qaws_mesh_3d_parameterize(uv, &r, m, NULL, 0, NULL);
		sprintf(msg, "flat grid, LSCM: %u chart, %u turned, angle distortion %.2e, area distortion %.2e", r.charts, r.flipped,
			(double)r.angle_distortion, (double)r.area_distortion);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && r.charts == 1 && r.flipped == 0 && r.angle_distortion < 1e-9 && r.area_distortion < 1e-9, msg);
		qaws_mesh_3d_destroy(m);
		free(uv);
	}
	{
		qaws_mesh_3d* m = dt_icosphere(1);
		qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * qaws_mesh_3d_get_triangle_count(m));
		qaws_param_report r[2];
		qaws_param_desc d;
		unsigned int k;
		qaws_status st[2];
		for (k = 0; k < 2; k++)
		{
			memset(&d, 0, sizeof d);
			d.method = k ? QAWS_PARAM_ARAP : QAWS_PARAM_LSCM;
			st[k] = qaws_mesh_3d_parameterize(uv, &r[k], m, NULL, 0, &d);
		}
		sprintf(msg, "a sphere's cap (%u triangles): LSCM %u turned, angles %.4f, areas %.4f; ARAP %u turned, angles %.4f, areas %.4f",
			qaws_mesh_3d_get_triangle_count(m), r[0].flipped, (double)r[0].angle_distortion, (double)r[0].area_distortion, r[1].flipped,
			(double)r[1].angle_distortion, (double)r[1].area_distortion);
		printf("  %s\n", msg);
		TEST_ASSERT(st[0] == QAWS_STATUS_OK && st[1] == QAWS_STATUS_OK && r[0].flipped == 0 && r[1].flipped == 0
			&& r[1].area_distortion < r[0].area_distortion && r[0].angle_distortion < r[1].angle_distortion, msg);
		qaws_mesh_3d_destroy(m);
		free(uv);
	}
	{
		/* the closed sphere: refused, then cut from its top vertex down to its bottom one along the lowest neighbours */
		qaws_mesh_3d* m = dt_icosphere(0);
		unsigned int nv = qaws_mesh_3d_get_vertex_count(m), top = 0, bot = 0, v, seam[64], ns = 0, ring[16], cnt, i;
		qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * qaws_mesh_3d_get_triangle_count(m));
		qaws_param_report r;
		qaws_status st;
		TEST_ASSERT(qaws_mesh_3d_parameterize(uv, &r, m, NULL, 0, NULL) == QAWS_STATUS_INVALID_TOPOLOGY, "a closed sphere without seams refused");
		for (v = 0; v < nv; v++)
		{
			if (qaws_mesh_3d_get_vertex(m, v).z > qaws_mesh_3d_get_vertex(m, top).z)
				top = v;
			if (qaws_mesh_3d_get_vertex(m, v).z < qaws_mesh_3d_get_vertex(m, bot).z)
				bot = v;
		}
		for (v = top; v != bot && ns < 32;)
		{
			unsigned int best = v;
			qaws_mesh_3d_get_vertex_ring(ring, &cnt, m, 16, v);
			for (i = 0; i < cnt; i++)
				if (qaws_mesh_3d_get_vertex(m, ring[i]).z < qaws_mesh_3d_get_vertex(m, best).z)
					best = ring[i];
			if (best == v)
				break;
			seam[2 * ns] = v;
			seam[2 * ns + 1] = best;
			ns++;
			v = best;
		}
		st = qaws_mesh_3d_parameterize(uv, &r, m, seam, ns, NULL);
		sprintf(msg, "the sphere cut from pole to pole (%u seam edges), LSCM: %u chart, %u turned", ns, r.charts, r.flipped);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && r.charts == 1 && r.flipped == 0, msg);
		seam[0] = top;
		seam[1] = bot;
		TEST_ASSERT(qaws_mesh_3d_parameterize(uv, &r, m, seam, 1, NULL) == QAWS_STATUS_INVALID_ARGUMENT, "a seam that is not an edge refused");
		{
			qaws_param_desc d;
			memset(&d, 0, sizeof d);
			d.method = QAWS_PARAM_BFF;
			TEST_ASSERT(qaws_mesh_3d_parameterize(uv, &r, m, seam, ns, &d) == QAWS_STATUS_UNSUPPORTED_OPERATION, "BFF reserved");
		}
		qaws_mesh_3d_destroy(m);
		free(uv);
	}
	{
		qaws_mesh_3d* m = dt_grid(4, 2);
		unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i;
		qaws_scalar* uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 6 * nt);
		qaws_param_report r;
		double max0 = -1e300, min1 = 1e300;
		qaws_status st = qaws_mesh_3d_parameterize(uv, &r, m, NULL, 0, NULL);
		for (i = 0; i < 3 * nt; i++)
		{
			if (i < 3 * nt / 2)
				max0 = fmax(max0, uv[2 * i]);
			else
				min1 = fmin(min1, uv[2 * i]);
		}
		sprintf(msg, "two grids: %u charts, the second from u %.3f past the first's end %.3f", r.charts, min1, max0);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && r.charts == 2 && min1 > max0, msg);
		qaws_mesh_3d_destroy(m);
		free(uv);
	}
}

int test_166_mesh_deform_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 166: deformation and parameterization\n");
	test_bend();
	test_turn();
	test_param();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
