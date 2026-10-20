/*
 * Test 151: adaptive surface tessellation (qaws_surface_tessellate_adaptive)
 *
 *   - a spherical cap at a chord tolerance: counts first, the same mesh on
 *     every call; the surface's distance to the mesh (at points inside each
 *     triangle, by its (u, v)) against the tolerance; no cracks: welded by (u, v),
 *     every edge off the domain's boundary is shared by two triangles, in
 *     opposite directions;
 *   - an angle tolerance: the normals at each triangle's corners;
 *   - a surface flat on one half, curved on the other: the triangles go to
 *     the curved half;
 *   - a short buffer, refusals.
 */

#include "test_common.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define TT_N 24
#define TT_VCAP 400000
#define TT_ICAP 1200000

static qaws_tessellation_vertex tt_v[TT_VCAP], tt_v2[TT_VCAP];
static unsigned int tt_i[TT_ICAP], tt_i2[TT_ICAP];

static qaws_surface* tt_surface(int half_flat)
{
	static qaws_vec3 p[TT_N * TT_N];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < TT_N; i++)
		for (j = 0; j < TT_N; j++)
		{
			double th = 0.4 + 2.2 * i / (TT_N - 1), ph = 2.5 * j / (TT_N - 1);
			if (half_flat)
			{
				p[i * TT_N + j].x = (qaws_scalar)i;
				p[i * TT_N + j].y = (qaws_scalar)j;
				p[i * TT_N + j].z = (qaws_scalar)(i < TT_N / 2 ? 0 : 1.5 * sin(1.3 * j) * cos(0.9 * i));
			}
			else
			{
				p[i * TT_N + j].x = (qaws_scalar)(3 * sin(th) * cos(ph));
				p[i * TT_N + j].y = (qaws_scalar)(3 * sin(th) * sin(ph));
				p[i * TT_N + j].z = (qaws_scalar)(3 * cos(th));
			}
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3; d.v_degree = 3; d.control_points = p; d.u_point_count = TT_N; d.v_point_count = TT_N;
	qaws_surface_create_bspline(&s, &d);
	return s;
}

/* the distance from p to the triangle abc, by projection on its plane and its edges */
static double tt_point_triangle(double const* p, qaws_vec3 a, qaws_vec3 b, qaws_vec3 c)
{
	double v[3][3] = { { a.x, a.y, a.z }, { b.x, b.y, b.z }, { c.x, c.y, c.z } }, n[3], e1[3], e2[3], l, h, best = 1e300, q[3];
	unsigned int k, m, inside = 1;
	for (k = 0; k < 3; k++)
	{
		e1[k] = v[1][k] - v[0][k];
		e2[k] = v[2][k] - v[0][k];
	}
	n[0] = e1[1] * e2[2] - e1[2] * e2[1];
	n[1] = e1[2] * e2[0] - e1[0] * e2[2];
	n[2] = e1[0] * e2[1] - e1[1] * e2[0];
	l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	h = ((p[0] - v[0][0]) * n[0] + (p[1] - v[0][1]) * n[1] + (p[2] - v[0][2]) * n[2]) / l;
	for (k = 0; k < 3; k++)
		q[k] = p[k] - h * n[k] / l;
	for (m = 0; m < 3; m++)
	{
		double const* s = v[m];
		double const* t = v[(m + 1) % 3];
		double ed[3] = { t[0] - s[0], t[1] - s[1], t[2] - s[2] }, w[3] = { q[0] - s[0], q[1] - s[1], q[2] - s[2] }, x[3], f, d;
		x[0] = ed[1] * w[2] - ed[2] * w[1];
		x[1] = ed[2] * w[0] - ed[0] * w[2];
		x[2] = ed[0] * w[1] - ed[1] * w[0];
		if (x[0] * n[0] + x[1] * n[1] + x[2] * n[2] < 0)
			inside = 0;
		f = (ed[0] * (p[0] - s[0]) + ed[1] * (p[1] - s[1]) + ed[2] * (p[2] - s[2])) / (ed[0] * ed[0] + ed[1] * ed[1] + ed[2] * ed[2]);
		f = f < 0 ? 0 : f > 1 ? 1 : f;
		d = sqrt((p[0] - s[0] - f * ed[0]) * (p[0] - s[0] - f * ed[0]) + (p[1] - s[1] - f * ed[1]) * (p[1] - s[1] - f * ed[1])
			+ (p[2] - s[2] - f * ed[2]) * (p[2] - s[2] - f * ed[2]));
		best = d < best ? d : best;
	}
	return inside ? fabs(h) : best;
}

#define TT_G 64

/* the largest distance from the surface to the mesh, at 10 points inside each triangle (by its (u, v)), each against the
   triangles near it in (u, v): grid cells list the triangles whose (u, v) box, widened by a cell, meets them */
static double tt_chord(qaws_surface const* s, unsigned int ni)
{
	static double const bary[10][3] = { { 1.0 / 3, 1.0 / 3, 1.0 / 3 }, { 0.5, 0.5, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 0.6, 0.2, 0.2 },
		{ 0.2, 0.6, 0.2 }, { 0.2, 0.2, 0.6 }, { 0.25, 0.75, 0 }, { 0, 0.25, 0.75 }, { 0.75, 0, 0.25 } };
	unsigned int* start = (unsigned int*)calloc(TT_G * TT_G + 1, sizeof(unsigned int));
	unsigned int* list = NULL;
	unsigned int pass, t, k, total = 0;
	double worst = 0;
	for (pass = 0; pass < 2; pass++)
	{
		for (t = 0; t < ni; t += 3)
		{
			double lo[2] = { 1, 1 }, hi[2] = { 0, 0 };
			int x0, x1, y0, y1, x, y;
			for (k = 0; k < 3; k++)
			{
				lo[0] = fmin(lo[0], tt_v[tt_i[t + k]].u); hi[0] = fmax(hi[0], tt_v[tt_i[t + k]].u);
				lo[1] = fmin(lo[1], tt_v[tt_i[t + k]].v); hi[1] = fmax(hi[1], tt_v[tt_i[t + k]].v);
			}
			x0 = (int)(lo[0] * TT_G) - 1; x1 = (int)(hi[0] * TT_G) + 1; y0 = (int)(lo[1] * TT_G) - 1; y1 = (int)(hi[1] * TT_G) + 1;
			for (x = x0 < 0 ? 0 : x0; x <= (x1 >= TT_G ? TT_G - 1 : x1); x++)
				for (y = y0 < 0 ? 0 : y0; y <= (y1 >= TT_G ? TT_G - 1 : y1); y++)
				{
					if (pass == 0)
						start[x * TT_G + y + 1]++;
					else
						list[start[x * TT_G + y]++] = t;
				}
		}
		if (pass == 0)
		{
			for (k = 1; k <= TT_G * TT_G; k++)
				start[k] += start[k - 1];
			total = start[TT_G * TT_G];
			list = (unsigned int*)malloc(sizeof(unsigned int) * (total ? total : 1));
		}
	}
	/* the fill moved each start to its cell's end: back by one cell */
	for (k = TT_G * TT_G; k > 0; k--)
		start[k] = start[k - 1];
	start[0] = 0;
	for (t = 0; t < ni; t += 3)
	{
		qaws_tessellation_vertex const* a = &tt_v[tt_i[t]];
		qaws_tessellation_vertex const* b = &tt_v[tt_i[t + 1]];
		qaws_tessellation_vertex const* c = &tt_v[tt_i[t + 2]];
		for (k = 0; k < 10; k++)
		{
			qaws_surface_eval_result r;
			double w0 = bary[k][0], w1 = bary[k][1], w2 = bary[k][2], p[3], u, v, best = 1e300;
			unsigned int cell, m;
			u = w0 * a->u + w1 * b->u + w2 * c->u;
			v = w0 * a->v + w1 * b->v + w2 * c->v;
			qaws_surface_evaluate(&r, s, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION);
			p[0] = r.position.x; p[1] = r.position.y; p[2] = r.position.z;
			cell = (unsigned int)(u * TT_G >= TT_G ? TT_G - 1 : u * TT_G) * TT_G + (unsigned int)(v * TT_G >= TT_G ? TT_G - 1 : v * TT_G);
			for (m = start[cell]; m < start[cell + 1]; m++)
			{
				unsigned int q = list[m];
				best = fmin(best, tt_point_triangle(p, tt_v[tt_i[q]].position, tt_v[tt_i[q + 1]].position, tt_v[tt_i[q + 2]].position));
			}
			worst = fmax(worst, best);
		}
	}
	free(start);
	free(list);
	return worst;
}

typedef struct tt_edge
{
	qaws_scalar a[2], b[2];             /* from, to: (u, v) */
	int forward;                        /* a < b in (u, v) order before the swap */
} tt_edge;

static int tt_less(qaws_scalar const* x, qaws_scalar const* y)
{
	return x[0] < y[0] || (x[0] == y[0] && x[1] < y[1]);
}

static int tt_cmp(void const* pa, void const* pb)
{
	tt_edge const* a = (tt_edge const*)pa;
	tt_edge const* b = (tt_edge const*)pb;
	if (tt_less(a->a, b->a)) return -1;
	if (tt_less(b->a, a->a)) return 1;
	if (tt_less(a->b, b->b)) return -1;
	if (tt_less(b->b, a->b)) return 1;
	return a->forward - b->forward;
}

/* edges off the domain's boundary not shared by exactly two triangles in opposite directions */
static unsigned int tt_cracks(qaws_surface const* s, unsigned int ni)
{
	static tt_edge e[TT_ICAP];
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	unsigned int n = 0, t, k, bad = 0;
	for (t = 0; t < ni; t += 3)
		for (k = 0; k < 3; k++)
		{
			qaws_tessellation_vertex const* p = &tt_v[tt_i[t + k]];
			qaws_tessellation_vertex const* q = &tt_v[tt_i[t + (k + 1) % 3]];
			qaws_scalar pa[2], pb[2];
			pa[0] = p->u; pa[1] = p->v; pb[0] = q->u; pb[1] = q->v;
			e[n].forward = tt_less(pa, pb);
			memcpy(e[n].a, e[n].forward ? pa : pb, sizeof pa);
			memcpy(e[n].b, e[n].forward ? pb : pa, sizeof pb);
			n++;
		}
	qsort(e, n, sizeof *e, tt_cmp);
	for (t = 0; t < n;)
	{
		unsigned int m = t + 1;
		int on_boundary;
		while (m < n && memcmp(e[m].a, e[t].a, sizeof e[t].a) == 0 && memcmp(e[m].b, e[t].b, sizeof e[t].b) == 0)
			m++;
		on_boundary = (e[t].a[0] == e[t].b[0] && (e[t].a[0] == ur.min_value || e[t].a[0] == ur.max_value))
			|| (e[t].a[1] == e[t].b[1] && (e[t].a[1] == vr.min_value || e[t].a[1] == vr.max_value));
		if (on_boundary ? m - t != 1 : (m - t != 2 || e[t].forward == e[t + 1].forward))
			bad++;
		t = m;
	}
	return bad;
}

static void test_chord_and_cracks(void)
{
	qaws_surface* s = tt_surface(0);
	qaws_surface_tessellation_desc d;
	unsigned int nv = 0, ni = 0, nv2 = 0, ni2 = 0, nv3, ni3, cracks;
	double err;
	char msg[240];
	memset(&d, 0, sizeof d);
	d.chord_tolerance = (qaws_scalar)2e-3;
	TEST_ASSERT(qaws_surface_tessellate_adaptive(s, &d, NULL, 0, &nv, NULL, 0, &ni) == QAWS_STATUS_OK && nv > 0 && ni > 0, "cap: counts first");
	TEST_ASSERT(nv <= TT_VCAP && ni <= TT_ICAP, "cap: within the test's buffers");
	TEST_ASSERT(qaws_surface_tessellate_adaptive(s, &d, tt_v, TT_VCAP, &nv2, tt_i, TT_ICAP, &ni2) == QAWS_STATUS_OK && nv2 == nv && ni2 == ni,
		"cap: the same counts filled");
	qaws_surface_tessellate_adaptive(s, &d, tt_v2, TT_VCAP, &nv3, tt_i2, TT_ICAP, &ni3);
	TEST_ASSERT(nv3 == nv && ni3 == ni && memcmp(tt_v, tt_v2, sizeof(qaws_tessellation_vertex) * nv) == 0
		&& memcmp(tt_i, tt_i2, sizeof(unsigned int) * ni) == 0, "cap: the same mesh every call");
	err = tt_chord(s, ni);
	cracks = tt_cracks(s, ni);
	snprintf(msg, sizeof msg, "cap at chord 2e-3: %u triangles, %u vertices, chord error %.3e, unmatched inner edges %u", ni / 3, nv, err, cracks);
	printf("  %s\n", msg);
	TEST_ASSERT(err <= 1.5 * 2e-3, msg);
	TEST_ASSERT(cracks == 0, "no cracks: every inner edge shared by two triangles in opposite directions");
	TEST_ASSERT(qaws_surface_tessellate_adaptive(s, &d, tt_v2, 10, &nv3, tt_i2, 30, &ni3) == QAWS_STATUS_BUFFER_TOO_SMALL && nv3 == nv && ni3 == ni,
		"short buffers: the counts");
	qaws_surface_destroy(s);
}

static void test_angle(void)
{
	qaws_surface* s = tt_surface(0);
	qaws_surface_tessellation_desc d;
	unsigned int nv = 0, ni = 0, t;
	double worst = 0;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.angle_tolerance = (qaws_scalar)0.15;
	qaws_surface_tessellate_adaptive(s, &d, tt_v, TT_VCAP, &nv, tt_i, TT_ICAP, &ni);
	for (t = 0; t < ni; t += 3)
	{
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 a = tt_v[tt_i[t + k]].normal, b = tt_v[tt_i[t + (k + 1) % 3]].normal;
			double x = (double)a.y * b.z - (double)a.z * b.y, y = (double)a.z * b.x - (double)a.x * b.z, z = (double)a.x * b.y - (double)a.y * b.x;
			worst = fmax(worst, atan2(sqrt(x * x + y * y + z * z), (double)a.x * b.x + (double)a.y * b.y + (double)a.z * b.z));
		}
	}
	snprintf(msg, sizeof msg, "cap at angle 0.15: %u triangles, largest angle between corner normals %.4f", ni / 3, worst);
	printf("  %s\n", msg);
	TEST_ASSERT(worst <= 0.15 && tt_cracks(s, ni) == 0, msg);
	qaws_surface_destroy(s);
}

static void test_adaptivity(void)
{
	qaws_surface* s = tt_surface(1);
	qaws_surface_tessellation_desc d;
	unsigned int nv = 0, ni = 0, t, flat = 0, curved = 0;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.chord_tolerance = (qaws_scalar)1e-2;
	qaws_surface_tessellate_adaptive(s, &d, tt_v, TT_VCAP, &nv, tt_i, TT_ICAP, &ni);
	for (t = 0; t < ni; t += 3)
	{
		double u = (tt_v[tt_i[t]].u + tt_v[tt_i[t + 1]].u + tt_v[tt_i[t + 2]].u) / 3;
		/* the flat half ends where the curved points start to act: u < 8 / 21 */
		if (u < 8.0 / 21)
			flat++;
		else
			curved++;
	}
	snprintf(msg, sizeof msg, "half flat: %u triangles on the flat part, %u on the curved part; cracks %u, chord error %.3e", flat, curved, tt_cracks(s, ni), tt_chord(s, ni));
	printf("  %s\n", msg);
	TEST_ASSERT(curved > 50 * flat && tt_cracks(s, ni) == 0 && tt_chord(s, ni) <= 1.5e-2, msg);
	qaws_surface_destroy(s);
}

static void test_refusals(void)
{
	qaws_surface* s = tt_surface(0);
	qaws_surface_tessellation_desc d;
	unsigned int nv, ni;
	memset(&d, 0, sizeof d);
	TEST_ASSERT(qaws_surface_tessellate_adaptive(s, &d, NULL, 0, &nv, NULL, 0, &ni) == QAWS_STATUS_INVALID_ARGUMENT, "no tolerance refused");
	d.angle_tolerance = -1;
	TEST_ASSERT(qaws_surface_tessellate_adaptive(s, &d, NULL, 0, &nv, NULL, 0, &ni) == QAWS_STATUS_INVALID_ARGUMENT, "a negative tolerance refused");
	d.angle_tolerance = (qaws_scalar)0.5;
	TEST_ASSERT(qaws_surface_tessellate_adaptive(s, &d, NULL, 0, NULL, NULL, 0, &ni) == QAWS_STATUS_INVALID_ARGUMENT, "no count refused");
	TEST_ASSERT(qaws_surface_tessellate_adaptive(NULL, &d, NULL, 0, &nv, NULL, 0, &ni) == QAWS_STATUS_INVALID_ARGUMENT, "no surface refused");
	qaws_surface_destroy(s);
}

int test_151_surface_tessellation_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 151: adaptive surface tessellation\n");
	test_chord_and_cracks();
	test_angle();
	test_adaptivity();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
