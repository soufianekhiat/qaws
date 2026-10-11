/*
 * Test 99: straight skeletons (qaws_skeleton)
 *
 *   Every shape is checked the same way:
 *     - every arc's ends lie at weighted distance w t from both its faces'
 *       edges (the defining property of the skeleton)
 *     - the faces tile the region: their areas sum to the region's
 *     - the wavefront at time t (a separate run of the events, stopped at
 *       t) encloses the area the faces leave unswept by t: each face clipped
 *       by its edge's line moved to w t
 *   and against closed forms: the square, rectangles (the ridge), regular
 *   polygons (every corner meeting at the centre at the inradius), a square
 *   with a centred square hole (every event at one time), weighted
 *   rectangles.
 */

#include "test_common.h"
#include "qaws_skeleton.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define SKT_PI 3.14159265358979323846
#define SKT_F32 (sizeof(qaws_scalar) == 4)

static int g_verbose_skt = 0;

typedef struct skt_shape
{
	char const* name;
	unsigned int nloops;
	unsigned int count[16];
	double pts[512];
} skt_shape;

/* the face's edge line, n . p = d, n inward */
static void skt_face_line(qaws_skeleton const* s, qaws_skeleton_face const* f, double* n, double* d)
{
	qaws_skeleton_node a, b;
	double l;
	qaws_skeleton_get_node(s, f->nodes[0], &a);
	qaws_skeleton_get_node(s, f->nodes[1], &b);
	l = hypot((double)b.position.x - a.position.x, (double)b.position.y - a.position.y);
	n[0] = -((double)b.position.y - a.position.y) / l;
	n[1] = ((double)b.position.x - a.position.x) / l;
	*d = n[0] * a.position.x + n[1] * a.position.y;
}

/* the area of the part of the polygon with n . p - d >= w t (unswept) */
static double skt_clip_area(double const* xy, unsigned int n, double const* nr, double d, double wt)
{
	double area = 0, buf[2 * 256], *q = buf;
	unsigned int i, m = 0;
	if (n > 120)
		q = (double*)malloc(sizeof(double) * 4 * n);
	for (i = 0; i < n; i++)
	{
		double const* a = xy + 2 * i;
		double const* b = xy + 2 * ((i + 1) % n);
		double fa = nr[0] * a[0] + nr[1] * a[1] - d - wt, fb = nr[0] * b[0] + nr[1] * b[1] - d - wt;
		if (fa >= 0) { q[2 * m] = a[0]; q[2 * m + 1] = a[1]; m++; }
		if ((fa >= 0) != (fb >= 0))
		{
			double s = fa / (fa - fb);
			q[2 * m] = a[0] + s * (b[0] - a[0]);
			q[2 * m + 1] = a[1] + s * (b[1] - a[1]);
			m++;
		}
	}
	for (i = 0; i < m; i++)
	{
		unsigned int j = (i + 1) % m;
		area += q[2 * i] * q[2 * j + 1] - q[2 * j] * q[2 * i + 1];
	}
	if (q != buf) free(q);
	return 0.5 * area;
}

static double skt_poly_area(double const* xy, unsigned int n)
{
	double a = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		unsigned int j = (i + 1) % n;
		a += xy[2 * i] * xy[2 * j + 1] - xy[2 * j] * xy[2 * i + 1];
	}
	return 0.5 * a;
}

static double skt_result_area(qaws_clip_result const* r)
{
	unsigned int i, n = r ? qaws_clip_result_get_path_count(r) : 0;
	double a = 0;
	for (i = 0; i < n; i++)
	{
		qaws_path_2d p;
		qaws_scalar v = 0;
		qaws_clip_result_get_path(r, i, &p);
		qaws_path_compute_area_2d(&p, &v);
		a += (double)v;
	}
	return a;
}

typedef struct skt_weights
{
	double w[64];
} skt_weights;

static qaws_scalar skt_weight(void* user, unsigned int path, unsigned int curve, unsigned int span)
{
	skt_weights const* w = (skt_weights const*)user;
	(void)path; (void)curve;
	return (qaws_scalar)w->w[span % 64];
}

/* builds the skeleton and checks it; returns it for closed-form checks */
static qaws_skeleton* skt_check(skt_shape const* sh, skt_weights const* wts)
{
	qaws_curve* cs[16];
	qaws_curve const* cv[16];
	qaws_path_2d p[16];
	qaws_skeleton_desc d;
	qaws_skeleton* s = NULL;
	qaws_status st;
	unsigned int l, i, off = 0, k;
	double area = 0, worst_dist = 0, faces = 0, worst_off = 0, size, tmax, tol;
	qaws_scalar buf[512];
	char msg[400];
	for (l = 0; l < sh->nloops; l++)
	{
		for (i = 0; i < 2 * sh->count[l]; i++) buf[i] = (qaws_scalar)sh->pts[2 * off + i];
		qaws_curve_create_polyline_2d(buf, sh->count[l], 1, &cs[l]);
		cv[l] = cs[l];
		p[l].curves = &cv[l]; p[l].curve_count = 1; p[l].closed = 1;
		area += skt_poly_area(sh->pts + 2 * off, sh->count[l]);
		off += sh->count[l];
	}
	{
		/* the region's area: the paths united */
		qaws_clip_result* u = NULL;
		if (qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, p, sh->nloops, NULL, 0, &u) == QAWS_STATUS_OK)
			area = skt_result_area(u);
		qaws_clip_result_destroy(u);
	}
	size = sqrt(fabs(area));
	tol = (SKT_F32 ? 2e-5 : 1e-9) * size;
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = sh->nloops; d.fill_rule = QAWS_FILL_NON_ZERO;
	if (wts) { d.weight_fn = skt_weight; d.weight_user = (void*)wts; }
	st = qaws_skeleton_create(&d, &s);
	sprintf(msg, "%s: created (status %d)", sh->name, (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && s, msg);
	if (!s)
	{
		for (l = 0; l < sh->nloops; l++) qaws_curve_destroy(cs[l]);
		return NULL;
	}
	tmax = (double)qaws_skeleton_get_max_time(s);
	/* arcs at weighted distance w t from both faces' edges */
	for (i = 0; i < qaws_skeleton_get_arc_count(s); i++)
	{
		qaws_skeleton_arc a;
		qaws_skeleton_get_arc(s, i, &a);
		for (k = 0; k < 4; k++)
		{
			qaws_skeleton_face f;
			qaws_skeleton_node nd;
			double n[2], dd, e;
			qaws_skeleton_get_face(s, (k & 1) ? a.face_b : a.face_a, &f);
			qaws_skeleton_get_node(s, (k & 2) ? a.node_b : a.node_a, &nd);
			skt_face_line(s, &f, n, &dd);
			e = fabs(n[0] * nd.position.x + n[1] * nd.position.y - dd - (double)f.weight * nd.time);
			if (e > worst_dist) worst_dist = e;
		}
	}
	/* faces tile the region */
	for (i = 0; i < qaws_skeleton_get_face_count(s); i++)
	{
		qaws_skeleton_face f;
		double xy[2 * 256];
		qaws_skeleton_get_face(s, i, &f);
		for (k = 0; k < f.node_count && k < 256; k++)
		{
			qaws_skeleton_node nd;
			qaws_skeleton_get_node(s, f.nodes[k], &nd);
			xy[2 * k] = nd.position.x; xy[2 * k + 1] = nd.position.y;
		}
		faces += skt_poly_area(xy, f.node_count);
	}
	/* the wavefront's area at five times against the faces' unswept parts */
	for (k = 1; k <= 5; k++)
	{
		double t = tmax * k / 6, unswept = 0, got;
		qaws_clip_result* r = NULL;
		for (i = 0; i < qaws_skeleton_get_face_count(s); i++)
		{
			qaws_skeleton_face f;
			double xy[2 * 256], n[2], dd;
			unsigned int j;
			qaws_skeleton_get_face(s, i, &f);
			for (j = 0; j < f.node_count && j < 256; j++)
			{
				qaws_skeleton_node nd;
				qaws_skeleton_get_node(s, f.nodes[j], &nd);
				xy[2 * j] = nd.position.x; xy[2 * j + 1] = nd.position.y;
			}
			skt_face_line(s, &f, n, &dd);
			unswept += skt_clip_area(xy, f.node_count, n, dd, (double)f.weight * t);
		}
		qaws_skeleton_offset(s, (qaws_scalar)t, &r);
		got = skt_result_area(r);
		if (fabs(got - unswept) > worst_off) worst_off = fabs(got - unswept);
		qaws_clip_result_destroy(r);
	}
	sprintf(msg, "%s: %u nodes, %u arcs, max time %.9f; arcs' distance error %.2g, faces' area %.12f (region %.12f), wavefront area error %.2g",
		sh->name, qaws_skeleton_get_node_count(s), qaws_skeleton_get_arc_count(s), tmax, worst_dist, faces, area, worst_off);
	TEST_ASSERT(worst_dist <= 10 * tol && fabs(faces - area) <= (SKT_F32 ? 1e-4 : 1e-9) * area &&
		worst_off <= (SKT_F32 ? 1e-4 : 1e-9) * area, msg);
	if (g_verbose_skt) printf("    %s\n", msg);
	for (l = 0; l < sh->nloops; l++) qaws_curve_destroy(cs[l]);
	return s;
}

static void skt_add(skt_shape* sh, double const* xy, unsigned int n)
{
	unsigned int i, off = 0;
	for (i = 0; i < sh->nloops; i++) off += sh->count[i];
	memcpy(sh->pts + 2 * off, xy, sizeof(double) * 2 * n);
	sh->count[sh->nloops++] = n;
}

static void skt_regular(skt_shape* sh, unsigned int n, double r, double cx, double cy, int cw)
{
	double xy[256];
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double a = 2 * SKT_PI * (cw ? n - i : i) / n;
		xy[2 * i] = cx + r * cos(a);
		xy[2 * i + 1] = cy + r * sin(a);
	}
	skt_add(sh, xy, n);
}

static void skt_free(qaws_skeleton* s)
{
	qaws_skeleton_destroy(s);
}

/* the longest arc with both ends at time t */
static double skt_ridge(qaws_skeleton const* s, double t)
{
	unsigned int i;
	double best = 0;
	for (i = 0; i < qaws_skeleton_get_arc_count(s); i++)
	{
		qaws_skeleton_arc a;
		qaws_skeleton_node x, y;
		qaws_skeleton_get_arc(s, i, &a);
		qaws_skeleton_get_node(s, a.node_a, &x);
		qaws_skeleton_get_node(s, a.node_b, &y);
		if (fabs((double)x.time - t) < 1e-6 && fabs((double)y.time - t) < 1e-6)
		{
			double l = hypot((double)x.position.x - y.position.x, (double)x.position.y - y.position.y);
			if (l > best) best = l;
		}
	}
	return best;
}

static void test_closed_forms(void)
{
	static double const square[8] = { 0, 0, 2, 0, 2, 2, 0, 2 };
	static double const rect[8] = { 0, 0, 3, 0, 3, 1, 0, 1 };
	static double const hole[8] = { 1, 1, 1, 3, 3, 3, 3, 1 };
	static double const outer[8] = { 0, 0, 4, 0, 4, 4, 0, 4 };
	skt_shape sh;
	qaws_skeleton* s;
	char msg[300];
	unsigned int n;
	double tol = SKT_F32 ? 1e-5 : 1e-12;
	memset(&sh, 0, sizeof(sh)); sh.name = "square"; skt_add(&sh, square, 4);
	s = skt_check(&sh, NULL);
	sprintf(msg, "square 2 x 2: max time %.15f (1), %u nodes (5), %u arcs (4)", s ? (double)qaws_skeleton_get_max_time(s) : -1,
		s ? qaws_skeleton_get_node_count(s) : 0, s ? qaws_skeleton_get_arc_count(s) : 0);
	TEST_ASSERT(s && fabs((double)qaws_skeleton_get_max_time(s) - 1) <= tol && qaws_skeleton_get_node_count(s) == 5 && qaws_skeleton_get_arc_count(s) == 4, msg);
	skt_free(s);
	memset(&sh, 0, sizeof(sh)); sh.name = "rectangle 3 x 1"; skt_add(&sh, rect, 4);
	s = skt_check(&sh, NULL);
	sprintf(msg, "rectangle 3 x 1: max time %.15f (0.5), ridge %.15f (2), %u nodes (6), %u arcs (5)", s ? (double)qaws_skeleton_get_max_time(s) : -1,
		s ? skt_ridge(s, 0.5) : -1, s ? qaws_skeleton_get_node_count(s) : 0, s ? qaws_skeleton_get_arc_count(s) : 0);
	TEST_ASSERT(s && fabs((double)qaws_skeleton_get_max_time(s) - 0.5) <= tol && fabs(skt_ridge(s, 0.5) - 2) <= 10 * tol &&
		qaws_skeleton_get_node_count(s) == 6 && qaws_skeleton_get_arc_count(s) == 5, msg);
	skt_free(s);
	for (n = 3; n <= 24; n += 7)
	{
		memset(&sh, 0, sizeof(sh)); sh.name = "regular polygon"; skt_regular(&sh, n, 1, 0.3, -0.2, 0);
		s = skt_check(&sh, NULL);
		sprintf(msg, "regular %u-gon: max time %.15f (the inradius %.15f), %u nodes (%u: all meet at the centre)", n,
			s ? (double)qaws_skeleton_get_max_time(s) : -1, cos(SKT_PI / n), s ? qaws_skeleton_get_node_count(s) : 0, n + 1);
		TEST_ASSERT(s && fabs((double)qaws_skeleton_get_max_time(s) - cos(SKT_PI / n)) <= tol && qaws_skeleton_get_node_count(s) == n + 1, msg);
		skt_free(s);
	}
	memset(&sh, 0, sizeof(sh)); sh.name = "square with a centred square hole"; skt_add(&sh, outer, 4); skt_add(&sh, hole, 4);
	s = skt_check(&sh, NULL);
	sprintf(msg, "4 x 4 square with a centred 2 x 2 hole: max time %.15f (0.5), every node at 0 or 0.5", s ? (double)qaws_skeleton_get_max_time(s) : -1);
	{
		int ok = s != NULL;
		unsigned int i;
		for (i = 0; s && i < qaws_skeleton_get_node_count(s); i++)
		{
			qaws_skeleton_node nd;
			qaws_skeleton_get_node(s, i, &nd);
			if (fabs((double)nd.time) > tol && fabs((double)nd.time - 0.5) > tol) ok = 0;
		}
		TEST_ASSERT(ok && fabs((double)qaws_skeleton_get_max_time(s) - 0.5) <= tol, msg);
	}
	skt_free(s);
	{
		/* weights 1 on the long sides, 2 on the short: the short sides'
		   faces reach 2 / sqrt(5) of ... the ridge ends where w t meets:
		   long sides meet at t = 0.5; the ridge is 3 - 2 * (0.5 * 2) = 1 */
		skt_weights w;
		memset(&w, 0, sizeof(w));
		w.w[0] = 1; w.w[1] = 2; w.w[2] = 1; w.w[3] = 2;
		memset(&sh, 0, sizeof(sh)); sh.name = "weighted rectangle"; skt_add(&sh, rect, 4);
		s = skt_check(&sh, &w);
		sprintf(msg, "rectangle 3 x 1, short sides at speed 2: max time %.15f (0.5), ridge %.15f (1)", s ? (double)qaws_skeleton_get_max_time(s) : -1,
			s ? skt_ridge(s, 0.5) : -1);
		TEST_ASSERT(s && fabs((double)qaws_skeleton_get_max_time(s) - 0.5) <= tol && fabs(skt_ridge(s, 0.5) - 1) <= 10 * tol, msg);
		skt_free(s);
	}
}

static void test_shapes(void)
{
	static double const ell[12] = { 0, 0, 4, 0, 4, 1, 1, 1, 1, 3, 0, 3 };
	static double const comb[24] = { 0, 0, 7, 0, 7, 3, 6, 3, 6, 1, 5, 1, 5, 3, 3, 3, 3, 1, 1, 1, 1, 3, 0, 3 };
	static double const plus[24] = { 1, 0, 2, 0, 2, 1, 3, 1, 3, 2, 2, 2, 2, 3, 1, 3, 1, 2, 0, 2, 0, 1, 1, 1 };
	static double const plan[32] = { 0, 0, 6, 0, 6, 2, 9, 2, 9, 6, 7, 6, 7, 8, 3, 8, 3, 5, 1, 5, 1, 4, 0, 4, 0, 2.5, -1.5, 2.5, -1.5, 1, 0, 1 };
	static double const spiral[28] = { 0, 0, 8, 0, 8, 7, 2, 7, 2, 3, 5, 3, 5, 4.2, 3.2, 4.2, 3.2, 5.8, 6.8, 5.8, 6.8, 1.2, 1.2, 1.2, 1.2, 7, 0, 7 };
	skt_shape sh;
	unsigned int i;
	memset(&sh, 0, sizeof(sh)); sh.name = "L"; skt_add(&sh, ell, 6); skt_free(skt_check(&sh, NULL));
	memset(&sh, 0, sizeof(sh)); sh.name = "comb"; skt_add(&sh, comb, 12); skt_free(skt_check(&sh, NULL));
	memset(&sh, 0, sizeof(sh)); sh.name = "plus (every event degenerate)"; skt_add(&sh, plus, 12); skt_free(skt_check(&sh, NULL));
	memset(&sh, 0, sizeof(sh)); sh.name = "floor plan"; skt_add(&sh, plan, 16); skt_free(skt_check(&sh, NULL));
	memset(&sh, 0, sizeof(sh)); sh.name = "orthogonal spiral"; skt_add(&sh, spiral, 14); skt_free(skt_check(&sh, NULL));
	{
		double xy[2 * 14];
		for (i = 0; i < 14; i++)
		{
			double a = 2 * SKT_PI * i / 14, r = (i & 1) ? 0.45 : 1.0;
			xy[2 * i] = r * cos(a); xy[2 * i + 1] = r * sin(a);
		}
		memset(&sh, 0, sizeof(sh)); sh.name = "7-point star"; skt_add(&sh, xy, 14); skt_free(skt_check(&sh, NULL));
	}
	{
		/* random star-shaped polygons, with a random hole */
		unsigned long long st = 12345;
		unsigned int r;
		for (r = 0; r < 6; r++)
		{
			double xy[2 * 60];
			unsigned int n = 20 + 8 * r;
			for (i = 0; i < n; i++)
			{
				double a = 2 * SKT_PI * i / n, rr;
				st = st * 6364136223846793005ull + 1442695040888963407ull;
				rr = 1.0 + 0.6 * ((st >> 11) * (1.0 / 9007199254740992.0));
				xy[2 * i] = rr * cos(a); xy[2 * i + 1] = rr * sin(a);
			}
			memset(&sh, 0, sizeof(sh)); sh.name = "random star-shaped polygon"; skt_add(&sh, xy, n);
			if (r & 1)
			{
				sh.name = "random star-shaped polygon with a hole";
				skt_regular(&sh, 5 + r, 0.4, 0.1 * r - 0.2, 0.05, 1);
			}
			skt_free(skt_check(&sh, NULL));
		}
	}
	{
		/* unions of random rectangles on an integer grid: orthogonal, nearly
		   every event degenerate (equal times, collinear edges, corridors) */
		unsigned long long st = 777;
		unsigned int r;
		for (r = 0; r < 12; r++)
		{
			unsigned int k, nr = 3 + r % 8;
			memset(&sh, 0, sizeof(sh));
			sh.name = "union of grid rectangles";
			for (k = 0; k < nr; k++)
			{
				double x0, y0, w, h, xy[8];
				st = st * 6364136223846793005ull + 1442695040888963407ull; x0 = (double)((st >> 33) % 8);
				st = st * 6364136223846793005ull + 1442695040888963407ull; y0 = (double)((st >> 33) % 8);
				st = st * 6364136223846793005ull + 1442695040888963407ull; w = 1 + (double)((st >> 33) % 5);
				st = st * 6364136223846793005ull + 1442695040888963407ull; h = 1 + (double)((st >> 33) % 5);
				xy[0] = x0; xy[1] = y0; xy[2] = x0 + w; xy[3] = y0; xy[4] = x0 + w; xy[5] = y0 + h; xy[6] = x0; xy[7] = y0 + h;
				skt_add(&sh, xy, 4);
			}
			skt_free(skt_check(&sh, NULL));
		}
	}
}

int test_99_skeleton_main(void)
{
	g_pass = 0;
	g_fail = 0;
	g_verbose_skt = getenv("QAWS_SKT_VERBOSE") != NULL;
	printf("Test 99: straight skeletons\n");
	test_closed_forms();
	test_shapes();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
