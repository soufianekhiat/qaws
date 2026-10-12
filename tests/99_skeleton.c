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

/* ------------------------------------------------------------------------ */
/*  Curved: lines and circular arcs                                         */
/* ------------------------------------------------------------------------ */

typedef struct skt_seg
{
	int arc;                    /* line: x0 y0 x1 y1; arc: cx cy r a0 a1 */
	double a[5];
} skt_seg;

typedef struct skt_cshape
{
	char const* name;
	unsigned int nloops;
	unsigned int count[8];
	skt_seg seg[64];
	double sigma;
	int round;                  /* the wavefront is the round offset (no reflex vertex kept) */
	double (*exact)(double t);  /* the wavefront's area in closed form, in place of the round offset */
	double exact_until;         /* the closed form holds up to this time */
} skt_cshape;

static void skt_cline(skt_cshape* sh, unsigned int* n, double x0, double y0, double x1, double y1)
{
	skt_seg* g = &sh->seg[(*n)++];
	g->arc = 0; g->a[0] = x0; g->a[1] = y0; g->a[2] = x1; g->a[3] = y1;
}

static void skt_carc(skt_cshape* sh, unsigned int* n, double cx, double cy, double r, double a0, double a1)
{
	skt_seg* g = &sh->seg[(*n)++];
	g->arc = 1; g->a[0] = cx; g->a[1] = cy; g->a[2] = r; g->a[3] = a0; g->a[4] = a1;
}

static qaws_curve* skt_seg_curve(skt_seg const* g)
{
	qaws_curve* c = NULL;
	if (!g->arc)
	{
		qaws_scalar cp[4];
		qaws_bezier_desc d;
		cp[0] = (qaws_scalar)g->a[0]; cp[1] = (qaws_scalar)g->a[1]; cp[2] = (qaws_scalar)g->a[2]; cp[3] = (qaws_scalar)g->a[3];
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 1; d.control_points = cp; d.control_point_count = 2;
		qaws_curve_create_bezier(&d, &c);
	}
	else
	{
		qaws_arc_segment s;
		qaws_arc_desc d;
		memset(&s, 0, sizeof(s));
		s.center[0] = (qaws_scalar)g->a[0]; s.center[1] = (qaws_scalar)g->a[1]; s.radius = (qaws_scalar)g->a[2];
		s.angle_start = (qaws_scalar)g->a[3]; s.angle_end = (qaws_scalar)g->a[4];
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = &s; d.segment_count = 1;
		qaws_curve_create_arc(&d, &c);
	}
	return c;
}

/* the support's residual at (p, t): 0 on it */
static double skt_face_resid(qaws_skeleton_face const* f, double x, double y, double t)
{
	if (f->circle)
		return f->side * ((double)f->radius - hypot(x - f->center.x, y - f->center.y)) - t;
	return (double)f->normal.x * x + (double)f->normal.y * y - (double)f->offset - (double)f->weight * t;
}

/* the distance from p to the nearest of the curves */
static double skt_boundary_dist(qaws_curve* const* cs, unsigned int n, double x, double y)
{
	double best = 1e300;
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		qaws_vec2 q;
		qaws_scalar t = 0;
		qaws_eval_result_2d e;
		double dd;
		q.x = (qaws_scalar)x; q.y = (qaws_scalar)y;
		qaws_curve_find_closest_parameter_2d(cs[i], q, &t);
		qaws_curve_evaluate_2d(cs[i], t, QAWS_EVAL_FLAG_POSITION, &e);
		dd = hypot((double)e.position.x - x, (double)e.position.y - y);
		if (dd < best) best = dd;
	}
	return best;
}

/* a face's area: its edge, then its arcs, sampled */
static double skt_face_area(qaws_skeleton const* s, unsigned int fi, unsigned int* missing)
{
	qaws_skeleton_face f;
	double area = 0, px = 0, py = 0;
	int have = 0;
	unsigned int k, na, i, m = 256, start;
	unsigned char* used = (unsigned char*)calloc(qaws_skeleton_get_arc_count(s) + 1, 1);
	qaws_skeleton_get_face(s, fi, &f);
	na = f.node_count;
#define SKT_PT(X, Y) do { double xx_ = (X), yy_ = (Y); if (have) area += px * yy_ - xx_ * py; px = xx_; py = yy_; have = 1; } while (0)
	{
		qaws_skeleton_node a;
		qaws_skeleton_get_node(s, f.nodes[0], &a);
		SKT_PT(a.position.x, a.position.y);
	}
	start = 0;
	if (f.edge)
	{
		qaws_range rg = qaws_curve_get_parameter_range(f.edge);
		for (k = 1; k <= m; k++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(f.edge, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * k / m), QAWS_EVAL_FLAG_POSITION, &e);
			SKT_PT(e.position.x, e.position.y);
		}
		{
			/* a closed edge ends where it starts: the arcs start at node 0 */
			qaws_eval_result_2d e0, e1;
			qaws_curve_evaluate_2d(f.edge, rg.min_value, QAWS_EVAL_FLAG_POSITION, &e0);
			qaws_curve_evaluate_2d(f.edge, rg.max_value, QAWS_EVAL_FLAG_POSITION, &e1);
			start = (na > 1 && hypot((double)e1.position.x - e0.position.x, (double)e1.position.y - e0.position.y) > (SKT_F32 ? 1e-5 : 1e-9)) ? 1 : 0;
		}
	}
	for (k = start; k < na; k++)
	{
		unsigned int a = f.nodes[k], b = f.nodes[(k + 1) % na], found = 0, pass;
		/* an unused arc first; a slit (a ray in and out) goes back along one */
		for (pass = 0; pass < 2 && !found; pass++)
		for (i = 0; i < qaws_skeleton_get_arc_count(s) && !found; i++)
		{
			qaws_skeleton_arc ar;
			qaws_skeleton_get_arc(s, i, &ar);
			if ((used[i] && pass == 0) || (ar.face_a != fi && ar.face_b != fi))
				continue;
			if ((ar.node_a == a && ar.node_b == b) || (ar.node_a == b && ar.node_b == a))
			{
				unsigned int j;
				int rev = ar.node_a != a;
				used[i] = 1;
				found = 1;
				for (j = 1; j <= 64; j++)
				{
					qaws_vec2 q;
					qaws_scalar tt;
					qaws_skeleton_evaluate_arc(s, i, (qaws_scalar)(rev ? 1 - j / 64.0 : j / 64.0), &q, &tt);
					SKT_PT(q.x, q.y);
				}
			}
		}
		if (!found) (*missing)++;
	}
#undef SKT_PT
	free(used);
	return 0.5 * area;
}

static void skt_check_curved(skt_cshape const* sh)
{
	qaws_curve* cs[64];
	qaws_curve const* cv[64];
	qaws_path_2d p[8];
	qaws_skeleton_desc d;
	qaws_skeleton* s = NULL;
	qaws_status st;
	unsigned int l, i, k, off = 0, ns = 0, missing = 0;
	double area = 0, size, worst_dist = 0, worst_axis = 0, worst_round = 0, faces = 0, tmax;
	char msg[500];
	for (l = 0; l < sh->nloops; l++)
	{
		for (i = 0; i < sh->count[l]; i++) { cs[off + i] = skt_seg_curve(&sh->seg[off + i]); cv[off + i] = cs[off + i]; }
		p[l].curves = cv + off; p[l].curve_count = sh->count[l]; p[l].closed = 1;
		off += sh->count[l];
	}
	ns = off;
	{
		qaws_clip_result* u = NULL;
		if (qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, p, sh->nloops, NULL, 0, &u) == QAWS_STATUS_OK)
			area = skt_result_area(u);
		qaws_clip_result_destroy(u);
	}
	size = sqrt(fabs(area));
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = sh->nloops; d.fill_rule = QAWS_FILL_NON_ZERO;
	d.splice_angle = (qaws_scalar)sh->sigma;
	st = qaws_skeleton_create(&d, &s);
	sprintf(msg, "%s: created (status %d)", sh->name, (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && s, msg);
	if (!s)
	{
		for (i = 0; i < ns; i++) qaws_curve_destroy(cs[i]);
		return;
	}
	tmax = (double)qaws_skeleton_get_max_time(s);
	/* every arc's points on both faces' supports at their times; at sigma =
	   pi, the non-smooth ones at distance t from the boundary */
	for (i = 0; i < qaws_skeleton_get_arc_count(s); i++)
	{
		qaws_skeleton_arc a;
		qaws_skeleton_face fa, fb;
		qaws_skeleton_get_arc(s, i, &a);
		qaws_skeleton_get_face(s, a.face_a, &fa);
		qaws_skeleton_get_face(s, a.face_b, &fb);
		for (k = 0; k <= 8; k++)
		{
			qaws_vec2 q;
			qaws_scalar t;
			double e1, e2;
			qaws_skeleton_evaluate_arc(s, i, (qaws_scalar)(k / 8.0), &q, &t);
			e1 = fabs(skt_face_resid(&fa, q.x, q.y, t));
			e2 = fabs(skt_face_resid(&fb, q.x, q.y, t));
			if (e1 > worst_dist) worst_dist = e1;
			if (e2 > worst_dist) worst_dist = e2;
			if (sh->sigma == SKT_PI && !a.smooth)
			{
				double e = fabs(skt_boundary_dist(cs, ns, q.x, q.y) - (double)t);
				if (e > worst_axis) worst_axis = e;
			}
		}
	}
	/* the faces tile the region */
	for (i = 0; i < qaws_skeleton_get_face_count(s); i++)
		faces += skt_face_area(s, i, &missing);
	/* the wavefront against the round offset */
	if (sh->round)
		for (k = 1; k <= 5; k++)
		{
			double t = tmax * k / 6, a1, a2;
			qaws_clip_result* r1 = NULL;
			qaws_clip_result* r2 = NULL;
			if (sh->exact && t > sh->exact_until)
				continue;
			qaws_skeleton_offset(s, (qaws_scalar)t, &r1);
			if (!sh->exact)
				qaws_offset_paths(p, sh->nloops, (qaws_scalar)-t, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &r2);
			a1 = skt_result_area(r1);
			a2 = sh->exact ? sh->exact(t) : skt_result_area(r2);
			if (fabs(a1 - a2) > worst_round) worst_round = fabs(a1 - a2);
			if (g_verbose_skt && fabs(a1 - a2) > 1e-3) printf("      t %.6f: skeleton offset %.9f, round offset %.9f\n", t, a1, a2);
			qaws_clip_result_destroy(r1);
			qaws_clip_result_destroy(r2);
		}
	sprintf(msg, "%s (sigma %.3f): %u nodes, %u arcs, max time %.9f; supports' error %.2g, medial axis' distance error %.2g, faces %.9f (region %.9f, %u gaps), round offset area error %.2g",
		sh->name, sh->sigma, qaws_skeleton_get_node_count(s), qaws_skeleton_get_arc_count(s), tmax, worst_dist, worst_axis, faces, area, missing, worst_round);
	TEST_ASSERT(worst_dist <= (SKT_F32 ? 1e-4 : 1e-8) * size && worst_axis <= (SKT_F32 ? 1e-4 : 1e-7) * size &&
		missing == 0 && fabs(faces - area) <= 2e-4 * area && worst_round <= (SKT_F32 ? 1e-3 : 1e-5) * area, msg);
	if (g_verbose_skt) printf("    %s\n", msg);
	qaws_skeleton_destroy(s);
	for (i = 0; i < ns; i++) qaws_curve_destroy(cs[i]);
}

/* the disc of radius 1 less the hole of radius 0.35, both moved by t, until
   the circles touch */
static double skt_holed_area(double t)
{
	return SKT_PI * ((1 - t) * (1 - t) - (0.35 + t) * (0.35 + t));
}

static void test_curved(void)
{
	skt_cshape sh;
	unsigned int n;
	double const pi = SKT_PI;
	/* a disc */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "disc"; sh.round = 1;
	skt_carc(&sh, &n, 0.2, -0.1, 1, 0, 2 * pi); sh.count[0] = n; sh.nloops = 1;
	skt_check_curved(&sh);
	/* a stadium */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "stadium"; sh.round = 1;
	skt_cline(&sh, &n, -1, -1, 1, -1); skt_carc(&sh, &n, 1, 0, 1, -pi / 2, pi / 2);
	skt_cline(&sh, &n, 1, 1, -1, 1); skt_carc(&sh, &n, -1, 0, 1, pi / 2, 3 * pi / 2);
	sh.count[0] = n; sh.nloops = 1;
	skt_check_curved(&sh);
	/* a rounded rectangle: smooth joints */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "rounded rectangle"; sh.round = 1;
	skt_cline(&sh, &n, 0.4, 0, 2.6, 0); skt_carc(&sh, &n, 2.6, 0.4, 0.4, -pi / 2, 0);
	skt_cline(&sh, &n, 3, 0.4, 3, 1.6); skt_carc(&sh, &n, 2.6, 1.6, 0.4, 0, pi / 2);
	skt_cline(&sh, &n, 2.6, 2, 0.4, 2); skt_carc(&sh, &n, 0.4, 1.6, 0.4, pi / 2, pi);
	skt_cline(&sh, &n, 0, 1.6, 0, 0.4); skt_carc(&sh, &n, 0.4, 0.4, 0.4, pi, 3 * pi / 2);
	sh.count[0] = n; sh.nloops = 1;
	skt_check_curved(&sh);
	/* a half disc */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "half disc"; sh.round = 1;
	skt_cline(&sh, &n, -1, 0, 1, 0); skt_carc(&sh, &n, 0, 0, 1, 0, pi); sh.count[0] = n; sh.nloops = 1;
	skt_check_curved(&sh);
	/* a lens: two convex arcs */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "lens"; sh.round = 1;
	skt_carc(&sh, &n, 0, 0.6, 1, pi + atan2(0.6, 0.8), 2 * pi - atan2(0.6, 0.8));
	skt_carc(&sh, &n, 0, -0.6, 1, atan2(0.6, 0.8), pi - atan2(0.6, 0.8));
	sh.count[0] = n; sh.nloops = 1;
	skt_check_curved(&sh);
	/* a disc with an off-centre hole: the circles squeeze */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "disc with an off-centre hole"; sh.round = 1;
	sh.exact = skt_holed_area;
	sh.exact_until = (1 - 0.35 - hypot(0.3, 0.1)) / 2;      /* the circles touch */
	skt_carc(&sh, &n, 0, 0, 1, 0, 2 * pi); sh.count[0] = 1;
	skt_carc(&sh, &n, 0.3, 0.1, 0.35, 2 * pi, 0); sh.count[1] = 1; sh.nloops = 2;
	skt_check_curved(&sh);
	/* a square with a round bite (concave arc, convex corners) */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "square with a round bite"; sh.round = 1;
	skt_cline(&sh, &n, 0, 0, 2, 0); skt_cline(&sh, &n, 2, 0, 2, 2); skt_cline(&sh, &n, 2, 2, 1.6, 2);
	skt_carc(&sh, &n, 1, 2, 0.6, 0, -pi); skt_cline(&sh, &n, 0.4, 2, 0, 2); skt_cline(&sh, &n, 0, 2, 0, 0);
	sh.count[0] = n; sh.nloops = 1;
	skt_check_curved(&sh);
	/* a square with a round bump: reflex corners, mitered then round */
	memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "square with a round bump, mitered";
	skt_cline(&sh, &n, 0, 0, 2, 0); skt_cline(&sh, &n, 2, 0, 2, 2); skt_cline(&sh, &n, 2, 2, 1.6, 2);
	skt_carc(&sh, &n, 1, 2, 0.6, 0, pi); skt_cline(&sh, &n, 0.4, 2, 0, 2); skt_cline(&sh, &n, 0, 2, 0, 0);
	sh.count[0] = n; sh.nloops = 1;
	skt_check_curved(&sh);
	sh.name = "square with a round bump, medial axis"; sh.sigma = pi; sh.round = 1;
	skt_check_curved(&sh);
	/* polygons' medial axes (sigma = pi): parabolic arcs at reflex corners */
	{
		static double const ell[12] = { 0, 0, 4, 0, 4, 1, 1, 1, 1, 3, 0, 3 };
		static double const comb[24] = { 0, 0, 7, 0, 7, 3, 6, 3, 6, 1, 5, 1, 5, 3, 3, 3, 3, 1, 1, 1, 1, 3, 0, 3 };
		unsigned int i;
		memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "L, medial axis"; sh.sigma = pi; sh.round = 1;
		for (i = 0; i < 6; i++) skt_cline(&sh, &n, ell[2 * i], ell[2 * i + 1], ell[2 * ((i + 1) % 6)], ell[2 * ((i + 1) % 6) + 1]);
		sh.count[0] = n; sh.nloops = 1;
		skt_check_curved(&sh);
		memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "comb, medial axis"; sh.sigma = pi; sh.round = 1;
		for (i = 0; i < 12; i++) skt_cline(&sh, &n, comb[2 * i], comb[2 * i + 1], comb[2 * ((i + 1) % 12)], comb[2 * ((i + 1) % 12) + 1]);
		sh.count[0] = n; sh.nloops = 1;
		skt_check_curved(&sh);
		memset(&sh, 0, sizeof(sh)); n = 0; sh.name = "7-point star, medial axis"; sh.sigma = pi; sh.round = 1;
		for (i = 0; i < 14; i++)
		{
			double a0 = 2 * pi * i / 14, a1 = 2 * pi * (i + 1) / 14, r0 = (i & 1) ? 0.45 : 1.0, r1 = ((i + 1) & 1) ? 0.45 : 1.0;
			skt_cline(&sh, &n, r0 * cos(a0), r0 * sin(a0), r1 * cos(a1), r1 * sin(a1));
		}
		sh.count[0] = n; sh.nloops = 1;
		skt_check_curved(&sh);
	}
}

/* ------------------------------------------------------------------------ */
/*  Free-form curves: through G1 arc splines                                */
/* ------------------------------------------------------------------------ */

static qaws_curve* skt_cubic(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
	qaws_scalar cp[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	cp[0] = (qaws_scalar)x0; cp[1] = (qaws_scalar)y0; cp[2] = (qaws_scalar)x1; cp[3] = (qaws_scalar)y1;
	cp[4] = (qaws_scalar)x2; cp[5] = (qaws_scalar)y2; cp[6] = (qaws_scalar)x3; cp[7] = (qaws_scalar)y3;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = 4;
	qaws_curve_create_bezier(&d, &c);
	return c;
}

/* every skeleton arc on its faces' supports; at sigma = pi each point at
   distance t from the true curves (within the fit); the faces tile the
   region of the spline */
static void skt_check_freeform(char const* name, qaws_curve* const* cs, unsigned int ncs, qaws_path_2d const* p, unsigned int np, double sigma)
{
	qaws_skeleton_desc d;
	qaws_skeleton* s = NULL;
	qaws_status st;
	unsigned int i, k, missing = 0;
	double area = 0, size, worst_dist = 0, worst_axis = 0, faces = 0, fit;
	char msg[400];
	{
		qaws_clip_result* u = NULL;
		if (qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, p, np, NULL, 0, &u) == QAWS_STATUS_OK)
			area = skt_result_area(u);
		qaws_clip_result_destroy(u);
	}
	size = sqrt(fabs(area));
	fit = 1e-4 * 2 * size;
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = np; d.fill_rule = QAWS_FILL_NON_ZERO;
	d.splice_angle = (qaws_scalar)sigma;
	d.curve_tolerance = (qaws_scalar)fit;
	st = qaws_skeleton_create(&d, &s);
	sprintf(msg, "%s: created (status %d)", name, (int)st);
	TEST_ASSERT(st == QAWS_STATUS_OK && s, msg);
	if (!s)
		return;
	for (i = 0; i < qaws_skeleton_get_arc_count(s); i++)
	{
		qaws_skeleton_arc a;
		qaws_skeleton_face fa, fb;
		qaws_skeleton_get_arc(s, i, &a);
		qaws_skeleton_get_face(s, a.face_a, &fa);
		qaws_skeleton_get_face(s, a.face_b, &fb);
		for (k = 0; k <= 8; k++)
		{
			qaws_vec2 q;
			qaws_scalar t;
			double e1, e2;
			qaws_skeleton_evaluate_arc(s, i, (qaws_scalar)(k / 8.0), &q, &t);
			e1 = fabs(skt_face_resid(&fa, q.x, q.y, t));
			e2 = fabs(skt_face_resid(&fb, q.x, q.y, t));
			if (e1 > worst_dist) worst_dist = e1;
			if (e2 > worst_dist) worst_dist = e2;
			if (sigma == SKT_PI && !a.smooth)
			{
				double e = fabs(skt_boundary_dist(cs, ncs, q.x, q.y) - (double)t);
				if (e > worst_axis) worst_axis = e;
			}
		}
	}
	for (i = 0; i < qaws_skeleton_get_face_count(s); i++)
		faces += skt_face_area(s, i, &missing);
	sprintf(msg, "%s (sigma %.3f, fit %.2g): %u faces, %u arcs, max time %.6f; supports' error %.2g, distance to the true curves %.2g, faces %.6f (region %.6f, %u gaps)",
		name, sigma, fit, qaws_skeleton_get_face_count(s), qaws_skeleton_get_arc_count(s), (double)qaws_skeleton_get_max_time(s),
		worst_dist, worst_axis, faces, area, missing);
	TEST_ASSERT(worst_dist <= (SKT_F32 ? 1e-4 : 1e-8) * size && worst_axis <= 2 * fit && missing == 0 &&
		fabs(faces - area) <= 1e-3 * area, msg);
	if (g_verbose_skt) printf("    %s\n", msg);
	qaws_skeleton_destroy(s);
}

static void test_freeform(void)
{
	qaws_curve* cs[8];
	qaws_curve const* cv[8];
	qaws_path_2d p[2];
	unsigned int i, k;
	double const pi = SKT_PI;
	for (k = 0; k < 2; k++)
	{
		double sigma = k ? pi : 0;
		/* an ellipse: one exact rational NURBS */
		qaws_curve_create_ellipse_2d(0.1, -0.2, 1.5, 0.8, 0.3, &cs[0]);
		cv[0] = cs[0];
		p[0].curves = cv; p[0].curve_count = 1; p[0].closed = 1;
		skt_check_freeform("ellipse (NURBS)", cs, 1, p, 1, sigma);
		/* an ellipse with an off-centre round hole */
		{
			qaws_arc_segment g;
			qaws_arc_desc d;
			memset(&g, 0, sizeof(g));
			g.center[0] = (qaws_scalar)0.5; g.center[1] = (qaws_scalar)-0.1; g.radius = (qaws_scalar)0.3;
			g.angle_start = (qaws_scalar)(2 * pi); g.angle_end = 0;
			memset(&d, 0, sizeof(d));
			d.dimension = QAWS_DIMENSION_2D; d.segments = &g; d.segment_count = 1;
			qaws_curve_create_arc(&d, &cs[1]);
			cv[1] = cs[1];
			p[1].curves = &cv[1]; p[1].curve_count = 1; p[1].closed = 1;
			skt_check_freeform("ellipse with a round hole", cs, 2, p, 2, sigma);
			qaws_curve_destroy(cs[1]);
		}
		qaws_curve_destroy(cs[0]);
		/* a blob of four cubic Beziers, with a dent */
		cs[0] = skt_cubic(0, 0, 1, -0.9, 2.4, -0.5, 2.2, 0.6);
		cs[1] = skt_cubic(2.2, 0.6, 2, 1.6, 1.4, 1.1, 1, 1);
		cs[2] = skt_cubic(1, 1, 0.6, 0.9, 0.2, 1.8, -0.5, 1.2);
		cs[3] = skt_cubic(-0.5, 1.2, -1.2, 0.6, -0.6, 0.2, 0, 0);
		for (i = 0; i < 4; i++) cv[i] = cs[i];
		p[0].curves = cv; p[0].curve_count = 4; p[0].closed = 1;
		skt_check_freeform("blob of four cubics", cs, 4, p, 1, sigma);
		for (i = 0; i < 4; i++) qaws_curve_destroy(cs[i]);
	}
}

/* ------------------------------------------------------------------------ */
/*  Random curved regions                                                   */
/* ------------------------------------------------------------------------ */

static unsigned long long g_skt_seed;

static double skt_rand(void)
{
	g_skt_seed = g_skt_seed * 6364136223846793005ull + 1442695040888963407ull;
	return (double)(g_skt_seed >> 11) * (1.0 / 9007199254740992.0);
}

static qaws_curve* skt_line_curve(double x0, double y0, double x1, double y1)
{
	skt_seg g;
	memset(&g, 0, sizeof(g));
	g.arc = 0; g.a[0] = x0; g.a[1] = y0; g.a[2] = x1; g.a[3] = y1;
	return skt_seg_curve(&g);
}

static qaws_curve* skt_arc_curve(double cx, double cy, double r, double a0, double a1)
{
	skt_seg g;
	g.arc = 1; g.a[0] = cx; g.a[1] = cy; g.a[2] = r; g.a[3] = a0; g.a[4] = a1;
	return skt_seg_curve(&g);
}

/* the arc from a to b with sagitta h, left of a to b when h > 0 */
static qaws_curve* skt_bulge(double ax, double ay, double bx, double by, double h)
{
	double mx = 0.5 * (ax + bx), my = 0.5 * (ay + by), dx = bx - ax, dy = by - ay, l = hypot(dx, dy), c = 0.5 * l;
	double r = (h * h + c * c) / (2 * fabs(h)), nx = -dy / l, ny = dx / l, sg = h > 0 ? 1 : -1;
	double cx = mx - sg * (r - fabs(h)) * nx, cy = my - sg * (r - fabs(h)) * ny;
	double a0 = atan2(ay - cy, ax - cx), a1 = atan2(by - cy, bx - cx), sw = a1 - a0;
	if (h > 0) { while (sw > 0) sw -= 2 * SKT_PI; }
	else { while (sw < 0) sw += 2 * SKT_PI; }
	return skt_arc_curve(cx, cy, r, a0, a0 + sw);
}

/* star-shaped polygons of 5 to 18 vertices whose sides are lines or arcs
   bulging either way, two in five with a round hole; both sigma. Each
   region that builds tiles with its faces and its arcs lie on their
   supports (to half the digits: events near tangencies are double roots);
   at sigma = pi the wavefront is the round offset */
static void test_random_curved(void)
{
	unsigned int const shapes = 40;
	unsigned int c, built = 0, runs = 0, bad = 0;
	double worst_dist = 0, worst_round = 0, worst_tile = 0;
	char msg[600], fails[300] = "";
	if (SKT_F32)
	{
		/* float vertices: the regions' corners rounded off their curves'
		   meetings; double builds only */
		printf("    random curved regions: double builds only\n");
		return;
	}
	g_skt_seed = 4242;
	for (c = 0; c < shapes; c++)
	{
		qaws_curve* cs[40];
		qaws_curve const* cv[40];
		qaws_path_2d p[2];
		double xy[2 * 20], area = 0, size;
		unsigned int n = 5 + (unsigned int)(skt_rand() * 14), i, np = 1, nc = 0, pass;
		for (i = 0; i < n; i++)
		{
			double a = 2 * SKT_PI * (i + 0.3 * skt_rand()) / n, rr = 1 + 0.6 * skt_rand();
			xy[2 * i] = rr * cos(a); xy[2 * i + 1] = rr * sin(a);
		}
		for (i = 0; i < n; i++)
		{
			double ax = xy[2 * i], ay = xy[2 * i + 1], bx = xy[2 * ((i + 1) % n)], by = xy[2 * ((i + 1) % n) + 1];
			double l = hypot(bx - ax, by - ay);
			cs[nc++] = skt_rand() < 0.35 ? skt_line_curve(ax, ay, bx, by) : skt_bulge(ax, ay, bx, by, (skt_rand() - 0.5) * 0.5 * l);
		}
		p[0].curves = cv; p[0].curve_count = nc; p[0].closed = 1;
		if (skt_rand() < 0.4)
		{
			double hx = 0.25 * (skt_rand() - 0.5), hy = 0.25 * (skt_rand() - 0.5);
			cs[nc++] = skt_arc_curve(hx, hy, 0.15 + 0.2 * skt_rand(), 2 * SKT_PI, 0);
			np = 2;
			p[1].curves = cv + nc - 1; p[1].curve_count = 1; p[1].closed = 1;
		}
		for (i = 0; i < nc; i++) cv[i] = cs[i];
		{
			qaws_clip_result* u = NULL;
			if (qaws_clip_boolean(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, p, np, NULL, 0, &u) == QAWS_STATUS_OK)
				area = skt_result_area(u);
			qaws_clip_result_destroy(u);
		}
		size = sqrt(fabs(area));
		for (pass = 0; pass < 2; pass++)
		{
			qaws_skeleton_desc d;
			qaws_skeleton* s = NULL;
			unsigned int missing = 0, k;
			double dist = 0, round = 0, faces = 0, tmax;
			int ok;
			memset(&d, 0, sizeof(d));
			d.paths = p; d.path_count = np; d.fill_rule = QAWS_FILL_NON_ZERO;
			d.splice_angle = (qaws_scalar)(pass ? SKT_PI : 0);
			runs++;
			if (qaws_skeleton_create(&d, &s) != QAWS_STATUS_OK)
			{
				if (strlen(fails) < sizeof(fails) - 40)
					sprintf(fails + strlen(fails), " %u%s", c, pass ? " (pi)" : "");
				continue;
			}
			built++;
			tmax = (double)qaws_skeleton_get_max_time(s);
			for (i = 0; i < qaws_skeleton_get_arc_count(s); i++)
			{
				qaws_skeleton_arc a;
				qaws_skeleton_face fa, fb;
				qaws_skeleton_get_arc(s, i, &a);
				qaws_skeleton_get_face(s, a.face_a, &fa);
				qaws_skeleton_get_face(s, a.face_b, &fb);
				for (k = 0; k <= 4; k++)
				{
					qaws_vec2 q;
					qaws_scalar t;
					double e1, e2;
					qaws_skeleton_evaluate_arc(s, i, (qaws_scalar)(k / 4.0), &q, &t);
					e1 = fabs(skt_face_resid(&fa, q.x, q.y, t));
					e2 = fabs(skt_face_resid(&fb, q.x, q.y, t));
					if (e1 > dist) dist = e1;
					if (e2 > dist) dist = e2;
				}
			}
			for (i = 0; i < qaws_skeleton_get_face_count(s); i++)
				faces += skt_face_area(s, i, &missing);
			if (pass)
				for (k = 1; k <= 4; k++)
				{
					double t = tmax * k / 5;
					qaws_clip_result* r1 = NULL;
					qaws_clip_result* r2 = NULL;
					qaws_skeleton_offset(s, (qaws_scalar)t, &r1);
					qaws_offset_paths(p, np, (qaws_scalar)-t, QAWS_JOIN_ROUND, QAWS_END_POLYGON, &r2);
					if (fabs(skt_result_area(r1) - skt_result_area(r2)) > round)
						round = fabs(skt_result_area(r1) - skt_result_area(r2));
					qaws_clip_result_destroy(r1);
					qaws_clip_result_destroy(r2);
				}
			ok = dist <= (SKT_F32 ? 1e-4 : 1e-5) * size && missing == 0 && fabs(faces - area) <= 2e-4 * area &&
				round <= (SKT_F32 ? 1e-3 : 1e-6) * area;
			if (!ok)
			{
				bad++;
				if (strlen(fails) < sizeof(fails) - 40)
					sprintf(fails + strlen(fails), " %u%s wrong", c, pass ? " (pi)" : "");
			}
			if (dist / size > worst_dist) worst_dist = dist / size;
			if (round / area > worst_round) worst_round = round / area;
			if (fabs(faces - area) / area > worst_tile) worst_tile = fabs(faces - area) / area;
			if (g_verbose_skt) printf("    random %u (sigma %s): %u faces, supports' error %.2g, faces %.6f of %.6f (%u gaps), round offset error %.2g\n",
				c, pass ? "pi" : "2 pi", qaws_skeleton_get_face_count(s), dist, faces, area, missing, round);
			qaws_skeleton_destroy(s);
		}
		for (i = 0; i < nc; i++) qaws_curve_destroy(cs[i]);
	}
	sprintf(msg, "random curved regions: %u of %u built, %u off; worst supports' error %.2g, tiling %.2g, round offset %.2g of the size or area; failed:%s",
		built, runs, bad, worst_dist, worst_tile, worst_round, fails[0] ? fails : " none");
	TEST_ASSERT(built + 1 >= runs && bad == 0, msg);
	if (g_verbose_skt) printf("    %s\n", msg);
}

int test_99_skeleton_main(void)
{
	g_pass = 0;
	g_fail = 0;
	g_verbose_skt = getenv("QAWS_SKT_VERBOSE") != NULL;
	printf("Test 99: straight skeletons\n");
	test_closed_forms();
	test_shapes();
	test_curved();
	test_freeform();
	test_random_curved();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
