/*
 * Test 91: triangulating regions of curves (qaws_triangulate)
 *
 *   - a square: 2 triangles, the area exact
 *   - a NURBS disc: boundary vertices on the circle, every free edge
 *     locally Delaunay, all triangles counter-clockwise, the area the
 *     inscribed polygon's
 *   - an annulus: no triangle in the hole, Euler characteristic 0
 *   - collinear and repeated input points, a bow-tie under even-odd and
 *     non-zero, two touching squares
 *   - refinement: minimum angle 30 degrees reached, maximum area kept, a
 *     size field kept; refined boundary vertices still on the curves
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

#define TRT_PI 3.14159265358979323846
#define TRT_F32 (sizeof(qaws_scalar) == 4)

typedef struct trt_stats
{
	unsigned int nv, nt, nb, cw, edges;
	double area, min_angle, max_area, max_edge, on_curve;
	unsigned int not_delaunay;
} trt_stats;

static double trt_cross(qaws_vec2 a, qaws_vec2 b, qaws_vec2 c)
{
	return ((double)b.x - a.x) * ((double)c.y - a.y) - ((double)b.y - a.y) * ((double)c.x - a.x);
}

/* in-circle in double (the test's inputs are tame) */
static double trt_incircle(qaws_vec2 a, qaws_vec2 b, qaws_vec2 c, qaws_vec2 d)
{
	double adx = (double)a.x - d.x, ady = (double)a.y - d.y, bdx = (double)b.x - d.x, bdy = (double)b.y - d.y, cdx = (double)c.x - d.x, cdy = (double)c.y - d.y;
	return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy) + (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
}

static void trt_measure(qaws_mesh_2d const* m, trt_stats* s)
{
	unsigned int i, k, j, nt = qaws_mesh_2d_get_triangle_count(m), nb = qaws_mesh_2d_get_boundary_edge_count(m);
	unsigned int* tri = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	memset(s, 0, sizeof(*s));
	s->nv = qaws_mesh_2d_get_vertex_count(m);
	s->nt = nt;
	s->nb = nb;
	s->min_angle = 180;
	for (i = 0; i < nt; i++)
	{
		qaws_vec2 p[3];
		double a, l[3];
		qaws_mesh_2d_get_triangle(m, i, &tri[3 * i]);
		for (k = 0; k < 3; k++) p[k] = qaws_mesh_2d_get_vertex(m, tri[3 * i + k]);
		a = 0.5 * trt_cross(p[0], p[1], p[2]);
		if (a <= 0) s->cw++;
		s->area += a;
		if (a > s->max_area) s->max_area = a;
		for (k = 0; k < 3; k++)
		{
			qaws_vec2 u = p[(k + 1) % 3], w = p[(k + 2) % 3];
			l[k] = hypot((double)w.x - u.x, (double)w.y - u.y);
			if (l[k] > s->max_edge) s->max_edge = l[k];
		}
		for (k = 0; k < 3; k++)
		{
			/* the angle at p[k] by the law of cosines */
			double b = l[(k + 1) % 3], c = l[(k + 2) % 3], cosv = (b * b + c * c - l[k] * l[k]) / (2 * b * c), ang;
			if (cosv > 1) cosv = 1;
			if (cosv < -1) cosv = -1;
			ang = acos(cosv) * 180 / TRT_PI;
			if (ang < s->min_angle) s->min_angle = ang;
		}
	}
	/* free edges locally Delaunay: the opposite vertex of the neighbour
	   outside the circumcircle (quadratic search: test meshes are small) */
	if (nt < 4000)
		for (i = 0; i < nt; i++)
			for (k = 0; k < 3; k++)
			{
				unsigned int a = tri[3 * i + (k + 1) % 3], b = tri[3 * i + (k + 2) % 3], boundary = 0;
				for (j = 0; j < nb && !boundary; j++)
				{
					unsigned int e[2];
					qaws_mesh_2d_get_boundary_edge(m, j, e);
					boundary = (e[0] == a && e[1] == b) || (e[0] == b && e[1] == a);
				}
				if (boundary) continue;
				for (j = 0; j < nt; j++)
				{
					unsigned int q;
					for (q = 0; q < 3; q++)
						if (tri[3 * j + (q + 1) % 3] == b && tri[3 * j + (q + 2) % 3] == a)
						{
							qaws_vec2 P0 = qaws_mesh_2d_get_vertex(m, tri[3 * i]), P1 = qaws_mesh_2d_get_vertex(m, tri[3 * i + 1]), P2 = qaws_mesh_2d_get_vertex(m, tri[3 * i + 2]);
							qaws_vec2 D = qaws_mesh_2d_get_vertex(m, tri[3 * j + q]);
							double scale = 0;
							{
								double ex = fabs((double)P0.x) + fabs((double)P1.x) + fabs((double)P2.x) + fabs((double)P0.y) + fabs((double)P1.y) + fabs((double)P2.y) + 1;
								scale = ex * ex * ex * ex;
							}
							if (trt_incircle(P0, P1, P2, D) > (TRT_F32 ? 1e-5 : 1e-11) * scale) s->not_delaunay++;
						}
				}
				s->edges++;
			}
	/* boundary vertices on their curves */
	for (i = 0; i < s->nv; i++)
	{
		unsigned int c;
		qaws_scalar t;
		if (qaws_mesh_2d_get_vertex_source(m, i, &c, &t))
		{
			qaws_eval_result_2d e;
			qaws_vec2 p = qaws_mesh_2d_get_vertex(m, i);
			double d;
			qaws_curve_evaluate_2d(&e, qaws_mesh_2d_get_boundary_curve(m, c), t, QAWS_EVAL_FLAG_POSITION);
			d = hypot((double)e.position.x - p.x, (double)e.position.y - p.y);
			if (d > s->on_curve) s->on_curve = d;
		}
	}
	free(tri);
}

static qaws_mesh_2d* trt_run(qaws_curve const** cs, unsigned int n, qaws_fill_rule fr, double tol, double min_angle, double max_area,
	qaws_size_fn size_fn)
{
	qaws_path_2d p[8];
	qaws_triangulate_desc d;
	qaws_mesh_2d* m = NULL;
	unsigned int i;
	for (i = 0; i < n; i++) { p[i].curves = &cs[i]; p[i].curve_count = 1; p[i].closed = 1; }
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = n; d.fill_rule = fr;
	d.tolerance = (qaws_scalar)tol; d.min_angle = (qaws_scalar)min_angle; d.max_area = (qaws_scalar)max_area;
	d.size_fn = size_fn;
	qaws_triangulate_region_2d(&d, &m);
	return m;
}

static qaws_curve* trt_poly(double const* xy, unsigned int n)
{
	qaws_scalar p[64];
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 2 * n; i++) p[i] = (qaws_scalar)xy[i];
	qaws_curve_create_polyline_2d(p, n, 1, &c);
	return c;
}

static void test_square(void)
{
	static double const sq[8] = { 0, 0, 2, 0, 2, 2, 0, 2 };
	qaws_curve* c = trt_poly(sq, 4);
	qaws_curve const* cc = c;
	qaws_mesh_2d* m = trt_run(&cc, 1, QAWS_FILL_NON_ZERO, 0, 0, 0, NULL);
	trt_stats s;
	char msg[200];
	TEST_ASSERT(m != NULL, "square: triangulated");
	if (m)
	{
		trt_measure(m, &s);
		sprintf(msg, "square: %u vertices, %u triangles, area %.12g, %u boundary edges", s.nv, s.nt, s.area, s.nb);
		TEST_ASSERT(s.nv == 4 && s.nt == 2 && fabs(s.area - 4) < 1e-12 && s.nb == 4 && s.cw == 0, msg);
	}
	qaws_mesh_2d_destroy(m);
	qaws_curve_destroy(c);
}

static void test_disc(void)
{
	qaws_curve* c = NULL;
	qaws_curve const* cc;
	qaws_mesh_2d* m;
	trt_stats s;
	char msg[300];
	qaws_curve_create_ellipse_2d(0, 0, 1, 1, 0, &c);
	cc = c;
	m = trt_run(&cc, 1, QAWS_FILL_NON_ZERO, 1e-3, 0, 0, NULL);
	TEST_ASSERT(m != NULL, "disc: triangulated");
	if (m)
	{
		double r = 0;
		unsigned int i;
		trt_measure(m, &s);
		for (i = 0; i < s.nv; i++)
		{
			qaws_vec2 p = qaws_mesh_2d_get_vertex(m, i);
			double d = fabs(hypot((double)p.x, (double)p.y) - 1);
			if (d > r) r = d;
		}
		/* the inscribed polygon's area: chord sag 1e-3 gives about pi - n sin(2pi/n)/2 */
		sprintf(msg, "disc: %u boundary vertices on the circle (%.2e), %u of %u free edges not Delaunay, %u clockwise, area %.6f (pi - %.2e)",
			s.nb, r, s.not_delaunay, s.edges, s.cw, s.area, TRT_PI - s.area);
		TEST_ASSERT(s.nb == s.nv && r < (TRT_F32 ? 1e-6 : 1e-12) && s.not_delaunay == 0 && s.cw == 0
			&& TRT_PI - s.area > 0 && TRT_PI - s.area < 2 * TRT_PI * 1e-3, msg);
		sprintf(msg, "disc: boundary vertices on their curves at their parameters (%.2e)", s.on_curve);
		TEST_ASSERT(s.on_curve < (TRT_F32 ? 1e-6 : 1e-12), msg);
	}
	qaws_mesh_2d_destroy(m);
	qaws_curve_destroy(c);
}

static void test_annulus(void)
{
	qaws_curve *o = NULL, *h = NULL, *hr = NULL;
	qaws_curve const* cs[2];
	qaws_mesh_2d* m;
	trt_stats s;
	char msg[300];
	qaws_curve_create_ellipse_2d(0, 0, 2, 2, 0, &o);
	qaws_curve_create_ellipse_2d((qaws_scalar)0.3, 0, 1, 1, 0, &h);
	qaws_curve_reverse(&hr, h);
	cs[0] = o; cs[1] = hr;
	m = trt_run(cs, 2, QAWS_FILL_NON_ZERO, 2e-3, 0, 0, NULL);
	TEST_ASSERT(m != NULL, "annulus: triangulated");
	if (m)
	{
		unsigned int i, inside = 0, nt;
		trt_measure(m, &s);
		nt = s.nt;
		for (i = 0; i < nt; i++)
		{
			unsigned int v[3];
			qaws_vec2 a, b, c;
			double cx, cy;
			qaws_mesh_2d_get_triangle(m, i, v);
			a = qaws_mesh_2d_get_vertex(m, v[0]); b = qaws_mesh_2d_get_vertex(m, v[1]); c = qaws_mesh_2d_get_vertex(m, v[2]);
			cx = ((double)a.x + b.x + c.x) / 3 - 0.3; cy = ((double)a.y + b.y + c.y) / 3;
			if (hypot(cx, cy) < 0.99) inside++;
		}
		/* V - E + F = 0 for one hole: E = (3F + boundary) / 2 */
		sprintf(msg, "annulus: %u triangles in the hole, Euler characteristic %d, area %.5f (3 pi = %.5f), %u not Delaunay",
			inside, (int)s.nv - (int)((3 * s.nt + s.nb) / 2) + (int)s.nt, s.area, 3 * TRT_PI, s.not_delaunay);
		TEST_ASSERT(inside == 0 && (int)s.nv - (int)((3 * s.nt + s.nb) / 2) + (int)s.nt == 0 && fabs(s.area - 3 * TRT_PI) < 0.02 && s.not_delaunay == 0, msg);
	}
	qaws_mesh_2d_destroy(m);
	qaws_curve_destroy(o); qaws_curve_destroy(h); qaws_curve_destroy(hr);
}

static void test_degenerate(void)
{
	char msg[300];
	/* collinear and repeated points on the boundary */
	{
		static double const q[18] = { 0, 0, 1, 0, 2, 0, 2, 0, 3, 0, 3, 3, 1.5, 3, 0, 3, 0, 1.5 };
		qaws_curve* c = trt_poly(q, 9);
		qaws_curve const* cc = c;
		qaws_mesh_2d* m = trt_run(&cc, 1, QAWS_FILL_NON_ZERO, 0, 0, 0, NULL);
		trt_stats s;
		if (m) trt_measure(m, &s);
		sprintf(msg, "collinear and repeated points: area %.12g of 9, %u clockwise", m ? s.area : 0, m ? s.cw : 0);
		TEST_ASSERT(m && fabs(s.area - 9) < 1e-9 && s.cw == 0, msg);
		qaws_mesh_2d_destroy(m);
		qaws_curve_destroy(c);
	}
	/* a bow-tie: two triangles of area 1 each under both rules */
	{
		static double const b[8] = { 0, 0, 2, 2, 2, 0, 0, 2 };
		qaws_curve* c = trt_poly(b, 4);
		qaws_curve const* cc = c;
		unsigned int r;
		for (r = 0; r < 2; r++)
		{
			qaws_mesh_2d* m = trt_run(&cc, 1, r ? QAWS_FILL_NON_ZERO : QAWS_FILL_EVEN_ODD, 0, 0, 0, NULL);
			trt_stats s;
			if (m) trt_measure(m, &s);
			sprintf(msg, "bow-tie (%s): area %.12g of 2, %u triangles", r ? "non-zero" : "even-odd", m ? s.area : 0, m ? s.nt : 0);
			TEST_ASSERT(m && fabs(s.area - 2) < 1e-9 && s.cw == 0, msg);
			qaws_mesh_2d_destroy(m);
		}
		qaws_curve_destroy(c);
	}
	/* two squares touching at a corner */
	{
		static double const a[8] = { 0, 0, 1, 0, 1, 1, 0, 1 }, b2[8] = { 1, 1, 2, 1, 2, 2, 1, 2 };
		qaws_curve *c1 = trt_poly(a, 4), *c2 = trt_poly(b2, 4);
		qaws_curve const* cs[2];
		qaws_mesh_2d* m;
		trt_stats s;
		cs[0] = c1; cs[1] = c2;
		m = trt_run(cs, 2, QAWS_FILL_NON_ZERO, 0, 0, 0, NULL);
		if (m) trt_measure(m, &s);
		sprintf(msg, "two squares touching at a corner: area %.12g of 2, %u triangles", m ? s.area : 0, m ? s.nt : 0);
		TEST_ASSERT(m && fabs(s.area - 2) < 1e-9 && s.nt == 4 && s.cw == 0, msg);
		qaws_mesh_2d_destroy(m);
		qaws_curve_destroy(c1); qaws_curve_destroy(c2);
	}
}

static qaws_scalar trt_size(void* user, qaws_vec2 p)
{
	(void)user;
	/* fine near the left, coarse to the right */
	return (qaws_scalar)(0.05 + 0.1 * ((double)p.x + 1));
}

static void test_refine(void)
{
	qaws_curve *o = NULL, *h = NULL, *hr = NULL;
	qaws_curve const* cs[2];
	char msg[300];
	qaws_curve_create_ellipse_2d(0, 0, 1, (qaws_scalar)0.7, 0, &o);
	qaws_curve_create_ellipse_2d((qaws_scalar)0.2, 0, (qaws_scalar)0.3, (qaws_scalar)0.2, 0, &h);
	qaws_curve_reverse(&hr, h);
	cs[0] = o; cs[1] = hr;
	{
		qaws_mesh_2d* m = trt_run(cs, 2, QAWS_FILL_NON_ZERO, 2e-3, 30, 0, NULL);
		trt_stats s;
		if (m) trt_measure(m, &s);
		sprintf(msg, "refined to 30 degrees: smallest angle %.3f, %u vertices, %u triangles, boundary on the curves (%.2e), %u not Delaunay",
			m ? s.min_angle : 0, m ? s.nv : 0, m ? s.nt : 0, m ? s.on_curve : 0, m ? s.not_delaunay : 0);
		TEST_ASSERT(m && s.nt > 0 && s.min_angle >= 30 - 1e-6 && s.cw == 0 && s.on_curve < (TRT_F32 ? 1e-6 : 1e-12) && s.not_delaunay == 0, msg);
		qaws_mesh_2d_destroy(m);
	}
	{
		qaws_mesh_2d* m = trt_run(cs, 2, QAWS_FILL_NON_ZERO, 2e-3, 20, 0.002, NULL);
		trt_stats s;
		if (m) trt_measure(m, &s);
		sprintf(msg, "maximum area 0.002: largest %.5f, smallest angle %.2f, %u triangles", m ? s.max_area : 0, m ? s.min_angle : 0, m ? s.nt : 0);
		TEST_ASSERT(m && s.nt > 0 && s.max_area <= 0.002 * (1 + 1e-9) && s.min_angle >= 20 - 1e-6 && s.cw == 0, msg);
		qaws_mesh_2d_destroy(m);
	}
	{
		qaws_mesh_2d* m = trt_run(cs, 2, QAWS_FILL_NON_ZERO, 2e-3, 0, 0, trt_size);
		trt_stats s;
		unsigned int i, over = 0;
		if (m)
		{
			trt_measure(m, &s);
			for (i = 0; i < s.nt; i++)
			{
				unsigned int v[3], k;
				qaws_vec2 p[3], c;
				double lmax = 0;
				qaws_mesh_2d_get_triangle(m, i, v);
				for (k = 0; k < 3; k++) p[k] = qaws_mesh_2d_get_vertex(m, v[k]);
				for (k = 0; k < 3; k++)
				{
					double l = hypot((double)p[(k + 1) % 3].x - p[k].x, (double)p[(k + 1) % 3].y - p[k].y);
					if (l > lmax) lmax = l;
				}
				c.x = (qaws_scalar)(((double)p[0].x + p[1].x + p[2].x) / 3);
				c.y = (qaws_scalar)(((double)p[0].y + p[1].y + p[2].y) / 3);
				if (lmax > (double)trt_size(NULL, c) * (1 + 1e-9)) over++;
			}
		}
		sprintf(msg, "size field: %u of %u triangles with an edge over the field", over, m ? s.nt : 0);
		TEST_ASSERT(m && s.nt > 0 && over == 0 && s.cw == 0, msg);
		qaws_mesh_2d_destroy(m);
	}
	qaws_curve_destroy(o); qaws_curve_destroy(h); qaws_curve_destroy(hr);
}

/* mean and largest edge length (each interior edge counted twice: the mean
   is what matters) */
static void trt_edges(qaws_mesh_2d const* m, double* mean, double* lmax)
{
	unsigned int i, k, nt = qaws_mesh_2d_get_triangle_count(m), n = 0;
	double s = 0;
	*lmax = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int v[3];
		qaws_mesh_2d_get_triangle(m, i, v);
		for (k = 0; k < 3; k++)
		{
			qaws_vec2 a = qaws_mesh_2d_get_vertex(m, v[k]), b = qaws_mesh_2d_get_vertex(m, v[(k + 1) % 3]);
			double l = hypot((double)b.x - a.x, (double)b.y - a.y);
			s += l;
			n++;
			if (l > *lmax) *lmax = l;
		}
	}
	*mean = n ? s / n : 0;
}

/* fine on the left of the ring, coarse on the right (positive on it) */
static qaws_scalar trt_ring_size(void* user, qaws_vec2 p)
{
	(void)user;
	return (qaws_scalar)(0.06 + 0.05 * ((double)p.x + 2));
}

static void test_remesh(void)
{
	char msg[300];
	/* a disc remeshed to edges of 0.1 */
	{
		qaws_curve* c = NULL;
		qaws_curve const* cc;
		qaws_mesh_2d* m;
		qaws_curve_create_ellipse_2d(0, 0, 1, 1, 0, &c);
		cc = c;
		m = trt_run(&cc, 1, QAWS_FILL_NON_ZERO, 1e-3, 0, 0, NULL);
		if (m && qaws_mesh_2d_remesh(m, (qaws_scalar)0.1, NULL, NULL, 6) == QAWS_STATUS_OK)
		{
			trt_stats s;
			double mean, lmax;
			trt_measure(m, &s);
			trt_edges(m, &mean, &lmax);
			sprintf(msg, "disc remeshed to 0.1: mean edge %.4f, longest %.4f, smallest angle %.1f, %u triangles, %u clockwise, area pi - %.2e, boundary on the circle (%.2e)",
				mean, lmax, s.min_angle, s.nt, s.cw, TRT_PI - s.area, s.on_curve);
			TEST_ASSERT(s.nt > 100 && mean > 0.08 && mean < 0.12 && lmax < 0.1 * 4.0 / 3.0 * 1.05 && s.min_angle > 25 && s.cw == 0
				&& fabs(TRT_PI - s.area) < 0.02 && s.on_curve < (TRT_F32 ? 1e-6 : 1e-12), msg);
		}
		else
			TEST_ASSERT(0, "disc remeshed");
		qaws_mesh_2d_destroy(m);
		qaws_curve_destroy(c);
	}
	/* an annulus under a size field: one hole kept */
	{
		qaws_curve *o = NULL, *h = NULL, *hr = NULL;
		qaws_curve const* cs[2];
		qaws_mesh_2d* m;
		qaws_curve_create_ellipse_2d(0, 0, 2, 2, 0, &o);
		qaws_curve_create_ellipse_2d((qaws_scalar)0.3, 0, 1, 1, 0, &h);
		qaws_curve_reverse(&hr, h);
		cs[0] = o; cs[1] = hr;
		m = trt_run(cs, 2, QAWS_FILL_NON_ZERO, 2e-3, 0, 0, NULL);
		if (m && qaws_mesh_2d_remesh(m, 0, trt_ring_size, NULL, 6) == QAWS_STATUS_OK)
		{
			trt_stats s;
			trt_measure(m, &s);
			sprintf(msg, "annulus remeshed under a size field: Euler characteristic %d, %u triangles, smallest angle %.1f, %u clockwise, area %.4f (3 pi %.4f), boundary on the curves (%.2e)",
				(int)s.nv - (int)((3 * s.nt + s.nb) / 2) + (int)s.nt, s.nt, s.min_angle, s.cw, s.area, 3 * TRT_PI, s.on_curve);
			TEST_ASSERT((int)s.nv - (int)((3 * s.nt + s.nb) / 2) + (int)s.nt == 0 && s.min_angle > 20 && s.cw == 0
				&& fabs(s.area - 3 * TRT_PI) < 0.05 && s.on_curve < (TRT_F32 ? 1e-6 : 1e-12), msg);
		}
		else
			TEST_ASSERT(0, "annulus remeshed");
		qaws_mesh_2d_destroy(m);
		qaws_curve_destroy(o); qaws_curve_destroy(h); qaws_curve_destroy(hr);
	}
}

/* faces straddling y = 1 (a vertex strictly above and one strictly below),
   over triangles and quads */
static unsigned int trt_straddle(qaws_mesh_2d const* m)
{
	unsigned int i, k, bad = 0;
	double const e = TRT_F32 ? 1e-5 : 1e-12;
	for (i = 0; i < qaws_mesh_2d_get_triangle_count(m); i++)
	{
		unsigned int v[3], up = 0, dn = 0;
		qaws_mesh_2d_get_triangle(m, i, v);
		for (k = 0; k < 3; k++)
		{
			double y = (double)qaws_mesh_2d_get_vertex(m, v[k]).y;
			up += y > 1 + e; dn += y < 1 - e;
		}
		bad += up && dn;
	}
	for (i = 0; i < qaws_mesh_2d_get_quad_count(m); i++)
	{
		unsigned int v[4], up = 0, dn = 0;
		qaws_mesh_2d_get_quad(m, i, v);
		for (k = 0; k < 4; k++)
		{
			double y = (double)qaws_mesh_2d_get_vertex(m, v[k]).y;
			up += y > 1 + e; dn += y < 1 - e;
		}
		bad += up && dn;
	}
	return bad;
}

/* the constraint edges on y = 1: their total length, and the worst vertex off the line */
static void trt_line(qaws_mesh_2d const* m, double* len, double* off)
{
	unsigned int i;
	*len = 0;
	*off = 0;
	for (i = 0; i < qaws_mesh_2d_get_constraint_edge_count(m); i++)
	{
		unsigned int e[2];
		qaws_vec2 a, b;
		qaws_mesh_2d_get_constraint_edge(m, i, e);
		a = qaws_mesh_2d_get_vertex(m, e[0]);
		b = qaws_mesh_2d_get_vertex(m, e[1]);
		if (fabs((double)a.y - 1) < 0.25 && fabs((double)b.y - 1) < 0.25)
		{
			*len += hypot((double)b.x - a.x, (double)b.y - a.y);
			if (fabs((double)a.y - 1) > *off) *off = fabs((double)a.y - 1);
			if (fabs((double)b.y - 1) > *off) *off = fabs((double)b.y - 1);
		}
	}
}

/* quads: convex and counter-clockwise, their total area */
static unsigned int trt_quads(qaws_mesh_2d const* m, double* area)
{
	unsigned int i, k, bad = 0;
	*area = 0;
	for (i = 0; i < qaws_mesh_2d_get_quad_count(m); i++)
	{
		unsigned int v[4];
		qaws_vec2 p[4];
		qaws_mesh_2d_get_quad(m, i, v);
		for (k = 0; k < 4; k++) p[k] = qaws_mesh_2d_get_vertex(m, v[k]);
		for (k = 0; k < 4; k++)
			if (!(trt_cross(p[k], p[(k + 1) % 4], p[(k + 2) % 4]) > 0)) { bad++; break; }
		*area += 0.5 * (trt_cross(p[0], p[1], p[2]) + trt_cross(p[0], p[2], p[3]));
	}
	return bad;
}

static void test_constraints_and_quads(void)
{
	static double const sq[8] = { 0, 0, 4, 0, 4, 4, 0, 4 };
	qaws_curve *box = trt_poly(sq, 4), *line = NULL, *circ = NULL;
	qaws_curve const *cb = box, *cl[2];
	qaws_path_2d p, cp[2];
	qaws_triangulate_desc d;
	char msg[400];
	unsigned int pass;
	{
		qaws_scalar l[4] = { 0, 1, 4, 1 };
		qaws_curve_create_polyline_2d(l, 2, 0, &line);
	}
	qaws_curve_create_ellipse_2d(2, (qaws_scalar)2.6, (qaws_scalar)0.8, (qaws_scalar)0.8, 0, &circ);
	cl[0] = line; cl[1] = circ;
	p.curves = &cb; p.curve_count = 1; p.closed = 1;
	cp[0].curves = &cl[0]; cp[0].curve_count = 1; cp[0].closed = 0;
	cp[1].curves = &cl[1]; cp[1].curve_count = 1; cp[1].closed = 1;
	for (pass = 0; pass < 3; pass++)
	{
		qaws_mesh_2d* m = NULL;
		trt_stats s;
		double len = 0, off = 0;
		unsigned int cross;
		memset(&d, 0, sizeof(d));
		d.paths = &p; d.path_count = 1; d.fill_rule = QAWS_FILL_NON_ZERO;
		d.constraints = cp; d.constraint_count = 2;
		d.tolerance = (qaws_scalar)2e-3;
		if (pass >= 1) d.min_angle = 28;
		if (qaws_triangulate_region_2d(&d, &m) != QAWS_STATUS_OK || !m) { TEST_ASSERT(0, "constrained square triangulated"); continue; }
		if (pass == 2 && qaws_mesh_2d_remesh(m, (qaws_scalar)0.2, NULL, NULL, 5) != QAWS_STATUS_OK) { TEST_ASSERT(0, "constrained square remeshed"); qaws_mesh_2d_destroy(m); continue; }
		trt_measure(m, &s);
		trt_line(m, &len, &off);
		cross = trt_straddle(m);
		sprintf(msg, "constrained square (%s): area %.10f of 16, %u triangles across the line y = 1, its constraint edges %.10f long (4), vertices on it within %.1e, %u constraint edges in all, smallest angle %.1f, %u clockwise",
			pass == 0 ? "plain" : (pass == 1 ? "refined to 28 degrees" : "remeshed to 0.2"), s.area, cross, len, off, qaws_mesh_2d_get_constraint_edge_count(m), s.min_angle, s.cw);
		TEST_ASSERT(fabs(s.area - 16) < 1e-9 && cross == 0 && fabs(len - 4) < 1e-9 && off < (TRT_F32 ? 1e-5 : 1e-12) && s.cw == 0
			&& qaws_mesh_2d_get_constraint_edge_count(m) > 20 && (pass == 0 || s.min_angle > 20), msg);
		if (pass == 2)
		{
			double qa = 0;
			unsigned int bad;
			if (qaws_mesh_2d_quadrangulate(m, 1, 4) != QAWS_STATUS_OK) TEST_ASSERT(0, "constrained square quadrangulated");
			bad = trt_quads(m, &qa);
			cross = trt_straddle(m);
			trt_line(m, &len, &off);
			sprintf(msg, "constrained square, all quads: %u quads, %u triangles, %u not convex, area %.10f of 16, %u across the line, its edges %.10f long",
				qaws_mesh_2d_get_quad_count(m), qaws_mesh_2d_get_triangle_count(m), bad, qa, cross, len);
			TEST_ASSERT(qaws_mesh_2d_get_quad_count(m) > 100 && qaws_mesh_2d_get_triangle_count(m) == 0 && bad == 0 && fabs(qa - 16) < 1e-9
				&& cross == 0 && fabs(len - 4) < 1e-9, msg);
		}
		qaws_mesh_2d_destroy(m);
	}
	/* crossing constraints: both diagonals of the square, refined */
	{
		qaws_scalar d1[4] = { 0, 0, 4, 4 }, d2[4] = { 0, 4, 4, 0 };
		qaws_curve *a = NULL, *b = NULL;
		qaws_curve const* cx[2];
		qaws_path_2d xp[2];
		qaws_mesh_2d* m = NULL;
		qaws_curve_create_polyline_2d(d1, 2, 0, &a);
		qaws_curve_create_polyline_2d(d2, 2, 0, &b);
		cx[0] = a; cx[1] = b;
		xp[0].curves = &cx[0]; xp[0].curve_count = 1; xp[0].closed = 0;
		xp[1].curves = &cx[1]; xp[1].curve_count = 1; xp[1].closed = 0;
		memset(&d, 0, sizeof(d));
		d.paths = &p; d.path_count = 1; d.fill_rule = QAWS_FILL_NON_ZERO;
		d.constraints = xp; d.constraint_count = 2; d.min_angle = 25;
		if (qaws_triangulate_region_2d(&d, &m) == QAWS_STATUS_OK && m)
		{
			unsigned int i, k, across = 0;
			double len = 0;
			trt_stats s;
			trt_measure(m, &s);
			for (i = 0; i < qaws_mesh_2d_get_triangle_count(m); i++)
			{
				unsigned int v[3], l1 = 0, r1 = 0, l2 = 0, r2 = 0;
				qaws_mesh_2d_get_triangle(m, i, v);
				for (k = 0; k < 3; k++)
				{
					qaws_vec2 q = qaws_mesh_2d_get_vertex(m, v[k]);
					double s1 = (double)q.y - q.x, s2 = (double)q.y + q.x - 4, e = TRT_F32 ? 1e-5 : 1e-12;
					l1 += s1 > e; r1 += s1 < -e; l2 += s2 > e; r2 += s2 < -e;
				}
				across += (l1 && r1) || (l2 && r2);
			}
			for (i = 0; i < qaws_mesh_2d_get_constraint_edge_count(m); i++)
			{
				unsigned int e[2];
				qaws_vec2 u, w;
				qaws_mesh_2d_get_constraint_edge(m, i, e);
				u = qaws_mesh_2d_get_vertex(m, e[0]); w = qaws_mesh_2d_get_vertex(m, e[1]);
				len += hypot((double)w.x - u.x, (double)w.y - u.y);
			}
			sprintf(msg, "crossing constraints (both diagonals), 25 degrees: %u triangles across a diagonal, constraint edges %.10f long (8 sqrt 2 = %.10f), area %.10f, smallest angle %.1f",
				across, len, 8 * sqrt(2.0), s.area, s.min_angle);
			TEST_ASSERT(across == 0 && fabs(len - 8 * sqrt(2.0)) < 1e-9 && fabs(s.area - 16) < 1e-9 && s.cw == 0, msg);
		}
		else
			TEST_ASSERT(0, "crossing constraints triangulated");
		qaws_mesh_2d_destroy(m);
		qaws_curve_destroy(a); qaws_curve_destroy(b);
	}
	/* a disc: quad-dominant, then all quads */
	{
		qaws_curve* c = NULL;
		qaws_curve const* cc;
		unsigned int all;
		qaws_curve_create_ellipse_2d(0, 0, 1, 1, 0, &c);
		cc = c;
		for (all = 0; all < 2; all++)
		{
			qaws_mesh_2d* m = trt_run(&cc, 1, QAWS_FILL_NON_ZERO, 1e-3, 0, 0, NULL);
			double qa = 0, ta = 0;
			unsigned int bad, nq, nt;
			trt_stats s;
			if (!m || qaws_mesh_2d_remesh(m, (qaws_scalar)0.12, NULL, NULL, 5) != QAWS_STATUS_OK || qaws_mesh_2d_quadrangulate(m, (int)all, 4) != QAWS_STATUS_OK)
			{
				TEST_ASSERT(0, "disc quadrangulated");
				qaws_mesh_2d_destroy(m);
				continue;
			}
			bad = trt_quads(m, &qa);
			trt_measure(m, &s);
			ta = s.area;
			nq = qaws_mesh_2d_get_quad_count(m);
			nt = qaws_mesh_2d_get_triangle_count(m);
			sprintf(msg, "disc, %s: %u quads, %u triangles left, %u quads not convex, area pi - %.2e, boundary on the circle (%.2e)",
				all ? "all quads" : "quad-dominant", nq, nt, bad, TRT_PI - (qa + ta), s.on_curve);
			TEST_ASSERT(nq > 50 && bad == 0 && s.cw == 0 && fabs(TRT_PI - (qa + ta)) < 0.01 && s.on_curve < (TRT_F32 ? 1e-6 : 1e-12)
				&& (all ? nt == 0 : nt * 5 < nq), msg);
			qaws_mesh_2d_destroy(m);
		}
		qaws_curve_destroy(c);
	}
	qaws_curve_destroy(box); qaws_curve_destroy(line); qaws_curve_destroy(circ);
}

int test_91_triangulate_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 91: triangulating regions of curves\n");
	test_square();
	test_disc();
	test_annulus();
	test_degenerate();
	test_refine();
	test_remesh();
	test_constraints_and_quads();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
