/*
 * Boolean operations on regions of curves: a planar arrangement of the input
 * curves, faces with subject / clip winding numbers, Clipper2's truth tables.
 *
 *   1. sources: every curve of every path (an unclosed closed-operand path
 *      gets a closing line)
 *   2. hits: one batched intersection of all sources, self-intersections
 *      included; crossings, touches and shared stretches
 *   3. vertices: curve ends and hits, clustered within the tolerance
 *   4. edges: each source cut at its vertices; edges of closed operands that
 *      run along each other merged, their winding steps summed; edges whose
 *      steps cancel dropped
 *   5. graph: half-edges sorted around each vertex by tangent (ties by the
 *      curves' bending); faces are the cycles of "next" (turn clockwise)
 *   6. faces: counter-clockwise cycles bound faces; each clockwise cycle is
 *      a component's outside and belongs to the face around it
 *   7. windings from the unbounded face across edges; inside by the clip
 *      type on the fill rule
 *   8. output: half-edges with the inside on the left walked into loops,
 *      pieces of one source fused, line runs joined into polylines; nesting
 *      from the face components
 *   9. open subjects: cut at the closed sources, each piece kept by where
 *      its middle lies
 */

#include "qaws_clip.h"
#include "qaws_curve.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "qaws_operations.h"
#include "qaws_curve_batch.h"
#include "qaws_bezier.h"
#include "qaws_rational_bezier.h"
#include "qaws_bspline.h"
#include "qaws_nurbs.h"
#include "qaws_composite.h"
#include "qaws_platform.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_flatten.h"
#include "internal/qaws_internal_path.h"
#include "internal/qaws_internal_broadphase.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if QAWS_SCALAR_IS_FLOAT
#define CL_TOL_REL 1e-4
#define CL_AREA_REL 1e-6
#else
#define CL_TOL_REL 6.4e-9
#define CL_AREA_REL 1e-12
#endif

#define CL_PI 3.14159265358979323846
#define CL_NONE QAWS_CLIP_NONE_INDEX

/* ======================================================================== */
/*  Data                                                                    */
/* ======================================================================== */

typedef struct cl_src
{
	qaws_curve const* c;
	unsigned int operand, path, curve;   /* operand 0 subject, 1 clip, 2 open subject */
	double base, scale;                  /* input curve parameter = base + scale t */
} cl_src;

typedef struct cl_point
{
	double p[2];
	unsigned int sa, sb;                 /* sources (sb CL_NONE: one) */
	double ta, tb;
	unsigned int parent;                 /* union-find */
} cl_point;

typedef struct cl_cut
{
	double t;
	unsigned int v;
} cl_cut;

typedef struct cl_vertex
{
	double p[2];
	unsigned int sa, sb;
	double ta, tb;
	unsigned int first, count;           /* outgoing half-edges, CCW (into vout) */
} cl_vertex;

typedef struct cl_edge
{
	unsigned int src;
	double t0, t1;                       /* along src, t0 < t1 */
	unsigned int v0, v1;
	int ds, dc;                          /* winding steps across it, right to left of t0 -> t1 */
	int alive;
} cl_edge;

typedef struct cl_half
{
	double ang;                          /* tangent angle leaving the origin */
	double dir[2];
	double len;                          /* chord length of the edge */
	unsigned int next, face, pos;
	int out, seen;
} cl_half;

typedef struct cl_loop_curve
{
	qaws_curve* c;
} cl_loop_curve;

struct qaws_clip_result
{
	qaws_curve** curves;                 /* owned */
	unsigned int curve_count, curve_cap;
	qaws_curve const** views;            /* curves as const for path views */
	struct cl_out_path
	{
		unsigned int first, count;       /* into curves */
		unsigned int vfirst, vcount;     /* into vertices */
		unsigned int parent, depth;
		int hole;
	} *paths, *open_paths;
	unsigned int path_count, path_cap, open_count, open_cap;
	qaws_clip_vertex* vertices;
	unsigned int vertex_count, vertex_cap;
};

typedef struct cl_ctx
{
	qaws_clip_desc const* d;
	cl_src* src;
	unsigned int nsrc, capsrc;
	qaws_curve** temp;                   /* lines we made (polyline segments, closing lines) */
	unsigned int ntemp, captemp;
	cl_point* pts;
	unsigned int npts, cappts;
	cl_vertex* v;
	unsigned int nv;
	cl_edge* e;
	unsigned int ne, cape;
	cl_edge* oe;                         /* open-subject edges */
	unsigned int noe, capoe;
	cl_half* h;
	unsigned int* vout;
	unsigned int* face_of_cycle;
	double tol, ext;
	qaws_clip_result* r;
} cl_ctx;

#define CL_GROW(ptr, n, cap, type) \
	do { if ((n) == (cap)) { unsigned int nc_ = (cap) ? 2 * (cap) : 64; type* np_ = (type*)realloc((ptr), sizeof(type) * nc_); \
		if (!np_) return QAWS_STATUS_ALLOCATION_FAILURE; (ptr) = np_; (cap) = nc_; } } while (0)

static qaws_status cl_eval(qaws_curve const* c, double t, double* p, double* d)
{
	qaws_eval_result_2d r;
	qaws_status s;
	if (c->kind == QAWS_CURVE_KIND_BEZIER && c->degree == 1)
	{
		/* lines in double whatever the scalar type */
		qaws_scalar const* cp = ((qaws_bezier_impl const*)c->impl)->control_points;
		double ax = cp[0], ay = cp[1], bx = cp[2], by = cp[3];
		p[0] = ax + t * (bx - ax);
		p[1] = ay + t * (by - ay);
		if (t <= 0) { p[0] = ax; p[1] = ay; }
		if (t >= 1) { p[0] = bx; p[1] = by; }
		if (d) { d[0] = bx - ax; d[1] = by - ay; }
		return QAWS_STATUS_OK;
	}
	s = qaws_curve_evaluate_2d(c, (qaws_scalar)t, d ? QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 : QAWS_EVAL_FLAG_POSITION, &r);
	p[0] = r.position.x; p[1] = r.position.y;
	if (d) { d[0] = r.d1.x; d[1] = r.d1.y; }
	return s;
}

static double cl_tmin(qaws_curve const* c) { return c->parameter_range.min_value; }
static double cl_tmax(qaws_curve const* c) { return c->parameter_range.max_value; }

/* ======================================================================== */
/*  1. Sources                                                              */
/* ======================================================================== */

static qaws_status cl_temp(cl_ctx* x, qaws_curve* c)
{
	CL_GROW(x->temp, x->ntemp, x->captemp, qaws_curve*);
	x->temp[x->ntemp++] = c;
	return QAWS_STATUS_OK;
}

static qaws_status cl_line(cl_ctx* x, double const* a, double const* b, qaws_curve** out)
{
	qaws_scalar cp[4];
	qaws_bezier_desc bd;
	qaws_status s;
	cp[0] = (qaws_scalar)a[0]; cp[1] = (qaws_scalar)a[1];
	cp[2] = (qaws_scalar)b[0]; cp[3] = (qaws_scalar)b[1];
	memset(&bd, 0, sizeof(bd));
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	s = qaws_curve_create_bezier(&bd, out);
	if (s == QAWS_STATUS_OK)
		s = cl_temp(x, *out);
	return s;
}

/* A source curve. Polylines (degree-1 B-splines) are split into one line
   per segment and composites into their segments, so every corner is a
   vertex and a path doubling back on itself meets itself as two sources;
   (base, scale) maps a source's parameter back to the input curve's. */
static qaws_status cl_add_curve(cl_ctx* x, qaws_curve const* c, unsigned int operand, unsigned int path, unsigned int curve,
	double base, double scale)
{
	if (c->kind == QAWS_CURVE_KIND_BSPLINE && c->degree == 1)
	{
		qaws_bspline_impl const* impl = (qaws_bspline_impl const*)c->impl;
		unsigned int k;
		for (k = 0; k + 1 < impl->control_point_count; k++)
		{
			double a[2], b[2], t0 = impl->knots[k + 1], t1 = impl->knots[k + 2];
			qaws_curve* line = NULL;
			qaws_status s;
			a[0] = impl->control_points[2 * k]; a[1] = impl->control_points[2 * k + 1];
			b[0] = impl->control_points[2 * k + 2]; b[1] = impl->control_points[2 * k + 3];
			if ((a[0] == b[0] && a[1] == b[1]) || !(t1 > t0))
				continue;
			s = cl_line(x, a, b, &line);
			if (s == QAWS_STATUS_OK)
				s = cl_add_curve(x, line, operand, path, curve, base + scale * t0, scale * (t1 - t0));
			if (s != QAWS_STATUS_OK)
				return s;
		}
		return QAWS_STATUS_OK;
	}
	if (c->kind == QAWS_CURVE_KIND_COMPOSITE)
	{
		qaws_composite_impl const* impl = (qaws_composite_impl const*)c->impl;
		unsigned int k;
		for (k = 0; k < impl->segment_count; k++)
		{
			qaws_curve const* g = impl->segments[k];
			double a = cl_tmin(g), b = cl_tmax(g);
			/* composite parameter k + (t - a) / (b - a) */
			qaws_status s = cl_add_curve(x, g, operand, path, curve, base + scale * (k - a / (b - a)), scale / (b - a));
			if (s != QAWS_STATUS_OK)
				return s;
		}
		return QAWS_STATUS_OK;
	}
	CL_GROW(x->src, x->nsrc, x->capsrc, cl_src);
	x->src[x->nsrc].c = c;
	x->src[x->nsrc].operand = operand;
	x->src[x->nsrc].path = path;
	x->src[x->nsrc].curve = curve;
	x->src[x->nsrc].base = base;
	x->src[x->nsrc].scale = scale;
	x->nsrc++;
	return QAWS_STATUS_OK;
}

static qaws_status cl_add_paths(cl_ctx* x, qaws_path_2d const* paths, unsigned int count, unsigned int operand)
{
	unsigned int i, j;
	for (i = 0; i < count; i++)
	{
		qaws_path_2d const* p = &paths[i];
		qaws_status s;
		if (p->curve_count && !p->curves)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (j = 0; j < p->curve_count; j++)
		{
			if (!p->curves[j])
				return QAWS_STATUS_INVALID_ARGUMENT;
			if (p->curves[j]->dimension != QAWS_DIMENSION_2D)
				return QAWS_STATUS_INVALID_DIMENSION;
			s = cl_add_curve(x, p->curves[j], operand, i, j, 0.0, 1.0);
			if (s != QAWS_STATUS_OK)
				return s;
		}
		/* a region's path that does not close gets the line back to its start */
		if (operand < 2 && p->curve_count)
		{
			qaws_curve const* f = p->curves[0];
			qaws_curve const* l = p->curves[p->curve_count - 1];
			double a[2], b[2];
			cl_eval(f, cl_tmin(f), a, NULL);
			cl_eval(l, cl_tmax(l), b, NULL);
			if (hypot(a[0] - b[0], a[1] - b[1]) > x->tol)
			{
				qaws_curve* line = NULL;
				s = cl_line(x, b, a, &line);
				if (s == QAWS_STATUS_OK)
					s = cl_add_curve(x, line, operand, i, p->curve_count, 0.0, 1.0);
				if (s != QAWS_STATUS_OK)
					return s;
			}
		}
	}
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  2-3. Hits and vertices                                                  */
/* ======================================================================== */

static qaws_status cl_point_add(cl_ctx* x, double const* p, unsigned int sa, double ta, unsigned int sb, double tb)
{
	cl_point* q;
	CL_GROW(x->pts, x->npts, x->cappts, cl_point);
	q = &x->pts[x->npts];
	q->p[0] = p[0]; q->p[1] = p[1];
	q->sa = sa; q->ta = ta; q->sb = sb; q->tb = tb;
	q->parent = x->npts;
	x->npts++;
	return QAWS_STATUS_OK;
}

static unsigned int cl_find(cl_point* pts, unsigned int i)
{
	while (pts[i].parent != i)
	{
		pts[i].parent = pts[pts[i].parent].parent;
		i = pts[i].parent;
	}
	return i;
}

static unsigned int cl_cell_hash(long long cx, long long cy, unsigned int mask)
{
	unsigned long long h = (unsigned long long)cx * 0x9E3779B97F4A7C15ull ^ (unsigned long long)cy * 0xC2B2AE3D27D4EB4Full;
	return (unsigned int)(h >> 32) & mask;
}

/* points within the tolerance are one vertex: a hash grid of tolerance-sized
   cells, each point joined with the points of the 9 cells around it */
static qaws_status cl_cluster(cl_ctx* x)
{
	unsigned int size = 1, mask, i, *head, *link;
	while (size < 2 * x->npts + 16) size <<= 1;
	mask = size - 1;
	head = (unsigned int*)malloc(sizeof(unsigned int) * (size + x->npts + 1));
	if (!head)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	link = head + size;
	for (i = 0; i < size; i++) head[i] = CL_NONE;
	for (i = 0; i < x->npts; i++)
	{
		long long cx = (long long)floor(x->pts[i].p[0] / x->tol), cy = (long long)floor(x->pts[i].p[1] / x->tol);
		int dx, dy;
		for (dx = -1; dx <= 1; dx++)
			for (dy = -1; dy <= 1; dy++)
			{
				unsigned int k = head[cl_cell_hash(cx + dx, cy + dy, mask)];
				for (; k != CL_NONE; k = link[k])
					if (hypot(x->pts[k].p[0] - x->pts[i].p[0], x->pts[k].p[1] - x->pts[i].p[1]) <= x->tol)
					{
						unsigned int a = cl_find(x->pts, k), b = cl_find(x->pts, i);
						if (a != b) x->pts[a > b ? a : b].parent = a < b ? a : b;
					}
			}
		{
			unsigned int hsh = cl_cell_hash(cx, cy, mask);
			link[i] = head[hsh];
			head[hsh] = i;
		}
	}
	free(head);
	/* compact roots into vertices: averaged position, a crossing's sources
	   preferred for the record */
	{
		unsigned int* id = (unsigned int*)malloc(sizeof(unsigned int) * (x->npts + 1));
		unsigned int* cnt;
		if (!id)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		x->nv = 0;
		for (i = 0; i < x->npts; i++)
			id[i] = cl_find(x->pts, i) == i ? x->nv++ : CL_NONE;
		/* every point straight to its root before parents are reused below */
		for (i = 0; i < x->npts; i++)
			x->pts[i].parent = cl_find(x->pts, i);
		x->v = (cl_vertex*)calloc(x->nv ? x->nv : 1, sizeof(cl_vertex));
		cnt = (unsigned int*)calloc(x->nv ? x->nv : 1, sizeof(unsigned int));
		if (!x->v || !cnt)
		{
			free(id); free(cnt);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		for (i = 0; i < x->nv; i++)
			x->v[i].sa = x->v[i].sb = CL_NONE;
		for (i = 0; i < x->npts; i++)
		{
			unsigned int r = id[x->pts[i].parent];
			cl_vertex* v = &x->v[r];
			v->p[0] += x->pts[i].p[0];
			v->p[1] += x->pts[i].p[1];
			cnt[r]++;
			if (v->sa == CL_NONE || (v->sb == CL_NONE && x->pts[i].sb != CL_NONE))
			{
				v->sa = x->pts[i].sa; v->ta = x->pts[i].ta;
				v->sb = x->pts[i].sb; v->tb = x->pts[i].tb;
			}
			x->pts[i].parent = r;   /* from here on: the vertex id */
		}
		for (i = 0; i < x->nv; i++)
		{
			x->v[i].p[0] /= cnt[i];
			x->v[i].p[1] /= cnt[i];
		}
		free(id);
		free(cnt);
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------------ */
/*  Line x line hits, in double: from float or integer inputs the cross     */
/*  products are exact, so each hit is one rounding away; collinear lines   */
/*  overlap exactly                                                         */
/* ------------------------------------------------------------------------ */

typedef struct cl_lseg
{
	double a[2], b[2], lo[2], hi[2];
	unsigned int src;
} cl_lseg;

static int cl_is_line_src(cl_ctx const* x, unsigned int i)
{
	qaws_curve const* c = x->src[i].c;
	return c->kind == QAWS_CURVE_KIND_BEZIER && c->degree == 1;
}

static double cl_cross(double ax, double ay, double bx, double by)
{
	return ax * by - ay * bx;
}

/* the line's parameter at its end points t = 0 and 1 (a Bezier: [0, 1]) */
static qaws_status cl_line_pair(cl_ctx* x, cl_lseg const* A, cl_lseg const* B)
{
	double d1[2] = { A->b[0] - A->a[0], A->b[1] - A->a[1] }, d2[2] = { B->b[0] - B->a[0], B->b[1] - B->a[1] };
	double w[2] = { B->a[0] - A->a[0], B->a[1] - A->a[1] };
	double den = cl_cross(d1[0], d1[1], d2[0], d2[1]);
	double l1 = d1[0] * d1[0] + d1[1] * d1[1];
	int a_open = x->src[A->src].operand == 2, b_open = x->src[B->src].operand == 2;
	if (a_open && b_open)
		return QAWS_STATUS_OK;
	if (den != 0)
	{
		double t = cl_cross(w[0], w[1], d2[0], d2[1]) / den, u = cl_cross(w[0], w[1], d1[0], d1[1]) / den, p[2];
		if (t < 0 || t > 1 || u < 0 || u > 1)
			return QAWS_STATUS_OK;
		/* both at their ends: a corner or a junction, already a point */
		if ((t == 0 || t == 1) && (u == 0 || u == 1))
			return QAWS_STATUS_OK;
		p[0] = A->a[0] + t * d1[0];
		p[1] = A->a[1] + t * d1[1];
		/* an end exactly on the other line stays exact */
		if (t == 0) { p[0] = A->a[0]; p[1] = A->a[1]; }
		else if (t == 1) { p[0] = A->b[0]; p[1] = A->b[1]; }
		else if (u == 0) { p[0] = B->a[0]; p[1] = B->a[1]; }
		else if (u == 1) { p[0] = B->b[0]; p[1] = B->b[1]; }
		return cl_point_add(x, p, A->src, t, B->src, u);
	}
	/* parallel: collinear lines overlap on the projections of their ends */
	if (cl_cross(w[0], w[1], d1[0], d1[1]) != 0 || !(l1 > 0))
		return QAWS_STATUS_OK;
	{
		double sc = (B->a[0] - A->a[0]) * d1[0] + (B->a[1] - A->a[1]) * d1[1];
		double sd = (B->b[0] - A->a[0]) * d1[0] + (B->b[1] - A->a[1]) * d1[1];
		double s0 = sc / l1, s1 = sd / l1, lo = s0 < s1 ? s0 : s1, hi = s0 < s1 ? s1 : s0;
		double ta, tb, ua, ub, pa[2], pb[2];
		qaws_status s;
		if (lo < 0) lo = 0;
		if (hi > 1) hi = 1;
		if (!(hi >= lo))
			return QAWS_STATUS_OK;
		/* B's parameters at the overlap's ends: u = (s - s0) / (s1 - s0) */
		ta = lo; tb = hi;
		ua = (lo - s0) / (s1 - s0); ub = (hi - s0) / (s1 - s0);
		pa[0] = A->a[0] + ta * d1[0]; pa[1] = A->a[1] + ta * d1[1];
		pb[0] = A->a[0] + tb * d1[0]; pb[1] = A->a[1] + tb * d1[1];
		if (ta == 0) { pa[0] = A->a[0]; pa[1] = A->a[1]; }
		if (tb == 1) { pb[0] = A->b[0]; pb[1] = A->b[1]; }
		if (ua == 0 || ua == 1) { pa[0] = ua == 0 ? B->a[0] : B->b[0]; pa[1] = ua == 0 ? B->a[1] : B->b[1]; }
		if (ub == 0 || ub == 1) { pb[0] = ub == 0 ? B->a[0] : B->b[0]; pb[1] = ub == 0 ? B->a[1] : B->b[1]; }
		s = cl_point_add(x, pa, A->src, ta, B->src, ua);
		if (s == QAWS_STATUS_OK && hi > lo)
			s = cl_point_add(x, pb, A->src, tb, B->src, ub);
		return s;
	}
}

typedef struct cl_line_ctx
{
	cl_ctx* x;
	cl_lseg const* L;
} cl_line_ctx;

static qaws_status cl_line_visit(void* user, unsigned int i, unsigned int j)
{
	cl_line_ctx* c = (cl_line_ctx*)user;
	cl_lseg const* A = &c->L[i];
	cl_lseg const* B = &c->L[j];
	return A->src < B->src ? cl_line_pair(c->x, A, B) : cl_line_pair(c->x, B, A);
}

static qaws_status cl_line_hits(cl_ctx* x)
{
	unsigned int n = 0, i, j;
	cl_lseg* L = (cl_lseg*)malloc(sizeof(cl_lseg) * (x->nsrc ? x->nsrc : 1));
	qaws_status s = QAWS_STATUS_OK;
	if (!L)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->nsrc; i++)
	{
		qaws_bezier_impl const* impl;
		cl_lseg* g;
		if (!cl_is_line_src(x, i))
			continue;
		impl = (qaws_bezier_impl const*)x->src[i].c->impl;
		g = &L[n++];
		g->a[0] = impl->control_points[0]; g->a[1] = impl->control_points[1];
		g->b[0] = impl->control_points[2]; g->b[1] = impl->control_points[3];
		g->lo[0] = g->a[0] < g->b[0] ? g->a[0] : g->b[0];
		g->lo[1] = g->a[1] < g->b[1] ? g->a[1] : g->b[1];
		g->hi[0] = g->a[0] < g->b[0] ? g->b[0] : g->a[0];
		g->hi[1] = g->a[1] < g->b[1] ? g->b[1] : g->a[1];
		g->src = i;
	}
	/* the pairs of overlapping boxes, from the shared uniform grid */
	if (n > 1)
	{
		qaws_bp_box* boxes = (qaws_bp_box*)malloc(sizeof(qaws_bp_box) * n);
		cl_line_ctx lc;
		if (!boxes) { free(L); return QAWS_STATUS_ALLOCATION_FAILURE; }
		for (i = 0; i < n; i++)
		{
			boxes[i].lo[0] = L[i].lo[0]; boxes[i].hi[0] = L[i].hi[0];
			boxes[i].lo[1] = L[i].lo[1]; boxes[i].hi[1] = L[i].hi[1];
			boxes[i].lo[2] = boxes[i].hi[2] = 0;
		}
		lc.x = x;
		lc.L = L;
		s = qaws_internal_broadphase(boxes, n, 2, NULL, cl_line_visit, &lc, NULL);
		free(boxes);
	}
	(void)j;
	free(L);
	return s;
}

static qaws_status cl_hits(cl_ctx* x)
{
	qaws_curve const** curves = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (x->nsrc ? x->nsrc : 1));
	unsigned int* fam = (unsigned int*)malloc(sizeof(unsigned int) * (x->nsrc ? x->nsrc : 1));
	qaws_curve_batch_desc bd;
	qaws_curve_batch_hit_2d* hits = NULL;
	unsigned int n = 0, cap = 0, i, ncurved = 0;
	qaws_status s = QAWS_STATUS_OK;
	if (!curves || !fam)
	{
		free((void*)curves); free(fam);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* lines meet lines below (exactly); the batch takes every pair with a
	   curve in it: all lines in one family, each curve in its own */
	for (i = 0; i < x->nsrc; i++)
	{
		curves[i] = x->src[i].c;
		fam[i] = cl_is_line_src(x, i) ? 0u : i + 1;
		ncurved += fam[i] != 0;
	}
	memset(&bd, 0, sizeof(bd));
	bd.curves = curves;
	bd.curve_count = x->nsrc;
	bd.families = fam;
	bd.flags = QAWS_CURVE_BATCH_SELF;
	bd.executor = x->d->executor;
	if (ncurved)
	{
		cap = 4 * x->nsrc + 64;
		hits = (qaws_curve_batch_hit_2d*)malloc(sizeof(qaws_curve_batch_hit_2d) * cap);
		if (!hits)
			s = QAWS_STATUS_ALLOCATION_FAILURE;
		else
			s = qaws_curve_batch_find_intersections_2d(&bd, hits, cap, &n, NULL);
		if (s == QAWS_STATUS_OK && n > cap)
		{
			/* more than guessed: once more with room for all */
			free(hits);
			cap = n;
			hits = (qaws_curve_batch_hit_2d*)malloc(sizeof(qaws_curve_batch_hit_2d) * cap);
			if (!hits)
				s = QAWS_STATUS_ALLOCATION_FAILURE;
			else
				s = qaws_curve_batch_find_intersections_2d(&bd, hits, cap, &n, NULL);
		}
		if (n > cap) n = cap;
	}
	if (s == QAWS_STATUS_OK)
		s = cl_line_hits(x);
	/* curve ends */
	for (i = 0; s == QAWS_STATUS_OK && i < x->nsrc; i++)
	{
		double p[2];
		cl_eval(x->src[i].c, cl_tmin(x->src[i].c), p, NULL);
		s = cl_point_add(x, p, i, cl_tmin(x->src[i].c), CL_NONE, 0);
		if (s == QAWS_STATUS_OK)
		{
			cl_eval(x->src[i].c, cl_tmax(x->src[i].c), p, NULL);
			s = cl_point_add(x, p, i, cl_tmax(x->src[i].c), CL_NONE, 0);
		}
	}
	for (i = 0; s == QAWS_STATUS_OK && i < n; i++)
	{
		qaws_curve_batch_hit_2d const* h = &hits[i];
		double p[2];
		/* open subjects are not cut by each other */
		if (x->src[h->curve_a].operand == 2 && x->src[h->curve_b].operand == 2)
			continue;
		/* two sources meeting at their ends (a corner, a junction): the ends
		   are points already */
		if (h->kind != QAWS_CURVE_HIT_OVERLAP &&
			(h->parameter_a <= cl_tmin(x->src[h->curve_a].c) || h->parameter_a >= cl_tmax(x->src[h->curve_a].c)) &&
			(h->parameter_b <= cl_tmin(x->src[h->curve_b].c) || h->parameter_b >= cl_tmax(x->src[h->curve_b].c)))
			continue;
		p[0] = h->position.x; p[1] = h->position.y;
		s = cl_point_add(x, p, h->curve_a, h->parameter_a, h->curve_b, h->parameter_b);
		if (s == QAWS_STATUS_OK && h->kind == QAWS_CURVE_HIT_OVERLAP)
		{
			cl_eval(x->src[h->curve_a].c, h->parameter_a_end, p, NULL);
			s = cl_point_add(x, p, h->curve_a, h->parameter_a_end, h->curve_b, h->parameter_b_end);
		}
	}
	free(hits);
	free((void*)curves);
	free(fam);
	return s == QAWS_STATUS_OK ? cl_cluster(x) : s;
}

/* ======================================================================== */
/*  4. Edges                                                                */
/* ======================================================================== */

static int cl_cmp_cut(void const* a, void const* b)
{
	double ta = ((cl_cut const*)a)->t, tb = ((cl_cut const*)b)->t;
	return ta < tb ? -1 : (ta > tb ? 1 : 0);
}

/* distance from q to the curve of an edge (Gauss-Newton on its range) */
static double cl_dist_to_edge(cl_ctx const* x, cl_edge const* e, double const* q)
{
	qaws_curve const* c = x->src[e->src].c;
	double t = 0.5 * (e->t0 + e->t1), p[2], d[2];
	unsigned int it;
	/* start from the best of a few samples */
	{
		double best = 1e300;
		unsigned int k;
		for (k = 0; k <= 8; k++)
		{
			double tk = e->t0 + (e->t1 - e->t0) * k / 8, dk;
			cl_eval(c, tk, p, NULL);
			dk = hypot(p[0] - q[0], p[1] - q[1]);
			if (dk < best) { best = dk; t = tk; }
		}
	}
	for (it = 0; it < 30; it++)
	{
		double dd, step;
		cl_eval(c, t, p, d);
		dd = d[0] * d[0] + d[1] * d[1];
		if (!(dd > 0)) break;
		step = ((p[0] - q[0]) * d[0] + (p[1] - q[1]) * d[1]) / dd;
		t -= step;
		if (t < e->t0) t = e->t0;
		if (t > e->t1) t = e->t1;
		if (fabs(step) <= (e->t1 - e->t0) * 1e-15) break;
	}
	cl_eval(c, t, p, NULL);
	return hypot(p[0] - q[0], p[1] - q[1]);
}

/* each source cut at its vertices into edges */
static qaws_status cl_edges(cl_ctx* x)
{
	cl_cut* cuts = NULL;
	unsigned int* count = (unsigned int*)calloc(x->nsrc + 1, sizeof(unsigned int));
	unsigned int* start = (unsigned int*)calloc(x->nsrc + 1, sizeof(unsigned int));
	unsigned int i, s;
	if (!count || !start)
	{
		free(count); free(start);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < x->npts; i++)
	{
		count[x->pts[i].sa]++;
		if (x->pts[i].sb != CL_NONE) count[x->pts[i].sb]++;
	}
	for (s = 0; s < x->nsrc; s++)
		start[s + 1] = start[s] + count[s];
	cuts = (cl_cut*)malloc(sizeof(cl_cut) * (start[x->nsrc] + 1));
	if (!cuts)
	{
		free(count); free(start);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memset(count, 0, sizeof(unsigned int) * x->nsrc);
	for (i = 0; i < x->npts; i++)
	{
		cl_point const* q = &x->pts[i];
		cuts[start[q->sa] + count[q->sa]].t = q->ta;
		cuts[start[q->sa] + count[q->sa]++].v = q->parent;
		if (q->sb != CL_NONE)
		{
			cuts[start[q->sb] + count[q->sb]].t = q->tb;
			cuts[start[q->sb] + count[q->sb]++].v = q->parent;
		}
	}
	for (s = 0; s < x->nsrc; s++)
	{
		cl_cut* c = &cuts[start[s]];
		unsigned int n = count[s], k, m = 0;
		int open = x->src[s].operand == 2;
		qsort(c, n, sizeof(cl_cut), cl_cmp_cut);
		/* cuts of one vertex next to each other are one cut, unless the curve
		   leaves the vertex between them (a closed curve's two ends) */
		for (k = 0; k < n; k++)
		{
			if (m && c[m - 1].v == c[k].v)
			{
				double mid[2];
				cl_eval(x->src[s].c, 0.5 * (c[m - 1].t + c[k].t), mid, NULL);
				if (hypot(mid[0] - x->v[c[k].v].p[0], mid[1] - x->v[c[k].v].p[1]) <= 2 * x->tol)
					continue;
			}
			c[m++] = c[k];
		}
		for (k = 0; k + 1 < m; k++)
		{
			cl_edge e;
			memset(&e, 0, sizeof(e));
			e.src = s;
			e.t0 = c[k].t; e.t1 = c[k + 1].t;
			e.v0 = c[k].v; e.v1 = c[k + 1].v;
			e.ds = x->src[s].operand == 0 ? 1 : 0;
			e.dc = x->src[s].operand == 1 ? 1 : 0;
			e.alive = 1;
			if (!(e.t1 > e.t0))
				continue;
			if (e.v0 == e.v1)
			{
				/* a loop only when it leaves the vertex */
				double mid[2];
				cl_eval(x->src[s].c, 0.5 * (e.t0 + e.t1), mid, NULL);
				if (hypot(mid[0] - x->v[e.v0].p[0], mid[1] - x->v[e.v0].p[1]) <= 2 * x->tol)
					continue;
			}
			if (open)
			{
				if (x->noe == x->capoe)
				{
					unsigned int nc = x->capoe ? 2 * x->capoe : 64;
					cl_edge* ne = (cl_edge*)realloc(x->oe, sizeof(cl_edge) * nc);
					if (!ne) { free(cuts); free(count); free(start); return QAWS_STATUS_ALLOCATION_FAILURE; }
					x->oe = ne; x->capoe = nc;
				}
				x->oe[x->noe++] = e;
			}
			else
			{
				if (x->ne == x->cape)
				{
					unsigned int nc = x->cape ? 2 * x->cape : 64;
					cl_edge* ne = (cl_edge*)realloc(x->e, sizeof(cl_edge) * nc);
					if (!ne) { free(cuts); free(count); free(start); return QAWS_STATUS_ALLOCATION_FAILURE; }
					x->e = ne; x->cape = nc;
				}
				x->e[x->ne++] = e;
			}
		}
	}
	free(cuts); free(count); free(start);
	return QAWS_STATUS_OK;
}

/* edges joining the same two vertices that run along each other become one
   edge; then edges whose winding steps cancel go */

typedef struct cl_pair_key { unsigned int lo, hi, e; } cl_pair_key;

static int cl_cmp_key(void const* a, void const* b)
{
	cl_pair_key const* p = (cl_pair_key const*)a;
	cl_pair_key const* q = (cl_pair_key const*)b;
	if (p->lo != q->lo) return p->lo < q->lo ? -1 : 1;
	if (p->hi != q->hi) return p->hi < q->hi ? -1 : 1;
	return p->e < q->e ? -1 : (p->e > q->e ? 1 : 0);
}

static qaws_status cl_merge_edges(cl_ctx* x)
{
	cl_pair_key* k = (cl_pair_key*)malloc(sizeof(cl_pair_key) * (x->ne ? x->ne : 1));
	unsigned int i, j, n = 0;
	if (!k)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < x->ne; i++)
	{
		k[i].lo = x->e[i].v0 < x->e[i].v1 ? x->e[i].v0 : x->e[i].v1;
		k[i].hi = x->e[i].v0 < x->e[i].v1 ? x->e[i].v1 : x->e[i].v0;
		k[i].e = i;
	}
	qsort(k, x->ne, sizeof(cl_pair_key), cl_cmp_key);
	for (i = 0; i < x->ne; i = j)
	{
		for (j = i + 1; j < x->ne && k[j].lo == k[i].lo && k[j].hi == k[i].hi; j++)
			;
		if (j - i < 2)
			continue;
		{
			unsigned int a, b;
			for (a = i; a < j; a++)
			{
				cl_edge* ea = &x->e[k[a].e];
				double ma[2];
				if (!ea->alive) continue;
				cl_eval(x->src[ea->src].c, 0.5 * (ea->t0 + ea->t1), ma, NULL);
				for (b = a + 1; b < j; b++)
				{
					cl_edge* eb = &x->e[k[b].e];
					double mb[2];
					int same;
					if (!eb->alive) continue;
					cl_eval(x->src[eb->src].c, 0.5 * (eb->t0 + eb->t1), mb, NULL);
					if (cl_dist_to_edge(x, ea, mb) > 4 * x->tol || cl_dist_to_edge(x, eb, ma) > 4 * x->tol)
						continue;
					/* direction: same end vertices, or for a loop the same tangent sense */
					if (ea->v0 != ea->v1)
						same = ea->v0 == eb->v0;
					else
					{
						double pa[2], da[2], pb[2], db[2];
						cl_eval(x->src[ea->src].c, ea->t0, pa, da);
						cl_eval(x->src[eb->src].c, eb->t0, pb, db);
						same = da[0] * db[0] + da[1] * db[1] >= 0;
					}
					ea->ds += same ? eb->ds : -eb->ds;
					ea->dc += same ? eb->dc : -eb->dc;
					eb->alive = 0;
				}
			}
		}
	}
	free(k);
	for (i = 0; i < x->ne; i++)
		if (x->e[i].alive && (x->e[i].ds || x->e[i].dc))
			x->e[n++] = x->e[i];
	x->ne = n;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  5. Graph                                                                */
/* ======================================================================== */

static unsigned int cl_origin(cl_ctx const* x, unsigned int h)
{
	cl_edge const* e = &x->e[h >> 1];
	return (h & 1) ? e->v1 : e->v0;
}

static unsigned int cl_target(cl_ctx const* x, unsigned int h)
{
	cl_edge const* e = &x->e[h >> 1];
	return (h & 1) ? e->v0 : e->v1;
}

/* a point at arc distance about s from the origin along half-edge h */
static void cl_probe(cl_ctx const* x, unsigned int h, double s, double* q)
{
	cl_edge const* e = &x->e[h >> 1];
	qaws_curve const* c = x->src[e->src].c;
	double t = (h & 1) ? e->t1 : e->t0, p[2], d[2], sp, dt;
	cl_eval(c, t, p, d);
	sp = hypot(d[0], d[1]);
	dt = sp > 0 ? s / sp : 1e-3 * (e->t1 - e->t0);
	if (dt > 0.5 * (e->t1 - e->t0)) dt = 0.5 * (e->t1 - e->t0);
	cl_eval(c, (h & 1) ? t - dt : t + dt, q, NULL);
}

static void cl_half_setup(cl_ctx* x, unsigned int h)
{
	cl_edge const* e = &x->e[h >> 1];
	qaws_curve const* c = x->src[e->src].c;
	double p[2], d[2], q[2], n;
	cl_half* hh = &x->h[h];
	cl_eval(c, (h & 1) ? e->t1 : e->t0, p, d);
	if (h & 1) { d[0] = -d[0]; d[1] = -d[1]; }
	cl_eval(c, (h & 1) ? e->t0 : e->t1, q, NULL);
	hh->len = hypot(q[0] - p[0], q[1] - p[1]);
	n = hypot(d[0], d[1]);
	if (!(n > 0))
	{
		/* a cusp: the direction of a nearby point */
		double r[2];
		cl_probe(x, h, 1e-6 * (hh->len > 0 ? hh->len : x->ext), r);
		d[0] = r[0] - p[0]; d[1] = r[1] - p[1];
		n = hypot(d[0], d[1]);
	}
	hh->dir[0] = n > 0 ? d[0] / n : 1;
	hh->dir[1] = n > 0 ? d[1] / n : 0;
	hh->ang = atan2(hh->dir[1], hh->dir[0]);
	hh->next = hh->face = CL_NONE;
	hh->out = hh->seen = 0;
}

/* CCW order of half-edges leaving one vertex: by tangent angle; tangent ties
   by where the curves are a short way along (the one bending left first) */
static int cl_before(cl_ctx const* x, unsigned int a, unsigned int b)
{
	cl_half const* ha = &x->h[a];
	cl_half const* hb = &x->h[b];
	double da = ha->ang - hb->ang;
	if (fabs(da) > 1e-9)
		return da < 0;
	{
		double s = 1e-4 * (ha->len < hb->len ? ha->len : hb->len), qa[2], qb[2], o[2], ca, cb;
		unsigned int v = cl_origin(x, a);
		o[0] = x->v[v].p[0]; o[1] = x->v[v].p[1];
		if (!(s > 0)) s = 1e-6 * x->ext;
		cl_probe(x, a, s, qa);
		cl_probe(x, b, s, qb);
		/* signed side of each probe from the common tangent */
		ca = ha->dir[0] * (qa[1] - o[1]) - ha->dir[1] * (qa[0] - o[0]);
		cb = ha->dir[0] * (qb[1] - o[1]) - ha->dir[1] * (qb[0] - o[0]);
		if (ca != cb)
			return ca < cb;
		return a < b;
	}
}

static qaws_status cl_graph(cl_ctx* x)
{
	unsigned int nh = 2 * x->ne, i, *cnt;
	x->h = (cl_half*)malloc(sizeof(cl_half) * (nh ? nh : 1));
	x->vout = (unsigned int*)malloc(sizeof(unsigned int) * (nh ? nh : 1));
	cnt = (unsigned int*)calloc(x->nv + 1, sizeof(unsigned int));
	if (!x->h || !x->vout || !cnt)
	{
		free(cnt);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nh; i++)
	{
		cl_half_setup(x, i);
		cnt[cl_origin(x, i)]++;
	}
	{
		unsigned int acc = 0;
		for (i = 0; i < x->nv; i++)
		{
			x->v[i].first = acc;
			acc += cnt[i];
			x->v[i].count = 0;
		}
	}
	for (i = 0; i < nh; i++)
	{
		cl_vertex* v = &x->v[cl_origin(x, i)];
		x->vout[v->first + v->count++] = i;
	}
	/* sort each vertex's fan (small: insertion sort) */
	for (i = 0; i < x->nv; i++)
	{
		unsigned int* f = &x->vout[x->v[i].first];
		unsigned int n = x->v[i].count, a, b;
		for (a = 1; a < n; a++)
		{
			unsigned int key = f[a];
			for (b = a; b > 0 && cl_before(x, key, f[b - 1]); b--)
				f[b] = f[b - 1];
			f[b] = key;
		}
		for (a = 0; a < n; a++)
			x->h[f[a]].pos = a;
	}
	/* next(h) = the half-edge clockwise of twin(h) around h's target */
	for (i = 0; i < nh; i++)
	{
		unsigned int tw = i ^ 1;
		cl_vertex const* v = &x->v[cl_origin(x, tw)];
		unsigned int p = x->h[tw].pos;
		x->h[i].next = x->vout[v->first + (p + v->count - 1) % v->count];
	}
	free(cnt);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  6-7. Faces and windings                                                 */
/* ======================================================================== */

/* integral of x dy along half-edge h */
static double cl_half_area(cl_ctx const* x, unsigned int h)
{
	cl_edge const* e = &x->e[h >> 1];
	double a = 0.0;
	if (cl_is_line_src(x, e->src))
	{
		/* a line from p to q: (x_p + x_q) / 2 (y_q - y_p), in double */
		double p[2], q[2];
		cl_eval(x->src[e->src].c, (h & 1) ? e->t1 : e->t0, p, NULL);
		cl_eval(x->src[e->src].c, (h & 1) ? e->t0 : e->t1, q, NULL);
		return 0.5 * (p[0] + q[0]) * (q[1] - p[1]);
	}
	qaws_internal_curve_area_2d(x->src[e->src].c, (h & 1) ? e->t1 : e->t0, (h & 1) ? e->t0 : e->t1,
		x->ext * x->ext * CL_AREA_REL, &a);
	return a;
}

/* winding of the cycle starting at half-edge h0 around point q: chords of
   pieces split until none can pass on the wrong side of q */
static double cl_piece_angle(cl_ctx const* x, qaws_curve const* c, double t0, double t1, double const* p0, double const* p1,
	double const* q, unsigned int depth)
{
	double tm = 0.5 * (t0 + t1), pm[2], ab[2], dev, dist, l2, u, cx, cy;
	cl_eval(c, tm, pm, NULL);
	ab[0] = p1[0] - p0[0]; ab[1] = p1[1] - p0[1];
	l2 = ab[0] * ab[0] + ab[1] * ab[1];
	/* deviation of the middle from the chord, distance of q from the chord */
	u = l2 > 0 ? ((pm[0] - p0[0]) * ab[0] + (pm[1] - p0[1]) * ab[1]) / l2 : 0;
	dev = hypot(pm[0] - p0[0] - u * ab[0], pm[1] - p0[1] - u * ab[1]);
	u = l2 > 0 ? ((q[0] - p0[0]) * ab[0] + (q[1] - p0[1]) * ab[1]) / l2 : 0;
	if (u < 0) u = 0;
	if (u > 1) u = 1;
	dist = hypot(q[0] - p0[0] - u * ab[0], q[1] - p0[1] - u * ab[1]);
	if (depth < 40 && (dist <= 8 * dev || depth < 2))
		return cl_piece_angle(x, c, t0, tm, p0, pm, q, depth + 1) + cl_piece_angle(x, c, tm, t1, pm, p1, q, depth + 1);
	cx = (p0[0] - q[0]) * (p1[1] - q[1]) - (p0[1] - q[1]) * (p1[0] - q[0]);
	cy = (p0[0] - q[0]) * (p1[0] - q[0]) + (p0[1] - q[1]) * (p1[1] - q[1]);
	return atan2(cx, cy);
}

static int cl_cycle_winding(cl_ctx const* x, unsigned int h0, double const* q)
{
	double angle = 0.0;
	unsigned int h = h0;
	do
	{
		cl_edge const* e = &x->e[h >> 1];
		qaws_curve const* c = x->src[e->src].c;
		double t0 = (h & 1) ? e->t1 : e->t0, t1 = (h & 1) ? e->t0 : e->t1;
		double const* p0 = x->v[cl_origin(x, h)].p;
		double const* p1 = x->v[cl_target(x, h)].p;
		angle += cl_piece_angle(x, c, t0, t1, p0, p1, q, 0);
		h = x->h[h].next;
	} while (h != h0);
	return (int)floor(angle / (2 * CL_PI) + 0.5);
}

typedef struct cl_cycle
{
	unsigned int h0;
	double area;
	double lo[2], hi[2];
	unsigned int comp;
} cl_cycle;

static unsigned int cl_vfind(unsigned int* p, unsigned int i)
{
	while (p[i] != i) { p[i] = p[p[i]]; i = p[i]; }
	return i;
}

/* faces: CCW cycles; each CW cycle belongs to the face around it (the
   smallest CCW cycle of another component that winds around it), or to the
   unbounded face 0. Writes x->h[].face and returns the face count. */
static qaws_status cl_faces(cl_ctx* x, unsigned int* out_faces, unsigned int** out_cycle_face, cl_cycle** out_cycles,
	unsigned int* out_ncycles)
{
	unsigned int nh = 2 * x->ne, i, ncyc = 0, capc = 0, nfaces = 1;
	cl_cycle* cyc = NULL;
	unsigned int* cface = NULL;
	unsigned int* vcomp = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv ? x->nv : 1));
	unsigned int* hcyc = (unsigned int*)malloc(sizeof(unsigned int) * (nh ? nh : 1));
	if (!vcomp || !hcyc)
	{
		free(vcomp); free(hcyc);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* connected components of the graph */
	for (i = 0; i < x->nv; i++) vcomp[i] = i;
	for (i = 0; i < x->ne; i++)
	{
		unsigned int a = cl_vfind(vcomp, x->e[i].v0), b = cl_vfind(vcomp, x->e[i].v1);
		if (a != b) vcomp[a] = b;
	}
	/* cycles */
	for (i = 0; i < nh; i++) hcyc[i] = CL_NONE;
	for (i = 0; i < nh; i++)
	{
		unsigned int h = i;
		cl_cycle c;
		if (hcyc[i] != CL_NONE) continue;
		c.h0 = i; c.area = 0.0;
		c.lo[0] = c.lo[1] = 1e300; c.hi[0] = c.hi[1] = -1e300;
		c.comp = cl_vfind(vcomp, cl_origin(x, i));
		do
		{
			double const* p = x->v[cl_origin(x, h)].p;
			hcyc[h] = ncyc;
			c.area += cl_half_area(x, h);
			if (p[0] < c.lo[0]) c.lo[0] = p[0];
			if (p[1] < c.lo[1]) c.lo[1] = p[1];
			if (p[0] > c.hi[0]) c.hi[0] = p[0];
			if (p[1] > c.hi[1]) c.hi[1] = p[1];
			/* the box of the curve between, roughly: its middle */
			{
				cl_edge const* e = &x->e[h >> 1];
				double m[2];
				unsigned int k;
				for (k = 1; k < 4; k++)
				{
					cl_eval(x->src[e->src].c, e->t0 + (e->t1 - e->t0) * k / 4, m, NULL);
					if (m[0] < c.lo[0]) c.lo[0] = m[0];
					if (m[1] < c.lo[1]) c.lo[1] = m[1];
					if (m[0] > c.hi[0]) c.hi[0] = m[0];
					if (m[1] > c.hi[1]) c.hi[1] = m[1];
				}
			}
			h = x->h[h].next;
		} while (h != i && hcyc[h] == CL_NONE);
		if (ncyc == capc)
		{
			unsigned int nc = capc ? 2 * capc : 64;
			cl_cycle* g = (cl_cycle*)realloc(cyc, sizeof(cl_cycle) * nc);
			if (!g) { free(cyc); free(vcomp); free(hcyc); return QAWS_STATUS_ALLOCATION_FAILURE; }
			cyc = g; capc = nc;
		}
		cyc[ncyc++] = c;
	}
	cface = (unsigned int*)malloc(sizeof(unsigned int) * (ncyc ? ncyc : 1));
	if (!cface)
	{
		free(cyc); free(vcomp); free(hcyc);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < ncyc; i++)
		cface[i] = cyc[i].area > 0 ? nfaces++ : CL_NONE;
	for (i = 0; i < ncyc; i++)
	{
		unsigned int j, best = CL_NONE;
		double q[2], best_area = 1e300;
		if (cface[i] != CL_NONE) continue;
		q[0] = x->v[cl_origin(x, cyc[i].h0)].p[0];
		q[1] = x->v[cl_origin(x, cyc[i].h0)].p[1];
		for (j = 0; j < ncyc; j++)
		{
			if (!(cyc[j].area > 0) || cyc[j].comp == cyc[i].comp || cyc[j].area >= best_area)
				continue;
			if (q[0] < cyc[j].lo[0] || q[0] > cyc[j].hi[0] || q[1] < cyc[j].lo[1] || q[1] > cyc[j].hi[1])
				continue;
			if (cl_cycle_winding(x, cyc[j].h0, q) != 0)
			{
				best = j;
				best_area = cyc[j].area;
			}
		}
		cface[i] = best == CL_NONE ? 0 : cface[best];
	}
	for (i = 0; i < nh; i++)
		x->h[i].face = cface[hcyc[i]];
	free(vcomp);
	free(hcyc);
	*out_faces = nfaces;
	*out_cycle_face = cface;
	*out_cycles = cyc;
	*out_ncycles = ncyc;
	return QAWS_STATUS_OK;
}

static int cl_inside(qaws_clip_type ct, qaws_fill_rule fr, int ws, int wc)
{
	int s = qaws_fill_rule_inside(fr, ws), c = qaws_fill_rule_inside(fr, wc);
	switch (ct)
	{
	case QAWS_CLIP_INTERSECTION: return s && c;
	case QAWS_CLIP_UNION: return s || c;
	case QAWS_CLIP_DIFFERENCE: return s && !c;
	case QAWS_CLIP_XOR: return s != c;
	default: return 0;
	}
}

/* windings from the unbounded face (0, 0) across edges: crossing h's edge
   from its right to its left adds its steps when h runs along the edge */
static qaws_status cl_windings(cl_ctx* x, unsigned int nfaces, int* ws, int* wc, int* known)
{
	unsigned int nh = 2 * x->ne, i, *queue, *first, *list, qh = 0, qt = 0;
	queue = (unsigned int*)malloc(sizeof(unsigned int) * (nfaces + 1));
	first = (unsigned int*)calloc(nfaces + 2, sizeof(unsigned int));
	list = (unsigned int*)malloc(sizeof(unsigned int) * (nh ? nh : 1));
	if (!queue || !first || !list)
	{
		free(queue); free(first); free(list);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nh; i++) first[x->h[i].face + 1]++;
	for (i = 0; i < nfaces; i++) first[i + 1] += first[i];
	{
		unsigned int* fill = (unsigned int*)calloc(nfaces + 1, sizeof(unsigned int));
		if (!fill) { free(queue); free(first); free(list); return QAWS_STATUS_ALLOCATION_FAILURE; }
		for (i = 0; i < nh; i++)
			list[first[x->h[i].face] + fill[x->h[i].face]++] = i;
		free(fill);
	}
	for (i = 0; i < nfaces; i++) known[i] = 0;
	ws[0] = wc[0] = 0;
	known[0] = 1;
	queue[qt++] = 0;
	while (qh < qt)
	{
		unsigned int f = queue[qh++], k;
		for (k = first[f]; k < first[f + 1]; k++)
		{
			unsigned int h = list[k], g = x->h[h ^ 1].face;
			cl_edge const* e = &x->e[h >> 1];
			int sg = (h & 1) ? -1 : 1;
			if (known[g]) continue;
			ws[g] = ws[f] - sg * e->ds;
			wc[g] = wc[f] - sg * e->dc;
			known[g] = 1;
			queue[qt++] = g;
		}
	}
	free(queue); free(first); free(list);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  8. Output                                                               */
/* ======================================================================== */

static qaws_status cl_res_curve(qaws_clip_result* r, qaws_curve* c)
{
	CL_GROW(r->curves, r->curve_count, r->curve_cap, qaws_curve*);
	r->curves[r->curve_count++] = c;
	return QAWS_STATUS_OK;
}

static qaws_status cl_res_vertex(qaws_clip_result* r, qaws_clip_vertex const* v)
{
	CL_GROW(r->vertices, r->vertex_count, r->vertex_cap, qaws_clip_vertex);
	r->vertices[r->vertex_count++] = *v;
	return QAWS_STATUS_OK;
}

static void cl_vertex_record(cl_ctx const* x, unsigned int vi, qaws_clip_vertex* o)
{
	cl_vertex const* v = &x->v[vi];
	memset(o, 0, sizeof(*o));
	o->position.x = (qaws_scalar)v->p[0];
	o->position.y = (qaws_scalar)v->p[1];
	o->operand_a = o->path_a = o->curve_a = CL_NONE;
	o->operand_b = o->path_b = o->curve_b = CL_NONE;
	if (v->sa != CL_NONE)
	{
		o->operand_a = x->src[v->sa].operand; o->path_a = x->src[v->sa].path; o->curve_a = x->src[v->sa].curve;
		o->parameter_a = (qaws_scalar)(x->src[v->sa].base + x->src[v->sa].scale * v->ta);
	}
	if (v->sb != CL_NONE)
	{
		o->operand_b = x->src[v->sb].operand; o->path_b = x->src[v->sb].path; o->curve_b = x->src[v->sb].curve;
		o->parameter_b = (qaws_scalar)(x->src[v->sb].base + x->src[v->sb].scale * v->tb);
		if (x->d->z_fn)
			o->z = x->d->z_fn(x->d->z_user, o);
	}
}

static int cl_is_line(qaws_curve const* c)
{
	return c->degree == 1 && (c->kind == QAWS_CURVE_KIND_BSPLINE || c->kind == QAWS_CURVE_KIND_BEZIER);
}

/* the new curve's end control points moved onto the vertices, so pieces meet
   exactly (curves interpolating their end control points only) */
static qaws_status cl_snap(qaws_curve** c, double const* p0, double const* p1)
{
	qaws_curve* old = *c;
	qaws_scalar buf[2 * 64];
	qaws_scalar* cp = buf;
	unsigned int n = 0;
	qaws_status s;
	qaws_curve_kind k = old->kind;
	if (k != QAWS_CURVE_KIND_BEZIER && k != QAWS_CURVE_KIND_BSPLINE && k != QAWS_CURVE_KIND_NURBS && k != QAWS_CURVE_KIND_RATIONAL_BEZIER)
		return QAWS_STATUS_OK;
	if (qaws_curve_get_control_points(old, buf, 64, &n) == QAWS_STATUS_BUFFER_TOO_SMALL)
	{
		cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
		qaws_curve_get_control_points(old, cp, n, &n);
	}
	cp[0] = (qaws_scalar)p0[0]; cp[1] = (qaws_scalar)p0[1];
	cp[2 * n - 2] = (qaws_scalar)p1[0]; cp[2 * n - 1] = (qaws_scalar)p1[1];
	switch (k)
	{
	case QAWS_CURVE_KIND_BEZIER:
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = old->degree; d.control_points = cp; d.control_point_count = n;
		s = qaws_curve_create_bezier(&d, c);
		break;
	}
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	{
		qaws_rational_bezier_impl const* impl = (qaws_rational_bezier_impl const*)old->impl;
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = old->degree; d.control_points = cp; d.control_point_count = n;
		d.weights = impl->weights; d.weight_count = n;
		s = qaws_curve_create_rational_bezier(&d, c);
		break;
	}
	case QAWS_CURVE_KIND_BSPLINE:
	{
		qaws_bspline_impl const* impl = (qaws_bspline_impl const*)old->impl;
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = old->degree; d.control_points = cp; d.control_point_count = n;
		d.knots = impl->knots; d.knot_count = impl->knot_count;
		s = qaws_curve_create_bspline(&d, c);
		break;
	}
	default:
	{
		qaws_nurbs_impl const* impl = (qaws_nurbs_impl const*)old->impl;
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = old->degree; d.control_points = cp; d.control_point_count = n;
		d.knots = impl->knots; d.knot_count = impl->knot_count; d.weights = impl->weights; d.weight_count = n;
		s = qaws_curve_create_nurbs(&d, c);
		break;
	}
	}
	if (cp != buf) free(cp);
	if (s == QAWS_STATUS_OK)
		qaws_curve_destroy(old);
	else
		*c = old;
	return s;
}

/* a piece of a source as an output curve: (src, ta -> tb), ends on vertices */
typedef struct cl_piece
{
	unsigned int src;
	double ta, tb;
	unsigned int va, vb;
	unsigned int h;                      /* the half-edge (closed loops) */
	double area;                         /* its integral of x dy */
} cl_piece;

/* b within tol of the line through a and c (a spike back along it too) */
static int cl_collinear(double const* a, double const* b, double const* c, double tol)
{
	double dx = c[0] - a[0], dy = c[1] - a[1], l = hypot(dx, dy);
	double cr = (b[0] - a[0]) * dy - (b[1] - a[1]) * dx;
	if (!(l > 0))
		return 1;
	return fabs(cr) <= tol * l;
}
/* corners of a degree-1 source strictly between parameters ta and tb, in order */
static unsigned int cl_line_corners(qaws_curve const* c, double ta, double tb, double* out, unsigned int cap)
{
	qaws_scalar buf[2 * 64];
	qaws_scalar* cp = buf;
	qaws_scalar const* kn = NULL;
	unsigned int n = 0, i, m = 0;
	if (c->kind != QAWS_CURVE_KIND_BSPLINE)
		return 0;
	if (qaws_curve_get_control_points(c, buf, 64, &n) == QAWS_STATUS_BUFFER_TOO_SMALL)
	{
		cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		if (!cp) return 0;
		qaws_curve_get_control_points(c, cp, n, &n);
	}
	kn = ((qaws_bspline_impl const*)c->impl)->knots;
	/* degree 1: control point i sits at knot i + 1 */
	if (ta < tb)
	{
		for (i = 0; i < n && m < cap; i++)
			if (kn[i + 1] > ta && kn[i + 1] < tb)
			{
				out[2 * m] = cp[2 * i]; out[2 * m + 1] = cp[2 * i + 1];
				m++;
			}
	}
	else
	{
		for (i = n; i-- > 0 && m < cap;)
			if (kn[i + 1] > tb && kn[i + 1] < ta)
			{
				out[2 * m] = cp[2 * i]; out[2 * m + 1] = cp[2 * i + 1];
				m++;
			}
	}
	if (cp != buf) free(cp);
	return m;
}

/* the pieces of one output path as curves: runs of one source with
   contiguous parameters fused into one extract; runs of lines joined into
   one polyline (corners kept, or collinear ones dropped) */
static qaws_status cl_emit_pieces(cl_ctx* x, cl_piece const* pc, unsigned int n, int closed, struct cl_out_path* op)
{
	qaws_clip_result* r = x->r;
	unsigned int i = 0;
	int preserve = (x->d->flags & QAWS_CLIP_PRESERVE_COLLINEAR) != 0;
	qaws_status s = QAWS_STATUS_OK;
	op->first = r->curve_count;
	op->vfirst = r->vertex_count;
	while (i < n && s == QAWS_STATUS_OK)
	{
		qaws_curve const* c = x->src[pc[i].src].c;
		if (cl_is_line(c))
		{
			/* a polyline over every consecutive line piece */
			double* pts = NULL;
			unsigned int vstart = r->vertex_count;
			unsigned int np = 0, capp = 0, j = i;
			while (j < n && cl_is_line(x->src[pc[j].src].c))
			{
				double corners[2 * 256];
				unsigned int nc, k;
				qaws_clip_vertex vr;
				if (np + 2 + 256 > capp)
				{
					unsigned int nc2 = capp ? 2 * capp + 512 : 1024;
					double* g = (double*)realloc(pts, sizeof(double) * 2 * nc2);
					if (!g) { free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }
					pts = g; capp = nc2;
				}
				pts[2 * np] = x->v[pc[j].va].p[0]; pts[2 * np + 1] = x->v[pc[j].va].p[1];
				np++;
				cl_vertex_record(x, pc[j].va, &vr);
				if ((s = cl_res_vertex(r, &vr)) != QAWS_STATUS_OK) { free(pts); return s; }
				nc = cl_line_corners(x->src[pc[j].src].c, pc[j].ta, pc[j].tb, corners, 256);
				for (k = 0; k < nc; k++)
				{
					if (np + 2 > capp)
					{
						double* g = (double*)realloc(pts, sizeof(double) * 2 * (2 * capp + 16));
						if (!g) { free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }
						pts = g; capp = 2 * capp + 16;
					}
					pts[2 * np] = corners[2 * k]; pts[2 * np + 1] = corners[2 * k + 1];
					np++;
					memset(&vr, 0, sizeof(vr));
					vr.position.x = (qaws_scalar)corners[2 * k]; vr.position.y = (qaws_scalar)corners[2 * k + 1];
					vr.operand_a = x->src[pc[j].src].operand; vr.path_a = x->src[pc[j].src].path; vr.curve_a = x->src[pc[j].src].curve;
					vr.operand_b = vr.path_b = vr.curve_b = CL_NONE;
					if ((s = cl_res_vertex(r, &vr)) != QAWS_STATUS_OK) { free(pts); return s; }
				}
				j++;
			}
			/* the end: the start of the next piece, or the loop's start */
			{
				unsigned int vend = pc[j - 1].vb;
				int whole = closed && i == 0 && j == n;
				qaws_clip_vertex endr;
				if (!whole)
				{
					pts[2 * np] = x->v[vend].p[0]; pts[2 * np + 1] = x->v[vend].p[1];
					np++;
					cl_vertex_record(x, vend, &endr);
					if ((s = cl_res_vertex(r, &endr)) != QAWS_STATUS_OK) { free(pts); return s; }
				}
				if (!preserve)
				{
					/* drop corners between collinear neighbours (and spikes) */
					unsigned int a, m = 0, changed = 1;
					while (changed && np >= 3)
					{
						changed = 0;
						m = 0;
						for (a = 0; a < np; a++)
						{
							int ends = !whole && (a == 0 || a == np - 1);
							double const* pp = whole ? &pts[2 * ((a + np - 1) % np)] : (a ? &pts[2 * (a - 1)] : NULL);
							double const* nn = whole ? &pts[2 * ((a + 1) % np)] : (a + 1 < np ? &pts[2 * (a + 1)] : NULL);
							if (!ends && pp && nn && cl_collinear(pp, &pts[2 * a], nn, x->tol))
							{
								changed = 1;
								continue;
							}
							pts[2 * m] = pts[2 * a]; pts[2 * m + 1] = pts[2 * a + 1];
							r->vertices[vstart + m] = r->vertices[vstart + a];
							m++;
						}
						np = m;
					}
				}
				/* the end vertex belongs to the next piece (or, open, to the caller) */
				r->vertex_count = vstart + np - (whole ? 0u : 1u);
				if (np >= 2)
				{
					qaws_scalar* sp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * np);
					qaws_curve* pl = NULL;
					unsigned int a;
					if (!sp) { free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }
					for (a = 0; a < 2 * np; a++) sp[a] = (qaws_scalar)pts[a];
					s = qaws_curve_create_polyline_2d(sp, np, whole, &pl);
					free(sp);
					if (s == QAWS_STATUS_OK)
						s = cl_res_curve(r, pl);
				}
			}
			free(pts);
			i = j;
			continue;
		}
		/* a curved run: one extract over contiguous pieces of one source */
		{
			unsigned int j = i + 1;
			qaws_curve* e = NULL;
			qaws_clip_vertex vr;
			while (j < n && pc[j].src == pc[i].src && pc[j].ta == pc[j - 1].tb &&
				((pc[i].tb > pc[i].ta) == (pc[j].tb > pc[j].ta)))
				j++;
			cl_vertex_record(x, pc[i].va, &vr);
			s = cl_res_vertex(r, &vr);
			if (s == QAWS_STATUS_OK)
				s = qaws_curve_extract(c, (qaws_scalar)pc[i].ta, (qaws_scalar)pc[j - 1].tb, &e);
			if (s == QAWS_STATUS_OK)
				s = cl_snap(&e, x->v[pc[i].va].p, x->v[pc[j - 1].vb].p);
			if (s == QAWS_STATUS_OK)
				s = cl_res_curve(r, e);
			else
				qaws_curve_destroy(e);
			i = j;
		}
	}
	op->count = r->curve_count - op->first;
	op->vcount = r->vertex_count - op->vfirst;
	return s;
}

static qaws_status cl_push_path(qaws_clip_result* r, int open, struct cl_out_path const* op)
{
	if (open)
	{
		CL_GROW(r->open_paths, r->open_count, r->open_cap, struct cl_out_path);
		r->open_paths[r->open_count++] = *op;
	}
	else
	{
		CL_GROW(r->paths, r->path_count, r->path_cap, struct cl_out_path);
		r->paths[r->path_count++] = *op;
	}
	return QAWS_STATUS_OK;
}

static qaws_status cl_output(cl_ctx* x, unsigned int nfaces, int const* inside)
{
	unsigned int nh = 2 * x->ne, i, nloop = 0;
	unsigned int* loop_of = (unsigned int*)malloc(sizeof(unsigned int) * (nh ? nh : 1));
	unsigned int* fcomp = (unsigned int*)malloc(sizeof(unsigned int) * (nfaces + 1));
	unsigned int* comp_outer = (unsigned int*)malloc(sizeof(unsigned int) * (nfaces + 1));
	unsigned int* comp_hole = (unsigned int*)malloc(sizeof(unsigned int) * (nfaces + 1));
	cl_piece* pc = NULL;
	unsigned int cappc = 0, capst = 0, *vpos = NULL;
	cl_piece* st = NULL;
	unsigned int* loop_left = NULL;
	unsigned int* loop_right = NULL;
	double* loop_area = NULL;
	unsigned int first_path = x->r->path_count;
	qaws_status s = QAWS_STATUS_OK;
	if (!loop_of || !fcomp || !comp_outer || !comp_hole)
	{
		s = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < nh; i++)
	{
		x->h[i].out = inside[x->h[i].face] && !inside[x->h[i ^ 1].face];
		loop_of[i] = CL_NONE;
	}
	/* face components: faces joined across edges that are not output */
	for (i = 0; i < nfaces; i++) fcomp[i] = i;
	for (i = 0; i < nh; i += 2)
		if (inside[x->h[i].face] == inside[x->h[i + 1].face])
		{
			unsigned int a = cl_vfind(fcomp, x->h[i].face), b = cl_vfind(fcomp, x->h[i + 1].face);
			if (a != b) fcomp[a] = b;
		}
	for (i = 0; i < nfaces; i++)
		comp_outer[i] = comp_hole[i] = CL_NONE;
	/* loops: from an output half-edge, the next output one clockwise round
	   its target from its twin; a loop through one vertex twice (two parts
	   touching at a point) is cut there into simple loops */
	vpos = (unsigned int*)malloc(sizeof(unsigned int) * (x->nv + 1));
	if (!vpos) { s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
	for (i = 0; i < x->nv; i++) vpos[i] = CL_NONE;
	for (i = 0; i < nh && s == QAWS_STATUS_OK; i++)
	{
		unsigned int h = i, np = 0, top = 0, k;
		if (!x->h[i].out || loop_of[i] != CL_NONE)
			continue;
		do
		{
			cl_edge const* e = &x->e[h >> 1];
			cl_vertex const* v;
			unsigned int tw, p, g = CL_NONE;
			loop_of[h] = 1;
			if (np == cappc)
			{
				unsigned int nc = cappc ? 2 * cappc : 64;
				cl_piece* gp = (cl_piece*)realloc(pc, sizeof(cl_piece) * nc);
				if (!gp) { s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
				pc = gp; cappc = nc;
			}
			pc[np].src = e->src;
			pc[np].ta = (h & 1) ? e->t1 : e->t0;
			pc[np].tb = (h & 1) ? e->t0 : e->t1;
			pc[np].va = cl_origin(x, h);
			pc[np].vb = cl_target(x, h);
			pc[np].h = h;
			pc[np].area = cl_half_area(x, h);
			np++;
			tw = h ^ 1;
			v = &x->v[cl_origin(x, tw)];
			p = x->h[tw].pos;
			for (k = 1; k <= v->count; k++)
			{
				unsigned int cand = x->vout[v->first + (p + v->count - k) % v->count];
				if (x->h[cand].out) { g = cand; break; }
			}
			if (g == CL_NONE) { s = QAWS_STATUS_INTERNAL_ERROR; goto done; }
			h = g;
		} while (h != i && loop_of[h] == CL_NONE);
		/* cut at repeated vertices: pieces are pushed in order; when one
		   ends at a vertex some earlier piece on the stack starts from, the
		   pieces from there on are a simple loop */
		if (np > capst)
		{
			cl_piece* g = (cl_piece*)realloc(st, sizeof(cl_piece) * np);
			if (!g) { s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
			st = g; capst = np;
		}
		for (k = 0; k < np && s == QAWS_STATUS_OK; k++)
		{
			unsigned int from;
			if (vpos[pc[k].va] == CL_NONE)
				vpos[pc[k].va] = top;
			st[top++] = pc[k];
			from = vpos[pc[k].vb];
			if (from == CL_NONE)
				continue;
			{
				unsigned int n = top - from, a;
				double area = 0.0;
				struct cl_out_path op;
				cl_piece* lp = &st[from];
				for (a = 0; a < n; a++)
				{
					area += lp[a].area;
					vpos[lp[a].va] = CL_NONE;
				}
				{
					void* g1 = realloc(loop_left, sizeof(unsigned int) * (nloop + 1));
					void* g2 = g1 ? realloc(loop_right, sizeof(unsigned int) * (nloop + 1)) : NULL;
					void* g3 = g2 ? realloc(loop_area, sizeof(double) * (nloop + 1)) : NULL;
					if (g1) loop_left = (unsigned int*)g1;
					if (g2) loop_right = (unsigned int*)g2;
					if (g3) loop_area = (double*)g3;
					if (!g3) { s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
				}
				loop_left[nloop] = cl_vfind(fcomp, x->h[lp[0].h].face);
				loop_right[nloop] = cl_vfind(fcomp, x->h[lp[0].h ^ 1].face);
				loop_area[nloop] = area;
				if (area > 0)
					comp_outer[loop_left[nloop]] = nloop;
				else
					comp_hole[loop_right[nloop]] = nloop;
				if (x->d->flags & QAWS_CLIP_REVERSE_SOLUTION)
				{
					/* traverse backwards */
					for (a = 0; a < n / 2; a++)
					{
						cl_piece t = lp[a];
						lp[a] = lp[n - 1 - a];
						lp[n - 1 - a] = t;
					}
					for (a = 0; a < n; a++)
					{
						double t = lp[a].ta; unsigned int vv = lp[a].va;
						lp[a].ta = lp[a].tb; lp[a].tb = t;
						lp[a].va = lp[a].vb; lp[a].vb = vv;
					}
				}
				memset(&op, 0, sizeof(op));
				op.hole = area < 0;
				s = cl_emit_pieces(x, lp, n, 1, &op);
				if (s == QAWS_STATUS_OK)
					s = cl_push_path(x->r, 0, &op);
				nloop++;
				top = from;
			}
		}
	}
	/* nesting: a hole's parent is the outer loop of the inside part on its
	   left; an outer loop's parent is the hole around the outside part on
	   its right */
	for (i = 0; s == QAWS_STATUS_OK && i < nloop; i++)
	{
		struct cl_out_path* op = &x->r->paths[first_path + i];
		unsigned int par = loop_area[i] > 0 ? comp_hole[loop_right[i]] : comp_outer[loop_left[i]];
		op->parent = par == CL_NONE ? CL_NONE : first_path + par;
	}
	for (i = 0; s == QAWS_STATUS_OK && i < nloop; i++)
	{
		struct cl_out_path* op = &x->r->paths[first_path + i];
		unsigned int d = 0, p = op->parent;
		while (p != CL_NONE && d <= nloop)
		{
			d++;
			p = x->r->paths[p].parent;
		}
		op->depth = d;
		op->hole = (d & 1) != 0;
	}
done:
	free(loop_of); free(fcomp); free(comp_outer); free(comp_hole);
	free(pc); free(loop_left); free(loop_right); free(loop_area); free(st); free(vpos);
	return s;
}

/* ======================================================================== */
/*  9. Open subjects                                                        */
/* ======================================================================== */

static qaws_status cl_open(cl_ctx* x)
{
	qaws_clip_desc const* d = x->d;
	unsigned int i = 0;
	cl_piece* pc = NULL;
	unsigned int np = 0, cappc = 0;
	qaws_status s = QAWS_STATUS_OK;
	unsigned int last_path = CL_NONE, last_v = CL_NONE;
	/* open edges are in source order, each source's in parameter order */
	for (i = 0; i <= x->noe && s == QAWS_STATUS_OK; i++)
	{
		int keep = 0;
		cl_edge const* e = i < x->noe ? &x->oe[i] : NULL;
		if (e)
		{
			double m[2];
			qaws_vec2 q;
			qaws_point_location ls = QAWS_POINT_OUTSIDE, lc = QAWS_POINT_OUTSIDE;
			int in_s, in_c;
			cl_eval(x->src[e->src].c, 0.5 * (e->t0 + e->t1), m, NULL);
			q.x = (qaws_scalar)m[0]; q.y = (qaws_scalar)m[1];
			if (d->subject_count)
				s = qaws_region_locate_point_2d(d->subjects, d->subject_count, d->fill_rule, q, (qaws_scalar)x->tol, &ls);
			if (s == QAWS_STATUS_OK && d->clip_count)
				s = qaws_region_locate_point_2d(d->clips, d->clip_count, d->fill_rule, q, (qaws_scalar)x->tol, &lc);
			in_s = ls == QAWS_POINT_INSIDE;
			in_c = lc == QAWS_POINT_INSIDE;
			switch (d->clip_type)
			{
			case QAWS_CLIP_INTERSECTION: keep = in_c; break;
			case QAWS_CLIP_UNION: keep = !in_s && !in_c; break;
			case QAWS_CLIP_NONE: keep = 0; break;
			default: keep = !in_c; break;
			}
		}
		/* close the current run when this piece does not continue it */
		if (np && (!keep || x->src[e->src].path != last_path || e->v0 != last_v))
		{
			struct cl_out_path op;
			memset(&op, 0, sizeof(op));
			op.parent = CL_NONE;
			s = cl_emit_pieces(x, pc, np, 0, &op);
			if (s == QAWS_STATUS_OK)
			{
				qaws_clip_vertex vr;
				/* the run's end vertex */
				cl_vertex_record(x, pc[np - 1].vb, &vr);
				s = cl_res_vertex(x->r, &vr);
				op.vcount++;
			}
			if (s == QAWS_STATUS_OK)
				s = cl_push_path(x->r, 1, &op);
			np = 0;
		}
		if (keep)
		{
			if (np == cappc)
			{
				unsigned int nc = cappc ? 2 * cappc : 64;
				cl_piece* g = (cl_piece*)realloc(pc, sizeof(cl_piece) * nc);
				if (!g) { free(pc); return QAWS_STATUS_ALLOCATION_FAILURE; }
				pc = g; cappc = nc;
			}
			pc[np].src = e->src; pc[np].ta = e->t0; pc[np].tb = e->t1; pc[np].va = e->v0; pc[np].vb = e->v1;
			np++;
			last_path = x->src[e->src].path;
			last_v = e->v1;
		}
	}
	free(pc);
	return s;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

static void cl_free_ctx(cl_ctx* x)
{
	unsigned int i;
	for (i = 0; i < x->ntemp; i++) qaws_curve_destroy(x->temp[i]);
	free(x->temp); free(x->src); free(x->pts); free(x->v); free(x->e); free(x->oe);
	free(x->h); free(x->vout);
}

qaws_status qaws_clip_execute(qaws_clip_desc const* desc, qaws_clip_result** out_result)
{
	cl_ctx x;
	qaws_status s;
	unsigned int nfaces = 0, *cface = NULL, ncyc = 0;
	cl_cycle* cyc = NULL;
	int *ws = NULL, *wc = NULL, *known = NULL, *inside = NULL;
	unsigned int i, npaths;
	if (!desc || !out_result || (desc->subject_count && !desc->subjects) || (desc->clip_count && !desc->clips) ||
		(desc->open_subject_count && !desc->open_subjects))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_result = NULL;
	memset(&x, 0, sizeof(x));
	x.d = desc;
	x.r = (qaws_clip_result*)calloc(1, sizeof(qaws_clip_result));
	npaths = desc->subject_count + desc->clip_count;
	(void)npaths;
	if (!x.r)
	{
		free(x.r);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* the extent and the tolerance, from the curves' sampling */
	{
		qaws_curve const** all;
		unsigned int n = 0, k, j;
		for (k = 0; k < desc->subject_count; k++) n += desc->subjects[k].curve_count;
		for (k = 0; k < desc->clip_count; k++) n += desc->clips[k].curve_count;
		for (k = 0; k < desc->open_subject_count; k++) n += desc->open_subjects[k].curve_count;
		all = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (n ? n : 1));
		if (!all) { cl_free_ctx(&x); free(x.r); return QAWS_STATUS_ALLOCATION_FAILURE; }
		n = 0;
		for (k = 0; k < desc->subject_count; k++)
			for (j = 0; j < desc->subjects[k].curve_count; j++) if (desc->subjects[k].curves[j]) all[n++] = desc->subjects[k].curves[j];
		for (k = 0; k < desc->clip_count; k++)
			for (j = 0; j < desc->clips[k].curve_count; j++) if (desc->clips[k].curves[j]) all[n++] = desc->clips[k].curves[j];
		for (k = 0; k < desc->open_subject_count; k++)
			for (j = 0; j < desc->open_subjects[k].curve_count; j++) if (desc->open_subjects[k].curves[j]) all[n++] = desc->open_subjects[k].curves[j];
		x.ext = n ? (double)qaws_internal_flatten_extent(all, n, 2, NULL, 0) : 1.0;
		free((void*)all);
		if (!(x.ext > 0)) x.ext = 1.0;
		x.tol = desc->tolerance > 0 ? (double)desc->tolerance : x.ext * CL_TOL_REL;
	}
	s = cl_add_paths(&x, desc->subjects, desc->subject_count, 0);
	if (s == QAWS_STATUS_OK) s = cl_add_paths(&x, desc->clips, desc->clip_count, 1);
	if (s == QAWS_STATUS_OK) s = cl_add_paths(&x, desc->open_subjects, desc->open_subject_count, 2);
	/* lines only: every hit is computed in double from the end points, one
	   rounding away, so vertices need only a rounding-sized tolerance */
	if (s == QAWS_STATUS_OK && !(desc->tolerance > 0))
	{
		unsigned int lines = 0;
		for (i = 0; i < x.nsrc; i++)
			lines += cl_is_line_src(&x, i);
		if (lines == x.nsrc)
			x.tol = x.ext * 1e-11;
	}
	if (s == QAWS_STATUS_OK && desc->clip_type != QAWS_CLIP_NONE && x.nsrc)
	{
		s = cl_hits(&x);
		if (s == QAWS_STATUS_OK) s = cl_edges(&x);
		if (s == QAWS_STATUS_OK) s = cl_merge_edges(&x);
		if (s == QAWS_STATUS_OK) s = cl_graph(&x);
		if (s == QAWS_STATUS_OK) s = cl_faces(&x, &nfaces, &cface, &cyc, &ncyc);
		if (s == QAWS_STATUS_OK)
		{
			ws = (int*)malloc(sizeof(int) * (nfaces + 1));
			wc = (int*)malloc(sizeof(int) * (nfaces + 1));
			known = (int*)malloc(sizeof(int) * (nfaces + 1));
			inside = (int*)malloc(sizeof(int) * (nfaces + 1));
			if (!ws || !wc || !known || !inside)
				s = QAWS_STATUS_ALLOCATION_FAILURE;
		}
		if (s == QAWS_STATUS_OK) s = cl_windings(&x, nfaces, ws, wc, known);
		if (s == QAWS_STATUS_OK)
		{
			for (i = 0; i < nfaces; i++)
				inside[i] = known[i] && cl_inside(desc->clip_type, desc->fill_rule, ws[i], wc[i]);
			s = cl_output(&x, nfaces, inside);
		}
		if (s == QAWS_STATUS_OK && x.noe)
			s = cl_open(&x);
	}
	free(cface); free(cyc); free(ws); free(wc); free(known); free(inside);
	cl_free_ctx(&x);
	if (s != QAWS_STATUS_OK)
	{
		qaws_clip_result_destroy(x.r);
		return s;
	}
	/* const views of the curves for the path views */
	x.r->views = (qaws_curve const**)malloc(sizeof(qaws_curve*) * (x.r->curve_count ? x.r->curve_count : 1));
	if (!x.r->views)
	{
		qaws_clip_result_destroy(x.r);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < x.r->curve_count; i++)
		x.r->views[i] = x.r->curves[i];
	*out_result = x.r;
	return QAWS_STATUS_OK;
}

void qaws_clip_result_destroy(qaws_clip_result* r)
{
	unsigned int i;
	if (!r)
		return;
	for (i = 0; i < r->curve_count; i++)
		qaws_curve_destroy(r->curves[i]);
	free(r->curves);
	free((void*)r->views);
	free(r->paths);
	free(r->open_paths);
	free(r->vertices);
	free(r);
}

unsigned int qaws_clip_result_get_path_count(qaws_clip_result const* r)
{
	return r ? r->path_count : 0;
}

qaws_status qaws_clip_result_get_path(qaws_clip_result const* r, unsigned int i, qaws_path_2d* out)
{
	if (!r || !out || i >= r->path_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out->curves = r->views + r->paths[i].first;
	out->curve_count = r->paths[i].count;
	out->closed = 1;
	return QAWS_STATUS_OK;
}

unsigned int qaws_clip_result_get_parent(qaws_clip_result const* r, unsigned int i)
{
	return r && i < r->path_count ? r->paths[i].parent : CL_NONE;
}

int qaws_clip_result_is_hole(qaws_clip_result const* r, unsigned int i)
{
	return r && i < r->path_count ? r->paths[i].hole : 0;
}

unsigned int qaws_clip_result_get_depth(qaws_clip_result const* r, unsigned int i)
{
	return r && i < r->path_count ? r->paths[i].depth : 0;
}

unsigned int qaws_clip_result_get_open_path_count(qaws_clip_result const* r)
{
	return r ? r->open_count : 0;
}

qaws_status qaws_clip_result_get_open_path(qaws_clip_result const* r, unsigned int i, qaws_path_2d* out)
{
	if (!r || !out || i >= r->open_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out->curves = r->views + r->open_paths[i].first;
	out->curve_count = r->open_paths[i].count;
	out->closed = 0;
	return QAWS_STATUS_OK;
}

qaws_status qaws_clip_result_get_vertices(qaws_clip_result const* r, int open, unsigned int i,
	qaws_clip_vertex const** out_vertices, unsigned int* out_count)
{
	struct cl_out_path const* p;
	if (!r || !out_vertices || !out_count || i >= (open ? r->open_count : r->path_count))
		return QAWS_STATUS_INVALID_ARGUMENT;
	p = open ? &r->open_paths[i] : &r->paths[i];
	*out_vertices = r->vertices + p->vfirst;
	*out_count = p->vcount;
	return QAWS_STATUS_OK;
}

qaws_status qaws_clip_boolean(qaws_clip_type clip_type, qaws_fill_rule fill_rule,
	qaws_path_2d const* subjects, unsigned int subject_count,
	qaws_path_2d const* clips, unsigned int clip_count,
	qaws_clip_result** out_result)
{
	qaws_clip_desc d;
	memset(&d, 0, sizeof(d));
	d.subjects = subjects;
	d.subject_count = subject_count;
	d.clips = clips;
	d.clip_count = clip_count;
	d.clip_type = clip_type;
	d.fill_rule = fill_rule;
	return qaws_clip_execute(&d, out_result);
}
