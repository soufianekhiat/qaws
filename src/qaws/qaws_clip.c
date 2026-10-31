/*
 * Boolean operations on regions of curves: a planar arrangement of the input
 * curves, faces with subject / clip winding numbers, Clipper2's truth tables.
 *
 *   1. sources: every curve of every path (an unclosed closed-operand path
 *      gets a closing line)
 *   2. hits: one batched intersection of all sources, self-intersections
 *      included; crossings, touches and shared stretches
 *   3. vertices: curve ends and hits, clustered within the tolerance; a
 *      vertex within the tolerance of a line cuts it there (snap rounding)
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
#include "internal/qaws_internal_memory.h"
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
#include "internal/qaws_internal_clip.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* split by tools/splitc.py: shared between qaws_clip.c and its moved part */
#define cl_eval QAWS_SFX(qaws_internal_clip_cl_eval)
#define cl_is_line_src QAWS_SFX(qaws_internal_clip_cl_is_line_src)
#define cl_origin QAWS_SFX(qaws_internal_clip_cl_origin)
#define cl_target QAWS_SFX(qaws_internal_clip_cl_target)
#define cl_half_area QAWS_SFX(qaws_internal_clip_cl_half_area)
#define cl_vfind QAWS_SFX(qaws_internal_clip_cl_vfind)
#define cl_res_vertex QAWS_SFX(qaws_internal_clip_cl_res_vertex)
#define cl_vertex_record QAWS_SFX(qaws_internal_clip_cl_vertex_record)
#define cl_emit_pieces QAWS_SFX(qaws_internal_clip_cl_emit_pieces)
#define cl_push_path QAWS_SFX(qaws_internal_clip_cl_push_path)
#define cl_output QAWS_SFX(qaws_internal_clip_cl_output)
#define cl_free_ctx QAWS_SFX(qaws_internal_clip_cl_free_ctx)

#if QAWS_SCALAR_IS_FLOAT
#define CL_TOL_REL 1e-4
#define CL_AREA_REL 1e-6
#define CL_ANG_TIE 5e-3     /* tangent angles this close are a tie: a tangent contact's noise, the root of 2e-5 */
#define CL_TANGENT_SIN 1e-2 /* hits this tangent may lie a tangent contact's reach off */
#define CL_PROBE_REL 0.25   /* tie probes this far along the shorter edge: bending beats rounding */
#else
#define CL_TOL_REL 6.4e-9
#define CL_AREA_REL 1e-12
#define CL_ANG_TIE 1e-5     /* a tangent contact's angles carry the square root of the hits' 1e-10 */
#define CL_TANGENT_SIN 1e-3
#define CL_PROBE_REL 1e-2
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
	double area;                         /* integral of x dy along it (NAN until computed) */
	unsigned int next, face, pos;
	int out, seen;
} cl_half;

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
	do { if ((n) >= 0x7FFFFFFFu) { return QAWS_STATUS_ALLOCATION_FAILURE; } /* doubling the capacity stays in unsigned int */ \
		if ((n) == (cap)) { unsigned int nc_ = (cap) ? 2 * (cap) : 64; type* np_ = (type*)qaws_internal_realloc((ptr), qaws_internal_mul2(sizeof(type), nc_)); \
		if (!np_) { return QAWS_STATUS_ALLOCATION_FAILURE; } (ptr) = np_; (cap) = nc_; } } while (0)

qaws_status cl_eval(qaws_curve const* c, double t, double* p, double* d)
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
	s = qaws_curve_evaluate_2d(&r, c, (qaws_scalar)t, d ? QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 : QAWS_EVAL_FLAG_POSITION);
	p[0] = r.position.x; p[1] = r.position.y;
	if (d) { d[0] = r.d1.x; d[1] = r.d1.y; }
	return s;
}

int cl_is_line_src(cl_ctx const* x, unsigned int i)
{
	qaws_curve const* c = x->src[i].c;
	return c->kind == QAWS_CURVE_KIND_BEZIER && c->degree == 1;
}

/* ======================================================================== */
/*  5. Graph                                                                */
/* ======================================================================== */

unsigned int cl_origin(cl_ctx const* x, unsigned int h)
{
	cl_edge const* e = &x->e[h >> 1];
	return (h & 1) ? e->v1 : e->v0;
}

unsigned int cl_target(cl_ctx const* x, unsigned int h)
{
	cl_edge const* e = &x->e[h >> 1];
	return (h & 1) ? e->v0 : e->v1;
}

/* ======================================================================== */
/*  6-7. Faces and windings                                                 */
/* ======================================================================== */

/* integral of x dy along half-edge h */
double cl_half_area(cl_ctx const* x, unsigned int h)
{
	cl_edge const* e = &x->e[h >> 1];
	double a = 0.0;
	if (cl_is_line_src(x, e->src))
	{
		/* a line from p to q: (x_p + x_q) / 2 (y_q - y_p), in double, on the
		   vertices (the points the graph is built on) */
		double const* p = x->v[cl_origin(x, h)].p;
		double const* q = x->v[cl_origin(x, h ^ 1)].p;
		return 0.5 * (p[0] + q[0]) * (q[1] - p[1]);
	}
	qaws_internal_curve_area_2d(&a, x->src[e->src].c, (h & 1) ? e->t1 : e->t0,
		(h & 1) ? e->t0 : e->t1, x->ext * x->ext * CL_AREA_REL);
	/* the curve's ends lie within the tolerance of the vertices: chords from
	   and to the vertices close it, else each gap weighs |x| times its size */
	{
		double const* p = x->v[cl_origin(x, h)].p;
		double const* q = x->v[cl_origin(x, h ^ 1)].p;
		double c0[2], c1[2];
		cl_eval(x->src[e->src].c, (h & 1) ? e->t1 : e->t0, c0, NULL);
		cl_eval(x->src[e->src].c, (h & 1) ? e->t0 : e->t1, c1, NULL);
		a += 0.5 * (p[0] + c0[0]) * (c0[1] - p[1]) + 0.5 * (c1[0] + q[0]) * (q[1] - c1[1]);
	}
	return a;
}

unsigned int cl_vfind(unsigned int* p, unsigned int i)
{
	while (p[i] != i) { p[i] = p[p[i]]; i = p[i]; }
	return i;
}

/* ======================================================================== */
/*  8. Output                                                               */
/* ======================================================================== */

static qaws_status cl_res_room(qaws_clip_result* r)
{
	CL_GROW(r->curves, r->curve_count, r->curve_cap, qaws_curve*);
	return QAWS_STATUS_OK;
}

/* c held by the result; destroyed when there is no room to hold it (it would leak) */
static qaws_status cl_res_curve(qaws_clip_result* r, qaws_curve* c)
{
	qaws_status s = cl_res_room(r);
	if (s != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(c);
		return s;
	}
	r->curves[r->curve_count++] = c;
	return QAWS_STATUS_OK;
}

qaws_status cl_res_vertex(qaws_clip_result* r, qaws_clip_vertex const* v)
{
	CL_GROW(r->vertices, r->vertex_count, r->vertex_cap, qaws_clip_vertex);
	r->vertices[r->vertex_count++] = *v;
	return QAWS_STATUS_OK;
}

void cl_vertex_record(cl_ctx const* x, unsigned int vi, qaws_clip_vertex* o)
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
	if (qaws_curve_get_control_points(buf, &n, old, 64) == QAWS_STATUS_BUFFER_TOO_SMALL)
	{
		cp = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 2, n));
		if (!cp) return QAWS_STATUS_ALLOCATION_FAILURE;
		qaws_curve_get_control_points(cp, &n, old, n);
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
		s = qaws_curve_create_bezier(c, &d);
		break;
	}
	case QAWS_CURVE_KIND_RATIONAL_BEZIER:
	{
		qaws_rational_bezier_impl const* impl = (qaws_rational_bezier_impl const*)old->impl;
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = old->degree; d.control_points = cp; d.control_point_count = n;
		d.weights = impl->weights; d.weight_count = n;
		s = qaws_curve_create_rational_bezier(c, &d);
		break;
	}
	case QAWS_CURVE_KIND_BSPLINE:
	{
		qaws_bspline_impl const* impl = (qaws_bspline_impl const*)old->impl;
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = old->degree; d.control_points = cp; d.control_point_count = n;
		d.knots = impl->knots; d.knot_count = impl->knot_count;
		s = qaws_curve_create_bspline(c, &d);
		break;
	}
	default:
	{
		qaws_nurbs_impl const* impl = (qaws_nurbs_impl const*)old->impl;
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = old->degree; d.control_points = cp; d.control_point_count = n;
		d.knots = impl->knots; d.knot_count = impl->knot_count; d.weights = impl->weights; d.weight_count = n;
		s = qaws_curve_create_nurbs(c, &d);
		break;
	}
	}
	if (cp != buf) qaws_internal_free(cp);
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
	if (qaws_curve_get_control_points(buf, &n, c, 64) == QAWS_STATUS_BUFFER_TOO_SMALL)
	{
		cp = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 2, n));
		if (!cp) return 0;
		qaws_curve_get_control_points(cp, &n, c, n);
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
	if (cp != buf) qaws_internal_free(cp);
	return m;
}

/* the pieces of one output path as curves: runs of one source with
   contiguous parameters fused into one extract; runs of lines joined into
   one polyline (corners kept, or collinear ones dropped) */
qaws_status cl_emit_pieces(cl_ctx* x, cl_piece const* pc, unsigned int n, int closed, struct cl_out_path* op)
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
					double* g = (double*)qaws_internal_realloc(pts, qaws_internal_mul3(sizeof(double), 2, nc2));
					if (!g) { qaws_internal_free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }
					pts = g; capp = nc2;
				}
				pts[2 * np] = x->v[pc[j].va].p[0]; pts[2 * np + 1] = x->v[pc[j].va].p[1];
				np++;
				cl_vertex_record(x, pc[j].va, &vr);
				if ((s = cl_res_vertex(r, &vr)) != QAWS_STATUS_OK) { qaws_internal_free(pts); return s; }
				nc = cl_line_corners(x->src[pc[j].src].c, pc[j].ta, pc[j].tb, corners, 256);
				for (k = 0; k < nc; k++)
				{
					if ((size_t)np + 2 > capp)
					{
						unsigned int ncap = qaws_internal_grow(capp, (size_t)np + 2, 16);
						double* g = ncap ? (double*)qaws_internal_realloc(pts, qaws_internal_mul3(sizeof(double), 2, ncap)) : NULL;
						if (!g) { qaws_internal_free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }
						pts = g; capp = ncap;
					}
					pts[2 * np] = corners[2 * k]; pts[2 * np + 1] = corners[2 * k + 1];
					np++;
					memset(&vr, 0, sizeof(vr));
					vr.position.x = (qaws_scalar)corners[2 * k]; vr.position.y = (qaws_scalar)corners[2 * k + 1];
					vr.operand_a = x->src[pc[j].src].operand; vr.path_a = x->src[pc[j].src].path; vr.curve_a = x->src[pc[j].src].curve;
					vr.operand_b = vr.path_b = vr.curve_b = CL_NONE;
					if ((s = cl_res_vertex(r, &vr)) != QAWS_STATUS_OK) { qaws_internal_free(pts); return s; }
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
					if ((s = cl_res_vertex(r, &endr)) != QAWS_STATUS_OK) { qaws_internal_free(pts); return s; }
				}
				if (!preserve)
				{
					/* drop corners between collinear neighbours (and spikes): each
					   with the run dropped since the last corner kept against the
					   chord from it, so a fine polygon of a curve thins out to the
					   tolerance instead of vanishing at once */
					unsigned int a, m = 0, changed = 1, run;
					while (changed && np >= 3)
					{
						changed = 0;
						m = 0;
						run = 0;
						for (a = 0; a < np; a++)
						{
							int ends = !whole && (a == 0 || a == np - 1), drop;
							double const* pp = m ? &pts[2 * (m - 1)] : (whole ? &pts[2 * (np - 1)] : NULL);
							double const* nn = whole ? &pts[2 * ((a + 1) % np)] : (a + 1 < np ? &pts[2 * (a + 1)] : NULL);
							unsigned int q;
							/* the run lies past the corners kept: not yet overwritten */
							drop = !ends && pp && nn && run < 256;
							for (q = a - run; drop && q <= a; q++)
								drop = cl_collinear(pp, &pts[2 * q], nn, x->tol);
							if (drop)
							{
								changed = 1;
								run++;
								continue;
							}
							run = 0;
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
					qaws_scalar* sp = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 2, np));
					qaws_curve* pl = NULL;
					unsigned int a;
					if (!sp) { qaws_internal_free(pts); return QAWS_STATUS_ALLOCATION_FAILURE; }
					for (a = 0; a < 2 * np; a++) sp[a] = (qaws_scalar)pts[a];
					s = qaws_curve_create_polyline_2d(&pl, sp, np, whole);
					qaws_internal_free(sp);
					if (s == QAWS_STATUS_OK)
						s = cl_res_curve(r, pl);
				}
			}
			qaws_internal_free(pts);
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
				s = qaws_curve_extract(&e, c, (qaws_scalar)pc[i].ta, (qaws_scalar)pc[j - 1].tb);
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

qaws_status cl_push_path(qaws_clip_result* r, int open, struct cl_out_path const* op)
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

qaws_status cl_output(cl_ctx* x, unsigned int nfaces, int const* inside)
{
	unsigned int nh = 2 * x->ne, i, nloop = 0;
	unsigned int* loop_of = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh ? nh : 1)));
	unsigned int* fcomp = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nfaces + 1)));
	unsigned int* comp_outer = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nfaces + 1)));
	unsigned int* comp_hole = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (nfaces + 1)));
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
	vpos = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (x->nv + 1)));
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
				unsigned int nc = qaws_internal_grow(cappc, (size_t)np + 1, 64);
				cl_piece* gp = nc ? (cl_piece*)qaws_internal_realloc(pc, qaws_internal_mul2(sizeof(cl_piece), nc)) : NULL;
				if (!gp) { s = QAWS_STATUS_ALLOCATION_FAILURE; goto done; }
				pc = gp; cappc = nc;
			}
			pc[np].src = e->src;
			pc[np].ta = (h & 1) ? e->t1 : e->t0;
			pc[np].tb = (h & 1) ? e->t0 : e->t1;
			pc[np].va = cl_origin(x, h);
			pc[np].vb = cl_target(x, h);
			pc[np].h = h;
			if (x->h[h].area != x->h[h].area)
				x->h[h].area = cl_half_area(x, h);
			pc[np].area = x->h[h].area;
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
			cl_piece* g = (cl_piece*)qaws_internal_realloc(st, qaws_internal_mul2(sizeof(cl_piece), np));
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
					void* g1 = qaws_internal_realloc(loop_left, qaws_internal_mul2(sizeof(unsigned int), (nloop + 1)));
					void* g2 = g1 ? qaws_internal_realloc(loop_right, qaws_internal_mul2(sizeof(unsigned int), (nloop + 1))) : NULL;
					void* g3 = g2 ? qaws_internal_realloc(loop_area, qaws_internal_mul2(sizeof(double), (nloop + 1))) : NULL;
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
	qaws_internal_free(loop_of); qaws_internal_free(fcomp); qaws_internal_free(comp_outer); qaws_internal_free(comp_hole);
	qaws_internal_free(pc); qaws_internal_free(loop_left); qaws_internal_free(loop_right); qaws_internal_free(loop_area); qaws_internal_free(st); qaws_internal_free(vpos);
	return s;
}

/* ======================================================================== */
/*  Public                                                                  */
/* ======================================================================== */

void cl_free_ctx(cl_ctx* x)
{
	unsigned int i;
	for (i = 0; i < x->ntemp; i++) qaws_curve_destroy(x->temp[i]);
	qaws_internal_free(x->temp); qaws_internal_free(x->src); qaws_internal_free(x->pts); qaws_internal_free(x->v); qaws_internal_free(x->e); qaws_internal_free(x->oe);
	qaws_internal_free(x->h); qaws_internal_free(x->vout);
}

void qaws_clip_result_destroy(qaws_clip_result* r)
{
	unsigned int i;
	if (!r)
		return;
	for (i = 0; i < r->curve_count; i++)
		qaws_curve_destroy(r->curves[i]);
	qaws_internal_free(r->curves);
	qaws_internal_free((void*)r->views);
	qaws_internal_free(r->paths);
	qaws_internal_free(r->open_paths);
	qaws_internal_free(r->vertices);
	qaws_internal_free(r);
}

unsigned int qaws_clip_result_get_path_count(qaws_clip_result const* r)
{
	return r ? r->path_count : 0;
}

qaws_status qaws_clip_result_get_path(qaws_path_2d* out, qaws_clip_result const* r, unsigned int i)
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

qaws_status qaws_clip_result_get_open_path(qaws_path_2d* out, qaws_clip_result const* r, unsigned int i)
{
	if (!r || !out || i >= r->open_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	out->curves = r->views + r->open_paths[i].first;
	out->curve_count = r->open_paths[i].count;
	out->closed = 0;
	return QAWS_STATUS_OK;
}

qaws_status qaws_clip_result_get_vertices(qaws_clip_vertex const** out_vertices, unsigned int* out_count, qaws_clip_result const* r,
	int open, unsigned int i)
{
	struct cl_out_path const* p;
	if (!r || !out_vertices || !out_count || i >= (open ? r->open_count : r->path_count))
		return QAWS_STATUS_INVALID_ARGUMENT;
	p = open ? &r->open_paths[i] : &r->paths[i];
	*out_vertices = r->vertices + p->vfirst;
	*out_count = p->vcount;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Internal: results put together                                          */
/* ======================================================================== */

qaws_status qaws_internal_clip_result_empty(qaws_clip_result** out)
{
	*out = (qaws_clip_result*)qaws_internal_zalloc(1, sizeof(qaws_clip_result));
	return *out ? QAWS_STATUS_OK : QAWS_STATUS_ALLOCATION_FAILURE;
}

qaws_status qaws_internal_clip_result_append(qaws_clip_result* dst, qaws_clip_result* src)
{
	unsigned int coff = dst->curve_count, voff = dst->vertex_count, poff = dst->path_count, i;
	qaws_curve const** views;
	/* room for every curve first: a failure leaves each curve owned by src alone */
	if (dst->curve_cap - dst->curve_count < src->curve_count)
	{
		unsigned int nc = qaws_internal_grow(0, qaws_internal_add2(dst->curve_count, src->curve_count), 1);
		qaws_curve** np = nc ? (qaws_curve**)qaws_internal_realloc(dst->curves, qaws_internal_mul2(sizeof(qaws_curve*), nc)) : NULL;
		if (!np)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		dst->curves = np;
		dst->curve_cap = nc;
	}
	for (i = 0; i < src->curve_count; i++)
		dst->curves[dst->curve_count++] = src->curves[i];
	src->curve_count = 0;   /* moved */
	for (i = 0; i < src->vertex_count; i++)
		if (cl_res_vertex(dst, &src->vertices[i]) != QAWS_STATUS_OK)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < src->path_count; i++)
	{
		struct cl_out_path p = src->paths[i];
		p.first += coff;
		p.vfirst += voff;
		if (p.parent != CL_NONE) p.parent += poff;
		if (cl_push_path(dst, 0, &p) != QAWS_STATUS_OK)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < src->open_count; i++)
	{
		struct cl_out_path p = src->open_paths[i];
		p.first += coff;
		p.vfirst += voff;
		if (cl_push_path(dst, 1, &p) != QAWS_STATUS_OK)
			return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	qaws_clip_result_destroy(src);
	views = (qaws_curve const**)qaws_internal_realloc((void*)dst->views, qaws_internal_mul2(sizeof(qaws_curve*), (dst->curve_count ? dst->curve_count : 1)));
	if (!views)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	dst->views = views;
	for (i = 0; i < dst->curve_count; i++)
		dst->views[i] = dst->curves[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_internal_clip_output_arrangement(qaws_internal_arrangement const* a, qaws_clip_result* result)
{
	cl_ctx x;
	unsigned int i, k, nh = 2 * a->edge_count;
	qaws_status s;
	memset(&x, 0, sizeof(x));
	x.d = a->desc;
	x.r = result;
	x.tol = a->tolerance;
	x.ext = a->extent;
	x.src = (cl_src*)qaws_internal_malloc(qaws_internal_mul2(sizeof(cl_src), (a->source_count + 1)));
	x.v = (cl_vertex*)qaws_internal_malloc(qaws_internal_mul2(sizeof(cl_vertex), (a->vertex_count + 1)));
	x.e = (cl_edge*)qaws_internal_malloc(qaws_internal_mul2(sizeof(cl_edge), (a->edge_count + 1)));
	x.h = (cl_half*)qaws_internal_zalloc(nh + 1, sizeof(cl_half));
	x.vout = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (nh + 1)));
	if (!x.src || !x.v || !x.e || !x.h || !x.vout || a->edge_count >= 0x7FFFFFFFu)   /* nh + 1 counts in unsigned int */
	{
		cl_free_ctx(&x);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < a->source_count; i++)
	{
		x.src[i].c = a->sources[i].curve;
		x.src[i].operand = a->sources[i].operand;
		x.src[i].path = a->sources[i].path;
		x.src[i].curve = a->sources[i].index;
		x.src[i].base = a->sources[i].base;
		x.src[i].scale = a->sources[i].scale;
	}
	x.nsrc = a->source_count;
	for (i = 0; i < a->vertex_count; i++)
	{
		cl_vertex* v = &x.v[i];
		v->p[0] = a->vertices[i].p[0];
		v->p[1] = a->vertices[i].p[1];
		v->sa = a->vertices[i].sa;
		v->sb = a->vertices[i].sb;
		v->ta = a->vertices[i].ta;
		v->tb = a->vertices[i].tb;
		v->first = a->vertices[i].first;
		v->count = a->vertices[i].count;
	}
	x.nv = a->vertex_count;
	for (i = 0; i < a->edge_count; i++)
	{
		cl_edge* e = &x.e[i];
		e->src = a->edges[i].src;
		e->t0 = a->edges[i].t0;
		e->t1 = a->edges[i].t1;
		e->v0 = a->edges[i].v0;
		e->v1 = a->edges[i].v1;
		e->ds = a->edges[i].ds;
		e->dc = a->edges[i].dc;
		e->alive = 1;
	}
	x.ne = a->edge_count;
	for (i = 0; i < nh; i++)
	{
		x.vout[i] = a->out[i];
		x.h[i].face = a->half_face[i];
		x.h[i].area = NAN;
	}
	for (i = 0; i < a->vertex_count; i++)
		for (k = 0; k < a->vertices[i].count; k++)
			x.h[x.vout[a->vertices[i].first + k]].pos = k;
	s = cl_output(&x, a->face_count, a->inside);
	cl_free_ctx(&x);
	/* const views of the curves for the path views */
	if (s == QAWS_STATUS_OK)
	{
		qaws_curve const** views = (qaws_curve const**)qaws_internal_realloc((void*)result->views, qaws_internal_mul2(sizeof(qaws_curve*), (result->curve_count + 1)));
		if (!views)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		result->views = views;
		for (i = 0; i < result->curve_count; i++)
			result->views[i] = result->curves[i];
	}
	return s;
}
