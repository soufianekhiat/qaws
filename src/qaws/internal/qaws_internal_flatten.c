#include "qaws_internal_flatten.h"
#include "qaws_internal_types.h"
#include "qaws_internal_surface.h"
#include "../qaws_eval.h"
#include "../qaws_surface.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define FL_CURVE_MAX_DEPTH 22
#define FL_SURF_MIN_DEPTH  2
#define FL_SURF_MAX_DEPTH  10

static qaws_scalar fl_dot(qaws_scalar const* a, qaws_scalar const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

qaws_status qaws_internal_curve_point(qaws_curve const* c, unsigned int dim, qaws_scalar t, unsigned int flags, qaws_scalar* p, qaws_scalar* d)
{
	qaws_status s;
	if (dim == 2)
	{
		qaws_eval_result_2d r;
		s = qaws_curve_evaluate_2d(c, t, flags, &r);
		p[0] = r.position.x; p[1] = r.position.y; p[2] = 0;
		if (d) { d[0] = r.d1.x; d[1] = r.d1.y; d[2] = 0; }
	}
	else
	{
		qaws_eval_result_3d r;
		s = qaws_curve_evaluate_3d(c, t, flags, &r);
		p[0] = r.position.x; p[1] = r.position.y; p[2] = r.position.z;
		if (d) { d[0] = r.d1.x; d[1] = r.d1.y; d[2] = r.d1.z; }
	}
	return s;
}

static qaws_status fl_surface_point(qaws_surface const* s, qaws_scalar u, qaws_scalar v, qaws_scalar* p)
{
	qaws_surface_eval_result r;
	qaws_status st = qaws_surface_evaluate(s, u, v, QAWS_SURFACE_EVAL_POSITION, &r);
	p[0] = r.position.x; p[1] = r.position.y; p[2] = r.position.z;
	return st;
}

/* ------------------------------------------------------------------ */
/*  Curves                                                             */
/* ------------------------------------------------------------------ */

typedef struct fl_curve
{
	qaws_curve const* c;
	unsigned int dim, owner, index;
	qaws_scalar flat;
	qaws_flat_seg** segs;
	unsigned int* count;
	unsigned int* cap;
} fl_curve;

/* distance from m to the segment [a, b] */
static qaws_scalar fl_point_segment(qaws_scalar const* m, qaws_scalar const* a, qaws_scalar const* b)
{
	qaws_scalar ab[3], am[3], l2, u, d[3];
	int k;
	for (k = 0; k < 3; k++) { ab[k] = b[k] - a[k]; am[k] = m[k] - a[k]; }
	l2 = fl_dot(ab, ab);
	u = l2 > 0 ? fl_dot(am, ab) / l2 : 0;
	if (u < 0) u = 0;
	if (u > 1) u = 1;
	for (k = 0; k < 3; k++) d[k] = am[k] - u * ab[k];
	return (qaws_scalar)sqrt(fl_dot(d, d));
}

static qaws_status fl_push_seg(fl_curve* f, qaws_scalar t0, qaws_scalar const* p0, qaws_scalar t1, qaws_scalar const* p1, qaws_scalar dev)
{
	qaws_flat_seg* s;
	int k;
	if (*f->count == *f->cap)
	{
		unsigned int cap = *f->cap ? *f->cap * 2 : 1024;
		qaws_flat_seg* g = (qaws_flat_seg*)realloc(*f->segs, cap * sizeof(qaws_flat_seg));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		*f->segs = g;
		*f->cap = cap;
	}
	s = &(*f->segs)[(*f->count)++];
	s->t0 = t0; s->t1 = t1;
	s->r = 2 * dev;
	for (k = 0; k < 3; k++)
	{
		s->p0[k] = p0[k];
		s->p1[k] = p1[k];
		s->lo[k] = (p0[k] < p1[k] ? p0[k] : p1[k]) - s->r;
		s->hi[k] = (p0[k] > p1[k] ? p0[k] : p1[k]) + s->r;
	}
	s->owner = f->owner;
	s->index = f->index++;
	return QAWS_STATUS_OK;
}

/* piece [t0, t1] with its known midpoint: the midpoint and both quarter
   points must lie within the flatness bound of the chord */
static qaws_status fl_curve_piece(fl_curve* f, qaws_scalar t0, qaws_scalar const* p0, qaws_scalar tm, qaws_scalar const* pm,
	qaws_scalar t1, qaws_scalar const* p1, unsigned int depth)
{
	qaws_scalar ta = (t0 + tm) / 2, tb = (tm + t1) / 2, pa[3], pb[3], dev, d;
	qaws_status s = qaws_internal_curve_point(f->c, f->dim, ta, QAWS_EVAL_FLAG_POSITION, pa, NULL);
	if (s == QAWS_STATUS_OK)
		s = qaws_internal_curve_point(f->c, f->dim, tb, QAWS_EVAL_FLAG_POSITION, pb, NULL);
	if (s != QAWS_STATUS_OK)
		return s;
	dev = fl_point_segment(pm, p0, p1);
	d = fl_point_segment(pa, p0, p1);
	if (d > dev) dev = d;
	d = fl_point_segment(pb, p0, p1);
	if (d > dev) dev = d;
	if (depth >= FL_CURVE_MAX_DEPTH || dev <= f->flat)
		return fl_push_seg(f, t0, p0, t1, p1, dev);
	s = fl_curve_piece(f, t0, p0, ta, pa, tm, pm, depth + 1);
	if (s != QAWS_STATUS_OK)
		return s;
	return fl_curve_piece(f, tm, pm, tb, pb, t1, p1, depth + 1);
}

qaws_status qaws_internal_flatten_curve(qaws_curve const* c, unsigned int dim, qaws_scalar flatness, unsigned int owner,
	qaws_flat_seg** segs, unsigned int* count, unsigned int* capacity)
{
	fl_curve f;
	qaws_scalar p0[3], p1[3], pm[3];
	unsigned int k, have = 0;
	f.c = c; f.dim = dim; f.owner = owner; f.index = 0; f.flat = flatness;
	f.segs = segs; f.count = count; f.cap = capacity;
	for (k = 0; k < c->span_count; k++)
	{
		qaws_scalar t0 = c->span_boundaries[k], t1 = c->span_boundaries[k + 1], tm = (t0 + t1) / 2;
		qaws_status s = QAWS_STATUS_OK;
		if (!(t1 > t0))
			continue;
		/* a span starts where the previous one ended */
		if (!have)
			s = qaws_internal_curve_point(c, dim, t0, QAWS_EVAL_FLAG_POSITION, p0, NULL);
		if (s == QAWS_STATUS_OK)
			s = qaws_internal_curve_point(c, dim, t1, QAWS_EVAL_FLAG_POSITION, p1, NULL);
		if (s == QAWS_STATUS_OK)
			s = qaws_internal_curve_point(c, dim, tm, QAWS_EVAL_FLAG_POSITION, pm, NULL);
		if (s == QAWS_STATUS_OK)
			s = fl_curve_piece(&f, t0, p0, tm, pm, t1, p1, 0);
		if (s != QAWS_STATUS_OK)
			return s;
		memcpy(p0, p1, sizeof(p0));
		have = 1;
	}
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Surfaces                                                           */
/* ------------------------------------------------------------------ */

typedef struct fl_surf
{
	qaws_surface const* s;
	unsigned int owner;
	qaws_scalar flat;
	qaws_flat_patch** patches;
	unsigned int* count;
	unsigned int* cap;
} fl_surf;

static qaws_scalar fl_dist(qaws_scalar const* a, qaws_scalar const* b)
{
	qaws_scalar d[3];
	d[0] = a[0] - b[0]; d[1] = a[1] - b[1]; d[2] = a[2] - b[2];
	return (qaws_scalar)sqrt(fl_dot(d, d));
}

/* patch with corners c[4] ((u0,v0) (u1,v0) (u0,v1) (u1,v1)) */
static qaws_status fl_surf_patch(fl_surf* f, qaws_scalar u0, qaws_scalar u1, qaws_scalar v0, qaws_scalar v1, qaws_scalar (*c)[3], unsigned int depth)
{
	qaws_scalar um = (u0 + u1) / 2, vm = (v0 + v1) / 2, m[5][3], b[3], dev = 0;
	qaws_status s;
	int k, i;
	/* centre, then edge midpoints: v0 edge, v1 edge, u0 edge, u1 edge */
	s = fl_surface_point(f->s, um, vm, m[0]);
	if (s == QAWS_STATUS_OK) s = fl_surface_point(f->s, um, v0, m[1]);
	if (s == QAWS_STATUS_OK) s = fl_surface_point(f->s, um, v1, m[2]);
	if (s == QAWS_STATUS_OK) s = fl_surface_point(f->s, u0, vm, m[3]);
	if (s == QAWS_STATUS_OK) s = fl_surface_point(f->s, u1, vm, m[4]);
	if (s != QAWS_STATUS_OK)
		return s;
	for (k = 0; k < 3; k++) b[k] = (c[0][k] + c[1][k] + c[2][k] + c[3][k]) / 4;
	dev = fl_dist(m[0], b);
	for (k = 0; k < 3; k++) b[k] = (c[0][k] + c[1][k]) / 2;
	if (fl_dist(m[1], b) > dev) dev = fl_dist(m[1], b);
	for (k = 0; k < 3; k++) b[k] = (c[2][k] + c[3][k]) / 2;
	if (fl_dist(m[2], b) > dev) dev = fl_dist(m[2], b);
	for (k = 0; k < 3; k++) b[k] = (c[0][k] + c[2][k]) / 2;
	if (fl_dist(m[3], b) > dev) dev = fl_dist(m[3], b);
	for (k = 0; k < 3; k++) b[k] = (c[1][k] + c[3][k]) / 2;
	if (fl_dist(m[4], b) > dev) dev = fl_dist(m[4], b);
	if (depth >= FL_SURF_MAX_DEPTH || (depth >= FL_SURF_MIN_DEPTH && dev <= f->flat))
	{
		qaws_flat_patch* p;
		if (*f->count == *f->cap)
		{
			unsigned int cap = *f->cap ? *f->cap * 2 : 1024;
			qaws_flat_patch* g = (qaws_flat_patch*)realloc(*f->patches, cap * sizeof(qaws_flat_patch));
			if (!g)
				return QAWS_STATUS_ALLOCATION_FAILURE;
			*f->patches = g;
			*f->cap = cap;
		}
		p = &(*f->patches)[(*f->count)++];
		p->u0 = u0; p->u1 = u1; p->v0 = v0; p->v1 = v1;
		p->r = 2 * dev;
		p->owner = f->owner;
		for (k = 0; k < 3; k++)
		{
			qaws_scalar lo = c[0][k], hi = c[0][k];
			for (i = 0; i < 4; i++)
			{
				p->p[i][k] = c[i][k];
				if (c[i][k] < lo) lo = c[i][k];
				if (c[i][k] > hi) hi = c[i][k];
			}
			for (i = 0; i < 5; i++)
			{
				if (m[i][k] < lo) lo = m[i][k];
				if (m[i][k] > hi) hi = m[i][k];
			}
			p->lo[k] = lo - p->r;
			p->hi[k] = hi + p->r;
		}
		return QAWS_STATUS_OK;
	}
	{
		/* children: (u0..um, v0..vm) (um..u1, v0..vm) (u0..um, vm..v1) (um..u1, vm..v1) */
		qaws_scalar q[4][3];
		for (k = 0; k < 3; k++) { q[0][k] = c[0][k]; q[1][k] = m[1][k]; q[2][k] = m[3][k]; q[3][k] = m[0][k]; }
		s = fl_surf_patch(f, u0, um, v0, vm, q, depth + 1);
		if (s != QAWS_STATUS_OK) return s;
		for (k = 0; k < 3; k++) { q[0][k] = m[1][k]; q[1][k] = c[1][k]; q[2][k] = m[0][k]; q[3][k] = m[4][k]; }
		s = fl_surf_patch(f, um, u1, v0, vm, q, depth + 1);
		if (s != QAWS_STATUS_OK) return s;
		for (k = 0; k < 3; k++) { q[0][k] = m[3][k]; q[1][k] = m[0][k]; q[2][k] = c[2][k]; q[3][k] = m[2][k]; }
		s = fl_surf_patch(f, u0, um, vm, v1, q, depth + 1);
		if (s != QAWS_STATUS_OK) return s;
		for (k = 0; k < 3; k++) { q[0][k] = m[0][k]; q[1][k] = m[4][k]; q[2][k] = m[2][k]; q[3][k] = c[3][k]; }
		return fl_surf_patch(f, um, u1, vm, v1, q, depth + 1);
	}
}

qaws_status qaws_internal_flatten_surface(qaws_surface const* s, qaws_scalar flatness, unsigned int owner,
	qaws_flat_patch** patches, unsigned int* count, unsigned int* capacity)
{
	fl_surf f;
	qaws_scalar c[4][3];
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	qaws_status st;
	f.s = s; f.owner = owner; f.flat = flatness;
	f.patches = patches; f.count = count; f.cap = capacity;
	st = fl_surface_point(s, ur.min_value, vr.min_value, c[0]);
	if (st == QAWS_STATUS_OK) st = fl_surface_point(s, ur.max_value, vr.min_value, c[1]);
	if (st == QAWS_STATUS_OK) st = fl_surface_point(s, ur.min_value, vr.max_value, c[2]);
	if (st == QAWS_STATUS_OK) st = fl_surface_point(s, ur.max_value, vr.max_value, c[3]);
	if (st != QAWS_STATUS_OK)
		return st;
	return fl_surf_patch(&f, ur.min_value, ur.max_value, vr.min_value, vr.max_value, c, 0);
}

qaws_scalar qaws_internal_flatten_extent(qaws_curve const* const* curves, unsigned int curve_count, unsigned int dim,
	qaws_surface const* const* surfaces, unsigned int surface_count)
{
	qaws_scalar lo[3], hi[3], ext = 0, p[3];
	unsigned int i, j, k, l;
	for (k = 0; k < 3; k++) { lo[k] = (qaws_scalar)HUGE_VAL; hi[k] = -(qaws_scalar)HUGE_VAL; }
#define FL_GROW(p) for (k = 0; k < 3; k++) { if ((p)[k] < lo[k]) lo[k] = (p)[k]; if ((p)[k] > hi[k]) hi[k] = (p)[k]; }
	for (i = 0; i < curve_count; i++)
	{
		qaws_range r = curves[i]->parameter_range;
		for (j = 0; j <= 16; j++)
			if (qaws_internal_curve_point(curves[i], dim, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)j / 16, QAWS_EVAL_FLAG_POSITION, p, NULL) == QAWS_STATUS_OK)
				FL_GROW(p)
	}
	for (i = 0; i < surface_count; i++)
	{
		qaws_range ur = qaws_surface_get_u_range(surfaces[i]), vr = qaws_surface_get_v_range(surfaces[i]);
		for (j = 0; j <= 8; j++)
			for (l = 0; l <= 8; l++)
				if (fl_surface_point(surfaces[i], ur.min_value + (ur.max_value - ur.min_value) * (qaws_scalar)j / 8,
					vr.min_value + (vr.max_value - vr.min_value) * (qaws_scalar)l / 8, p) == QAWS_STATUS_OK)
					FL_GROW(p)
	}
#undef FL_GROW
	for (k = 0; k < 3; k++)
		if (hi[k] - lo[k] > ext)
			ext = hi[k] - lo[k];
	return ext > 0 ? ext : 1;
}
