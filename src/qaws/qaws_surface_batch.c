#include "qaws_surface_batch.h"
#include "qaws_eval.h"
#include "qaws_surface.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_broadphase.h"
#include "internal/qaws_internal_flatten.h"
#include "internal/qaws_internal_batch.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
 * Batched curve / surface intersection.
 *
 *   1. flatten every curve into chord segments and every surface into
 *      (u, v) patches (qaws_internal_flatten_*), boxes inflated by twice the
 *      measured deviation;
 *   2. one grid over all pieces (qaws_internal_broadphase); only segment /
 *      patch pairs are kept;
 *   3. the segment against the patch's two triangles, with a slack scaled by
 *      both inflations, seeds (t, u, v);
 *   4. Newton on S(u, v) - C(t) = 0, polished to rounding level;
 *   5. coincident hits of one curve / surface pair are kept once.
 */

#define SB_NEWTON_ITERS 24

#if QAWS_SCALAR_IS_FLOAT
#define SB_POS_REL ((qaws_scalar)2e-5)
#define SB_PAR_REL ((qaws_scalar)1e-4)
#else
#define SB_POS_REL ((qaws_scalar)1e-10)
#define SB_PAR_REL ((qaws_scalar)1e-7)
#endif

typedef struct sb_hit
{
	unsigned int curve, surface;
	qaws_scalar t, u, v;
	qaws_scalar p[3];
} sb_hit;

typedef struct sb_ctx
{
	qaws_curve_surface_batch_desc const* desc;
	qaws_flat_seg* segs;
	unsigned int nseg, capseg;
	qaws_flat_patch* patches;
	unsigned int npatch, cappatch;
	sb_hit* hits;
	unsigned int nhit, caphit;
	qaws_scalar flat, pos_tol;
	qaws_surface_batch_stats stats;
} sb_ctx;

static qaws_scalar sb_dot(qaws_scalar const* a, qaws_scalar const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void sb_cross(qaws_scalar const* a, qaws_scalar const* b, qaws_scalar* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

static qaws_scalar sb_det3(qaws_scalar const* a, qaws_scalar const* b, qaws_scalar const* c)
{
	qaws_scalar x[3];
	sb_cross(b, c, x);
	return sb_dot(a, x);
}

/* Newton on F(t, u, v) = S(u, v) - C(t), J = [-C', Su, Sv] */
static int sb_newton(sb_ctx* x, qaws_curve const* c, qaws_surface const* s, qaws_scalar* t, qaws_scalar* u, qaws_scalar* v, qaws_scalar* pos)
{
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	qaws_scalar tmin = c->parameter_range.min_value, tmax = c->parameter_range.max_value;
	unsigned int it, polish = 0;
	for (it = 0; it < SB_NEWTON_ITERS; it++)
	{
		qaws_scalar pc[3], dc[3], f[3], a[3], su[3], sv[3], det;
		qaws_surface_eval_result r;
		int k;
		if (qaws_internal_curve_point(c, 3, *t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1, pc, dc) != QAWS_STATUS_OK
			|| qaws_surface_evaluate(s, *u, *v, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r) != QAWS_STATUS_OK)
			return 0;
		f[0] = r.position.x - pc[0]; f[1] = r.position.y - pc[1]; f[2] = r.position.z - pc[2];
		su[0] = r.du.x; su[1] = r.du.y; su[2] = r.du.z;
		sv[0] = r.dv.x; sv[1] = r.dv.y; sv[2] = r.dv.z;
		for (k = 0; k < 3; k++)
			a[k] = -dc[k];
		/* converged: two more steps take it to rounding level */
		if (sqrt(sb_dot(f, f)) <= x->pos_tol && polish++ == 2)
		{
			pos[0] = (r.position.x + pc[0]) / 2; pos[1] = (r.position.y + pc[1]) / 2; pos[2] = (r.position.z + pc[2]) / 2;
			return 1;
		}
		det = sb_det3(a, su, sv);
		if (!(fabs(det) > (qaws_scalar)1e-12 * sqrt(sb_dot(a, a) * sb_dot(su, su) * sb_dot(sv, sv))))
		{
			if (!polish)
				return 0;
			pos[0] = (r.position.x + pc[0]) / 2; pos[1] = (r.position.y + pc[1]) / 2; pos[2] = (r.position.z + pc[2]) / 2;
			return 1;
		}
		/* Cramer: J d = f */
		*t -= sb_det3(f, su, sv) / det;
		*u -= sb_det3(a, f, sv) / det;
		*v -= sb_det3(a, su, f) / det;
		if (*t < tmin) *t = tmin;
		if (*t > tmax) *t = tmax;
		if (*u < ur.min_value) *u = ur.min_value;
		if (*u > ur.max_value) *u = ur.max_value;
		if (*v < vr.min_value) *v = vr.min_value;
		if (*v > vr.max_value) *v = vr.max_value;
	}
	return 0;
}

static qaws_status sb_push(sb_ctx* x, unsigned int curve, unsigned int surface, qaws_scalar t, qaws_scalar u, qaws_scalar v, qaws_scalar const* p)
{
	sb_hit* h;
	if (x->nhit == x->caphit)
	{
		unsigned int cap = x->caphit ? x->caphit * 2 : 256;
		sb_hit* g = (sb_hit*)realloc(x->hits, cap * sizeof(sb_hit));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		x->hits = g;
		x->caphit = cap;
	}
	h = &x->hits[x->nhit++];
	h->curve = curve; h->surface = surface;
	h->t = t; h->u = u; h->v = v;
	h->p[0] = p[0]; h->p[1] = p[1]; h->p[2] = p[2];
	return QAWS_STATUS_OK;
}

/* segment [a, b] against triangle (q0, q1, q2): line parameter s and
   barycentrics (b1, b2), with slack; 0 when parallel or outside */
static int sb_seg_tri(qaws_scalar const* a, qaws_scalar const* b, qaws_scalar const* q0, qaws_scalar const* q1, qaws_scalar const* q2,
	qaws_scalar eps, qaws_scalar delta, qaws_scalar* s, qaws_scalar* b1, qaws_scalar* b2)
{
	qaws_scalar d[3], e1[3], e2[3], h[3], w[3], q[3], det;
	int k;
	for (k = 0; k < 3; k++)
	{
		d[k] = b[k] - a[k];
		e1[k] = q1[k] - q0[k];
		e2[k] = q2[k] - q0[k];
		w[k] = a[k] - q0[k];
	}
	sb_cross(d, e2, h);
	det = sb_dot(e1, h);
	if (!(fabs(det) > (qaws_scalar)1e-14 * sqrt(sb_dot(d, d) * sb_dot(e1, e1) * sb_dot(e2, e2))))
		return 0;
	*b1 = sb_dot(w, h) / det;
	sb_cross(w, e1, q);
	*b2 = sb_dot(d, q) / det;
	*s = sb_dot(e2, q) / det;
	return *b1 >= -eps && *b2 >= -eps && *b1 + *b2 <= 1 + eps && *s >= -delta && *s <= 1 + delta;
}

static qaws_status sb_pair(sb_ctx* x, qaws_flat_seg const* S, qaws_flat_patch const* P)
{
	qaws_curve const* c = x->desc->curves[S->owner];
	qaws_surface const* sf = x->desc->surfaces[P->owner];
	/* triangles (00, 10, 11) and (00, 11, 01) with their (u, v) corners */
	static int const tri[2][3] = { { 0, 1, 3 }, { 0, 3, 2 } };
	qaws_scalar uc[4], vc[4], slack = P->r + S->r + x->pos_tol, size, len, d[3];
	int k, i;
	uc[0] = P->u0; vc[0] = P->v0; uc[1] = P->u1; vc[1] = P->v0;
	uc[2] = P->u0; vc[2] = P->v1; uc[3] = P->u1; vc[3] = P->v1;
	for (i = 0; i < 3; i++)
		d[i] = P->p[3][i] - P->p[0][i];
	size = (qaws_scalar)sqrt(sb_dot(d, d));
	for (i = 0; i < 3; i++)
		d[i] = S->p1[i] - S->p0[i];
	len = (qaws_scalar)sqrt(sb_dot(d, d));
	for (k = 0; k < 2; k++)
	{
		qaws_scalar s, b1, b2, t, u, v, pos[3];
		qaws_scalar eps = size > 0 ? 2 * slack / size : 1, delta = len > 0 ? 2 * slack / len : 1;
		if (!sb_seg_tri(S->p0, S->p1, P->p[tri[k][0]], P->p[tri[k][1]], P->p[tri[k][2]], eps, delta, &s, &b1, &b2))
			continue;
		x->stats.newton_count++;
		if (s < 0) s = 0;
		if (s > 1) s = 1;
		t = S->t0 + s * (S->t1 - S->t0);
		u = (1 - b1 - b2) * uc[tri[k][0]] + b1 * uc[tri[k][1]] + b2 * uc[tri[k][2]];
		v = (1 - b1 - b2) * vc[tri[k][0]] + b1 * vc[tri[k][1]] + b2 * vc[tri[k][2]];
		if (sb_newton(x, c, sf, &t, &u, &v, pos))
		{
			qaws_status st = sb_push(x, S->owner, P->owner, t, u, v, pos);
			if (st != QAWS_STATUS_OK)
				return st;
		}
	}
	return QAWS_STATUS_OK;
}

/* pieces 0..nseg-1 are segments, then patches: only mixed pairs */
static int sb_accept(void* user, unsigned int i, unsigned int j)
{
	sb_ctx const* x = (sb_ctx const*)user;
	return i < x->nseg && j >= x->nseg;
}

static qaws_status sb_visit(void* user, unsigned int i, unsigned int j)
{
	sb_ctx* x = (sb_ctx*)user;
	return sb_pair(x, &x->segs[i], &x->patches[j - x->nseg]);
}

static int sb_cmp(void const* p, void const* q)
{
	sb_hit const* a = (sb_hit const*)p;
	sb_hit const* b = (sb_hit const*)q;
	if (a->curve != b->curve) return a->curve < b->curve ? -1 : 1;
	if (a->surface != b->surface) return a->surface < b->surface ? -1 : 1;
	if (a->t != b->t) return a->t < b->t ? -1 : 1;
	return 0;
}

static void sb_merge(sb_ctx* x)
{
	unsigned int i, n = 0, group = 0;
	qsort(x->hits, x->nhit, sizeof(sb_hit), sb_cmp);
	for (i = 0; i < x->nhit; i++)
	{
		sb_hit const* h = &x->hits[i];
		qaws_curve const* c = x->desc->curves[h->curve];
		qaws_scalar pt = (c->parameter_range.max_value - c->parameter_range.min_value) * SB_PAR_REL;
		unsigned int j;
		int dup = 0;
		if (n == 0 || x->hits[n - 1].curve != h->curve || x->hits[n - 1].surface != h->surface)
			group = n;
		for (j = group; j < n && !dup; j++)
		{
			sb_hit const* k = &x->hits[j];
			qaws_scalar d[3];
			d[0] = k->p[0] - h->p[0]; d[1] = k->p[1] - h->p[1]; d[2] = k->p[2] - h->p[2];
			dup = fabs(k->t - h->t) <= pt || sqrt(sb_dot(d, d)) <= 16 * x->pos_tol;
		}
		if (!dup)
			x->hits[n++] = *h;
	}
	x->nhit = n;
}

/* the grid over the prepared segments and patches, refinement, merge, output */
static qaws_status sb_solve(sb_ctx* x, qaws_curve_surface_batch_hit* out_hits, unsigned int hit_capacity, unsigned int* out_count)
{
	qaws_bp_box* boxes = NULL;
	qaws_bp_stats bs;
	unsigned int i, k, n;
	qaws_status s = QAWS_STATUS_OK;
	x->stats.segment_count = x->nseg;
	x->stats.patch_count = x->npatch;
	if (x->nseg && x->npatch)
	{
		boxes = (qaws_bp_box*)malloc((x->nseg + x->npatch) * sizeof(qaws_bp_box));
		if (!boxes)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < x->nseg + x->npatch; i++)
			for (k = 0; k < 3; k++)
			{
				boxes[i].lo[k] = (double)(i < x->nseg ? x->segs[i].lo[k] : x->patches[i - x->nseg].lo[k]);
				boxes[i].hi[k] = (double)(i < x->nseg ? x->segs[i].hi[k] : x->patches[i - x->nseg].hi[k]);
			}
		s = qaws_internal_broadphase(boxes, x->nseg + x->npatch, 3, sb_accept, sb_visit, x, &bs);
		x->stats.cell_count = bs.cell_count;
		x->stats.candidate_count = bs.candidate_count;
		free(boxes);
	}
	if (s != QAWS_STATUS_OK)
		return s;
	sb_merge(x);
	for (n = 0; n < x->nhit && n < hit_capacity; n++)
	{
		sb_hit const* h = &x->hits[n];
		qaws_curve_surface_batch_hit* o = &out_hits[n];
		o->curve = h->curve; o->surface = h->surface;
		o->t = h->t; o->u = h->u; o->v = h->v;
		o->position.x = h->p[0]; o->position.y = h->p[1]; o->position.z = h->p[2];
	}
	*out_count = x->nhit;
	x->stats.hit_count = x->nhit;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_surface_batch_find_intersections(
	qaws_curve_surface_batch_desc const* desc,
	qaws_curve_surface_batch_hit* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_surface_batch_stats* out_stats)
{
	sb_ctx x;
	qaws_scalar ext;
	unsigned int i;
	qaws_status s = QAWS_STATUS_OK;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (!desc || !out_count || (!out_hits && hit_capacity) || (desc->curve_count && !desc->curves) || (desc->surface_count && !desc->surfaces))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	for (i = 0; i < desc->curve_count; i++)
	{
		if (!desc->curves[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (desc->curves[i]->dimension != QAWS_DIMENSION_3D)
			return QAWS_STATUS_INVALID_DIMENSION;
	}
	for (i = 0; i < desc->surface_count; i++)
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&x, 0, sizeof(x));
	x.desc = desc;
	ext = qaws_internal_flatten_extent(desc->curves, desc->curve_count, 3, desc->surfaces, desc->surface_count);
	x.flat = desc->flatness > 0 ? desc->flatness : ext / 512;
	x.pos_tol = ext * SB_POS_REL;
	for (i = 0; i < desc->curve_count && s == QAWS_STATUS_OK; i++)
		s = qaws_internal_flatten_curve(desc->curves[i], 3, x.flat, i, &x.segs, &x.nseg, &x.capseg);
	for (i = 0; i < desc->surface_count && s == QAWS_STATUS_OK; i++)
		s = qaws_internal_flatten_surface(desc->surfaces[i], x.flat, i, &x.patches, &x.npatch, &x.cappatch);
	if (s == QAWS_STATUS_OK)
		s = sb_solve(&x, out_hits, hit_capacity, out_count);
	if (out_stats)
		*out_stats = x.stats;
	free(x.segs);
	free(x.patches);
	free(x.hits);
	return s;
}

/* ------------------------------------------------------------------ */
/*  Prepared surface sets                                              */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_set_create(qaws_surface_batch_desc const* desc, qaws_surface_set** out_set)
{
	qaws_surface_set* set;
	unsigned int i, n, cap = 0;
	qaws_status s = QAWS_STATUS_OK;
	if (!desc || !out_set || (desc->surface_count && !desc->surfaces))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_set = NULL;
	n = desc->surface_count;
	for (i = 0; i < n; i++)
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	set = (qaws_surface_set*)calloc(1, sizeof(qaws_surface_set));
	if (!set)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	set->desc = *desc;
	set->desc.surfaces = (qaws_surface const* const*)malloc((n ? n : 1) * sizeof(qaws_surface*));
	set->desc.families = desc->families ? (unsigned int const*)malloc((n ? n : 1) * sizeof(unsigned int)) : NULL;
	if (!set->desc.surfaces || (desc->families && !set->desc.families))
	{
		qaws_surface_set_destroy(set);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy((void*)set->desc.surfaces, desc->surfaces, n * sizeof(qaws_surface*));
	if (desc->families)
		memcpy((void*)set->desc.families, desc->families, n * sizeof(unsigned int));
	set->ext = qaws_internal_flatten_extent(NULL, 0, 3, desc->surfaces, n);
	set->flat = desc->flatness > 0 ? desc->flatness : set->ext / 512;
	for (i = 0; i < n && s == QAWS_STATUS_OK; i++)
		s = qaws_internal_flatten_surface(desc->surfaces[i], set->flat, i, &set->patches, &set->npatch, &cap);
	if (s != QAWS_STATUS_OK)
	{
		qaws_surface_set_destroy(set);
		return s;
	}
	*out_set = set;
	return QAWS_STATUS_OK;
}

void qaws_surface_set_destroy(qaws_surface_set* set)
{
	if (!set)
		return;
	free((void*)set->desc.surfaces);
	free((void*)set->desc.families);
	free(set->patches);
	free(set);
}

unsigned int qaws_surface_set_get_patch_count(qaws_surface_set const* set)
{
	return set ? set->npatch : 0;
}

qaws_status qaws_curve_set_find_surface_intersections(
	qaws_curve_set const* curves,
	qaws_surface_set const* surfaces,
	qaws_curve_surface_batch_hit* out_hits,
	unsigned int hit_capacity,
	unsigned int* out_count,
	qaws_surface_batch_stats* out_stats)
{
	sb_ctx x;
	qaws_curve_surface_batch_desc d;
	qaws_status s;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (!curves || !surfaces || !out_count || (!out_hits && hit_capacity))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (curves->desc.curve_count && curves->dim != 3)
		return QAWS_STATUS_INVALID_DIMENSION;
	memset(&d, 0, sizeof(d));
	d.curves = curves->desc.curves;
	d.curve_count = curves->desc.curve_count;
	d.surfaces = surfaces->desc.surfaces;
	d.surface_count = surfaces->desc.surface_count;
	memset(&x, 0, sizeof(x));
	x.desc = &d;
	/* the prepared pieces as they are */
	x.segs = curves->segs;
	x.nseg = curves->nseg;
	x.patches = surfaces->patches;
	x.npatch = surfaces->npatch;
	x.flat = curves->flat > surfaces->flat ? curves->flat : surfaces->flat;
	x.pos_tol = (curves->ext > surfaces->ext ? curves->ext : surfaces->ext) * SB_POS_REL;
	s = sb_solve(&x, out_hits, hit_capacity, out_count);
	if (out_stats)
		*out_stats = x.stats;
	free(x.hits);
	return s;
}

/* ------------------------------------------------------------------ */
/*  Surfaces x surfaces                                                */
/* ------------------------------------------------------------------ */

typedef struct ss_seg
{
	unsigned int a, b;          /* surfaces, a < b */
	qaws_ssi_point p[2];
} ss_seg;

typedef struct ss_ctx
{
	qaws_surface_batch_desc const* desc;
	qaws_flat_patch* patches;
	unsigned int npatch, cappatch;
	ss_seg* segs;
	unsigned int nseg, capseg;
	qaws_scalar flat, pos_tol;
	qaws_surface_batch_stats stats;
} ss_ctx;

static qaws_scalar ss_det3x3(qaws_scalar const m[3][3])
{
	return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
		+ m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

static void ss_finish(qaws_ssi_point* pt, qaws_surface_eval_result const* ra, qaws_surface_eval_result const* rb)
{
	pt->position.x = (ra->position.x + rb->position.x) / 2;
	pt->position.y = (ra->position.y + rb->position.y) / 2;
	pt->position.z = (ra->position.z + rb->position.z) / 2;
}

/* minimum-norm Newton on F = S_a(u1, v1) - S_b(u2, v2): dx = J^T (J J^T)^-1 F */
static int ss_newton(ss_ctx* x, qaws_surface const* sa, qaws_surface const* sb, qaws_ssi_point* pt)
{
	qaws_range ua = qaws_surface_get_u_range(sa), va = qaws_surface_get_v_range(sa);
	qaws_range ub = qaws_surface_get_u_range(sb), vb = qaws_surface_get_v_range(sb);
	unsigned int it, polish = 0, flags = QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV;
	for (it = 0; it < SB_NEWTON_ITERS; it++)
	{
		qaws_surface_eval_result ra, rb;
		qaws_scalar f[3], J[4][3], M[3][3], y[3], det;
		int r, c, k;
		if (qaws_surface_evaluate(sa, pt->u1, pt->v1, flags, &ra) != QAWS_STATUS_OK || qaws_surface_evaluate(sb, pt->u2, pt->v2, flags, &rb) != QAWS_STATUS_OK)
			return 0;
		f[0] = ra.position.x - rb.position.x; f[1] = ra.position.y - rb.position.y; f[2] = ra.position.z - rb.position.z;
		/* converged: two more steps take it to rounding level */
		if (sqrt(sb_dot(f, f)) <= x->pos_tol && polish++ == 2)
		{
			ss_finish(pt, &ra, &rb);
			return 1;
		}
		J[0][0] = ra.du.x; J[0][1] = ra.du.y; J[0][2] = ra.du.z;
		J[1][0] = ra.dv.x; J[1][1] = ra.dv.y; J[1][2] = ra.dv.z;
		J[2][0] = -rb.du.x; J[2][1] = -rb.du.y; J[2][2] = -rb.du.z;
		J[3][0] = -rb.dv.x; J[3][1] = -rb.dv.y; J[3][2] = -rb.dv.z;
		for (r = 0; r < 3; r++)
			for (c = 0; c < 3; c++)
			{
				M[r][c] = 0;
				for (k = 0; k < 4; k++)
					M[r][c] += J[k][r] * J[k][c];
			}
		det = ss_det3x3(M);
		if (!(fabs(det) > (qaws_scalar)1e-14 * (M[0][0] * M[1][1] * M[2][2])))
		{
			if (!polish)
				return 0;
			ss_finish(pt, &ra, &rb);
			return 1;
		}
		for (c = 0; c < 3; c++)
		{
			qaws_scalar Mc[3][3];
			for (r = 0; r < 3; r++)
				for (k = 0; k < 3; k++)
					Mc[r][k] = k == c ? f[r] : M[r][k];
			y[c] = ss_det3x3(Mc) / det;
		}
		pt->u1 -= sb_dot(J[0], y);
		pt->v1 -= sb_dot(J[1], y);
		pt->u2 -= sb_dot(J[2], y);
		pt->v2 -= sb_dot(J[3], y);
		if (pt->u1 < ua.min_value) pt->u1 = ua.min_value;
		if (pt->u1 > ua.max_value) pt->u1 = ua.max_value;
		if (pt->v1 < va.min_value) pt->v1 = va.min_value;
		if (pt->v1 > va.max_value) pt->v1 = va.max_value;
		if (pt->u2 < ub.min_value) pt->u2 = ub.min_value;
		if (pt->u2 > ub.max_value) pt->u2 = ub.max_value;
		if (pt->v2 < vb.min_value) pt->v2 = vb.min_value;
		if (pt->v2 > vb.max_value) pt->v2 = vb.max_value;
	}
	return 0;
}

/* barycentrics of p in triangle (q0, q1, q2), projected on its plane */
static void ss_bary(qaws_scalar const* p, qaws_scalar const* q0, qaws_scalar const* q1, qaws_scalar const* q2, qaws_scalar* w)
{
	qaws_scalar v0[3], v1[3], v2[3], d00, d01, d11, d20, d21, den;
	int k;
	for (k = 0; k < 3; k++)
	{
		v0[k] = q1[k] - q0[k];
		v1[k] = q2[k] - q0[k];
		v2[k] = p[k] - q0[k];
	}
	d00 = sb_dot(v0, v0); d01 = sb_dot(v0, v1); d11 = sb_dot(v1, v1);
	d20 = sb_dot(v2, v0); d21 = sb_dot(v2, v1);
	den = d00 * d11 - d01 * d01;
	if (!(den > 0))
	{
		w[0] = 1; w[1] = w[2] = 0;
		return;
	}
	w[1] = (d11 * d20 - d01 * d21) / den;
	w[2] = (d00 * d21 - d01 * d20) / den;
	w[0] = 1 - w[1] - w[2];
}

typedef struct ss_cand
{
	qaws_scalar p[3];
	qaws_scalar wa[3], wb[3];   /* barycentrics in the triangles of A and B */
} ss_cand;

/* points where the edges of triangle e cross the plane of triangle g (or its
   vertices lie exactly on it) inside g,
   appended to out[n..4); e_is_a: e is the triangle of surface A */
static unsigned int ss_edges(qaws_scalar const* const* e, qaws_scalar const* const* g, int e_is_a, ss_cand* out, unsigned int n)
{
	qaws_scalar g1[3], g2[3], nrm[3], d[3], eps = (qaws_scalar)1e-9;
	int i, k;
	for (k = 0; k < 3; k++)
	{
		g1[k] = g[1][k] - g[0][k];
		g2[k] = g[2][k] - g[0][k];
	}
	sb_cross(g1, g2, nrm);
	for (i = 0; i < 3; i++)
	{
		qaws_scalar w[3];
		for (k = 0; k < 3; k++)
			w[k] = e[i][k] - g[0][k];
		d[i] = sb_dot(nrm, w);
	}
	for (i = 0; i < 3 && n < 4; i++)
	{
		int j = (i + 1) % 3;
		qaws_scalar s, wo[3], we[3];
		ss_cand* c;
		/* a strict crossing of edge (i, j), or vertex i exactly on the plane */
		if (d[i] == 0)
			s = 0;
		else if ((d[i] < 0 && d[j] > 0) || (d[i] > 0 && d[j] < 0))
			s = d[i] / (d[i] - d[j]);
		else
			continue;
		c = &out[n];
		for (k = 0; k < 3; k++)
			c->p[k] = e[i][k] + s * (e[j][k] - e[i][k]);
		ss_bary(c->p, g[0], g[1], g[2], wo);
		if (wo[0] < -eps || wo[1] < -eps || wo[2] < -eps)
			continue;
		we[0] = we[1] = we[2] = 0;
		we[i] = 1 - s;
		we[j] = s;
		for (k = 0; k < 3; k++)
		{
			c->wa[k] = e_is_a ? we[k] : wo[k];
			c->wb[k] = e_is_a ? wo[k] : we[k];
		}
		n++;
	}
	return n;
}

static qaws_status ss_pair(ss_ctx* x, qaws_flat_patch const* A, qaws_flat_patch const* B)
{
	static int const tri[2][3] = { { 0, 1, 3 }, { 0, 3, 2 } };
	qaws_surface const* sa = x->desc->surfaces[A->owner];
	qaws_surface const* sb = x->desc->surfaces[B->owner];
	qaws_scalar ua[4], va[4], ub[4], vb[4];
	int i, j, k;
	ua[0] = A->u0; va[0] = A->v0; ua[1] = A->u1; va[1] = A->v0; ua[2] = A->u0; va[2] = A->v1; ua[3] = A->u1; va[3] = A->v1;
	ub[0] = B->u0; vb[0] = B->v0; ub[1] = B->u1; vb[1] = B->v0; ub[2] = B->u0; vb[2] = B->v1; ub[3] = B->u1; vb[3] = B->v1;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
		{
			qaws_scalar const* ta[3];
			qaws_scalar const* tb[3];
			ss_cand c[4];
			unsigned int n, p0 = 0, p1 = 0, q, r;
			qaws_scalar best = -1;
			ss_seg sg;
			for (k = 0; k < 3; k++)
			{
				ta[k] = A->p[tri[i][k]];
				tb[k] = B->p[tri[j][k]];
			}
			n = ss_edges(ta, tb, 1, c, 0);
			n = ss_edges(tb, ta, 0, c, n);
			if (n < 2)
				continue;
			/* the two candidates farthest apart end the segment */
			for (q = 0; q < n; q++)
				for (r = q + 1; r < n; r++)
				{
					qaws_scalar dd[3], l;
					for (k = 0; k < 3; k++)
						dd[k] = c[q].p[k] - c[r].p[k];
					l = sb_dot(dd, dd);
					if (l > best)
					{
						best = l;
						p0 = q;
						p1 = r;
					}
				}
			if (!(best > 0))
				continue;
			sg.a = A->owner;
			sg.b = B->owner;
			for (q = 0; q < 2; q++)
			{
				ss_cand const* e = &c[q ? p1 : p0];
				qaws_ssi_point* pt = &sg.p[q];
				pt->u1 = pt->v1 = pt->u2 = pt->v2 = 0;
				for (k = 0; k < 3; k++)
				{
					pt->u1 += e->wa[k] * ua[tri[i][k]];
					pt->v1 += e->wa[k] * va[tri[i][k]];
					pt->u2 += e->wb[k] * ub[tri[j][k]];
					pt->v2 += e->wb[k] * vb[tri[j][k]];
				}
				x->stats.newton_count++;
				if (!ss_newton(x, sa, sb, pt))
					break;
			}
			if (q < 2)
				continue;
			if (x->nseg == x->capseg)
			{
				unsigned int cap = x->capseg ? x->capseg * 2 : 256;
				ss_seg* g = (ss_seg*)realloc(x->segs, cap * sizeof(ss_seg));
				if (!g)
					return QAWS_STATUS_ALLOCATION_FAILURE;
				x->segs = g;
				x->capseg = cap;
			}
			x->segs[x->nseg++] = sg;
		}
	return QAWS_STATUS_OK;
}

static int ss_accept(void* user, unsigned int i, unsigned int j)
{
	ss_ctx const* x = (ss_ctx const*)user;
	unsigned int a = x->patches[i].owner, b = x->patches[j].owner;
	return a != b && (!x->desc->families || x->desc->families[a] != x->desc->families[b]);
}

static qaws_status ss_visit(void* user, unsigned int i, unsigned int j)
{
	ss_ctx* x = (ss_ctx*)user;
	return ss_pair(x, &x->patches[i], &x->patches[j]);
}

static int ss_cmp_seg(void const* p, void const* q)
{
	ss_seg const* a = (ss_seg const*)p;
	ss_seg const* b = (ss_seg const*)q;
	if (a->a != b->a) return a->a < b->a ? -1 : 1;
	if (a->b != b->b) return a->b < b->b ? -1 : 1;
	return 0;
}

typedef struct ss_end
{
	qaws_scalar x;
	unsigned int id;            /* 2 segment + end */
} ss_end;

static int ss_cmp_end(void const* p, void const* q)
{
	qaws_scalar a = ((ss_end const*)p)->x, b = ((ss_end const*)q)->x;
	return a < b ? -1 : (a > b ? 1 : 0);
}

static qaws_scalar ss_dist2(qaws_ssi_point const* p, qaws_ssi_point const* q)
{
	qaws_scalar dx = p->position.x - q->position.x, dy = p->position.y - q->position.y, dz = p->position.z - q->position.z;
	return dx * dx + dy * dy + dz * dz;
}

typedef struct ss_bridge
{
	unsigned int a, b;          /* dangling root nodes */
	qaws_scalar d2;
} ss_bridge;

static int ss_cmp_bridge(void const* p, void const* q)
{
	qaws_scalar a = ((ss_bridge const*)p)->d2, b = ((ss_bridge const*)q)->d2;
	return a < b ? -1 : (a > b ? 1 : 0);
}

static unsigned int ss_find(unsigned int* parent, unsigned int i)
{
	while (parent[i] != i)
	{
		parent[i] = parent[parent[i]];
		i = parent[i];
	}
	return i;
}

typedef struct ss_out
{
	qaws_surface_batch_curve* curves;
	unsigned int curve_cap, ncurve;
	qaws_ssi_point* points;
	unsigned int point_cap, npoint;
} ss_out;

/* one polyline; written only when the whole curve fits */
static void ss_emit(ss_out* o, unsigned int a, unsigned int b, qaws_ssi_point const* const* pts, unsigned int n, int closed)
{
	unsigned int k;
	if (o->ncurve < o->curve_cap && o->npoint + n <= o->point_cap)
	{
		qaws_surface_batch_curve* c = &o->curves[o->ncurve];
		c->surface_a = a;
		c->surface_b = b;
		c->first = o->npoint;
		c->count = n;
		c->closed = closed;
		for (k = 0; k < n; k++)
			o->points[o->npoint + k] = *pts[k];
	}
	o->ncurve++;
	o->npoint += n;
}

/* Nodes of the segment graph in parent (union-find over the 2 m end points):
   end points that coincide are joined, then each dangling end is bridged
   to its nearest dangling partner within the bridge distance, one to one,
   so runs of short segments never collapse. deg receives the final degree
   of every root. */
static qaws_status ss_nodes(ss_ctx const* x, ss_seg const* segs, unsigned int m, ss_end* ends, unsigned int* parent, unsigned int* deg)
{
	unsigned int ne = 2 * m, i, j;
	qaws_scalar tight = 64 * x->pos_tol, bridge = 4 * x->flat + 16 * x->pos_tol;
	qaws_status st = QAWS_STATUS_OK;
#define SS_PT(e) (&segs[(e) / 2].p[(e) % 2])
	memset(deg, 0, ne * sizeof(unsigned int));
	/* 1. end points that coincide (neighbouring triangle pairs compute the
	   same edge / plane point): a tight union, along a sweep in x */
	for (i = 0; i < ne; i++)
	{
		ends[i].x = SS_PT(i)->position.x;
		ends[i].id = i;
		parent[i] = i;
	}
	qsort(ends, ne, sizeof(ss_end), ss_cmp_end);
	for (i = 0; i < ne; i++)
		for (j = i + 1; j < ne && ends[j].x - ends[i].x <= tight; j++)
			if (ss_dist2(SS_PT(ends[i].id), SS_PT(ends[j].id)) <= tight * tight)
			{
				unsigned int ri = ss_find(parent, ends[i].id), rj = ss_find(parent, ends[j].id);
				if (ri != rj)
					parent[ri < rj ? rj : ri] = ri < rj ? ri : rj;
			}
	/* 2. gaps (quadtree T-junctions): each dangling end joins its nearest
	   dangling partner within the bridge distance, one to one, so runs of
	   short segments never collapse */
	for (i = 0; i < m; i++)
	{
		unsigned int u = ss_find(parent, 2 * i), v = ss_find(parent, 2 * i + 1);
		if (u != v)
		{
			deg[u]++;
			deg[v]++;
		}
	}
	{
		ss_bridge* br = NULL;
		unsigned int nbr = 0, capbr = 0;
		for (i = 0; i < ne; i++)
		{
			unsigned int ri = ss_find(parent, ends[i].id);
			if (deg[ri] != 1)
				continue;
			for (j = i + 1; j < ne && ends[j].x - ends[i].x <= bridge; j++)
			{
				unsigned int rj = ss_find(parent, ends[j].id);
				qaws_scalar d2;
				if (rj == ri || deg[rj] != 1)
					continue;
				d2 = ss_dist2(SS_PT(ends[i].id), SS_PT(ends[j].id));
				if (d2 > bridge * bridge)
					continue;
				if (nbr == capbr)
				{
					ss_bridge* g;
					capbr = capbr ? capbr * 2 : 64;
					g = (ss_bridge*)realloc(br, capbr * sizeof(ss_bridge));
					if (!g)
					{
						free(br);
						st = QAWS_STATUS_ALLOCATION_FAILURE;
						goto done;
					}
					br = g;
				}
				br[nbr].a = ri;
				br[nbr].b = rj;
				br[nbr].d2 = d2;
				nbr++;
			}
		}
		qsort(br, nbr, sizeof(ss_bridge), ss_cmp_bridge);
		for (i = 0; i < nbr; i++)
		{
			unsigned int a = br[i].a, b = br[i].b;
			if (deg[a] != 1 || deg[b] != 1 || ss_find(parent, a) != a || ss_find(parent, b) != b)
				continue;
			parent[a < b ? b : a] = a < b ? a : b;
			deg[a < b ? a : b] = 2;
			deg[a < b ? b : a] = 0;
		}
		free(br);
	}
	/* final degrees */
	memset(deg, 0, ne * sizeof(unsigned int));
	for (i = 0; i < m; i++)
	{
		unsigned int u = ss_find(parent, 2 * i), v = ss_find(parent, 2 * i + 1);
		if (u != v)
		{
			deg[u]++;
			deg[v]++;
		}
	}
done:
#undef SS_PT
	return st;
}

/* is the point on a parameter boundary of either surface (a true open end)? */
static int ss_on_boundary(qaws_surface const* sa, qaws_surface const* sb, qaws_ssi_point const* p)
{
	qaws_range r[4];
	qaws_scalar v[4];
	int k;
	r[0] = qaws_surface_get_u_range(sa); r[1] = qaws_surface_get_v_range(sa);
	r[2] = qaws_surface_get_u_range(sb); r[3] = qaws_surface_get_v_range(sb);
	v[0] = p->u1; v[1] = p->v1; v[2] = p->u2; v[3] = p->v2;
	for (k = 0; k < 4; k++)
	{
		qaws_scalar e = (r[k].max_value - r[k].min_value) * SB_PAR_REL;
		if (v[k] <= r[k].min_value + e || v[k] >= r[k].max_value - e)
			return 1;
	}
	return 0;
}

static qaws_status ss_append(ss_seg** L, unsigned int* m, unsigned int* cap, unsigned int a, unsigned int b, qaws_ssi_point const* p, qaws_ssi_point const* q)
{
	if (*m == *cap)
	{
		unsigned int c = *cap ? *cap * 2 : 64;
		ss_seg* g = (ss_seg*)realloc(*L, c * sizeof(ss_seg));
		if (!g)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		*L = g;
		*cap = c;
	}
	(*L)[*m].a = a;
	(*L)[*m].b = b;
	(*L)[*m].p[0] = *p;
	(*L)[*m].p[1] = *q;
	(*m)++;
	return QAWS_STATUS_OK;
}

/* one predictor step of length h along T = n_a x n_b (oriented by dir),
   corrected onto both surfaces */
static int ss_step(ss_ctx* x, qaws_surface const* sa, qaws_surface const* sb, qaws_ssi_point* p, qaws_scalar* dir, qaws_scalar h)
{
	unsigned int flags = QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV;
	qaws_surface_eval_result ra, rb;
	qaws_scalar au[3], av[3], bu[3], bv[3], na[3], nb[3], t[3], l, d[3];
	qaws_surface_eval_result const* rr[2];
	qaws_scalar const* uu[2];
	qaws_scalar const* vv[2];
	qaws_scalar delta[2][2];
	int k, s;
	if (qaws_surface_evaluate(sa, p->u1, p->v1, flags, &ra) != QAWS_STATUS_OK || qaws_surface_evaluate(sb, p->u2, p->v2, flags, &rb) != QAWS_STATUS_OK)
		return 0;
	au[0] = ra.du.x; au[1] = ra.du.y; au[2] = ra.du.z; av[0] = ra.dv.x; av[1] = ra.dv.y; av[2] = ra.dv.z;
	bu[0] = rb.du.x; bu[1] = rb.du.y; bu[2] = rb.du.z; bv[0] = rb.dv.x; bv[1] = rb.dv.y; bv[2] = rb.dv.z;
	sb_cross(au, av, na);
	sb_cross(bu, bv, nb);
	sb_cross(na, nb, t);
	l = (qaws_scalar)sqrt(sb_dot(t, t));
	if (!(l > 0))
		return 0;
	if (sb_dot(t, dir) < 0)
		l = -l;
	for (k = 0; k < 3; k++)
		d[k] = h * t[k] / l;
	/* parameter moves of each surface for the move d (least squares) */
	rr[0] = &ra; rr[1] = &rb;
	uu[0] = au; vv[0] = av; uu[1] = bu; vv[1] = bv;
	for (s = 0; s < 2; s++)
	{
		qaws_scalar m00 = sb_dot(uu[s], uu[s]), m01 = sb_dot(uu[s], vv[s]), m11 = sb_dot(vv[s], vv[s]);
		qaws_scalar g0 = sb_dot(uu[s], d), g1 = sb_dot(vv[s], d), det = m00 * m11 - m01 * m01;
		if (!(det > 0))
			return 0;
		delta[s][0] = (m11 * g0 - m01 * g1) / det;
		delta[s][1] = (m00 * g1 - m01 * g0) / det;
	}
	(void)rr;
	p->u1 += delta[0][0]; p->v1 += delta[0][1];
	p->u2 += delta[1][0]; p->v2 += delta[1][1];
	for (k = 0; k < 3; k++)
		dir[k] = d[k];
	return ss_newton(x, sa, sb, p);
}

#define SS_MARCH_STEPS 96

/* Gaps where the flattened triangles missed (a shallow crossing angle):
   from every dangling end off the surface boundaries, march along the
   intersection curve until another dangling end is within reach; the new
   segments are appended to *L. */
static qaws_status ss_close_gaps(ss_ctx* x, ss_seg** L, unsigned int* m, unsigned int* cap)
{
	unsigned int n0 = *m, ne = 2 * n0, i;
	qaws_scalar h = 2 * x->flat;
	qaws_surface const* sa = x->desc->surfaces[(*L)[0].a];
	qaws_surface const* sb = x->desc->surfaces[(*L)[0].b];
	unsigned int a = (*L)[0].a, b = (*L)[0].b;
	ss_end* ends = (ss_end*)malloc(ne * sizeof(ss_end));
	unsigned int* parent = (unsigned int*)malloc(ne * sizeof(unsigned int));
	unsigned int* deg = (unsigned int*)calloc(ne, sizeof(unsigned int));
	unsigned char* dangling = (unsigned char*)calloc(ne, 1);
	qaws_status st = QAWS_STATUS_OK;
#define SG_PT(e) (&(*L)[(e) / 2].p[(e) % 2])
	if (!ends || !parent || !deg || !dangling)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* the nodes after joining and bridging: only real gaps remain */
	st = ss_nodes(x, *L, n0, ends, parent, deg);
	if (st != QAWS_STATUS_OK)
		goto done;
	/* dangling: a node of degree 1 off the boundaries, marked on the end
	   point whose segment leaves the node (it gives the outward direction) */
	for (i = 0; i < ne; i++)
	{
		unsigned int r = ss_find(parent, i);
		if (deg[r] == 1 && ss_find(parent, i ^ 1u) != r && !ss_on_boundary(sa, sb, SG_PT(i)))
			dangling[i] = 1;
	}
	for (i = 0; i < ne && st == QAWS_STATUS_OK; i++)
	{
		qaws_ssi_point cur, nxt;
		qaws_scalar dir[3];
		unsigned int k, e = i;
		if (!dangling[i])
			continue;
		dangling[i] = 0;
		/* the end of the segment holding the root point i */
		cur = *SG_PT(e);
		{
			qaws_ssi_point const* q = &(*L)[e / 2].p[1 - e % 2];
			dir[0] = cur.position.x - q->position.x;
			dir[1] = cur.position.y - q->position.y;
			dir[2] = cur.position.z - q->position.z;
		}
		for (k = 0; k < SS_MARCH_STEPS && st == QAWS_STATUS_OK; k++)
		{
			unsigned int t, hit = ne;
			nxt = cur;
			if (!ss_step(x, sa, sb, &nxt, dir, h))
				break;
			/* another dangling end within reach closes the gap */
			for (t = 0; t < ne; t++)
				if (dangling[t] && ss_dist2(&nxt, SG_PT(t)) <= (qaws_scalar)2.25 * h * h)
				{
					hit = t;
					break;
				}
			if (hit < ne)
			{
				qaws_ssi_point target = *SG_PT(hit);
				st = ss_append(L, m, cap, a, b, &cur, &nxt);
				if (st == QAWS_STATUS_OK)
					st = ss_append(L, m, cap, a, b, &nxt, &target);
				dangling[hit] = 0;
				break;
			}
			st = ss_append(L, m, cap, a, b, &cur, &nxt);
			if (ss_on_boundary(sa, sb, &nxt))
				break;
			cur = nxt;
		}
	}
#undef SG_PT
done:
	free(ends);
	free(parent);
	free(deg);
	free(dangling);
	return st;
}

/* chains the m segments of one surface pair into polylines */
static qaws_status ss_chain(ss_ctx* x, ss_seg const* segs, unsigned int m, ss_out* o)
{
	unsigned int ne = 2 * m, i;
	ss_end* ends = (ss_end*)malloc(ne * sizeof(ss_end));
	unsigned int* parent = (unsigned int*)malloc(ne * sizeof(unsigned int));
	unsigned int* deg = (unsigned int*)calloc(ne, sizeof(unsigned int));
	unsigned int* adj_start = (unsigned int*)calloc(ne + 1, sizeof(unsigned int));
	unsigned int* adj = (unsigned int*)malloc(ne * sizeof(unsigned int));
	unsigned char* used = (unsigned char*)calloc(m ? m : 1, 1);
	qaws_ssi_point const** line = (qaws_ssi_point const**)malloc((m + 2) * sizeof(qaws_ssi_point*));
	qaws_status st = QAWS_STATUS_OK;
#define SS_PT(e) (&segs[(e) / 2].p[(e) % 2])
	if (!ends || !parent || !deg || !adj_start || !adj || !used || !line)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	st = ss_nodes(x, segs, m, ends, parent, deg);
	if (st != QAWS_STATUS_OK)
		goto done;
	memset(deg, 0, ne * sizeof(unsigned int));
	/* graph on the root nodes: segment s joins root(2s) and root(2s + 1) */
	for (i = 0; i < m; i++)
	{
		unsigned int u = ss_find(parent, 2 * i), v = ss_find(parent, 2 * i + 1);
		if (u == v)
		{
			used[i] = 1;   /* shorter than the merge distance */
			continue;
		}
		deg[u]++;
		deg[v]++;
	}
	for (i = 0; i < ne; i++)
		adj_start[i + 1] = adj_start[i] + deg[i];
	memset(deg, 0, ne * sizeof(unsigned int));
	for (i = 0; i < m; i++)
		if (!used[i])
		{
			unsigned int u = ss_find(parent, 2 * i), v = ss_find(parent, 2 * i + 1);
			adj[adj_start[u] + deg[u]++] = i;
			adj[adj_start[v] + deg[v]++] = i;
		}
	/* walks: from ends and junctions first, then the loops */
	{
		int pass;
		for (pass = 0; pass < 2; pass++)
			for (i = 0; i < ne; i++)
			{
				unsigned int e;
				if (ss_find(parent, i) != i || deg[i] == 0 || (pass == 0 && deg[i] == 2))
					continue;
				for (e = adj_start[i]; e < adj_start[i] + deg[i]; e++)
				{
					unsigned int node = i, s = adj[e], n = 0;
					int closed = 0;
					if (used[s])
						continue;
					line[n++] = SS_PT(i);
					for (;;)
					{
						unsigned int u = ss_find(parent, 2 * s), f;
						unsigned int next = u == node ? ss_find(parent, 2 * s + 1) : u;
						used[s] = 1;
						node = next;
						if (node == i)
						{
							closed = 1;
							break;
						}
						line[n++] = SS_PT(node);
						if (deg[node] != 2)
							break;
						s = (unsigned int)-1;
						for (f = adj_start[node]; f < adj_start[node] + deg[node]; f++)
							if (!used[adj[f]])
								s = adj[f];
						if (s == (unsigned int)-1)
							break;
					}
					ss_emit(o, segs[0].a, segs[0].b, line, n, closed);
				}
			}
	}
#undef SS_PT
done:
	free(ends);
	free(parent);
	free(deg);
	free(adj_start);
	free(adj);
	free(used);
	free(line);
	return st;
}

/* the grid over the prepared patches, segments, gap closing, chaining */
static qaws_status ss_solve(ss_ctx* x, ss_out* o)
{
	qaws_bp_box* boxes;
	qaws_bp_stats bs;
	unsigned int i, k, g = 0;
	qaws_status s = QAWS_STATUS_OK;
	x->stats.patch_count = x->npatch;
	if (x->npatch > 1)
	{
		boxes = (qaws_bp_box*)malloc(x->npatch * sizeof(qaws_bp_box));
		if (!boxes)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		for (i = 0; i < x->npatch; i++)
			for (k = 0; k < 3; k++)
			{
				boxes[i].lo[k] = (double)x->patches[i].lo[k];
				boxes[i].hi[k] = (double)x->patches[i].hi[k];
			}
		s = qaws_internal_broadphase(boxes, x->npatch, 3, ss_accept, ss_visit, x, &bs);
		x->stats.cell_count = bs.cell_count;
		x->stats.candidate_count = bs.candidate_count;
		free(boxes);
	}
	if (s != QAWS_STATUS_OK)
		return s;
	qsort(x->segs, x->nseg, sizeof(ss_seg), ss_cmp_seg);
	for (i = 1; i <= x->nseg && s == QAWS_STATUS_OK; i++)
		if (i == x->nseg || x->segs[i].a != x->segs[g].a || x->segs[i].b != x->segs[g].b)
		{
			/* the group's segments, gaps closed by marching, then chained */
			unsigned int m = i - g, cap = m;
			ss_seg* L = (ss_seg*)malloc(m * sizeof(ss_seg));
			if (!L)
				return QAWS_STATUS_ALLOCATION_FAILURE;
			memcpy(L, x->segs + g, m * sizeof(ss_seg));
			s = ss_close_gaps(x, &L, &m, &cap);
			if (s == QAWS_STATUS_OK)
				s = ss_chain(x, L, m, o);
			free(L);
			g = i;
		}
	x->stats.hit_count = o->npoint;
	return s;
}

static qaws_status ss_check_out(unsigned int* out_curve_count, unsigned int* out_point_count, qaws_surface_batch_curve* out_curves,
	unsigned int curve_capacity, qaws_ssi_point* out_points, unsigned int point_capacity)
{
	if (!out_curve_count || !out_point_count || (!out_curves && curve_capacity) || (!out_points && point_capacity))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve_count = 0;
	*out_point_count = 0;
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_batch_find_intersections(
	qaws_surface_batch_desc const* desc,
	qaws_surface_batch_curve* out_curves,
	unsigned int curve_capacity,
	unsigned int* out_curve_count,
	qaws_ssi_point* out_points,
	unsigned int point_capacity,
	unsigned int* out_point_count,
	qaws_surface_batch_stats* out_stats)
{
	ss_ctx x;
	ss_out o;
	qaws_scalar ext;
	unsigned int i;
	qaws_status s;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	s = ss_check_out(out_curve_count, out_point_count, out_curves, curve_capacity, out_points, point_capacity);
	if (s != QAWS_STATUS_OK || !desc || (desc->surface_count && !desc->surfaces))
		return s != QAWS_STATUS_OK ? s : QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->surface_count; i++)
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&x, 0, sizeof(x));
	memset(&o, 0, sizeof(o));
	x.desc = desc;
	ext = qaws_internal_flatten_extent(NULL, 0, 3, desc->surfaces, desc->surface_count);
	x.flat = desc->flatness > 0 ? desc->flatness : ext / 512;
	x.pos_tol = ext * SB_POS_REL;
	for (i = 0; i < desc->surface_count && s == QAWS_STATUS_OK; i++)
		s = qaws_internal_flatten_surface(desc->surfaces[i], x.flat, i, &x.patches, &x.npatch, &x.cappatch);
	o.curves = out_curves;
	o.curve_cap = curve_capacity;
	o.points = out_points;
	o.point_cap = point_capacity;
	if (s == QAWS_STATUS_OK)
		s = ss_solve(&x, &o);
	if (s == QAWS_STATUS_OK)
	{
		*out_curve_count = o.ncurve;
		*out_point_count = o.npoint;
	}
	if (out_stats)
		*out_stats = x.stats;
	free(x.patches);
	free(x.segs);
	return s;
}

qaws_status qaws_surface_set_find_intersections(
	qaws_surface_set const* set,
	qaws_surface_set const* other,
	qaws_surface_batch_curve* out_curves,
	unsigned int curve_capacity,
	unsigned int* out_curve_count,
	qaws_ssi_point* out_points,
	unsigned int point_capacity,
	unsigned int* out_point_count,
	qaws_surface_batch_stats* out_stats)
{
	ss_ctx x;
	ss_out o;
	qaws_surface_batch_desc d;
	qaws_surface const** surfs = NULL;
	unsigned int* fam = NULL;
	qaws_flat_patch* patches = NULL;
	unsigned int na, nb, i;
	qaws_status s;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	s = ss_check_out(out_curve_count, out_point_count, out_curves, curve_capacity, out_points, point_capacity);
	if (s != QAWS_STATUS_OK || !set)
		return s != QAWS_STATUS_OK ? s : QAWS_STATUS_INVALID_ARGUMENT;
	memset(&x, 0, sizeof(x));
	memset(&o, 0, sizeof(o));
	na = set->desc.surface_count;
	nb = other ? other->desc.surface_count : 0;
	if (!other)
	{
		/* the prepared patches as they are */
		x.desc = &set->desc;
		x.patches = set->patches;
		x.npatch = set->npatch;
		x.flat = set->flat;
		x.pos_tol = set->ext * SB_POS_REL;
	}
	else
	{
		/* both sets side by side: other's surfaces numbered after set's, in another family */
		surfs = (qaws_surface const**)malloc((na + nb ? na + nb : 1) * sizeof(qaws_surface*));
		fam = (unsigned int*)malloc((na + nb ? na + nb : 1) * sizeof(unsigned int));
		patches = (qaws_flat_patch*)malloc((set->npatch + other->npatch ? set->npatch + other->npatch : 1) * sizeof(qaws_flat_patch));
		if (!surfs || !fam || !patches)
		{
			s = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		memcpy((void*)surfs, set->desc.surfaces, na * sizeof(qaws_surface*));
		memcpy((void*)(surfs + na), other->desc.surfaces, nb * sizeof(qaws_surface*));
		for (i = 0; i < na + nb; i++)
			fam[i] = i >= na;
		memcpy(patches, set->patches, set->npatch * sizeof(qaws_flat_patch));
		memcpy(patches + set->npatch, other->patches, other->npatch * sizeof(qaws_flat_patch));
		for (i = set->npatch; i < set->npatch + other->npatch; i++)
			patches[i].owner += na;
		memset(&d, 0, sizeof(d));
		d.surfaces = surfs;
		d.surface_count = na + nb;
		d.families = fam;
		x.desc = &d;
		x.patches = patches;
		x.npatch = set->npatch + other->npatch;
		x.flat = set->flat > other->flat ? set->flat : other->flat;
		x.pos_tol = (set->ext > other->ext ? set->ext : other->ext) * SB_POS_REL;
	}
	o.curves = out_curves;
	o.curve_cap = curve_capacity;
	o.points = out_points;
	o.point_cap = point_capacity;
	s = ss_solve(&x, &o);
	if (s == QAWS_STATUS_OK)
	{
		/* across two sets, surface_b indexes other */
		if (other)
			for (i = 0; i < o.ncurve && i < curve_capacity; i++)
				out_curves[i].surface_b -= na;
		*out_curve_count = o.ncurve;
		*out_point_count = o.npoint;
	}
done:
	if (out_stats)
		*out_stats = x.stats;
	free(x.segs);
	free((void*)surfs);
	free(fam);
	free(patches);
	return s;
}
