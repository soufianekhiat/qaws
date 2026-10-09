#include "qaws_exact_surface.h"
#include "qaws_exact_solve.h"
#include "../internal/qaws_internal_broadphase.h"
#include <stdlib.h>
#include <string.h>

/*
 * Certified intersections of N exact curves.
 *
 * Every span of every curve gets a sound box from its control points
 * (qaws_exact_span_box); the shared grid (qaws_internal_broadphase) gives the
 * span pairs of different curves whose boxes overlap, and each one is
 * certified by qaws_exact_span_pair_hits. A curve pair with a span pair that
 * fails certification is dropped as a whole and counted.
 */

typedef struct eb_span
{
	unsigned int curve, index;
} eb_span;

typedef struct eb_ctx
{
	qaws_exact_batch_desc const* desc;
	eb_span* spans;
	qaws_exact_batch_hit* hits;
	unsigned int nhit, caphit;
	uint64_t* failed;   /* curve_a * curve_count + curve_b */
	unsigned int nfail, capfail;
	qaws_exact_pair* tmp;
} eb_ctx;

static int eb_grow(void** p, unsigned int* cap, unsigned int need, size_t elem)
{
	unsigned int c;
	void* q;
	if (need <= *cap)
		return 1;
	c = *cap ? *cap : 64;
	while (c < need)
		c *= 2;
	q = realloc(*p, (size_t)c * elem);
	if (!q)
		return 0;
	*p = q;
	*cap = c;
	return 1;
}

static qaws_status eb_fail(eb_ctx* x, unsigned int a, unsigned int b)
{
	if (!eb_grow((void**)&x->failed, &x->capfail, x->nfail + 1, sizeof(uint64_t)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	x->failed[x->nfail++] = (uint64_t)a * x->desc->curve_count + b;
	return QAWS_STATUS_OK;
}

static qaws_status eb_add(eb_ctx* x, unsigned int a, unsigned int b, qaws_exact_pair const* p, unsigned int n)
{
	unsigned int k;
	if (!eb_grow((void**)&x->hits, &x->caphit, x->nhit + n, sizeof(qaws_exact_batch_hit)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (k = 0; k < n; k++)
	{
		qaws_exact_batch_hit* h = &x->hits[x->nhit++];
		h->curve_a = a;
		h->curve_b = b;
		h->pair = p[k];
	}
	return QAWS_STATUS_OK;
}

static int eb_accept(void* user, unsigned int i, unsigned int j)
{
	eb_ctx const* x = (eb_ctx const*)user;
	unsigned int a = x->spans[i].curve, b = x->spans[j].curve;
	return a != b && (!x->desc->families || x->desc->families[a] != x->desc->families[b]);
}

static qaws_status eb_visit(void* user, unsigned int i, unsigned int j)
{
	eb_ctx* x = (eb_ctx*)user;
	eb_span const* A = &x->spans[i];
	eb_span const* B = &x->spans[j];
	unsigned int n = 0;
	int common = 0;
	qaws_status st = qaws_exact_span_pair_hits(x->desc->curves[A->curve], A->index, x->desc->curves[B->curve], B->index,
		x->tmp, QAWS_EXACT_SOLVE_MAX_ROOTS, &n, &common);
	if (st == QAWS_STATUS_OK)
		return eb_add(x, A->curve, B->curve, x->tmp, n);
	if (st == QAWS_STATUS_ALLOCATION_FAILURE)
		return st;
	return eb_fail(x, A->curve, B->curve);
}

static int eb_cmp_u64(void const* p, void const* q)
{
	uint64_t a = *(uint64_t const*)p, b = *(uint64_t const*)q;
	return a < b ? -1 : (a > b ? 1 : 0);
}

static int eb_cmp_hit(void const* p, void const* q)
{
	qaws_exact_batch_hit const* a = (qaws_exact_batch_hit const*)p;
	qaws_exact_batch_hit const* b = (qaws_exact_batch_hit const*)q;
	if (a->curve_a != b->curve_a) return a->curve_a < b->curve_a ? -1 : 1;
	if (a->curve_b != b->curve_b) return a->curve_b < b->curve_b ? -1 : 1;
	if (a->pair.a_lo != b->pair.a_lo) return a->pair.a_lo < b->pair.a_lo ? -1 : 1;
	if (a->pair.b_lo != b->pair.b_lo) return a->pair.b_lo < b->pair.b_lo ? -1 : 1;
	return 0;
}

/* each curve against itself */
static qaws_status eb_self(eb_ctx* x, unsigned int i)
{
	unsigned int cap = 64, n = 0;
	qaws_exact_pair* buf = NULL;
	qaws_status st;
	for (;;)
	{
		qaws_exact_pair* g = (qaws_exact_pair*)realloc(buf, cap * sizeof(qaws_exact_pair));
		if (!g)
		{
			free(buf);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		buf = g;
		st = qaws_exact_curve_self_hits(x->desc->curves[i], buf, cap, &n);
		if (st != QAWS_STATUS_BUFFER_TOO_SMALL)
			break;
		cap *= 4;
	}
	if (st == QAWS_STATUS_OK)
		st = eb_add(x, i, i, buf, n);
	else if (st != QAWS_STATUS_ALLOCATION_FAILURE)
		st = eb_fail(x, i, i);
	free(buf);
	return st;
}

qaws_status qaws_exact_curve_batch_hits(qaws_exact_batch_desc const* desc, qaws_exact_batch_hit* out_hits, unsigned int capacity,
	unsigned int* out_count, qaws_exact_batch_stats* out_stats)
{
	eb_ctx x;
	qaws_bp_box* boxes = NULL;
	qaws_bp_stats bs;
	qaws_exact_batch_stats stats;
	unsigned int i, k, nspan = 0, dim, n;
	qaws_status st = QAWS_STATUS_OK;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (!desc || !out_count || (!out_hits && capacity) || (desc->curve_count && !desc->curves))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (desc->curve_count == 0)
		return QAWS_STATUS_OK;
	for (i = 0; i < desc->curve_count; i++)
	{
		qaws_exact_curve const* c = desc->curves[i];
		if (!c || (c->dimension != 2 && c->dimension != 3) || c->dimension != desc->curves[0]->dimension)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (c->space_exp2 != desc->curves[0]->space_exp2)
			return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
		nspan += c->span_count;
	}
	dim = (unsigned int)desc->curves[0]->dimension;
	memset(&x, 0, sizeof(x));
	memset(&stats, 0, sizeof(stats));
	x.desc = desc;
	x.spans = (eb_span*)malloc((nspan ? nspan : 1) * sizeof(eb_span));
	boxes = (qaws_bp_box*)malloc((nspan ? nspan : 1) * sizeof(qaws_bp_box));
	x.tmp = (qaws_exact_pair*)malloc(QAWS_EXACT_SOLVE_MAX_ROOTS * sizeof(qaws_exact_pair));
	if (!x.spans || !boxes || !x.tmp)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	nspan = 0;
	for (i = 0; i < desc->curve_count; i++)
		for (k = 0; k < desc->curves[i]->span_count; k++, nspan++)
		{
			x.spans[nspan].curve = i;
			x.spans[nspan].index = k;
			qaws_exact_span_box(&desc->curves[i]->spans[k], dim, boxes[nspan].lo, boxes[nspan].hi);
		}
	stats.span_count = nspan;

	st = qaws_internal_broadphase(boxes, nspan, dim, eb_accept, eb_visit, &x, &bs);
	stats.cell_count = bs.cell_count;
	stats.candidate_count = bs.candidate_count;
	if (desc->flags & QAWS_EXACT_BATCH_SELF)
		for (i = 0; i < desc->curve_count && st == QAWS_STATUS_OK; i++)
			st = eb_self(&x, i);
	if (st != QAWS_STATUS_OK)
		goto done;

	/* drop the hits of curve pairs that failed, then sort */
	qsort(x.failed, x.nfail, sizeof(uint64_t), eb_cmp_u64);
	for (i = 0, k = 0; i < x.nfail; i++)
		if (k == 0 || x.failed[k - 1] != x.failed[i])
			x.failed[k++] = x.failed[i];
	x.nfail = k;
	stats.uncertified_count = x.nfail;
	for (i = 0, n = 0; i < x.nhit; i++)
	{
		uint64_t key = (uint64_t)x.hits[i].curve_a * desc->curve_count + x.hits[i].curve_b;
		if (!x.nfail || !bsearch(&key, x.failed, x.nfail, sizeof(uint64_t), eb_cmp_u64))
			x.hits[n++] = x.hits[i];
	}
	x.nhit = n;
	qsort(x.hits, x.nhit, sizeof(qaws_exact_batch_hit), eb_cmp_hit);
	for (i = 0; i < x.nhit && i < capacity; i++)
		out_hits[i] = x.hits[i];
	*out_count = x.nhit;
	stats.hit_count = x.nhit;
	st = x.nhit > capacity ? QAWS_STATUS_BUFFER_TOO_SMALL : (x.nfail ? QAWS_STATUS_CERTIFICATION_FAILED : QAWS_STATUS_OK);
done:
	if (out_stats)
		*out_stats = stats;
	free(x.spans);
	free(x.hits);
	free(x.failed);
	free(x.tmp);
	free(boxes);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Curves x surfaces                                                  */
/* ------------------------------------------------------------------ */

typedef struct es_piece
{
	unsigned int owner, a, b;   /* span: curve, span index; patch: surface, iu, iv */
} es_piece;

typedef struct es_ctx
{
	qaws_exact_surface_batch_desc const* desc;
	es_piece* pieces;
	unsigned int nspan;         /* pieces [0, nspan) are spans, then patches */
	qaws_exact_curve_surface_batch_hit* hits;
	unsigned int nhit, caphit;
	uint64_t* failed;           /* curve * surface_count + surface */
	unsigned int nfail, capfail;
	qaws_exact_curve_surface_hit tmp[64];
} es_ctx;

static int es_accept(void* user, unsigned int i, unsigned int j)
{
	es_ctx const* x = (es_ctx const*)user;
	return i < x->nspan && j >= x->nspan;
}

static qaws_status es_visit(void* user, unsigned int i, unsigned int j)
{
	es_ctx* x = (es_ctx*)user;
	es_piece const* S = &x->pieces[i];
	es_piece const* P = &x->pieces[j];
	unsigned int n = 0, k;
	qaws_status st = qaws_exact_span_patch_hits(x->desc->curves[S->owner], S->a, x->desc->surfaces[P->owner], P->a, P->b, x->tmp, 64, &n);
	if (st == QAWS_STATUS_ALLOCATION_FAILURE)
		return st;
	if (st != QAWS_STATUS_OK)
	{
		if (!eb_grow((void**)&x->failed, &x->capfail, x->nfail + 1, sizeof(uint64_t)))
			return QAWS_STATUS_ALLOCATION_FAILURE;
		x->failed[x->nfail++] = (uint64_t)S->owner * x->desc->surface_count + P->owner;
		return QAWS_STATUS_OK;
	}
	if (!eb_grow((void**)&x->hits, &x->caphit, x->nhit + n, sizeof(qaws_exact_curve_surface_batch_hit)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (k = 0; k < n; k++)
	{
		qaws_exact_curve_surface_batch_hit* h = &x->hits[x->nhit++];
		h->curve = S->owner;
		h->surface = P->owner;
		h->hit = x->tmp[k];
	}
	return QAWS_STATUS_OK;
}

static int es_cmp_hit(void const* p, void const* q)
{
	qaws_exact_curve_surface_batch_hit const* a = (qaws_exact_curve_surface_batch_hit const*)p;
	qaws_exact_curve_surface_batch_hit const* b = (qaws_exact_curve_surface_batch_hit const*)q;
	if (a->curve != b->curve) return a->curve < b->curve ? -1 : 1;
	if (a->surface != b->surface) return a->surface < b->surface ? -1 : 1;
	if (a->hit.t_lo != b->hit.t_lo) return a->hit.t_lo < b->hit.t_lo ? -1 : 1;
	return 0;
}

qaws_status qaws_exact_curve_surface_batch_hits(qaws_exact_surface_batch_desc const* desc, qaws_exact_curve_surface_batch_hit* out_hits,
	unsigned int capacity, unsigned int* out_count, qaws_exact_batch_stats* out_stats)
{
	es_ctx* x;
	qaws_bp_box* boxes = NULL;
	qaws_bp_stats bs;
	qaws_exact_batch_stats stats;
	unsigned int i, k, iu, iv, npiece = 0, nspan = 0, n, group;
	int space = 0;
	qaws_status st = QAWS_STATUS_OK;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (!desc || !out_count || (!out_hits && capacity) || (desc->curve_count && !desc->curves) || (desc->surface_count && !desc->surfaces))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	memset(&stats, 0, sizeof(stats));
	for (i = 0; i < desc->curve_count; i++)
		if (!desc->curves[i] || desc->curves[i]->dimension != 3)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < desc->surface_count; i++)
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	/* one exact space for everything */
	space = desc->curve_count ? desc->curves[0]->space_exp2 : (desc->surface_count ? desc->surfaces[0]->space_exp2 : 0);
	for (i = 0; i < desc->curve_count; i++)
	{
		if (desc->curves[i]->space_exp2 != space)
			return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
		nspan += desc->curves[i]->span_count;
	}
	npiece = nspan;
	for (i = 0; i < desc->surface_count; i++)
	{
		if (desc->surfaces[i]->space_exp2 != space)
			return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
		npiece += desc->surfaces[i]->nu * desc->surfaces[i]->nv;
	}
	stats.span_count = nspan;
	stats.patch_count = npiece - nspan;
	if (!nspan || npiece == nspan)
	{
		if (out_stats)
			*out_stats = stats;
		return QAWS_STATUS_OK;
	}
	x = (es_ctx*)calloc(1, sizeof(es_ctx));
	boxes = (qaws_bp_box*)malloc(npiece * sizeof(qaws_bp_box));
	if (!x || !boxes)
	{
		free(x);
		free(boxes);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	x->desc = desc;
	x->nspan = nspan;
	x->pieces = (es_piece*)malloc(npiece * sizeof(es_piece));
	if (!x->pieces)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	npiece = 0;
	for (i = 0; i < desc->curve_count; i++)
		for (k = 0; k < desc->curves[i]->span_count; k++, npiece++)
		{
			x->pieces[npiece].owner = i;
			x->pieces[npiece].a = k;
			x->pieces[npiece].b = 0;
			qaws_exact_span_box(&desc->curves[i]->spans[k], 3, boxes[npiece].lo, boxes[npiece].hi);
		}
	for (i = 0; i < desc->surface_count; i++)
		for (iu = 0; iu < desc->surfaces[i]->nu; iu++)
			for (iv = 0; iv < desc->surfaces[i]->nv; iv++, npiece++)
			{
				x->pieces[npiece].owner = i;
				x->pieces[npiece].a = iu;
				x->pieces[npiece].b = iv;
				qaws_exact_patch_box(desc->surfaces[i], iu, iv, boxes[npiece].lo, boxes[npiece].hi);
			}
	st = qaws_internal_broadphase(boxes, npiece, 3, es_accept, es_visit, x, &bs);
	stats.cell_count = bs.cell_count;
	stats.candidate_count = bs.candidate_count;
	if (st != QAWS_STATUS_OK)
		goto done;

	/* drop failed pairs; sort; a root on a span or patch edge once */
	qsort(x->failed, x->nfail, sizeof(uint64_t), eb_cmp_u64);
	for (i = 0, k = 0; i < x->nfail; i++)
		if (k == 0 || x->failed[k - 1] != x->failed[i])
			x->failed[k++] = x->failed[i];
	x->nfail = k;
	stats.uncertified_count = x->nfail;
	for (i = 0, n = 0; i < x->nhit; i++)
	{
		uint64_t key = (uint64_t)x->hits[i].curve * desc->surface_count + x->hits[i].surface;
		if (!x->nfail || !bsearch(&key, x->failed, x->nfail, sizeof(uint64_t), eb_cmp_u64))
			x->hits[n++] = x->hits[i];
	}
	x->nhit = n;
	qsort(x->hits, x->nhit, sizeof(qaws_exact_curve_surface_batch_hit), es_cmp_hit);
	for (i = 0, n = 0, group = 0; i < x->nhit; i++)
	{
		unsigned int r;
		int dup = 0;
		if (n == 0 || x->hits[n - 1].curve != x->hits[i].curve || x->hits[n - 1].surface != x->hits[i].surface)
			group = n;
		for (r = group; r < n && !dup; r++)
			dup = qaws_exact_curve_surface_hits_overlap(&x->hits[r].hit, &x->hits[i].hit);
		if (!dup)
			x->hits[n++] = x->hits[i];
	}
	x->nhit = n;
	for (i = 0; i < x->nhit && i < capacity; i++)
		out_hits[i] = x->hits[i];
	*out_count = x->nhit;
	stats.hit_count = x->nhit;
	st = x->nhit > capacity ? QAWS_STATUS_BUFFER_TOO_SMALL : (x->nfail ? QAWS_STATUS_CERTIFICATION_FAILED : QAWS_STATUS_OK);
done:
	if (out_stats)
		*out_stats = stats;
	free(x->pieces);
	free(x->hits);
	free(x->failed);
	free(x);
	free(boxes);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Surfaces x surfaces                                                */
/* ------------------------------------------------------------------ */

typedef struct ess_pair
{
	unsigned int a, b;          /* surfaces, a < b */
	unsigned int q[4];          /* iu1, iv1, iu2, iv2 */
} ess_pair;

typedef struct ess_ctx
{
	qaws_exact_ssi_batch_desc const* desc;
	es_piece* pieces;
	ess_pair* pairs;
	unsigned int npair, cappair;
} ess_ctx;

static int ess_accept(void* user, unsigned int i, unsigned int j)
{
	ess_ctx const* x = (ess_ctx const*)user;
	unsigned int a = x->pieces[i].owner, b = x->pieces[j].owner;
	return a != b && (!x->desc->families || x->desc->families[a] != x->desc->families[b]);
}

static qaws_status ess_visit(void* user, unsigned int i, unsigned int j)
{
	ess_ctx* x = (ess_ctx*)user;
	es_piece const* P = &x->pieces[i];
	es_piece const* Q = &x->pieces[j];
	ess_pair* p;
	if (!eb_grow((void**)&x->pairs, &x->cappair, x->npair + 1, sizeof(ess_pair)))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* pieces are numbered surface by surface: i < j gives P's surface first */
	p = &x->pairs[x->npair++];
	p->a = P->owner;
	p->b = Q->owner;
	p->q[0] = P->a;
	p->q[1] = P->b;
	p->q[2] = Q->a;
	p->q[3] = Q->b;
	return QAWS_STATUS_OK;
}

static int ess_cmp(void const* p, void const* q)
{
	ess_pair const* a = (ess_pair const*)p;
	ess_pair const* b = (ess_pair const*)q;
	unsigned int k;
	if (a->a != b->a) return a->a < b->a ? -1 : 1;
	if (a->b != b->b) return a->b < b->b ? -1 : 1;
	for (k = 0; k < 4; k++)
		if (a->q[k] != b->q[k])
			return a->q[k] < b->q[k] ? -1 : 1;
	return 0;
}

qaws_status qaws_exact_surface_batch_hits(qaws_exact_ssi_batch_desc const* desc, qaws_exact_ssi_point* out_points, unsigned int point_capacity,
	unsigned int* out_point_count, qaws_exact_ssi_batch_branch* out_branches, unsigned int branch_capacity, unsigned int* out_branch_count,
	qaws_exact_batch_stats* out_stats)
{
	ess_ctx x;
	qaws_bp_box* boxes = NULL;
	qaws_bp_stats bs;
	qaws_exact_batch_stats stats;
	unsigned int* list = NULL;
	qaws_exact_ssi_branch* tmp = NULL;
	unsigned int i, k, iu, iv, npiece = 0, np = 0, nb = 0, g;
	int space;
	qaws_status st = QAWS_STATUS_OK;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (!desc || !out_point_count || !out_branch_count || (!out_points && point_capacity) || (!out_branches && branch_capacity)
		|| (desc->surface_count && !desc->surfaces))
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_point_count = 0;
	*out_branch_count = 0;
	memset(&stats, 0, sizeof(stats));
	for (i = 0; i < desc->surface_count; i++)
		if (!desc->surfaces[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	space = desc->surface_count ? desc->surfaces[0]->space_exp2 : 0;
	for (i = 0; i < desc->surface_count; i++)
	{
		if (desc->surfaces[i]->space_exp2 != space)
			return QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE;
		npiece += desc->surfaces[i]->nu * desc->surfaces[i]->nv;
	}
	stats.patch_count = npiece;
	memset(&x, 0, sizeof(x));
	x.desc = desc;
	x.pieces = (es_piece*)malloc((npiece ? npiece : 1) * sizeof(es_piece));
	boxes = (qaws_bp_box*)malloc((npiece ? npiece : 1) * sizeof(qaws_bp_box));
	tmp = (qaws_exact_ssi_branch*)malloc((branch_capacity ? branch_capacity : 1) * sizeof(qaws_exact_ssi_branch));
	if (!x.pieces || !boxes || !tmp)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	npiece = 0;
	for (i = 0; i < desc->surface_count; i++)
		for (iu = 0; iu < desc->surfaces[i]->nu; iu++)
			for (iv = 0; iv < desc->surfaces[i]->nv; iv++, npiece++)
			{
				x.pieces[npiece].owner = i;
				x.pieces[npiece].a = iu;
				x.pieces[npiece].b = iv;
				qaws_exact_patch_box(desc->surfaces[i], iu, iv, boxes[npiece].lo, boxes[npiece].hi);
			}
	st = qaws_internal_broadphase(boxes, npiece, 3, ess_accept, ess_visit, &x, &bs);
	stats.cell_count = bs.cell_count;
	stats.candidate_count = bs.candidate_count;
	if (st != QAWS_STATUS_OK)
		goto done;
	qsort(x.pairs, x.npair, sizeof(ess_pair), ess_cmp);
	list = (unsigned int*)malloc((x.npair ? x.npair : 1) * 4 * sizeof(unsigned int));
	if (!list)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* one certified solve per surface pair, over its candidate patch pairs */
	for (g = 0, i = 1; i <= x.npair && st == QAWS_STATUS_OK; i++)
	{
		unsigned int n = 0, pc = 0, bc = 0, a, b;
		qaws_status s2;
		if (i < x.npair && x.pairs[i].a == x.pairs[g].a && x.pairs[i].b == x.pairs[g].b)
			continue;
		a = x.pairs[g].a;
		b = x.pairs[g].b;
		for (k = g; k < i; k++, n++)
			memcpy(&list[4 * n], x.pairs[k].q, 4 * sizeof(unsigned int));
		g = i;
		s2 = qaws_exact_ssi_solve(desc->surfaces[a], desc->surfaces[b], desc->min_depth, list, n, out_points + np, point_capacity - np, &pc,
			tmp, branch_capacity - nb, &bc);
		if (s2 == QAWS_STATUS_BUFFER_TOO_SMALL || s2 == QAWS_STATUS_ALLOCATION_FAILURE)
		{
			st = s2;
			break;
		}
		if (s2 != QAWS_STATUS_OK)
		{
			stats.uncertified_count++;
			continue;
		}
		for (k = 0; k < bc; k++)
		{
			out_branches[nb + k].surface_a = a;
			out_branches[nb + k].surface_b = b;
			out_branches[nb + k].branch = tmp[k];
			out_branches[nb + k].branch.first += np;
		}
		np += pc;
		nb += bc;
	}
	*out_point_count = np;
	*out_branch_count = nb;
	stats.hit_count = np;
	if (st == QAWS_STATUS_OK && stats.uncertified_count)
		st = QAWS_STATUS_CERTIFICATION_FAILED;
done:
	if (out_stats)
		*out_stats = stats;
	free(x.pieces);
	free(x.pairs);
	free(boxes);
	free(list);
	free(tmp);
	return st;
}
