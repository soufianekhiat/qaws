#include "qaws_exact_curve.h"
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
