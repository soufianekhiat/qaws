#include "qaws_internal_broadphase.h"
#include "qaws_internal_memory.h"
#include "qaws_internal_parallel.h"
#include "../qaws_platform.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define BP_MAX_CELLS (1u << 22)
#if QAWS_SCALAR_IS_FLOAT
#define BP_FLAT 1e-5          /* an axis this thin against the widest is flat */
#else
#define BP_FLAT 1e-9
#endif
#define BP_MAX_AXIS  4096u

typedef struct bp_grid
{
	double lo[3], inv[3];
	unsigned int g[3];
} bp_grid;

static unsigned int bp_cell(bp_grid const* G, double v, unsigned int k)
{
	double c = (v - G->lo[k]) * G->inv[k];
	if (!(c > 0))
		return 0;
	if (c >= (double)G->g[k])
		return G->g[k] - 1;
	return (unsigned int)c;
}

static int bp_unbounded(qaws_bp_box const* b, unsigned int dim)
{
	unsigned int k;
	for (k = 0; k < dim; k++)
		if (!(b->lo[k] <= b->hi[k]))
			return 1;
	return 0;
}

/* cell ranges of box b on each axis */
static void bp_range(bp_grid const* G, qaws_bp_box const* b, unsigned int dim, unsigned int* c0, unsigned int* c1)
{
	unsigned int k;
	int unb = bp_unbounded(b, dim);
	for (k = 0; k < 3; k++)
	{
		if (unb || k >= dim)
		{
			c0[k] = 0;
			c1[k] = G->g[k] - 1;
		}
		else
		{
			c0[k] = bp_cell(G, b->lo[k], k);
			c1[k] = bp_cell(G, b->hi[k], k);
		}
	}
}

/* A built grid: cell ranges of every box filled by counting sort. */
struct qaws_bp_grid
{
	bp_grid G;
	unsigned int dim, ncell;
	unsigned int count;          /* boxes */
	unsigned int* start;        /* ncell + 1 */
	unsigned int* items;
	double hmin;                /* smallest cell side over the axes with extent */
	int fit;                    /* cells sized by the boxes (no sweep needed) */
	unsigned int nchunk;        /* the pair walk's chunks, about equal work */
	unsigned int* chunk;        /* nchunk + 1 starts: cell, row (a crowded cell's rows split among chunks) */
};

#define BP_PAIR_GRAIN 64   /* cells per chunk */
#define BP_BOX_GRAIN 16384 /* boxes per chunk of the build */
#define BP_SLABS 16        /* layers of cells filled apart through an executor */
#define BP_FIT_CELLS 16    /* a fitted grid's cells per box at most */

/* the pair walk cut into nchunk runs of about equal work: a cell weighs one, and each of its rows (its boxes in order)
   one per box after it; rows of g0 cells (row_work: their weights) skipped whole when no run starts in them */
static unsigned int* bp_chunks(unsigned int const* start, unsigned int g0, unsigned int nrow, double const* row_work, unsigned int nchunk,
	double total)
{
	unsigned int* at = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nchunk + 1));
	double step = total / (nchunk ? nchunk : 1), acc = 0;
	unsigned int r, c, e, k = 1;
	if (!at)
		return NULL;
	at[0] = at[1] = 0;
	for (r = 0; r < nrow && k < nchunk; r++)
	{
		if (acc + row_work[r] < k * step)
		{
			acc += row_work[r];
			continue;
		}
		for (c = r * g0; c < (r + 1) * g0; c++)
		{
			unsigned int n = start[c + 1] - start[c];
			for (; k < nchunk && acc >= k * step; k++)
			{
				at[2 * k] = c;
				at[2 * k + 1] = 0;
			}
			acc += 1;
			for (e = 0; e < n; e++)
			{
				for (; k < nchunk && acc >= k * step; k++)
				{
					at[2 * k] = c;
					at[2 * k + 1] = e;
				}
				acc += n - 1 - e;
			}
		}
	}
	for (; k <= nchunk; k++)
	{
		at[2 * k] = nrow * g0;
		at[2 * k + 1] = 0;
	}
	return at;
}

/* the build's state: boxes in chunks for the bounds and cell ranges, then cells in slabs of whole layers (each slab
   counts and fills its own cells, its boxes in order: the same grid whatever the slabs) */
typedef struct bp_make
{
	qaws_bp_box const* boxes;
	unsigned int count, dim, nslab;
	bp_grid G;
	double* part;                      /* per box chunk: lo[3], hi[3], widest sides summed, bounded boxes */
	unsigned short* range;             /* per box: first and last cells on x, y, z */
	unsigned int* start;
	unsigned int* items;
	unsigned int* slab;                /* per slab: its items, then the items before it */
	double* row_work;                  /* per row of cells along x: the pair walk's work */
} bp_make;

static qaws_status bp_bounds_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	bp_make const* M = (bp_make const*)user;
	double* p = M->part + 8 * (size_t)chunk;
	unsigned int i, k;
	for (k = 0; k < 3; k++)
	{
		p[k] = HUGE_VAL;
		p[3 + k] = -HUGE_VAL;
	}
	p[6] = p[7] = 0;
	for (i = begin; i < end; i++)
		if (!bp_unbounded(&M->boxes[i], M->dim))
		{
			double w = 0;
			for (k = 0; k < M->dim; k++)
			{
				if (M->boxes[i].lo[k] < p[k]) p[k] = M->boxes[i].lo[k];
				if (M->boxes[i].hi[k] > p[3 + k]) p[3 + k] = M->boxes[i].hi[k];
				w = fmax(w, M->boxes[i].hi[k] - M->boxes[i].lo[k]);
			}
			p[6] += w;
			p[7] += 1;
		}
	return QAWS_STATUS_OK;
}

static qaws_status bp_range_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	bp_make const* M = (bp_make const*)user;
	unsigned int i, k;
	(void)chunk;
	for (i = begin; i < end; i++)
	{
		unsigned int c0[3], c1[3];
		bp_range(&M->G, &M->boxes[i], M->dim, c0, c1);
		for (k = 0; k < 3; k++)
		{
			M->range[6 * (size_t)i + 2 * k] = (unsigned short)c0[k];
			M->range[6 * (size_t)i + 2 * k + 1] = (unsigned short)c1[k];
		}
	}
	return QAWS_STATUS_OK;
}

/* slab s: its layers [z0, z1) */
static void bp_slab(bp_make const* M, unsigned int s, unsigned int* z0, unsigned int* z1)
{
	*z0 = (unsigned int)((unsigned long long)M->G.g[2] * s / M->nslab);
	*z1 = (unsigned int)((unsigned long long)M->G.g[2] * (s + 1) / M->nslab);
}

/* a slab's cells counted, then their starts from the slab's first item, and its rows' work */
static qaws_status bp_count_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	bp_make const* M = (bp_make const*)user;
	unsigned int s, i, a, b, c, z0, z1, g0 = M->G.g[0], g1 = M->G.g[1];
	(void)chunk;
	for (s = begin; s < end; s++)
	{
		unsigned int sum = 0;
		bp_slab(M, s, &z0, &z1);
		for (i = 0; i < M->count; i++)
		{
			unsigned short const* r = M->range + 6 * (size_t)i;
			unsigned int lo = r[4] > z0 ? r[4] : z0, hi = r[5] + 1u < z1 ? r[5] + 1u : z1;
			for (c = lo; c < hi; c++)
				for (b = r[2]; b <= r[3]; b++)
					for (a = r[0]; a <= r[1]; a++)
						M->start[(c * g1 + b) * g0 + a]++;
		}
		for (b = z0 * g1; b < z1 * g1; b++)
		{
			double w = 0;
			for (a = b * g0; a < (b + 1) * g0; a++)
			{
				unsigned int n = M->start[a];
				M->start[a] = sum;
				sum += n;
				w += 1 + 0.5 * n * ((double)n - 1);
			}
			M->row_work[b] = w;
		}
		M->slab[s] = sum;
	}
	return QAWS_STATUS_OK;
}

/* a slab's items, its boxes in order (each cell's start moved on as it fills, then put back) */
static qaws_status bp_fill_job(void* user, unsigned int chunk, unsigned int begin, unsigned int end)
{
	bp_make const* M = (bp_make const*)user;
	unsigned int s, i, a, b, c, z0, z1, g0 = M->G.g[0], g1 = M->G.g[1];
	(void)chunk;
	for (s = begin; s < end; s++)
	{
		bp_slab(M, s, &z0, &z1);
		if (z0 == z1)
			continue;
		for (c = z0 * g1 * g0; c < z1 * g1 * g0; c++)
			M->start[c] += M->slab[s];
		for (i = 0; i < M->count; i++)
		{
			unsigned short const* r = M->range + 6 * (size_t)i;
			unsigned int lo = r[4] > z0 ? r[4] : z0, hi = r[5] + 1u < z1 ? r[5] + 1u : z1;
			for (c = lo; c < hi; c++)
				for (b = r[2]; b <= r[3]; b++)
					for (a = r[0]; a <= r[1]; a++)
						M->items[M->start[(c * g1 + b) * g0 + a]++] = i;
		}
		for (c = z1 * g1 * g0 - 1; c > z0 * g1 * g0; c--)
			M->start[c] = M->start[c - 1];
		M->start[z0 * g1 * g0] = M->slab[s];
	}
	return QAWS_STATUS_OK;
}

/* fit: cells about the boxes' mean size (surfaces: the volume rule makes cells far bigger than their boxes) */
static qaws_status bp_build(qaws_bp_box const* boxes, unsigned int count, unsigned int dim, struct qaws_bp_grid* B, int fit,
	qaws_batch_executor const* executor)
{
	bp_make M;
	double hi[3], h = 0, work = 0, mean = 0;
	unsigned int i, k, s, ncell, total = 0, bounded = 0, nch = qaws_internal_chunk_count(count, BP_BOX_GRAIN);
	unsigned int* chunk = NULL;
	qaws_status st;
	memset(B, 0, sizeof(*B));
	memset(&M, 0, sizeof M);
	M.boxes = boxes;
	M.count = count;
	M.dim = dim;
	/* scene box of the bounded boxes (and their mean widest side), by chunks summed in order */
	M.part = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 8, (size_t)nch + 1));
	if (!M.part)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = qaws_internal_parallel(count, BP_BOX_GRAIN, bp_bounds_job, &M, executor);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (k = 0; k < 3; k++)
	{
		M.G.lo[k] = HUGE_VAL;
		hi[k] = -HUGE_VAL;
	}
	for (i = 0; i < nch; i++)
	{
		double const* p = M.part + 8 * (size_t)i;
		for (k = 0; k < dim; k++)
		{
			if (p[k] < M.G.lo[k]) M.G.lo[k] = p[k];
			if (p[3 + k] > hi[k]) hi[k] = p[3 + k];
		}
		mean += p[6];
		bounded += (unsigned int)p[7];
	}
	if (!bounded)
		for (k = 0; k < 3; k++)
		{
			M.G.lo[k] = 0;
			hi[k] = 0;
		}
	for (k = dim; k < 3; k++)
	{
		M.G.lo[k] = 0;
		hi[k] = 0;
	}

	/* cell size: about one bounded box per cell over the axes with extent */
	{
		double vol = 1, ext = 0;
		unsigned int axes = 0;
		for (k = 0; k < dim; k++)
			if (hi[k] - M.G.lo[k] > ext)
				ext = hi[k] - M.G.lo[k];
		/* (an extent within the boxes' rounding is flat: float boxes of a
		   curve on a line are a rounding thick) */
		for (k = 0; k < dim; k++)
			if (hi[k] - M.G.lo[k] > ext * BP_FLAT)
			{
				vol *= hi[k] - M.G.lo[k];
				axes++;
			}
		h = axes ? pow(vol / (bounded ? bounded : 1), 1.0 / axes) : 0;
		if (fit && bounded)
		{
			/* the boxes' mean widest side, made a little bigger: a box spans about two cells an axis */
			double cells = 1;
			mean = 1.5 * mean / bounded;
			if (mean > 0 && mean < h)
				h = mean;
			/* boxes along a curve fill few cells: not past BP_FIT_CELLS a box (the empty ones walked for nothing) */
			for (k = 0; k < dim; k++)
				if (hi[k] - M.G.lo[k] > ext * BP_FLAT)
					cells *= (hi[k] - M.G.lo[k]) / h;
			if (h > 0 && cells > BP_FIT_CELLS * (double)bounded)
				h *= pow(cells / (BP_FIT_CELLS * (double)bounded), 1.0 / axes);
		}
	}
	for (k = 0; k < 3; k++)
	{
		double n = k < dim && h > 0 ? ceil((hi[k] - M.G.lo[k]) / h) : 1;
		M.G.g[k] = n < 1 ? 1 : (n > BP_MAX_AXIS ? BP_MAX_AXIS : (unsigned int)n);
	}
	while ((double)M.G.g[0] * M.G.g[1] * M.G.g[2] > BP_MAX_CELLS)
		for (k = 0; k < 3; k++)
			M.G.g[k] = (M.G.g[k] + 1) / 2;
	for (k = 0; k < 3; k++)
		M.G.inv[k] = hi[k] > M.G.lo[k] ? (double)M.G.g[k] / (hi[k] - M.G.lo[k]) : 0;
	ncell = M.G.g[0] * M.G.g[1] * M.G.g[2];

	/* count, prefix, fill: whole layers per slab */
	M.nslab = executor && M.G.g[2] > 1 ? (M.G.g[2] < BP_SLABS ? M.G.g[2] : BP_SLABS) : 1;
	M.range = (unsigned short*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned short), 6, (size_t)count + 1));
	M.start = (unsigned int*)qaws_internal_zalloc((size_t)ncell + 1, sizeof(unsigned int));
	M.slab = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * M.nslab);
	M.row_work = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)M.G.g[1] * M.G.g[2]));
	st = !M.range || !M.start || !M.slab || !M.row_work ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_parallel(count, BP_BOX_GRAIN, bp_range_job, &M, executor);
	if (st == QAWS_STATUS_OK)
		st = qaws_internal_parallel(M.nslab, 1, bp_count_job, &M, executor);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (s = 0; s < M.nslab; s++)
	{
		unsigned int n = M.slab[s];
		M.slab[s] = total;
		total += n;
	}
	M.start[ncell] = total;
	for (i = 0; i < M.G.g[1] * M.G.g[2]; i++)
		work += M.row_work[i];
	M.items = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(((size_t)(total ? total : 1)), sizeof(unsigned int)));
	st = M.items ? qaws_internal_parallel(M.nslab, 1, bp_fill_job, &M, executor) : QAWS_STATUS_ALLOCATION_FAILURE;
	if (st != QAWS_STATUS_OK)
		goto done;
	B->nchunk = qaws_internal_chunk_count(ncell, BP_PAIR_GRAIN);
	chunk = bp_chunks(M.start, M.G.g[0], M.G.g[1] * M.G.g[2], M.row_work, B->nchunk, work);
	if (!chunk)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	B->G = M.G;
	B->dim = dim;
	B->count = count;
	B->ncell = ncell;
	B->start = M.start;
	B->items = M.items;
	B->chunk = chunk;
	B->fit = fit;
	M.start = M.items = NULL;
	/* (over the axes the rings advance along: a flat axis of one cell, its
	   side a rounding, would make the ring bound crawl) */
	B->hmin = HUGE_VAL;
	for (k = 0; k < dim; k++)
		if (M.G.inv[k] > 0 && M.G.g[k] > 1 && 1 / M.G.inv[k] < B->hmin)
			B->hmin = 1 / M.G.inv[k];
	if (!(B->hmin < HUGE_VAL))
		B->hmin = 0;
done:
	qaws_internal_free(M.part);
	qaws_internal_free(M.range);
	qaws_internal_free(M.start);
	qaws_internal_free(M.items);
	qaws_internal_free(M.slab);
	qaws_internal_free(M.row_work);
	return st;
}

#define BP_LOCAL 32

/* the pairs from row `from` to row `to` (cell, row: a cell's rows are its boxes in order, each paired with those after
   it), each in the cell holding the low corner of the overlap */
static qaws_status bp_cells(struct qaws_bp_grid const* GR, qaws_bp_box const* boxes, qaws_bp_accept accept, void* accept_user,
	qaws_status (*visit)(void* user, unsigned int chunk, unsigned int i, unsigned int j), void* user, unsigned int chunk,
	unsigned int const* from, unsigned int const* to, unsigned int* candidates)
{
	unsigned int i, k, dim = GR->dim, c1 = to[0] + (to[1] > 0);
	qaws_status s = QAWS_STATUS_OK;
	qaws_bp_box stack_box[BP_LOCAL], *local = stack_box;   /* a cell's boxes read once, not once per pair */
	unsigned int stack_cell[2 * BP_LOCAL], *lcell = stack_cell, *pick, most = 0;
	for (i = from[0]; i < c1; i++)
		if (GR->start[i + 1] - GR->start[i] > most)
			most = GR->start[i + 1] - GR->start[i];
	if (most > BP_LOCAL)
	{
		/* crowded cells (fans round a vertex): one buffer for the largest */
		local = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), most));
		lcell = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, most));
		if (!local || !lcell)
		{
			qaws_internal_free(local == stack_box ? NULL : local);
			qaws_internal_free(lcell == stack_cell ? NULL : lcell);
			local = stack_box;
			lcell = stack_cell;
			most = BP_LOCAL;   /* (the general path for the larger cells) */
		}
	}
	else
		most = BP_LOCAL;
	pick = lcell + most;
	for (i = from[0]; i < c1 && s == QAWS_STATUS_OK; i++)
	{
		unsigned int e, f, b0 = GR->start[i], n = GR->start[i + 1] - b0, here[3];
		unsigned int r0 = i == from[0] ? from[1] : 0, r1 = i == to[0] ? to[1] : n;
		int near = n <= most;
		if (n < 2 || r0 >= r1)
			continue;
		here[0] = i % GR->G.g[0];
		here[1] = i / GR->G.g[0] % GR->G.g[1];
		here[2] = i / GR->G.g[0] / GR->G.g[1];
		for (e = 0; near && e < n; e++)
		{
			local[e] = boxes[GR->items[b0 + e]];
			if (bp_unbounded(&local[e], dim))
				near = 0;   /* (the general path below) */
			else
				for (k = 0, lcell[e] = 0; k < 3; k++)
					lcell[e] |= (unsigned int)((k < dim ? bp_cell(&GR->G, local[e].lo[k], k) : 0) == here[k]) << k;
		}
		if (near)
		{
			/* bounded boxes: the overlap's low corner is in the cell of the larger low cells (cells grow with the
			   coordinate, a box's low cell is at most this one): on each axis, one of the two starts here (the
			   axes it does, bit by bit in lcell). The partners of e picked without branches, then visited in order. */
			for (e = r0; e < r1 && s == QAWS_STATUS_OK; e++)
			{
				qaws_bp_box const* A = &local[e];
				unsigned int m = 0, j;
				for (f = e + 1; f < n; f++)
				{
					qaws_bp_box const* B = &local[f];
					int over = (lcell[e] | lcell[f]) == 7;
					for (k = 0; k < dim; k++)
						over &= (A->lo[k] <= B->hi[k]) & (B->lo[k] <= A->hi[k]);
					pick[m] = f;
					m += (unsigned int)over;
				}
				for (j = 0; j < m && s == QAWS_STATUS_OK; j++)
				{
					unsigned int p = GR->items[b0 + e], q = GR->items[b0 + pick[j]], lo_i, hi_i;
					lo_i = p < q ? p : q;
					hi_i = p < q ? q : p;
					if (accept && !accept(accept_user, lo_i, hi_i))
						continue;
					(*candidates)++;
					s = visit(user, chunk, lo_i, hi_i);
				}
			}
			continue;
		}
		for (e = b0 + r0; e < b0 + r1 && s == QAWS_STATUS_OK; e++)
			for (f = e + 1; f < GR->start[i + 1] && s == QAWS_STATUS_OK; f++)
			{
				unsigned int p = GR->items[e], q = GR->items[f], lo_i = p < q ? p : q, hi_i = p < q ? q : p;
				qaws_bp_box const* A = &boxes[p];
				qaws_bp_box const* B = &boxes[q];
				int ua = bp_unbounded(A, dim), ub = bp_unbounded(B, dim), over = 1;
				unsigned int cell[3];
				for (k = 0; k < 3; k++)
				{
					double c;
					if (k >= dim)
					{
						cell[k] = 0;
						continue;
					}
					if (!ua && !ub)
					{
						if (!(A->lo[k] <= B->hi[k] && B->lo[k] <= A->hi[k]))
						{
							over = 0;
							break;
						}
						c = A->lo[k] > B->lo[k] ? A->lo[k] : B->lo[k];
					}
					else if (ua && ub)
						c = GR->G.lo[k];
					else
						c = ua ? B->lo[k] : A->lo[k];
					cell[k] = bp_cell(&GR->G, c, k);
				}
				if (!over || (cell[2] * GR->G.g[1] + cell[1]) * GR->G.g[0] + cell[0] != i)
					continue;
				if (accept && !accept(accept_user, lo_i, hi_i))
					continue;
				(*candidates)++;
				s = visit(user, chunk, lo_i, hi_i);
			}
	}
	if (local != stack_box)
		qaws_internal_free(local);
	if (lcell != stack_cell)
		qaws_internal_free(lcell);
	return s;
}

/* ------------------------------------------------------------------ */
/*  Sweep and prune: better than the grid when boxes span many cells   */
/* ------------------------------------------------------------------ */

typedef struct bp_key
{
	double lo;
	unsigned int i;
} bp_key;

static int bp_cmp_key(void const* a, void const* b)
{
	double x = ((bp_key const*)a)->lo, y = ((bp_key const*)b)->lo;
	return x < y ? -1 : (x > y ? 1 : 0);
}

/* the axis with the largest extent of the bounded boxes */
static unsigned int bp_axis(qaws_bp_box const* boxes, unsigned int count, unsigned int dim)
{
	double lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, best = -1;
	unsigned int i, k, axis = 0;
	for (i = 0; i < count; i++)
		if (!bp_unbounded(&boxes[i], dim))
			for (k = 0; k < dim; k++)
			{
				if (boxes[i].lo[k] < lo[k]) lo[k] = boxes[i].lo[k];
				if (boxes[i].hi[k] > hi[k]) hi[k] = boxes[i].hi[k];
			}
	for (k = 0; k < dim; k++)
		if (hi[k] - lo[k] > best)
		{
			best = hi[k] - lo[k];
			axis = k;
		}
	return axis;
}

/* bounded boxes sorted by their low end on the axis (unbounded ones counted apart) */
static bp_key* bp_sorted(qaws_bp_box const* boxes, unsigned int count, unsigned int dim, unsigned int axis, unsigned int* out_n)
{
	bp_key* keys = (bp_key*)qaws_internal_malloc(qaws_internal_mul2(sizeof(bp_key), (count ? count : 1)));
	unsigned int i, n = 0;
	if (!keys)
		return NULL;
	for (i = 0; i < count; i++)
		if (!bp_unbounded(&boxes[i], dim))
		{
			keys[n].lo = boxes[i].lo[axis];
			keys[n].i = i;
			n++;
		}
	qsort(keys, n, sizeof(bp_key), bp_cmp_key);
	*out_n = n;
	return keys;
}

/* pairs overlapping on the axis (what the sweep tests), by binary search */
static double bp_sweep_cost(qaws_bp_box const* boxes, bp_key const* keys, unsigned int n, unsigned int axis)
{
	double cost = 0;
	unsigned int a;
	for (a = 0; a < n; a++)
	{
		double h = boxes[keys[a].i].hi[axis];
		unsigned int lo = a + 1, hi = n;
		while (lo < hi)
		{
			unsigned int m = (lo + hi) / 2;
			if (keys[m].lo <= h) lo = m + 1; else hi = m;
		}
		cost += lo - a - 1;
	}
	return cost;
}

/* pair tests the grid would make: every cell's pairs */
static double bp_grid_cost(struct qaws_bp_grid const* GR)
{
	double cost = 0;
	unsigned int c;
	for (c = 0; c < GR->ncell; c++)
	{
		double m = GR->start[c + 1] - GR->start[c];
		cost += 0.5 * m * (m - 1);
	}
	return cost;
}

static qaws_status bp_sweep(qaws_bp_box const* boxes, unsigned int count, unsigned int dim, unsigned int axis, bp_key const* keys, unsigned int n,
	qaws_bp_accept accept, void* accept_user, qaws_status (*visit)(void* user, unsigned int chunk, unsigned int i, unsigned int j), void* user,
	unsigned int* candidates)
{
	unsigned int a, b, k;
	qaws_status s = QAWS_STATUS_OK;
	for (a = 0; a < n && s == QAWS_STATUS_OK; a++)
	{
		qaws_bp_box const* A = &boxes[keys[a].i];
		double h = A->hi[axis];
		for (b = a + 1; b < n && keys[b].lo <= h && s == QAWS_STATUS_OK; b++)
		{
			qaws_bp_box const* B = &boxes[keys[b].i];
			unsigned int p = keys[a].i, q = keys[b].i;
			int over = 1;
			for (k = 0; k < dim; k++)
				if (k != axis && !(A->lo[k] <= B->hi[k] && B->lo[k] <= A->hi[k]))
				{
					over = 0;
					break;
				}
			if (!over)
				continue;
			if (p > q) { unsigned int t = p; p = q; q = t; }
			if (accept && !accept(accept_user, p, q))
				continue;
			(*candidates)++;
			s = visit(user, 0, p, q);
		}
	}
	/* unbounded boxes meet everything */
	for (a = 0; a < count && s == QAWS_STATUS_OK; a++)
		if (bp_unbounded(&boxes[a], dim))
			for (b = 0; b < count && s == QAWS_STATUS_OK; b++)
			{
				unsigned int p = a < b ? a : b, q = a < b ? b : a;
				if (b == a || (bp_unbounded(&boxes[b], dim) && b < a))
					continue;
				if (accept && !accept(accept_user, p, q))
					continue;
				(*candidates)++;
				s = visit(user, 0, p, q);
			}
	return s;
}

typedef struct bp_serial
{
	qaws_bp_visit visit;
	void* user;
} bp_serial;

static qaws_status bp_serial_visit(void* user, unsigned int chunk, unsigned int i, unsigned int j)
{
	bp_serial const* b = (bp_serial const*)user;
	(void)chunk;
	return b->visit(b->user, i, j);
}


qaws_status qaws_internal_broadphase(
	qaws_bp_stats* out_stats,
	qaws_bp_box const* boxes,
	unsigned int count,
	unsigned int dim,
	qaws_bp_accept accept,
	qaws_bp_visit visit,
	void* user)
{
	struct qaws_bp_grid GR;
	bp_serial b;
	unsigned int candidates = 0;
	qaws_status s;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (count < 2)
		return QAWS_STATUS_OK;
	s = bp_build(boxes, count, dim, &GR, 0, NULL);
	if (s != QAWS_STATUS_OK)
		return s;
	if (out_stats)
		out_stats->cell_count = GR.ncell;
	b.visit = visit;
	b.user = user;
	{
		/* the sweep when it tests fewer pairs than the grid (boxes spanning many cells) */
		unsigned int axis = bp_axis(boxes, count, dim), n = 0;
		bp_key* keys = bp_sorted(boxes, count, dim, axis, &n);
		if (keys && bp_sweep_cost(boxes, keys, n, axis) < bp_grid_cost(&GR))
			s = bp_sweep(boxes, count, dim, axis, keys, n, accept, user, bp_serial_visit, &b, &candidates);
		else
			s = bp_cells(&GR, boxes, accept, user, bp_serial_visit, &b, 0, GR.chunk, GR.chunk + 2 * GR.nchunk, &candidates);
		qaws_internal_free(keys);
	}
	if (out_stats)
		out_stats->candidate_count = candidates;
	qaws_internal_free(GR.start);
	qaws_internal_free(GR.items);
	qaws_internal_free(GR.chunk);
	return s;
}

qaws_status qaws_internal_grid_create(qaws_bp_grid** out_grid, qaws_bp_box const* boxes, unsigned int count, unsigned int dim)
{
	struct qaws_bp_grid* g;
	qaws_status s;
	*out_grid = NULL;
	g = (struct qaws_bp_grid*)qaws_internal_malloc(sizeof(struct qaws_bp_grid));
	if (!g)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s = bp_build(boxes, count, dim, g, 0, NULL);
	if (s != QAWS_STATUS_OK)
	{
		qaws_internal_free(g);
		return s;
	}
	*out_grid = g;
	return QAWS_STATUS_OK;
}

qaws_status qaws_internal_grid_create_fit(qaws_bp_grid** out_grid, qaws_bp_box const* boxes, unsigned int count, unsigned int dim,
	qaws_batch_executor const* executor)
{
	struct qaws_bp_grid* g;
	qaws_status s;
	*out_grid = NULL;
	g = (struct qaws_bp_grid*)qaws_internal_malloc(sizeof(struct qaws_bp_grid));
	if (!g)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s = bp_build(boxes, count, dim, g, 1, executor);
	if (s != QAWS_STATUS_OK)
	{
		qaws_internal_free(g);
		return s;
	}
	*out_grid = g;
	return QAWS_STATUS_OK;
}

void qaws_internal_grid_destroy(qaws_bp_grid* grid)
{
	if (!grid)
		return;
	qaws_internal_free(grid->start);
	qaws_internal_free(grid->items);
	qaws_internal_free(grid->chunk);
	qaws_internal_free(grid);
}

double qaws_internal_grid_ring(qaws_bp_grid const* g, double const* p, unsigned int r, void (*visit)(void* user, unsigned int item), void* user)
{
	int c[3], lo[3], hi[3], a, b, d, inside = 0;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		c[k] = k < g->dim ? (int)bp_cell(&g->G, p[k], k) : 0;
		lo[k] = k < g->dim ? c[k] - (int)r : 0;
		hi[k] = k < g->dim ? c[k] + (int)r : 0;
		/* the ring has cells in the grid when one of its faces does */
		if (k < g->dim && ((lo[k] >= 0) || (hi[k] <= (int)g->G.g[k] - 1)))
			inside = 1;
	}
	if (r > 0 && !inside)
		return HUGE_VAL;
	for (d = lo[2] < 0 ? 0 : lo[2]; d <= hi[2] && d < (int)g->G.g[2]; d++)
		for (b = lo[1] < 0 ? 0 : lo[1]; b <= hi[1] && b < (int)g->G.g[1]; b++)
			for (a = lo[0] < 0 ? 0 : lo[0]; a <= hi[0] && a < (int)g->G.g[0]; a++)
			{
				unsigned int cell, e;
				/* only the shell: Chebyshev distance r from c */
				int m = abs(a - c[0]);
				if (abs(b - c[1]) > m) m = abs(b - c[1]);
				if (abs(d - c[2]) > m) m = abs(d - c[2]);
				if (m != (int)r)
					continue;
				cell = ((unsigned int)d * g->G.g[1] + (unsigned int)b) * g->G.g[0] + (unsigned int)a;
				for (e = g->start[cell]; e < g->start[cell + 1]; e++)
					visit(user, g->items[e]);
			}
	/* every point of a cell r rings away is at least (r - 1) cells from p */
	return r > 1 ? (r - 1) * g->hmin : 0;
}

int qaws_internal_grid_ray(qaws_bp_grid const* g, double const* o, double const* d, double tmax,
	int (*visit)(void* user, unsigned int const* items, unsigned int count, double t_exit), void* user)
{
	double t0 = 0, t1 = tmax, tnext[3], tstep[3], hi;
	int cell[3], step[3], n[3];
	unsigned int k;
	/* clip to the grid's box (axes past dim have one cell) */
	for (k = 0; k < 3; k++)
	{
		n[k] = (int)g->G.g[k];
		if (k >= g->dim || g->G.inv[k] <= 0)
			continue;
		hi = g->G.lo[k] + n[k] / g->G.inv[k];
		if (d[k] == 0)
		{
			if (o[k] < g->G.lo[k] || o[k] > hi)
				return 0;
			continue;
		}
		{
			double ta = (g->G.lo[k] - o[k]) / d[k], tb = (hi - o[k]) / d[k];
			if (ta > tb) { double s = ta; ta = tb; tb = s; }
			if (ta > t0) t0 = ta;
			if (tb < t1) t1 = tb;
		}
	}
	if (!(t0 <= t1))
		return 0;
	/* the entry cell and the steps */
	for (k = 0; k < 3; k++)
	{
		if (k >= g->dim || g->G.inv[k] <= 0)
		{
			cell[k] = 0;
			step[k] = 0;
			tnext[k] = HUGE_VAL;
			tstep[k] = HUGE_VAL;
			continue;
		}
		cell[k] = (int)bp_cell(&g->G, o[k] + t0 * d[k], k);
		if (d[k] > 0)
		{
			step[k] = 1;
			tnext[k] = (g->G.lo[k] + (cell[k] + 1) / g->G.inv[k] - o[k]) / d[k];
			tstep[k] = 1 / (g->G.inv[k] * d[k]);
		}
		else if (d[k] < 0)
		{
			step[k] = -1;
			tnext[k] = (g->G.lo[k] + cell[k] / g->G.inv[k] - o[k]) / d[k];
			tstep[k] = -1 / (g->G.inv[k] * d[k]);
		}
		else
		{
			step[k] = 0;
			tnext[k] = HUGE_VAL;
			tstep[k] = HUGE_VAL;
		}
	}
	for (;;)
	{
		unsigned int c = ((unsigned int)cell[2] * g->G.g[1] + (unsigned int)cell[1]) * g->G.g[0] + (unsigned int)cell[0];
		double texit = tnext[0];
		unsigned int a = 0;
		for (k = 1; k < 3; k++)
			if (tnext[k] < texit)
			{
				texit = tnext[k];
				a = k;
			}
		if (texit > t1)
			texit = t1;
		if (g->start[c + 1] > g->start[c] && visit(user, g->items + g->start[c], g->start[c + 1] - g->start[c], texit))
			return 1;
		if (texit >= t1)
			return 0;
		cell[a] += step[a];
		if (cell[a] < 0 || cell[a] >= n[a])
			return 0;
		tnext[a] += tstep[a];
	}
}

unsigned int qaws_internal_grid_pair_chunks(qaws_bp_grid const* grid)
{
	return grid->nchunk;
}

typedef struct bp_par
{
	struct qaws_bp_grid const* GR;
	qaws_bp_box const* boxes;
	qaws_bp_accept accept;
	void* accept_user;
	qaws_status (*visit)(void* user, unsigned int chunk, unsigned int i, unsigned int j);
	void* visit_user;
	unsigned int* candidates;    /* one per chunk */
} bp_par;

/* one chunk of the grid's (one item each) */
static qaws_status bp_par_job(void* ctx, unsigned int chunk, unsigned int begin, unsigned int end)
{
	bp_par const* P = (bp_par const*)ctx;
	unsigned int n = 0;   /* (counted here: the chunks' counters share cache lines) */
	qaws_status s = bp_cells(P->GR, P->boxes, P->accept, P->accept_user, P->visit, P->visit_user, chunk, P->GR->chunk + 2 * (size_t)begin,
		P->GR->chunk + 2 * (size_t)end, &n);
	P->candidates[chunk] = n;
	return s;
}

qaws_status qaws_internal_grid_pairs(unsigned int* out_candidates, qaws_bp_grid const* grid, qaws_bp_box const* boxes, qaws_bp_accept accept,
	void* accept_user, qaws_status (*visit)(void* user, unsigned int chunk, unsigned int i, unsigned int j), void* visit_user,
	qaws_batch_executor const* executor)
{
	bp_par P;
	unsigned int k, n = qaws_internal_grid_pair_chunks(grid);
	qaws_status s;
	if (out_candidates)
		*out_candidates = 0;
	if (n == 0)
		return QAWS_STATUS_OK;
	P.GR = grid;
	P.boxes = boxes;
	P.accept = accept;
	P.accept_user = accept_user;
	P.visit = visit;
	P.visit_user = visit_user;
	P.candidates = (unsigned int*)qaws_internal_zalloc(n, sizeof(unsigned int));
	if (!P.candidates)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	{
		/* the sweep (serial, chunk 0) when it tests fewer pairs than the grid */
		unsigned int axis = grid->fit ? 0 : bp_axis(boxes, grid->count, grid->dim), nk = 0;
		bp_key* keys = grid->fit ? NULL : bp_sorted(boxes, grid->count, grid->dim, axis, &nk);
		if (keys && bp_sweep_cost(boxes, keys, nk, axis) < bp_grid_cost(grid))
			s = bp_sweep(boxes, grid->count, grid->dim, axis, keys, nk, accept, accept_user, visit, visit_user, &P.candidates[0]);
		else
			s = qaws_internal_parallel(n, 1, bp_par_job, &P, executor);
		qaws_internal_free(keys);
	}
	for (k = 0; k < n && out_candidates; k++)
		*out_candidates += P.candidates[k];
	qaws_internal_free(P.candidates);
	return s;
}

unsigned int qaws_internal_grid_cells(qaws_bp_grid const* grid)
{
	return grid->ncell;
}
