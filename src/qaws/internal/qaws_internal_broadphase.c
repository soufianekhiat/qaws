#include "qaws_internal_broadphase.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define BP_MAX_CELLS (1u << 22)
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
	unsigned int* start;        /* ncell + 1 */
	unsigned int* items;
	double hmin;                /* smallest cell side over the axes with extent */
};

static qaws_status bp_build(qaws_bp_box const* boxes, unsigned int count, unsigned int dim, struct qaws_bp_grid* B)
{
	bp_grid G;
	double hi[3], h = 0;
	unsigned int i, k, ncell, total = 0, bounded = 0;
	unsigned int *start = NULL, *fill = NULL, *items = NULL;
	memset(B, 0, sizeof(*B));
	/* scene box of the bounded boxes */
	for (k = 0; k < 3; k++)
	{
		G.lo[k] = HUGE_VAL;
		hi[k] = -HUGE_VAL;
	}
	for (i = 0; i < count; i++)
		if (!bp_unbounded(&boxes[i], dim))
		{
			bounded++;
			for (k = 0; k < dim; k++)
			{
				if (boxes[i].lo[k] < G.lo[k]) G.lo[k] = boxes[i].lo[k];
				if (boxes[i].hi[k] > hi[k]) hi[k] = boxes[i].hi[k];
			}
		}
	if (!bounded)
		for (k = 0; k < 3; k++)
		{
			G.lo[k] = 0;
			hi[k] = 0;
		}
	for (k = dim; k < 3; k++)
	{
		G.lo[k] = 0;
		hi[k] = 0;
	}

	/* cell size: about one bounded box per cell over the axes with extent */
	{
		double vol = 1, ext = 0;
		unsigned int axes = 0;
		for (k = 0; k < dim; k++)
			if (hi[k] - G.lo[k] > ext)
				ext = hi[k] - G.lo[k];
		for (k = 0; k < dim; k++)
			if (hi[k] - G.lo[k] > ext * 1e-9)
			{
				vol *= hi[k] - G.lo[k];
				axes++;
			}
		h = axes ? pow(vol / (bounded ? bounded : 1), 1.0 / axes) : 0;
	}
	for (k = 0; k < 3; k++)
	{
		double n = k < dim && h > 0 ? ceil((hi[k] - G.lo[k]) / h) : 1;
		G.g[k] = n < 1 ? 1 : (n > BP_MAX_AXIS ? BP_MAX_AXIS : (unsigned int)n);
	}
	while ((double)G.g[0] * G.g[1] * G.g[2] > BP_MAX_CELLS)
		for (k = 0; k < 3; k++)
			G.g[k] = (G.g[k] + 1) / 2;
	for (k = 0; k < 3; k++)
		G.inv[k] = hi[k] > G.lo[k] ? (double)G.g[k] / (hi[k] - G.lo[k]) : 0;
	ncell = G.g[0] * G.g[1] * G.g[2];

	/* count, prefix, fill */
	start = (unsigned int*)calloc((size_t)ncell + 1, sizeof(unsigned int));
	fill = (unsigned int*)malloc(((size_t)ncell + 1) * sizeof(unsigned int));
	if (!start || !fill)
	{
		goto done;
	}
	for (i = 0; i < count; i++)
	{
		unsigned int c0[3], c1[3], a, b, c;
		bp_range(&G, &boxes[i], dim, c0, c1);
		for (c = c0[2]; c <= c1[2]; c++)
			for (b = c0[1]; b <= c1[1]; b++)
				for (a = c0[0]; a <= c1[0]; a++)
					start[(c * G.g[1] + b) * G.g[0] + a]++;
	}
	for (i = 0; i < ncell; i++)
	{
		unsigned int n = start[i];
		start[i] = total;
		total += n;
	}
	start[ncell] = total;
	memcpy(fill, start, ((size_t)ncell + 1) * sizeof(unsigned int));
	items = (unsigned int*)malloc((size_t)(total ? total : 1) * sizeof(unsigned int));
	if (!items)
	{
		goto done;
	}
	for (i = 0; i < count; i++)
	{
		unsigned int c0[3], c1[3], a, b, c;
		bp_range(&G, &boxes[i], dim, c0, c1);
		for (c = c0[2]; c <= c1[2]; c++)
			for (b = c0[1]; b <= c1[1]; b++)
				for (a = c0[0]; a <= c1[0]; a++)
					items[fill[(c * G.g[1] + b) * G.g[0] + a]++] = i;
	}
	free(fill);
	B->G = G;
	B->dim = dim;
	B->ncell = ncell;
	B->start = start;
	B->items = items;
	B->hmin = HUGE_VAL;
	for (k = 0; k < dim; k++)
		if (G.inv[k] > 0 && 1 / G.inv[k] < B->hmin)
			B->hmin = 1 / G.inv[k];
	if (!(B->hmin < HUGE_VAL))
		B->hmin = 0;
	return QAWS_STATUS_OK;
done:
	free(start);
	free(fill);
	free(items);
	return QAWS_STATUS_ALLOCATION_FAILURE;
}

qaws_status qaws_internal_broadphase(
	qaws_bp_box const* boxes,
	unsigned int count,
	unsigned int dim,
	qaws_bp_accept accept,
	qaws_bp_visit visit,
	void* user,
	qaws_bp_stats* out_stats)
{
	struct qaws_bp_grid GR;
	unsigned int i, k;
	qaws_status s;
	if (out_stats)
		memset(out_stats, 0, sizeof(*out_stats));
	if (count < 2)
		return QAWS_STATUS_OK;
	s = bp_build(boxes, count, dim, &GR);
	if (s != QAWS_STATUS_OK)
		return s;
	if (out_stats)
		out_stats->cell_count = GR.ncell;
	/* pairs, each in the cell holding the low corner of the overlap */
	for (i = 0; i < GR.ncell && s == QAWS_STATUS_OK; i++)
	{
		unsigned int e, f;
		for (e = GR.start[i]; e < GR.start[i + 1] && s == QAWS_STATUS_OK; e++)
			for (f = e + 1; f < GR.start[i + 1] && s == QAWS_STATUS_OK; f++)
			{
				unsigned int p = GR.items[e], q = GR.items[f], lo_i = p < q ? p : q, hi_i = p < q ? q : p;
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
						c = GR.G.lo[k];
					else
						c = ua ? B->lo[k] : A->lo[k];
					cell[k] = bp_cell(&GR.G, c, k);
				}
				if (!over || (cell[2] * GR.G.g[1] + cell[1]) * GR.G.g[0] + cell[0] != i)
					continue;
				if (accept && !accept(user, lo_i, hi_i))
					continue;
				if (out_stats)
					out_stats->candidate_count++;
				s = visit(user, lo_i, hi_i);
			}
	}
	free(GR.start);
	free(GR.items);
	return s;
}

qaws_status qaws_internal_grid_create(qaws_bp_box const* boxes, unsigned int count, unsigned int dim, qaws_bp_grid** out_grid)
{
	struct qaws_bp_grid* g;
	qaws_status s;
	*out_grid = NULL;
	g = (struct qaws_bp_grid*)malloc(sizeof(struct qaws_bp_grid));
	if (!g)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	s = bp_build(boxes, count, dim, g);
	if (s != QAWS_STATUS_OK)
	{
		free(g);
		return s;
	}
	*out_grid = g;
	return QAWS_STATUS_OK;
}

void qaws_internal_grid_destroy(qaws_bp_grid* grid)
{
	if (!grid)
		return;
	free(grid->start);
	free(grid->items);
	free(grid);
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
