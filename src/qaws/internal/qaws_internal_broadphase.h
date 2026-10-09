#ifndef QAWS_INTERNAL_BROADPHASE_H
#define QAWS_INTERNAL_BROADPHASE_H

#include "../qaws_status.h"

/*
 * Shared broad phase of the batched operations: a uniform grid over boxes
 * (about one box per cell, filled by counting sort). Every pair i < j of
 * boxes that overlap and that `accept` allows is visited exactly once, in
 * the cell holding the low corner of the two boxes' overlap.
 *
 * A box with lo > hi on some axis (e.g. lo = +inf, hi = -inf) is unbounded:
 * it overlaps everything and is placed in every cell.
 */

typedef struct qaws_bp_box
{
	double lo[3], hi[3];
} qaws_bp_box;

/* may i and j (i < j) be a pair? (NULL: every pair) */
typedef int (*qaws_bp_accept)(void* user, unsigned int i, unsigned int j);
/* one overlapping pair; anything but OK stops the traversal and is returned */
typedef qaws_status (*qaws_bp_visit)(void* user, unsigned int i, unsigned int j);

typedef struct qaws_bp_stats
{
	unsigned int cell_count;
	unsigned int candidate_count;   /* pairs visited */
} qaws_bp_stats;

/* dim 2 or 3 (dim 2 ignores the z of the boxes). out_stats may be NULL. */
qaws_status qaws_internal_broadphase(
	qaws_bp_box const* boxes,
	unsigned int count,
	unsigned int dim,
	qaws_bp_accept accept,
	qaws_bp_visit visit,
	void* user,
	qaws_bp_stats* out_stats);

/*
 * The same grid kept for point queries: qaws_internal_grid_ring visits the
 * items of every cell at Chebyshev distance r from the cell of p (an item
 * spanning several cells is visited once per cell) and returns a lower
 * bound on the distance from p to those cells, or HUGE_VAL when ring r has
 * no cell in the grid. Rings 0, 1, 2, ... search outward.
 */
typedef struct qaws_bp_grid qaws_bp_grid;

qaws_status qaws_internal_grid_create(qaws_bp_box const* boxes, unsigned int count, unsigned int dim, qaws_bp_grid** out_grid);
void qaws_internal_grid_destroy(qaws_bp_grid* grid);
double qaws_internal_grid_ring(qaws_bp_grid const* grid, double const* p, unsigned int r, void (*visit)(void* user, unsigned int item), void* user);

/*
 * The cells crossed by the ray o + t d, t in [0, tmax], in order (clipped to
 * the grid's box): visit gets each non-empty cell's items and the t where
 * the ray leaves the cell, and returns non-zero to stop. Returns 1 when
 * stopped by visit.
 */
int qaws_internal_grid_ray(qaws_bp_grid const* grid, double const* o, double const* d, double tmax,
	int (*visit)(void* user, unsigned int const* items, unsigned int count, double t_exit), void* user);

/*
 * The pair walk of qaws_internal_broadphase over a built grid, its cells in
 * chunks run through the executor (any order, any thread): visit gets the
 * chunk of the cell it runs in, one of qaws_internal_grid_pair_chunks(grid),
 * for per-chunk storage. out_candidates (may be NULL): the pairs visited.
 */
struct qaws_batch_executor;
unsigned int qaws_internal_grid_cells(qaws_bp_grid const* grid);
unsigned int qaws_internal_grid_pair_chunks(qaws_bp_grid const* grid);
qaws_status qaws_internal_grid_pairs(qaws_bp_grid const* grid, qaws_bp_box const* boxes, qaws_bp_accept accept, void* accept_user,
	qaws_status (*visit)(void* user, unsigned int chunk, unsigned int i, unsigned int j), void* visit_user, struct qaws_batch_executor const* executor,
	unsigned int* out_candidates);

#endif /* QAWS_INTERNAL_BROADPHASE_H */
