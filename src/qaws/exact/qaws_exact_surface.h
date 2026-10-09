#ifndef QAWS_EXACT_SURFACE_H
#define QAWS_EXACT_SURFACE_H

/* Internal: exact surface layout and the exact rational evaluation. */

#include "qaws_exact_curve.h"

/* Patches are integer homogeneous Bezier nets over [ub[iu], ub[iu+1]] x
   [vb[iv], vb[iv+1]] on the parameter lattices u = U 2^-u_shift,
   v = V 2^-v_shift. */
struct qaws_exact_surface
{
	unsigned int p, q;                /* degrees in u and v */
	unsigned int nu, nv;              /* patch counts */
	int64_t* ub;                      /* nu + 1 breaks */
	int64_t* vb;                      /* nv + 1 breaks */
	qaws_exact_int** patch;           /* nu nv nets of (p + 1)(q + 1) (wx, wy, wz, w), u index major */
	int u_shift, v_shift;
	int space_exp2;
};

/*
 * d^(a+b) S / du^a dv^b at (U 2^-u_shift, V 2^-v_shift) exactly, a + b <= 2:
 * num[c] / den in lattice units (den > 0), the patch chosen as for curves
 * (lower end inclusive, the last patch at the domain end).
 */
qaws_status qaws_exact_surface_eval_rational(qaws_exact_surface const* s, int64_t U, int64_t V, unsigned int a, unsigned int b,
	qaws_exact_int* num, qaws_exact_int* den);

/* Su x Sv exactly: num[c] / den in squared lattice units. */
qaws_status qaws_exact_surface_normal_rational(qaws_exact_surface const* s, int64_t U, int64_t V, qaws_exact_int* num,
	qaws_exact_int* den);

/* Box of patch (iu, iv) in lattice units from its control net: sound (two
   ulps outward) when every weight is positive, otherwise unbounded. */
void qaws_exact_patch_box(qaws_exact_surface const* s, unsigned int iu, unsigned int iv, double lo[3], double hi[3]);

/* Certified intersections of span ks of a 3D curve with patch (iu, iv),
   appended to out without deduplication (a root on a span or patch edge is
   found by both neighbours). */
qaws_status qaws_exact_span_patch_hits(qaws_exact_curve const* curve, unsigned int ks, qaws_exact_surface const* surface, unsigned int iu,
	unsigned int iv, qaws_exact_curve_surface_hit* out, unsigned int capacity, unsigned int* count);

/* do two hits' enclosures overlap in t, u and v (the same root)? */
int qaws_exact_curve_surface_hits_overlap(qaws_exact_curve_surface_hit const* a, qaws_exact_curve_surface_hit const* b);

/* qaws_exact_surface_surface_hits over the listed patch pairs only
   (pairs[4 k .. 4 k + 3] = iu1, iv1, iu2, iv2); the others must be proven
   apart by the caller. */
qaws_status qaws_exact_ssi_solve(qaws_exact_surface const* a, qaws_exact_surface const* b, unsigned int min_depth, unsigned int const* pairs,
	unsigned int pair_count, qaws_exact_ssi_point* out_points, unsigned int point_capacity, unsigned int* out_point_count,
	qaws_exact_ssi_branch* out_branches, unsigned int branch_capacity, unsigned int* out_branch_count);

#endif /* QAWS_EXACT_SURFACE_H */
