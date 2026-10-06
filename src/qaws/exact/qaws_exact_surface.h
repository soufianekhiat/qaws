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

#endif /* QAWS_EXACT_SURFACE_H */
