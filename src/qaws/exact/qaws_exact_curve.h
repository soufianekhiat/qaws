#ifndef QAWS_EXACT_CURVE_H
#define QAWS_EXACT_CURVE_H

/* Internal: exact curve layout and the exact rational evaluation. */

#include "../qaws_exact.h"
#include "qaws_exact_int.h"

#define QAWS_EXACT_MAX_DEGREE 16

/* One span: an integer homogeneous Bezier over the parameter lattice
   interval [a, b] (t = T 2^-param_shift). */
typedef struct qaws_exact_span
{
	unsigned int degree;
	int64_t a, b;
	qaws_exact_int* h;    /* (degree + 1) points of (w x, w y, [w z], w) */
} qaws_exact_span;

struct qaws_exact_curve
{
	int dimension;                 /* 2 or 3 */
	unsigned int span_count;
	qaws_exact_span* spans;
	int param_shift;               /* t = T 2^-param_shift */
	int space_exp2;
	qaws_exact_desc desc;
};

/*
 * C^(k)(T 2^-param_shift) exactly: num[c] / den per component c (den > 0 is
 * shared), k <= 3, T clamped to the domain. The span holding T is the one
 * with a <= T < b, the last one at the end of the domain. Values are in
 * lattice units (multiply by 2^space_exp2 for world units); derivatives
 * are with respect to the public parameter t.
 */
qaws_status qaws_exact_curve_eval_rational(qaws_exact_curve const* curve, int64_t T, unsigned int k, qaws_exact_int* num,
	qaws_exact_int* den);

/* ------------------------------------------------------------------
 * Shared with the exact surfaces
 * ------------------------------------------------------------------ */

#define QAWS_EXACT_MAX_KNOTS 256

/* A homogeneous point as fractions num[0..D) / den, den > 0. */
typedef struct qaws_exact_hfrac
{
	qaws_exact_int num[4];
	qaws_exact_int den;
} qaws_exact_hfrac;

/* x -> nearest lattice integer i (x ~ i 2^exp2, ties to even), |i| < 2^bits. */
qaws_status qaws_exact_quantize(double x, int exp2, unsigned int bits, int64_t* out, double* err);

/* Parameter shift of a domain whose largest magnitude is tmax. */
int qaws_exact_param_shift_for(double tmax, unsigned int param_bits);

/* Bezier control point j of the degree-p B-spline span s (knots K, local
   points local[0..p] = P_(s-p)..P_s), exactly. */
qaws_status qaws_exact_blossom(qaws_exact_hfrac const* local, int64_t const* K, unsigned int s, unsigned int p, unsigned int j,
	unsigned int D, qaws_exact_hfrac* out);

/* b[0..p] over their lcm into integers h[(p + 1) D], divided by the common gcd. */
qaws_status qaws_exact_clear_denominators(qaws_exact_hfrac* b, unsigned int p, unsigned int D, qaws_exact_int* h);

/* Certified winding number of a closed loop of 2D exact curves around the rational point (P0 / P2, P1 / P2)
   in lattice units (P2 > 0); the loop is assumed already validated. */
qaws_status qaws_exact_winding_2d_hom(qaws_exact_curve const* const* pieces, unsigned int count, qaws_exact_int const* P, int* out_winding);

/* t = (a 2^depth + (b - a) index) 2^-(depth + shift): the nearest double and whether it is exact. */
qaws_status qaws_exact_span_param_to_double(qaws_exact_span const* sp, int shift, uint64_t index, int depth, double* out, int* exact);

/* d^j H / ds^j of the span at s = x / (b - a), scaled by (b - a)^degree. */
qaws_status qaws_exact_homogeneous_derivative(qaws_exact_span const* sp, unsigned int D, int64_t x, unsigned int j,
	qaws_exact_int* out);

/* Box of a span in lattice units, from its control points: sound (two ulps
   outward) when every weight is positive; otherwise unbounded (lo = +inf,
   hi = -inf) on every axis. z is 0 for 2D. */
void qaws_exact_span_box(qaws_exact_span const* sp, unsigned int dim, double lo[3], double hi[3]);

/* Certified intersections of span ia of a with span ib of b, appended to out
   (curve a's parameter first); a point exactly at a knot is reported only by
   the span it ends. *common is set when the spans failed on a common
   component. Returns BUFFER_TOO_SMALL past capacity. */
qaws_status qaws_exact_span_pair_hits(qaws_exact_curve const* a, unsigned int ia, qaws_exact_curve const* b, unsigned int ib,
	qaws_exact_pair* out, unsigned int capacity, unsigned int* count, int* common);

#endif /* QAWS_EXACT_CURVE_H */
