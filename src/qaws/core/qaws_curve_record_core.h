#ifndef QAWS_CURVE_RECORD_CORE_H
#define QAWS_CURVE_RECORD_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"

/* ===================================================================
 * B-spline and NURBS curves read from buffers (the shader ABI)
 *
 * Each curve is a qaws_curve_record: its degree, point count, dimension
 * and offsets into knot, point and weight buffers that any number of
 * curves share. The arrays may be any length: the span search reads the
 * knot buffer (a bounded binary search, no window copied first), and only
 * the de Boor window (degree + 1 points, 2 (degree + 1) knots) is loaded.
 * The degree is at most QAWS_CORE_MAX_POINTS - 1.
 *
 * The code including this header defines first:
 *   QAWS_BUFFER_KNOT(i)    knot i
 *   QAWS_BUFFER_POINT(i)   float i of the points (packed, dimension floats each)
 *   QAWS_BUFFER_WEIGHT(i)  weight i (rational records only)
 * e.g. in HLSL, with StructuredBuffer<float> knots:
 *   #define QAWS_BUFFER_KNOT(i) knots[i]
 * The knots are those of qaws_curve_create_bspline / _nurbs (count +
 * degree + 1 values, the domain [knot[degree], knot[count]]); the results
 * are theirs (test 144; on a GPU, examples/eval_gpu_parity.c).
 * =================================================================== */

QAWS_STRUCT(qaws_curve_record) {
    int degree;
    int count;          /* control points */
    int dimension;      /* 2 or 3 */
    int rational;       /* non-zero: weights at weight_offset */
    int knot_offset;    /* count + degree + 1 knots */
    int point_offset;   /* count * dimension floats */
    int weight_offset;  /* count weights */
    int reserved;
} QAWS_STRUCT_END(qaws_curve_record);

#define QAWS_RECORD_SEARCH_STEPS 32

/* the knot span of t in [knot[degree], knot[count]]: knot[span] <= t < knot[span + 1], the last non-empty span at the end */
QAWS_INLINE int qaws_record_find_span(qaws_curve_record r, qaws_scalar t)
{
    int lo = r.degree;
    int hi = r.count;
    int it;
    if (t >= QAWS_BUFFER_KNOT(r.knot_offset + r.count)) {
        /* a knot repeated at the end leaves empty spans before it (at most degree of them) */
        lo = r.count - 1;
        for (it = 0; it < QAWS_CORE_MAX_POINTS; it++)
            if (lo > r.degree && QAWS_BUFFER_KNOT(r.knot_offset + lo) >= QAWS_BUFFER_KNOT(r.knot_offset + r.count))
                lo = lo - 1;
        return lo;
    }
    /* invariant knot[lo] <= t < knot[hi] */
    for (it = 0; it < QAWS_RECORD_SEARCH_STEPS; it++) {
        if (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (t < QAWS_BUFFER_KNOT(r.knot_offset + mid))
                hi = mid;
            else
                lo = mid;
        }
    }
    return lo;
}

/* Position, D1 and D2 (eval_flags 0x2, 0x4) at t, clamped to the domain; 2D records leave z at zero. Rational records by the
   quotient rule on the weighted sums. */
QAWS_INLINE qaws_eval_3d qaws_record_eval_3d(qaws_curve_record r, qaws_scalar t, int eval_flags)
{
    qaws_scalar lk[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar A[9];
    qaws_scalar W[3];
    qaws_scalar t0, t1, inv;
    int d = r.degree;
    int span, j, k, c, max_deriv = 0;
    qaws_eval_3d result;

    t0 = QAWS_BUFFER_KNOT(r.knot_offset + d);
    t1 = QAWS_BUFFER_KNOT(r.knot_offset + r.count);
    t = qaws_clamp(t, t0, t1);
    span = qaws_record_find_span(r, t);
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
        lk[j] = j < 2 * (d + 1) ? QAWS_BUFFER_KNOT(r.knot_offset + span - d + j) : QAWS_ZERO;
    if ((eval_flags & 0x4) != 0) max_deriv = 2;
    else if ((eval_flags & 0x2) != 0) max_deriv = 1;
    if (max_deriv > d) max_deriv = d;
    qaws_bspline_basis_derivs(lk, d, d, t, max_deriv, ders);

    for (k = 0; k < 9; k++) A[k] = QAWS_ZERO;
    for (k = 0; k < 3; k++) W[k] = QAWS_ZERO;
    for (j = 0; j <= d; j++) {
        int idx = span - d + j;
        qaws_scalar x[3];
        qaws_scalar w = QAWS_ONE;
        for (c = 0; c < 3; c++)
            x[c] = c < r.dimension ? QAWS_BUFFER_POINT(r.point_offset + idx * r.dimension + c) : QAWS_ZERO;
        if (r.rational != 0)
            w = QAWS_BUFFER_WEIGHT(r.weight_offset + idx);
        for (k = 0; k <= max_deriv; k++) {
            qaws_scalar bw = ders[k * (d + 1) + j] * w;
            W[k] += bw;
            for (c = 0; c < 3; c++)
                A[k * 3 + c] += bw * x[c];
        }
    }
    if (r.rational != 0) {
        inv = QAWS_ONE / W[0];
        for (c = 0; c < 3; c++) {
            qaws_scalar c0 = A[c] * inv;
            qaws_scalar c1 = (A[3 + c] - c0 * W[1]) * inv;
            A[6 + c] = (A[6 + c] - QAWS_LITERAL(2.0) * c1 * W[1] - c0 * W[2]) * inv;
            A[c] = c0;
            A[3 + c] = c1;
        }
    }
    result.position.x = A[0]; result.position.y = A[1]; result.position.z = A[2];
    result.d1.x = A[3]; result.d1.y = A[4]; result.d1.z = A[5];
    result.d2.x = A[6]; result.d2.y = A[7]; result.d2.z = A[8];
    return result;
}

#endif /* QAWS_CURVE_RECORD_CORE_H */
