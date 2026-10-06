#ifndef QAWS_BSPLINE_BASIS_CORE_H
#define QAWS_BSPLINE_BASIS_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

#define QAWS_CORE_MAX_DERIV 3

/* ===================================================================
 * Knot span search — available on all backends
 *
 * Finds index i such that knots[i] <= t < knots[i+1].
 * =================================================================== */

#if QAWS_BACKEND == QAWS_BACKEND_HALIDE

/*
 * Halide: branchless linear scan returning a symbolic integer expression.
 *
 * When t is a Halide Var or RDom the returned Expr evaluates the
 * correct span per sample at runtime.  O(num_cp) selects, but
 * num_cp <= QAWS_CORE_MAX_POINTS (17) so this is always bounded.
 */
QAWS_INLINE Halide::Expr qaws_find_span(
    qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int num_cp,
    qaws_scalar t)
{
    int n = num_cp - 1;
    Halide::Expr span = Halide::Expr(degree);
    int i;

    for (i = degree; i < n; i++) {
        /* spans starting at the domain end are empty */
        span = Halide::select(t >= knots[i + 1] && knots[i + 1] < knots[n + 1],
                              Halide::Expr(i + 1), span);
    }

    return span;
}

#else /* C, HLSL, GLSL */

/*
 * Bounded binary search (shader-safe).
 * 20 iterations suffice for any practical knot vector.
 */
QAWS_INLINE int qaws_find_span(
    qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int num_cp,
    qaws_scalar t)
{
    int n = num_cp - 1;
    int low, high, mid, iter, s;

    /* at the domain ends, the nearest non-empty span (a knot repeated at
       the end of an unclamped vector leaves empty spans there) */
    if (t >= knots[n + 1]) {
        s = n;
        for (iter = 0; iter < QAWS_CORE_MAX_POINTS; iter++)
            if (s > degree && knots[s] >= knots[n + 1])
                s--;
        return s;
    }
    if (t <= knots[degree]) {
        s = degree;
        for (iter = 0; iter < QAWS_CORE_MAX_POINTS; iter++)
            if (s < n && knots[s + 1] <= knots[degree])
                s++;
        return s;
    }

    low = degree;
    high = n;

    for (iter = 0; iter < 20; iter++) {
        if (low >= high)
            break;
        mid = (low + high) / 2;
        if (t < knots[mid])
            high = mid;
        else if (t >= knots[mid + 1])
            low = mid + 1;
        else
            return mid;
    }

    return low;
}

#endif /* QAWS_BACKEND == QAWS_BACKEND_HALIDE */

/* ===================================================================
 * B-spline basis functions (no derivatives)
 *
 * Matches qaws_internal_bspline_basis() but with int types
 * for shader compatibility.
 * =================================================================== */

QAWS_INLINE void qaws_bspline_basis(
    qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    QAWS_OUT qaws_scalar out_N[QAWS_CORE_MAX_POINTS])
{
    qaws_scalar left[QAWS_CORE_MAX_POINTS];
    qaws_scalar right[QAWS_CORE_MAX_POINTS];
    int j, r;
    qaws_scalar saved, denom, temp;

    out_N[0] = QAWS_ONE;

    for (j = 1; j <= degree; j++) {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        saved = QAWS_ZERO;

        for (r = 0; r < j; r++) {
            denom = right[r + 1] + left[j - r];
            temp = QAWS_SELECT(denom != QAWS_ZERO, out_N[r] / denom, QAWS_ZERO);
            out_N[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }

        out_N[j] = saved;
    }
}

/* ===================================================================
 * B-spline basis functions + derivatives (stack buffer only)
 *
 * Matches qaws_internal_bspline_basis_derivs() but uses
 * fixed-size stack buffers instead of malloc.
 *
 * out_ders: (max_k+1) x (degree+1) row-major
 * =================================================================== */

QAWS_INLINE void qaws_bspline_basis_derivs(
    qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    int max_k,
    QAWS_OUT qaws_scalar out_ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS])
{
    int p = degree;
    int order = p + 1;
    qaws_scalar ndu[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS];
    qaws_scalar left_arr[QAWS_CORE_MAX_POINTS];
    qaws_scalar right_arr[QAWS_CORE_MAX_POINTS];
    qaws_scalar a0[QAWS_CORE_MAX_POINTS];
    qaws_scalar a1[QAWS_CORE_MAX_POINTS];
    int j, r, d_idx, s1, s2, rj, pk, j1, j2;
    qaws_scalar saved, temp, dd, factor;

    ndu[0] = QAWS_ONE;

    for (j = 1; j <= p; j++) {
        left_arr[j] = t - knots[span + 1 - j];
        right_arr[j] = knots[span + j] - t;
        saved = QAWS_ZERO;

        for (r = 0; r < j; r++) {
            ndu[j * order + r] = right_arr[r + 1] + left_arr[j - r];
            temp = ndu[r * order + (j - 1)] / ndu[j * order + r];
            ndu[r * order + j] = saved + right_arr[r + 1] * temp;
            saved = left_arr[j - r] * temp;
        }

        ndu[j * order + j] = saved;
    }

    for (j = 0; j <= p; j++) {
        out_ders[0 * order + j] = ndu[j * order + p];
    }

    for (r = 0; r <= p; r++) {
        s1 = 0;
        s2 = 1;
        a0[0] = QAWS_ONE;

        for (d_idx = 1; d_idx <= max_k; d_idx++) {
            dd = QAWS_ZERO;
            rj = r - d_idx;
            pk = p - d_idx;

            if (rj >= 0) {
                if (s1 == 0)
                    a1[0] = a0[0] / ndu[(pk + 1) * order + rj];
                else
                    a1[0] = a1[0]; /* placeholder: swapped arrays handle this */
                /* Re-derive using swap logic */
            }

            j1 = (rj >= -1) ? 1 : (-rj);
            j2 = (r - 1 <= pk) ? (d_idx - 1) : (p - r);

            /* Use a[s1] for reading, a[s2] for writing.
               We use a0/a1 and swap indices. */
            if (rj >= 0) {
                if (s1 == 0)
                    a1[0] = a0[0] / ndu[(pk + 1) * order + rj];
                else
                    a0[0] = a1[0] / ndu[(pk + 1) * order + rj];
                if (s1 == 0)
                    dd = a1[0] * ndu[rj * order + pk];
                else
                    dd = a0[0] * ndu[rj * order + pk];
            }

            for (j = j1; j <= j2; j++) {
                if (s1 == 0) {
                    a1[j] = (a0[j] - a0[j - 1]) / ndu[(pk + 1) * order + (rj + j)];
                    dd += a1[j] * ndu[(rj + j) * order + pk];
                } else {
                    a0[j] = (a1[j] - a1[j - 1]) / ndu[(pk + 1) * order + (rj + j)];
                    dd += a0[j] * ndu[(rj + j) * order + pk];
                }
            }

            if (r <= pk) {
                if (s1 == 0) {
                    a1[d_idx] = -a0[d_idx - 1] / ndu[(pk + 1) * order + r];
                    dd += a1[d_idx] * ndu[r * order + pk];
                } else {
                    a0[d_idx] = -a1[d_idx - 1] / ndu[(pk + 1) * order + r];
                    dd += a0[d_idx] * ndu[r * order + pk];
                }
            }

            out_ders[d_idx * order + r] = dd;

            /* Swap s1, s2 */
            { int tmp = s1; s1 = s2; s2 = tmp; }
        }
    }

    factor = (qaws_scalar)p;
    for (d_idx = 1; d_idx <= max_k; d_idx++) {
        for (j = 0; j <= p; j++) {
            out_ders[d_idx * order + j] *= factor;
        }
        factor *= (qaws_scalar)(p - d_idx);
    }
}

#endif /* QAWS_BSPLINE_BASIS_CORE_H */
