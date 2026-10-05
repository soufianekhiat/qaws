#ifndef QAWS_BSPLINE_DIFF_CORE_H
#define QAWS_BSPLINE_DIFF_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"

/*
 * B-spline derivative kernels for every backend (C, HLSL, GLSL, Halide).
 *
 * Same conventions as qaws_bspline_eval_3d: local_cp holds the degree + 1
 * control points of the span, local_knots the knot window around it and
 * span the span index inside that window. Jets are qaws_eval_3d
 * (position, d1, d2).
 *
 *   qaws_bspline_tangent_3d    tangent jet along control point tangents
 *                              and a parameter tangent
 *   qaws_bspline_adjoint_cp_3d contribution of one sample's jet adjoint to
 *                              local control point j: the per-element term
 *                              of a gather kernel (one thread per control
 *                              point, looping over the samples it supports;
 *                              see qaws_curve_build_support_index)
 *   qaws_bspline_adjoint_t_3d  parameter adjoint of one sample
 *
 * They reproduce qaws_curve_eval_batch_tangent_3d / _adjoint_3d of the C
 * runtime for B-spline curves.
 */

/* Basis rows 0..3 with rows above the degree set to zero. */
QAWS_INLINE void qaws_bspline_diff_rows(
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    QAWS_OUT qaws_scalar out_ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS])
{
    int stride = degree + 1;
    int max_k = degree < QAWS_CORE_MAX_DERIV ? degree : QAWS_CORE_MAX_DERIV;
    int k, j;
    qaws_bspline_basis_derivs(local_knots, degree, span, t, max_k, out_ders);
    for (k = max_k + 1; k <= QAWS_CORE_MAX_DERIV; k++)
        for (j = 0; j < stride; j++)
            out_ders[k * stride + j] = QAWS_ZERO;
}

/* sum_j N_j^(k) values[j] */
QAWS_INLINE qaws_vec3 qaws_bspline_diff_row_sum(
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS],
    qaws_scalar values[QAWS_CORE_MAX_POINTS * 3],
    int degree,
    int k)
{
    qaws_vec3 r;
    int j;
    qaws_scalar n;
    r.x = QAWS_ZERO; r.y = QAWS_ZERO; r.z = QAWS_ZERO;
    for (j = 0; j <= degree; j++) {
        n = ders[k * (degree + 1) + j];
        r.x += n * values[j * 3 + 0];
        r.y += n * values[j * 3 + 1];
        r.z += n * values[j * 3 + 2];
    }
    return r;
}

/* d C^(k) = sum N^(k) P_dot + t_dot sum N^(k+1) P,  k = 0, 1, 2 */
QAWS_INLINE qaws_eval_3d qaws_bspline_tangent_3d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    qaws_scalar t_dot)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_eval_3d r;
    qaws_vec3 a, b;

    qaws_bspline_diff_rows(local_knots, degree, span, t, ders);

    a = qaws_bspline_diff_row_sum(ders, local_cp_dot, degree, 0);
    b = qaws_bspline_diff_row_sum(ders, local_cp, degree, 1);
    r.position.x = a.x + t_dot * b.x;
    r.position.y = a.y + t_dot * b.y;
    r.position.z = a.z + t_dot * b.z;

    a = qaws_bspline_diff_row_sum(ders, local_cp_dot, degree, 1);
    b = qaws_bspline_diff_row_sum(ders, local_cp, degree, 2);
    r.d1.x = a.x + t_dot * b.x;
    r.d1.y = a.y + t_dot * b.y;
    r.d1.z = a.z + t_dot * b.z;

    a = qaws_bspline_diff_row_sum(ders, local_cp_dot, degree, 2);
    b = qaws_bspline_diff_row_sum(ders, local_cp, degree, 3);
    r.d2.x = a.x + t_dot * b.x;
    r.d2.y = a.y + t_dot * b.y;
    r.d2.z = a.z + t_dot * b.z;
    return r;
}

/* sum_k N_j^(k)(t) ybar_k for local control point j (0..degree) */
QAWS_INLINE qaws_vec3 qaws_bspline_adjoint_cp_3d(
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    int j,
    qaws_eval_3d ybar)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar n0, n1, n2;
    qaws_vec3 g;
    int stride = degree + 1;

    qaws_bspline_diff_rows(local_knots, degree, span, t, ders);
    n0 = ders[0 * stride + j];
    n1 = ders[1 * stride + j];
    n2 = ders[2 * stride + j];
    g.x = n0 * ybar.position.x + n1 * ybar.d1.x + n2 * ybar.d2.x;
    g.y = n0 * ybar.position.y + n1 * ybar.d1.y + n2 * ybar.d2.y;
    g.z = n0 * ybar.position.z + n1 * ybar.d1.z + n2 * ybar.d2.z;
    return g;
}

/* sum_k ybar_k . C^(k+1)(t) */
QAWS_INLINE qaws_scalar qaws_bspline_adjoint_t_3d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    qaws_eval_3d ybar)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_vec3 c1, c2, c3;

    qaws_bspline_diff_rows(local_knots, degree, span, t, ders);
    c1 = qaws_bspline_diff_row_sum(ders, local_cp, degree, 1);
    c2 = qaws_bspline_diff_row_sum(ders, local_cp, degree, 2);
    c3 = qaws_bspline_diff_row_sum(ders, local_cp, degree, 3);
    return ybar.position.x * c1.x + ybar.position.y * c1.y + ybar.position.z * c1.z
         + ybar.d1.x * c2.x + ybar.d1.y * c2.y + ybar.d1.z * c2.z
         + ybar.d2.x * c3.x + ybar.d2.y * c3.y + ybar.d2.z * c3.z;
}

#endif /* QAWS_BSPLINE_DIFF_CORE_H */
