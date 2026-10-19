#ifndef QAWS_NURBS_EVAL_CORE_H
#define QAWS_NURBS_EVAL_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"

/* ===================================================================
 * NURBS evaluation (2D)
 *
 * Weighted rational evaluation using A^(k)/W^(k) accumulation
 * and quotient rule for derivatives.
 *
 * local_cp:      (degree+1) control points [x0,y0, x1,y1, ...]
 * local_knots:   knot window
 * local_weights: (degree+1) weights
 * =================================================================== */

QAWS_INLINE qaws_eval_2d qaws_nurbs_eval_2d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_weights[QAWS_CORE_MAX_POINTS],
    int degree,
    int span,
    qaws_scalar t,
    int eval_flags)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int max_deriv = 0;
    int order = degree + 1;
    int k, j;
    qaws_scalar W_deriv[4];
    qaws_scalar Ax[4], Ay[4];
    qaws_scalar Nkj, w, Nw;
    qaws_scalar inv_w0;
    qaws_scalar Cx, Cy, C1x, C1y, C2x, C2y;
    qaws_eval_2d result;

    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO;

    if ((eval_flags & 0x2) != 0) max_deriv = 1;
    if ((eval_flags & 0x4) != 0) max_deriv = 2;
    if ((eval_flags & 0x8) != 0) max_deriv = 3;

    qaws_bspline_basis_derivs(local_knots, degree, span, t, max_deriv, ders);

    W_deriv[0] = QAWS_ZERO; W_deriv[1] = QAWS_ZERO;
    W_deriv[2] = QAWS_ZERO; W_deriv[3] = QAWS_ZERO;
    Ax[0] = QAWS_ZERO; Ax[1] = QAWS_ZERO; Ax[2] = QAWS_ZERO; Ax[3] = QAWS_ZERO;
    Ay[0] = QAWS_ZERO; Ay[1] = QAWS_ZERO; Ay[2] = QAWS_ZERO; Ay[3] = QAWS_ZERO;

    for (k = 0; k <= max_deriv; k++) {
        for (j = 0; j < order; j++) {
            Nkj = ders[k * order + j];
            w = local_weights[j];
            Nw = Nkj * w;
            W_deriv[k] += Nw;
            Ax[k] += Nw * local_cp[j * 2 + 0];
            Ay[k] += Nw * local_cp[j * 2 + 1];
        }
    }

    inv_w0 = QAWS_ONE / W_deriv[0];
    Cx = Ax[0] * inv_w0;
    Cy = Ay[0] * inv_w0;

    result.position.x = Cx;
    result.position.y = Cy;

    C1x = QAWS_ZERO; C1y = QAWS_ZERO;
    C2x = QAWS_ZERO; C2y = QAWS_ZERO;

    if (max_deriv >= 1) {
        C1x = (Ax[1] - Cx * W_deriv[1]) * inv_w0;
        C1y = (Ay[1] - Cy * W_deriv[1]) * inv_w0;
        result.d1.x = C1x;
        result.d1.y = C1y;
    }

    if (max_deriv >= 2) {
        C2x = (Ax[2] - QAWS_LITERAL(2.0) * C1x * W_deriv[1] - Cx * W_deriv[2]) * inv_w0;
        C2y = (Ay[2] - QAWS_LITERAL(2.0) * C1y * W_deriv[1] - Cy * W_deriv[2]) * inv_w0;
        result.d2.x = C2x;
        result.d2.y = C2y;
    }

    return result;
}

/* ===================================================================
 * NURBS evaluation (3D)
 * =================================================================== */

QAWS_INLINE qaws_eval_3d qaws_nurbs_eval_3d(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_weights[QAWS_CORE_MAX_POINTS],
    int degree,
    int span,
    qaws_scalar t,
    int eval_flags)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int max_deriv = 0;
    int order = degree + 1;
    int k, j;
    qaws_scalar W_deriv[4];
    qaws_scalar Ax[4], Ay[4], Az[4];
    qaws_scalar Nkj, w, Nw;
    qaws_scalar inv_w0;
    qaws_scalar Cx, Cy, Cz, C1x, C1y, C1z, C2x, C2y, C2z;
    qaws_eval_3d result;

    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO; result.position.z = QAWS_ZERO;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO; result.d1.z = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO; result.d2.z = QAWS_ZERO;

    if ((eval_flags & 0x2) != 0) max_deriv = 1;
    if ((eval_flags & 0x4) != 0) max_deriv = 2;
    if ((eval_flags & 0x8) != 0) max_deriv = 3;

    qaws_bspline_basis_derivs(local_knots, degree, span, t, max_deriv, ders);

    W_deriv[0] = QAWS_ZERO; W_deriv[1] = QAWS_ZERO;
    W_deriv[2] = QAWS_ZERO; W_deriv[3] = QAWS_ZERO;
    Ax[0] = QAWS_ZERO; Ax[1] = QAWS_ZERO; Ax[2] = QAWS_ZERO; Ax[3] = QAWS_ZERO;
    Ay[0] = QAWS_ZERO; Ay[1] = QAWS_ZERO; Ay[2] = QAWS_ZERO; Ay[3] = QAWS_ZERO;
    Az[0] = QAWS_ZERO; Az[1] = QAWS_ZERO; Az[2] = QAWS_ZERO; Az[3] = QAWS_ZERO;

    for (k = 0; k <= max_deriv; k++) {
        for (j = 0; j < order; j++) {
            Nkj = ders[k * order + j];
            w = local_weights[j];
            Nw = Nkj * w;
            W_deriv[k] += Nw;
            Ax[k] += Nw * local_cp[j * 3 + 0];
            Ay[k] += Nw * local_cp[j * 3 + 1];
            Az[k] += Nw * local_cp[j * 3 + 2];
        }
    }

    inv_w0 = QAWS_ONE / W_deriv[0];
    Cx = Ax[0] * inv_w0;
    Cy = Ay[0] * inv_w0;
    Cz = Az[0] * inv_w0;

    result.position.x = Cx;
    result.position.y = Cy;
    result.position.z = Cz;

    C1x = QAWS_ZERO; C1y = QAWS_ZERO; C1z = QAWS_ZERO;
    C2x = QAWS_ZERO; C2y = QAWS_ZERO; C2z = QAWS_ZERO;

    if (max_deriv >= 1) {
        C1x = (Ax[1] - Cx * W_deriv[1]) * inv_w0;
        C1y = (Ay[1] - Cy * W_deriv[1]) * inv_w0;
        C1z = (Az[1] - Cz * W_deriv[1]) * inv_w0;
        result.d1.x = C1x;
        result.d1.y = C1y;
        result.d1.z = C1z;
    }

    if (max_deriv >= 2) {
        C2x = (Ax[2] - QAWS_LITERAL(2.0) * C1x * W_deriv[1] - Cx * W_deriv[2]) * inv_w0;
        C2y = (Ay[2] - QAWS_LITERAL(2.0) * C1y * W_deriv[1] - Cy * W_deriv[2]) * inv_w0;
        C2z = (Az[2] - QAWS_LITERAL(2.0) * C1z * W_deriv[1] - Cz * W_deriv[2]) * inv_w0;
        result.d2.x = C2x;
        result.d2.y = C2y;
        result.d2.z = C2z;
    }

    return result;
}

#endif /* QAWS_NURBS_EVAL_CORE_H */
