#ifndef QAWS_RATIONAL_BEZIER_CORE_H
#define QAWS_RATIONAL_BEZIER_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_decasteljau_core.h"

/* ===================================================================
 * Rational Bezier evaluation (2D)
 *
 * weighted_cp: homogeneous control points (w*x, w*y, w) per point,
 *              i.e. (degree+1) * 3 scalars.
 *
 * Uses De Casteljau on homogeneous coords + quotient rule.
 * =================================================================== */

QAWS_INLINE qaws_eval_2d qaws_rational_bezier_eval_2d(
    qaws_scalar weighted_cp[QAWS_CORE_MAX_POINTS * 3],
    int degree,
    qaws_scalar t,
    int eval_flags)
{
    qaws_eval_2d result;
    qaws_scalar d1p[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar d2p[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar w0, inv_w0;
    qaws_scalar Cx, Cy, C1x, C1y, C2x, C2y;
    qaws_vec3 h_pos, h_d1, h_d2;

    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO;

    /* Evaluate homogeneous curve (3 components: wx, wy, w) */
    h_pos = qaws_decasteljau_3d(weighted_cp, degree, t);
    w0 = h_pos.z;
    inv_w0 = QAWS_ONE / w0;
    Cx = h_pos.x * inv_w0;
    Cy = h_pos.y * inv_w0;

    result.position.x = Cx;
    result.position.y = Cy;

    C1x = QAWS_ZERO; C1y = QAWS_ZERO;
    C2x = QAWS_ZERO; C2y = QAWS_ZERO;

    h_d1.x = QAWS_ZERO; h_d1.y = QAWS_ZERO; h_d1.z = QAWS_ZERO;
    h_d2.x = QAWS_ZERO; h_d2.y = QAWS_ZERO; h_d2.z = QAWS_ZERO;

    /* D1 (also for D2: the quotient rule reads it) */
    if ((eval_flags & 0x6) && degree >= 1) {
        qaws_bezier_deriv_points_3d(weighted_cp, degree, d1p);
        h_d1 = qaws_decasteljau_3d(d1p, degree - 1, t);
        C1x = (h_d1.x - Cx * h_d1.z) * inv_w0;
        C1y = (h_d1.y - Cy * h_d1.z) * inv_w0;
        if (eval_flags & 0x2) {
            result.d1.x = C1x;
            result.d1.y = C1y;
        }
    }

    /* D2 (the homogeneous second derivative vanishes below degree 2, the quotient's does not) */
    if ((eval_flags & 0x4) && degree >= 1) {
        if (degree >= 2) {
            qaws_bezier_deriv_points_3d(d1p, degree - 1, d2p);
            h_d2 = qaws_decasteljau_3d(d2p, degree - 2, t);
        }
        C2x = (h_d2.x - QAWS_LITERAL(2.0) * C1x * h_d1.z - Cx * h_d2.z) * inv_w0;
        C2y = (h_d2.y - QAWS_LITERAL(2.0) * C1y * h_d1.z - Cy * h_d2.z) * inv_w0;
        result.d2.x = C2x;
        result.d2.y = C2y;
    }

    return result;
}

/* ===================================================================
 * Rational Bezier evaluation (3D)
 *
 * weighted_cp: homogeneous control points (w*x, w*y, w*z, w),
 *              i.e. (degree+1) * 4 scalars.
 *
 * We use a flat-array De Casteljau with hdim=4 here since
 * qaws_decasteljau_3d only handles 3 components.
 * =================================================================== */

QAWS_INLINE qaws_eval_3d qaws_rational_bezier_eval_3d(
    qaws_scalar weighted_cp[QAWS_CORE_MAX_POINTS * 4],
    int degree,
    qaws_scalar t,
    int eval_flags)
{
    qaws_eval_3d result;
    /* For 3D rational, hdim=4, we do manual De Casteljau on 4 components */
    qaws_scalar work[QAWS_CORE_MAX_POINTS * 4];
    qaws_scalar d1p[QAWS_CORE_MAX_POINTS * 4];
    qaws_scalar d2p[QAWS_CORE_MAX_POINTS * 4];
    qaws_scalar hbuf[4], hd1[4], hd2[4];
    qaws_scalar w0, inv_w0;
    qaws_scalar Cx, Cy, Cz, C1x, C1y, C1z, C2x, C2y, C2z;
    qaws_scalar u = QAWS_ONE - t;
    qaws_scalar deg;
    int r, i, c;

    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO; result.position.z = QAWS_ZERO;
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO; result.d1.z = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO; result.d2.z = QAWS_ZERO;

    /* De Casteljau on 4-component homogeneous points */
    for (i = 0; i <= degree; i++)
        for (c = 0; c < 4; c++)
            work[i * 4 + c] = weighted_cp[i * 4 + c];

    for (r = 1; r <= degree; r++)
        for (i = 0; i <= degree - r; i++)
            for (c = 0; c < 4; c++)
                work[i * 4 + c] = u * work[i * 4 + c] + t * work[(i + 1) * 4 + c];

    for (c = 0; c < 4; c++) hbuf[c] = work[c];

    w0 = hbuf[3];
    inv_w0 = QAWS_ONE / w0;
    Cx = hbuf[0] * inv_w0;
    Cy = hbuf[1] * inv_w0;
    Cz = hbuf[2] * inv_w0;

    result.position.x = Cx;
    result.position.y = Cy;
    result.position.z = Cz;

    C1x = QAWS_ZERO; C1y = QAWS_ZERO; C1z = QAWS_ZERO;
    C2x = QAWS_ZERO; C2y = QAWS_ZERO; C2z = QAWS_ZERO;

    for (c = 0; c < 4; c++) { hd1[c] = QAWS_ZERO; hd2[c] = QAWS_ZERO; }

    /* D1 (also for D2: the quotient rule reads it) */
    if ((eval_flags & 0x6) && degree >= 1) {
        deg = (qaws_scalar)degree;
        for (i = 0; i < degree; i++)
            for (c = 0; c < 4; c++)
                d1p[i * 4 + c] = deg * (weighted_cp[(i + 1) * 4 + c] - weighted_cp[i * 4 + c]);

        /* De Casteljau on d1p */
        for (i = 0; i <= degree - 1; i++)
            for (c = 0; c < 4; c++)
                work[i * 4 + c] = d1p[i * 4 + c];
        for (r = 1; r <= degree - 1; r++)
            for (i = 0; i <= degree - 1 - r; i++)
                for (c = 0; c < 4; c++)
                    work[i * 4 + c] = u * work[i * 4 + c] + t * work[(i + 1) * 4 + c];
        for (c = 0; c < 4; c++) hd1[c] = work[c];

        C1x = (hd1[0] - Cx * hd1[3]) * inv_w0;
        C1y = (hd1[1] - Cy * hd1[3]) * inv_w0;
        C1z = (hd1[2] - Cz * hd1[3]) * inv_w0;
        if (eval_flags & 0x2) {
            result.d1.x = C1x;
            result.d1.y = C1y;
            result.d1.z = C1z;
        }
    }

    /* D2 (the homogeneous second derivative vanishes below degree 2, the quotient's does not) */
    if ((eval_flags & 0x4) && degree >= 1) {
        if (degree >= 2) {
            deg = (qaws_scalar)(degree - 1);
            for (i = 0; i < degree - 1; i++)
                for (c = 0; c < 4; c++)
                    d2p[i * 4 + c] = deg * (d1p[(i + 1) * 4 + c] - d1p[i * 4 + c]);

            for (i = 0; i <= degree - 2; i++)
                for (c = 0; c < 4; c++)
                    work[i * 4 + c] = d2p[i * 4 + c];
            for (r = 1; r <= degree - 2; r++)
                for (i = 0; i <= degree - 2 - r; i++)
                    for (c = 0; c < 4; c++)
                        work[i * 4 + c] = u * work[i * 4 + c] + t * work[(i + 1) * 4 + c];
            for (c = 0; c < 4; c++) hd2[c] = work[c];
        }

        C2x = (hd2[0] - QAWS_LITERAL(2.0) * C1x * hd1[3] - Cx * hd2[3]) * inv_w0;
        C2y = (hd2[1] - QAWS_LITERAL(2.0) * C1y * hd1[3] - Cy * hd2[3]) * inv_w0;
        C2z = (hd2[2] - QAWS_LITERAL(2.0) * C1z * hd1[3] - Cz * hd2[3]) * inv_w0;
        result.d2.x = C2x;
        result.d2.y = C2y;
        result.d2.z = C2z;
    }

    return result;
}

#endif /* QAWS_RATIONAL_BEZIER_CORE_H */
