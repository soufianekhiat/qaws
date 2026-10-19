#ifndef QAWS_HORNER_CORE_H
#define QAWS_HORNER_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ===================================================================
 * Horner polynomial evaluation (2D)
 *
 * Stride-based: coefficients stored as [x0,y0, x1,y1, ...] for 2D.
 * P(t) = c0 + c1*t + c2*t^2 + ... + cn*t^n
 * =================================================================== */

QAWS_INLINE qaws_vec2 qaws_horner_2d(
    qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    qaws_scalar t)
{
    qaws_vec2 result;
    int i;
    result.x = coeffs[degree * 2 + 0];
    result.y = coeffs[degree * 2 + 1];
    for (i = degree - 1; i >= 0; i--) {
        result.x = result.x * t + coeffs[i * 2 + 0];
        result.y = result.y * t + coeffs[i * 2 + 1];
    }
    return result;
}

QAWS_INLINE qaws_vec3 qaws_horner_3d(
    qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 3],
    int degree,
    qaws_scalar t)
{
    qaws_vec3 result;
    int i;
    result.x = coeffs[degree * 3 + 0];
    result.y = coeffs[degree * 3 + 1];
    result.z = coeffs[degree * 3 + 2];
    for (i = degree - 1; i >= 0; i--) {
        result.x = result.x * t + coeffs[i * 3 + 0];
        result.y = result.y * t + coeffs[i * 3 + 1];
        result.z = result.z * t + coeffs[i * 3 + 2];
    }
    return result;
}

/* ===================================================================
 * Horner with derivatives (2D)
 *
 * P'(t) = sum_{i=1}^{degree} i * c_i * t^{i-1}
 * P''(t) = sum_{i=2}^{degree} i*(i-1) * c_i * t^{i-2}
 * =================================================================== */

QAWS_INLINE qaws_eval_2d qaws_horner_eval_2d(
    qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    qaws_scalar t,
    int eval_flags)
{
    qaws_eval_2d result;
    int i;

    result.position = qaws_horner_2d(coeffs, degree, t);
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO;

    if (((eval_flags & 0x2) != 0) && degree >= 1) {
        qaws_scalar d1x = QAWS_ITOF(degree) * coeffs[degree * 2 + 0];
        qaws_scalar d1y = QAWS_ITOF(degree) * coeffs[degree * 2 + 1];
        for (i = degree - 1; i >= 1; i--) {
            d1x = d1x * t + QAWS_ITOF(i) * coeffs[i * 2 + 0];
            d1y = d1y * t + QAWS_ITOF(i) * coeffs[i * 2 + 1];
        }
        result.d1.x = d1x;
        result.d1.y = d1y;
    }

    if (((eval_flags & 0x4) != 0) && degree >= 2) {
        qaws_scalar d2x = QAWS_ITOF(degree * (degree - 1)) * coeffs[degree * 2 + 0];
        qaws_scalar d2y = QAWS_ITOF(degree * (degree - 1)) * coeffs[degree * 2 + 1];
        for (i = degree - 1; i >= 2; i--) {
            d2x = d2x * t + QAWS_ITOF(i * (i - 1)) * coeffs[i * 2 + 0];
            d2y = d2y * t + QAWS_ITOF(i * (i - 1)) * coeffs[i * 2 + 1];
        }
        result.d2.x = d2x;
        result.d2.y = d2y;
    }

    return result;
}

/* ===================================================================
 * Horner with derivatives (3D)
 * =================================================================== */

QAWS_INLINE qaws_eval_3d qaws_horner_eval_3d(
    qaws_scalar coeffs[QAWS_CORE_MAX_POINTS * 3],
    int degree,
    qaws_scalar t,
    int eval_flags)
{
    qaws_eval_3d result;
    int i;

    result.position = qaws_horner_3d(coeffs, degree, t);
    result.d1.x = QAWS_ZERO; result.d1.y = QAWS_ZERO; result.d1.z = QAWS_ZERO;
    result.d2.x = QAWS_ZERO; result.d2.y = QAWS_ZERO; result.d2.z = QAWS_ZERO;

    if (((eval_flags & 0x2) != 0) && degree >= 1) {
        qaws_scalar d1x = QAWS_ITOF(degree) * coeffs[degree * 3 + 0];
        qaws_scalar d1y = QAWS_ITOF(degree) * coeffs[degree * 3 + 1];
        qaws_scalar d1z = QAWS_ITOF(degree) * coeffs[degree * 3 + 2];
        for (i = degree - 1; i >= 1; i--) {
            d1x = d1x * t + QAWS_ITOF(i) * coeffs[i * 3 + 0];
            d1y = d1y * t + QAWS_ITOF(i) * coeffs[i * 3 + 1];
            d1z = d1z * t + QAWS_ITOF(i) * coeffs[i * 3 + 2];
        }
        result.d1.x = d1x;
        result.d1.y = d1y;
        result.d1.z = d1z;
    }

    if (((eval_flags & 0x4) != 0) && degree >= 2) {
        qaws_scalar d2x = QAWS_ITOF(degree * (degree - 1)) * coeffs[degree * 3 + 0];
        qaws_scalar d2y = QAWS_ITOF(degree * (degree - 1)) * coeffs[degree * 3 + 1];
        qaws_scalar d2z = QAWS_ITOF(degree * (degree - 1)) * coeffs[degree * 3 + 2];
        for (i = degree - 1; i >= 2; i--) {
            d2x = d2x * t + QAWS_ITOF(i * (i - 1)) * coeffs[i * 3 + 0];
            d2y = d2y * t + QAWS_ITOF(i * (i - 1)) * coeffs[i * 3 + 1];
            d2z = d2z * t + QAWS_ITOF(i * (i - 1)) * coeffs[i * 3 + 2];
        }
        result.d2.x = d2x;
        result.d2.y = d2y;
        result.d2.z = d2z;
    }

    return result;
}

#endif /* QAWS_HORNER_CORE_H */
