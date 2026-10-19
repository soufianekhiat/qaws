#ifndef QAWS_FRAME_CORE_H
#define QAWS_FRAME_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ===================================================================
 * Rotation-minimizing frames: one double-reflection step (Wang, Juttler,
 * Zheng and Liu, "Computation of rotation minimizing frames", ACM TOG
 * 2008). From the frame vector r0 (normal to the unit tangent t0) at the
 * point x0, the vector r1 at x1, unit tangent t1: reflect r0 and t0 in the
 * plane bisecting x0 x1, then in the plane bisecting the reflected tangent
 * and t1. Fourth-order accurate in the step; the frame (t, r, t x r)
 * follows the curve without twisting about it.
 * =================================================================== */

QAWS_INLINE qaws_scalar qaws_frame_dot(qaws_vec3 a, qaws_vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

QAWS_INLINE qaws_vec3 qaws_frame_reflect(qaws_vec3 v, qaws_vec3 n, qaws_scalar nn)
{
    /* v reflected in the plane through the origin normal to n (nn = n.n > 0) */
    qaws_scalar k = QAWS_LITERAL(2.0) * qaws_frame_dot(n, v) / nn;
    qaws_vec3 r;
    r.x = v.x - k * n.x;
    r.y = v.y - k * n.y;
    r.z = v.z - k * n.z;
    return r;
}

QAWS_INLINE qaws_vec3 qaws_double_reflection(qaws_vec3 x0, qaws_vec3 t0, qaws_vec3 r0, qaws_vec3 x1, qaws_vec3 t1)
{
    qaws_vec3 v1, v2, rl, tl;
    qaws_scalar c1, c2;
    v1.x = x1.x - x0.x;
    v1.y = x1.y - x0.y;
    v1.z = x1.z - x0.z;
    c1 = qaws_frame_dot(v1, v1);
    rl = r0;
    tl = t0;
    if (c1 > QAWS_ZERO) {
        rl = qaws_frame_reflect(r0, v1, c1);
        tl = qaws_frame_reflect(t0, v1, c1);
    }
    v2.x = t1.x - tl.x;
    v2.y = t1.y - tl.y;
    v2.z = t1.z - tl.z;
    c2 = qaws_frame_dot(v2, v2);
    if (c2 > QAWS_ZERO)
        rl = qaws_frame_reflect(rl, v2, c2);
    return rl;
}

#endif /* QAWS_FRAME_CORE_H */
