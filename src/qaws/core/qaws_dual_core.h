#ifndef QAWS_DUAL_CORE_H
#define QAWS_DUAL_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"

/* ===================================================================
 * Second-order forward numbers and their adjoint rules.
 *
 * A dual value carries a primal value v, its directional derivative t
 * along one tangent, and the second directional derivative tt along the
 * same tangent (not divided by 2):
 *
 *     x(e) = v + e * t + e^2/2 * tt + O(e^3)
 *
 * These kernels are the building blocks of every differential rule built
 * on top of spatial jets (normals, curvature, frames, rational projection,
 * offsets). They are backend-neutral: no pointers, no allocation.
 *
 * Adjoint kernels take primal values and an output adjoint and return the
 * input adjoint contributions (to be added by the caller).
 * =================================================================== */

QAWS_TYPE_DEF struct {
    qaws_scalar v, t, tt;
} qaws_dual1;

QAWS_TYPE_DEF struct {
    qaws_vec3 v, t, tt;
} qaws_dual3;

QAWS_TYPE_DEF struct {
    qaws_vec3 a, b;
} qaws_vec3_pair;

QAWS_TYPE_DEF struct {
    qaws_vec3 x;
    qaws_scalar s;
} qaws_vec3_scalar;

/* -------------------------------------------------------------------
 * Plain vec3 helpers
 * ------------------------------------------------------------------- */

QAWS_INLINE qaws_vec3 qaws_v3(qaws_scalar x, qaws_scalar y, qaws_scalar z)
{
    qaws_vec3 r; r.x = x; r.y = y; r.z = z; return r;
}

QAWS_INLINE qaws_vec3 qaws_v3_add(qaws_vec3 a, qaws_vec3 b)
{
    return qaws_v3(a.x + b.x, a.y + b.y, a.z + b.z);
}

QAWS_INLINE qaws_vec3 qaws_v3_sub(qaws_vec3 a, qaws_vec3 b)
{
    return qaws_v3(a.x - b.x, a.y - b.y, a.z - b.z);
}

QAWS_INLINE qaws_vec3 qaws_v3_scale(qaws_vec3 a, qaws_scalar s)
{
    return qaws_v3(a.x * s, a.y * s, a.z * s);
}

/* a + b * s */
QAWS_INLINE qaws_vec3 qaws_v3_axpy(qaws_vec3 a, qaws_vec3 b, qaws_scalar s)
{
    return qaws_v3(a.x + b.x * s, a.y + b.y * s, a.z + b.z * s);
}

QAWS_INLINE qaws_scalar qaws_v3_dot(qaws_vec3 a, qaws_vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

QAWS_INLINE qaws_vec3 qaws_v3_cross(qaws_vec3 a, qaws_vec3 b)
{
    return qaws_v3(a.y * b.z - a.z * b.y,
                   a.z * b.x - a.x * b.z,
                   a.x * b.y - a.y * b.x);
}

QAWS_INLINE qaws_vec3 qaws_v3_zero(void)
{
    return qaws_v3(QAWS_ZERO, QAWS_ZERO, QAWS_ZERO);
}

/* -------------------------------------------------------------------
 * Construction
 * ------------------------------------------------------------------- */

QAWS_INLINE qaws_dual1 qaws_dual1_make(qaws_scalar v, qaws_scalar t, qaws_scalar tt)
{
    qaws_dual1 r; r.v = v; r.t = t; r.tt = tt; return r;
}

QAWS_INLINE qaws_dual1 qaws_dual1_const(qaws_scalar v)
{
    return qaws_dual1_make(v, QAWS_ZERO, QAWS_ZERO);
}

QAWS_INLINE qaws_dual3 qaws_dual3_make(qaws_vec3 v, qaws_vec3 t, qaws_vec3 tt)
{
    qaws_dual3 r; r.v = v; r.t = t; r.tt = tt; return r;
}

QAWS_INLINE qaws_dual3 qaws_dual3_const(qaws_vec3 v)
{
    return qaws_dual3_make(v, qaws_v3_zero(), qaws_v3_zero());
}

/* -------------------------------------------------------------------
 * Scalar rules
 * ------------------------------------------------------------------- */

QAWS_INLINE qaws_dual1 qaws_dual1_add(qaws_dual1 a, qaws_dual1 b)
{
    return qaws_dual1_make(a.v + b.v, a.t + b.t, a.tt + b.tt);
}

QAWS_INLINE qaws_dual1 qaws_dual1_sub(qaws_dual1 a, qaws_dual1 b)
{
    return qaws_dual1_make(a.v - b.v, a.t - b.t, a.tt - b.tt);
}

QAWS_INLINE qaws_dual1 qaws_dual1_mul(qaws_dual1 a, qaws_dual1 b)
{
    return qaws_dual1_make(
        a.v * b.v,
        a.t * b.v + a.v * b.t,
        a.tt * b.v + QAWS_LITERAL(2.0) * a.t * b.t + a.v * b.tt);
}

/* q = a / b:  q' = (a' - q b') / b,  q'' = (a'' - 2 q' b' - q b'') / b */
QAWS_INLINE qaws_dual1 qaws_dual1_div(qaws_dual1 a, qaws_dual1 b)
{
    qaws_dual1 q;
    q.v = a.v / b.v;
    q.t = (a.t - q.v * b.t) / b.v;
    q.tt = (a.tt - QAWS_LITERAL(2.0) * q.t * b.t - q.v * b.tt) / b.v;
    return q;
}

/* s = sqrt(a):  s' = a' / (2 s),  s'' = (a'' - 2 s'^2) / (2 s).  Guarded at 0. */
QAWS_INLINE qaws_dual1 qaws_dual1_sqrt(qaws_dual1 a)
{
    qaws_dual1 s;
    qaws_scalar two_s;
    s.v = QAWS_SQRT(a.v);
    two_s = QAWS_LITERAL(2.0) * s.v;
    s.t = QAWS_SELECT(two_s > QAWS_ZERO, a.t / two_s, QAWS_ZERO);
    s.tt = QAWS_SELECT(two_s > QAWS_ZERO,
                       (a.tt - QAWS_LITERAL(2.0) * s.t * s.t) / two_s, QAWS_ZERO);
    return s;
}

/* -------------------------------------------------------------------
 * Vector rules
 * ------------------------------------------------------------------- */

QAWS_INLINE qaws_dual3 qaws_dual3_add(qaws_dual3 a, qaws_dual3 b)
{
    return qaws_dual3_make(qaws_v3_add(a.v, b.v), qaws_v3_add(a.t, b.t), qaws_v3_add(a.tt, b.tt));
}

QAWS_INLINE qaws_dual3 qaws_dual3_sub(qaws_dual3 a, qaws_dual3 b)
{
    return qaws_dual3_make(qaws_v3_sub(a.v, b.v), qaws_v3_sub(a.t, b.t), qaws_v3_sub(a.tt, b.tt));
}

/* x * s for a dual scalar s */
QAWS_INLINE qaws_dual3 qaws_dual3_scale(qaws_dual3 x, qaws_dual1 s)
{
    qaws_dual3 r;
    r.v = qaws_v3_scale(x.v, s.v);
    r.t = qaws_v3_add(qaws_v3_scale(x.t, s.v), qaws_v3_scale(x.v, s.t));
    r.tt = qaws_v3_add(qaws_v3_add(qaws_v3_scale(x.tt, s.v),
                                   qaws_v3_scale(x.t, QAWS_LITERAL(2.0) * s.t)),
                       qaws_v3_scale(x.v, s.tt));
    return r;
}

QAWS_INLINE qaws_dual1 qaws_dual3_dot(qaws_dual3 a, qaws_dual3 b)
{
    return qaws_dual1_make(
        qaws_v3_dot(a.v, b.v),
        qaws_v3_dot(a.t, b.v) + qaws_v3_dot(a.v, b.t),
        qaws_v3_dot(a.tt, b.v) + QAWS_LITERAL(2.0) * qaws_v3_dot(a.t, b.t) + qaws_v3_dot(a.v, b.tt));
}

QAWS_INLINE qaws_dual3 qaws_dual3_cross(qaws_dual3 a, qaws_dual3 b)
{
    qaws_dual3 r;
    r.v = qaws_v3_cross(a.v, b.v);
    r.t = qaws_v3_add(qaws_v3_cross(a.t, b.v), qaws_v3_cross(a.v, b.t));
    r.tt = qaws_v3_add(qaws_v3_add(qaws_v3_cross(a.tt, b.v),
                                   qaws_v3_scale(qaws_v3_cross(a.t, b.t), QAWS_LITERAL(2.0))),
                       qaws_v3_cross(a.v, b.tt));
    return r;
}

/* x / w for a dual scalar w:  p' = (x' - p w') / w,  p'' = (x'' - 2 p' w' - p w'') / w */
QAWS_INLINE qaws_dual3 qaws_dual3_div(qaws_dual3 x, qaws_dual1 w)
{
    qaws_dual3 p;
    qaws_scalar inv = QAWS_ONE / w.v;
    p.v = qaws_v3_scale(x.v, inv);
    p.t = qaws_v3_scale(qaws_v3_sub(x.t, qaws_v3_scale(p.v, w.t)), inv);
    p.tt = qaws_v3_scale(
        qaws_v3_sub(qaws_v3_sub(x.tt, qaws_v3_scale(p.t, QAWS_LITERAL(2.0) * w.t)),
                    qaws_v3_scale(p.v, w.tt)),
        inv);
    return p;
}

QAWS_INLINE qaws_dual1 qaws_dual3_length(qaws_dual3 a)
{
    return qaws_dual1_sqrt(qaws_dual3_dot(a, a));
}

/* Multiply every coefficient by a constant. */
QAWS_INLINE qaws_dual3 qaws_dual3_mul_const(qaws_dual3 a, qaws_scalar s)
{
    return qaws_dual3_make(qaws_v3_scale(a.v, s), qaws_v3_scale(a.t, s), qaws_v3_scale(a.tt, s));
}

/* Unit vector. Degenerate input returns a zero dual (callers choose a fallback). */
QAWS_INLINE qaws_dual3 qaws_dual3_normalize(qaws_dual3 a)
{
    qaws_dual1 len = qaws_dual3_length(a);
    qaws_scalar ok = QAWS_SELECT(len.v > QAWS_LITERAL(1e-30), QAWS_ONE, QAWS_ZERO);
    len.v = QAWS_SELECT(len.v > QAWS_LITERAL(1e-30), len.v, QAWS_ONE);
    return qaws_dual3_mul_const(qaws_dual3_div(a, len), ok);
}

/* -------------------------------------------------------------------
 * Adjoint rules (first order, primal values given)
 * ------------------------------------------------------------------- */

/* s = a . b */
QAWS_INLINE qaws_vec3_pair qaws_dot_adjoint(qaws_vec3 a, qaws_vec3 b, qaws_scalar s_bar)
{
    qaws_vec3_pair r;
    r.a = qaws_v3_scale(b, s_bar);
    r.b = qaws_v3_scale(a, s_bar);
    return r;
}

/* c = a x b:  a_bar = b x c_bar,  b_bar = c_bar x a */
QAWS_INLINE qaws_vec3_pair qaws_cross_adjoint(qaws_vec3 a, qaws_vec3 b, qaws_vec3 c_bar)
{
    qaws_vec3_pair r;
    r.a = qaws_v3_cross(b, c_bar);
    r.b = qaws_v3_cross(c_bar, a);
    return r;
}

/* l = |a| */
QAWS_INLINE qaws_vec3 qaws_length_adjoint(qaws_vec3 a, qaws_scalar l_bar)
{
    qaws_scalar len = QAWS_SQRT(qaws_v3_dot(a, a));
    qaws_scalar inv = QAWS_SELECT(len > QAWS_LITERAL(1e-30), QAWS_ONE / len, QAWS_ZERO);
    return qaws_v3_scale(a, l_bar * inv);
}

/* n = a / |a|:  a_bar = (n_bar - n (n . n_bar)) / |a| */
QAWS_INLINE qaws_vec3 qaws_normalize_adjoint(qaws_vec3 a, qaws_vec3 n_bar)
{
    qaws_scalar len = QAWS_SQRT(qaws_v3_dot(a, a));
    qaws_scalar inv = QAWS_SELECT(len > QAWS_LITERAL(1e-30), QAWS_ONE / len, QAWS_ZERO);
    qaws_vec3 n = qaws_v3_scale(a, inv);
    return qaws_v3_scale(qaws_v3_sub(n_bar, qaws_v3_scale(n, qaws_v3_dot(n, n_bar))), inv);
}

/* p = x / w:  x_bar = p_bar / w,  w_bar = -(p_bar . p) / w */
QAWS_INLINE qaws_vec3_scalar qaws_div_adjoint(qaws_vec3 x, qaws_scalar w, qaws_vec3 p_bar)
{
    qaws_vec3_scalar r;
    qaws_scalar inv = QAWS_ONE / w;
    qaws_vec3 p = qaws_v3_scale(x, inv);
    r.x = qaws_v3_scale(p_bar, inv);
    r.s = -qaws_v3_dot(p_bar, p) * inv;
    return r;
}

/* y = x * s:  x_bar = s y_bar,  s_bar = x . y_bar */
QAWS_INLINE qaws_vec3_scalar qaws_scale_adjoint(qaws_vec3 x, qaws_scalar s, qaws_vec3 y_bar)
{
    qaws_vec3_scalar r;
    r.x = qaws_v3_scale(y_bar, s);
    r.s = qaws_v3_dot(x, y_bar);
    return r;
}

#endif /* QAWS_DUAL_CORE_H */
