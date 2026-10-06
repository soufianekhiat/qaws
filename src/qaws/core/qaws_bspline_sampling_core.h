#ifndef QAWS_BSPLINE_SAMPLING_CORE_H
#define QAWS_BSPLINE_SAMPLING_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"
#include "qaws_bspline_diff_core.h"
#include "qaws_dual_core.h"

/*
 * Inverse-CDF sampling and integral functionals of B-spline curves, as
 * kernels for every backend (C, HLSL, GLSL). Same window conventions as
 * qaws_bspline_diff_core.h: local_cp holds the degree + 1 control points of
 * the span, local_knots the knot window around it, span the span index in
 * that window; local_cp_dot holds the control point tangents (zeros for
 * none). Results match qaws_curve_cdf_sample_* and qaws_curve_functional_*
 * with the default quadrature: every span is split into
 * QAWS_SAMPLING_PIECES pieces of 8 Gauss points.
 *
 * Integrands f(C', C''), on second-order dual numbers along the control
 * point tangents:
 *   QAWS_CORE_INTEGRAND_SPEED         |C'|: arc-length measure, length
 *   QAWS_CORE_INTEGRAND_CURVATURE     sqrt(floor^2 |C'|^2 + |C' x C''|^2 / |C'|^4)
 *                                     (curvature measure)
 *   QAWS_CORE_INTEGRAND_BENDING       |C''|^2
 *   QAWS_CORE_INTEGRAND_CURVATURE_SQ  |C' x C''|^2 / |C'|^5
 *
 * Kernels (one thread per span, per sample or per control point):
 *   qaws_bspline_integrate        value, tangent and second tangent of the
 *                                 integral over [a, b] inside one span; the
 *                                 span totals give the prefix measures
 *   qaws_bspline_cdf_solve        parameter at which the measure from a
 *                                 reaches `remaining` (safeguarded Newton,
 *                                 fixed iteration count, branch free)
 *   qaws_bspline_cdf_tangent      first and second tangents of a sample from
 *                                 the target and measure rates (implicit
 *                                 function theorem)
 *   qaws_bspline_cdf_lambda       adjoint multiplier (t_bar + p_bar . C') / m
 *   qaws_bspline_node_adjoint_cp  c (N_j' df/dC' + N_j'' df/dC'') at one
 *                                 quadrature node: the per-element term of a
 *                                 gather over the nodes a control point
 *                                 supports (functional gradients, and the
 *                                 measure terms of the sampling adjoint)
 */

#define QAWS_SAMPLING_PIECES 8

#define QAWS_CORE_INTEGRAND_SPEED        0
#define QAWS_CORE_INTEGRAND_CURVATURE    1
#define QAWS_CORE_INTEGRAND_BENDING      2
#define QAWS_CORE_INTEGRAND_CURVATURE_SQ 3

/* 8-point Gauss-Legendre rule on [-1, 1]: xw[0] node, xw[1] weight. */
QAWS_INLINE void qaws_gauss8(int k, QAWS_OUT qaws_scalar xw[2])
{
    int m = k < 4 ? 3 - k : k - 4;
    qaws_scalar s = k < 4 ? -QAWS_ONE : QAWS_ONE;
    qaws_scalar x = QAWS_LITERAL(0.1834346424956498);
    qaws_scalar w = QAWS_LITERAL(0.3626837833783620);
    if (m == 1) { x = QAWS_LITERAL(0.5255324099163290); w = QAWS_LITERAL(0.3137066458778873); }
    if (m == 2) { x = QAWS_LITERAL(0.7966664774136267); w = QAWS_LITERAL(0.2223810344533745); }
    if (m == 3) { x = QAWS_LITERAL(0.9602898564975363); w = QAWS_LITERAL(0.1012285362903763); }
    xw[0] = s * x;
    xw[1] = w;
}

/* f(C', C'') on dual jets. */
QAWS_INLINE qaws_dual1 qaws_core_integrand(int id, qaws_scalar curvature_floor, qaws_dual3 d1, qaws_dual3 d2)
{
    qaws_dual1 s = qaws_dual3_length(d1);
    qaws_dual3 x = qaws_dual3_cross(d1, d2);
    qaws_dual1 xx = qaws_dual3_dot(x, x);
    qaws_dual1 s2 = qaws_dual1_mul(s, s);
    qaws_dual1 s4 = qaws_dual1_mul(s2, s2);
    qaws_dual1 r = s;
    if (id == QAWS_CORE_INTEGRAND_CURVATURE) {
        qaws_scalar f2 = curvature_floor * curvature_floor;
        r = qaws_dual1_sqrt(qaws_dual1_add(qaws_dual1_make(f2 * s2.v, f2 * s2.t, f2 * s2.tt), qaws_dual1_div(xx, s4)));
    }
    if (id == QAWS_CORE_INTEGRAND_BENDING)
        r = qaws_dual3_dot(d2, d2);
    if (id == QAWS_CORE_INTEGRAND_CURVATURE_SQ)
        r = qaws_dual1_div(xx, qaws_dual1_mul(s4, s));
    return r;
}

/* Jets (C', C'', C''') at t and the first two along the control point
   tangents. */
QAWS_INLINE void qaws_bspline_jet_rows(
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    QAWS_OUT qaws_vec3 c[4],
    QAWS_OUT qaws_vec3 c_dot[3])
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int k;
    qaws_bspline_diff_rows(local_knots, degree, span, t, ders);
    for (k = 0; k < 4; k++)
        c[k] = qaws_bspline_diff_row_sum(ders, local_cp, degree, k);
    for (k = 0; k < 3; k++)
        c_dot[k] = qaws_bspline_diff_row_sum(ders, local_cp_dot, degree, k);
}

/* f at t along the control point tangents (v, t, tt). */
QAWS_INLINE qaws_dual1 qaws_bspline_integrand_at(
    int id,
    qaws_scalar curvature_floor,
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t)
{
    qaws_vec3 c[4];
    qaws_vec3 c_dot[3];
    qaws_bspline_jet_rows(local_cp, local_cp_dot, local_knots, degree, span, t, c, c_dot);
    return qaws_core_integrand(id, curvature_floor,
        qaws_dual3_make(c[1], c_dot[1], qaws_v3_zero()), qaws_dual3_make(c[2], c_dot[2], qaws_v3_zero()));
}

/* Integral of f over [a, b] inside one span, with its first and second
   derivatives along the control point tangents (zero when b <= a). */
QAWS_INLINE qaws_dual1 qaws_bspline_integrate(
    int id,
    qaws_scalar curvature_floor,
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar a,
    qaws_scalar b)
{
    qaws_dual1 sum = qaws_dual1_const(QAWS_ZERO);
    qaws_scalar h = (b - a) / QAWS_LITERAL(8.0);
    qaws_scalar lo = a;
    qaws_scalar xw[2];
    qaws_dual1 f;
    qaws_scalar W;
    int piece, k;
    for (piece = 0; piece < QAWS_SAMPLING_PIECES; piece++) {
        for (k = 0; k < 8; k++) {
            qaws_gauss8(k, xw);
            W = QAWS_LITERAL(0.5) * h * xw[1];
            f = qaws_bspline_integrand_at(id, curvature_floor, local_cp, local_cp_dot, local_knots, degree, span,
                lo + QAWS_LITERAL(0.5) * h * (QAWS_ONE + xw[0]));
            sum = qaws_dual1_make(sum.v + W * f.v, sum.t + W * f.t, sum.tt + W * f.tt);
        }
        lo = lo + h;
    }
    sum.v = QAWS_SELECT(b > a, sum.v, QAWS_ZERO);
    sum.t = QAWS_SELECT(b > a, sum.t, QAWS_ZERO);
    sum.tt = QAWS_SELECT(b > a, sum.tt, QAWS_ZERO);
    return sum;
}

/* Parameter in [a, b] (one span) at which the measure from a reaches
   `remaining`: safeguarded Newton, fixed iteration count. */
QAWS_INLINE qaws_scalar qaws_bspline_cdf_solve(
    int id,
    qaws_scalar curvature_floor,
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar a,
    qaws_scalar b,
    qaws_scalar remaining)
{
    qaws_scalar zero_dot[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar lo = a, hi = b, t, f, m, next, len, total;
    int it;
    for (it = 0; it < QAWS_CORE_MAX_POINTS * 3; it++)
        zero_dot[it] = QAWS_ZERO;
    total = qaws_bspline_integrate(id, curvature_floor, local_cp, zero_dot, local_knots, degree, span, a, b).v;
    len = QAWS_SELECT(total > QAWS_ZERO, remaining / total, QAWS_ZERO);
    t = a + (b - a) * QAWS_SELECT(len < QAWS_ZERO, QAWS_ZERO, QAWS_SELECT(len > QAWS_ONE, QAWS_ONE, len));
    for (it = 0; it < 40; it++) {
        f = qaws_bspline_integrate(id, curvature_floor, local_cp, zero_dot, local_knots, degree, span, a, t).v - remaining;
        hi = QAWS_SELECT(f > QAWS_ZERO, t, hi);
        lo = QAWS_SELECT(f > QAWS_ZERO, lo, t);
        m = qaws_bspline_integrand_at(id, curvature_floor, local_cp, zero_dot, local_knots, degree, span, t).v;
        next = QAWS_SELECT(m > QAWS_ZERO, t - f / QAWS_SELECT(m > QAWS_ZERO, m, QAWS_ONE), QAWS_LITERAL(0.5) * (lo + hi));
        next = QAWS_SELECT(next > lo && next < hi, next, QAWS_LITERAL(0.5) * (lo + hi));
        t = QAWS_SELECT(f == QAWS_ZERO, t, next);
    }
    return t;
}

/*
 * Tangents of the sample at t: with sigma the target measure and M(t) the
 * measure up to t, m = M_t, the implicit function theorem gives
 *   t'  = (sigma' - M') / m,
 *   t'' = (sigma'' - M'' - 2 m_e t' - m_t t'^2) / m
 * (M', M'' at fixed t: the prefix of the full spans plus
 * qaws_bspline_integrate over [a_k, t]; m_e the rate of m along the
 * control point tangents, m_t its t-derivative). out_t = (t', t''),
 * out_p = (position tangent, second tangent).
 */
QAWS_INLINE void qaws_bspline_cdf_tangent(
    int id,
    qaws_scalar curvature_floor,
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_cp_dot[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    qaws_scalar sigma_dot,
    qaws_scalar sigma_dot2,
    qaws_scalar measure_dot,
    qaws_scalar measure_dot2,
    QAWS_OUT qaws_scalar out_t[2],
    QAWS_OUT qaws_vec3 out_p[2])
{
    qaws_vec3 c[4];
    qaws_vec3 c_dot[3];
    qaws_dual1 me, mt;
    qaws_scalar td, tdd, m;
    qaws_bspline_jet_rows(local_cp, local_cp_dot, local_knots, degree, span, t, c, c_dot);
    me = qaws_core_integrand(id, curvature_floor,
        qaws_dual3_make(c[1], c_dot[1], qaws_v3_zero()), qaws_dual3_make(c[2], c_dot[2], qaws_v3_zero()));
    mt = qaws_core_integrand(id, curvature_floor,
        qaws_dual3_make(c[1], c[2], qaws_v3_zero()), qaws_dual3_make(c[2], c[3], qaws_v3_zero()));
    m = QAWS_SELECT(me.v > QAWS_ZERO, me.v, QAWS_ONE);
    td = (sigma_dot - measure_dot) / m;
    tdd = (sigma_dot2 - measure_dot2 - QAWS_LITERAL(2.0) * me.t * td - mt.t * td * td) / m;
    out_t[0] = QAWS_SELECT(me.v > QAWS_ZERO, td, QAWS_ZERO);
    out_t[1] = QAWS_SELECT(me.v > QAWS_ZERO, tdd, QAWS_ZERO);
    /* p' = C_dot + C' t',  p'' = C'' t'^2 + 2 C'_dot t' + C' t'' */
    out_p[0] = qaws_v3_axpy(c_dot[0], c[1], out_t[0]);
    out_p[1] = qaws_v3_axpy(qaws_v3_axpy(qaws_v3_scale(c[2], out_t[0] * out_t[0]), c_dot[1], QAWS_LITERAL(2.0) * out_t[0]),
        c[1], out_t[1]);
}

/* lambda = (t_bar + p_bar . C') / m of one sample (0 where m vanishes). */
QAWS_INLINE qaws_scalar qaws_bspline_cdf_lambda(
    int id,
    qaws_scalar curvature_floor,
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    qaws_scalar t_bar,
    qaws_vec3 p_bar)
{
    qaws_vec3 c[4];
    qaws_vec3 c_dot[3];
    qaws_scalar m;
    qaws_bspline_jet_rows(local_cp, local_cp, local_knots, degree, span, t, c, c_dot);
    m = qaws_core_integrand(id, curvature_floor, qaws_dual3_const(c[1]), qaws_dual3_const(c[2])).v;
    return QAWS_SELECT(m > QAWS_ZERO, (t_bar + qaws_v3_dot(p_bar, c[1])) / QAWS_SELECT(m > QAWS_ZERO, m, QAWS_ONE), QAWS_ZERO);
}

/* c (N_j' df/dC' + N_j'' df/dC'') at the node t for local control point j. */
QAWS_INLINE qaws_vec3 qaws_bspline_node_adjoint_cp(
    int id,
    qaws_scalar curvature_floor,
    qaws_scalar local_cp[QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    int span,
    qaws_scalar t,
    int j,
    qaws_scalar c)
{
    qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_vec3 c1, c2, e, g;
    qaws_scalar g1[3], g2[3];
    int i;
    qaws_bspline_diff_rows(local_knots, degree, span, t, ders);
    c1 = qaws_bspline_diff_row_sum(ders, local_cp, degree, 1);
    c2 = qaws_bspline_diff_row_sum(ders, local_cp, degree, 2);
    /* gradient by seeding each component of C' and C'' */
    for (i = 0; i < 3; i++) {
        e = qaws_v3(i == 0 ? QAWS_ONE : QAWS_ZERO, i == 1 ? QAWS_ONE : QAWS_ZERO, i == 2 ? QAWS_ONE : QAWS_ZERO);
        g1[i] = qaws_core_integrand(id, curvature_floor, qaws_dual3_make(c1, e, qaws_v3_zero()), qaws_dual3_const(c2)).t;
        g2[i] = qaws_core_integrand(id, curvature_floor, qaws_dual3_const(c1), qaws_dual3_make(c2, e, qaws_v3_zero())).t;
    }
    i = degree + 1;
    g.x = c * (ders[1 * i + j] * g1[0] + ders[2 * i + j] * g2[0]);
    g.y = c * (ders[1 * i + j] * g1[1] + ders[2 * i + j] * g2[1]);
    g.z = c * (ders[1 * i + j] * g1[2] + ders[2 * i + j] * g2[2]);
    return g;
}

#endif /* QAWS_BSPLINE_SAMPLING_CORE_H */
