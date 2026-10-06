#ifndef QAWS_BSPLINE_SURFACE_SAMPLING_CORE_H
#define QAWS_BSPLINE_SURFACE_SAMPLING_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"
#include "qaws_dual_core.h"
#include "qaws_bspline_sampling_core.h"

/*
 * Inverse-CDF warping of the unit square onto a B-spline surface by its
 * area measure w = |S_u x S_v|, as kernels for every backend (C, HLSL,
 * GLSL). They match qaws_surface_cdf_sample_* of the runtime with the area
 * measure, `cells` cells per direction and 8 Gauss points per cell:
 *   A(u) = xi_u A_total,  B(v; u) = xi_v B(v1; u),
 *   A(u) = integral_{u0}^{u} B(v1; .),  B(b; u) = integral_{v0}^{b} w(u, .).
 *
 * The whole control grid is passed (a v-line crosses several knot spans):
 * cp[(i * v_count + j) * 3 + c], i along u, at most QAWS_CORE_MAX_POINTS
 * per direction; cp_dot the control point tangents; dir scales them (0:
 * value and parameter rates only). Duals carry first and second rates
 * along one path: control points moving at cp_dot, and u, v at their own
 * rates.
 *
 * Kernels (one thread per u cell, per sample, per control point):
 *   qaws_patch_cell_mass      measure of u cell c with its rates; an
 *                             exclusive prefix sum gives prefix_v/_t/_tt
 *   qaws_patch_cdf_solve      (u, v) of a point xi (fixed iteration count)
 *   qaws_patch_cdf_tangent    first and second tangents of (u, v) and of
 *                             the position along cp_dot and xi_dot
 *   qaws_patch_cdf_multipliers  adjoint multipliers (lambda, mu) and the
 *                             xi adjoints of one sample
 *   qaws_patch_line_adjoint_cp  c times the gradient of B(b; u) for one
 *                             control point: the per-line term of the
 *                             gather (conditional lines at a sample, the
 *                             lines of its partial u cell, the lines of
 *                             the full u cells)
 */

#define QAWS_CORE_SURFACE_GRID (QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3)
#define QAWS_CORE_MAX_CELLS 32

#define QAWS_PATCH_PARAMS \
    qaws_scalar cp[QAWS_CORE_SURFACE_GRID], qaws_scalar cp_dot[QAWS_CORE_SURFACE_GRID], \
    qaws_scalar u_knots[QAWS_CORE_MAX_POINTS * 2], qaws_scalar v_knots[QAWS_CORE_MAX_POINTS * 2], \
    int u_degree, int v_degree, int u_count, int v_count
#define QAWS_PATCH_ARGS cp, cp_dot, u_knots, v_knots, u_degree, v_degree, u_count, v_count

/* sum_ij N_i^(a)(u) N_j^(b)(v) P_ij over the spans (us, vs). */
QAWS_INLINE qaws_vec3 qaws_patch_sum(
    qaws_scalar grid[QAWS_CORE_SURFACE_GRID],
    qaws_scalar nu[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS],
    qaws_scalar nv[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS],
    int u_degree, int v_degree, int v_count, int us, int vs, int a, int b)
{
    qaws_vec3 r = qaws_v3_zero();
    int i, j, at, uo = u_degree + 1, vo = v_degree + 1;
    qaws_scalar w;
    for (i = 0; i < uo; i++)
        for (j = 0; j < vo; j++) {
            w = nu[a * uo + i] * nv[b * vo + j];
            at = ((us - u_degree + i) * v_count + (vs - v_degree + j)) * 3;
            r = qaws_v3_axpy(r, qaws_v3(grid[at], grid[at + 1], grid[at + 2]), w);
        }
    return r;
}

/* Basis rows 0..3 at t (rows above the degree zero). */
QAWS_INLINE void qaws_patch_rows(
    qaws_scalar knots[QAWS_CORE_MAX_POINTS * 2], int degree, int span, qaws_scalar t,
    QAWS_OUT qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS])
{
    int max_k = degree < QAWS_CORE_MAX_DERIV ? degree : QAWS_CORE_MAX_DERIV;
    int k, j;
    qaws_bspline_basis_derivs(knots, degree, span, t, max_k, ders);
    for (k = max_k + 1; k <= QAWS_CORE_MAX_DERIV; k++)
        for (j = 0; j <= degree; j++)
            ders[k * (degree + 1) + j] = QAWS_ZERO;
}

/* Dual (S, S_u, S_v) at the dual point (u, v): jets along the path. */
QAWS_INLINE void qaws_patch_jets(QAWS_PATCH_PARAMS, qaws_scalar dir, qaws_dual1 u, qaws_dual1 v, QAWS_OUT qaws_dual3 y[3])
{
    qaws_scalar nu[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar nv[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int us = qaws_find_span(u_knots, u_degree, u_count, u.v);
    int vs = qaws_find_span(v_knots, v_degree, v_count, v.v);
    qaws_vec3 S00, S10, S01, S20, S11, S02, S30, S21, S12, S03, D00, D10, D01, D20, D11, D02;
    qaws_scalar ut = u.t, vt = v.t;
    qaws_patch_rows(u_knots, u_degree, us, u.v, nu);
    qaws_patch_rows(v_knots, v_degree, vs, v.v, nv);
    S00 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 0);
    S10 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 1, 0);
    S01 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 1);
    S20 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 2, 0);
    S11 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 1, 1);
    S02 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 2);
    S30 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 3, 0);
    S21 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 2, 1);
    S12 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 1, 2);
    S03 = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 3);
    D00 = qaws_v3_scale(qaws_patch_sum(cp_dot, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 0), dir);
    D10 = qaws_v3_scale(qaws_patch_sum(cp_dot, nu, nv, u_degree, v_degree, v_count, us, vs, 1, 0), dir);
    D01 = qaws_v3_scale(qaws_patch_sum(cp_dot, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 1), dir);
    D20 = qaws_v3_scale(qaws_patch_sum(cp_dot, nu, nv, u_degree, v_degree, v_count, us, vs, 2, 0), dir);
    D11 = qaws_v3_scale(qaws_patch_sum(cp_dot, nu, nv, u_degree, v_degree, v_count, us, vs, 1, 1), dir);
    D02 = qaws_v3_scale(qaws_patch_sum(cp_dot, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 2), dir);
    /* x = X(u, v) + e X_dot: x' = X_dot + X_u u' + X_v v',
       x'' = 2 (X_dot,u u' + X_dot,v v') + X_uu u'^2 + 2 X_uv u' v' + X_vv v'^2 + X_u u'' + X_v v'' */
    y[0].v = S00;
    y[0].t = qaws_v3_axpy(qaws_v3_axpy(D00, S10, ut), S01, vt);
    y[0].tt = qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_zero(),
        D10, QAWS_LITERAL(2.0) * ut), D01, QAWS_LITERAL(2.0) * vt), S20, ut * ut), S11, QAWS_LITERAL(2.0) * ut * vt), S02, vt * vt),
        S10, u.tt), S01, v.tt);
    y[1].v = S10;
    y[1].t = qaws_v3_axpy(qaws_v3_axpy(D10, S20, ut), S11, vt);
    y[1].tt = qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_zero(),
        D20, QAWS_LITERAL(2.0) * ut), D11, QAWS_LITERAL(2.0) * vt), S30, ut * ut), S21, QAWS_LITERAL(2.0) * ut * vt), S12, vt * vt),
        S20, u.tt), S11, v.tt);
    y[2].v = S01;
    y[2].t = qaws_v3_axpy(qaws_v3_axpy(D01, S11, ut), S02, vt);
    y[2].tt = qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_axpy(qaws_v3_zero(),
        D11, QAWS_LITERAL(2.0) * ut), D02, QAWS_LITERAL(2.0) * vt), S21, ut * ut), S12, QAWS_LITERAL(2.0) * ut * vt), S03, vt * vt),
        S11, u.tt), S02, v.tt);
}

/* w = |S_u x S_v| at the dual point. */
QAWS_INLINE qaws_dual1 qaws_patch_area_at(QAWS_PATCH_PARAMS, qaws_scalar dir, qaws_dual1 u, qaws_dual1 v)
{
    qaws_dual3 y[3];
    qaws_patch_jets(QAWS_PATCH_ARGS, dir, u, v, y);
    return qaws_dual3_length(qaws_dual3_cross(y[1], y[2]));
}

QAWS_INLINE qaws_dual1 qaws_dual_affine(qaws_scalar a, qaws_dual1 x, qaws_scalar b)
{
    return qaws_dual1_make(a * x.v + b, a * x.t, a * x.tt);
}

/* B(b; u) = integral_{v0}^{b} w(u, .): full v cells, then the partial cell
   whose nodes and weights move with b. */
QAWS_INLINE qaws_dual1 qaws_patch_line(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, qaws_dual1 u, qaws_dual1 b)
{
    qaws_scalar v0 = v_knots[v_degree], v1 = v_knots[v_count];
    qaws_scalar hv = (v1 - v0) / QAWS_ITOF(cells);
    qaws_scalar xw[2];
    qaws_scalar vk, base = v0;
    qaws_dual1 sum = qaws_dual1_const(QAWS_ZERO), len, f;
    int c, q, k;
    k = QAWS_FTOI((b.v - v0) / hv);
    k = k < 0 ? 0 : (k >= cells ? cells - 1 : k);
    for (c = 0; c < QAWS_CORE_MAX_CELLS; c++) {
        if (c < k) {
            for (q = 0; q < 8; q++) {
                qaws_gauss8(q, xw);
                f = qaws_patch_area_at(QAWS_PATCH_ARGS, dir, u, qaws_dual1_const(base + QAWS_LITERAL(0.5) * hv * (QAWS_ONE + xw[0])));
                sum = qaws_dual1_add(sum, qaws_dual_affine(QAWS_LITERAL(0.5) * hv * xw[1], f, QAWS_ZERO));
            }
        }
        base = base + QAWS_SELECT(c < k, hv, QAWS_ZERO);
    }
    vk = base;
    len = qaws_dual_affine(QAWS_ONE, b, -vk);
    for (q = 0; q < 8; q++) {
        qaws_gauss8(q, xw);
        f = qaws_patch_area_at(QAWS_PATCH_ARGS, dir, u, qaws_dual_affine(QAWS_LITERAL(0.5) * (QAWS_ONE + xw[0]), len, vk));
        sum = qaws_dual1_add(sum, qaws_dual1_mul(qaws_dual_affine(QAWS_LITERAL(0.5) * xw[1], len, QAWS_ZERO), f));
    }
    return sum;
}

/* Measure of u cell c: integral over the cell of B(v1; .). */
QAWS_INLINE qaws_dual1 qaws_patch_cell_mass(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, int c)
{
    qaws_scalar u0 = u_knots[u_degree], u1 = u_knots[u_count], v1 = v_knots[v_count];
    qaws_scalar hu = (u1 - u0) / QAWS_ITOF(cells);
    qaws_scalar lo = u0 + hu * QAWS_ITOF(c);
    qaws_scalar xw[2];
    qaws_dual1 sum = qaws_dual1_const(QAWS_ZERO), b;
    int q;
    for (q = 0; q < 8; q++) {
        qaws_gauss8(q, xw);
        b = qaws_patch_line(QAWS_PATCH_ARGS, dir, cells, qaws_dual1_const(lo + QAWS_LITERAL(0.5) * hu * (QAWS_ONE + xw[0])),
            qaws_dual1_const(v1));
        sum = qaws_dual1_add(sum, qaws_dual_affine(QAWS_LITERAL(0.5) * hu * xw[1], b, QAWS_ZERO));
    }
    return sum;
}

/* A(u) = prefix[k] + the partial u cell moving with u. */
QAWS_INLINE qaws_dual1 qaws_patch_marginal(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells,
    qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar prefix_t[QAWS_CORE_MAX_CELLS + 1],
    qaws_scalar prefix_tt[QAWS_CORE_MAX_CELLS + 1], qaws_dual1 u)
{
    qaws_scalar u0 = u_knots[u_degree], u1 = u_knots[u_count], v1 = v_knots[v_count];
    qaws_scalar hu = (u1 - u0) / QAWS_ITOF(cells);
    qaws_scalar xw[2];
    qaws_scalar uk;
    qaws_dual1 sum, len, b;
    int k, q;
    k = QAWS_FTOI((u.v - u0) / hu);
    k = k < 0 ? 0 : (k >= cells ? cells - 1 : k);
    uk = u0 + hu * QAWS_ITOF(k);
    sum = qaws_dual1_make(prefix_v[k], dir * prefix_t[k], dir * prefix_tt[k]);
    len = qaws_dual_affine(QAWS_ONE, u, -uk);
    for (q = 0; q < 8; q++) {
        qaws_gauss8(q, xw);
        b = qaws_patch_line(QAWS_PATCH_ARGS, dir, cells, qaws_dual_affine(QAWS_LITERAL(0.5) * (QAWS_ONE + xw[0]), len, uk),
            qaws_dual1_const(v1));
        sum = qaws_dual1_add(sum, qaws_dual1_mul(qaws_dual_affine(QAWS_LITERAL(0.5) * xw[1], len, QAWS_ZERO), b));
    }
    return sum;
}

/* F(u) = A(u) - xi_u A_total. */
QAWS_INLINE qaws_dual1 qaws_patch_F(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells,
    qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar prefix_t[QAWS_CORE_MAX_CELLS + 1],
    qaws_scalar prefix_tt[QAWS_CORE_MAX_CELLS + 1], qaws_dual1 u, qaws_dual1 xi)
{
    qaws_dual1 total = qaws_dual1_make(prefix_v[cells], dir * prefix_t[cells], dir * prefix_tt[cells]);
    return qaws_dual1_sub(qaws_patch_marginal(QAWS_PATCH_ARGS, dir, cells, prefix_v, prefix_t, prefix_tt, u),
        qaws_dual1_mul(xi, total));
}

/* G(v; u) = B(v; u) - xi_v B(v1; u). */
QAWS_INLINE qaws_dual1 qaws_patch_G(QAWS_PATCH_PARAMS, qaws_scalar dir, int cells, qaws_dual1 u, qaws_dual1 v, qaws_dual1 xi)
{
    qaws_scalar v1 = v_knots[v_count];
    return qaws_dual1_sub(qaws_patch_line(QAWS_PATCH_ARGS, dir, cells, u, v),
        qaws_dual1_mul(xi, qaws_patch_line(QAWS_PATCH_ARGS, dir, cells, u, qaws_dual1_const(v1))));
}

/* (u, v) of the point xi: safeguarded Newton in the cell holding the
   target, fixed iteration count. out_uv = (u, v). */
QAWS_INLINE void qaws_patch_cdf_solve(QAWS_PATCH_PARAMS, int cells,
    qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar xu, qaws_scalar xv, QAWS_OUT qaws_scalar out_uv[2])
{
    qaws_scalar u0 = u_knots[u_degree], u1 = u_knots[u_count], v0 = v_knots[v_degree], v1 = v_knots[v_count];
    qaws_scalar hu = (u1 - u0) / QAWS_ITOF(cells), hv = (v1 - v0) / QAWS_ITOF(cells);
    qaws_scalar total = prefix_v[cells], target, lo, hi, x, f, s, b1, bk;
    qaws_scalar zeros[QAWS_CORE_MAX_CELLS + 1];
    qaws_dual1 r;
    int k = 0, it, c;
    for (c = 0; c <= QAWS_CORE_MAX_CELLS; c++)
        zeros[c] = QAWS_ZERO;
    xu = xu < QAWS_ZERO ? QAWS_ZERO : (xu > QAWS_ONE ? QAWS_ONE : xu);
    xv = xv < QAWS_ZERO ? QAWS_ZERO : (xv > QAWS_ONE ? QAWS_ONE : xv);
    target = xu * total;
    for (c = 0; c + 1 < QAWS_CORE_MAX_CELLS; c++)
        if (c + 1 < cells && c == k && prefix_v[c + 1] <= target)
            k = c + 1;
    lo = u0 + hu * QAWS_ITOF(k);
    hi = lo + hu;
    x = QAWS_LITERAL(0.5) * (lo + hi);
    for (it = 0; it < 40; it++) {
        r = qaws_patch_F(QAWS_PATCH_ARGS, QAWS_ZERO, cells, prefix_v, zeros, zeros, qaws_dual1_make(x, QAWS_ONE, QAWS_ZERO),
            qaws_dual1_const(xu));
        f = r.v;
        s = r.t;
        hi = QAWS_SELECT(f > QAWS_ZERO, x, hi);
        lo = QAWS_SELECT(f > QAWS_ZERO, lo, x);
        x = QAWS_SELECT(s > QAWS_ZERO, x - f / QAWS_SELECT(s > QAWS_ZERO, s, QAWS_ONE), QAWS_LITERAL(0.5) * (lo + hi));
        x = QAWS_SELECT(x > lo && x < hi, x, QAWS_LITERAL(0.5) * (lo + hi));
    }
    out_uv[0] = x;
    /* the v cell holding xi_v B(v1; u) */
    b1 = qaws_patch_line(QAWS_PATCH_ARGS, QAWS_ZERO, cells, qaws_dual1_const(x), qaws_dual1_const(v1)).v;
    k = 0;
    for (c = 0; c + 1 < QAWS_CORE_MAX_CELLS; c++)
        if (c + 1 < cells && c == k) {
            bk = qaws_patch_line(QAWS_PATCH_ARGS, QAWS_ZERO, cells, qaws_dual1_const(x),
                qaws_dual1_const(v0 + hv * QAWS_ITOF(c + 1))).v;
            if (!(bk > xv * b1))
                k = c + 1;
        }
    lo = v0 + hv * QAWS_ITOF(k);
    hi = lo + hv;
    s = QAWS_LITERAL(0.5) * (lo + hi);
    for (it = 0; it < 40; it++) {
        r = qaws_patch_G(QAWS_PATCH_ARGS, QAWS_ZERO, cells, qaws_dual1_const(x), qaws_dual1_make(s, QAWS_ONE, QAWS_ZERO),
            qaws_dual1_const(xv));
        f = r.v;
        hi = QAWS_SELECT(f > QAWS_ZERO, s, hi);
        lo = QAWS_SELECT(f > QAWS_ZERO, lo, s);
        s = QAWS_SELECT(r.t > QAWS_ZERO, s - f / QAWS_SELECT(r.t > QAWS_ZERO, r.t, QAWS_ONE), QAWS_LITERAL(0.5) * (lo + hi));
        s = QAWS_SELECT(s > lo && s < hi, s, QAWS_LITERAL(0.5) * (lo + hi));
    }
    out_uv[1] = s;
}

/*
 * Tangents of the sample (u, v) of xi along cp_dot and xi_dot: one dual
 * Newton step per order on the discrete equations,
 *   u' = -F'/F_u,  u'' = -F''(u moving at u')/F_u, then v likewise with u
 * moving. result = (u', v', u'', v''), out_p = (position tangent, second).
 */
QAWS_INLINE void qaws_patch_cdf_tangent(QAWS_PATCH_PARAMS, int cells,
    qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar prefix_t[QAWS_CORE_MAX_CELLS + 1],
    qaws_scalar prefix_tt[QAWS_CORE_MAX_CELLS + 1], qaws_scalar xu, qaws_scalar xv, qaws_scalar xu_dot, qaws_scalar xv_dot,
    qaws_scalar u, qaws_scalar v, QAWS_OUT qaws_scalar result[4], QAWS_OUT qaws_vec3 out_p[2])
{
    qaws_dual1 xiu = qaws_dual1_make(xu, xu_dot, QAWS_ZERO), xiv = qaws_dual1_make(xv, xv_dot, QAWS_ZERO);
    qaws_dual1 F, G, ud, vd;
    qaws_dual3 y[3];
    qaws_scalar Fu, Gv, u1, u2, v1, v2;
    F = qaws_patch_F(QAWS_PATCH_ARGS, QAWS_ZERO, cells, prefix_v, prefix_t, prefix_tt, qaws_dual1_make(u, QAWS_ONE, QAWS_ZERO),
        qaws_dual1_const(xu));
    Fu = QAWS_SELECT(F.t > QAWS_ZERO, F.t, QAWS_ONE);
    G = qaws_patch_G(QAWS_PATCH_ARGS, QAWS_ZERO, cells, qaws_dual1_const(u), qaws_dual1_make(v, QAWS_ONE, QAWS_ZERO), qaws_dual1_const(xv));
    Gv = QAWS_SELECT(G.t > QAWS_ZERO, G.t, QAWS_ONE);
    F = qaws_patch_F(QAWS_PATCH_ARGS, QAWS_ONE, cells, prefix_v, prefix_t, prefix_tt, qaws_dual1_const(u), xiu);
    u1 = -F.t / Fu;
    F = qaws_patch_F(QAWS_PATCH_ARGS, QAWS_ONE, cells, prefix_v, prefix_t, prefix_tt, qaws_dual1_make(u, u1, QAWS_ZERO), xiu);
    u2 = -F.tt / Fu;
    ud = qaws_dual1_make(u, u1, u2);
    G = qaws_patch_G(QAWS_PATCH_ARGS, QAWS_ONE, cells, ud, qaws_dual1_const(v), xiv);
    v1 = -G.t / Gv;
    G = qaws_patch_G(QAWS_PATCH_ARGS, QAWS_ONE, cells, ud, qaws_dual1_make(v, v1, QAWS_ZERO), xiv);
    v2 = -G.tt / Gv;
    vd = qaws_dual1_make(v, v1, v2);
    qaws_patch_jets(QAWS_PATCH_ARGS, QAWS_ONE, ud, vd, y);
    result[0] = u1;
    result[1] = v1;
    result[2] = u2;
    result[3] = v2;
    out_p[0] = y[0].t;
    out_p[1] = y[0].tt;
}

/*
 * Adjoint multipliers of one sample with adjoint (t_bar_u, t_bar_v, p_bar):
 *   v_bar = t_bar_v + p_bar . S_v,  mu = v_bar / G_v,
 *   u_bar = t_bar_u + p_bar . S_u - mu G_u,  lambda = u_bar / F_u.
 * result = (lambda, mu, xi_u adjoint lambda A_total, xi_v adjoint mu B(v1; u)).
 */
QAWS_INLINE void qaws_patch_cdf_multipliers(QAWS_PATCH_PARAMS, int cells,
    qaws_scalar prefix_v[QAWS_CORE_MAX_CELLS + 1], qaws_scalar xu, qaws_scalar xv, qaws_scalar u, qaws_scalar v,
    qaws_scalar u_bar, qaws_scalar v_bar, qaws_vec3 p_bar, QAWS_OUT qaws_scalar result[4])
{
    qaws_scalar zeros[QAWS_CORE_MAX_CELLS + 1];
    qaws_scalar v1 = v_knots[v_count];
    qaws_dual1 F, Gv, Gu;
    qaws_dual3 y[3];
    qaws_scalar ub, vb, mu, lambda;
    int c;
    for (c = 0; c <= QAWS_CORE_MAX_CELLS; c++)
        zeros[c] = QAWS_ZERO;
    F = qaws_patch_F(QAWS_PATCH_ARGS, QAWS_ZERO, cells, prefix_v, zeros, zeros, qaws_dual1_make(u, QAWS_ONE, QAWS_ZERO),
        qaws_dual1_const(xu));
    Gv = qaws_patch_G(QAWS_PATCH_ARGS, QAWS_ZERO, cells, qaws_dual1_const(u), qaws_dual1_make(v, QAWS_ONE, QAWS_ZERO), qaws_dual1_const(xv));
    Gu = qaws_patch_G(QAWS_PATCH_ARGS, QAWS_ZERO, cells, qaws_dual1_make(u, QAWS_ONE, QAWS_ZERO), qaws_dual1_const(v), qaws_dual1_const(xv));
    qaws_patch_jets(QAWS_PATCH_ARGS, QAWS_ZERO, qaws_dual1_const(u), qaws_dual1_const(v), y);
    vb = v_bar + qaws_v3_dot(p_bar, y[2].v);
    ub = u_bar + qaws_v3_dot(p_bar, y[1].v);
    mu = QAWS_SELECT(Gv.t > QAWS_ZERO, vb / QAWS_SELECT(Gv.t > QAWS_ZERO, Gv.t, QAWS_ONE), QAWS_ZERO);
    ub = ub - mu * Gu.t;
    lambda = QAWS_SELECT(F.t > QAWS_ZERO, ub / QAWS_SELECT(F.t > QAWS_ZERO, F.t, QAWS_ONE), QAWS_ZERO);
    result[0] = lambda;
    result[1] = mu;
    result[2] = lambda * prefix_v[cells];
    result[3] = mu * qaws_patch_line(QAWS_PATCH_ARGS, QAWS_ZERO, cells, qaws_dual1_const(u), qaws_dual1_const(v1)).v;
}

/* c (N_i'(u) N_j(v) S_v x n + N_i(u) N_j'(v) n x S_u) at a node, for global
   control point (i, j); zero outside the node's spans. */
QAWS_INLINE qaws_vec3 qaws_patch_area_adjoint_cp(QAWS_PATCH_PARAMS, qaws_scalar u, qaws_scalar v, int i, int j, qaws_scalar c)
{
    qaws_scalar nu[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar nv[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int us = qaws_find_span(u_knots, u_degree, u_count, u);
    int vs = qaws_find_span(v_knots, v_degree, v_count, v);
    int li = i - (us - u_degree), lj = j - (vs - v_degree);
    qaws_vec3 su, sv, n, gu, gv, g;
    qaws_scalar len, a, b;
    g = qaws_v3_zero();
    if (li < 0 || li > u_degree || lj < 0 || lj > v_degree)
        return g;
    qaws_patch_rows(u_knots, u_degree, us, u, nu);
    qaws_patch_rows(v_knots, v_degree, vs, v, nv);
    su = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 1, 0);
    sv = qaws_patch_sum(cp, nu, nv, u_degree, v_degree, v_count, us, vs, 0, 1);
    n = qaws_v3_cross(su, sv);
    len = QAWS_SQRT(qaws_v3_dot(n, n));
    n = qaws_v3_scale(n, QAWS_SELECT(len > QAWS_ZERO, QAWS_ONE / QAWS_SELECT(len > QAWS_ZERO, len, QAWS_ONE), QAWS_ZERO));
    gu = qaws_v3_cross(sv, n);
    gv = qaws_v3_cross(n, su);
    a = nu[1 * (u_degree + 1) + li] * nv[0 * (v_degree + 1) + lj];
    b = nu[0 * (u_degree + 1) + li] * nv[1 * (v_degree + 1) + lj];
    g = qaws_v3_scale(qaws_v3_add(qaws_v3_scale(gu, a), qaws_v3_scale(gv, b)), c);
    return g;
}

/* c times the gradient of B(b; u) for global control point (i, j). */
QAWS_INLINE qaws_vec3 qaws_patch_line_adjoint_cp(QAWS_PATCH_PARAMS, int cells, qaws_scalar u, qaws_scalar b, int i, int j, qaws_scalar c)
{
    qaws_scalar v0 = v_knots[v_degree], v1 = v_knots[v_count];
    qaws_scalar hv = (v1 - v0) / QAWS_ITOF(cells);
    qaws_scalar xw[2];
    qaws_scalar base = v0, len;
    qaws_vec3 g = qaws_v3_zero();
    int cc, q, k;
    k = QAWS_FTOI((b - v0) / hv);
    k = k < 0 ? 0 : (k >= cells ? cells - 1 : k);
    for (cc = 0; cc < QAWS_CORE_MAX_CELLS; cc++) {
        if (cc < k) {
            for (q = 0; q < 8; q++) {
                qaws_gauss8(q, xw);
                g = qaws_v3_add(g, qaws_patch_area_adjoint_cp(QAWS_PATCH_ARGS, u, base + QAWS_LITERAL(0.5) * hv * (QAWS_ONE + xw[0]), i, j,
                    c * QAWS_LITERAL(0.5) * hv * xw[1]));
            }
        }
        base = base + QAWS_SELECT(cc < k, hv, QAWS_ZERO);
    }
    len = b - base;
    for (q = 0; q < 8; q++) {
        qaws_gauss8(q, xw);
        g = qaws_v3_add(g, qaws_patch_area_adjoint_cp(QAWS_PATCH_ARGS, u, base + QAWS_LITERAL(0.5) * len * (QAWS_ONE + xw[0]), i, j,
            c * QAWS_LITERAL(0.5) * len * xw[1]));
    }
    return g;
}

#endif /* QAWS_BSPLINE_SURFACE_SAMPLING_CORE_H */
