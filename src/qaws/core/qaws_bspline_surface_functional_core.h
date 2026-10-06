#ifndef QAWS_BSPLINE_SURFACE_FUNCTIONAL_CORE_H
#define QAWS_BSPLINE_SURFACE_FUNCTIONAL_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"
#include "qaws_dual_core.h"

/*
 * Integral functionals of B-spline surfaces, as kernels for every backend
 * (C, HLSL, GLSL). They match qaws_surface_functional_* of the runtime:
 * the domain is split into cells x cells cells of 4 x 4 Gauss points, and
 * every node is evaluated on the knot window of its own spans (cells need
 * not align with knots).
 *
 * Window of a node: local_cp holds the (u_degree + 1) x (v_degree + 1)
 * control points of its spans, local_cp[(i * (v_degree + 1) + j) * 3 + c]
 * for the global point (u_span - u_degree + i, v_span - v_degree + j);
 * local_u_knots / local_v_knots the knot windows around the spans (span
 * index u_degree / v_degree inside them); local_cp_dot the control point
 * tangents (zeros for none).
 *
 * Integrands f(S_u, S_v, S_uu, S_uv, S_vv), on second-order dual numbers
 * along the control point tangents:
 *   QAWS_CORE_SURFACE_AREA        |S_u x S_v|
 *   QAWS_CORE_SURFACE_THIN_PLATE  |S_uu|^2 + 2 |S_uv|^2 + |S_vv|^2
 *   QAWS_CORE_SURFACE_WILLMORE    H^2 |S_u x S_v|
 *
 * Kernels:
 *   qaws_surface_cell_node             parameters and weight of node q of a cell
 *   qaws_bspline_surface_node          weighted integrand at a node with its
 *                                      first and second rates (one thread per
 *                                      node, then a sum)
 *   qaws_bspline_surface_node_adjoint_cp
 *                                      c sum_k N^(k)_ij df/dy_k at a node for
 *                                      local control point (i, j): the
 *                                      per-element term of a gather over the
 *                                      nodes a control point supports
 */

#define QAWS_CORE_SURFACE_AREA       0
#define QAWS_CORE_SURFACE_THIN_PLATE 1
#define QAWS_CORE_SURFACE_WILLMORE   2

#define QAWS_CORE_SURFACE_WINDOW (QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3)

/* 4-point Gauss-Legendre rule on [-1, 1]: xw[0] node, xw[1] weight. */
QAWS_INLINE void qaws_gauss4(int k, QAWS_OUT qaws_scalar xw[2])
{
    qaws_scalar s = k < 2 ? -QAWS_ONE : QAWS_ONE;
    int m = k < 2 ? 1 - k : k - 2;
    xw[0] = s * QAWS_SELECT(m == 0, QAWS_LITERAL(0.3399810435848563), QAWS_LITERAL(0.8611363115940526));
    xw[1] = QAWS_SELECT(m == 0, QAWS_LITERAL(0.6521451548625461), QAWS_LITERAL(0.3478548451374538));
}

/* Node q (0..15) of cell (ci, cj) of a cells x cells grid over
   [u0, u1] x [v0, v1]: node = (u, v, weight). */
QAWS_INLINE void qaws_surface_cell_node(
    qaws_scalar u0, qaws_scalar u1, qaws_scalar v0, qaws_scalar v1,
    qaws_scalar cell_u, qaws_scalar cell_v, qaws_scalar cells, int q,
    QAWS_OUT qaws_scalar node[3])
{
    qaws_scalar xa[2];
    qaws_scalar xb[2];
    qaws_scalar du = (u1 - u0) / cells, dv = (v1 - v0) / cells;
    qaws_gauss4(q / 4, xa);
    qaws_gauss4(q - 4 * (q / 4), xb);
    node[0] = u0 + du * (cell_u + QAWS_LITERAL(0.5) + QAWS_LITERAL(0.5) * xa[0]);
    node[1] = v0 + dv * (cell_v + QAWS_LITERAL(0.5) + QAWS_LITERAL(0.5) * xb[0]);
    node[2] = QAWS_LITERAL(0.25) * du * dv * xa[1] * xb[1];
}

/* f on dual jets y = (S_u, S_v, S_uu, S_uv, S_vv). */
QAWS_INLINE qaws_dual1 qaws_core_surface_integrand(int id, qaws_dual3 su, qaws_dual3 sv, qaws_dual3 suu, qaws_dual3 suv, qaws_dual3 svv)
{
    qaws_dual3 nr = qaws_dual3_cross(su, sv);
    qaws_dual1 nl = qaws_dual3_length(nr);
    qaws_dual1 r = nl;
    if (id == QAWS_CORE_SURFACE_THIN_PLATE) {
        qaws_dual1 a = qaws_dual3_dot(suu, suu);
        qaws_dual1 b = qaws_dual3_dot(suv, suv);
        qaws_dual1 c = qaws_dual3_dot(svv, svv);
        r = qaws_dual1_add(qaws_dual1_add(a, qaws_dual1_make(QAWS_LITERAL(2.0) * b.v, QAWS_LITERAL(2.0) * b.t, QAWS_LITERAL(2.0) * b.tt)), c);
    }
    if (id == QAWS_CORE_SURFACE_WILLMORE) {
        qaws_dual3 N = qaws_dual3_div(nr, nl);
        qaws_dual1 E = qaws_dual3_dot(su, su);
        qaws_dual1 F = qaws_dual3_dot(su, sv);
        qaws_dual1 G = qaws_dual3_dot(sv, sv);
        qaws_dual1 L = qaws_dual3_dot(suu, N);
        qaws_dual1 M = qaws_dual3_dot(suv, N);
        qaws_dual1 Nn = qaws_dual3_dot(svv, N);
        qaws_dual1 det = qaws_dual1_sub(qaws_dual1_mul(E, G), qaws_dual1_mul(F, F));
        qaws_dual1 num = qaws_dual1_sub(qaws_dual1_add(qaws_dual1_mul(E, Nn), qaws_dual1_mul(G, L)),
            qaws_dual1_mul(qaws_dual1_make(QAWS_LITERAL(2.0) * F.v, QAWS_LITERAL(2.0) * F.t, QAWS_LITERAL(2.0) * F.tt), M));
        qaws_dual1 H = qaws_dual1_div(num, qaws_dual1_make(QAWS_LITERAL(2.0) * det.v, QAWS_LITERAL(2.0) * det.t, QAWS_LITERAL(2.0) * det.tt));
        r = qaws_dual1_mul(qaws_dual1_mul(H, H), nl);
    }
    return r;
}

/* Basis rows 0..2 of a window, rows above the degree zero. */
QAWS_INLINE void qaws_surface_rows(
    qaws_scalar local_knots[QAWS_CORE_MAX_POINTS * 2],
    int degree,
    qaws_scalar t,
    QAWS_OUT qaws_scalar ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS])
{
    int max_k = degree < 2 ? degree : 2;
    int k, j;
    qaws_bspline_basis_derivs(local_knots, degree, degree, t, max_k, ders);
    for (k = max_k + 1; k <= 2; k++)
        for (j = 0; j <= degree; j++)
            ders[k * (degree + 1) + j] = QAWS_ZERO;
}

/* Jet entries (S_u, S_v, S_uu, S_uv, S_vv) of the window and of its
   tangents: y[k] for the value, y_dot[k] along local_cp_dot. */
QAWS_INLINE void qaws_bspline_surface_jet_rows(
    qaws_scalar local_cp[QAWS_CORE_SURFACE_WINDOW],
    qaws_scalar local_cp_dot[QAWS_CORE_SURFACE_WINDOW],
    qaws_scalar local_u_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_v_knots[QAWS_CORE_MAX_POINTS * 2],
    int u_degree,
    int v_degree,
    qaws_scalar u,
    qaws_scalar v,
    QAWS_OUT qaws_vec3 y[5],
    QAWS_OUT qaws_vec3 y_dot[5])
{
    qaws_scalar nu[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar nv[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    int ku_of[5];
    int kv_of[5];
    int k, i, j, uo = u_degree + 1, vo = v_degree + 1, at;
    qaws_scalar w;
    ku_of[0] = 1; kv_of[0] = 0;
    ku_of[1] = 0; kv_of[1] = 1;
    ku_of[2] = 2; kv_of[2] = 0;
    ku_of[3] = 1; kv_of[3] = 1;
    ku_of[4] = 0; kv_of[4] = 2;
    qaws_surface_rows(local_u_knots, u_degree, u, nu);
    qaws_surface_rows(local_v_knots, v_degree, v, nv);
    for (k = 0; k < 5; k++) {
        y[k] = qaws_v3_zero();
        y_dot[k] = qaws_v3_zero();
        for (i = 0; i < uo; i++)
            for (j = 0; j < vo; j++) {
                w = nu[ku_of[k] * uo + i] * nv[kv_of[k] * vo + j];
                at = (i * vo + j) * 3;
                y[k] = qaws_v3_axpy(y[k], qaws_v3(local_cp[at], local_cp[at + 1], local_cp[at + 2]), w);
                y_dot[k] = qaws_v3_axpy(y_dot[k], qaws_v3(local_cp_dot[at], local_cp_dot[at + 1], local_cp_dot[at + 2]), w);
            }
    }
}

/* weight f at (u, v) with its first and second rates along local_cp_dot. */
QAWS_INLINE qaws_dual1 qaws_bspline_surface_node(
    int id,
    qaws_scalar local_cp[QAWS_CORE_SURFACE_WINDOW],
    qaws_scalar local_cp_dot[QAWS_CORE_SURFACE_WINDOW],
    qaws_scalar local_u_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_v_knots[QAWS_CORE_MAX_POINTS * 2],
    int u_degree,
    int v_degree,
    qaws_scalar u,
    qaws_scalar v,
    qaws_scalar weight)
{
    qaws_vec3 y[5];
    qaws_vec3 yd[5];
    qaws_dual1 f;
    qaws_bspline_surface_jet_rows(local_cp, local_cp_dot, local_u_knots, local_v_knots, u_degree, v_degree, u, v, y, yd);
    f = qaws_core_surface_integrand(id,
        qaws_dual3_make(y[0], yd[0], qaws_v3_zero()), qaws_dual3_make(y[1], yd[1], qaws_v3_zero()),
        qaws_dual3_make(y[2], yd[2], qaws_v3_zero()), qaws_dual3_make(y[3], yd[3], qaws_v3_zero()),
        qaws_dual3_make(y[4], yd[4], qaws_v3_zero()));
    return qaws_dual1_make(weight * f.v, weight * f.t, weight * f.tt);
}

/* c sum_k N_i^(ku_k)(u) N_j^(kv_k)(v) df/dy_k at (u, v) for local control
   point (i, j). */
QAWS_INLINE qaws_vec3 qaws_bspline_surface_node_adjoint_cp(
    int id,
    qaws_scalar local_cp[QAWS_CORE_SURFACE_WINDOW],
    qaws_scalar local_u_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar local_v_knots[QAWS_CORE_MAX_POINTS * 2],
    int u_degree,
    int v_degree,
    qaws_scalar u,
    qaws_scalar v,
    int i,
    int j,
    qaws_scalar c)
{
    qaws_scalar nu[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar nv[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_vec3 y[5];
    qaws_vec3 yd[5];
    qaws_dual3 s[5];
    qaws_vec3 g = qaws_v3_zero();
    qaws_vec3 e;
    qaws_scalar n, df;
    int ku_of[5];
    int kv_of[5];
    int k, m, comp, uo = u_degree + 1, vo = v_degree + 1;
    ku_of[0] = 1; kv_of[0] = 0;
    ku_of[1] = 0; kv_of[1] = 1;
    ku_of[2] = 2; kv_of[2] = 0;
    ku_of[3] = 1; kv_of[3] = 1;
    ku_of[4] = 0; kv_of[4] = 2;
    qaws_bspline_surface_jet_rows(local_cp, local_cp, local_u_knots, local_v_knots, u_degree, v_degree, u, v, y, yd);
    qaws_surface_rows(local_u_knots, u_degree, u, nu);
    qaws_surface_rows(local_v_knots, v_degree, v, nv);
    for (k = 0; k < 5; k++) {
        n = nu[ku_of[k] * uo + i] * nv[kv_of[k] * vo + j];
        for (comp = 0; comp < 3; comp++) {
            e = qaws_v3(comp == 0 ? QAWS_ONE : QAWS_ZERO, comp == 1 ? QAWS_ONE : QAWS_ZERO, comp == 2 ? QAWS_ONE : QAWS_ZERO);
            for (m = 0; m < 5; m++)
                s[m] = qaws_dual3_make(y[m], qaws_v3_scale(e, m == k ? QAWS_ONE : QAWS_ZERO), qaws_v3_zero());
            df = qaws_core_surface_integrand(id, s[0], s[1], s[2], s[3], s[4]).t;
            if (comp == 0) g.x += c * n * df;
            if (comp == 1) g.y += c * n * df;
            if (comp == 2) g.z += c * n * df;
        }
    }
    return g;
}

#endif /* QAWS_BSPLINE_SURFACE_FUNCTIONAL_CORE_H */
