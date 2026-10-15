#ifndef QAWS_REPULSIVE_CORE_H
#define QAWS_REPULSIVE_CORE_H

#include "../qaws_platform.h"

/*
 * Tangent-point energies (qaws_repulsive.h), the pair terms as kernels for
 * every backend (C, HLSL, GLSL): closed forms of the direct evaluator's
 * terms, with first derivatives, and the metric's pair weights. The loops
 * over pairs belong to the caller (a thread per element I looping over
 * every J, then a gather per vertex over the elements holding it: no
 * atomics, deterministic).
 *
 * Curves (Repulsive Curves, Eq. 17), edges A = (a0, a1), B = (b0, b1), 3D
 * (a 2D curve has z = 0):
 *   e(A, B) = l_A^(1 - alpha) l_B / 4 sum_ij |e_A x (a_i - b_j)|^alpha |a_i - b_j|^-beta
 * for edges sharing no vertex (the caller skips the others).
 * Meshes (Repulsive Surfaces, Eq. 10), faces A, B (corners 3 x 3):
 *   e(A, B) = 1/4 |m_A|^(1 - p) |m_B| |m_A . (X_A - X_B)|^p |X_A - X_B|^-2p
 * for every A != B, m the face's cross product, X its barycenter.
 *
 * The gradient of the energy at a vertex is the sum, over the elements I
 * holding it, of I's block: sum over J of d/dI [e(I, J) + e(J, I)]
 * (qaws_repulsive_*_pair with which 0 and 1).
 *
 * The metric A = B + B0 (the fractional Sobolev inner product) applied to
 * a vertex function u: element I gathers sum_J kB (g_I - g_J) and
 * sum_J kB0 (u_I - u_J), g_I its derivative (coefficients from
 * qaws_repulsive_*_coefficients), u_I its average; each vertex then
 * collects c_Ia . accB_I + acc0_I / k from its elements.
 *
 * Status: the pair functions return 0, 1 (an element of no length or
 * area) or 2 (coincident points), and write nothing else on failure.
 */

/* x^y for x >= 0. When 2 y is a whole number (alpha / 2, beta / 2 and p
   usually are) by products and one square root: a GPU's pow is
   exp2(y log2 x), about 5e-6 off in single precision, and the curves'
   gradient is a small difference of large pair terms that magnifies it */
QAWS_INLINE qaws_scalar qaws_rp_pow(qaws_scalar x, qaws_scalar y)
{
    qaws_scalar y2 = y + y, r = QAWS_ONE;
    int n, i;
    if (!(QAWS_FLOOR(y2) == y2 && QAWS_FABS(y2) <= QAWS_LITERAL(64.0)))
        return QAWS_POW(x, y);
    n = QAWS_FTOI(QAWS_FABS(y2));
    for (i = 0; i < n / 2; i++) r = r * x;
    if (n - 2 * (n / 2) == 1)
        r = r * QAWS_SQRT(x);
    return y < QAWS_ZERO ? QAWS_ONE / r : r;
}

QAWS_INLINE qaws_scalar qaws_rp_dot3(qaws_scalar a[3], qaws_scalar b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

QAWS_INLINE void qaws_rp_cross3(qaws_scalar a[3], qaws_scalar b[3], QAWS_OUT qaws_scalar r[3])
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

/* e(A, B) for edges A = (p[0..2], p[3..5]), B = (q[0..2], q[3..5]); with
   which 0 the gradient in A's two vertices (g[0..5]), with which 1 in B's */
QAWS_INLINE int qaws_repulsive_curve_pair(qaws_scalar p[6], qaws_scalar q[6], qaws_scalar alpha, qaws_scalar beta, int which,
    QAWS_OUT qaws_scalar e_out[1], QAWS_OUT qaws_scalar g[6])
{
    qaws_scalar e[3];
    qaws_scalar f[3];
    qaws_scalar ge[3];
    qaws_scalar la = QAWS_ZERO, lb = QAWS_ZERO, W, S = QAWS_ZERO, a2 = QAWS_LITERAL(0.5) * alpha, b2 = QAWS_LITERAL(0.5) * beta;
    int i, j, c;
    for (c = 0; c < 3; c++) {
        e[c] = p[3 + c] - p[c];
        f[c] = q[3 + c] - q[c];
        la += e[c] * e[c];
        lb += f[c] * f[c];
        ge[c] = QAWS_ZERO;
    }
    for (c = 0; c < 6; c++) g[c] = QAWS_ZERO;
    e_out[0] = QAWS_ZERO;
    if (!(la > QAWS_ZERO && lb > QAWS_ZERO))
        return 1;
    W = QAWS_LITERAL(0.25) * qaws_rp_pow(la, QAWS_LITERAL(0.5) * (QAWS_ONE - alpha)) * QAWS_SQRT(lb);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++) {
            qaws_scalar d[3];
            qaws_scalar w[3];
            qaws_scalar dw[3];
            qaws_scalar we[3];
            qaws_scalar r2 = QAWS_ZERO, c2, k, kc, kr;
            for (c = 0; c < 3; c++) {
                d[c] = p[3 * i + c] - q[3 * j + c];
                r2 += d[c] * d[c];
            }
            if (!(r2 > QAWS_ZERO)) {
                for (c = 0; c < 6; c++) g[c] = QAWS_ZERO;   /* (nothing partial on failure) */
                return 2;
            }
            /* |e x d|^2 from the cross product itself: la r2 - (e.d)^2 loses
               every digit in single precision for neighbors nearly in line */
            qaws_rp_cross3(e, d, w);
            c2 = qaws_rp_dot3(w, w);
            k = qaws_rp_pow(c2, a2) * qaws_rp_pow(r2, -b2);
            S += k;
            /* dk = kc dc2 + kr dr2; dc2 = 2 (d x w) . de + 2 (w x e) . dd */
            kc = c2 > QAWS_ZERO ? a2 * qaws_rp_pow(c2, a2 - QAWS_ONE) * qaws_rp_pow(r2, -b2) : QAWS_SELECT(alpha == QAWS_LITERAL(2.0), qaws_rp_pow(r2, -b2), QAWS_ZERO);
            kr = -b2 * k / r2;
            qaws_rp_cross3(d, w, dw);
            qaws_rp_cross3(w, e, we);
            for (c = 0; c < 3; c++) {
                qaws_scalar gd = QAWS_LITERAL(2.0) * kc * we[c] + QAWS_LITERAL(2.0) * kr * d[c];
                if (which == 0) {
                    ge[c] += QAWS_LITERAL(2.0) * kc * dw[c];
                    g[3 * i + c] += W * gd;
                } else {
                    g[3 * j + c] -= W * gd;
                }
            }
        }
    e_out[0] = W * S;
    /* W's own: (1 - alpha) e / la on A's edge, f / lb on B's */
    for (c = 0; c < 3; c++) {
        qaws_scalar t = which == 0 ? W * ge[c] + (QAWS_ONE - alpha) * W * S * e[c] / la : W * S * f[c] / lb;
        g[3 + c] += t;
        g[c] -= t;
    }
    return 0;
}

/* d (u . m) / d corner of a face, m = e1 x e2 (e1 = x1 - x0, e2 = x2 - x0):
   e2 x u, u x e1, and minus their sum at corner 0 */
QAWS_INLINE void qaws_rp_dm(qaws_scalar e1[3], qaws_scalar e2[3], qaws_scalar u[3], QAWS_OUT qaws_scalar r[9])
{
    qaws_scalar a[3];
    qaws_scalar b[3];
    int c;
    qaws_rp_cross3(e2, u, a);
    qaws_rp_cross3(u, e1, b);
    for (c = 0; c < 3; c++) {
        r[c] = -a[c] - b[c];
        r[3 + c] = a[c];
        r[6 + c] = b[c];
    }
}

/* e(A, B) for faces A = p (corners 3 x 3), B = q; with which 0 the
   gradient in A's corners (g[0..8]), with which 1 in B's */
QAWS_INLINE int qaws_repulsive_mesh_pair(qaws_scalar p[9], qaws_scalar q[9], qaws_scalar pw, int which,
    QAWS_OUT qaws_scalar e_out[1], QAWS_OUT qaws_scalar g[9])
{
    qaws_scalar e1[3];
    qaws_scalar e2[3];
    qaws_scalar f1[3];
    qaws_scalar f2[3];
    qaws_scalar ma[3];
    qaws_scalar mb[3];
    qaws_scalar d[3];
    qaws_scalar u[3];
    qaws_scalar J[9];
    qaws_scalar a, b, s = QAWS_ZERO, r2 = QAWS_ZERO, e, ks, b2 = pw;
    int c, k;
    for (c = 0; c < 3; c++) {
        e1[c] = p[3 + c] - p[c];
        e2[c] = p[6 + c] - p[c];
        f1[c] = q[3 + c] - q[c];
        f2[c] = q[6 + c] - q[c];
    }
    for (c = 0; c < 9; c++) g[c] = QAWS_ZERO;
    e_out[0] = QAWS_ZERO;
    qaws_rp_cross3(e1, e2, ma);
    qaws_rp_cross3(f1, f2, mb);
    a = QAWS_SQRT(qaws_rp_dot3(ma, ma));
    b = QAWS_SQRT(qaws_rp_dot3(mb, mb));
    if (!(a > QAWS_ZERO && b > QAWS_ZERO))
        return 1;
    for (c = 0; c < 3; c++) {
        d[c] = (p[c] + p[3 + c] + p[6 + c] - q[c] - q[3 + c] - q[6 + c]) / QAWS_LITERAL(3.0);
        s += ma[c] * d[c];
        r2 += d[c] * d[c];
    }
    if (!(r2 > QAWS_ZERO))
        return 2;
    e = QAWS_LITERAL(0.25) * qaws_rp_pow(a, QAWS_ONE - pw) * b * qaws_rp_pow(QAWS_FABS(s), pw) * qaws_rp_pow(r2, -b2);
    e_out[0] = e;
    /* de = e ((1 - p) da / a + db / b - p dr2 / r2) + ks ds, ks = 1/4 a^(1 - p) b p |s|^(p - 1) sgn(s) r2^-p */
    ks = QAWS_LITERAL(0.25) * qaws_rp_pow(a, QAWS_ONE - pw) * b * pw * qaws_rp_pow(QAWS_FABS(s), pw - QAWS_ONE) * QAWS_SELECT(s < QAWS_ZERO, -QAWS_ONE, QAWS_ONE) * qaws_rp_pow(r2, -b2);
    if (which == 0) {
        for (c = 0; c < 3; c++) u[c] = e * (QAWS_ONE - pw) * ma[c] / (a * a) + ks * d[c];
        qaws_rp_dm(e1, e2, u, J);
        for (k = 0; k < 3; k++)
            for (c = 0; c < 3; c++) g[3 * k + c] = J[3 * k + c] + (ks * ma[c] - QAWS_LITERAL(2.0) * b2 * e * d[c] / r2) / QAWS_LITERAL(3.0);
    } else {
        for (c = 0; c < 3; c++) u[c] = e * mb[c] / (b * b);
        qaws_rp_dm(f1, f2, u, J);
        for (k = 0; k < 3; k++)
            for (c = 0; c < 3; c++) g[3 * k + c] = J[3 * k + c] - (ks * ma[c] - QAWS_LITERAL(2.0) * b2 * e * d[c] / r2) / QAWS_LITERAL(3.0);
    }
    return 0;
}

/* ---- the metric ---- */

/* an edge's derivative coefficients (g_I = u_1 c[3..5] + u_0 c[0..2]:
   (u_1 - u_0) / l T) and weight l */
QAWS_INLINE qaws_scalar qaws_repulsive_curve_coefficients(qaws_scalar p[6], QAWS_OUT qaws_scalar c[6])
{
    qaws_scalar e[3];
    qaws_scalar l2 = QAWS_ZERO;
    int k;
    for (k = 0; k < 3; k++) {
        e[k] = p[3 + k] - p[k];
        l2 += e[k] * e[k];
    }
    for (k = 0; k < 3; k++) {
        c[k] = l2 > QAWS_ZERO ? -e[k] / l2 : QAWS_ZERO;
        c[3 + k] = -c[k];
    }
    return QAWS_SQRT(l2);
}

/* a face's hat gradients N x e_opposite / (2 a) (c[3 k .. 3 k + 2] for
   corner k) and weight a */
QAWS_INLINE qaws_scalar qaws_repulsive_mesh_coefficients(qaws_scalar p[9], QAWS_OUT qaws_scalar c[9])
{
    qaws_scalar e1[3];
    qaws_scalar e2[3];
    qaws_scalar m[3];
    qaws_scalar n[3];
    qaws_scalar l;
    int k, q;
    for (q = 0; q < 3; q++) {
        e1[q] = p[3 + q] - p[q];
        e2[q] = p[6 + q] - p[q];
    }
    qaws_rp_cross3(e1, e2, m);
    l = QAWS_SQRT(qaws_rp_dot3(m, m));
    for (q = 0; q < 3; q++) n[q] = l > QAWS_ZERO ? m[q] / l : QAWS_ZERO;
    for (k = 0; k < 3; k++) {
        qaws_scalar o[3];
        qaws_scalar r[3];
        int i0 = k == 0 ? 1 : (k == 1 ? 2 : 0), i1 = k == 0 ? 2 : (k == 1 ? 0 : 1);
        for (q = 0; q < 3; q++) o[q] = p[3 * i1 + q] - p[3 * i0 + q];
        qaws_rp_cross3(n, o, r);
        for (q = 0; q < 3; q++) c[3 * k + q] = l > QAWS_ZERO ? r[q] / l : QAWS_ZERO;
    }
    return QAWS_LITERAL(0.5) * l;
}

/* the curves' weights for the pair of edges (A, B) sharing no vertex, as
   row A gathers them: kB = 2 w, w = 1/4 l_A l_B sum_ij r^-ex, and kB0 the
   low-order weights of both orders, |T x d|^2 / r^4 r^-ex; ex = 2 sigma + 1,
   sigma = (beta - 1) / alpha - 1. kw[0] = kB, kw[1] = kB0. */
QAWS_INLINE int qaws_repulsive_curve_metric_pair(qaws_scalar p[6], qaws_scalar q[6], qaws_scalar ex, QAWS_OUT qaws_scalar kw[2])
{
    qaws_scalar e[3];
    qaws_scalar f[3];
    qaws_scalar la = QAWS_ZERO, lb = QAWS_ZERO, w = QAWS_ZERO, w0 = QAWS_ZERO;
    int i, j, c;
    for (c = 0; c < 3; c++) {
        e[c] = p[3 + c] - p[c];
        f[c] = q[3 + c] - q[c];
        la += e[c] * e[c];
        lb += f[c] * f[c];
    }
    kw[0] = QAWS_ZERO;
    kw[1] = QAWS_ZERO;
    if (!(la > QAWS_ZERO && lb > QAWS_ZERO))
        return 1;
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++) {
            qaws_scalar d[3];
            qaws_scalar xa[3];
            qaws_scalar xb[3];
            qaws_scalar r2 = QAWS_ZERO, inv;
            for (c = 0; c < 3; c++) {
                d[c] = p[3 * i + c] - q[3 * j + c];
                r2 += d[c] * d[c];
            }
            if (!(r2 > QAWS_ZERO))
                return 2;
            inv = qaws_rp_pow(r2, QAWS_LITERAL(-0.5) * ex);
            /* |T x d|^2 by the cross product (stable in single precision) */
            qaws_rp_cross3(e, d, xa);
            qaws_rp_cross3(f, d, xb);
            w += inv;
            w0 += (qaws_rp_dot3(xa, xa) / la + qaws_rp_dot3(xb, xb) / lb) / (r2 * r2) * inv;
        }
    kw[0] = QAWS_LITERAL(0.5) * QAWS_SQRT(la * lb) * w;
    kw[1] = QAWS_LITERAL(0.25) * QAWS_SQRT(la * lb) * w0;
    return 0;
}

/* the meshes' weights for faces (A, B), A != B: w = a_A a_B r^-ex,
   kB = 2 w, kB0 = w ((N_A . d)^2 + (N_B . d)^2) / r^4; ex = 2 sigma + 2,
   sigma = 1 - 2 / p */
QAWS_INLINE int qaws_repulsive_mesh_metric_pair(qaws_scalar p[9], qaws_scalar q[9], qaws_scalar ex, QAWS_OUT qaws_scalar kw[2])
{
    qaws_scalar e1[3];
    qaws_scalar e2[3];
    qaws_scalar f1[3];
    qaws_scalar f2[3];
    qaws_scalar ma[3];
    qaws_scalar mb[3];
    qaws_scalar d[3];
    qaws_scalar a, b, r2 = QAWS_ZERO, na = QAWS_ZERO, nb = QAWS_ZERO, w;
    int c;
    for (c = 0; c < 3; c++) {
        e1[c] = p[3 + c] - p[c];
        e2[c] = p[6 + c] - p[c];
        f1[c] = q[3 + c] - q[c];
        f2[c] = q[6 + c] - q[c];
    }
    kw[0] = QAWS_ZERO;
    kw[1] = QAWS_ZERO;
    qaws_rp_cross3(e1, e2, ma);
    qaws_rp_cross3(f1, f2, mb);
    a = QAWS_SQRT(qaws_rp_dot3(ma, ma));
    b = QAWS_SQRT(qaws_rp_dot3(mb, mb));
    if (!(a > QAWS_ZERO && b > QAWS_ZERO))
        return 1;
    for (c = 0; c < 3; c++) {
        d[c] = (p[c] + p[3 + c] + p[6 + c] - q[c] - q[3 + c] - q[6 + c]) / QAWS_LITERAL(3.0);
        r2 += d[c] * d[c];
        na += ma[c] * d[c];
        nb += mb[c] * d[c];
    }
    if (!(r2 > QAWS_ZERO))
        return 2;
    na /= a;
    nb /= b;
    w = QAWS_LITERAL(0.25) * a * b * qaws_rp_pow(r2, QAWS_LITERAL(-0.5) * ex);
    kw[0] = QAWS_LITERAL(2.0) * w;
    kw[1] = w * (na * na + nb * nb) / (r2 * r2);
    return 0;
}

#endif /* QAWS_REPULSIVE_CORE_H */
