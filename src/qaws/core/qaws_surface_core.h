#ifndef QAWS_SURFACE_CORE_H
#define QAWS_SURFACE_CORE_H

#include "../qaws_platform.h"
#include "../qaws_core_types.h"
#include "qaws_bspline_basis_core.h"

/* ===================================================================
 * Surface eval flag bits
 * =================================================================== */

/* the values of QAWS_SURFACE_EVAL_* (qaws_surface_types.h): one flag word for the runtime and the cores */
#define QAWS_SURFACE_CORE_POSITION  0x1
#define QAWS_SURFACE_CORE_DU        0x2
#define QAWS_SURFACE_CORE_DV        0x4
#define QAWS_SURFACE_CORE_NORMAL    0x8
#define QAWS_SURFACE_CORE_DUU       0x10
#define QAWS_SURFACE_CORE_DUV       0x20
#define QAWS_SURFACE_CORE_DVV       0x40

/* ===================================================================
 * Surface evaluation result type
 * =================================================================== */

#ifndef QAWS_SURFACE_EVAL_3D_DEFINED
#define QAWS_SURFACE_EVAL_3D_DEFINED

QAWS_TYPE_DEF struct {
    qaws_vec3 position;
    qaws_vec3 du, dv;
    qaws_vec3 duu, dvv, duv;
    qaws_vec3 normal;
} qaws_surface_eval_3d;

#endif /* QAWS_SURFACE_EVAL_3D_DEFINED */

/* ===================================================================
 * Tensor-product Bezier surface evaluation (3D)
 *
 * cp: flat array of 3D control points, row-major:
 *     cp[(i * v_count + j) * 3 + k] for row i, column j, component k
 * u_degree, v_degree: polynomial degrees in u and v
 * v_count: number of control points per u-row (= v_degree + 1)
 * u, v: parameter values in [0, 1]
 * eval_flags: bitmask of QAWS_SURFACE_CORE_* flags
 * =================================================================== */

QAWS_INLINE qaws_surface_eval_3d qaws_surface_bezier_eval_3d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3],
    int u_degree,
    int v_degree,
    int v_count,
    qaws_scalar u,
    qaws_scalar v,
    int eval_flags)
{
    qaws_surface_eval_3d result;
    qaws_scalar work[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar col[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar dv_col[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar dvv_col[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar dcol[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar ddcol[QAWS_CORE_MAX_POINTS * 3];
    int u_count = u_degree + 1;
    int row_stride = v_count * 3;
    int i, r, k;
    qaws_scalar s, nx, ny, nz, len;

    /* Zero the result */
    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO; result.position.z = QAWS_ZERO;
    result.du.x = QAWS_ZERO; result.du.y = QAWS_ZERO; result.du.z = QAWS_ZERO;
    result.dv.x = QAWS_ZERO; result.dv.y = QAWS_ZERO; result.dv.z = QAWS_ZERO;
    result.duu.x = QAWS_ZERO; result.duu.y = QAWS_ZERO; result.duu.z = QAWS_ZERO;
    result.dvv.x = QAWS_ZERO; result.dvv.y = QAWS_ZERO; result.dvv.z = QAWS_ZERO;
    result.duv.x = QAWS_ZERO; result.duv.y = QAWS_ZERO; result.duv.z = QAWS_ZERO;
    result.normal.x = QAWS_ZERO; result.normal.y = QAWS_ZERO; result.normal.z = QAWS_ZERO;

    s = QAWS_ONE - v;

    /* Step 1: For each u-row, evaluate at v using de Casteljau to get a column */
    for (i = 0; i < u_count; i++) {
        for (k = 0; k <= v_degree; k++) {
            work[k * 3 + 0] = cp[i * row_stride + k * 3 + 0];
            work[k * 3 + 1] = cp[i * row_stride + k * 3 + 1];
            work[k * 3 + 2] = cp[i * row_stride + k * 3 + 2];
        }
        for (r = 1; r <= v_degree; r++) {
            for (k = 0; k <= v_degree - r; k++) {
                work[k * 3 + 0] = s * work[k * 3 + 0] + v * work[(k + 1) * 3 + 0];
                work[k * 3 + 1] = s * work[k * 3 + 1] + v * work[(k + 1) * 3 + 1];
                work[k * 3 + 2] = s * work[k * 3 + 2] + v * work[(k + 1) * 3 + 2];
            }
        }
        col[i * 3 + 0] = work[0];
        col[i * 3 + 1] = work[1];
        col[i * 3 + 2] = work[2];
    }

    /* Position: evaluate column at u */
    if (eval_flags & QAWS_SURFACE_CORE_POSITION) {
        for (k = 0; k < u_count; k++) {
            work[k * 3 + 0] = col[k * 3 + 0];
            work[k * 3 + 1] = col[k * 3 + 1];
            work[k * 3 + 2] = col[k * 3 + 2];
        }
        s = QAWS_ONE - u;
        for (r = 1; r <= u_degree; r++) {
            for (k = 0; k <= u_degree - r; k++) {
                work[k * 3 + 0] = s * work[k * 3 + 0] + u * work[(k + 1) * 3 + 0];
                work[k * 3 + 1] = s * work[k * 3 + 1] + u * work[(k + 1) * 3 + 1];
                work[k * 3 + 2] = s * work[k * 3 + 2] + u * work[(k + 1) * 3 + 2];
            }
        }
        result.position.x = work[0];
        result.position.y = work[1];
        result.position.z = work[2];
    }

    /* du: derivative control points in u, then evaluate at u */
    if (eval_flags & (QAWS_SURFACE_CORE_DU | QAWS_SURFACE_CORE_DUU | QAWS_SURFACE_CORE_DUV | QAWS_SURFACE_CORE_NORMAL)) {
        qaws_scalar deg_u = (qaws_scalar)u_degree;

        /* D(col)/du: first-difference the column */
        for (k = 0; k < u_degree; k++) {
            dcol[k * 3 + 0] = deg_u * (col[(k + 1) * 3 + 0] - col[k * 3 + 0]);
            dcol[k * 3 + 1] = deg_u * (col[(k + 1) * 3 + 1] - col[k * 3 + 1]);
            dcol[k * 3 + 2] = deg_u * (col[(k + 1) * 3 + 2] - col[k * 3 + 2]);
        }

        /* Evaluate dcol at u */
        for (k = 0; k < u_degree; k++) {
            work[k * 3 + 0] = dcol[k * 3 + 0];
            work[k * 3 + 1] = dcol[k * 3 + 1];
            work[k * 3 + 2] = dcol[k * 3 + 2];
        }
        s = QAWS_ONE - u;
        for (r = 1; r < u_degree; r++) {
            for (k = 0; k < u_degree - r; k++) {
                work[k * 3 + 0] = s * work[k * 3 + 0] + u * work[(k + 1) * 3 + 0];
                work[k * 3 + 1] = s * work[k * 3 + 1] + u * work[(k + 1) * 3 + 1];
                work[k * 3 + 2] = s * work[k * 3 + 2] + u * work[(k + 1) * 3 + 2];
            }
        }
        result.du.x = work[0];
        result.du.y = work[1];
        result.du.z = work[2];

        /* duu: second derivative in u */
        if ((eval_flags & QAWS_SURFACE_CORE_DUU) && u_degree >= 2) {
            qaws_scalar deg_u1 = (qaws_scalar)(u_degree - 1);
            for (k = 0; k < u_degree - 1; k++) {
                ddcol[k * 3 + 0] = deg_u1 * (dcol[(k + 1) * 3 + 0] - dcol[k * 3 + 0]);
                ddcol[k * 3 + 1] = deg_u1 * (dcol[(k + 1) * 3 + 1] - dcol[k * 3 + 1]);
                ddcol[k * 3 + 2] = deg_u1 * (dcol[(k + 1) * 3 + 2] - dcol[k * 3 + 2]);
            }
            for (k = 0; k < u_degree - 1; k++) {
                work[k * 3 + 0] = ddcol[k * 3 + 0];
                work[k * 3 + 1] = ddcol[k * 3 + 1];
                work[k * 3 + 2] = ddcol[k * 3 + 2];
            }
            s = QAWS_ONE - u;
            for (r = 1; r < u_degree - 1; r++) {
                for (k = 0; k < u_degree - 1 - r; k++) {
                    work[k * 3 + 0] = s * work[k * 3 + 0] + u * work[(k + 1) * 3 + 0];
                    work[k * 3 + 1] = s * work[k * 3 + 1] + u * work[(k + 1) * 3 + 1];
                    work[k * 3 + 2] = s * work[k * 3 + 2] + u * work[(k + 1) * 3 + 2];
                }
            }
            result.duu.x = work[0];
            result.duu.y = work[1];
            result.duu.z = work[2];
        }
    }

    /* dv: differentiate each row in v, evaluate at v, then evaluate column at u */
    if (eval_flags & (QAWS_SURFACE_CORE_DV | QAWS_SURFACE_CORE_DVV | QAWS_SURFACE_CORE_DUV | QAWS_SURFACE_CORE_NORMAL)) {
        qaws_scalar deg_v = (qaws_scalar)v_degree;

        for (i = 0; i < u_count; i++) {
            /* Derivative CPs for this row in v */
            qaws_scalar dv_row[QAWS_CORE_MAX_POINTS * 3];
            for (k = 0; k < v_degree; k++) {
                dv_row[k * 3 + 0] = deg_v * (cp[i * row_stride + (k + 1) * 3 + 0] - cp[i * row_stride + k * 3 + 0]);
                dv_row[k * 3 + 1] = deg_v * (cp[i * row_stride + (k + 1) * 3 + 1] - cp[i * row_stride + k * 3 + 1]);
                dv_row[k * 3 + 2] = deg_v * (cp[i * row_stride + (k + 1) * 3 + 2] - cp[i * row_stride + k * 3 + 2]);
            }
            /* Evaluate dv_row at v using de Casteljau */
            for (k = 0; k < v_degree; k++) {
                work[k * 3 + 0] = dv_row[k * 3 + 0];
                work[k * 3 + 1] = dv_row[k * 3 + 1];
                work[k * 3 + 2] = dv_row[k * 3 + 2];
            }
            s = QAWS_ONE - v;
            for (r = 1; r < v_degree; r++) {
                for (k = 0; k < v_degree - r; k++) {
                    work[k * 3 + 0] = s * work[k * 3 + 0] + v * work[(k + 1) * 3 + 0];
                    work[k * 3 + 1] = s * work[k * 3 + 1] + v * work[(k + 1) * 3 + 1];
                    work[k * 3 + 2] = s * work[k * 3 + 2] + v * work[(k + 1) * 3 + 2];
                }
            }
            dv_col[i * 3 + 0] = work[0];
            dv_col[i * 3 + 1] = work[1];
            dv_col[i * 3 + 2] = work[2];
        }

        /* Evaluate dv_col at u */
        for (k = 0; k < u_count; k++) {
            work[k * 3 + 0] = dv_col[k * 3 + 0];
            work[k * 3 + 1] = dv_col[k * 3 + 1];
            work[k * 3 + 2] = dv_col[k * 3 + 2];
        }
        s = QAWS_ONE - u;
        for (r = 1; r <= u_degree; r++) {
            for (k = 0; k <= u_degree - r; k++) {
                work[k * 3 + 0] = s * work[k * 3 + 0] + u * work[(k + 1) * 3 + 0];
                work[k * 3 + 1] = s * work[k * 3 + 1] + u * work[(k + 1) * 3 + 1];
                work[k * 3 + 2] = s * work[k * 3 + 2] + u * work[(k + 1) * 3 + 2];
            }
        }
        result.dv.x = work[0];
        result.dv.y = work[1];
        result.dv.z = work[2];

        /* dvv: second derivative in v */
        if ((eval_flags & QAWS_SURFACE_CORE_DVV) && v_degree >= 2) {
            qaws_scalar deg_v1 = (qaws_scalar)(v_degree - 1);

            for (i = 0; i < u_count; i++) {
                qaws_scalar dv_row[QAWS_CORE_MAX_POINTS * 3];
                qaws_scalar dvv_row[QAWS_CORE_MAX_POINTS * 3];

                /* First v-derivative CPs */
                for (k = 0; k < v_degree; k++) {
                    dv_row[k * 3 + 0] = deg_v * (cp[i * row_stride + (k + 1) * 3 + 0] - cp[i * row_stride + k * 3 + 0]);
                    dv_row[k * 3 + 1] = deg_v * (cp[i * row_stride + (k + 1) * 3 + 1] - cp[i * row_stride + k * 3 + 1]);
                    dv_row[k * 3 + 2] = deg_v * (cp[i * row_stride + (k + 1) * 3 + 2] - cp[i * row_stride + k * 3 + 2]);
                }
                /* Second v-derivative CPs */
                for (k = 0; k < v_degree - 1; k++) {
                    dvv_row[k * 3 + 0] = deg_v1 * (dv_row[(k + 1) * 3 + 0] - dv_row[k * 3 + 0]);
                    dvv_row[k * 3 + 1] = deg_v1 * (dv_row[(k + 1) * 3 + 1] - dv_row[k * 3 + 1]);
                    dvv_row[k * 3 + 2] = deg_v1 * (dv_row[(k + 1) * 3 + 2] - dv_row[k * 3 + 2]);
                }
                /* Evaluate dvv_row at v */
                for (k = 0; k < v_degree - 1; k++) {
                    work[k * 3 + 0] = dvv_row[k * 3 + 0];
                    work[k * 3 + 1] = dvv_row[k * 3 + 1];
                    work[k * 3 + 2] = dvv_row[k * 3 + 2];
                }
                s = QAWS_ONE - v;
                for (r = 1; r < v_degree - 1; r++) {
                    for (k = 0; k < v_degree - 1 - r; k++) {
                        work[k * 3 + 0] = s * work[k * 3 + 0] + v * work[(k + 1) * 3 + 0];
                        work[k * 3 + 1] = s * work[k * 3 + 1] + v * work[(k + 1) * 3 + 1];
                        work[k * 3 + 2] = s * work[k * 3 + 2] + v * work[(k + 1) * 3 + 2];
                    }
                }
                dvv_col[i * 3 + 0] = work[0];
                dvv_col[i * 3 + 1] = work[1];
                dvv_col[i * 3 + 2] = work[2];
            }

            /* Evaluate dvv_col at u */
            for (k = 0; k < u_count; k++) {
                work[k * 3 + 0] = dvv_col[k * 3 + 0];
                work[k * 3 + 1] = dvv_col[k * 3 + 1];
                work[k * 3 + 2] = dvv_col[k * 3 + 2];
            }
            s = QAWS_ONE - u;
            for (r = 1; r <= u_degree; r++) {
                for (k = 0; k <= u_degree - r; k++) {
                    work[k * 3 + 0] = s * work[k * 3 + 0] + u * work[(k + 1) * 3 + 0];
                    work[k * 3 + 1] = s * work[k * 3 + 1] + u * work[(k + 1) * 3 + 1];
                    work[k * 3 + 2] = s * work[k * 3 + 2] + u * work[(k + 1) * 3 + 2];
                }
            }
            result.dvv.x = work[0];
            result.dvv.y = work[1];
            result.dvv.z = work[2];
        }

        /* duv: differentiate the dv_col in u, then evaluate at u */
        if (eval_flags & QAWS_SURFACE_CORE_DUV) {
            qaws_scalar duv_col[QAWS_CORE_MAX_POINTS * 3];
            qaws_scalar deg_u = (qaws_scalar)u_degree;

            for (k = 0; k < u_degree; k++) {
                duv_col[k * 3 + 0] = deg_u * (dv_col[(k + 1) * 3 + 0] - dv_col[k * 3 + 0]);
                duv_col[k * 3 + 1] = deg_u * (dv_col[(k + 1) * 3 + 1] - dv_col[k * 3 + 1]);
                duv_col[k * 3 + 2] = deg_u * (dv_col[(k + 1) * 3 + 2] - dv_col[k * 3 + 2]);
            }

            for (k = 0; k < u_degree; k++) {
                work[k * 3 + 0] = duv_col[k * 3 + 0];
                work[k * 3 + 1] = duv_col[k * 3 + 1];
                work[k * 3 + 2] = duv_col[k * 3 + 2];
            }
            s = QAWS_ONE - u;
            for (r = 1; r < u_degree; r++) {
                for (k = 0; k < u_degree - r; k++) {
                    work[k * 3 + 0] = s * work[k * 3 + 0] + u * work[(k + 1) * 3 + 0];
                    work[k * 3 + 1] = s * work[k * 3 + 1] + u * work[(k + 1) * 3 + 1];
                    work[k * 3 + 2] = s * work[k * 3 + 2] + u * work[(k + 1) * 3 + 2];
                }
            }
            result.duv.x = work[0];
            result.duv.y = work[1];
            result.duv.z = work[2];
        }
    }

    /* Normal: cross product of du x dv, normalized */
    if (eval_flags & QAWS_SURFACE_CORE_NORMAL) {
        nx = result.du.y * result.dv.z - result.du.z * result.dv.y;
        ny = result.du.z * result.dv.x - result.du.x * result.dv.z;
        nz = result.du.x * result.dv.y - result.du.y * result.dv.x;
        len = QAWS_SQRT(nx * nx + ny * ny + nz * nz);
        if (len > QAWS_LITERAL(1e-12)) {
            result.normal.x = nx / len;
            result.normal.y = ny / len;
            result.normal.z = nz / len;
        } else {
            result.normal.x = QAWS_ZERO;
            result.normal.y = QAWS_ZERO;
            result.normal.z = QAWS_ONE;
        }
    }

    return result;
}

/* ===================================================================
 * Tensor-product B-spline surface evaluation (3D)
 *
 * Uses basis function derivatives from qaws_bspline_basis_core.h.
 *
 * cp: flat array of 3D control points, row-major:
 *     cp[(i * v_count + j) * 3 + k]
 * u_knots, v_knots: knot vectors
 * u_degree, v_degree: polynomial degrees
 * u_count, v_count: number of control points in u and v
 * u_span, v_span: precomputed knot span indices
 * u, v: parameter values
 * eval_flags: bitmask of QAWS_SURFACE_CORE_* flags
 * =================================================================== */

QAWS_INLINE qaws_surface_eval_3d qaws_surface_bspline_eval_3d(
    qaws_scalar cp[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3],
    qaws_scalar u_knots[QAWS_CORE_MAX_POINTS * 2],
    qaws_scalar v_knots[QAWS_CORE_MAX_POINTS * 2],
    int u_degree,
    int v_degree,
    int u_count,
    int v_count,
    int u_span,
    int v_span,
    qaws_scalar u,
    qaws_scalar v,
    int eval_flags)
{
    qaws_surface_eval_3d result;
    int u_order = u_degree + 1;
    int v_order = v_degree + 1;
    int max_u_deriv = 0;
    int max_v_deriv = 0;
    qaws_scalar u_ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar v_ders[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar S[3][3][3]; /* S[ku][kv][xyz] */
    int ku, kv, i, j;
    qaws_scalar Nu, Nv, w, nx, ny, nz, len;
    int ci, cj, cp_idx;

    /* Zero the result */
    result.position.x = QAWS_ZERO; result.position.y = QAWS_ZERO; result.position.z = QAWS_ZERO;
    result.du.x = QAWS_ZERO; result.du.y = QAWS_ZERO; result.du.z = QAWS_ZERO;
    result.dv.x = QAWS_ZERO; result.dv.y = QAWS_ZERO; result.dv.z = QAWS_ZERO;
    result.duu.x = QAWS_ZERO; result.duu.y = QAWS_ZERO; result.duu.z = QAWS_ZERO;
    result.dvv.x = QAWS_ZERO; result.dvv.y = QAWS_ZERO; result.dvv.z = QAWS_ZERO;
    result.duv.x = QAWS_ZERO; result.duv.y = QAWS_ZERO; result.duv.z = QAWS_ZERO;
    result.normal.x = QAWS_ZERO; result.normal.y = QAWS_ZERO; result.normal.z = QAWS_ZERO;

    /* Determine derivative orders needed */
    if (eval_flags & QAWS_SURFACE_CORE_DU) max_u_deriv = 1;
    if (eval_flags & QAWS_SURFACE_CORE_DUU) max_u_deriv = 2;
    if (eval_flags & QAWS_SURFACE_CORE_DV) max_v_deriv = 1;
    if (eval_flags & QAWS_SURFACE_CORE_DVV) max_v_deriv = 2;
    if (eval_flags & QAWS_SURFACE_CORE_DUV) {
        if (max_u_deriv < 1) max_u_deriv = 1;
        if (max_v_deriv < 1) max_v_deriv = 1;
    }
    if (eval_flags & QAWS_SURFACE_CORE_NORMAL) {
        if (max_u_deriv < 1) max_u_deriv = 1;
        if (max_v_deriv < 1) max_v_deriv = 1;
    }
    if (max_u_deriv > u_degree) max_u_deriv = u_degree;
    if (max_v_deriv > v_degree) max_v_deriv = v_degree;

    /* Compute basis function derivatives */
    qaws_bspline_basis_derivs(u_knots, u_degree, u_span, u, max_u_deriv, u_ders);
    qaws_bspline_basis_derivs(v_knots, v_degree, v_span, v, max_v_deriv, v_ders);

    /* Zero S accumulator */
    for (ku = 0; ku < 3; ku++)
        for (kv = 0; kv < 3; kv++) {
            S[ku][kv][0] = QAWS_ZERO;
            S[ku][kv][1] = QAWS_ZERO;
            S[ku][kv][2] = QAWS_ZERO;
        }

    /* Compute S[ku][kv] = sum_i sum_j N_i^(ku)(u) * N_j^(kv)(v) * P[i][j] */
    for (ku = 0; ku <= max_u_deriv; ku++) {
        for (kv = 0; kv <= max_v_deriv; kv++) {
            for (i = 0; i < u_order; i++) {
                ci = u_span - u_degree + i;
                Nu = u_ders[ku * u_order + i];
                for (j = 0; j < v_order; j++) {
                    cj = v_span - v_degree + j;
                    Nv = v_ders[kv * v_order + j];
                    w = Nu * Nv;
                    cp_idx = (ci * v_count + cj) * 3;
                    S[ku][kv][0] += w * cp[cp_idx + 0];
                    S[ku][kv][1] += w * cp[cp_idx + 1];
                    S[ku][kv][2] += w * cp[cp_idx + 2];
                }
            }
        }
    }

    /* Position */
    if (eval_flags & QAWS_SURFACE_CORE_POSITION) {
        result.position.x = S[0][0][0];
        result.position.y = S[0][0][1];
        result.position.z = S[0][0][2];
    }

    /* du */
    if (eval_flags & (QAWS_SURFACE_CORE_DU | QAWS_SURFACE_CORE_NORMAL)) {
        result.du.x = S[1][0][0];
        result.du.y = S[1][0][1];
        result.du.z = S[1][0][2];
    }

    /* dv */
    if (eval_flags & (QAWS_SURFACE_CORE_DV | QAWS_SURFACE_CORE_NORMAL)) {
        result.dv.x = S[0][1][0];
        result.dv.y = S[0][1][1];
        result.dv.z = S[0][1][2];
    }

    /* duu */
    if ((eval_flags & QAWS_SURFACE_CORE_DUU) && max_u_deriv >= 2) {
        result.duu.x = S[2][0][0];
        result.duu.y = S[2][0][1];
        result.duu.z = S[2][0][2];
    }

    /* dvv */
    if ((eval_flags & QAWS_SURFACE_CORE_DVV) && max_v_deriv >= 2) {
        result.dvv.x = S[0][2][0];
        result.dvv.y = S[0][2][1];
        result.dvv.z = S[0][2][2];
    }

    /* duv */
    if ((eval_flags & QAWS_SURFACE_CORE_DUV) && max_u_deriv >= 1 && max_v_deriv >= 1) {
        result.duv.x = S[1][1][0];
        result.duv.y = S[1][1][1];
        result.duv.z = S[1][1][2];
    }

    /* Normal: cross product of du x dv, normalized */
    if (eval_flags & QAWS_SURFACE_CORE_NORMAL) {
        nx = result.du.y * result.dv.z - result.du.z * result.dv.y;
        ny = result.du.z * result.dv.x - result.du.x * result.dv.z;
        nz = result.du.x * result.dv.y - result.du.y * result.dv.x;
        len = QAWS_SQRT(nx * nx + ny * ny + nz * nz);
        if (len > QAWS_LITERAL(1e-12)) {
            result.normal.x = nx / len;
            result.normal.y = ny / len;
            result.normal.z = nz / len;
        } else {
            result.normal.x = QAWS_ZERO;
            result.normal.y = QAWS_ZERO;
            result.normal.z = QAWS_ONE;
        }
    }

    return result;
}

#endif /* QAWS_SURFACE_CORE_H */
