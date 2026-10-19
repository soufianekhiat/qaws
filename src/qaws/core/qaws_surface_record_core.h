#ifndef QAWS_SURFACE_RECORD_CORE_H
#define QAWS_SURFACE_RECORD_CORE_H

#include "qaws_curve_record_core.h"
#include "qaws_surface_core.h"

/* ===================================================================
 * B-spline and NURBS surfaces read from buffers (the shader ABI)
 *
 * As qaws_curve_record_core.h for surfaces: a qaws_surface_record per
 * surface (degrees, point counts, offsets into the knot, point and weight
 * buffers), the same accessor macros, a bounded span search in each knot
 * vector and only the (u_degree + 1) x (v_degree + 1) window of points
 * loaded, so any number of points; degrees at most QAWS_CORE_MAX_POINTS - 1.
 * Point (i, j), i < u_count, j < v_count, is point i * v_count + j (3
 * floats), as qaws_surface_create_bspline / _nurbs; the results are theirs
 * (test 144; on a GPU, examples/eval_gpu_parity.c). eval_flags are the
 * QAWS_SURFACE_CORE_* (= QAWS_SURFACE_EVAL_*) values.
 * =================================================================== */

QAWS_STRUCT(qaws_surface_record) {
    int u_degree;
    int v_degree;
    int u_count;
    int v_count;
    int rational;       /* non-zero: weights at weight_offset, one per point */
    int u_knot_offset;  /* u_count + u_degree + 1 knots */
    int v_knot_offset;  /* v_count + v_degree + 1 knots */
    int point_offset;   /* u_count * v_count * 3 floats */
    int weight_offset;
    int reserved0;
    int reserved1;
    int reserved2;
} QAWS_STRUCT_END(qaws_surface_record);

QAWS_INLINE qaws_surface_eval_3d qaws_surface_record_eval(qaws_surface_record r, qaws_scalar u, qaws_scalar v, int eval_flags)
{
    qaws_scalar lku[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar lkv[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar nu[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar nv[(QAWS_CORE_MAX_DERIV + 1) * QAWS_CORE_MAX_POINTS];
    qaws_scalar A[27];
    qaws_scalar W[9];
    qaws_scalar S[27];
    qaws_scalar inv, len;
    int ud = r.u_degree;
    int vd = r.v_degree;
    int us, vs, i, j, ku, kv, c, k, mu = 0, mv = 0, bu, bv;
    qaws_vec3 n;
    qaws_surface_eval_3d result;

    u = qaws_clamp(u, QAWS_BUFFER_KNOT(r.u_knot_offset + ud), QAWS_BUFFER_KNOT(r.u_knot_offset + r.u_count));
    v = qaws_clamp(v, QAWS_BUFFER_KNOT(r.v_knot_offset + vd), QAWS_BUFFER_KNOT(r.v_knot_offset + r.v_count));
    us = qaws_buffer_find_span(r.u_knot_offset, ud, r.u_count, u);
    vs = qaws_buffer_find_span(r.v_knot_offset, vd, r.v_count, v);
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++) {
        lku[j] = j < 2 * (ud + 1) ? QAWS_BUFFER_KNOT(r.u_knot_offset + us - ud + j) : QAWS_ZERO;
        lkv[j] = j < 2 * (vd + 1) ? QAWS_BUFFER_KNOT(r.v_knot_offset + vs - vd + j) : QAWS_ZERO;
    }
    if ((eval_flags & (QAWS_SURFACE_CORE_DU | QAWS_SURFACE_CORE_DUV | QAWS_SURFACE_CORE_NORMAL)) != 0) mu = 1;
    if ((eval_flags & QAWS_SURFACE_CORE_DUU) != 0) mu = 2;
    if ((eval_flags & (QAWS_SURFACE_CORE_DV | QAWS_SURFACE_CORE_DUV | QAWS_SURFACE_CORE_NORMAL)) != 0) mv = 1;
    if ((eval_flags & QAWS_SURFACE_CORE_DVV) != 0) mv = 2;
    bu = mu < ud ? mu : ud;
    bv = mv < vd ? mv : vd;
    qaws_bspline_basis_derivs(lku, ud, ud, u, bu, nu);
    qaws_bspline_basis_derivs(lkv, vd, vd, v, bv, nv);

    for (k = 0; k < 27; k++) A[k] = QAWS_ZERO;
    for (k = 0; k < 9; k++) W[k] = QAWS_ZERO;
    for (i = 0; i <= ud; i++) {
        for (j = 0; j <= vd; j++) {
            int idx = (us - ud + i) * r.v_count + (vs - vd + j);
            qaws_scalar x[3];
            qaws_scalar w = QAWS_ONE;
            for (c = 0; c < 3; c++)
                x[c] = QAWS_BUFFER_POINT(r.point_offset + idx * 3 + c);
            if (r.rational != 0)
                w = QAWS_BUFFER_WEIGHT(r.weight_offset + idx);
            for (ku = 0; ku <= bu; ku++) {
                for (kv = 0; kv <= bv; kv++) {
                    qaws_scalar bw = nu[ku * (ud + 1) + i] * nv[kv * (vd + 1) + j] * w;
                    W[ku * 3 + kv] += bw;
                    for (c = 0; c < 3; c++)
                        A[(ku * 3 + kv) * 3 + c] += bw * x[c];
                }
            }
        }
    }
    for (k = 0; k < 27; k++) S[k] = A[k];
    if (r.rational != 0) {
        /* the quotient rule: S = A / W and its partials to second order */
        inv = QAWS_ONE / W[0];
        for (c = 0; c < 3; c++) {
            qaws_scalar s00 = A[c] * inv;
            qaws_scalar s10 = (A[9 + c] - s00 * W[3]) * inv;
            qaws_scalar s01 = (A[3 + c] - s00 * W[1]) * inv;
            S[c] = s00;
            S[9 + c] = s10;
            S[3 + c] = s01;
            S[18 + c] = (A[18 + c] - QAWS_LITERAL(2.0) * s10 * W[3] - s00 * W[6]) * inv;
            S[6 + c] = (A[6 + c] - QAWS_LITERAL(2.0) * s01 * W[1] - s00 * W[2]) * inv;
            S[12 + c] = (A[12 + c] - s10 * W[1] - s01 * W[3] - s00 * W[4]) * inv;
        }
    }
    /* (k, l) at (k * 3 + l) * 3 */
    result.position.x = S[0];  result.position.y = S[1];  result.position.z = S[2];
    result.du.x = S[9];  result.du.y = S[10]; result.du.z = S[11];
    result.dv.x = S[3];  result.dv.y = S[4];  result.dv.z = S[5];
    result.duu.x = S[18]; result.duu.y = S[19]; result.duu.z = S[20];
    result.duv.x = S[12]; result.duv.y = S[13]; result.duv.z = S[14];
    result.dvv.x = S[6]; result.dvv.y = S[7]; result.dvv.z = S[8];
    n.x = result.du.y * result.dv.z - result.du.z * result.dv.y;
    n.y = result.du.z * result.dv.x - result.du.x * result.dv.z;
    n.z = result.du.x * result.dv.y - result.du.y * result.dv.x;
    len = QAWS_SQRT(n.x * n.x + n.y * n.y + n.z * n.z);
    if (len > QAWS_LITERAL(1e-12)) {
        result.normal.x = n.x / len; result.normal.y = n.y / len; result.normal.z = n.z / len;
    } else {
        result.normal.x = QAWS_ZERO; result.normal.y = QAWS_ZERO; result.normal.z = QAWS_ONE;
    }
    return result;
}

#endif /* QAWS_SURFACE_RECORD_CORE_H */
