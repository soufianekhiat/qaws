#include "../qaws_platform.h"
#include "../qaws_types.h"
#include "../qaws_core_types.h"
#include "../internal/qaws_internal_basis.h"

/* ===================================================================
 * Batch Bezier tensor-product surface evaluation
 *
 * Evaluates a Bezier surface at a grid of (u,v) parameter pairs.
 * Uses the De Casteljau algorithm in two passes: first reduces along
 * u for each row of control points, then along v on the intermediate
 * results.
 *
 * control_points: flat row-major array of 3D control points,
 *                 (u_degree+1) * (v_degree+1) points,
 *                 [x00,y00,z00, x01,y01,z01, ..., xMN,yMN,zMN]
 *                 Row i contains v_degree+1 points.
 * u_degree:       degree in u direction
 * v_degree:       degree in v direction
 * uv_values:      interleaved parameter pairs [u0,v0, u1,v1, ...]
 * pair_count:     number of (u,v) pairs
 * out_positions:  receives pair_count * 3 scalars [x0,y0,z0, ...]
 * =================================================================== */

void qaws_batch_surface_bezier_eval(
    qaws_scalar const* control_points,
    unsigned int u_degree,
    unsigned int v_degree,
    qaws_scalar const* uv_values,
    unsigned int pair_count,
    qaws_scalar* out_positions)
{
    unsigned int u_count = u_degree + 1;
    unsigned int v_count = v_degree + 1;
    unsigned int idx;

    for (idx = 0; idx < pair_count; idx++) {
        qaws_scalar u = uv_values[idx * 2];
        qaws_scalar v = uv_values[idx * 2 + 1];
        qaws_scalar one_u = QAWS_ONE - u;
        qaws_scalar one_v = QAWS_ONE - v;

        /* First pass: for each row i in u, reduce along v using De Casteljau.
         * This yields u_count intermediate 3D points. */
        qaws_scalar tmp_x[QAWS_CORE_MAX_POINTS];
        qaws_scalar tmp_y[QAWS_CORE_MAX_POINTS];
        qaws_scalar tmp_z[QAWS_CORE_MAX_POINTS];
        unsigned int i, j, r;

        for (i = 0; i < u_count; i++) {
            /* Copy row i control points */
            qaws_scalar vx[QAWS_CORE_MAX_POINTS];
            qaws_scalar vy[QAWS_CORE_MAX_POINTS];
            qaws_scalar vz[QAWS_CORE_MAX_POINTS];

            for (j = 0; j < v_count; j++) {
                unsigned int ci = (i * v_count + j) * 3;
                vx[j] = control_points[ci];
                vy[j] = control_points[ci + 1];
                vz[j] = control_points[ci + 2];
            }

            /* De Casteljau along v */
            for (r = 1; r <= v_degree; r++) {
                for (j = 0; j <= v_degree - r; j++) {
                    vx[j] = one_v * vx[j] + v * vx[j + 1];
                    vy[j] = one_v * vy[j] + v * vy[j + 1];
                    vz[j] = one_v * vz[j] + v * vz[j + 1];
                }
            }

            tmp_x[i] = vx[0];
            tmp_y[i] = vy[0];
            tmp_z[i] = vz[0];
        }

        /* Second pass: De Casteljau along u on the intermediate points */
        for (r = 1; r <= u_degree; r++) {
            for (i = 0; i <= u_degree - r; i++) {
                tmp_x[i] = one_u * tmp_x[i] + u * tmp_x[i + 1];
                tmp_y[i] = one_u * tmp_y[i] + u * tmp_y[i + 1];
                tmp_z[i] = one_u * tmp_z[i] + u * tmp_z[i + 1];
            }
        }

        out_positions[idx * 3]     = tmp_x[0];
        out_positions[idx * 3 + 1] = tmp_y[0];
        out_positions[idx * 3 + 2] = tmp_z[0];
    }
}

/* ===================================================================
 * Batch B-spline tensor-product surface evaluation
 *
 * Evaluates a B-spline surface at a set of (u,v) parameter pairs.
 * Uses basis function evaluation in each direction and combines
 * with a tensor-product summation.
 *
 * control_points:  row-major 3D CPs, u_cp_count * v_cp_count points
 * u_knots/v_knots: full knot vectors
 * u/v_knot_count:  number of knots in each direction
 * u/v_cp_count:    number of CPs in each direction
 * u/v_degree:      degree in each direction
 * uv_values:       interleaved [u0,v0, u1,v1, ...]
 * pair_count:      number of (u,v) pairs
 * out_positions:   receives pair_count * 3 scalars
 * =================================================================== */

void qaws_batch_surface_bspline_eval(
    qaws_scalar const* control_points,
    qaws_scalar const* u_knots,
    unsigned int u_knot_count,
    qaws_scalar const* v_knots,
    unsigned int v_knot_count,
    unsigned int u_cp_count,
    unsigned int v_cp_count,
    unsigned int u_degree,
    unsigned int v_degree,
    qaws_scalar const* uv_values,
    unsigned int pair_count,
    qaws_scalar* out_positions)
{
    unsigned int u_order = u_degree + 1;
    unsigned int v_order = v_degree + 1;
    unsigned int idx;

    for (idx = 0; idx < pair_count; idx++) {
        qaws_scalar u = uv_values[idx * 2];
        qaws_scalar v = uv_values[idx * 2 + 1];
        unsigned int u_span, v_span;
        qaws_scalar u_basis[QAWS_CORE_MAX_POINTS];
        qaws_scalar v_basis[QAWS_CORE_MAX_POINTS];
        qaws_scalar px, py, pz;
        unsigned int ui, vi;

        u_span = qaws_internal_find_knot_span(u_knots, u_knot_count,
                                               u_degree, u_cp_count, u);
        v_span = qaws_internal_find_knot_span(v_knots, v_knot_count,
                                               v_degree, v_cp_count, v);

        qaws_internal_bspline_basis(u_basis, u_knots, u_knot_count,
                                     u_degree, u_span, u);
        qaws_internal_bspline_basis(v_basis, v_knots, v_knot_count,
                                     v_degree, v_span, v);

        /* Tensor-product summation:
         * S(u,v) = sum_i sum_j N_i(u) * N_j(v) * P_{i,j} */
        px = QAWS_ZERO;
        py = QAWS_ZERO;
        pz = QAWS_ZERO;

        for (ui = 0; ui < u_order; ui++) {
            unsigned int row = u_span - u_degree + ui;
            qaws_scalar nu = u_basis[ui];

            for (vi = 0; vi < v_order; vi++) {
                unsigned int col = v_span - v_degree + vi;
                unsigned int ci = (row * v_cp_count + col) * 3;
                qaws_scalar nv = v_basis[vi];
                qaws_scalar w = nu * nv;

                px += w * control_points[ci];
                py += w * control_points[ci + 1];
                pz += w * control_points[ci + 2];
            }
        }

        out_positions[idx * 3]     = px;
        out_positions[idx * 3 + 1] = py;
        out_positions[idx * 3 + 2] = pz;
    }
}

/* ===================================================================
 * Batch NURBS tensor-product surface evaluation
 *
 * Same as B-spline surface but applies the rational quotient rule:
 * S(u,v) = sum(N_i*N_j*w_{ij}*P_{ij}) / sum(N_i*N_j*w_{ij})
 *
 * weights: per-control-point weights, row-major, u_cp_count*v_cp_count
 * =================================================================== */

void qaws_batch_surface_nurbs_eval(
    qaws_scalar const* control_points,
    qaws_scalar const* weights,
    qaws_scalar const* u_knots,
    unsigned int u_knot_count,
    qaws_scalar const* v_knots,
    unsigned int v_knot_count,
    unsigned int u_cp_count,
    unsigned int v_cp_count,
    unsigned int u_degree,
    unsigned int v_degree,
    qaws_scalar const* uv_values,
    unsigned int pair_count,
    qaws_scalar* out_positions)
{
    unsigned int u_order = u_degree + 1;
    unsigned int v_order = v_degree + 1;
    unsigned int idx;

    for (idx = 0; idx < pair_count; idx++) {
        qaws_scalar u = uv_values[idx * 2];
        qaws_scalar v = uv_values[idx * 2 + 1];
        unsigned int u_span, v_span;
        qaws_scalar u_basis[QAWS_CORE_MAX_POINTS];
        qaws_scalar v_basis[QAWS_CORE_MAX_POINTS];
        qaws_scalar ax, ay, az, w_sum;
        unsigned int ui, vi;

        u_span = qaws_internal_find_knot_span(u_knots, u_knot_count,
                                               u_degree, u_cp_count, u);
        v_span = qaws_internal_find_knot_span(v_knots, v_knot_count,
                                               v_degree, v_cp_count, v);

        qaws_internal_bspline_basis(u_basis, u_knots, u_knot_count,
                                     u_degree, u_span, u);
        qaws_internal_bspline_basis(v_basis, v_knots, v_knot_count,
                                     v_degree, v_span, v);

        ax = QAWS_ZERO;
        ay = QAWS_ZERO;
        az = QAWS_ZERO;
        w_sum = QAWS_ZERO;

        for (ui = 0; ui < u_order; ui++) {
            unsigned int row = u_span - u_degree + ui;
            qaws_scalar nu = u_basis[ui];

            for (vi = 0; vi < v_order; vi++) {
                unsigned int col = v_span - v_degree + vi;
                unsigned int flat = row * v_cp_count + col;
                unsigned int ci = flat * 3;
                qaws_scalar nv = v_basis[vi];
                qaws_scalar nw = nu * nv * weights[flat];

                w_sum += nw;
                ax += nw * control_points[ci];
                ay += nw * control_points[ci + 1];
                az += nw * control_points[ci + 2];
            }
        }

        {
            qaws_scalar inv_w = QAWS_ONE / w_sum;
            out_positions[idx * 3]     = ax * inv_w;
            out_positions[idx * 3 + 1] = ay * inv_w;
            out_positions[idx * 3 + 2] = az * inv_w;
        }
    }
}
