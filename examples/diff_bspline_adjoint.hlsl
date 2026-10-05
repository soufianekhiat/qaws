/*
 * B-spline derivatives on the GPU with the backend-neutral kernels of
 * core/qaws_bspline_diff_core.h.
 *
 * CSTangent        one thread per sample: tangent jet along control point
 *                  tangents and a parameter tangent
 * CSGather         one thread per control point: sums the jet adjoints of
 *                  the samples it supports (CSR lists from
 *                  qaws_curve_build_support_index). No atomics, the result
 *                  is deterministic.
 * CSParamAdjoint   one thread per sample: adjoint of its parameter
 *
 * Jets and jet adjoints are 9 floats per sample: position, d1, d2.
 */

#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_bspline_basis_core.h"
#include "../src/qaws/core/qaws_bspline_diff_core.h"

StructuredBuffer<float> knots : register(t0);
StructuredBuffer<float> control_points : register(t1);    /* 3 per point */
StructuredBuffer<float> params : register(t2);            /* t per sample */
StructuredBuffer<float> jet_adjoints : register(t3);      /* 9 per sample */
StructuredBuffer<uint> support_offsets : register(t4);    /* control_point_count + 1 */
StructuredBuffer<uint> support_samples : register(t5);
StructuredBuffer<float> cp_tangents : register(t6);       /* 3 per point */
StructuredBuffer<float> param_tangents : register(t7);    /* per sample */

RWStructuredBuffer<float> cp_adjoints : register(u0);     /* 3 per point */
RWStructuredBuffer<float> param_adjoints : register(u1);  /* per sample */
RWStructuredBuffer<float> tangents : register(u2);        /* 9 per sample */

cbuffer Params : register(b0) {
    uint degree;
    uint num_cp;
    uint sample_count;
};

/* Span of t and the local window around it. */
int load_window(float t, out float local_knots[QAWS_CORE_MAX_POINTS * 2], out float local_cp[QAWS_CORE_MAX_POINTS * 3])
{
    float all_knots[QAWS_CORE_MAX_POINTS * 2];
    int span, j;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
        all_knots[j] = j < (int)(num_cp + degree + 1) ? knots[j] : 0.0;
    span = qaws_find_span(all_knots, (int)degree, (int)num_cp, t);
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
        local_knots[j] = j < (int)(2 * (degree + 1)) ? knots[span - (int)degree + j] : 0.0;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 3; j++)
        local_cp[j] = 0.0;
    for (j = 0; j <= (int)degree; j++) {
        uint c = (uint)(span - (int)degree + j) * 3;
        local_cp[j * 3 + 0] = control_points[c + 0];
        local_cp[j * 3 + 1] = control_points[c + 1];
        local_cp[j * 3 + 2] = control_points[c + 2];
    }
    return span;
}

qaws_eval_3d load_jet_adjoint(uint s)
{
    qaws_eval_3d y;
    y.position.x = jet_adjoints[s * 9 + 0]; y.position.y = jet_adjoints[s * 9 + 1]; y.position.z = jet_adjoints[s * 9 + 2];
    y.d1.x = jet_adjoints[s * 9 + 3]; y.d1.y = jet_adjoints[s * 9 + 4]; y.d1.z = jet_adjoints[s * 9 + 5];
    y.d2.x = jet_adjoints[s * 9 + 6]; y.d2.y = jet_adjoints[s * 9 + 7]; y.d2.z = jet_adjoints[s * 9 + 8];
    return y;
}

[numthreads(64, 1, 1)]
void CSTangent(uint3 id : SV_DispatchThreadID)
{
    float local_knots[QAWS_CORE_MAX_POINTS * 2];
    float local_cp[QAWS_CORE_MAX_POINTS * 3];
    float local_dot[QAWS_CORE_MAX_POINTS * 3];
    qaws_eval_3d r;
    int span, j;
    uint s = id.x;
    if (s >= sample_count)
        return;
    span = load_window(params[s], local_knots, local_cp);
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 3; j++)
        local_dot[j] = 0.0;
    for (j = 0; j <= (int)degree; j++) {
        uint c = (uint)(span - (int)degree + j) * 3;
        local_dot[j * 3 + 0] = cp_tangents[c + 0];
        local_dot[j * 3 + 1] = cp_tangents[c + 1];
        local_dot[j * 3 + 2] = cp_tangents[c + 2];
    }
    r = qaws_bspline_tangent_3d(local_cp, local_dot, local_knots, (int)degree, (int)degree, params[s], param_tangents[s]);
    tangents[s * 9 + 0] = r.position.x; tangents[s * 9 + 1] = r.position.y; tangents[s * 9 + 2] = r.position.z;
    tangents[s * 9 + 3] = r.d1.x; tangents[s * 9 + 4] = r.d1.y; tangents[s * 9 + 5] = r.d1.z;
    tangents[s * 9 + 6] = r.d2.x; tangents[s * 9 + 7] = r.d2.y; tangents[s * 9 + 8] = r.d2.z;
}

[numthreads(64, 1, 1)]
void CSGather(uint3 id : SV_DispatchThreadID)
{
    float local_knots[QAWS_CORE_MAX_POINTS * 2];
    float local_cp[QAWS_CORE_MAX_POINTS * 3];
    float3 acc = float3(0.0, 0.0, 0.0);
    uint e = id.x, k;
    if (e >= num_cp)
        return;
    for (k = support_offsets[e]; k < support_offsets[e + 1]; k++) {
        uint s = support_samples[k];
        int span = load_window(params[s], local_knots, local_cp);
        qaws_vec3 g = qaws_bspline_adjoint_cp_3d(local_knots, (int)degree, (int)degree, params[s],
            (int)e - (span - (int)degree), load_jet_adjoint(s));
        acc += float3(g.x, g.y, g.z);
    }
    cp_adjoints[e * 3 + 0] += acc.x;
    cp_adjoints[e * 3 + 1] += acc.y;
    cp_adjoints[e * 3 + 2] += acc.z;
}

[numthreads(64, 1, 1)]
void CSParamAdjoint(uint3 id : SV_DispatchThreadID)
{
    float local_knots[QAWS_CORE_MAX_POINTS * 2];
    float local_cp[QAWS_CORE_MAX_POINTS * 3];
    uint s = id.x;
    if (s >= sample_count)
        return;
    load_window(params[s], local_knots, local_cp);
    param_adjoints[s] += qaws_bspline_adjoint_t_3d(local_cp, local_knots, (int)degree, (int)degree, params[s],
        load_jet_adjoint(s));
}
