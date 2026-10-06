/*
 * Inverse-CDF warping of the unit square onto a B-spline surface by its
 * area measure on the GPU, with the backend-neutral kernels of
 * core/qaws_bspline_surface_sampling_core.h.
 *
 * CSCellMass     one thread per u cell: measure of the cell and its rates
 *                along the control point tangents (an exclusive prefix sum
 *                gives prefix_v / prefix_t / prefix_tt)
 * CSSample       one thread per point xi: (u, v), the position and the
 *                first and second tangents along cp_tangents and xi_tangents
 * CSMultipliers  one thread per sample: adjoint multipliers (lambda, mu)
 *                and the xi adjoints
 *
 * The control point gradient is a gather over lines (see
 * qaws_patch_line_adjoint_cp and test 62_diff_surface_sampling_core).
 * Control points are row-major, (i * v_count + j) * 3, i along u.
 */

#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_bspline_surface_sampling_core.h"

StructuredBuffer<float> u_knot_buffer : register(t0);
StructuredBuffer<float> v_knot_buffer : register(t1);
StructuredBuffer<float> control_points : register(t2);   /* 3 per point */
StructuredBuffer<float> cp_tangents : register(t3);      /* 3 per point */
StructuredBuffer<float> prefix : register(t4);           /* 3 per cell + 1: value, rate, second rate */
StructuredBuffer<float> points : register(t5);           /* 4 per sample: xi_u, xi_v, xi_u rate, xi_v rate */
StructuredBuffer<float> sample_adjoints : register(t6);  /* 5 per sample: u_bar, v_bar, p_bar */

RWStructuredBuffer<float> cell_mass : register(u0);      /* 3 per cell */
RWStructuredBuffer<float> samples : register(u1);        /* 14 per sample: u, v, u', v', u'', v'', p', p'' */
RWStructuredBuffer<float> multipliers : register(u2);    /* 4 per sample: lambda, mu, xi_u adjoint, xi_v adjoint */

cbuffer Params : register(b0) {
    uint u_degree;
    uint v_degree;
    uint u_count;
    uint v_count;
    uint cells;
    uint sample_count;
};

struct patch {
    float cp[QAWS_CORE_SURFACE_GRID];
    float dot[QAWS_CORE_SURFACE_GRID];
    float uk[QAWS_CORE_MAX_POINTS * 2];
    float vk[QAWS_CORE_MAX_POINTS * 2];
};

patch load_patch()
{
    patch p;
    int i;
    for (i = 0; i < QAWS_CORE_SURFACE_GRID; i++) {
        p.cp[i] = i < (int)(u_count * v_count * 3) ? control_points[i] : 0.0;
        p.dot[i] = i < (int)(u_count * v_count * 3) ? cp_tangents[i] : 0.0;
    }
    for (i = 0; i < QAWS_CORE_MAX_POINTS * 2; i++) {
        p.uk[i] = i < (int)(u_count + u_degree + 1) ? u_knot_buffer[i] : 0.0;
        p.vk[i] = i < (int)(v_count + v_degree + 1) ? v_knot_buffer[i] : 0.0;
    }
    return p;
}

#define PATCH_ARGS(p) p.cp, p.dot, p.uk, p.vk, (int)u_degree, (int)v_degree, (int)u_count, (int)v_count

[numthreads(32, 1, 1)]
void CSCellMass(uint3 id : SV_DispatchThreadID)
{
    patch p;
    qaws_dual1 m;
    if (id.x >= cells)
        return;
    p = load_patch();
    m = qaws_patch_cell_mass(PATCH_ARGS(p), 1.0, (int)cells, (int)id.x);
    cell_mass[id.x * 3 + 0] = m.v;
    cell_mass[id.x * 3 + 1] = m.t;
    cell_mass[id.x * 3 + 2] = m.tt;
}

[numthreads(32, 1, 1)]
void CSSample(uint3 id : SV_DispatchThreadID)
{
    uint n = id.x;
    patch p;
    float pv[QAWS_CORE_MAX_CELLS + 1];
    float pt[QAWS_CORE_MAX_CELLS + 1];
    float ptt[QAWS_CORE_MAX_CELLS + 1];
    float uv[2];
    float tg[4];
    qaws_vec3 tp[2];
    int c;
    if (n >= sample_count)
        return;
    p = load_patch();
    for (c = 0; c <= QAWS_CORE_MAX_CELLS; c++) {
        pv[c] = c <= (int)cells ? prefix[c * 3 + 0] : 0.0;
        pt[c] = c <= (int)cells ? prefix[c * 3 + 1] : 0.0;
        ptt[c] = c <= (int)cells ? prefix[c * 3 + 2] : 0.0;
    }
    qaws_patch_cdf_solve(PATCH_ARGS(p), (int)cells, pv, points[n * 4 + 0], points[n * 4 + 1], uv);
    qaws_patch_cdf_tangent(PATCH_ARGS(p), (int)cells, pv, pt, ptt, points[n * 4 + 0], points[n * 4 + 1], points[n * 4 + 2],
        points[n * 4 + 3], uv[0], uv[1], tg, tp);
    samples[n * 14 + 0] = uv[0];
    samples[n * 14 + 1] = uv[1];
    for (c = 0; c < 4; c++)
        samples[n * 14 + 2 + c] = tg[c];
    samples[n * 14 + 6] = tp[0].x; samples[n * 14 + 7] = tp[0].y; samples[n * 14 + 8] = tp[0].z;
    samples[n * 14 + 9] = tp[1].x; samples[n * 14 + 10] = tp[1].y; samples[n * 14 + 11] = tp[1].z;
}

[numthreads(32, 1, 1)]
void CSMultipliers(uint3 id : SV_DispatchThreadID)
{
    uint n = id.x;
    patch p;
    float pv[QAWS_CORE_MAX_CELLS + 1];
    float out_m[4];
    int c;
    if (n >= sample_count)
        return;
    p = load_patch();
    for (c = 0; c <= QAWS_CORE_MAX_CELLS; c++)
        pv[c] = c <= (int)cells ? prefix[c * 3 + 0] : 0.0;
    qaws_patch_cdf_multipliers(PATCH_ARGS(p), (int)cells, pv, points[n * 4 + 0], points[n * 4 + 1], samples[n * 14 + 0],
        samples[n * 14 + 1], sample_adjoints[n * 5 + 0], sample_adjoints[n * 5 + 1],
        qaws_v3(sample_adjoints[n * 5 + 2], sample_adjoints[n * 5 + 3], sample_adjoints[n * 5 + 4]), out_m);
    for (c = 0; c < 4; c++)
        multipliers[n * 4 + c] = out_m[c];
}
