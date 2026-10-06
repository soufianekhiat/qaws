/*
 * Integral functionals of a B-spline surface on the GPU (area, thin plate,
 * Willmore), with the backend-neutral kernels of
 * core/qaws_bspline_surface_functional_core.h.
 *
 * CSNode    one thread per quadrature node (cells x cells cells of 4 x 4
 *           Gauss points): weighted integrand with its first and second
 *           rates along the control point tangents; a sum of node_values
 *           gives the functional and its tangents
 * CSGather  one thread per control point: gradient over the nodes whose
 *           knot window holds it (no atomics, deterministic)
 *
 * Control points are row-major, (i * v_count + j) * 3, i along u.
 */

#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_bspline_surface_functional_core.h"

StructuredBuffer<float> u_knots : register(t0);
StructuredBuffer<float> v_knots : register(t1);
StructuredBuffer<float> control_points : register(t2);   /* 3 per point */
StructuredBuffer<float> cp_tangents : register(t3);      /* 3 per point */

RWStructuredBuffer<float> node_values : register(u0);    /* 3 per node: value, rate, second rate */
RWStructuredBuffer<float> cp_gradient : register(u1);    /* 3 per point */

cbuffer Params : register(b0) {
    uint u_degree;
    uint v_degree;
    uint u_count;
    uint v_count;
    uint cells;
    int functional;          /* QAWS_CORE_SURFACE_AREA, _THIN_PLATE, _WILLMORE */
};

struct window {
    float cp[QAWS_CORE_SURFACE_WINDOW];
    float dot[QAWS_CORE_SURFACE_WINDOW];
    float uk[QAWS_CORE_MAX_POINTS * 2];
    float vk[QAWS_CORE_MAX_POINTS * 2];
    int ufirst, vfirst;
};

window load_window(float u, float v)
{
    window w;
    float all[QAWS_CORE_MAX_POINTS * 2];
    int us, vs, i, j;
    for (i = 0; i < QAWS_CORE_MAX_POINTS * 2; i++)
        all[i] = i < (int)(u_count + u_degree + 1) ? u_knots[i] : 0.0;
    us = qaws_find_span(all, (int)u_degree, (int)u_count, u);
    for (i = 0; i < QAWS_CORE_MAX_POINTS * 2; i++)
        all[i] = i < (int)(v_count + v_degree + 1) ? v_knots[i] : 0.0;
    vs = qaws_find_span(all, (int)v_degree, (int)v_count, v);
    for (i = 0; i < QAWS_CORE_MAX_POINTS * 2; i++) {
        w.uk[i] = i < (int)(2 * (u_degree + 1)) ? u_knots[us - (int)u_degree + i] : 0.0;
        w.vk[i] = i < (int)(2 * (v_degree + 1)) ? v_knots[vs - (int)v_degree + i] : 0.0;
    }
    w.ufirst = us - (int)u_degree;
    w.vfirst = vs - (int)v_degree;
    for (i = 0; i < QAWS_CORE_SURFACE_WINDOW; i++) {
        w.cp[i] = 0.0;
        w.dot[i] = 0.0;
    }
    for (i = 0; i <= (int)u_degree; i++)
        for (j = 0; j <= (int)v_degree; j++) {
            uint g = ((uint)(w.ufirst + i) * v_count + (uint)(w.vfirst + j)) * 3;
            int l = (i * (int)(v_degree + 1) + j) * 3;
            w.cp[l + 0] = control_points[g + 0];
            w.cp[l + 1] = control_points[g + 1];
            w.cp[l + 2] = control_points[g + 2];
            w.dot[l + 0] = cp_tangents[g + 0];
            w.dot[l + 1] = cp_tangents[g + 1];
            w.dot[l + 2] = cp_tangents[g + 2];
        }
    return w;
}

float3 node_of(uint n)
{
    float out_node[3];
    uint cell = n / 16, q = n % 16;
    qaws_surface_cell_node(0.0, 1.0, 0.0, 1.0, (float)(cell / cells), (float)(cell % cells), (float)cells, (int)q, out_node);
    return float3(out_node[0], out_node[1], out_node[2]);
}

[numthreads(64, 1, 1)]
void CSNode(uint3 id : SV_DispatchThreadID)
{
    uint n = id.x;
    float3 node;
    window w;
    qaws_dual1 f;
    if (n >= cells * cells * 16)
        return;
    node = node_of(n);
    w = load_window(node.x, node.y);
    f = qaws_bspline_surface_node(functional, w.cp, w.dot, w.uk, w.vk, (int)u_degree, (int)v_degree, node.x, node.y, node.z);
    node_values[n * 3 + 0] = f.v;
    node_values[n * 3 + 1] = f.t;
    node_values[n * 3 + 2] = f.tt;
}

[numthreads(64, 1, 1)]
void CSGather(uint3 id : SV_DispatchThreadID)
{
    uint e = id.x, n;
    int i, j;
    float3 acc = float3(0.0, 0.0, 0.0);
    if (e >= u_count * v_count)
        return;
    i = (int)(e / v_count);
    j = (int)(e % v_count);
    for (n = 0; n < cells * cells * 16; n++) {
        float3 node = node_of(n);
        window w = load_window(node.x, node.y);
        qaws_vec3 g;
        if (i < w.ufirst || i > w.ufirst + (int)u_degree || j < w.vfirst || j > w.vfirst + (int)v_degree)
            continue;
        g = qaws_bspline_surface_node_adjoint_cp(functional, w.cp, w.uk, w.vk, (int)u_degree, (int)v_degree, node.x, node.y,
            i - w.ufirst, j - w.vfirst, node.z);
        acc += float3(g.x, g.y, g.z);
    }
    cp_gradient[e * 3 + 0] = acc.x;
    cp_gradient[e * 3 + 1] = acc.y;
    cp_gradient[e * 3 + 2] = acc.z;
}
