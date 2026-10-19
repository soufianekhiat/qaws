/*
 * The evaluation cores on the GPU, one family per dispatch, for
 * examples/eval_gpu_parity.c: thread i evaluates parameter i (curves: t;
 * surfaces: u, v) and writes 9 values (position, d1, d2) or 21 (position,
 * du, dv, normal, duu, duv, dvv). The control data is one float buffer read
 * at the offsets given.
 */

#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_arc_core.h"
#include "../src/qaws/core/qaws_bspline_basis_core.h"
#include "../src/qaws/core/qaws_bspline_eval_core.h"
#include "../src/qaws/core/qaws_cubic_poly_core.h"
#include "../src/qaws/core/qaws_decasteljau_core.h"
#include "../src/qaws/core/qaws_horner_core.h"
#include "../src/qaws/core/qaws_nurbs_eval_core.h"
#include "../src/qaws/core/qaws_rational_bezier_core.h"
#include "../src/qaws/core/qaws_surface_core.h"

StructuredBuffer<float> data : register(t0);
StructuredBuffer<float> params : register(t1);
RWStructuredBuffer<float> results : register(u0);

cbuffer Params : register(b0)
{
    uint family;      /* 0 Bezier, 1 rational Bezier, 2 B-spline, 3 NURBS, 4 cubic, 5 polynomial, 6 arc, 7 Bezier surface, 8 B-spline surface */
    uint count;
    uint degree;      /* u degree for surfaces */
    uint n;           /* points (u count for surfaces) */
    uint v_degree;
    uint v_count;
    uint off_knots;   /* (u) knots */
    uint off_extra;   /* weights, or the v knots */
};

void put3(uint base, qaws_vec3 v)
{
    results[base] = v.x;
    results[base + 1] = v.y;
    results[base + 2] = v.z;
}

/* the span of t in knots[off .. off + cnt + deg] and its window */
int window(uint off, int deg, int cnt, float t, out float lk[QAWS_CORE_MAX_POINTS * 2])
{
    float all[QAWS_CORE_MAX_POINTS * 2];
    int j;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
        all[j] = j < cnt + deg + 1 ? data[off + j] : 0.0;
    int span = qaws_find_span(all, deg, cnt, t);
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
        lk[j] = j < 2 * (deg + 1) ? all[span - deg + j] : 0.0;
    return span;
}

/* (one entry per surface family and one for the curves: shader model 5 caps the registers of one entry at 4096) */
void put_surface(uint i, qaws_surface_eval_3d s)
{
    put3(21 * i, s.position); put3(21 * i + 3, s.du); put3(21 * i + 6, s.dv); put3(21 * i + 9, s.normal);
    put3(21 * i + 12, s.duu); put3(21 * i + 15, s.duv); put3(21 * i + 18, s.dvv);
}

[numthreads(64, 1, 1)]
void CSBezierSurface(uint3 id : SV_DispatchThreadID)
{
    uint i = id.x;
    int j;
    float grid[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3];
    if (i >= count)
        return;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3; j++)
        grid[j] = j < (int)(n * v_count * 3) ? data[j] : 0.0;
    put_surface(i, qaws_surface_bezier_eval_3d(grid, (int)degree, (int)v_degree, (int)v_count, params[2 * i], params[2 * i + 1], 127));
}

[numthreads(64, 1, 1)]
void CSBsplineSurface(uint3 id : SV_DispatchThreadID)
{
    uint i = id.x;
    int j;
    float grid[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3];
    float uk[QAWS_CORE_MAX_POINTS * 2], vk[QAWS_CORE_MAX_POINTS * 2];
    float u, v;
    if (i >= count)
        return;
    u = params[2 * i];
    v = params[2 * i + 1];
    for (j = 0; j < QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3; j++)
        grid[j] = j < (int)(n * v_count * 3) ? data[j] : 0.0;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
    {
        uk[j] = j < (int)(n + degree + 1) ? data[off_knots + j] : 0.0;
        vk[j] = j < (int)(v_count + v_degree + 1) ? data[off_extra + j] : 0.0;
    }
    put_surface(i, qaws_surface_bspline_eval_3d(grid, uk, vk, (int)degree, (int)v_degree, (int)n, (int)v_count, qaws_find_span(uk, (int)degree, (int)n, u),
        qaws_find_span(vk, (int)v_degree, (int)v_count, v), u, v, 127));
}

[numthreads(64, 1, 1)]
void CSCurve(uint3 id : SV_DispatchThreadID)
{
    uint i = id.x;
    if (i >= count)
        return;
    int d = (int)degree, j;
    float t = params[i];
    float cp[QAWS_CORE_MAX_POINTS * 4];
    float cp3[QAWS_CORE_MAX_POINTS * 3];
    float lk[QAWS_CORE_MAX_POINTS * 2];
    float lw[QAWS_CORE_MAX_POINTS];
    qaws_eval_3d e;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 4; j++)
        cp[j] = j < (int)(n * 4) ? data[j] : 0.0;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 3; j++)
        cp3[j] = cp[j];
    if (family == 0)
        e = qaws_decasteljau_eval_3d(cp3, d, t, 7);
    else if (family == 1)
        e = qaws_rational_bezier_eval_3d(cp, d, t, 7);
    else if (family == 2 || family == 3)
    {
        int span = window(off_knots, d, (int)n, t, lk);
        float wcp[QAWS_CORE_MAX_POINTS * 3];
        for (j = 0; j < QAWS_CORE_MAX_POINTS * 3; j++)
            wcp[j] = j < 3 * (d + 1) ? data[(span - d) * 3 + j] : 0.0;
        for (j = 0; j < QAWS_CORE_MAX_POINTS; j++)
            lw[j] = family == 3 && j <= d ? data[off_extra + span - d + j] : 1.0;
        if (family == 2)
            e = qaws_bspline_eval_3d(wcp, lk, d, d, t, 7);
        else
            e = qaws_nurbs_eval_3d(wcp, lk, lw, d, d, t, 7);
    }
    else if (family == 4)
    {
        qaws_vec3 a, b, c, dd;
        a.x = cp[0]; a.y = cp[1]; a.z = cp[2];
        b.x = cp[3]; b.y = cp[4]; b.z = cp[5];
        c.x = cp[6]; c.y = cp[7]; c.z = cp[8];
        dd.x = cp[9]; dd.y = cp[10]; dd.z = cp[11];
        e = qaws_cubic_eval_3d(a, b, c, dd, t, 7);
    }
    else if (family == 5)
        e = qaws_horner_eval_3d(cp3, d, t, 7);
    else
    {
        qaws_vec3 center, au, av;
        center.x = cp[0]; center.y = cp[1]; center.z = cp[2];
        au.x = cp[3]; au.y = cp[4]; au.z = cp[5];
        av.x = cp[6]; av.y = cp[7]; av.z = cp[8];
        e = qaws_arc_eval_full_3d(center, au, av, cp[9], cp[10], t, 7);
    }
    put3(9 * i, e.position); put3(9 * i + 3, e.d1); put3(9 * i + 6, e.d2);
}
