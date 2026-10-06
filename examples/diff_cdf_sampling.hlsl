/*
 * Inverse-CDF sampling of a B-spline curve on the GPU, first and second
 * order forward and the adjoint, with the backend-neutral kernels of
 * core/qaws_bspline_sampling_core.h.
 *
 * CSSpanMeasure  one thread per span: measure of the span and its first and
 *                second rates along the control point tangents (an
 *                exclusive prefix sum of these, on the host or by a scan
 *                pass, gives `prefix`)
 * CSSample       one thread per sample: target sigma = distance + fraction
 *                total, its span, the solved parameter, the position and the
 *                first and second tangents
 * CSLambda       one thread per sample: adjoint multiplier
 *                lambda = (t_bar + p_bar . C') / m (also the distance adjoint)
 * CSGather       one thread per control point: the gradient of
 *                sum(t_bar t + p_bar . p) over the quadrature nodes and
 *                samples of the spans it supports, with the span weights
 *                F - sum of lambda beyond the span (F = sum lambda fraction).
 *                No atomics, deterministic.
 *
 * Spans are the non-empty knot intervals; span_knot[k] is the knot index
 * of the start of span k. Samples are listed per span in CSR form
 * (span_offsets / span_samples) for the gather.
 */

#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_bspline_sampling_core.h"

StructuredBuffer<float> knots : register(t0);
StructuredBuffer<float> control_points : register(t1);   /* 3 per point */
StructuredBuffer<float> cp_tangents : register(t2);      /* 3 per point */
StructuredBuffer<uint> span_knot : register(t3);         /* span_count */
StructuredBuffer<float> prefix : register(t4);           /* 3 per span + 1: measure, rate, second rate */
StructuredBuffer<float> targets : register(t5);          /* 4 per sample: distance, fraction, distance rate, second rate */
StructuredBuffer<float> sample_adjoints : register(t6);  /* 4 per sample: t_bar, p_bar */
StructuredBuffer<uint> span_offsets : register(t7);      /* span_count + 1 */
StructuredBuffer<uint> span_samples : register(t8);
StructuredBuffer<float> span_weights : register(t9);     /* per span: F - sum of lambda beyond it */
StructuredBuffer<float> lambdas : register(t10);         /* per sample */

RWStructuredBuffer<float> span_measures : register(u0);  /* 3 per span */
RWStructuredBuffer<float> samples : register(u1);        /* 12 per sample: t, p, t', p', t'', p'' */
RWStructuredBuffer<uint> sample_spans : register(u2);    /* per sample */
RWStructuredBuffer<float> lambda_out : register(u3);     /* per sample */
RWStructuredBuffer<float> cp_adjoints : register(u4);    /* 3 per point */

cbuffer Params : register(b0) {
    uint degree;
    uint num_cp;
    uint span_count;
    uint sample_count;
    int integrand;          /* QAWS_CORE_INTEGRAND_SPEED or _CURVATURE */
    float curvature_floor;
};

struct window {
    float lk[QAWS_CORE_MAX_POINTS * 2];
    float lcp[QAWS_CORE_MAX_POINTS * 3];
    float ldot[QAWS_CORE_MAX_POINTS * 3];
    float a, b;
    int first;
};

window load_span(uint k)
{
    window w;
    int s = (int)span_knot[k], j;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
        w.lk[j] = j < (int)(2 * (degree + 1)) ? knots[s - (int)degree + j] : 0.0;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 3; j++) {
        w.lcp[j] = 0.0;
        w.ldot[j] = 0.0;
    }
    for (j = 0; j <= (int)degree; j++) {
        uint c = (uint)(s - (int)degree + j) * 3;
        w.lcp[j * 3 + 0] = control_points[c + 0];
        w.lcp[j * 3 + 1] = control_points[c + 1];
        w.lcp[j * 3 + 2] = control_points[c + 2];
        w.ldot[j * 3 + 0] = cp_tangents[c + 0];
        w.ldot[j * 3 + 1] = cp_tangents[c + 1];
        w.ldot[j * 3 + 2] = cp_tangents[c + 2];
    }
    w.a = knots[s];
    w.b = knots[s + 1];
    w.first = s - (int)degree;
    return w;
}

[numthreads(64, 1, 1)]
void CSSpanMeasure(uint3 id : SV_DispatchThreadID)
{
    uint k = id.x;
    window w;
    qaws_dual1 m;
    if (k >= span_count)
        return;
    w = load_span(k);
    m = qaws_bspline_integrate(integrand, curvature_floor, w.lcp, w.ldot, w.lk, (int)degree, (int)degree, w.a, w.b);
    span_measures[k * 3 + 0] = m.v;
    span_measures[k * 3 + 1] = m.t;
    span_measures[k * 3 + 2] = m.tt;
}

[numthreads(64, 1, 1)]
void CSSample(uint3 id : SV_DispatchThreadID)
{
    uint i = id.x, k = 0;
    float total, sigma, t;
    float ot[2];
    qaws_vec3 op[2];
    qaws_vec3 c[4];
    qaws_vec3 cd[3];
    qaws_dual1 part;
    window w;
    if (i >= sample_count)
        return;
    total = prefix[span_count * 3];
    sigma = clamp(targets[i * 4 + 0] + targets[i * 4 + 1] * total, 0.0, total);
    while (k + 1 < span_count && prefix[(k + 1) * 3] <= sigma)
        k++;
    w = load_span(k);
    t = qaws_bspline_cdf_solve(integrand, curvature_floor, w.lcp, w.lk, (int)degree, (int)degree, w.a, w.b, sigma - prefix[k * 3]);
    part = qaws_bspline_integrate(integrand, curvature_floor, w.lcp, w.ldot, w.lk, (int)degree, (int)degree, w.a, t);
    qaws_bspline_cdf_tangent(integrand, curvature_floor, w.lcp, w.ldot, w.lk, (int)degree, (int)degree, t,
        targets[i * 4 + 2] + targets[i * 4 + 1] * prefix[span_count * 3 + 1],
        targets[i * 4 + 3] + targets[i * 4 + 1] * prefix[span_count * 3 + 2],
        prefix[k * 3 + 1] + part.t, prefix[k * 3 + 2] + part.tt, ot, op);
    qaws_bspline_jet_rows(w.lcp, w.ldot, w.lk, (int)degree, (int)degree, t, c, cd);
    samples[i * 12 + 0] = t;
    samples[i * 12 + 1] = c[0].x; samples[i * 12 + 2] = c[0].y; samples[i * 12 + 3] = c[0].z;
    samples[i * 12 + 4] = ot[0];
    samples[i * 12 + 5] = op[0].x; samples[i * 12 + 6] = op[0].y; samples[i * 12 + 7] = op[0].z;
    samples[i * 12 + 8] = ot[1];
    samples[i * 12 + 9] = op[1].x; samples[i * 12 + 10] = op[1].y; samples[i * 12 + 11] = op[1].z;
    sample_spans[i] = k;
}

[numthreads(64, 1, 1)]
void CSLambda(uint3 id : SV_DispatchThreadID)
{
    uint i = id.x;
    window w;
    qaws_vec3 pb;
    if (i >= sample_count)
        return;
    w = load_span(sample_spans[i]);
    pb = qaws_v3(sample_adjoints[i * 4 + 1], sample_adjoints[i * 4 + 2], sample_adjoints[i * 4 + 3]);
    lambda_out[i] = qaws_bspline_cdf_lambda(integrand, curvature_floor, w.lcp, w.lk, (int)degree, (int)degree,
        samples[i * 12], sample_adjoints[i * 4], pb);
}

[numthreads(64, 1, 1)]
void CSGather(uint3 id : SV_DispatchThreadID)
{
    uint e = id.x, k, n, q;
    float3 acc = float3(0.0, 0.0, 0.0);
    float xw[2];
    if (e >= num_cp)
        return;
    for (k = 0; k < span_count; k++) {
        window w = load_span(k);
        int j = (int)e - w.first;
        float h = (w.b - w.a) / QAWS_SAMPLING_PIECES;
        qaws_vec3 g;
        if (j < 0 || j > (int)degree)
            continue;
        /* full span nodes */
        for (q = 0; q < 8 * QAWS_SAMPLING_PIECES; q++) {
            qaws_gauss8((int)(q % 8), xw);
            g = qaws_bspline_node_adjoint_cp(integrand, curvature_floor, w.lcp, w.lk, (int)degree, (int)degree,
                w.a + h * (float)(q / 8) + 0.5 * h * (1.0 + xw[0]), j, span_weights[k] * 0.5 * h * xw[1]);
            acc += float3(g.x, g.y, g.z);
        }
        /* samples in the span: partial integrals and positions */
        for (n = span_offsets[k]; n < span_offsets[k + 1]; n++) {
            uint i = span_samples[n];
            float t = samples[i * 12];
            float hp = (t - w.a) / QAWS_SAMPLING_PIECES;
            qaws_eval_3d yb;
            for (q = 0; q < 8 * QAWS_SAMPLING_PIECES; q++) {
                qaws_gauss8((int)(q % 8), xw);
                g = qaws_bspline_node_adjoint_cp(integrand, curvature_floor, w.lcp, w.lk, (int)degree, (int)degree,
                    w.a + hp * (float)(q / 8) + 0.5 * hp * (1.0 + xw[0]), j, -lambdas[i] * 0.5 * hp * xw[1]);
                acc += float3(g.x, g.y, g.z);
            }
            yb.position = qaws_v3(sample_adjoints[i * 4 + 1], sample_adjoints[i * 4 + 2], sample_adjoints[i * 4 + 3]);
            yb.d1 = qaws_v3_zero();
            yb.d2 = qaws_v3_zero();
            g = qaws_bspline_adjoint_cp_3d(w.lk, (int)degree, (int)degree, t, j, yb);
            acc += float3(g.x, g.y, g.z);
        }
    }
    cp_adjoints[e * 3 + 0] += acc.x;
    cp_adjoints[e * 3 + 1] += acc.y;
    cp_adjoints[e * 3 + 2] += acc.z;
}
