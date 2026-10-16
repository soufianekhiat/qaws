/*
 * Tangent-point energies on the GPU (qaws_repulsive.h, R6), with the
 * backend-neutral pair kernels of core/qaws_repulsive_core.h: the direct
 * all-pairs energy, its gradient, and the fractional Sobolev metric's
 * product, for curves (k = 2: edges) or meshes (k = 3: faces). Storage is
 * O(elements); nothing of size N^2.
 *
 * CSElement  one thread per element I, looping over every element J: its
 *            energy share sum_J e(I, J), its gradient block
 *            sum_J d/dI [e(I, J) + e(J, I)] (3 per corner), and with u its
 *            metric row: accB = sum kB (g_I - g_J), acc0 = sum kB0
 *            (u_I - u_J), and its coefficients
 * CSGather   one thread per vertex: the gradient and y = A u from the
 *            elements holding it (a CSR adjacency; no atomics, so the
 *            result does not depend on scheduling)
 * The energy is the sum of the shares (on the host, in double).
 */

#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_repulsive_core.h"

StructuredBuffer<float> positions : register(t0);       /* 3 per vertex */
StructuredBuffer<uint> indices : register(t1);          /* k per element */
StructuredBuffer<uint> adj_offsets : register(t2);      /* vertex_count + 1 */
StructuredBuffer<uint> adj_entries : register(t3);      /* element * k + corner */
StructuredBuffer<float> u_values : register(t4);        /* 1 per vertex */

RWStructuredBuffer<float> energy_shares : register(u0); /* 1 per element */
RWStructuredBuffer<float> grad_blocks : register(u1);   /* 9 per element */
RWStructuredBuffer<float> metric_rows : register(u2);   /* 4 per element: accB, acc0 */
RWStructuredBuffer<float> coefficients : register(u3);  /* 9 per element */
RWStructuredBuffer<float> gradient : register(u4);      /* 3 per vertex */
RWStructuredBuffer<float> product : register(u5);       /* 1 per vertex */

cbuffer Params : register(b0) {
    uint element_count;
    uint vertex_count;
    uint k;            /* 2 curves, 3 meshes */
    uint pad0;
    float alpha;       /* curves: alpha; meshes: p */
    float beta;        /* curves: beta */
    float ex;          /* the metric's exponent */
    float pad1;
};

void load_element(uint e, out float p[9])
{
    uint a, c;
    for (a = 0; a < 9; a++) p[a] = 0;
    for (a = 0; a < k; a++)
        for (c = 0; c < 3; c++) p[3 * a + c] = positions[3 * indices[k * e + a] + c];
}

[numthreads(64, 1, 1)]
void CSElement(uint3 tid : SV_DispatchThreadID)
{
    uint I = tid.x, J, a, c;
    float p[9];
    float q[9];
    float cI[9];
    float cJ[9];
    float blk[9];
    float gI[3];
    float acc[4];
    float share = 0, uI = 0;
    if (I >= element_count)
        return;
    load_element(I, p);
    for (c = 0; c < 9; c++) { blk[c] = 0; cI[c] = 0; }
    for (c = 0; c < 4; c++) acc[c] = 0;
    if (k == 2) {
        float p6[6]; float c6[6];
        for (c = 0; c < 6; c++) p6[c] = p[c];
        qaws_repulsive_curve_coefficients(p6, c6);
        for (c = 0; c < 6; c++) cI[c] = c6[c];
    } else {
        qaws_repulsive_mesh_coefficients(p, cI);
    }
    for (c = 0; c < 3; c++) gI[c] = 0;
    for (a = 0; a < k; a++) {
        float ua = u_values[indices[k * I + a]];
        uI += ua / (float)k;
        for (c = 0; c < 3; c++) gI[c] += ua * cI[3 * a + c];
    }
    for (J = 0; J < element_count; J++) {
        float e1[1];
        float g[9];
        float kw[2];
        float gJ[3];
        float uJ = 0;
        uint touch = 0, b;
        if (J == I)
            continue;
        if (k == 2) {
            for (a = 0; a < 2; a++)
                for (b = 0; b < 2; b++) touch |= indices[2 * I + a] == indices[2 * J + b] ? 1u : 0u;
            if (touch != 0)
                continue;   /* (curves: edges sharing a vertex do not interact) */
        }
        load_element(J, q);
        for (c = 0; c < 9; c++) cJ[c] = 0;
        if (k == 2) {
            float p6[6]; float q6[6]; float g6[6]; float c6[6];
            for (c = 0; c < 6; c++) { p6[c] = p[c]; q6[c] = q[c]; }
            qaws_repulsive_curve_pair(e1, g6, p6, q6, alpha, beta, 0);
            share += e1[0];
            for (c = 0; c < 6; c++) blk[c] += g6[c];
            qaws_repulsive_curve_pair(e1, g6, q6, p6, alpha, beta, 1);
            for (c = 0; c < 6; c++) blk[c] += g6[c];
            qaws_repulsive_curve_metric_pair(kw, p6, q6, ex);
            qaws_repulsive_curve_coefficients(q6, c6);
            for (c = 0; c < 6; c++) cJ[c] = c6[c];
        } else {
            qaws_repulsive_mesh_pair(e1, g, p, q, alpha, 0);
            share += e1[0];
            for (c = 0; c < 9; c++) blk[c] += g[c];
            qaws_repulsive_mesh_pair(e1, g, q, p, alpha, 1);
            for (c = 0; c < 9; c++) blk[c] += g[c];
            qaws_repulsive_mesh_metric_pair(kw, p, q, ex);
            qaws_repulsive_mesh_coefficients(q, cJ);
        }
        for (c = 0; c < 3; c++) gJ[c] = 0;
        for (a = 0; a < k; a++) {
            float ua = u_values[indices[k * J + a]];
            uJ += ua / (float)k;
            for (c = 0; c < 3; c++) gJ[c] += ua * cJ[3 * a + c];
        }
        for (c = 0; c < 3; c++) acc[c] += kw[0] * (gI[c] - gJ[c]);
        acc[3] += kw[1] * (uI - uJ);
    }
    energy_shares[I] = share;
    for (c = 0; c < 9; c++) { grad_blocks[9 * I + c] = blk[c]; coefficients[9 * I + c] = cI[c]; }
    for (c = 0; c < 4; c++) metric_rows[4 * I + c] = acc[c];
}

[numthreads(64, 1, 1)]
void CSGather(uint3 tid : SV_DispatchThreadID)
{
    uint v = tid.x, a, c;
    float g[3];
    float y = 0;
    if (v >= vertex_count)
        return;
    for (c = 0; c < 3; c++) g[c] = 0;
    for (a = adj_offsets[v]; a < adj_offsets[v + 1]; a++) {
        uint el = adj_entries[a] / k, slot = adj_entries[a] - k * (adj_entries[a] / k);
        for (c = 0; c < 3; c++) {
            g[c] += grad_blocks[9 * el + 3 * slot + c];
            y += coefficients[9 * el + 3 * slot + c] * metric_rows[4 * el + c];
        }
        y += metric_rows[4 * el + 3] / (float)k;
    }
    for (c = 0; c < 3; c++) gradient[3 * v + c] = g[c];
    product[v] = y;
}
