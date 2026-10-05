#include "../src/qaws/qaws_hlsl.h"
#include "../src/qaws/core/qaws_bspline_basis_core.h"
#include "../src/qaws/core/qaws_bspline_eval_core.h"

StructuredBuffer<float> knots : register(t0);
StructuredBuffer<float> control_points : register(t1);
RWStructuredBuffer<float> results : register(u0);

cbuffer Params : register(b0) {
    uint degree;
    uint num_cp;
    uint sample_count;
};

[numthreads(256, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID)
{
    if (id.x >= sample_count)
        return;

    /* Map thread to parameter value */
    float t_min = knots[degree];
    float t_max = knots[num_cp];
    float t = t_min + (float)id.x / (float)(sample_count - 1) * (t_max - t_min);

    /* Find knot span on the full knot vector */
    float all_knots[QAWS_CORE_MAX_POINTS * 2];
    float local_knots[QAWS_CORE_MAX_POINTS * 2];
    float local_cp[QAWS_CORE_MAX_POINTS * 2];
    int j;
    for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
        all_knots[j] = j < (int)(num_cp + degree + 1) ? knots[j] : 0.0;
    int span = qaws_find_span(all_knots, (int)degree, (int)num_cp, t);

    /* Extract local knot window and control points */
    for (j = 0; j < (int)(2 * (degree + 1)); j++)
        local_knots[j] = knots[span - (int)degree + j];

    for (j = 0; j <= (int)degree; j++) {
        uint cp_idx = (span - degree + (uint)j) * 2;
        local_cp[j * 2 + 0] = control_points[cp_idx + 0];
        local_cp[j * 2 + 1] = control_points[cp_idx + 1];
    }

    qaws_eval_2d result = qaws_bspline_eval_2d(
        local_cp, local_knots, (int)degree, (int)degree, t, 1);

    results[id.x * 2]     = result.position.x;
    results[id.x * 2 + 1] = result.position.y;
}
