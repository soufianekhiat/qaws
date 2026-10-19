#include "../../src/qaws/qaws_hlsl.h"

StructuredBuffer<float> input_data : register(t0);
RWStructuredBuffer<float> results : register(u0);

#include "cores_check.h"

[numthreads(64, 1, 1)]
void CSMain(uint3 id : SV_DispatchThreadID)
{
    results[id.x] = qaws_check_cores((float)id.x / 64.0);
}
