/*
 * Every evaluation core called once, for the shader compile checks
 * (tests/shaders/check_shaders.ps1): cores_check.hlsl and cores_check.comp
 * include this after declaring input_data, a buffer of at least
 * QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3 values. Compiling is what
 * is checked; the values are not meaningful.
 */

#include "../../src/qaws/core/qaws_arc_core.h"
#include "../../src/qaws/core/qaws_bspline_basis_core.h"
#include "../../src/qaws/core/qaws_bspline_eval_core.h"
#include "../../src/qaws/core/qaws_cubic_poly_core.h"
#include "../../src/qaws/core/qaws_decasteljau_core.h"
#include "../../src/qaws/core/qaws_horner_core.h"
#include "../../src/qaws/core/qaws_nurbs_eval_core.h"
#include "../../src/qaws/core/qaws_rational_bezier_core.h"
#include "../../src/qaws/core/qaws_surface_core.h"

#define QAWS_BUFFER_KNOT(i) input_data[i]
#define QAWS_BUFFER_POINT(i) input_data[i]
#define QAWS_BUFFER_WEIGHT(i) input_data[i]
#include "../../src/qaws/core/qaws_curve_record_core.h"
#include "../../src/qaws/core/qaws_surface_record_core.h"

qaws_scalar qaws_check_cores(qaws_scalar t)
{
    qaws_scalar cp2[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar cp3[QAWS_CORE_MAX_POINTS * 3];
    qaws_scalar cp4[QAWS_CORE_MAX_POINTS * 4];
    qaws_scalar kn[QAWS_CORE_MAX_POINTS * 2];
    qaws_scalar w[QAWS_CORE_MAX_POINTS];
    qaws_scalar grid[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3];
    qaws_vec2 a2;
    qaws_vec2 b2;
    qaws_vec3 a3;
    qaws_vec3 b3;
    qaws_eval_2d e2;
    qaws_eval_3d e3;
    qaws_surface_eval_3d es;
    qaws_curve_record rec;
    qaws_surface_record srec;
    qaws_scalar s = QAWS_ZERO;
    int i;
    int span;

    for (i = 0; i < QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3; i++)
        grid[i] = input_data[i];
    for (i = 0; i < QAWS_CORE_MAX_POINTS * 4; i++)
        cp4[i] = grid[i];
    for (i = 0; i < QAWS_CORE_MAX_POINTS * 3; i++)
        cp3[i] = grid[i];
    for (i = 0; i < QAWS_CORE_MAX_POINTS * 2; i++)
    {
        cp2[i] = grid[i];
        kn[i] = QAWS_ITOF(i);
    }
    for (i = 0; i < QAWS_CORE_MAX_POINTS; i++)
        w[i] = QAWS_ONE;
    a2.x = t; a2.y = QAWS_ONE;
    b2.x = QAWS_ONE; b2.y = t;
    a3.x = t; a3.y = QAWS_ONE; a3.z = QAWS_ZERO;
    b3.x = QAWS_ONE; b3.y = t; b3.z = QAWS_ONE;

    e2 = qaws_arc_eval_full_2d(a2, b2, a2, QAWS_ZERO, QAWS_ONE, t, 7);
    s += e2.position.x + e2.d1.y + e2.d2.x;
    e3 = qaws_arc_eval_full_3d(a3, b3, a3, QAWS_ZERO, QAWS_ONE, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;

    span = qaws_find_span(kn, 3, 6, t);
    e2 = qaws_bspline_eval_2d(cp2, kn, 3, span, t, 7);
    s += e2.position.x + e2.d1.y + e2.d2.x;
    e3 = qaws_bspline_eval_3d(cp3, kn, 3, span, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;

    e2 = qaws_cubic_eval_2d(a2, b2, a2, b2, t, 7);
    s += e2.position.x + e2.d1.y + e2.d2.x;
    e3 = qaws_cubic_eval_3d(a3, b3, a3, b3, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;

    e2 = qaws_decasteljau_eval_2d(cp2, 3, t, 7);
    s += e2.position.x + e2.d1.y + e2.d2.x;
    e3 = qaws_decasteljau_eval_3d(cp3, 3, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;

    e2 = qaws_horner_eval_2d(cp2, 3, t, 7);
    s += e2.position.x + e2.d1.y + e2.d2.x;
    e3 = qaws_horner_eval_3d(cp3, 3, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;

    e2 = qaws_nurbs_eval_2d(cp2, kn, w, 3, span, t, 7);
    s += e2.position.x + e2.d1.y + e2.d2.x;
    e3 = qaws_nurbs_eval_3d(cp3, kn, w, 3, span, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;

    e2 = qaws_rational_bezier_eval_2d(cp3, 3, t, 7);
    s += e2.position.x + e2.d1.y + e2.d2.x;
    e3 = qaws_rational_bezier_eval_3d(cp4, 3, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;

    es = qaws_surface_bezier_eval_3d(grid, 3, 3, 4, t, t, 127);
    s += es.position.x + es.du.y + es.dv.z + es.normal.x + es.duu.y + es.duv.z + es.dvv.x;
    es = qaws_surface_bspline_eval_3d(grid, kn, kn, 3, 3, 6, 6, 3, 3, t, t, 127);
    s += es.position.x + es.du.y + es.dv.z + es.normal.x + es.duu.y + es.duv.z + es.dvv.x;
    rec.degree = 3; rec.count = 40; rec.dimension = 3; rec.rational = 1;
    rec.knot_offset = 0; rec.point_offset = 64; rec.weight_offset = 256; rec.reserved = 0;
    e3 = qaws_record_eval_3d(rec, t, 7);
    s += e3.position.x + e3.d1.y + e3.d2.z;
    srec.u_degree = 3; srec.v_degree = 2; srec.u_count = 8; srec.v_count = 6; srec.rational = 1;
    srec.u_knot_offset = 0; srec.v_knot_offset = 16; srec.point_offset = 64; srec.weight_offset = 256;
    srec.reserved0 = 0; srec.reserved1 = 0; srec.reserved2 = 0;
    es = qaws_surface_record_eval(srec, t, t, 127);
    s += es.position.x + es.du.y + es.dv.z + es.normal.x + es.duu.y + es.duv.z + es.dvv.x;
    return s;
}
