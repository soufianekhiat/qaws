/*
 * The C side of examples/eval_record_halide.cpp: the record cores compiled for C (as test 144 runs them against the
 * runtime) over the same buffers, so the Halide pipelines can be compared with them. f64.
 */

#define QAWS_SCALAR_IS_FLOAT 0
#include "../src/qaws/qaws_platform.h"
#include "../src/qaws/qaws_core_types.h"

static double const* er_knots;
static double const* er_points;
static double const* er_weights;
#define QAWS_BUFFER_KNOT(i) er_knots[i]
#define QAWS_BUFFER_POINT(i) er_points[i]
#define QAWS_BUFFER_WEIGHT(i) er_weights[i]
#include "../src/qaws/core/qaws_curve_record_core.h"
#include "../src/qaws/core/qaws_surface_record_core.h"

/* rec: degree, count, dimension, rational, knot_offset, point_offset, weight_offset; out: 9 per t (position, D1, D2) */
void eval_record_ref_curve(double* out, int const* rec, double const* knots, double const* points, double const* weights,
	double const* t, int n)
{
	qaws_curve_record r;
	int i;
	er_knots = knots; er_points = points; er_weights = weights;
	r.degree = rec[0]; r.count = rec[1]; r.dimension = rec[2]; r.rational = rec[3];
	r.knot_offset = rec[4]; r.point_offset = rec[5]; r.weight_offset = rec[6]; r.reserved = 0;
	for (i = 0; i < n; i++)
	{
		qaws_eval_3d e = qaws_record_eval_3d(r, t[i], 7);
		double* o = out + 9 * i;
		o[0] = e.position.x; o[1] = e.position.y; o[2] = e.position.z;
		o[3] = e.d1.x; o[4] = e.d1.y; o[5] = e.d1.z;
		o[6] = e.d2.x; o[7] = e.d2.y; o[8] = e.d2.z;
	}
}

/* rec: u_degree, v_degree, u_count, v_count, rational, u_knot_offset, v_knot_offset, point_offset, weight_offset;
   out: 21 per (u, v) (position, du, dv, duu, duv, dvv, normal) */
void eval_record_ref_surface(double* out, int const* rec, double const* knots, double const* points, double const* weights,
	double const* uv, int n)
{
	qaws_surface_record r;
	int i;
	er_knots = knots; er_points = points; er_weights = weights;
	r.u_degree = rec[0]; r.v_degree = rec[1]; r.u_count = rec[2]; r.v_count = rec[3]; r.rational = rec[4];
	r.u_knot_offset = rec[5]; r.v_knot_offset = rec[6]; r.point_offset = rec[7]; r.weight_offset = rec[8];
	r.reserved0 = r.reserved1 = r.reserved2 = 0;
	for (i = 0; i < n; i++)
	{
		qaws_surface_eval_3d e = qaws_surface_record_eval(r, uv[2 * i], uv[2 * i + 1], 127);
		qaws_vec3 const* v[7] = { &e.position, &e.du, &e.dv, &e.duu, &e.duv, &e.dvv, &e.normal };
		int k;
		for (k = 0; k < 7; k++)
		{
			out[21 * i + 3 * k] = v[k]->x;
			out[21 * i + 3 * k + 1] = v[k]->y;
			out[21 * i + 3 * k + 2] = v[k]->z;
		}
	}
}
