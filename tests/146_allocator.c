/*
 * Test 146: the allocator through every curve and surface constructor
 *
 * Each constructor's _ex variant runs with an allocator failing at its k-th
 * allocation, for k = 0, 1, ... until it succeeds. Every failure must give
 * QAWS_STATUS_ALLOCATION_FAILURE, no curve, and nothing left allocated;
 * the success, once destroyed, must leave nothing either. The allocator
 * counts what it hands out and gets back, so memory taken from malloc but
 * released through it (or the reverse) shows as a count that does not
 * return to zero.
 */

#include "test_common.h"
#include "qaws_alloc.h"
#include "qaws_view.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef struct al_state
{
	int calls, fail_at, live;
} al_state;

static void* al_alloc(unsigned long size, void* u)
{
	al_state* s = (al_state*)u;
	if (s->calls++ == s->fail_at)
		return NULL;
	s->live++;
	return malloc(size ? size : 1);
}

static void al_dealloc(void* ptr, void* u)
{
	if (ptr)
		((al_state*)u)->live--;
	free(ptr);
}

static qaws_vec3 const al_pts[6] = { { 0, 0, 0 }, { 1, 2, 0.5f }, { 2.5f, 1.5f, -0.5f }, { 3, -1, 1 }, { 4.5f, 0, 0.2f }, { 5, 1, 0 } };
static qaws_vec3 const al_tan[6] = { { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 0, -1, 1 }, { 1, 0, -1 }, { 1, 1, 1 } };

/* curve constructors: case c with the allocator into *out */
static qaws_status al_curve(unsigned int c, qaws_allocator const* al, qaws_curve** out)
{
	*out = NULL;
	switch (c)
	{
	case 0:
	{
		qaws_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = al_pts; d.control_point_count = 4;
		return qaws_curve_create_bezier_ex(out, &d, al);
	}
	case 1:
	{
		qaws_hermite_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.points = al_pts; d.derivatives = al_tan; d.point_count = 6; d.derivative_count = 6;
		return qaws_curve_create_hermite_ex(out, &d, al);
	}
	case 2:
	{
		qaws_catmull_rom_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.control_points = al_pts; d.control_point_count = 6;
		return qaws_curve_create_catmull_rom_ex(out, &d, al);
	}
	case 3:
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = al_pts; d.control_point_count = 6;
		return qaws_curve_create_bspline_ex(out, &d, al);
	}
	case 4:
	{
		static qaws_scalar const w[6] = { 1, 2, 0.5f, 1, 1.5f, 1 }, k[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = al_pts; d.control_point_count = 6; d.weights = w; d.weight_count = 6;
		d.knots = k; d.knot_count = 10;
		return qaws_curve_create_nurbs_ex(out, &d, al);
	}
	case 5:
	{
		static qaws_scalar const times[6] = { 0, 1, 2, 3.5f, 4, 5 };
		qaws_trajectory_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.key_positions = al_pts; d.key_count = 6; d.key_times = times; d.key_time_count = 6;
		return qaws_curve_create_trajectory_ex(out, &d, al);
	}
	case 6:
	{
		static qaws_scalar const w[4] = { 1, 0.7f, 1.3f, 1 };
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = al_pts; d.control_point_count = 4; d.weights = w; d.weight_count = 4;
		return qaws_curve_create_rational_bezier_ex(&d, al, out);
	}
	case 7:
	{
		qaws_yuksel_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.control_points = al_pts; d.control_point_count = 6; d.mode = QAWS_YUKSEL_MODE_HYBRID;
		return qaws_curve_create_yuksel_ex(&d, al, out);
	}
	case 8:
	{
		qaws_arc_segment seg[2];
		qaws_arc_desc d;
		memset(seg, 0, sizeof seg);
		seg[0].radius = 1; seg[0].angle_end = 1.5f;
		seg[1].center[0] = 0.2f; seg[1].radius = 0.8f; seg[1].angle_start = 1.5f; seg[1].angle_end = 3;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_2D; d.segments = seg; d.segment_count = 2;
		return qaws_curve_create_arc_ex(&d, al, out);
	}
	case 9:
	{
		qaws_polynomial_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.coefficients = al_pts; d.coefficient_count = 4;
		return qaws_curve_create_polynomial_ex(&d, al, out);
	}
	case 10:
	{
		qaws_clothoid_desc d;
		memset(&d, 0, sizeof d);
		d.start_curvature = 0.1f; d.end_curvature = 1.2f; d.length = 3;
		return qaws_curve_create_clothoid_ex(&d, al, out);
	}
	case 11:
	{
		/* two segments made with the same allocator; ownership moves to the composite only on success */
		qaws_curve* seg[2] = { NULL, NULL };
		qaws_composite_desc d;
		qaws_status st;
		if ((st = al_curve(0, al, &seg[0])) != QAWS_STATUS_OK)
			return st;
		if ((st = al_curve(9, al, &seg[1])) != QAWS_STATUS_OK)
		{
			qaws_curve_destroy(seg[0]);
			return st;
		}
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.segments = seg; d.segment_count = 2;
		if ((st = qaws_curve_create_composite_ex(&d, al, out)) != QAWS_STATUS_OK)
		{
			qaws_curve_destroy(seg[0]);
			qaws_curve_destroy(seg[1]);
		}
		return st;
	}
	case 12:
	{
		qaws_subdivision_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_3; d.control_points = al_pts; d.control_point_count = 6;
		d.refinement_levels = 3;
		return qaws_curve_create_subdivision_ex(&d, al, out);
	}
	default:
		return QAWS_STATUS_INVALID_ARGUMENT;
	}
}

static qaws_status al_surface(unsigned int c, qaws_allocator const* al, qaws_surface** out)
{
	static qaws_vec3 grid[16];
	static qaws_scalar const w[16] = { 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1 };
	unsigned int i;
	for (i = 0; i < 16; i++)
	{
		grid[i].x = (qaws_scalar)(i / 4);
		grid[i].y = (qaws_scalar)(i % 4);
		grid[i].z = (qaws_scalar)(0.3 * sin(1.7 * i));
	}
	*out = NULL;
	if (c == 0)
	{
		qaws_surface_bezier_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 3; d.v_degree = 3; d.control_points = grid; d.u_point_count = 4; d.v_point_count = 4;
		return qaws_surface_create_bezier_ex(out, &d, al);
	}
	if (c == 1)
	{
		qaws_surface_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 2; d.v_degree = 3; d.control_points = grid; d.u_point_count = 4; d.v_point_count = 4;
		return qaws_surface_create_bspline_ex(out, &d, al);
	}
	{
		qaws_surface_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 2; d.v_degree = 2; d.control_points = grid; d.u_point_count = 4; d.v_point_count = 4; d.weights = w;
		return qaws_surface_create_nurbs_ex(out, &d, al);
	}
}

static void test_constructors(void)
{
	static char const* const names[16] = { "Bezier", "Hermite", "Catmull-Rom", "B-spline", "NURBS", "trajectory", "rational Bezier", "Yuksel", "arc",
		"polynomial", "clothoid", "composite", "subdivision", "Bezier surface", "B-spline surface", "NURBS surface" };
	char report[1024] = "";
	unsigned int c, bad = 0;
	for (c = 0; c < 16; c++)
	{
		al_state s;
		qaws_allocator al;
		int k, ok = 0, wrong = 0, leaks = 0;
		al.alloc = al_alloc;
		al.dealloc = al_dealloc;
		al.user_data = &s;
		for (k = 0; k < 256 && !ok; k++)
		{
			qaws_curve* curve = NULL;
			qaws_surface* surface = NULL;
			qaws_status st;
			s.calls = 0;
			s.fail_at = k;
			s.live = 0;
			st = c < 13 ? al_curve(c, &al, &curve) : al_surface(c - 13, &al, &surface);
			if (st == QAWS_STATUS_OK)
			{
				ok = k + 1;
				/* (an evaluation, so a curve's lazily allocated state, if any, is counted too) */
				if (curve)
				{
					qaws_eval_result_3d r3;
					qaws_eval_result_2d r2;
					qaws_range r = qaws_curve_get_parameter_range(curve);
					if (qaws_curve_get_dimension(curve) == QAWS_DIMENSION_3D)
						qaws_curve_evaluate_3d(&r3, curve, (r.min_value + r.max_value) / 2, QAWS_EVAL_FLAG_POSITION);
					else
						qaws_curve_evaluate_2d(&r2, curve, (r.min_value + r.max_value) / 2, QAWS_EVAL_FLAG_POSITION);
				}
				qaws_curve_destroy(curve);
				qaws_surface_destroy(surface);
			}
			else
				wrong += st != QAWS_STATUS_ALLOCATION_FAILURE || curve || surface;
			leaks += s.live != 0;
		}
		if (!ok || wrong || leaks)
		{
			bad++;
			sprintf(report + strlen(report), " %s (succeeded %d, wrong %d, leaks %d);", names[c], ok, wrong, leaks);
		}
	}
	TEST_ASSERT(bad == 0, bad ? report : "every constructor: each allocation failure reported and cleaned up, nothing left after destroy");
	if (bad)
		printf("  %s\n", report);
}

/* inputs the surfaces borrow, made once without the allocator */
typedef struct al_inputs
{
	qaws_curve *c0, *c1, *d0, *d1, *mid, *path, *profile2, *uv[4];
	qaws_surface *plane_a, *plane_b;
} al_inputs;

static qaws_curve* al_line(qaws_dimension dim, double ax, double ay, double az, double bx, double by, double bz, double bump)
{
	qaws_vec3 p[4];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 4; i++)
	{
		double t = i / 3.0;
		p[i].x = (qaws_scalar)(ax + (bx - ax) * t);
		p[i].y = (qaws_scalar)(ay + (by - ay) * t);
		p[i].z = (qaws_scalar)(az + (bz - az) * t + (i == 1 || i == 2 ? bump : 0));
	}
	if (dim == QAWS_DIMENSION_2D)
		for (i = 0; i < 4; i++)
		{
			((qaws_vec2*)p)[i].x = p[i].x;
			((qaws_vec2*)p)[i].y = p[i].y;
		}
	memset(&d, 0, sizeof d);
	d.dimension = dim; d.degree = 3; d.control_points = p; d.control_point_count = 4;
	qaws_curve_create_bezier(&c, &d);
	return c;
}

static int al_inputs_make(al_inputs* in)
{
	qaws_surface_bilinear_desc a, b;
	memset(in, 0, sizeof *in);
	in->c0 = al_line(QAWS_DIMENSION_3D, 0, 0, 0, 1, 0, 0, 0.2);
	in->c1 = al_line(QAWS_DIMENSION_3D, 0, 1, 0, 1, 1, 0, -0.1);
	in->d0 = al_line(QAWS_DIMENSION_3D, 0, 0, 0, 0, 1, 0, 0.15);
	in->d1 = al_line(QAWS_DIMENSION_3D, 1, 0, 0, 1, 1, 0, -0.05);
	in->mid = al_line(QAWS_DIMENSION_3D, 0, 0.5, 0.3, 1, 0.5, 0.3, 0.1);
	in->path = al_line(QAWS_DIMENSION_3D, 0, 0, 0, 0, 0, 3, 0.5);
	in->profile2 = al_line(QAWS_DIMENSION_2D, 0.5, 0, 0, 1, 1, 0, 0);
	/* a closed square in the unit uv square */
	in->uv[0] = al_line(QAWS_DIMENSION_2D, 0.2, 0.2, 0, 0.8, 0.2, 0, 0);
	in->uv[1] = al_line(QAWS_DIMENSION_2D, 0.8, 0.2, 0, 0.8, 0.8, 0, 0);
	in->uv[2] = al_line(QAWS_DIMENSION_2D, 0.8, 0.8, 0, 0.2, 0.8, 0, 0);
	in->uv[3] = al_line(QAWS_DIMENSION_2D, 0.2, 0.8, 0, 0.2, 0.2, 0, 0);
	/* two planes crossing through each other's interiors (z = 0 and x = 0), so the fillet has an intersection to follow */
	memset(&a, 0, sizeof a);
	a.p00.x = -1; a.p00.y = -1; a.p10.x = 1; a.p10.y = -1; a.p01.x = -1; a.p01.y = 1; a.p11.x = 1; a.p11.y = 1;
	memset(&b, 0, sizeof b);
	b.p00.y = -1; b.p00.z = -1; b.p10.y = 1; b.p10.z = -1; b.p01.y = -1; b.p01.z = 1; b.p11.y = 1; b.p11.z = 1;
	/* (tilted, off the sampling grid: x = 0.137 + 0.21 z) */
	b.p00.x = b.p10.x = (qaws_scalar)(0.137 - 0.21);
	b.p01.x = b.p11.x = (qaws_scalar)(0.137 + 0.21);
	qaws_surface_create_bilinear(&in->plane_a, &a);
	qaws_surface_create_bilinear(&in->plane_b, &b);
	return in->c0 && in->c1 && in->d0 && in->d1 && in->mid && in->path && in->profile2 && in->uv[0] && in->uv[1] && in->uv[2] && in->uv[3]
		&& in->plane_a && in->plane_b;
}

static void al_inputs_free(al_inputs* in)
{
	unsigned int i;
	qaws_curve_destroy(in->c0); qaws_curve_destroy(in->c1); qaws_curve_destroy(in->d0); qaws_curve_destroy(in->d1);
	qaws_curve_destroy(in->mid); qaws_curve_destroy(in->path); qaws_curve_destroy(in->profile2);
	for (i = 0; i < 4; i++) qaws_curve_destroy(in->uv[i]);
	qaws_surface_destroy(in->plane_a);
	qaws_surface_destroy(in->plane_b);
}

#define AL_SURFACE_CASES 16

/* the other surface constructors; case 15 (the fillet with trims) also makes two curves */
static qaws_status al_surface2(unsigned int c, al_inputs const* in, qaws_allocator const* al, qaws_surface** out, qaws_curve** ta, qaws_curve** tb)
{
	*out = NULL;
	switch (c)
	{
	case 0:
	{
		qaws_surface_bilinear_desc d;
		memset(&d, 0, sizeof d);
		d.p10.x = 1; d.p01.y = 1; d.p11.x = 1; d.p11.y = 1; d.p11.z = 0.5f;
		return qaws_surface_create_bilinear_ex(&d, al, out);
	}
	case 1:
	{
		qaws_surface_biquadratic_desc d;
		unsigned int i;
		for (i = 0; i < 9; i++)
		{
			d.control_points[i].x = (qaws_scalar)(i / 3);
			d.control_points[i].y = (qaws_scalar)(i % 3);
			d.control_points[i].z = (qaws_scalar)(i == 4 ? 1 : 0);
		}
		return qaws_surface_create_biquadratic_ex(&d, al, out);
	}
	case 2:
	{
		qaws_surface_coons_desc d;
		d.c0 = in->c0; d.c1 = in->c1; d.d0 = in->d0; d.d1 = in->d1;
		return qaws_surface_create_coons_ex(&d, al, out);
	}
	case 3:
	{
		qaws_surface_extrusion_desc d;
		memset(&d, 0, sizeof d);
		d.profile = in->c0; d.direction.z = 1; d.length = 2;
		return qaws_surface_create_extrusion_ex(&d, al, out);
	}
	case 4:
	{
		qaws_surface_offset_desc d;
		d.base = in->plane_a; d.distance = (qaws_scalar)0.1;
		return qaws_surface_create_offset_ex(&d, al, out);
	}
	case 5:
	{
		qaws_surface_pipe_desc d;
		memset(&d, 0, sizeof d);
		d.path = in->path; d.radius_x = (qaws_scalar)0.2;
		return qaws_surface_create_pipe_ex(&d, al, out);
	}
	case 6:
	{
		qaws_surface_revolution_desc d;
		memset(&d, 0, sizeof d);
		d.profile = in->profile2; d.axis_direction.y = 1;
		return qaws_surface_create_revolution_ex(&d, al, out);
	}
	case 7:
	{
		qaws_surface_ruled_desc d;
		d.curve_a = in->c0; d.curve_b = in->c1;
		return qaws_surface_create_ruled_ex(&d, al, out);
	}
	case 8:
	{
		qaws_surface_swept_desc d;
		memset(&d, 0, sizeof d);
		d.path = in->path; d.profile = in->profile2;
		return qaws_surface_create_swept_ex(&d, al, out);
	}
	case 9:
	{
		qaws_curve const* sections[3];
		qaws_surface_loft_desc d;
		sections[0] = in->c0; sections[1] = in->mid; sections[2] = in->c1;
		memset(&d, 0, sizeof d);
		d.sections = sections; d.section_count = 3;
		return qaws_surface_create_loft_ex(&d, al, out);
	}
	case 10:
	{
		static qaws_scalar const p01[2] = { 0, 1 };
		qaws_curve const* uc[2];
		qaws_curve const* vc[2];
		qaws_surface_gordon_desc d;
		uc[0] = in->c0; uc[1] = in->c1; vc[0] = in->d0; vc[1] = in->d1;
		d.u_curves = uc; d.u_curve_count = 2; d.v_params = p01; d.v_curves = vc; d.v_curve_count = 2; d.u_params = p01;
		return qaws_surface_create_gordon_ex(&d, al, out);
	}
	case 11:
	{
		qaws_curve const* loop_curves[4];
		qaws_trim_loop loop;
		qaws_surface_trim_desc d;
		unsigned int i;
		for (i = 0; i < 4; i++) loop_curves[i] = in->uv[i];
		loop.curves = loop_curves; loop.curve_count = 4; loop.is_outer = 1;
		d.base = in->plane_a; d.loops = &loop; d.loop_count = 1;
		return qaws_surface_create_trimmed_ex(&d, al, out);
	}
	case 12:
	{
		static qaws_scalar const k[4] = { 0, 0, 1, 1 };
		qaws_tspline_control_point cps[4];
		qaws_surface_tspline_desc d;
		unsigned int i;
		for (i = 0; i < 4; i++)
		{
			cps[i].position.x = (qaws_scalar)(i & 1); cps[i].position.y = (qaws_scalar)(i >> 1); cps[i].position.z = (qaws_scalar)(i == 3);
			cps[i].weight = 1; cps[i].u_knots = k; cps[i].v_knots = k;
		}
		memset(&d, 0, sizeof d);
		d.control_points = cps; d.control_point_count = 4; d.u_degree = 1; d.v_degree = 1;
		d.u_range.max_value = 1; d.v_range.max_value = 1;
		return qaws_surface_create_tspline_ex(&d, al, out);
	}
	case 13:
	{
		/* a cube of 6 quads */
		static qaws_vec3 const v[8] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
		static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
		static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
		qaws_surface_subdiv_desc d;
		memset(&d, 0, sizeof d);
		d.vertices = v; d.vertex_count = 8; d.face_indices = f; d.face_sizes = fs; d.face_count = 6;
		d.scheme = QAWS_SUBDIV_CATMULL_CLARK; d.subdivision_level = 1;
		return qaws_surface_create_subdiv_ex(&d, al, out);
	}
	default:
	{
		qaws_surface_fillet_desc d;
		memset(&d, 0, sizeof d);
		d.surface_a = in->plane_a; d.surface_b = in->plane_b; d.radius = (qaws_scalar)0.3; d.sample_count = 32; d.arc_segments = 8;
		if (c == 14)
			return qaws_surface_create_fillet_ex(&d, al, out);
		return qaws_surface_create_fillet_with_trims_ex(&d, al, out, ta, tb);
	}
	}
}

static void test_surface_constructors(void)
{
	static char const* const names[AL_SURFACE_CASES] = { "bilinear", "biquadratic", "Coons", "extrusion", "offset", "pipe", "revolution", "ruled",
		"swept", "loft", "Gordon", "trimmed", "T-spline", "subdivision", "fillet", "fillet with trims" };
	char report[1024] = "";
	al_inputs in;
	unsigned int c, bad = 0;
	if (!al_inputs_make(&in))
	{
		TEST_ASSERT(0, "surface inputs made");
		al_inputs_free(&in);
		return;
	}
	for (c = 0; c < AL_SURFACE_CASES; c++)
	{
		al_state s;
		qaws_allocator al;
		int k, ok = 0, wrong = 0, leaks = 0;
		qaws_status last_wrong = QAWS_STATUS_OK;
		al.alloc = al_alloc;
		al.dealloc = al_dealloc;
		al.user_data = &s;
		for (k = 0; k < 4096 && !ok; k++)
		{
			qaws_surface* surface = NULL;
			qaws_curve *ta = NULL, *tb = NULL;
			qaws_status st;
			s.calls = 0;
			s.fail_at = k;
			s.live = 0;
			st = al_surface2(c, &in, &al, &surface, &ta, &tb);
			if (st == QAWS_STATUS_OK)
			{
				qaws_surface_eval_result r;
				ok = k + 1;
				/* an evaluation, so scratch past the stack buffers would be counted too */
				qaws_surface_evaluate(&r, surface, (qaws_scalar)0.4, (qaws_scalar)0.6, QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL);
				qaws_surface_destroy(surface);
				qaws_curve_destroy(ta);
				qaws_curve_destroy(tb);
			}
			else
			{
				wrong += st != QAWS_STATUS_ALLOCATION_FAILURE || surface || ta || tb;
				if (st != QAWS_STATUS_ALLOCATION_FAILURE)
					last_wrong = st;
			}
			leaks += s.live != 0;
		}
		if (!ok || wrong || leaks)
		{
			bad++;
			sprintf(report + strlen(report), " %s (succeeded %d, wrong %d (%s), leaks %d);", names[c], ok, wrong, qaws_status_to_string(last_wrong), leaks);
		}
	}
	TEST_ASSERT(bad == 0, bad ? report : "every surface constructor: each allocation failure reported and cleaned up, nothing left after destroy");
	if (bad)
		printf("  %s\n", report);
	al_inputs_free(&in);
}

/* derived curves: an operation's results come from its source curve's allocator */
static qaws_status al_derive(unsigned int c, qaws_curve const* src, qaws_curve const* src2, qaws_curve const* bez, qaws_curve** a,
	qaws_curve** b)
{
	switch (c)
	{
	case 0: return qaws_curve_reverse(a, src);
	case 1: return qaws_curve_extract(a, src, (qaws_scalar)0.4, (qaws_scalar)2.1);
	case 2: return qaws_curve_split(a, b, src, (qaws_scalar)1.3);
	case 3: return qaws_curve_elevate_degree(a, bez);
	case 4: return qaws_curve_convert_bspline_to_nurbs(a, src);
	case 5: return qaws_curve_reparameterize_arc_length(a, src, 64);
	default: return qaws_curve_match_arc_length(a, src2, src, 64);
	}
}

static void test_derived(void)
{
	static char const* const names[7] = { "reverse", "extract", "split", "degree elevation", "B-spline to NURBS", "arc-length curve",
		"length-matched curve" };
	char report[1024] = "";
	al_state s;
	qaws_allocator al;
	qaws_curve *src = NULL, *src2 = NULL, *bez = NULL;
	unsigned int c, bad = 0;
	int base;
	al.alloc = al_alloc;
	al.dealloc = al_dealloc;
	al.user_data = &s;
	s.calls = 0;
	s.fail_at = -1;
	s.live = 0;
	if (al_curve(3, &al, &src) != QAWS_STATUS_OK || al_curve(0, NULL, &src2) != QAWS_STATUS_OK || al_curve(0, &al, &bez) != QAWS_STATUS_OK)
	{
		TEST_ASSERT(0, "derived: sources made");
		qaws_curve_destroy(src);
		qaws_curve_destroy(src2);
		qaws_curve_destroy(bez);
		return;
	}
	base = s.live;
	for (c = 0; c < 7; c++)
	{
		int k, ok = 0, wrong = 0, leaks = 0, counted = 1;
		for (k = 0; k < 256 && !ok; k++)
		{
			qaws_curve *a = NULL, *b = NULL;
			qaws_status st;
			s.calls = 0;
			s.fail_at = k;
			st = al_derive(c, src, src2, bez, &a, &b);
			if (st == QAWS_STATUS_OK)
			{
				ok = k + 1;
				counted = s.live > base;      /* the results live in the source's allocator */
				qaws_curve_destroy(a);
				qaws_curve_destroy(b);
			}
			else
				wrong += st != QAWS_STATUS_ALLOCATION_FAILURE || a || b;
			leaks += s.live != base;
		}
		if (!ok || wrong || leaks || !counted)
		{
			bad++;
			sprintf(report + strlen(report), " %s (succeeded %d, wrong %d, leaks %d, results in the allocator %d);", names[c], ok, wrong, leaks, counted);
		}
	}
	s.fail_at = -1;
	qaws_curve_destroy(src);
	qaws_curve_destroy(src2);
	qaws_curve_destroy(bez);
	TEST_ASSERT(bad == 0 && s.live == 0, bad ? report : "derived curves in the source's allocator; failures reported, nothing left");
	if (bad)
		printf("  %s\n", report);
}

int test_146_allocator_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 146: the allocator through every constructor\n");
	test_constructors();
	test_surface_constructors();
	test_derived();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
