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
		return qaws_curve_create_bezier_ex(&d, al, out);
	}
	case 1:
	{
		qaws_hermite_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.points = al_pts; d.derivatives = al_tan; d.point_count = 6; d.derivative_count = 6;
		return qaws_curve_create_hermite_ex(&d, al, out);
	}
	case 2:
	{
		qaws_catmull_rom_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.control_points = al_pts; d.control_point_count = 6;
		return qaws_curve_create_catmull_rom_ex(&d, al, out);
	}
	case 3:
	{
		qaws_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = al_pts; d.control_point_count = 6;
		return qaws_curve_create_bspline_ex(&d, al, out);
	}
	case 4:
	{
		static qaws_scalar const w[6] = { 1, 2, 0.5f, 1, 1.5f, 1 }, k[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };
		qaws_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.control_points = al_pts; d.control_point_count = 6; d.weights = w; d.weight_count = 6;
		d.knots = k; d.knot_count = 10;
		return qaws_curve_create_nurbs_ex(&d, al, out);
	}
	case 5:
	{
		static qaws_scalar const times[6] = { 0, 1, 2, 3.5f, 4, 5 };
		qaws_trajectory_desc d;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_3D; d.degree = 3; d.key_positions = al_pts; d.key_count = 6; d.key_times = times; d.key_time_count = 6;
		return qaws_curve_create_trajectory_ex(&d, al, out);
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
		return qaws_surface_create_bezier_ex(&d, al, out);
	}
	if (c == 1)
	{
		qaws_surface_bspline_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 2; d.v_degree = 3; d.control_points = grid; d.u_point_count = 4; d.v_point_count = 4;
		return qaws_surface_create_bspline_ex(&d, al, out);
	}
	{
		qaws_surface_nurbs_desc d;
		memset(&d, 0, sizeof d);
		d.u_degree = 2; d.v_degree = 2; d.control_points = grid; d.u_point_count = 4; d.v_point_count = 4; d.weights = w;
		return qaws_surface_create_nurbs_ex(&d, al, out);
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
						qaws_curve_evaluate_3d(curve, (r.min_value + r.max_value) / 2, QAWS_EVAL_FLAG_POSITION, &r3);
					else
						qaws_curve_evaluate_2d(curve, (r.min_value + r.max_value) / 2, QAWS_EVAL_FLAG_POSITION, &r2);
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

int test_146_allocator_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 146: the allocator through every constructor\n");
	test_constructors();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
