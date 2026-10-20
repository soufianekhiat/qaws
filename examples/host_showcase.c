/*
 * host_showcase.c - features for host applications
 *
 * Writes SVG figures into ./showcase/:
 *   host1_frames.svg   a ribbon swept along an S-shaped space curve with a
 *                      straight middle: by the Frenet frame (undefined on the
 *                      straight piece, flipping at the inflection) and by
 *                      rotation-minimizing frames (qaws_curve_compute_
 *                      rotation_minimizing_frames_3d: continuous, twist-free)
 *   host2_arc_length.svg  points at uniform parameters against points at equal
 *                      arc length (qaws_curve_compute_parameter_at_arc_length)
 *   host3_closest.svg  the point of a spiked circle closest to its center: the
 *                      nearest of 32 samples against the certified search
 *                      (qaws_curve_find_closest_point_certified_2d)
 *   host4_tessellation.svg  an S whose middle lies on its chord: the adaptive
 *                      sampler's midpoint test against qaws_curve_tessellate_2d
 *                      at chord and angle tolerances
 *   host5_surface_closest.svg  the point of a dented spherical cap closest to
 *                      its center: the nearest of an 8 x 8 grid against
 *                      qaws_surface_find_closest_point_certified
 *   host6_surface_tessellation.svg  a B-spline flat on half its net: the
 *                      quadtree with an absolute flatness threshold against
 *                      qaws_surface_tessellate_adaptive (chord and angle)
 *   host7_pipe.svg     a pipe around a twisted path: rotation-minimizing frames,
 *                      lit by the normal of its analytic derivatives
 *   host8_parameterization.svg  interpolating fits through uneven points by
 *                      uniform, chordal and centripetal parameters
 */

#include "qaws.h"
#include "example_svg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HS_SAMPLES 160

static void hs_quad(svg* s, double const* xy, char const* fill, double opacity)
{
	fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" fill-opacity=\"%.2f\" stroke=\"%s\" "
		"stroke-width=\"0.4\"/>\n", xy[0], xy[1], xy[2], xy[3], xy[4], xy[5], xy[6], xy[7], fill, opacity, fill);
}

/* the ribbon P +- w N over the curve's samples: frames from Frenet (ok[i] 0 where undefined) or the rotation-minimizing ones */
static void hs_ribbon(svg* s, projection const* pr, qaws_curve const* c, int frenet, double w, char const* front, char const* back)
{
	qaws_scalar t[HS_SAMPLES], fr[HS_SAMPLES * 9];
	qaws_vec3 pos[HS_SAMPLES], nrm[HS_SAMPLES];
	int ok[HS_SAMPLES], i;
	double line[HS_SAMPLES * 2];
	qaws_range r = qaws_curve_get_parameter_range(c);
	for (i = 0; i < HS_SAMPLES; i++)
		t[i] = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / (HS_SAMPLES - 1);
	if (!frenet)
		qaws_curve_compute_rotation_minimizing_frames_3d(c, t, HS_SAMPLES, NULL, fr);
	for (i = 0; i < HS_SAMPLES; i++)
	{
		qaws_eval_result_3d e;
		qaws_curve_evaluate_3d(&e, c, t[i], QAWS_EVAL_FLAG_POSITION);
		pos[i] = e.position;
		if (frenet)
		{
			qaws_vec3 tg, bn;
			ok[i] = qaws_curve_compute_frenet_frame_3d(&tg, &nrm[i], &bn, c, t[i]) == QAWS_STATUS_OK
				&& isfinite((double)nrm[i].x) && (double)nrm[i].x * nrm[i].x + (double)nrm[i].y * nrm[i].y + (double)nrm[i].z * nrm[i].z > 0.5;
		}
		else
		{
			nrm[i].x = fr[9 * i + 3]; nrm[i].y = fr[9 * i + 4]; nrm[i].z = fr[9 * i + 5];
			ok[i] = 1;
		}
		project(pr, pos[i].x, pos[i].y, pos[i].z, &line[2 * i], &line[2 * i + 1]);
	}
	for (i = 0; i + 1 < HS_SAMPLES; i++)
	{
		double q[8], up;
		if (!ok[i] || !ok[i + 1])
			continue;
		project(pr, pos[i].x + w * nrm[i].x, pos[i].y + w * nrm[i].y, pos[i].z + w * nrm[i].z, &q[0], &q[1]);
		project(pr, pos[i + 1].x + w * nrm[i + 1].x, pos[i + 1].y + w * nrm[i + 1].y, pos[i + 1].z + w * nrm[i + 1].z, &q[2], &q[3]);
		project(pr, pos[i + 1].x - w * nrm[i + 1].x, pos[i + 1].y - w * nrm[i + 1].y, pos[i + 1].z - w * nrm[i + 1].z, &q[4], &q[5]);
		project(pr, pos[i].x - w * nrm[i].x, pos[i].y - w * nrm[i].y, pos[i].z - w * nrm[i].z, &q[6], &q[7]);
		/* which face shows: the screen orientation of the quad */
		up = (q[2] - q[0]) * (q[7] - q[1]) - (q[3] - q[1]) * (q[6] - q[0]);
		hs_quad(s, q, up > 0 ? front : back, 0.85);
	}
	svg_polyline(s, line, HS_SAMPLES, "#1b1f24", 1.6, 1, 0);
}

static void figure_frames(void)
{
	/* an S in space with a straight middle: two cubic arcs joined by a line */
	qaws_vec3 p[10] = { { -4, -1, 0 }, { -3.5f, 1, 0.3f }, { -2, 1.5f, 0.6f }, { -1, 0.5f, 0.9f }, { 0, 0, 1 }, { 1, -0.5f, 1.1f },
		{ 2, -1.5f, 1.4f }, { 3.5f, -1, 1.7f }, { 4, 1, 2 }, { 3, 2, 2.5f } };
	qaws_scalar knots[14] = { 0, 0, 0, 0, 1, 2, 2, 2, 3, 4, 4, 4, 4, 4 };
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	projection left = { 330, 400, 66, 1.0 }, right = { 930, 400, 66, 1.0 };
	viewport vl = { 24, 80, 576, 620, 0, 1, 0, 1 }, vr = { 624, 80, 576, 620, 0, 1, 0, 1 };
	svg s;
	/* collinear points 3..5 with a triple knot make the middle straight */
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 10;
	d.knots = knots;
	d.knot_count = 14;
	if (qaws_curve_create_bspline(&c, &d) != QAWS_STATUS_OK)
	{
		printf("host1: curve not created\n");
		return;
	}
	if (!svg_open(&s, "showcase/host1_frames.svg", 1224, 730, "Frames along a curve",
		"a ribbon swept by the normal: Frenet (left) is undefined on the straight piece and flips at the inflection; "
		"rotation-minimizing frames (right) follow without twisting"))
		return;
	svg_panel(&s, &vl, "Frenet frame");
	svg_panel(&s, &vr, "rotation-minimizing frames (double reflection)");
	hs_ribbon(&s, &left, c, 1, 0.35, "#f4a261", "#9c4f1a");
	hs_ribbon(&s, &right, c, 0, 0.35, "#2a9d8f", "#145c54");
	svg_text(&s, 40, 690, 13, "#57606a", "start", "gaps: no Frenet normal (zero curvature); colour change: the ribbon turned over");
	svg_text(&s, 640, 690, 13, "#57606a", "start", "qaws_curve_compute_rotation_minimizing_frames_3d, 160 samples, fourth-order accurate");
	svg_close(&s);
	qaws_curve_destroy(c);
	printf("wrote showcase/host1_frames.svg\n");
}

#define HS_DOTS 25

/* the curve and HS_DOTS points on it: at uniform parameters, or at equal lengths (qaws_curve_compute_parameter_at_arc_length) */
static void hs_dots(svg* s, viewport const* v, qaws_curve const* c, int by_length, char const* color)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	qaws_scalar total = 0, t, t_prev = r.min_value;
	double line[HS_SAMPLES * 2], xy[HS_DOTS * 2], dmin = 1e300, dmax = 0;
	char note[160];
	int i;
	for (i = 0; i < HS_SAMPLES; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(&e, c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / (HS_SAMPLES - 1), QAWS_EVAL_FLAG_POSITION);
		line[2 * i] = vx(v, e.position.x);
		line[2 * i + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, line, HS_SAMPLES, "#1b1f24", 1.6, 1, 0);
	qaws_curve_compute_arc_length(&total, c, r.min_value, r.max_value);
	for (i = 0; i < HS_DOTS; i++)
	{
		qaws_eval_result_2d e;
		t = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / (HS_DOTS - 1);
		if (by_length)
			qaws_curve_compute_parameter_at_arc_length(c, r.min_value, total * (qaws_scalar)i / (HS_DOTS - 1), &t);
		qaws_curve_evaluate_2d(&e, c, t, QAWS_EVAL_FLAG_POSITION);
		xy[2 * i] = e.position.x;
		xy[2 * i + 1] = e.position.y;
		if (i > 0)
		{
			qaws_scalar arc = 0;
			double d;
			qaws_curve_compute_arc_length(&arc, c, t_prev, t);
			d = arc;
			dmin = d < dmin ? d : dmin;
			dmax = d > dmax ? d : dmax;
		}
		t_prev = t;
	}
	for (i = 0; i < HS_DOTS; i++)
		svg_circle(s, vx(v, xy[2 * i]), vy(v, xy[2 * i + 1]), 5, color, "#1b1f24");
	snprintf(note, sizeof note, "arc lengths between neighbours: %.6f to %.6f", dmin, dmax);
	svg_text(s, v->x0 + 16, v->y0 + v->h - 16, 13, "#57606a", "start", note);
}

static void figure_arc_length(void)
{
	/* a cubic B-spline whose knots crowd the parameter into the first turn */
	qaws_vec2 p[8] = { { 0, 0 }, { 1, 3 }, { 3, 3.5f }, { 4, 0.5f }, { 6, -0.5f }, { 7.5f, 2.5f }, { 9, 3 }, { 10, 0 } };
	qaws_scalar knots[12] = { 0, 0, 0, 0, 0.05f, 0.1f, 0.6f, 0.95f, 1, 1, 1, 1 };
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	viewport vl = { 24, 80, 576, 560, -0.8, 10.8, -1.9, 4.6 }, vr = { 624, 80, 576, 560, -0.8, 10.8, -1.9, 4.6 };
	svg s;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 8;
	d.knots = knots;
	d.knot_count = 12;
	if (qaws_curve_create_bspline(&c, &d) != QAWS_STATUS_OK)
	{
		printf("host2: curve not created\n");
		return;
	}
	if (!svg_open(&s, "showcase/host2_arc_length.svg", 1224, 670, "Points at equal arc length",
		"25 points on a B-spline with uneven parameter speed: at uniform parameters (left) they bunch and spread; "
		"at equal lengths (right) they are evenly spaced"))
		return;
	svg_panel(&s, &vl, "uniform parameters");
	svg_panel(&s, &vr, "equal arc length (safeguarded Newton)");
	hs_dots(&s, &vl, c, 0, "#f4a261");
	hs_dots(&s, &vr, c, 1, "#2a9d8f");
	svg_text(&s, 640, 662, 13, "#57606a", "start", "qaws_curve_compute_parameter_at_arc_length over adaptive Gauss-Kronrod per span");
	svg_close(&s);
	qaws_curve_destroy(c);
	printf("wrote showcase/host2_arc_length.svg\n");
}

#define HS_RING 200

static void hs_ring_curve(svg* s, viewport const* v, qaws_curve const* c)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double line[1200 * 2];
	int i;
	for (i = 0; i < 1200; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(&e, c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 1199, QAWS_EVAL_FLAG_POSITION);
		line[2 * i] = vx(v, e.position.x);
		line[2 * i + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, line, 1200, "#1b1f24", 1.6, 1, 0);
	svg_circle(s, vx(v, 0), vy(v, 0), 5, "#d1242f", "#1b1f24");
}

static void figure_closest(void)
{
	/* a circle of radius 5 as a cubic B-spline, one point pulled to radius 0.5: a spike between the 32 samples */
	static qaws_vec2 p[HS_RING];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	qaws_range r;
	qaws_vec2 q = { 0, 0 };
	qaws_scalar t = 0, dist = 0, gap = 0, ts = 0;
	double best = 1e300;
	viewport vl = { 24, 80, 576, 576, -6, 6, -6, 6 }, vr = { 624, 80, 576, 576, -6, 6, -6, 6 };
	qaws_eval_result_2d e;
	char note[200];
	svg s;
	int i;
	for (i = 0; i < HS_RING; i++)
	{
		double rad = i == 67 ? 0.5 : 5.0, a = 6.283185307179586 * i / (HS_RING - 1);
		p[i].x = (qaws_scalar)(rad * cos(a));
		p[i].y = (qaws_scalar)(rad * sin(a));
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = HS_RING;
	if (qaws_curve_create_bspline(&c, &d) != QAWS_STATUS_OK)
	{
		printf("host3: curve not created\n");
		return;
	}
	if (!svg_open(&s, "showcase/host3_closest.svg", 1224, 700, "Certified closest point",
		"the point of the curve closest to the circle's center (red): a spike lies between the 32 samples a sampled search starts from")) 
		return;
	svg_panel(&s, &vl, "nearest of 32 samples (Newton stays on the circle: equal distances)");
	svg_panel(&s, &vr, "certified: branch and bound over the Bezier spans");
	hs_ring_curve(&s, &vl, c);
	hs_ring_curve(&s, &vr, c);
	r = qaws_curve_get_parameter_range(c);
	for (i = 0; i < 32; i++)
	{
		double dd;
		qaws_curve_evaluate_2d(&e, c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 31, QAWS_EVAL_FLAG_POSITION);
		dd = hypot(e.position.x, e.position.y);
		svg_circle(&s, vx(&vl, e.position.x), vy(&vl, e.position.y), 3.5, "#f4a261", "#1b1f24");
		if (dd < best)
		{
			best = dd;
			ts = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 31;
		}
	}
	qaws_curve_evaluate_2d(&e, c, ts, QAWS_EVAL_FLAG_POSITION);
	svg_line(&s, vx(&vl, 0), vy(&vl, 0), vx(&vl, e.position.x), vy(&vl, e.position.y), "#f4a261", 2.5, 1);
	snprintf(note, sizeof note, "distance %.6f", best);
	svg_text(&s, vl.x0 + 16, vl.y0 + vl.h - 16, 13, "#57606a", "start", note);
	qaws_curve_find_closest_point_certified_2d(c, q, &t, &dist, &gap);
	qaws_curve_evaluate_2d(&e, c, t, QAWS_EVAL_FLAG_POSITION);
	svg_line(&s, vx(&vr, 0), vy(&vr, 0), vx(&vr, e.position.x), vy(&vr, e.position.y), "#2a9d8f", 2.5, 1);
	svg_circle(&s, vx(&vr, e.position.x), vy(&vr, e.position.y), 5, "#2a9d8f", "#1b1f24");
	snprintf(note, sizeof note, "distance %.12f, gap %.1e (the minimum is in [d - gap, d])", (double)dist, (double)gap);
	svg_text(&s, vr.x0 + 16, vr.y0 + vr.h - 16, 13, "#57606a", "start", note);
	svg_text(&s, 640, 692, 13, "#57606a", "start", "qaws_curve_find_closest_point_certified_2d: convex-hull and fat-line bounds, de Casteljau halving");
	svg_close(&s);
	qaws_curve_destroy(c);
	printf("wrote showcase/host3_closest.svg\n");
}

#define HS_TESS 512

/* the curve, its polyline and vertices */
static void hs_polyline(svg* s, viewport const* v, qaws_curve const* c, qaws_vec2 const* pts, unsigned int n, char const* color, char const* label)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	double line[600 * 2], poly[HS_TESS * 2];
	char note[160];
	unsigned int i;
	for (i = 0; i < 600; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(&e, c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 599, QAWS_EVAL_FLAG_POSITION);
		line[2 * i] = vx(v, e.position.x);
		line[2 * i + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, line, 600, "#8c959f", 6, 0.35, 0);
	for (i = 0; i < n; i++)
	{
		poly[2 * i] = vx(v, pts[i].x);
		poly[2 * i + 1] = vy(v, pts[i].y);
	}
	svg_polyline(s, poly, (int)n, color, 2, 1, 0);
	for (i = 0; i < n; i++)
		svg_circle(s, poly[2 * i], poly[2 * i + 1], 3.5, color, "#1b1f24");
	snprintf(note, sizeof note, "%s: %u vertices", label, n);
	svg_text(s, v->x0 + 16, v->y0 + v->h - 16, 13, "#57606a", "start", note);
}

static void figure_tessellation(void)
{
	/* an S whose middle (4.5, 0) lies on its chord */
	qaws_vec2 p[4] = { { 0, 0 }, { 3, 4 }, { 6, -4 }, { 9, 0 } };
	static qaws_scalar t[HS_TESS];
	static qaws_vec2 pts[HS_TESS], old_pts[HS_TESS];
	qaws_bezier_desc d;
	qaws_sampling_desc sd;
	qaws_curve_tessellation_desc td;
	qaws_curve* c = NULL;
	unsigned int n = 0, old_n = 0;
	viewport vl = { 24, 80, 576, 520, -0.8, 9.8, -2.4, 2.8 }, vr = { 624, 80, 576, 520, -0.8, 9.8, -2.4, 2.8 };
	svg s;
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 4;
	if (qaws_curve_create_bezier(&c, &d) != QAWS_STATUS_OK)
	{
		printf("host4: curve not created\n");
		return;
	}
	/* the adaptive sampler: a midpoint test from one segment */
	memset(&sd, 0, sizeof sd);
	sd.sample_count = 2;
	sd.error_tolerance = (qaws_scalar)0.02;
	sd.use_adaptive_sampling = 1;
	qaws_curve_sample_2d(old_pts, &old_n, c, HS_TESS, &sd);
	memset(&td, 0, sizeof td);
	td.chord_tolerance = (qaws_scalar)0.02;
	td.angle_tolerance = (qaws_scalar)0.35;
	qaws_curve_tessellate_2d(c, &td, NULL, NULL, 0, &n);
	if (n > HS_TESS || qaws_curve_tessellate_2d(c, &td, t, pts, HS_TESS, &n) != QAWS_STATUS_OK)
	{
		printf("host4: tessellation failed\n");
		qaws_curve_destroy(c);
		return;
	}
	if (!svg_open(&s, "showcase/host4_tessellation.svg", 1224, 640, "Adaptive curve tessellation",
		"the same chord tolerance 0.02 on a cubic Bezier S (grey): a midpoint test accepts a segment whose middle lies on its chord"))
		return;
	svg_panel(&s, &vl, "adaptive sampling (midpoint test, from one segment)");
	svg_panel(&s, &vr, "qaws_curve_tessellate_2d (chord 0.02, angle 0.35 rad)");
	hs_polyline(&s, &vl, c, old_pts, old_n, "#f4a261", "qaws_curve_sample_2d, use_adaptive_sampling");
	hs_polyline(&s, &vr, c, pts, n, "#2a9d8f", "tests at 1/4, 1/2, 3/4; counts first, caller buffers");
	svg_close(&s);
	qaws_curve_destroy(c);
	printf("wrote showcase/host4_tessellation.svg\n");
}

#define HS_CAP 40

static double hs_radius(void const* user, qaws_scalar u, qaws_scalar v)
{
	qaws_surface_eval_result r;
	qaws_surface_evaluate(&r, (qaws_surface const*)user, u, v, QAWS_SURFACE_EVAL_POSITION);
	return sqrt((double)r.position.x * r.position.x + (double)r.position.y * r.position.y + (double)r.position.z * r.position.z);
}

static void hs_mark(svg* s, projection const* pr, qaws_vec3 p, char const* color, char const* label)
{
	double x0, y0, x1, y1;
	project(pr, 0, 0, 0, &x0, &y0);
	project(pr, p.x, p.y, p.z, &x1, &y1);
	svg_line(s, x0, y0, x1, y1, color, 2.5, 1);
	svg_circle(s, x1, y1, 5, color, "#1b1f24");
	svg_circle(s, x0, y0, 5, "#d1242f", "#1b1f24");
	svg_text(s, x1 + 8, y1 - 8, 13, "#1b1f24", "start", label);
}

static void figure_surface_closest(void)
{
	/* a spherical cap of radius 5 around the origin, one net point pulled to radius 0.5: a dent between the 8 x 8 samples */
	static qaws_vec3 p[HS_CAP * HS_CAP];
	qaws_surface_bspline_desc d;
	qaws_surface* c = NULL;
	qaws_vec3 q = { 0, 0, 0 }, cl, near_pt = { 0, 0, 0 };
	qaws_scalar u, v, dist, gap;
	projection left = { 312, 400, 37, 1.0 }, right = { 912, 400, 37, 1.0 };
	viewport vl = { 24, 80, 576, 620, 0, 1, 0, 1 }, vr = { 624, 80, 576, 620, 0, 1, 0, 1 };
	double best = 1e300, lx, ly;
	char note[200];
	svg s;
	unsigned int i, j;
	for (i = 0; i < HS_CAP; i++)
		for (j = 0; j < HS_CAP; j++)
		{
			double th = 0.3 + 2.5 * i / (HS_CAP - 1), ph = 3.0 * j / (HS_CAP - 1), r = i == 14 && j == 14 ? 0.5 : 5.0;
			p[i * HS_CAP + j].x = (qaws_scalar)(r * sin(th) * cos(ph));
			p[i * HS_CAP + j].y = (qaws_scalar)(r * sin(th) * sin(ph));
			p[i * HS_CAP + j].z = (qaws_scalar)(r * cos(th));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = p;
	d.u_point_count = HS_CAP;
	d.v_point_count = HS_CAP;
	if (qaws_surface_create_bspline(&c, &d) != QAWS_STATUS_OK)
	{
		printf("host5: surface not created\n");
		return;
	}
	if (!svg_open(&s, "showcase/host5_surface_closest.svg", 1224, 760, "Certified closest point on a surface",
		"a bicubic B-spline cap of radius 5 around the red point, coloured by distance to it, with a dent between the samples of an 8 x 8 grid"))
		return;
	svg_panel(&s, &vl, "nearest of the 8 x 8 grid (Newton stays on the sphere)");
	svg_panel(&s, &vr, "certified: branch and bound over the Bezier patches");
	draw_quads(&s, &left, c, 110, 110, hs_radius, c, 2.9, 5.0);
	draw_quads(&s, &right, c, 110, 110, hs_radius, c, 2.9, 5.0);
	for (i = 0; i < 8; i++)
		for (j = 0; j < 8; j++)
		{
			qaws_surface_eval_result r;
			double dd;
			qaws_surface_evaluate(&r, c, (qaws_scalar)(i / 7.0), (qaws_scalar)(j / 7.0), QAWS_SURFACE_EVAL_POSITION);
			dd = sqrt((double)r.position.x * r.position.x + (double)r.position.y * r.position.y + (double)r.position.z * r.position.z);
			project(&left, r.position.x, r.position.y, r.position.z, &lx, &ly);
			svg_circle(&s, lx, ly, 3, "#f4a261", "#1b1f24");
			if (dd < best)
			{
				best = dd;
				near_pt = r.position;
			}
		}
	snprintf(note, sizeof note, "%.6f", best);
	hs_mark(&s, &left, near_pt, "#f4a261", note);
	qaws_surface_find_closest_point_certified(c, q, &u, &v, &cl, &dist, &gap);
	snprintf(note, sizeof note, "%.9f (gap %.0e)", (double)dist, (double)gap);
	hs_mark(&s, &right, cl, "#2a9d8f", note);
	svg_colorbar(&s, 640, 716, 300, 12, "2.9", "5.0");
	svg_text(&s, 960, 727, 13, "#57606a", "start", "distance to the red point");
	svg_close(&s);
	qaws_surface_destroy(c);
	printf("wrote showcase/host5_surface_closest.svg\n");
}

#define HS_TV 120000
#define HS_TI 360000

/* triangles shaded by their normal against a light, back to front */
static void hs_mesh(svg* s, projection const* pr, qaws_tessellation_vertex const* v, unsigned int const* idx, unsigned int ni)
{
	static unsigned int order[HS_TI / 3];
	static double depth[HS_TI / 3];
	unsigned int n = ni / 3, t, k;
	for (t = 0; t < n; t++)
	{
		order[t] = t;
		depth[t] = 0;
		for (k = 0; k < 3; k++)
			depth[t] += v[idx[3 * t + k]].position.x + v[idx[3 * t + k]].position.y;
	}
	/* shell sort by depth, far first */
	for (k = n / 2; k > 0; k /= 2)
		for (t = k; t < n; t++)
		{
			unsigned int o = order[t], m = t;
			while (m >= k && depth[order[m - k]] > depth[o])
			{
				order[m] = order[m - k];
				m -= k;
			}
			order[m] = o;
		}
	for (t = 0; t < n; t++)
	{
		qaws_vec3 const* p[3];
		double x[3], y[3], e1[3], e2[3], nx, ny, nz, l, lit;
		int g;
		for (k = 0; k < 3; k++)
		{
			p[k] = &v[idx[3 * order[t] + k]].position;
			project(pr, p[k]->x, p[k]->y, p[k]->z, &x[k], &y[k]);
		}
		e1[0] = p[1]->x - p[0]->x; e1[1] = p[1]->y - p[0]->y; e1[2] = p[1]->z - p[0]->z;
		e2[0] = p[2]->x - p[0]->x; e2[1] = p[2]->y - p[0]->y; e2[2] = p[2]->z - p[0]->z;
		nx = e1[1] * e2[2] - e1[2] * e2[1];
		ny = e1[2] * e2[0] - e1[0] * e2[2];
		nz = e1[0] * e2[1] - e1[1] * e2[0];
		l = sqrt(nx * nx + ny * ny + nz * nz);
		lit = l > 0 ? fabs(-0.6 * nx + 0.35 * ny + 0.72 * nz) / l : 1;
		g = (int)(40 + 215 * lit * lit);
		fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"rgb(%d,%d,%d)\" stroke=\"#1b1f24\" stroke-width=\"0.25\" "
			"stroke-opacity=\"0.7\"/>\n", x[0], y[0], x[1], y[1], x[2], y[2], g / 3, g * 2 / 3, g);
	}
}

static void figure_tessellation_surface(void)
{
	/* flat for the first half of the net, waves on the other */
	static qaws_vec3 p[24 * 24];
	static qaws_tessellation_vertex v[HS_TV];
	static unsigned int idx[HS_TI];
	qaws_surface_bspline_desc d;
	qaws_tessellation_desc od;
	qaws_surface_tessellation_desc td;
	qaws_surface* c = NULL;
	unsigned int i, j, nv = 0, ni = 0;
	projection left = { 312, 230, 12.5, 1.0 }, right = { 912, 230, 12.5, 1.0 };
	viewport vl = { 24, 80, 576, 520, 0, 1, 0, 1 }, vr = { 624, 80, 576, 520, 0, 1, 0, 1 };
	char note[200];
	svg s;
	for (i = 0; i < 24; i++)
		for (j = 0; j < 24; j++)
		{
			p[i * 24 + j].x = (qaws_scalar)i;
			p[i * 24 + j].y = (qaws_scalar)j;
			p[i * 24 + j].z = (qaws_scalar)(i < 12 ? 0 : 1.5 * sin(1.3 * j) * cos(0.9 * i));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = p;
	d.u_point_count = 24;
	d.v_point_count = 24;
	if (qaws_surface_create_bspline(&c, &d) != QAWS_STATUS_OK)
	{
		printf("host6: surface not created\n");
		return;
	}
	if (!svg_open(&s, "showcase/host6_surface_tessellation.svg", 1224, 640, "Adaptive surface tessellation",
		"a bicubic B-spline flat on half its net: a quadtree with an absolute flatness threshold against chord and angle tolerances"))
		return;
	svg_panel(&s, &vl, "qaws_surface_tessellate (flatness 0.02, depth 6)");
	svg_panel(&s, &vr, "qaws_surface_tessellate_adaptive (chord 0.02, angle 0.5)");
	memset(&od, 0, sizeof od);
	od.max_depth = 6;
	od.curvature_threshold = (qaws_scalar)0.02;
	qaws_surface_tessellate(v, &nv, idx, &ni, c, HS_TV, HS_TI, &od);
	hs_mesh(&s, &left, v, idx, ni);
	snprintf(note, sizeof note, "%u triangles", ni / 3);
	svg_text(&s, vl.x0 + 16, vl.y0 + vl.h - 16, 13, "#57606a", "start", note);
	memset(&td, 0, sizeof td);
	td.chord_tolerance = (qaws_scalar)0.02;
	td.angle_tolerance = (qaws_scalar)0.5;
	qaws_surface_tessellate_adaptive(c, &td, NULL, 0, &nv, NULL, 0, &ni);
	if (nv <= HS_TV && ni <= HS_TI && qaws_surface_tessellate_adaptive(c, &td, v, HS_TV, &nv, idx, HS_TI, &ni) == QAWS_STATUS_OK)
	{
		hs_mesh(&s, &right, v, idx, ni);
		snprintf(note, sizeof note, "%u triangles, no cracks, counts first into caller buffers", ni / 3);
		svg_text(&s, vr.x0 + 16, vr.y0 + vr.h - 16, 13, "#57606a", "start", note);
	}
	svg_close(&s);
	qaws_surface_destroy(c);
	printf("wrote showcase/host6_surface_tessellation.svg\n");
}

#define HS_PU 260
#define HS_PV 24

/* a pipe as quads back to front, lit by its analytic normal, striped along v (a twisting frame would wind the stripes) */
static void hs_pipe(svg* s, projection const* pr, qaws_surface const* pipe)
{
	static int order[HS_PU * HS_PV];
	static double depth[HS_PU * HS_PV];
	int n = HS_PU * HS_PV, k, m, gap;
	for (k = 0; k < n; k++)
	{
		qaws_surface_eval_result r;
		qaws_surface_evaluate(&r, pipe, (qaws_scalar)((k / HS_PV + 0.5) / HS_PU), (qaws_scalar)((k % HS_PV + 0.5) / HS_PV), QAWS_SURFACE_EVAL_POSITION);
		depth[k] = r.position.x + r.position.y - 0.4 * r.position.z;
		order[k] = k;
	}
	for (gap = n / 2; gap > 0; gap /= 2)
		for (k = gap; k < n; k++)
		{
			int o = order[k];
			for (m = k; m >= gap && depth[order[m - gap]] > depth[o]; m -= gap)
				order[m] = order[m - gap];
			order[m] = o;
		}
	for (k = 0; k < n; k++)
	{
		int i = order[k] / HS_PV, j = order[k] % HS_PV, c;
		double x[4], y[4], lit;
		qaws_surface_eval_result r[4], mid;
		int du[4] = { 0, 1, 1, 0 }, dv[4] = { 0, 0, 1, 1 };
		for (c = 0; c < 4; c++)
		{
			qaws_surface_evaluate(&r[c], pipe, (qaws_scalar)((double)(i + du[c]) / HS_PU), (qaws_scalar)((double)(j + dv[c]) / HS_PV), QAWS_SURFACE_EVAL_POSITION);
			project(pr, r[c].position.x, r[c].position.y, r[c].position.z, &x[c], &y[c]);
		}
		qaws_surface_evaluate(&mid, pipe, (qaws_scalar)((i + 0.5) / HS_PU), (qaws_scalar)((j + 0.5) / HS_PV), QAWS_SURFACE_EVAL_NORMAL);
		lit = fabs(-0.5 * mid.normal.x + 0.3 * mid.normal.y + 0.81 * mid.normal.z);
		lit = 0.25 + 0.75 * lit;
		fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"rgb(%d,%d,%d)\" stroke=\"none\"/>\n",
			x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3],
			(int)(lit * (j / 3 % 2 ? 42 : 244)), (int)(lit * (j / 3 % 2 ? 157 : 162)), (int)(lit * (j / 3 % 2 ? 143 : 97)));
	}
}

static void figure_pipe(void)
{
	/* a twisted path: a helix with wobbles */
	qaws_vec3 p[12];
	qaws_bspline_desc d;
	qaws_surface_pipe_desc pd;
	qaws_curve* c = NULL;
	qaws_surface* pipe = NULL;
	projection pr = { 612, 420, 52, 1.0 };
	viewport v = { 24, 80, 1176, 620, 0, 1, 0, 1 };
	svg s;
	unsigned int i;
	for (i = 0; i < 12; i++)
	{
		double a = 0.75 * i;
		p[i].x = (qaws_scalar)(2.2 * cos(a) + 0.5 * sin(2.1 * a));
		p[i].y = (qaws_scalar)(2.2 * sin(a));
		p[i].z = (qaws_scalar)(0.55 * i + 0.6 * cos(1.7 * a) - 3);
	}
	memset(&d, 0, sizeof d);
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 12;
	if (qaws_curve_create_bspline(&c, &d) != QAWS_STATUS_OK)
		return;
	memset(&pd, 0, sizeof pd);
	pd.path = c;
	pd.radius_x = (qaws_scalar)0.45;
	if (qaws_surface_create_pipe(&pipe, &pd) != QAWS_STATUS_OK)
	{
		qaws_curve_destroy(c);
		return;
	}
	if (svg_open(&s, "showcase/host7_pipe.svg", 1224, 730, "Pipe with analytic derivatives",
		"a pipe around a twisted cubic path: rotation-minimizing frames by double reflection, lit by the normal of the analytic du and dv"))
	{
		svg_panel(&s, &v, "qaws_surface_create_pipe: stripes along v stay straight (no twist); smooth shading from the analytic normal");
		hs_pipe(&s, &pr, pipe);
		svg_close(&s);
		printf("wrote showcase/host7_pipe.svg\n");
	}
	qaws_surface_destroy(pipe);
	qaws_curve_destroy(c);
}

static void figure_parameterization(void)
{
	/* uneven points: a tight cluster, then long jumps and a sharp turn */
	qaws_vec2 p[8] = { { 0, 0 }, { 0.3f, 0.25f }, { 0.5f, 0.3f }, { 4, 3.2f }, { 6.2f, 3.4f }, { 6.5f, 1.2f }, { 9.5f, 0.6f }, { 9.8f, 3.6f } };
	static char const* names[3] = { "uniform", "chordal (default)", "centripetal" };
	static char const* colors[3] = { "#d1242f", "#0969da", "#2a9d8f" };
	qaws_parameterization par[3] = { QAWS_PARAMETERIZATION_UNIFORM, QAWS_PARAMETERIZATION_CHORDAL, QAWS_PARAMETERIZATION_CENTRIPETAL };
	viewport v = { 24, 80, 1176, 600, -1.5, 11.5, -2.2, 5.6 };
	double line[400 * 2];
	svg s;
	unsigned int k, i;
	if (!svg_open(&s, "showcase/host8_parameterization.svg", 1224, 720, "Parameterizations of fitting",
		"cubic B-splines through the same eight unevenly spaced points (as many control points as points), by qaws_bspline_fit_desc.parameterization"))
		return;
	svg_panel(&s, &v, NULL);
	for (k = 0; k < 3; k++)
	{
		qaws_bspline_fit_desc d;
		qaws_curve* c = NULL;
		memset(&d, 0, sizeof d);
		d.dimension = QAWS_DIMENSION_2D;
		d.data_points = p;
		d.data_point_count = 8;
		d.degree = 3;
		d.control_point_count = 8;
		d.parameterization = par[k];
		if (qaws_curve_fit_bspline(&c, &d) != QAWS_STATUS_OK)
			continue;
		for (i = 0; i < 400; i++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(&e, c, (qaws_scalar)(i / 399.0), QAWS_EVAL_FLAG_POSITION);
			line[2 * i] = vx(&v, e.position.x);
			line[2 * i + 1] = vy(&v, e.position.y);
		}
		svg_polyline(&s, line, 400, colors[k], 2.4, 0.9, 0);
		svg_line(&s, 60, 112 + 22 * k, 92, 112 + 22 * k, colors[k], 3, 1);
		svg_text(&s, 100, 117 + 22 * k, 14, "#1b1f24", "start", names[k]);
		qaws_curve_destroy(c);
	}
	for (i = 0; i < 8; i++)
		svg_circle(&s, vx(&v, p[i].x), vy(&v, p[i].y), 5, "#ffffff", "#1b1f24");
	svg_text(&s, 40, 705, 13, "#57606a", "start", "uniform steps kink inside the cluster, chord lengths swing wide after it; their square roots (centripetal) follow the points");
	svg_close(&s);
	printf("wrote showcase/host8_parameterization.svg\n");
}

int main(void)
{
	figure_parameterization();
	figure_pipe();
	figure_tessellation_surface();
	figure_surface_closest();
	figure_frames();
	figure_arc_length();
	figure_closest();
	figure_tessellation();
	return 0;
}
