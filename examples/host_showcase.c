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
		qaws_curve_evaluate_3d(c, t[i], QAWS_EVAL_FLAG_POSITION, &e);
		pos[i] = e.position;
		if (frenet)
		{
			qaws_vec3 tg, bn;
			ok[i] = qaws_curve_compute_frenet_frame_3d(c, t[i], &tg, &nrm[i], &bn) == QAWS_STATUS_OK
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
	if (qaws_curve_create_bspline(&d, &c) != QAWS_STATUS_OK)
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
		qaws_curve_evaluate_2d(c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / (HS_SAMPLES - 1), QAWS_EVAL_FLAG_POSITION, &e);
		line[2 * i] = vx(v, e.position.x);
		line[2 * i + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, line, HS_SAMPLES, "#1b1f24", 1.6, 1, 0);
	qaws_curve_compute_arc_length(c, r.min_value, r.max_value, &total);
	for (i = 0; i < HS_DOTS; i++)
	{
		qaws_eval_result_2d e;
		t = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / (HS_DOTS - 1);
		if (by_length)
			qaws_curve_compute_parameter_at_arc_length(c, r.min_value, total * (qaws_scalar)i / (HS_DOTS - 1), &t);
		qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION, &e);
		xy[2 * i] = e.position.x;
		xy[2 * i + 1] = e.position.y;
		if (i > 0)
		{
			qaws_scalar arc = 0;
			double d;
			qaws_curve_compute_arc_length(c, t_prev, t, &arc);
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
	if (qaws_curve_create_bspline(&d, &c) != QAWS_STATUS_OK)
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
		qaws_curve_evaluate_2d(c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 1199, QAWS_EVAL_FLAG_POSITION, &e);
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
	if (qaws_curve_create_bspline(&d, &c) != QAWS_STATUS_OK)
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
		qaws_curve_evaluate_2d(c, r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 31, QAWS_EVAL_FLAG_POSITION, &e);
		dd = hypot(e.position.x, e.position.y);
		svg_circle(&s, vx(&vl, e.position.x), vy(&vl, e.position.y), 3.5, "#f4a261", "#1b1f24");
		if (dd < best)
		{
			best = dd;
			ts = r.min_value + (r.max_value - r.min_value) * (qaws_scalar)i / 31;
		}
	}
	qaws_curve_evaluate_2d(c, ts, QAWS_EVAL_FLAG_POSITION, &e);
	svg_line(&s, vx(&vl, 0), vy(&vl, 0), vx(&vl, e.position.x), vy(&vl, e.position.y), "#f4a261", 2.5, 1);
	snprintf(note, sizeof note, "distance %.6f", best);
	svg_text(&s, vl.x0 + 16, vl.y0 + vl.h - 16, 13, "#57606a", "start", note);
	qaws_curve_find_closest_point_certified_2d(c, q, &t, &dist, &gap);
	qaws_curve_evaluate_2d(c, t, QAWS_EVAL_FLAG_POSITION, &e);
	svg_line(&s, vx(&vr, 0), vy(&vr, 0), vx(&vr, e.position.x), vy(&vr, e.position.y), "#2a9d8f", 2.5, 1);
	svg_circle(&s, vx(&vr, e.position.x), vy(&vr, e.position.y), 5, "#2a9d8f", "#1b1f24");
	snprintf(note, sizeof note, "distance %.12f, gap %.1e (the minimum is in [d - gap, d])", (double)dist, (double)gap);
	svg_text(&s, vr.x0 + 16, vr.y0 + vr.h - 16, 13, "#57606a", "start", note);
	svg_text(&s, 640, 692, 13, "#57606a", "start", "qaws_curve_find_closest_point_certified_2d: convex-hull and fat-line bounds, de Casteljau halving");
	svg_close(&s);
	qaws_curve_destroy(c);
	printf("wrote showcase/host3_closest.svg\n");
}

int main(void)
{
	figure_frames();
	figure_arc_length();
	figure_closest();
	return 0;
}
