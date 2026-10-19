/*
 * host_showcase.c - features for host applications
 *
 * Writes SVG figures into ./showcase/:
 *   host1_frames.svg   a ribbon swept along an S-shaped space curve with a
 *                      straight middle: by the Frenet frame (undefined on the
 *                      straight piece, flipping at the inflection) and by
 *                      rotation-minimizing frames (qaws_curve_compute_
 *                      rotation_minimizing_frames_3d: continuous, twist-free)
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

int main(void)
{
	figure_frames();
	return 0;
}
