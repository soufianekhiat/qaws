/* Figure 1: qaws_curve_extract on a curve of every kind. The source in grey,
   the piece between 20% and 75% of its domain as the new curve in colour,
   with that curve's control polygon when it has one, and the kind it came
   out as. */

#define EX_SAMPLES 300

static qaws_curve* ex_make(int which, char const** name)
{
	qaws_curve* c = NULL;
	switch (which)
	{
	case 0:
	{
		static qaws_scalar cp[] = { 0, 0, 1, 3, 3, -1, 4, 2, 5, 0, 6, 1 };
		qaws_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 5; d.control_points = cp; d.control_point_count = 6;
		qaws_curve_create_bezier(&d, &c);
		*name = "Bezier (degree 5)";
		break;
	}
	case 1:
	{
		static qaws_scalar cp[] = { 2, 0, 2, 2, 0, 2 }, w[] = { 1, (qaws_scalar)0.70710678118654752, 1 };
		qaws_rational_bezier_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = cp; d.control_point_count = 3;
		d.weights = w; d.weight_count = 3;
		qaws_curve_create_rational_bezier(&d, &c);
		*name = "rational Bezier (quarter circle)";
		break;
	}
	case 2:
	{
		static qaws_scalar cp[] = { 0, 0, 1, 2, 2, -1, 3, 1, 4, 3, 5, 0, 6, 2 };
		static qaws_scalar kn[] = { 0, 0, 0, 0, (qaws_scalar)0.3, (qaws_scalar)0.45, (qaws_scalar)0.8, 1, 1, 1, 1 };
		qaws_bspline_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = 7;
		d.knots = kn; d.knot_count = 11;
		qaws_curve_create_bspline(&d, &c);
		*name = "cubic B-spline";
		break;
	}
	case 3:
	{
		static double const ux[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, uy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
		static qaws_scalar kn[12] = { 0, 0, 0, (qaws_scalar)0.25, (qaws_scalar)0.25, (qaws_scalar)0.5,
			(qaws_scalar)0.5, (qaws_scalar)0.75, (qaws_scalar)0.75, 1, 1, 1 };
		qaws_scalar cp[18], w[9];
		qaws_nurbs_desc d;
		unsigned int i;
		for (i = 0; i < 9; i++)
		{
			cp[2 * i] = (qaws_scalar)(3 * ux[i]);
			cp[2 * i + 1] = (qaws_scalar)(2 * uy[i]);
			w[i] = (qaws_scalar)(i % 2 ? sqrt(0.5) : 1);
		}
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 2; d.control_points = cp; d.control_point_count = 9;
		d.knots = kn; d.knot_count = 12; d.weights = w; d.weight_count = 9;
		qaws_curve_create_nurbs(&d, &c);
		*name = "NURBS ellipse";
		break;
	}
	case 4:
	{
		static qaws_scalar p[] = { 0, 0, 2, 1, 3, -1, 5, 0 }, v[] = { 1, 2, 2, 0, 1, -2, 1, 1 };
		qaws_hermite_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.points = p; d.derivatives = v;
		d.point_count = 4; d.derivative_count = 4;
		qaws_curve_create_hermite(&d, &c);
		*name = "Hermite";
		break;
	}
	case 5:
	case 6:
	{
		static qaws_scalar p[] = { 0, 0, 1, 2, 3, 2, 4, 0, 2, -1, 1, (qaws_scalar)-0.5 };
		qaws_catmull_rom_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.control_points = p; d.control_point_count = 6;
		d.parameterization = which == 5 ? QAWS_PARAMETERIZATION_CENTRIPETAL : QAWS_PARAMETERIZATION_CHORDAL;
		d.closed = which == 6;
		qaws_curve_create_catmull_rom(&d, &c);
		*name = which == 5 ? "Catmull-Rom (centripetal)" : "Catmull-Rom (chordal, closed)";
		break;
	}
	case 7:
	{
		static qaws_scalar p[] = { 0, 0, 2, 3, 5, 1, 6, 4 }, t[] = { 0, 1, (qaws_scalar)2.5, 3 };
		qaws_trajectory_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.key_positions = p; d.key_count = 4;
		d.key_times = t; d.key_time_count = 4;
		qaws_curve_create_trajectory(&d, &c);
		*name = "trajectory";
		break;
	}
	case 8:
	{
		static qaws_scalar p[] = { 0, 0, 1, 2, 3, 2, 4, 0, 6, 1 };
		qaws_yuksel_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.control_points = p; d.control_point_count = 5;
		d.mode = QAWS_YUKSEL_MODE_CIRCULAR;
		qaws_curve_create_yuksel(&d, &c);
		*name = "Yuksel (circular, fitted)";
		break;
	}
	case 9:
	{
		static qaws_scalar co[] = { 1, 0, 2, 1, -1, 3, (qaws_scalar)0.5, -2, (qaws_scalar)0.25, (qaws_scalar)0.5 };
		qaws_polynomial_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.degree = 4; d.coefficients = co; d.coefficient_count = 5;
		d.t_min = -1; d.t_max = 2;
		qaws_curve_create_polynomial(&d, &c);
		*name = "polynomial (degree 4)";
		break;
	}
	case 10:
	{
		qaws_arc_segment s[2];
		qaws_arc_desc d;
		memset(s, 0, sizeof(s));
		s[0].radius = 2; s[0].angle_start = 0; s[0].angle_end = (qaws_scalar)1.5;
		s[1].center[0] = (qaws_scalar)cos(1.5); s[1].center[1] = (qaws_scalar)sin(1.5);
		s[1].radius = 1; s[1].angle_start = (qaws_scalar)1.5; s[1].angle_end = (qaws_scalar)4.5;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = s; d.segment_count = 2;
		qaws_curve_create_arc(&d, &c);
		*name = "arc (two segments)";
		break;
	}
	case 11:
	{
		qaws_clothoid_desc d;
		memset(&d, 0, sizeof(d));
		d.start_angle = (qaws_scalar)0.3; d.start_curvature = (qaws_scalar)-0.2;
		d.end_curvature = (qaws_scalar)1.1; d.length = 6;
		qaws_curve_create_clothoid(&d, &c);
		*name = "clothoid";
		break;
	}
	case 12:
	{
		static qaws_scalar p[] = { 0, 0, 3, 0, 4, 2, 2, 4, -1, 2 };
		qaws_subdivision_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.scheme = QAWS_SUBDIVISION_LANE_RIESENFELD_3;
		d.control_points = p; d.control_point_count = 5; d.closed = 1; d.refinement_levels = 3;
		qaws_curve_create_subdivision(&d, &c);
		*name = "subdivision (closed)";
		break;
	}
	default:
	{
		static qaws_scalar l[] = { 0, 0, 2, 0 };
		static qaws_scalar b[] = { 2, 0, 3, 0, 4, 1, 4, 2 };
		qaws_arc_segment s;
		qaws_bezier_desc bd;
		qaws_arc_desc ad;
		qaws_composite_desc d;
		qaws_curve* seg[3] = { NULL, NULL, NULL };
		memset(&bd, 0, sizeof(bd));
		bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = l; bd.control_point_count = 2;
		qaws_curve_create_bezier(&bd, &seg[0]);
		bd.degree = 3; bd.control_points = b; bd.control_point_count = 4;
		qaws_curve_create_bezier(&bd, &seg[1]);
		memset(&s, 0, sizeof(s));
		s.center[0] = 3; s.center[1] = 2; s.radius = 1; s.angle_start = 0; s.angle_end = (qaws_scalar)PI;
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = &s; ad.segment_count = 1;
		qaws_curve_create_arc(&ad, &seg[2]);
		memset(&d, 0, sizeof(d));
		d.dimension = QAWS_DIMENSION_2D; d.segments = seg; d.segment_count = 3;
		qaws_curve_create_composite(&d, &c);
		*name = "composite (line, Bezier, arc)";
		break;
	}
	}
	return c;
}

static char const* ex_kind_name(qaws_curve_kind k)
{
	switch (k)
	{
	case QAWS_CURVE_KIND_BEZIER: return "Bezier";
	case QAWS_CURVE_KIND_RATIONAL_BEZIER: return "rational Bezier";
	case QAWS_CURVE_KIND_BSPLINE: return "B-spline";
	case QAWS_CURVE_KIND_NURBS: return "NURBS";
	case QAWS_CURVE_KIND_POLYNOMIAL: return "polynomial";
	case QAWS_CURVE_KIND_ARC: return "arc";
	case QAWS_CURVE_KIND_CLOTHOID: return "clothoid";
	case QAWS_CURVE_KIND_COMPOSITE: return "composite";
	default: return "?";
	}
}

static void ex_sample(qaws_curve const* c, double t0, double t1, double* xy, int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		qaws_eval_result_2d r;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(t0 + (t1 - t0) * i / (n - 1)), QAWS_EVAL_FLAG_POSITION, &r);
		xy[2 * i] = r.position.x;
		xy[2 * i + 1] = r.position.y;
	}
}

static void demo_extract(void)
{
	int const cols = 5, rows = 3, pw = 300, ph = 230, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	static double src[2 * EX_SAMPLES], pc[2 * EX_SAMPLES], scr[2 * EX_SAMPLES];
	double worst = 0.0;
	svg s;
	int k;

	if (!svg_open(&s, "showcase/b2d1_extract.svg", W, H, "Exact pieces of curves, as new curves",
		"qaws_curve_extract between 20% and 75% of the domain: source in grey, the new curve in colour with its control polygon; the piece lies on the source to rounding"))
		return;
	for (k = 0; k < 14; k++)
	{
		char const* name = "";
		qaws_curve* c = ex_make(k, &name);
		qaws_curve* e = NULL;
		qaws_range r, re;
		viewport v;
		double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, t0, t1, cx, cy, half, err = 0.0;
		char label[160];
		int i, col = k % cols, row = k / cols;
		if (!c) continue;
		r = qaws_curve_get_parameter_range(c);
		t0 = r.min_value + 0.2 * (r.max_value - r.min_value);
		t1 = r.min_value + 0.75 * (r.max_value - r.min_value);
		qaws_curve_extract(c, (qaws_scalar)t0, (qaws_scalar)t1, &e);
		ex_sample(c, r.min_value, r.max_value, src, EX_SAMPLES);
		for (i = 0; i < EX_SAMPLES; i++)
		{
			if (src[2 * i] < lo[0]) lo[0] = src[2 * i];
			if (src[2 * i] > hi[0]) hi[0] = src[2 * i];
			if (src[2 * i + 1] < lo[1]) lo[1] = src[2 * i + 1];
			if (src[2 * i + 1] > hi[1]) hi[1] = src[2 * i + 1];
		}
		cx = 0.5 * (lo[0] + hi[0]); cy = 0.5 * (lo[1] + hi[1]);
		half = 0.5 * (hi[0] - lo[0] > (hi[1] - lo[1]) * pw / ph ? hi[0] - lo[0] : (hi[1] - lo[1]) * pw / ph) * 1.25;
		v.x0 = gap + col * (pw + gap); v.y0 = top + row * (ph + gap); v.w = pw; v.h = ph;
		v.xmin = cx - half; v.xmax = cx + half;
		v.ymin = cy - half * ph / pw - 0.06 * half; v.ymax = cy + half * ph / pw - 0.06 * half;
		v.ymin -= 0.1 * half; v.ymax += 0.1 * half;
		svg_panel(&s, &v, name);
		for (i = 0; i < EX_SAMPLES; i++)
		{
			scr[2 * i] = vx(&v, src[2 * i]);
			scr[2 * i + 1] = vy(&v, src[2 * i + 1]);
		}
		svg_polyline(&s, scr, EX_SAMPLES, "#8c959f", 2.0, 0.7, 0);
		if (e)
		{
			qaws_curve_kind kind = qaws_curve_get_kind(e);
			re = qaws_curve_get_parameter_range(e);
			ex_sample(e, re.min_value, re.max_value, pc, EX_SAMPLES);
			/* deviation from the source at the same parameters, where they map linearly */
			if (kind == QAWS_CURVE_KIND_BSPLINE || kind == QAWS_CURVE_KIND_NURBS || kind == QAWS_CURVE_KIND_BEZIER ||
				kind == QAWS_CURVE_KIND_POLYNOMIAL || kind == QAWS_CURVE_KIND_RATIONAL_BEZIER)
			{
				for (i = 0; i < EX_SAMPLES; i += 7)
				{
					qaws_eval_result_2d a;
					double t = t0 + (t1 - t0) * i / (EX_SAMPLES - 1), d;
					qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION, &a);
					d = hypot(a.position.x - pc[2 * i], a.position.y - pc[2 * i + 1]);
					if (k != 8 && d > err) err = d;
				}
			}
			if (err > worst) worst = err;
			/* control polygon of the new curve */
			if (kind == QAWS_CURVE_KIND_BEZIER || kind == QAWS_CURVE_KIND_BSPLINE ||
				kind == QAWS_CURVE_KIND_NURBS || kind == QAWS_CURVE_KIND_RATIONAL_BEZIER)
			{
				qaws_scalar cps[3 * 256];
				unsigned int n = 0, j;
				if (qaws_curve_get_control_points(e, cps, 256, &n) == QAWS_STATUS_OK && n > 1)
				{
					double poly[2 * 256];
					for (j = 0; j < n && j < 256; j++)
					{
						poly[2 * j] = vx(&v, cps[2 * j]);
						poly[2 * j + 1] = vy(&v, cps[2 * j + 1]);
					}
					svg_polyline(&s, poly, (int)n, "#cf222e", 1.0, 0.55, 1);
					for (j = 0; j < n && j < 256; j++)
						svg_circle(&s, poly[2 * j], poly[2 * j + 1], 2.2, "#ffffff", "#cf222e");
				}
			}
			for (i = 0; i < EX_SAMPLES; i++)
			{
				scr[2 * i] = vx(&v, pc[2 * i]);
				scr[2 * i + 1] = vy(&v, pc[2 * i + 1]);
			}
			svg_polyline(&s, scr, EX_SAMPLES, "#0969da", 3.2, 0.95, 0);
			svg_circle(&s, scr[0], scr[1], 4, "#1a7f37", "#ffffff");
			svg_circle(&s, scr[2 * EX_SAMPLES - 2], scr[2 * EX_SAMPLES - 1], 4, "#bf8700", "#ffffff");
			sprintf(label, "piece: %s", ex_kind_name(kind));
			svg_text(&s, v.x0 + 10, v.y0 + ph - 12, 12, "#0969da", "start", label);
			qaws_curve_destroy(e);
		}
		qaws_curve_destroy(c);
	}
	/* legend in the last panel slot */
	{
		double x = gap + 4 * (pw + gap) + 16, y = top + 2 * (ph + gap) + 40;
		char line[160];
		svg_line(&s, x, y, x + 40, y, "#8c959f", 2.0, 0.7);
		svg_text(&s, x + 50, y + 4, 13, "#24292f", "start", "source curve");
		svg_line(&s, x, y + 28, x + 40, y + 28, "#0969da", 3.2, 0.95);
		svg_text(&s, x + 50, y + 32, 13, "#24292f", "start", "extracted piece (new curve)");
		svg_line(&s, x, y + 56, x + 40, y + 56, "#cf222e", 1.0, 0.55);
		svg_text(&s, x + 50, y + 60, 13, "#24292f", "start", "its control polygon");
		svg_circle(&s, x + 6, y + 84, 4, "#1a7f37", "#ffffff");
		svg_circle(&s, x + 30, y + 84, 4, "#bf8700", "#ffffff");
		svg_text(&s, x + 50, y + 88, 13, "#24292f", "start", "C(t0), C(t1)");
		sprintf(line, "worst deviation: %.1e", worst);
		svg_text(&s, x, y + 120, 13, "#57606a", "start", line);
	}
	svg_close(&s);
	printf("    extract: 14 kinds, worst deviation of a piece from its source %.1e\n", worst);
}
