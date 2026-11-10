/* Figure 2: qaws_curve_create_interpolating_spline with each method through the
   same points, open and closed, with a curvature comb (a hair at each sample,
   along the normal, as long as the curvature, capped), and qaws_curve_offset_cubic_2d
   of each closed spline at four distances. */

#define SP_SAMPLES 400
#define SP_COMB 160

static qaws_vec2 const sp_points[7] = { { 0, 0 }, { 2, 1 }, { 4, (qaws_scalar)0.5 }, { 5, 3 }, { (qaws_scalar)3.5, (qaws_scalar)4.5 },
	{ (qaws_scalar)1.5, (qaws_scalar)3.2 }, { (qaws_scalar)0.3, (qaws_scalar)2.4 } };

static qaws_curve* sp_make(qaws_interpolating_method m, int closed)
{
	qaws_interpolating_spline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof d);
	d.points = sp_points;
	d.point_count = 7;
	d.method = m;
	d.closed = closed;
	return qaws_curve_create_interpolating_spline(&c, &d) == QAWS_STATUS_OK ? c : NULL;
}

static void sp_draw(svg* s, viewport const* v, qaws_curve const* c, char const* color, double width)
{
	static double xy[2 * SP_SAMPLES];
	qaws_range r = qaws_curve_get_parameter_range(c);
	int i;
	for (i = 0; i < SP_SAMPLES; i++)
	{
		qaws_eval_result_2d e;
		qaws_curve_evaluate_2d(&e, c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * i / (SP_SAMPLES - 1)), QAWS_EVAL_FLAG_POSITION);
		xy[2 * i] = vx(v, e.position.x);
		xy[2 * i + 1] = vy(v, e.position.y);
	}
	svg_polyline(s, xy, SP_SAMPLES, color, width, 0.95, 0);
}

static void sp_comb(svg* s, viewport const* v, qaws_curve const* c, double scale)
{
	static double tip[2 * SP_COMB];
	qaws_range r = qaws_curve_get_parameter_range(c);
	int i;
	for (i = 0; i < SP_COMB; i++)
	{
		qaws_eval_result_2d e;
		double dx, dy, len, k, x, y;
		qaws_curve_evaluate_2d(&e, c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * i / (SP_COMB - 1)),
			QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2);
		dx = e.d1.x; dy = e.d1.y;
		len = sqrt(dx * dx + dy * dy);
		k = len > 0 ? scale * (dx * e.d2.y - dy * e.d2.x) / (len * len * len) : 0;
		/* capped so a tight turn stays in the panel */
		k = k > 1.2 ? 1.2 : k < -1.2 ? -1.2 : k;
		/* away from the centre of curvature */
		x = e.position.x + k * dy / (len > 0 ? len : 1);
		y = e.position.y - k * dx / (len > 0 ? len : 1);
		tip[2 * i] = vx(v, x);
		tip[2 * i + 1] = vy(v, y);
		svg_line(s, vx(v, e.position.x), vy(v, e.position.y), tip[2 * i], tip[2 * i + 1], "#bf8700", 0.8, 0.6);
	}
	svg_polyline(s, tip, SP_COMB, "#bf8700", 1.0, 0.8, 0);
}

static void demo_splines(void)
{
	static qaws_interpolating_method const methods[5] = { QAWS_INTERPOLATING_CUBIC_C2, QAWS_INTERPOLATING_HOBBY,
		QAWS_INTERPOLATING_KAPPA, QAWS_INTERPOLATING_CUBIC_G2, QAWS_INTERPOLATING_S9 };
	static char const* const names[5] = { "uniform C2 cubic", "Hobby", "kappa curves", "cubic G2", "S9" };
	static double const dist[4] = { 0.6, 0.3, -0.3, -0.6 };
	static char const* const dcolor[4] = { "#8250df", "#0969da", "#1a7f37", "#cf222e" };
	int const cols = 5, rows = 3, pw = 300, ph = 260, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 40;
	svg s;
	int m, row;

	if (!svg_open(&s, "shots/b2d26_splines.svg", W, H, "Splines through points, and their offsets as cubics",
		"qaws_curve_create_interpolating_spline: open, closed with its curvature comb; qaws_curve_offset_cubic_2d of the closed one at +-0.3, +-0.6 (tolerance 1e-3)"))
		return;
	for (row = 0; row < rows; row++)
	{
		for (m = 0; m < 5; m++)
		{
			qaws_curve* c = sp_make(methods[m], row > 0);
			viewport v;
			char label[160];
			int i;
			v.x0 = gap + m * (pw + gap); v.y0 = top + row * (ph + gap); v.w = pw; v.h = ph;
			v.xmin = -1.3; v.xmax = 6.3;
			v.ymin = 2.25 - 0.5 * 7.6 * ph / pw; v.ymax = 2.25 + 0.5 * 7.6 * ph / pw;
			sprintf(label, "%s, %s", names[m], row == 0 ? "open" : row == 1 ? "closed, curvature comb" : "closed, offsets");
			svg_panel(&s, &v, label);
			if (!c)
			{
				svg_text(&s, v.x0 + 10, v.y0 + ph / 2, 12, "#cf222e", "start", "not made");
				continue;
			}
			if (row == 2)
			{
				unsigned int spans = 0;
				int d;
				for (d = 0; d < 4; d++)
				{
					qaws_curve* o = NULL;
					if (qaws_curve_offset_cubic_2d(&o, c, (qaws_scalar)dist[d], (qaws_scalar)1e-3) == QAWS_STATUS_OK && o)
					{
						spans += qaws_curve_get_span_count(o);
						sp_draw(&s, &v, o, dcolor[d], 1.4);
						qaws_curve_destroy(o);
					}
				}
				sp_draw(&s, &v, c, "#57606a", 2.0);
				sprintf(label, "%u offset cubics in all", spans);
				svg_text(&s, v.x0 + 10, v.y0 + ph - 12, 12, "#57606a", "start", label);
			}
			else
			{
				if (row == 1)
					sp_comb(&s, &v, c, 0.45);
				sp_draw(&s, &v, c, "#0969da", 2.4);
				sprintf(label, "%u cubic or quadratic spans", qaws_curve_get_span_count(c));
				svg_text(&s, v.x0 + 10, v.y0 + ph - 12, 12, "#57606a", "start", label);
			}
			for (i = 0; i < 7; i++)
				svg_circle(&s, vx(&v, sp_points[i].x), vy(&v, sp_points[i].y), 3.5, "#ffffff", "#24292f");
			qaws_curve_destroy(c);
		}
	}
	{
		double x = gap + 10, y = H - 22;
		int d;
		char line[64];
		for (d = 0; d < 4; d++)
		{
			svg_line(&s, x, y, x + 30, y, dcolor[d], 1.4, 0.95);
			sprintf(line, "offset %+.1f (left positive)", dist[d]);
			svg_text(&s, x + 36, y + 4, 12, "#24292f", "start", line);
			x += 230;
		}
		svg_line(&s, x, y, x + 30, y, "#bf8700", 1.0, 0.8);
		svg_text(&s, x + 36, y + 4, 12, "#24292f", "start", "curvature comb");
	}
	svg_close(&s);
	printf("    splines: 5 methods open and closed, offsets at 4 distances\n");
}
