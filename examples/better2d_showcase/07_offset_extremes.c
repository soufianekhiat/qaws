/* Figure 7: the hard cases of stroke expansion (as in Levien and Uguray's
   paper), stroked with round ends: a cusp of the curve, collinear control
   points doubling back, a loop, a near-cusp, a tiny curve under a wide
   stroke, a wide inflection. Dashed: the plain parallel curves on both sides
   (swallowtails, jumps); thin red and blue: the raw chains (G1 cubics, the
   evolute through backward runs, round joins at cusps of the curve); filled:
   the stroke after the positive union, its joints dotted. Every one is
   checked against a distance oracle in test 90. */

typedef struct ox_case
{
	char const* name;
	unsigned int n;           /* control points (quadratic or cubic Bezier) */
	double xy[8];
	double delta;
} ox_case;

static void demo_offset_extremes(void)
{
	static ox_case const cases[8] = {
		{ "a cusp (C' = 0), stroked 0.5", 4, { 0, 0, 3, 3, 0, 3, 3, 0 }, 0.5 },
		{ "the same cusp, stroked 1.6", 4, { 0, 0, 3, 3, 0, 3, 3, 0 }, 1.6 },
		{ "collinear points doubling back", 4, { 0, 0, 3, 0, -1, 0, 2, 0 }, 0.5 },
		{ "a loop, stroked 0.8", 4, { 0, 0, 4, 3, -1, 3, 3, 0 }, 0.8 },
		{ "a near-cusp (radius 0.01), stroked 1", 4, { 0, 0, 3, 3, 0.05, 3, 3, 0 }, 1.0 },
		{ "a tight parabola (radius 0.025), stroked 0.6", 3, { -0.2, 0.4, 0, -0.4, 0.2, 0.4 }, 0.6 },
		{ "a tiny S stroked 3 times its size", 4, { 0, 0, 0.1, 0.2, 0.2, -0.2, 0.3, 0 }, 1.0 },
		{ "an inflection, stroked 2", 4, { 0, 0, 1, 2, 2, -2, 3, 0 }, 2.0 },
	};
	int const cols = 4, rows = 2, pw = 320, ph = 300, gap = 16, top = 84;
	int const W = cols * pw + (cols + 1) * gap, H = top + rows * (ph + gap) + 8;
	svg s;
	int k;
	double t0;
	char label[200];
	if (!svg_open(&s, "showcase/b2d7_offset_extremes.svg", W, H, "Offsets of curves: the hard cases",
		"Strokes with round ends. Dashed: plain parallel curves; thin red / blue: raw chains on each side; filled: the stroke after the union, joints dotted. All pass the distance oracle (test 90)"))
		return;
	t0 = (double)clock() / CLOCKS_PER_SEC;
	for (k = 0; k < 8; k++)
	{
		ox_case const* cs = &cases[k];
		qaws_curve* c = co_bezier(cs->xy, cs->n);
		qaws_curve const* cc = c;
		qaws_path_2d p;
		qaws_clip_result* r = NULL;
		qaws_vec2 lo, hi;
		viewport v;
		unsigned int side, cub = 0, chain[2] = { 0, 0 };
		double span, w, h;
		int x0 = gap + (k % cols) * (pw + gap), y0 = top + (k / cols) * (ph + gap);
		p.curves = &cc; p.curve_count = 1; p.closed = 0;
		qaws_path_compute_bounds_2d(&p, &lo, &hi);
		w = (double)(hi.x - lo.x) + 2 * cs->delta;
		h = (double)(hi.y - lo.y) + 2 * cs->delta;
		span = 1.12 * (w > h * pw / (ph - 50) ? w : h * pw / (ph - 50));
		co_view(&v, x0, y0, pw, ph, 0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y) - 0.04 * span, span);
		svg_panel(&s, &v, cs->name);
		if (qaws_offset_paths(&p, 1, (qaws_scalar)cs->delta, QAWS_JOIN_ROUND, QAWS_END_ROUND, &r) == QAWS_STATUS_OK)
			co_result(&s, &v, r, "#c8e1ff", "#0969da", &cub, 0);
		for (side = 0; side < 2; side++)
		{
			qaws_curve* ch[64];
			double d = side ? -cs->delta : cs->delta;
			unsigned int nk = 0, q;
			co_parallel(&s, &v, c, d);
			qaws_internal_offset_curve(c, d, 1e-6 * span, ch, 64, &nk, NULL, 0, NULL);
			for (q = 0; q < nk && q < 64; q++)
			{
				co_curve(&s, &v, ch[q], side ? "#0550ae" : "#cf222e", 1.0, 0.9, 0);
				qaws_curve_destroy(ch[q]);
			}
			chain[side] = nk;
		}
		if (r)
		{
			unsigned int i, j;
			for (i = 0; i < qaws_clip_result_get_path_count(r); i++)
			{
				qaws_path_2d q;
				qaws_clip_result_get_path(r, i, &q);
				for (j = 0; j < q.curve_count; j++)
					co_joints(&s, &v, q.curves[j], "#0969da", 2.2);
			}
		}
		co_curve(&s, &v, c, "#24292f", 1.8, 0.9, 0);
		sprintf(label, "stroke: %u cubics and arcs, %u region%s", cub, r ? qaws_clip_result_get_path_count(r) : 0,
			r && qaws_clip_result_get_path_count(r) == 1 ? "" : "s");
		svg_text(&s, x0 + 10, y0 + ph - 12, 11, "#57606a", "start", label);
		qaws_clip_result_destroy(r);
		qaws_curve_destroy(c);
	}
	svg_close(&s);
	printf("    offset extremes: 8 panels in %.3f s\n", (double)clock() / CLOCKS_PER_SEC - t0);
}
