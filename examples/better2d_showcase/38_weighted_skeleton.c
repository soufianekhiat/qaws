/* Figure 38: weighted skeletons (qaws_skeleton's weight_fn on lines and
   circular arcs). Each shape unweighted (top) and weighted (bottom): the
   edges drawn orange move at their weight (labelled), arcs and the circles
   opened at reflex corners too (at the faster neighbour's speed), so the
   medial axis (sigma = pi) is the weighted one, each edge's distance over
   its weight. */

typedef struct wsk_user
{
	csk_shape const* sh;
	double const* w;
} wsk_user;

static qaws_scalar wsk_weight(void* user, unsigned int path, unsigned int curve, unsigned int span)
{
	wsk_user const* u = (wsk_user const*)user;
	unsigned int l, off = 0;
	(void)span;
	for (l = 0; l < path; l++) off += u->sh->count[l];
	return (qaws_scalar)u->w[off + curve];
}

static void demo_weighted_skeleton(void)
{
	static csk_shape sh[4];
	static double const wr[4] = { 2, 1, 1, 1 }, wl[6] = { 2, 1, 1, 1.5, 1, 1 }, ws[4] = { 1, 2, 1, 0.5 }, wb[6] = { 1, 1, 1, 0.5, 1, 1 };
	static double const* const wts[4] = { wr, wl, ws, wb };
	static double const sig[4] = { 0, PI, 0, PI };
	unsigned int n, i;
	double const pi = PI;
	svg s;
	double pw = 330, ph = 250;
	clock_t t0 = clock();
	memset(sh, 0, sizeof(sh));
	n = 0; sh[0].name = "rectangle, mitered";
	csk_line(&sh[0], &n, 0, 0, 4, 0); csk_line(&sh[0], &n, 4, 0, 4, 2); csk_line(&sh[0], &n, 4, 2, 0, 2); csk_line(&sh[0], &n, 0, 2, 0, 0);
	sh[0].count[0] = n; sh[0].nloops = 1; sh[0].box[0] = -0.1; sh[0].box[1] = 4.1; sh[0].box[2] = -0.1; sh[0].box[3] = 2.1;
	n = 0; sh[1].name = "L, medial axis";
	csk_line(&sh[1], &n, 0, 0, 4, 0); csk_line(&sh[1], &n, 4, 0, 4, 1); csk_line(&sh[1], &n, 4, 1, 1, 1);
	csk_line(&sh[1], &n, 1, 1, 1, 3); csk_line(&sh[1], &n, 1, 3, 0, 3); csk_line(&sh[1], &n, 0, 3, 0, 0);
	sh[1].count[0] = n; sh[1].nloops = 1; sh[1].box[0] = -0.1; sh[1].box[1] = 4.1; sh[1].box[2] = -0.1; sh[1].box[3] = 3.1;
	n = 0; sh[2].name = "stadium, mitered";
	csk_line(&sh[2], &n, -1, -1, 1, -1); csk_arc(&sh[2], &n, 1, 0, 1, -pi / 2, pi / 2);
	csk_line(&sh[2], &n, 1, 1, -1, 1); csk_arc(&sh[2], &n, -1, 0, 1, pi / 2, 3 * pi / 2);
	sh[2].count[0] = n; sh[2].nloops = 1; sh[2].box[0] = -2.05; sh[2].box[1] = 2.05; sh[2].box[2] = -1.05; sh[2].box[3] = 1.05;
	n = 0; sh[3].name = "square with a round bump, medial axis";
	csk_line(&sh[3], &n, 0, 0, 2, 0); csk_line(&sh[3], &n, 2, 0, 2, 2); csk_line(&sh[3], &n, 2, 2, 1.6, 2);
	csk_arc(&sh[3], &n, 1, 2, 0.6, 0, pi); csk_line(&sh[3], &n, 0.4, 2, 0, 2); csk_line(&sh[3], &n, 0, 2, 0, 0);
	sh[3].count[0] = n; sh[3].nloops = 1; sh[3].box[0] = -0.1; sh[3].box[1] = 2.1; sh[3].box[2] = -0.1; sh[3].box[3] = 2.7;
	if (!svg_open(&s, "showcase/b2d38_weighted_skeleton.svg", (int)(40 + 4 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Weighted skeletons: every edge at its own speed, lines and arcs",
		"qaws_skeleton with weight_fn: top unweighted, bottom the orange edges at their weights; the circles opened at reflex corners move at the faster neighbour's speed"))
		return;
	for (i = 0; i < 4; i++)
	{
		unsigned int c;
		for (c = 0; c < 2; c++)
		{
			viewport v;
			wsk_user u;
			double bx = sh[i].box[1] - sh[i].box[0], by = sh[i].box[3] - sh[i].box[2], sc;
			unsigned int k;
			char label[160];
			sc = (pw - 20) / bx < (ph - 30) / by ? (pw - 20) / bx : (ph - 30) / by;
			v.w = bx * sc; v.h = by * sc;
			v.x0 = 24 + i * (pw + 12) + (pw - v.w) / 2;
			v.y0 = 104 + c * (ph + 44) + 10 + (ph - 30 - v.h) / 2;
			v.xmin = sh[i].box[0]; v.xmax = sh[i].box[1]; v.ymin = sh[i].box[2]; v.ymax = sh[i].box[3];
			sprintf(label, "%s%s", sh[i].name, c ? ", weighted" : "");
			svg_text(&s, 24 + i * (pw + 12) + pw / 2, 104 + c * (ph + 44), 12, "#24292f", "middle", label);
			u.sh = &sh[i];
			u.w = wts[i];
			csk_draw(&s, &v, &sh[i], sig[i], c ? wsk_weight : NULL, &u);
			/* the weighted edges over the boundary, with their weights */
			for (k = 0; c && k < sh[i].count[0]; k++)
				if (wts[i][k] != 1)
				{
					csk_seg const* g = &sh[i].seg[k];
					qaws_curve* cv = csk_curve(g);
					qaws_range rg = qaws_curve_get_parameter_range(cv);
					qaws_eval_result_2d e;
					unsigned int j;
					char wtxt[32];
					fprintf(s.f, "<polyline fill=\"none\" stroke=\"#bc4c00\" stroke-width=\"3\" points=\"");
					for (j = 0; j <= 32; j++)
					{
						qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(rg.min_value + (rg.max_value - rg.min_value) * j / 32), QAWS_EVAL_FLAG_POSITION);
						fprintf(s.f, "%.2f,%.2f ", vx(&v, e.position.x), vy(&v, e.position.y));
					}
					fprintf(s.f, "\"/>\n");
					qaws_curve_evaluate_2d(&e, cv, (qaws_scalar)(0.5 * (rg.min_value + rg.max_value)), QAWS_EVAL_FLAG_POSITION);
					sprintf(wtxt, "w %g", wts[i][k]);
					svg_text(&s, vx(&v, e.position.x), vy(&v, e.position.y) - 6, 11, "#bc4c00", "middle", wtxt);
					qaws_curve_destroy(cv);
				}
		}
	}
	svg_close(&s);
	printf("figure 38: weighted skeletons, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
