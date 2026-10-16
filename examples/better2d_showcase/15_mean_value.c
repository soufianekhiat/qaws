/* Figure 15: mean value coordinates over regions bounded by curves
   (qaws_mean_value). Top: a region of cubic Beziers with a round hole, its
   handles (the curves' control points) at rest and moved: a grid inside is
   carried along. Bottom: the transfinite interpolant of data given on the
   boundary, over a concave polygon with a hole and over a NURBS ellipse. */

static qaws_curve* mvf_cubic(double const* xy)
{
	qaws_scalar cp[8];
	qaws_bezier_desc d;
	qaws_curve* c = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++) cp[i] = (qaws_scalar)xy[i];
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D; d.degree = 3; d.control_points = cp; d.control_point_count = 4;
	qaws_curve_create_bezier(&c, &d);
	return c;
}

static void mvf_points(svg* s, viewport const* v, qaws_vec2 const* p, unsigned int n, char const* color, double width)
{
	double* xy = (double*)malloc(sizeof(double) * 2 * (n + 1));
	unsigned int k;
	for (k = 0; k < n; k++) { xy[2 * k] = vx(v, p[k].x); xy[2 * k + 1] = vy(v, p[k].y); }
	svg_polyline(s, xy, (int)n, color, width, 1.0, 0);
	free(xy);
}

/* the region's grid lines (cut by the region), each mapped by the handles hh */
static void mvf_grid(svg* s, viewport const* v, qaws_mean_value const* m, qaws_region_sampler const* rs, qaws_vec2 const* hh,
	double const* box, unsigned int lines)
{
	unsigned int line, k;
	for (line = 0; line < 2 * lines; line++)
	{
		int horiz = line < lines;
		double c = horiz ? box[2] + (box[3] - box[2]) * ((line % lines) + 0.5) / lines : box[0] + (box[1] - box[0]) * ((line % lines) + 0.5) / lines;
		qaws_vec2 run[401], out[401];
		unsigned int nr = 0;
		for (k = 0; k <= 401; k++)
		{
			qaws_vec2 q;
			int in = 0;
			if (horiz) { q.x = (qaws_scalar)(box[0] + (box[1] - box[0]) * k / 400.0); q.y = (qaws_scalar)c; }
			else { q.x = (qaws_scalar)c; q.y = (qaws_scalar)(box[2] + (box[3] - box[2]) * k / 400.0); }
			if (k <= 400) qaws_region_sampler_contains(rs, q, &in);
			if (in) run[nr++] = q;
			else
			{
				if (nr > 1 && qaws_mean_value_deform(m, hh, run, nr, out) == QAWS_STATUS_OK)
					mvf_points(s, v, out, nr, horiz ? "#2f6fbf" : "#cf6a2a", 0.9);
				nr = 0;
			}
		}
	}
}

/* the boundary curves mapped by the handles hh (on the boundary the
   coordinates are the curves' basis: the moved curves themselves) */
static void mvf_boundary(svg* s, viewport const* v, qaws_mean_value const* m, qaws_curve* const* cs, unsigned int nc, qaws_vec2 const* hh)
{
	unsigned int i, k;
	for (i = 0; i < nc; i++)
	{
		qaws_range rg = qaws_curve_get_parameter_range(cs[i]);
		qaws_vec2 bp[161], bo[161];
		for (k = 0; k <= 160; k++)
		{
			qaws_eval_result_2d e;
			qaws_curve_evaluate_2d(&e, cs[i], rg.min_value + (rg.max_value - rg.min_value) * k / 160, QAWS_EVAL_FLAG_POSITION);
			bp[k] = e.position;
		}
		if (qaws_mean_value_deform(m, hh, bp, 161, bo) == QAWS_STATUS_OK)
			mvf_points(s, v, bo, 161, "#1b1f24", 2.2);
	}
}

/* the interpolant of f over the region as cells, by the heat ramp */
typedef struct mvf_data
{
	int kind;
} mvf_data;

static void mvf_fn(void* user, unsigned int path, unsigned int curve, qaws_scalar t, qaws_vec2 p, qaws_scalar* out)
{
	mvf_data const* d = (mvf_data const*)user;
	double x = p.x, y = p.y;
	(void)curve; (void)t;
	if (d->kind == 0)
		out[0] = (qaws_scalar)(path ? 1.0 : 0.5 + 0.5 * sin(3 * x) * cos(2 * y));     /* the hole hot */
	else
		out[0] = (qaws_scalar)(0.5 + 0.5 * sin(5 * atan2(y + 0.2, x - 0.1)));          /* five waves round the ellipse */
}

static void mvf_field(svg* s, viewport const* v, qaws_mean_value const* m, qaws_region_sampler const* rs, double const* box, int kind)
{
	unsigned int i, j, nx = 120, ny = 80;
	mvf_data d;
	d.kind = kind;
	for (j = 0; j < ny; j++)
		for (i = 0; i < nx; i++)
		{
			qaws_vec2 q;
			int in = 0;
			qaws_scalar val = 0;
			char col[32];
			double x0 = box[0] + (box[1] - box[0]) * i / nx, x1 = box[0] + (box[1] - box[0]) * (i + 1) / nx;
			double y0 = box[2] + (box[3] - box[2]) * j / ny, y1 = box[2] + (box[3] - box[2]) * (j + 1) / ny;
			q.x = (qaws_scalar)(0.5 * (x0 + x1)); q.y = (qaws_scalar)(0.5 * (y0 + y1));
			qaws_region_sampler_contains(rs, q, &in);
			if (!in || qaws_mean_value_interpolate(m, q, mvf_fn, &d, 1, &val) != QAWS_STATUS_OK)
				continue;
			heat((double)val, col);
			fprintf(s->f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\"/>\n",
				vx(v, x0), vy(v, y1), vx(v, x1) - vx(v, x0) + 0.4, vy(v, y0) - vy(v, y1) + 0.4, col);
		}
}

static void demo_mean_value(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 560, ph = 360;
	static double const blob[4][8] = {
		{ 0, 0, 1, -0.9, 2.4, -0.5, 2.2, 0.6 }, { 2.2, 0.6, 2, 1.6, 1.4, 1.1, 1, 1 },
		{ 1, 1, 0.6, 0.9, 0.2, 1.8, -0.5, 1.2 }, { -0.5, 1.2, -1.2, 0.6, -0.6, 0.2, 0, 0 } };
	qaws_curve* cs[5];
	qaws_curve const* cv[5];
	qaws_path_2d p[2];
	qaws_mean_value_desc d;
	qaws_mean_value* m = NULL;
	qaws_region_sampler_desc rd;
	qaws_region_sampler* rs = NULL;
	unsigned int i, nh = 0, panel;
	qaws_vec2 *h = NULL, *mv = NULL;
	double gbox[4] = { -1.0, 2.4, -0.7, 1.6 };
	if (!svg_open(&s, "showcase/b2d15_mean_value.svg", (int)(48 + 2 * (pw + 12)), (int)(110 + 2 * (ph + 44)),
		"Mean value coordinates over regions bounded by curves",
		"qaws_mean_value: the handles are the boundary curves' control points; coordinates by boundary integrals of their basis functions"))
		return;
	/* the blob with a round hole */
	for (i = 0; i < 4; i++) { cs[i] = mvf_cubic(blob[i]); cv[i] = cs[i]; }
	{
		qaws_arc_segment g;
		qaws_arc_desc ad;
		memset(&g, 0, sizeof(g));
		g.center[0] = (qaws_scalar)1.15; g.center[1] = (qaws_scalar)0.35; g.radius = (qaws_scalar)0.28;
		g.angle_start = (qaws_scalar)(2 * 3.14159265358979323846); g.angle_end = 0;
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
		qaws_curve_create_arc(&cs[4], &ad);
		cv[4] = cs[4];
	}
	p[0].curves = cv; p[0].curve_count = 4; p[0].closed = 1;
	p[1].curves = cv + 4; p[1].curve_count = 1; p[1].closed = 1;
	memset(&d, 0, sizeof(d));
	d.paths = p; d.path_count = 2;
	memset(&rd, 0, sizeof(rd));
	rd.paths = p; rd.path_count = 2; rd.fill_rule = QAWS_FILL_EVEN_ODD;
	if (qaws_mean_value_create(&d, &m) == QAWS_STATUS_OK && qaws_region_sampler_create(&rd, &rs) == QAWS_STATUS_OK)
	{
		nh = qaws_mean_value_get_handle_count(m);
		h = (qaws_vec2*)malloc(sizeof(qaws_vec2) * nh);
		mv = (qaws_vec2*)malloc(sizeof(qaws_vec2) * nh);
		for (i = 0; i < nh; i++)
		{
			qaws_mean_value_get_handle(m, i, &h[i]);
			mv[i].x = (qaws_scalar)(h[i].x + 0.15 * h[i].y * h[i].y);
			mv[i].y = (qaws_scalar)(h[i].y + 0.16 * h[i].x * h[i].x - 0.1 * h[i].x);
			if (i >= 12) mv[i].x -= (qaws_scalar)0.12;
		}
		for (panel = 0; panel < 2; panel++)
		{
			viewport v;
			qaws_vec2 const* hh = panel ? mv : h;
			unsigned int k;
			v.x0 = 24 + panel * (pw + 12); v.y0 = 104; v.w = pw; v.h = ph;
			v.xmin = -1.5; v.xmax = 3.1; v.ymin = -1.0; v.ymax = -1.0 + 4.6 * ph / pw;
			svg_panel(&s, &v, panel ? "handles moved: the inside carried along" : "at rest: 4 cubic Beziers and an arc (as rational quadratics), 20 handles");
			mvf_grid(&s, &v, m, rs, hh, gbox, 22);
			mvf_boundary(&s, &v, m, cs, 5, hh);
			for (k = 0; k < 2; k++)
			{
				double xy[2 * 13];
				unsigned int a = k ? 12 : 0, b = k ? nh : 12, j, n = 0;
				for (j = a; j < b; j++) { xy[2 * n] = vx(&v, hh[j].x); xy[2 * n + 1] = vy(&v, hh[j].y); n++; }
				xy[2 * n] = vx(&v, hh[a].x); xy[2 * n + 1] = vy(&v, hh[a].y); n++;
				svg_polyline(&s, xy, (int)n, "#8c959f", 1.0, 1.0, 1);
			}
			for (k = 0; k < nh; k++)
				svg_circle(&s, vx(&v, hh[k].x), vy(&v, hh[k].y), 4.2, k < 12 ? "#1a7f37" : "#8250df", "#ffffff");
		}
	}
	free(h); free(mv);
	qaws_region_sampler_destroy(rs); rs = NULL;
	qaws_mean_value_destroy(m); m = NULL;
	for (i = 0; i < 5; i++) qaws_curve_destroy(cs[i]);
	/* the transfinite interpolant: an L with a hole (the hole hot), and an ellipse */
	for (panel = 0; panel < 2; panel++)
	{
		viewport v;
		qaws_curve *c0 = NULL, *c1 = NULL;
		qaws_curve const* pv[2];
		unsigned int np = 1;
		double box[4];
		char const* label;
		if (panel == 0)
		{
			qaws_scalar L[12] = { 0, 0, 2, 0, 2, (qaws_scalar)0.8, (qaws_scalar)0.9, (qaws_scalar)0.8, (qaws_scalar)0.9, (qaws_scalar)1.6, 0, (qaws_scalar)1.6 };
			qaws_scalar H[8] = { (qaws_scalar)0.25, (qaws_scalar)0.25, (qaws_scalar)0.6, (qaws_scalar)0.25, (qaws_scalar)0.6, (qaws_scalar)0.55, (qaws_scalar)0.25, (qaws_scalar)0.55 };
			qaws_curve_create_polyline_2d(L, 6, 1, &c0);
			qaws_curve_create_polyline_2d(H, 4, 1, &c1);
			np = 2;
			box[0] = -0.05; box[1] = 2.05; box[2] = -0.05; box[3] = 1.65;
			label = "a concave polygon with a hole: 1 on the hole, a wave on the outside";
		}
		else
		{
			qaws_curve_create_ellipse_2d((qaws_scalar)0.1, (qaws_scalar)-0.2, (qaws_scalar)1.5, (qaws_scalar)0.8, (qaws_scalar)0.3, &c0);
			box[0] = -1.5; box[1] = 1.7; box[2] = -1.15; box[3] = 0.75;
			label = "a NURBS ellipse: five waves round its boundary";
		}
		pv[0] = c0; pv[1] = c1;
		p[0].curves = pv; p[0].curve_count = 1; p[0].closed = 1;
		p[1].curves = pv + 1; p[1].curve_count = 1; p[1].closed = 1;
		d.paths = p; d.path_count = np;
		rd.paths = p; rd.path_count = np;
		v.x0 = 24 + panel * (pw + 12); v.y0 = 104 + ph + 44; v.w = pw; v.h = ph;
		{
			double sx = (box[1] - box[0]) / (pw - 20), sy = (box[3] - box[2]) / (ph - 30), sc = sx > sy ? sx : sy;
			double cx = 0.5 * (box[0] + box[1]), cy = 0.5 * (box[2] + box[3]);
			v.xmin = cx - 0.5 * pw * sc; v.xmax = cx + 0.5 * pw * sc; v.ymin = cy - 0.5 * ph * sc; v.ymax = cy + 0.5 * ph * sc;
		}
		svg_panel(&s, &v, label);
		if (qaws_mean_value_create(&d, &m) == QAWS_STATUS_OK && qaws_region_sampler_create(&rd, &rs) == QAWS_STATUS_OK)
			mvf_field(&s, &v, m, rs, box, panel);
		for (i = 0; i < np; i++)
		{
			qaws_range rg = qaws_curve_get_parameter_range(pv[i]);
			qaws_vec2 bp[241];
			unsigned int k;
			for (k = 0; k <= 240; k++)
			{
				qaws_eval_result_2d e;
				qaws_curve_evaluate_2d(&e, pv[i], rg.min_value + (rg.max_value - rg.min_value) * k / 240, QAWS_EVAL_FLAG_POSITION);
				bp[k] = e.position;
			}
			mvf_points(&s, &v, bp, 241, "#1b1f24", 1.6);
		}
		qaws_region_sampler_destroy(rs); rs = NULL;
		qaws_mean_value_destroy(m); m = NULL;
		qaws_curve_destroy(c0);
		qaws_curve_destroy(c1);
	}
	svg_close(&s);
	printf("figure 15: mean value coordinates, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
