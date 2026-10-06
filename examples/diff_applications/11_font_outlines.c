/* ================================================================== */
/* 11. Font outline fitting: glyph contours to cubic B-splines         */
/* ================================================================== */

/*
 * A rendered glyph (photos/glyph_*.ppm, from examples/glyph_to_ppm.ps1) is
 * traced at its half-ink isoline; every contour becomes one closed cubic
 * B-spline. Corners are found on the trace and become triple knots (C0);
 * between corners the knots are spaced by arc length. The control points
 * descend sum |C(t_j) - q_j|^2 over the traced points q_j through the batch
 * adjoint, with the parameters t_j re-projected onto the curve at every
 * step (parameter correction).
 *
 * Closed curves: the knot gaps repeat around the loop and the first three
 * control points are repeated at the end, so the spline is periodic; the
 * adjoints of a repeated point are summed into its one variable.
 */

#define FO_MAX_CONTOURS 16
#define FO_MAX_CP 160
#define FO_MAX_Q 2400
#define FO_SPAN_PX 42.0
#define FO_CORNER_DEG 38.0

typedef struct fo_loop
{
	int n;                            /* unique control points (= knot gaps) */
	qaws_scalar p[FO_MAX_CP * 2];
	qaws_scalar knots[FO_MAX_CP + 7];
	int corner[FO_MAX_CP];            /* control point i interpolates a corner */
	int nq;
	double q[FO_MAX_Q * 2];           /* traced points */
	double t[FO_MAX_Q];               /* their curve parameters */
	double length;
	int corners;
} fo_loop;

static qaws_curve* fo_curve(fo_loop const* L)
{
	qaws_scalar ext[(FO_MAX_CP + 3) * 2];
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	int i;
	for (i = 0; i < L->n + 3; i++)
	{
		ext[2 * i] = L->p[2 * (i % L->n)];
		ext[2 * i + 1] = L->p[2 * (i % L->n) + 1];
	}
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = ext;
	d.control_point_count = (unsigned int)(L->n + 3);
	d.knots = L->knots;
	d.knot_count = (unsigned int)(L->n + 7);
	qaws_curve_create_bspline(&d, &c);
	return c;
}

/* Point of the closed polyline at arc length s (wrapped). */
static void fo_poly_at(double const* xy, double const* cum, int n, double total, double s, double* x, double* y)
{
	int lo = 0, hi = n;
	double f;
	s = fmod(s, total);
	if (s < 0) s += total;
	while (hi - lo > 1)
	{
		int mid = (lo + hi) / 2;
		if (cum[mid] <= s) lo = mid; else hi = mid;
	}
	{
		int a = lo, b = (lo + 1) % n;
		double seg = cum[lo + 1] - cum[lo];
		f = seg > 0 ? (s - cum[lo]) / seg : 0;
		*x = xy[2 * a] + f * (xy[2 * b] - xy[2 * a]);
		*y = xy[2 * a + 1] + f * (xy[2 * b + 1] - xy[2 * a + 1]);
	}
}

/* Builds the loop from a closed trace: corners, knot gaps, initial control
   points at the Greville abscissae, and the traced points to fit. */
static int fo_build(contour const* c, fo_loop* L)
{
	int n = c->n, i, k, start = 0, nc = 0, cidx[64];
	double* cum = (double*)malloc(sizeof(double) * (size_t)(n + 1));
	double* ang = (double*)calloc((size_t)n, sizeof(double));
	double* xy = (double*)malloc(sizeof(double) * 2 * (size_t)n);
	double total, gaps[FO_MAX_CP];
	int ng = 0;
	memset(L, 0, sizeof(*L));
	/* turning angle over about 5 px on each side */
	for (i = 0; i < n; i++)
	{
		int a = (i - 5 + n) % n, b = (i + 5) % n;
		double ux = c->xy[2 * i] - c->xy[2 * a], uy = c->xy[2 * i + 1] - c->xy[2 * a + 1];
		double vx_ = c->xy[2 * b] - c->xy[2 * i], vy_ = c->xy[2 * b + 1] - c->xy[2 * i + 1];
		double cr = ux * vy_ - uy * vx_, dt = ux * vx_ + uy * vy_;
		ang[i] = fabs(atan2(cr, dt)) * 180 / 3.14159265358979;
	}
	for (i = 0; i < n && nc < 64; i++)
	{
		int is_max = ang[i] > FO_CORNER_DEG;
		for (k = -6; k <= 6 && is_max; k++)
			if (k && (ang[(i + k + n) % n] > ang[i] || (ang[(i + k + n) % n] == ang[i] && k < 0)))
				is_max = 0;
		if (is_max)
			cidx[nc++] = i;
	}
	/* rotate the trace so it starts at the first corner */
	if (nc)
		start = cidx[0];
	for (i = 0; i < n; i++)
	{
		xy[2 * i] = c->xy[2 * ((i + start) % n)];
		xy[2 * i + 1] = c->xy[2 * ((i + start) % n) + 1];
	}
	for (k = 0; k < nc; k++)
		cidx[k] = (cidx[k] - start + n) % n;
	cum[0] = 0;
	for (i = 0; i < n; i++)
		cum[i + 1] = cum[i] + hypot(xy[2 * ((i + 1) % n)] - xy[2 * i], xy[2 * ((i + 1) % n) + 1] - xy[2 * i + 1]);
	total = cum[n];
	L->length = total;
	L->corners = nc;
	/* knot gaps: segments between corners at about FO_SPAN_PX each, two zero
	   gaps at every corner (a triple knot) */
	{
		int segs = nc ? nc : 1;
		for (k = 0; k < segs; k++)
		{
			double s0 = nc ? cum[cidx[k]] : 0, s1 = nc ? (k + 1 < nc ? cum[cidx[k + 1]] : total) : total;
			int ns = (int)((s1 - s0) / FO_SPAN_PX + 0.5), j;
			if (ns < 1) ns = 1;
			if (nc)
			{
				gaps[ng++] = 0;
				gaps[ng++] = 0;
			}
			for (j = 0; j < ns && ng < FO_MAX_CP - 2; j++)
				gaps[ng++] = (s1 - s0) / ns;
		}
	}
	L->n = ng;
	/* knots: the gap after u_i is gaps[(i - 3) mod n], so u_3 is arc length 0 */
	L->knots[0] = 0;
	for (i = 0; i < ng + 6; i++)
		L->knots[i + 1] = (qaws_scalar)(L->knots[i] + gaps[((i - 3) % ng + ng) % ng]);
	/* Greville abscissae */
	for (i = 0; i < ng; i++)
	{
		double g = (L->knots[i + 1] + L->knots[i + 2] + L->knots[i + 3]) / 3.0 - L->knots[3], x, y;
		fo_poly_at(xy, cum, n, total, g, &x, &y);
		L->p[2 * i] = (qaws_scalar)x;
		L->p[2 * i + 1] = (qaws_scalar)y;
		L->corner[i] = nc && L->knots[i + 1] == L->knots[i + 2] && L->knots[i + 2] == L->knots[i + 3];
	}
	/* traced points every ~1.5 px, with their arc-length parameters */
	{
		int m = (int)(total / 1.5);
		if (m > FO_MAX_Q) m = FO_MAX_Q;
		for (i = 0; i < m; i++)
		{
			double s = total * i / m;
			fo_poly_at(xy, cum, n, total, s, &L->q[2 * i], &L->q[2 * i + 1]);
			L->t[i] = L->knots[3] + s;
		}
		L->nq = m;
	}
	free(cum);
	free(ang);
	free(xy);
	return ng >= 4;
}

/* Parameter correction: a few Newton steps of |C(t) - q|^2 per point. */
static void fo_project(qaws_curve const* c, fo_loop* L)
{
	double lo = L->knots[3], hi = L->knots[L->n + 3], per = hi - lo;
	int i, it;
	for (i = 0; i < L->nq; i++)
		for (it = 0; it < 3; it++)
		{
			qaws_eval_result_2d e;
			double dx, dy, g, H, t = L->t[i];
			qaws_curve_evaluate_2d(c, (qaws_scalar)t, QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2, &e);
			dx = e.position.x - L->q[2 * i];
			dy = e.position.y - L->q[2 * i + 1];
			g = dx * e.d1.x + dy * e.d1.y;
			H = e.d1.x * e.d1.x + e.d1.y * e.d1.y + dx * e.d2.x + dy * e.d2.y;
			if (H > 1e-9)
				t -= g / H;
			t = lo + fmod(fmod(t - lo, per) + per, per);
			L->t[i] = t;
		}
}

/* Mean squared distance and its gradient on the unique control points. */
static double fo_energy(qaws_curve const* c, fo_loop const* L, qaws_scalar* grad, double* max_err)
{
	static qaws_scalar ts[FO_MAX_Q];
	static qaws_curve_jet_2d prim[FO_MAX_Q], tan[FO_MAX_Q], bar[FO_MAX_Q];
	qaws_scalar ext[(FO_MAX_CP + 3) * 2];
	double e = 0, mx = 0;
	int i;
	for (i = 0; i < L->nq; i++)
		ts[i] = (qaws_scalar)L->t[i];
	qaws_curve_eval_batch_tangent_2d(NULL, c, ts, NULL, (unsigned int)L->nq, QAWS_EVAL_FLAG_POSITION, NULL, prim, tan);
	for (i = 0; i < L->nq; i++)
	{
		double dx = prim[i].d[0].x - L->q[2 * i], dy = prim[i].d[0].y - L->q[2 * i + 1], d2 = dx * dx + dy * dy;
		e += d2 / L->nq;
		if (d2 > mx) mx = d2;
		memset(&bar[i], 0, sizeof(bar[i]));
		bar[i].d[0].x = (qaws_scalar)(2 * dx / L->nq);
		bar[i].d[0].y = (qaws_scalar)(2 * dy / L->nq);
		bar[i].channels = QAWS_EVAL_FLAG_POSITION;
	}
	if (max_err)
		*max_err = sqrt(mx);
	if (grad)
	{
		qaws_field_view fv;
		qaws_diff_views views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, ext, (unsigned int)(L->n + 3), 2);
		memset(ext, 0, sizeof(qaws_scalar) * 2 * (size_t)(L->n + 3));
		qaws_curve_eval_batch_adjoint_2d(NULL, c, ts, (unsigned int)L->nq, QAWS_EVAL_FLAG_POSITION, bar, &views, NULL);
		memset(grad, 0, sizeof(qaws_scalar) * 2 * (size_t)L->n);
		for (i = 0; i < L->n + 3; i++)
		{
			grad[2 * (i % L->n)] += ext[2 * i];
			grad[2 * (i % L->n) + 1] += ext[2 * i + 1];
		}
	}
	return e;
}

static void fo_fit(fo_loop* L, double* out_mean, double* out_max)
{
	qaws_scalar g[FO_MAX_CP * 2];
	static adam opt;
	double e = 0, mx = 0;
	int it;
	memset(&opt, 0, sizeof(opt));
	for (it = 0; it < 400; it++)
	{
		qaws_curve* c = fo_curve(L);
		if (it % 4 == 0)
			fo_project(c, L);
		e = fo_energy(c, L, g, &mx);
		qaws_curve_destroy(c);
		adam_step(&opt, L->p, g, 2 * L->n, it < 250 ? 0.08 : 0.02);
	}
	{
		qaws_curve* c = fo_curve(L);
		fo_project(c, L);
		e = fo_energy(c, L, NULL, &mx);
		qaws_curve_destroy(c);
	}
	*out_mean = sqrt(e);
	*out_max = mx;
}

static void fo_curve_path(svg* s, imgmap const* m, qaws_curve const* c, int start)
{
	qaws_range r = qaws_curve_get_parameter_range(c);
	int k, n = 400;
	for (k = 0; k <= n; k++)
	{
		qaws_eval_result_2d e;
		double x, y;
		qaws_curve_evaluate_2d(c, (qaws_scalar)(r.min_value + (r.max_value - r.min_value) * k / n), QAWS_EVAL_FLAG_POSITION, &e);
		map_pt(m, e.position.x, e.position.y, &x, &y);
		fprintf(s->f, "%s%.2f %.2f", (k == 0 && start) ? "M" : (k == 0 ? "M" : "L"), x, y);
	}
	fprintf(s->f, "Z");
}

static void app_font_outlines(void)
{
	static char const* const names[6] = { "amp", "g", "R", "S", "at", "Q" };
	static char const* const labels[6] = { "&amp;", "g", "R", "S", "@", "Q" };
	static contour cs[FO_MAX_CONTOURS];
	static fo_loop loops[FO_MAX_CONTOURS];
	svg s;
	int gi, total_q = 0, total_cp = 0;
	double worst = 0;
	svg_open(&s, "showcase/app11_font_outlines.svg", 1500, 900, "Font outlines: glyph traces fitted by closed cubic B-splines",
		"Georgia glyphs rendered at 400 px, traced at the half-ink isoline; corners become triple knots, knots elsewhere follow arc length. Control points descend the squared distance through the batch adjoint, with parameter correction.");
	for (gi = 0; gi < 6; gi++)
	{
		char path[512], buf[200];
		image img;
		int nc, i, ncp = 0, nq = 0, ncorner = 0;
		double gmean = 0, gmax = 0;
		viewport v;
		imgmap m;
		sprintf(path, "%s/glyph_%s.ppm", g_photos, names[gi]);
		if (!image_load_ppm(path, &img))
		{
			printf("11_font_outlines: %s missing (examples/glyph_to_ppm.ps1)\n", path);
			continue;
		}
		blur(img.lum, img.w, img.h, 0.8);
		nc = extract_contours(img.lum, img.w, img.h, 0.5, 30, cs, FO_MAX_CONTOURS);
		v.x0 = 20 + (gi % 3) * 490;
		v.y0 = 70 + (gi / 3) * 415;
		v.w = 470;
		v.h = 400;
		v.xmin = 0; v.xmax = 1; v.ymin = 0; v.ymax = 1;
		sprintf(buf, "\"%s\"", labels[gi]);
		svg_panel(&s, &v, buf);
		{
			/* panel framed on the glyph's bounding box */
			double x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
			int k;
			for (i = 0; i < nc; i++)
				for (k = 0; k < cs[i].n; k++)
				{
					if (cs[i].xy[2 * k] < x0) x0 = cs[i].xy[2 * k];
					if (cs[i].xy[2 * k] > x1) x1 = cs[i].xy[2 * k];
					if (cs[i].xy[2 * k + 1] < y0) y0 = cs[i].xy[2 * k + 1];
					if (cs[i].xy[2 * k + 1] > y1) y1 = cs[i].xy[2 * k + 1];
				}
			m.scale = 320.0 / (y1 - y0);
			if ((x1 - x0) * m.scale > 430) m.scale = 430.0 / (x1 - x0);
			m.ox = x0;
			m.oy = y0;
			m.x0 = v.x0 + (v.w - (x1 - x0) * m.scale) / 2;
			m.y0 = v.y0 + 32 + (330 - (y1 - y0) * m.scale) / 2;
		}
		/* fit every closed contour */
		fprintf(s.f, "<path fill=\"#d0d7de\" fill-rule=\"evenodd\" stroke=\"#24292f\" stroke-width=\"1.4\" d=\"");
		for (i = 0; i < nc; i++)
		{
			double mean, mx;
			qaws_curve* c;
			if (!cs[i].closed || !fo_build(&cs[i], &loops[i]))
			{
				loops[i].n = 0;
				continue;
			}
			fo_fit(&loops[i], &mean, &mx);
			c = fo_curve(&loops[i]);
			fo_curve_path(&s, &m, c, i == 0);
			qaws_curve_destroy(c);
			ncp += loops[i].n;
			nq += cs[i].n;
			ncorner += loops[i].corners;
			gmean += mean * mean * loops[i].nq;
			if (mx > gmax) gmax = mx;
		}
		fprintf(s.f, "\"/>\n");
		/* control polygons and corners */
		for (i = 0; i < nc; i++)
		{
			fo_loop const* L = &loops[i];
			double xy[2 * (FO_MAX_CP + 1)];
			int k, tot = 0;
			if (!L->n)
				continue;
			for (k = 0; k <= L->n; k++)
				map_pt(&m, L->p[2 * (k % L->n)], L->p[2 * (k % L->n) + 1], &xy[2 * k], &xy[2 * k + 1]);
			svg_polyline(&s, xy, L->n + 1, "#0969da", 0.8, 0.7, 0);
			for (k = 0; k < L->n; k++)
				if (L->corner[k])
					fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"6\" height=\"6\" fill=\"#cf222e\"/>\n", xy[2 * k] - 3, xy[2 * k + 1] - 3);
				else
					svg_circle(&s, xy[2 * k], xy[2 * k + 1], 2.4, "#ffffff", "#0969da");
			tot += L->nq;
			(void)tot;
		}
		{
			int nqs = 0;
			for (i = 0; i < nc; i++)
				nqs += loops[i].n ? loops[i].nq : 0;
			gmean = nqs ? sqrt(gmean / nqs) : 0;
		}
		sprintf(buf, "%d trace points -> %d control points, %d corners; rms %.2f px, max %.2f px", nq, ncp, ncorner, gmean, gmax);
		svg_text(&s, v.x0 + 10, v.y0 + v.h - 10, 11, "#57606a", "start", buf);
		printf("11_font_outlines: %s: %s\n", names[gi], buf);
		total_q += nq;
		total_cp += ncp;
		if (gmax > worst) worst = gmax;
		for (i = 0; i < nc; i++)
			free(cs[i].xy);
		image_free(&img);
	}
	svg_close(&s);
	printf("11_font_outlines: %d trace points -> %d control points, worst max error %.2f px\n", total_q, total_cp, worst);
}
