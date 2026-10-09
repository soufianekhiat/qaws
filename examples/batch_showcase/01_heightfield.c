/* Figure 1: contour lines and gradient lines of a heightfield, every
   crossing in one batched call; the time of the batch against every
   contour / gradient pair as the number of curves grows. */

#define HF_BUMPS 6

static double const hf_bump[HF_BUMPS][4] = {
	/* cx, cy, height, width */
	{ -0.45, -0.35, 1.0, 0.42 },
	{ 0.50, 0.40, 0.85, 0.36 },
	{ 0.35, -0.55, -0.7, 0.30 },
	{ -0.50, 0.55, 0.6, 0.28 },
	{ 0.05, 0.05, -0.45, 0.25 },
	{ 0.75, -0.10, 0.4, 0.22 },
};

static double hf_height(double x, double y, double* gx, double* gy)
{
	double h = 0;
	int i;
	if (gx) *gx = 0;
	if (gy) *gy = 0;
	for (i = 0; i < HF_BUMPS; i++)
	{
		double dx = x - hf_bump[i][0], dy = y - hf_bump[i][1], s2 = hf_bump[i][3] * hf_bump[i][3];
		double e = hf_bump[i][2] * exp(-(dx * dx + dy * dy) / s2);
		h += e;
		if (gx) *gx += -2 * dx / s2 * e;
		if (gy) *gy += -2 * dy / s2 * e;
	}
	return h;
}

typedef struct hf_scene
{
	qaws_curve** curves;
	unsigned int* family;   /* 0 contour, 1 gradient */
	unsigned int count, cap;
} hf_scene;

/* cubic B-spline with the polyline points (every `step`-th) as control points */
static void hf_add(hf_scene* sc, double const* pts, unsigned int n, unsigned int step, unsigned int family)
{
	qaws_scalar* cps;
	qaws_scalar* kn;
	unsigned int m = 0, i, nk;
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	if (n < 2)
		return;
	cps = (qaws_scalar*)malloc((n / step + 3) * 2 * sizeof(qaws_scalar));
	for (i = 0; i < n; i += step, m++)
	{
		cps[2 * m] = (qaws_scalar)pts[2 * i];
		cps[2 * m + 1] = (qaws_scalar)pts[2 * i + 1];
	}
	if ((n - 1) % step)
	{
		cps[2 * m] = (qaws_scalar)pts[2 * (n - 1)];
		cps[2 * m + 1] = (qaws_scalar)pts[2 * (n - 1) + 1];
		m++;
	}
	if (m < 4)
	{
		free(cps);
		return;
	}
	nk = m + 4;
	kn = (qaws_scalar*)malloc(nk * sizeof(qaws_scalar));
	for (i = 0; i < nk; i++)
		kn[i] = (qaws_scalar)(i < 4 ? 0 : (i >= m ? m - 3 : i - 3));
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_2D;
	d.degree = 3;
	d.control_points = cps;
	d.control_point_count = m;
	d.knots = kn;
	d.knot_count = nk;
	if (qaws_curve_create_bspline(&d, &c) == QAWS_STATUS_OK)
	{
		if (sc->count == sc->cap)
		{
			sc->cap = sc->cap ? sc->cap * 2 : 64;
			sc->curves = (qaws_curve**)realloc(sc->curves, sc->cap * sizeof(qaws_curve*));
			sc->family = (unsigned int*)realloc(sc->family, sc->cap * sizeof(unsigned int));
		}
		sc->curves[sc->count] = c;
		sc->family[sc->count++] = family;
	}
	free(cps);
	free(kn);
}

#define HF_GRID 160

/* marching squares at level L; crossed grid edges are chained through the
   cells they bound, so every contour comes out as one polyline */
static void hf_contours(hf_scene* sc, double const* f, double L)
{
	unsigned int const G = HF_GRID, HN = (G + 1) * G, NE = 2 * HN;
	int* link = (int*)malloc(NE * 2 * sizeof(int));
	unsigned char* used = (unsigned char*)calloc(NE, 1);
	double* pts = (double*)malloc(NE * 2 * sizeof(double));
	unsigned int i, j, e;
#define HF_F(i, j) (f[(j) * (G + 1) + (i)] - L)
#define HF_H(i, j) ((j) * G + (i))
#define HF_V(i, j) (HN + (j) * (G + 1) + (i))
	for (e = 0; e < NE * 2; e++)
		link[e] = -1;
	for (j = 0; j < G; j++)
		for (i = 0; i < G; i++)
		{
			/* corners 0 (i,j) 1 (i+1,j) 2 (i+1,j+1) 3 (i,j+1); edges 0 bottom 1 right 2 top 3 left */
			double v[4];
			int in[4], ed[4], cr[4], nc = 0, k, pair[4];
			v[0] = HF_F(i, j); v[1] = HF_F(i + 1, j); v[2] = HF_F(i + 1, j + 1); v[3] = HF_F(i, j + 1);
			ed[0] = HF_H(i, j); ed[1] = HF_V(i + 1, j); ed[2] = HF_H(i, j + 1); ed[3] = HF_V(i, j);
			for (k = 0; k < 4; k++)
				in[k] = v[k] >= 0;
			for (k = 0; k < 4; k++)
				if (in[k] != in[(k + 1) & 3])
					cr[nc++] = k;
			if (nc == 2)
			{
				pair[0] = cr[0];
				pair[1] = cr[1];
			}
			else if (nc == 4)
			{
				/* saddle: the centre decides which corners connect */
				int c_in = (v[0] + v[1] + v[2] + v[3]) / 4 >= 0;
				if (c_in == in[0]) { pair[0] = 0; pair[1] = 1; pair[2] = 2; pair[3] = 3; }
				else { pair[0] = 3; pair[1] = 0; pair[2] = 1; pair[3] = 2; }
			}
			for (k = 0; k + 1 < nc; k += 2)
			{
				int a = ed[pair[k]], b = ed[pair[k + 1]];
				link[2 * a + (link[2 * a] >= 0)] = b;
				link[2 * b + (link[2 * b] >= 0)] = a;
			}
		}
	/* chain: open contours first (an end on the border), then loops */
	{
		int pass;
		double* poly = (double*)malloc(NE * 2 * sizeof(double));
		for (pass = 0; pass < 2; pass++)
			for (e = 0; e < NE; e++)
			{
				int cur, prev = -1, n = 0;
				if (used[e] || link[2 * e] < 0 || (pass == 0 && link[2 * e + 1] >= 0))
					continue;
				cur = (int)e;
				while (cur >= 0 && !used[cur])
				{
					int nx;
					unsigned int ci, cj, di, dj;
					double a, b, t;
					used[cur] = 1;
					if ((unsigned int)cur < HN) { ci = cur % G; cj = cur / G; di = ci + 1; dj = cj; }
					else { ci = (cur - HN) % (G + 1); cj = (cur - HN) / (G + 1); di = ci; dj = cj + 1; }
					a = HF_F(ci, cj);
					b = HF_F(di, dj);
					t = a / (a - b);
					poly[2 * n] = -1 + 2.0 * (ci + t * ((double)di - ci)) / G;
					poly[2 * n + 1] = -1 + 2.0 * (cj + t * ((double)dj - cj)) / G;
					n++;
					nx = link[2 * cur] != prev ? link[2 * cur] : link[2 * cur + 1];
					prev = cur;
					cur = nx;
				}
				if (cur == (int)e)   /* closed: repeat the first point */
				{
					poly[2 * n] = poly[0];
					poly[2 * n + 1] = poly[1];
					n++;
				}
				if (n >= 12)
					hf_add(sc, poly, (unsigned int)n, 4, 0);
			}
		free(poly);
	}
#undef HF_F
#undef HF_H
#undef HF_V
	free(link);
	free(used);
	free(pts);
}

/* gradient line through (x, y): RK4 along +/- grad h / |grad h| */
static void hf_gradient(hf_scene* sc, double x, double y)
{
	enum { STEPS = 500 };
	static double buf[2][2 * STEPS + 2], line[4 * STEPS + 4];
	unsigned int cnt[2], dir, i, n = 0;
	double const ds = 0.008;
	for (dir = 0; dir < 2; dir++)
	{
		double px = x, py = y, sg = dir ? -1 : 1;
		cnt[dir] = 0;
		for (i = 0; i < STEPS; i++)
		{
			double kx[4], ky[4], gx, gy, l, qx = px, qy = py;
			int s;
			for (s = 0; s < 4; s++)
			{
				hf_height(qx, qy, &gx, &gy);
				l = sqrt(gx * gx + gy * gy);
				if (l < 2e-3)
					break;
				kx[s] = sg * gx / l;
				ky[s] = sg * gy / l;
				qx = px + (s == 2 ? ds : ds / 2) * kx[s];
				qy = py + (s == 2 ? ds : ds / 2) * ky[s];
			}
			if (s < 4)
				break;
			px += ds / 6 * (kx[0] + 2 * kx[1] + 2 * kx[2] + kx[3]);
			py += ds / 6 * (ky[0] + 2 * ky[1] + 2 * ky[2] + ky[3]);
			if (fabs(px) > 1 || fabs(py) > 1)
				break;
			buf[dir][2 * cnt[dir]] = px;
			buf[dir][2 * cnt[dir] + 1] = py;
			cnt[dir]++;
		}
	}
	for (i = cnt[1]; i-- > 0;)
	{
		line[2 * n] = buf[1][2 * i];
		line[2 * n + 1] = buf[1][2 * i + 1];
		n++;
	}
	line[2 * n] = x;
	line[2 * n + 1] = y;
	n++;
	for (i = 0; i < cnt[0]; i++)
	{
		line[2 * n] = buf[0][2 * i];
		line[2 * n + 1] = buf[0][2 * i + 1];
		n++;
	}
	if (n >= 12)
		hf_add(sc, line, n, 3, 1);
}

/* n contour levels and about n gradient lines (seeds on a jittered grid) */
static void hf_build(hf_scene* sc, unsigned int n, double const* f)
{
	unsigned int k, s = (unsigned int)ceil(sqrt((double)n)), i, j;
	memset(sc, 0, sizeof(*sc));
	for (k = 0; k < n; k++)
		hf_contours(sc, f, -0.65 + 1.55 * (k + 0.5) / n);
	for (j = 0; j < s; j++)
		for (i = 0; i < s; i++)
		{
			double jx = 0.3 * sin(12.9898 * (i + 1) + 78.233 * (j + 1)), jy = 0.3 * cos(39.3468 * (i + 1) + 11.135 * (j + 1));
			hf_gradient(sc, -0.92 + 1.84 * (i + 0.5 + jx) / s, -0.92 + 1.84 * (j + 0.5 + jy) / s);
		}
}

static void hf_free(hf_scene* sc)
{
	unsigned int i;
	for (i = 0; i < sc->count; i++)
		qaws_curve_destroy(sc->curves[i]);
	free(sc->curves);
	free(sc->family);
}

static double hf_now(void)
{
	return (double)clock() / CLOCKS_PER_SEC;
}

/* every crossing by the batch; and by the pairwise call over all contour /
   gradient pairs (pairwise skipped when `pairwise` is 0) */
static void hf_time(hf_scene const* sc, int pairwise, double* tb, double* tp, unsigned int* nb, unsigned int* np,
	qaws_curve_batch_hit_2d** hits, qaws_curve_batch_stats* st)
{
	qaws_curve_batch_desc d;
	unsigned int cap = 1 << 16, i, j;
	double t0;
	*hits = (qaws_curve_batch_hit_2d*)malloc(cap * sizeof(qaws_curve_batch_hit_2d));
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_curve const* const*)sc->curves;
	d.curve_count = sc->count;
	d.families = sc->family;
	/* repeated for at least 50 ms: clock() ticks in milliseconds */
	{
		unsigned int reps = 0;
		t0 = hf_now();
		do
		{
			qaws_curve_batch_find_intersections_2d(&d, *hits, cap, nb, st);
			reps++;
		} while (hf_now() - t0 < 0.05);
		*tb = (hf_now() - t0) / reps;
	}
	*np = 0;
	*tp = 0;
	if (!pairwise)
		return;
	t0 = hf_now();
	for (i = 0; i < sc->count; i++)
		for (j = i + 1; j < sc->count; j++)
			if (sc->family[i] != sc->family[j])
			{
				qaws_intersection_2d buf[64];
				unsigned int k = 0;
				qaws_curve_find_intersections_2d(sc->curves[i], sc->curves[j], buf, 64, &k);
				*np += k;
			}
	*tp = hf_now() - t0;
}

static void demo_heightfield(void)
{
	static double f[(HF_GRID + 1) * (HF_GRID + 1)];
	unsigned int const sizes[5] = { 6, 12, 24, 48, 96 };
	double tb[5], tp[5], fig_tb, fig_tp;
	unsigned int curves[5], i, j, k;
	hf_scene sc;
	qaws_curve_batch_hit_2d* hits = NULL;
	qaws_curve_batch_stats st;
	unsigned int nb = 0, np = 0;
	svg s;
	viewport v, p;
	char sub[320], buf[160];
	double* xy = (double*)malloc(2 * 400 * sizeof(double));

	for (j = 0; j <= HF_GRID; j++)
		for (i = 0; i <= HF_GRID; i++)
			f[j * (HF_GRID + 1) + i] = hf_height(-1 + 2.0 * i / HF_GRID, -1 + 2.0 * j / HF_GRID, NULL, NULL);

	/* scaling: the pairwise loop only up to 48, it is the slow side */
	for (k = 0; k < 5; k++)
	{
		qaws_curve_batch_hit_2d* h = NULL;
		qaws_curve_batch_stats st2;
		hf_build(&sc, sizes[k], f);
		hf_time(&sc, sizes[k] <= 48, &tb[k], &tp[k], &nb, &np, &h, &st2);
		curves[k] = sc.count;
		printf("    heightfield %2u levels: %4u curves, batch %6u hits %.5f s (%u segments, %u candidates)", sizes[k], sc.count, nb, tb[k],
			st2.segment_count, st2.candidate_count);
		if (sizes[k] <= 48)
			printf(", pairwise %6u hits %.3f s, x%.0f", np, tp[k], tp[k] / (tb[k] > 1e-4 ? tb[k] : 1e-4));
		printf("\n");
		free(h);
		hf_free(&sc);
	}

	/* the figure: 24 levels */
	hf_build(&sc, 24, f);
	hf_time(&sc, 0, &fig_tb, &fig_tp, &nb, &np, &hits, &st);
	sprintf(sub, "%u contour lines and gradient lines of a heightfield as cubic B-splines; %u crossings in one call (%.1f ms, %u chord segments, %u candidate pairs)",
		sc.count, nb, fig_tb * 1000, st.segment_count, st.candidate_count);
	if (!svg_open(&s, "showcase/batch1_heightfield.svg", 1400, 760, "Batched intersections: contours x gradient lines", sub))
		return;
	v.x0 = 30; v.y0 = 80; v.w = 660; v.h = 660; v.xmin = -1; v.xmax = 1; v.ymin = -1; v.ymax = 1;
	svg_panel(&s, &v, NULL);
	for (j = 0; j < 66; j++)
		for (i = 0; i < 66; i++)
		{
			char col[16];
			double x = -1 + 2.0 * (i + 0.5) / 66, y = -1 + 2.0 * (j + 0.5) / 66, t = (hf_height(x, y, NULL, NULL) + 0.75) / 1.75;
			heat(t < 0 ? 0 : (t > 1 ? 1 : t), col);
			fprintf(s.f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"10.3\" height=\"10.3\" fill=\"%s\" fill-opacity=\"0.28\"/>\n",
				vx(&v, -1 + 2.0 * i / 66), vy(&v, -1 + 2.0 * (j + 1) / 66), col);
		}
	for (k = 0; k < sc.count; k++)
	{
		curve_polyline(sc.curves[k], &v, xy, 400);
		if (sc.family[k] == 0)
			svg_polyline(&s, xy, 400, "#24292f", 1.1, 0.85, 0);
		else
			svg_polyline(&s, xy, 400, "#0969da", 1.0, 0.75, 0);
	}
	for (k = 0; k < nb && k < (1u << 16); k++)
		svg_circle(&s, vx(&v, hits[k].position.x), vy(&v, hits[k].position.y), 2.3, "#cf222e", "#ffffff");
	svg_text(&s, v.x0 + 10, v.y0 + 20, 13, "#24292f", "start", "contours (black, family 0) x gradient lines (blue, family 1): red = crossing");

	/* scaling plot, log-log */
	p.x0 = 760; p.y0 = 80; p.w = 610; p.h = 620;
	p.xmin = log10(10.0); p.xmax = log10(2000.0); p.ymin = -4; p.ymax = 2;
	svg_panel(&s, &p, "time vs number of curves (log-log)");
	for (k = 0; k <= 6; k++)
	{
		double y = p.ymin + k;
		svg_line(&s, p.x0, vy(&p, y), p.x0 + p.w, vy(&p, y), "#d0d7de", 1, 1);
		sprintf(buf, "%g s", pow(10, y));
		svg_text(&s, p.x0 + 6, vy(&p, y) - 4, 11, "#57606a", "start", buf);
	}
	for (k = 1; k <= 3; k++)
	{
		double x = (double)k;
		svg_line(&s, vx(&p, x), p.y0, vx(&p, x), p.y0 + p.h, "#d0d7de", 1, 1);
		sprintf(buf, "%g curves", pow(10, x));
		svg_text(&s, vx(&p, x), p.y0 + p.h + 16, 11, "#57606a", "middle", buf);
	}
	{
		double pb[10], pp[8];
		unsigned int npp = 0;
		for (k = 0; k < 5; k++)
		{
			pb[2 * k] = vx(&p, log10((double)curves[k]));
			pb[2 * k + 1] = vy(&p, log10(tb[k] > 1e-4 ? tb[k] : 1e-4));
			if (sizes[k] <= 48)
			{
				pp[2 * npp] = pb[2 * k];
				pp[2 * npp + 1] = vy(&p, log10(tp[k]));
				npp++;
			}
		}
		svg_polyline(&s, pp, (int)npp, "#8250df", 2.2, 1, 0);
		svg_polyline(&s, pb, 5, "#1a7f37", 2.2, 1, 0);
		for (k = 0; k < npp; k++)
			svg_circle(&s, pp[2 * k], pp[2 * k + 1], 4, "#8250df", "#ffffff");
		for (k = 0; k < 5; k++)
			svg_circle(&s, pb[2 * k], pb[2 * k + 1], 4, "#1a7f37", "#ffffff");
		svg_text(&s, p.x0 + 20, p.y0 + 46, 13, "#8250df", "start", "pairwise: every contour / gradient pair");
		svg_text(&s, p.x0 + 20, p.y0 + 66, 13, "#1a7f37", "start", "batch: one flattening + one grid for all curves");
		k = 3;
		sprintf(buf, "%u curves: x%.0f faster", curves[k], tp[k] / (tb[k] > 1e-4 ? tb[k] : 1e-4));
		svg_text(&s, pb[2 * k] + 8, pb[2 * k + 1] + 18, 12, "#1a7f37", "start", buf);
	}
	svg_close(&s);
	printf("  -> showcase/batch1_heightfield.svg\n");
	free(hits);
	free(xy);
	hf_free(&sc);
}
