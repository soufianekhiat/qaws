/* Figure 35: the hierarchical mode of the fractional Sobolev descent
   (qaws_repulsive_curves_solve with theta > 0, R4): Barnes-Hut energy and
   gradient, the metric applied from Barnes-Hut lists, GMRES with the
   Laplacian-sandwich preconditioner. A lumpy (2, 5) torus knot of 600
   edges relaxed with its total length held; the time per step against
   size, dense against hierarchical.  */

static void s35_knot(qaws_scalar* x, unsigned int* e, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * 3.14159265358979323846 * i / n, r = 1 + 0.45 * cos(5 * t) + 0.12 * sin(13 * t);
		x[3 * i] = (qaws_scalar)(r * cos(2 * t)); x[3 * i + 1] = (qaws_scalar)(r * sin(2 * t)); x[3 * i + 2] = (qaws_scalar)(0.45 * sin(5 * t) + 0.08 * cos(9 * t));
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % n;
	}
}

/* the closed curve (vertices in order) seen from above, y tilted by z,
   fitted to the panel: runs of 8 edges drawn far to near, the near ones
   darker, each with a rim so that crossings read */
static void s35_draw(svg* s, viewport const* v, qaws_scalar const* x, unsigned int n, int const* rgb)
{
	enum { RUN = 8 };
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, sc, zlo = 1e300, zhi = -1e300, xy[2 * (RUN + 1)];
	unsigned int nr = (n + RUN - 1) / RUN, i, k, *order = (unsigned int*)malloc(sizeof(unsigned int) * nr);
	double* depth = (double*)malloc(sizeof(double) * nr);
	for (i = 0; i < n; i++)
	{
		double py = (double)x[3 * i + 1] * 0.94 - 0.34 * (double)x[3 * i + 2];
		lo[0] = fmin(lo[0], (double)x[3 * i]); hi[0] = fmax(hi[0], (double)x[3 * i]); lo[1] = fmin(lo[1], py); hi[1] = fmax(hi[1], py);
		zlo = fmin(zlo, (double)x[3 * i + 2]); zhi = fmax(zhi, (double)x[3 * i + 2]);
	}
	sc = fmin((v->w - 40) / (hi[0] - lo[0]), (v->h - 60) / (hi[1] - lo[1]));
	for (i = 0; i < nr; i++)
	{
		depth[i] = 0;
		for (k = 0; k <= RUN; k++) depth[i] += (double)x[3 * ((RUN * i + k) % n) + 2];
		order[i] = i;
	}
	for (i = 1; i < nr; i++)
	{
		unsigned int key = order[i];
		k = i;
		while (k > 0 && depth[order[k - 1]] > depth[key]) { order[k] = order[k - 1]; k--; }
		order[k] = key;
	}
	for (i = 0; i < nr; i++)
	{
		unsigned int r = order[i], m = 0;
		double t = zhi > zlo ? (depth[r] / (RUN + 1) - zlo) / (zhi - zlo) : 1, sh = 0.55 + 0.45 * t;
		char col[32];
		for (k = 0; k <= RUN && RUN * r + k <= n; k++, m++)
		{
			qaws_scalar const* p = x + 3 * ((RUN * r + k) % n);
			xy[2 * m] = v->x0 + v->w / 2 + ((double)p[0] - 0.5 * (lo[0] + hi[0])) * sc;
			xy[2 * m + 1] = v->y0 + 20 + v->h / 2 - (((double)p[1] * 0.94 - 0.34 * (double)p[2]) - 0.5 * (lo[1] + hi[1])) * sc;
		}
		sprintf(col, "rgb(%d,%d,%d)", (int)(rgb[0] * sh), (int)(rgb[1] * sh), (int)(rgb[2] * sh));
		if (m > 3) svg_polyline(s, xy + 2, (int)m - 2, "#f6f8fa", 6, 1, 0);   /* (the rim inside the run: its joins stay drawn) */
		svg_polyline(s, xy, (int)m, col, 2.2, 1, 0);
	}
	free(order); free(depth);
}

/* seconds per step at n edges (two steps) */
static double s35_step(unsigned int n, double theta)
{
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * n), it;
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solver* s = NULL;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_report r;
	clock_t c0;
	s35_knot(x, e, n);
	memset(&d, 0, sizeof(d));
	d.dimension = 3; d.vertex_count = n; d.edges = e; d.edge_count = n;
	qaws_repulsive_curves_create(&d, &c);
	memset(&sd, 0, sizeof(sd));
	sd.hold_total_length = 1; sd.theta = (qaws_scalar)theta;
	qaws_repulsive_curves_solver_create(c, &sd, x, &s);
	c0 = clock();
	for (it = 0; it < 2; it++) qaws_repulsive_curves_solver_step(s, x, &r);
	qaws_repulsive_curves_solver_destroy(s);
	qaws_repulsive_curves_destroy(c);
	free(x); free(e);
	return (double)(clock() - c0) / CLOCKS_PER_SEC / 2;
}

static void demo_repulsive_hsolve(void)
{
	enum { N = 600 };
	clock_t t0 = clock();
	svg s;
	double const pw = 330, ph = 330;
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * N), *x0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * N);
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * N);
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_report r;
	viewport v;
	static int const blue[3] = { 9, 105, 218 }, grey[3] = { 120, 128, 138 };
	char label[200];
	unsigned int k;
	if (!svg_open(&s, "showcase/b2d35_repulsive_hsolve.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + ph + 60),
		"Repulsive curves, hierarchical",
		"theta > 0: Barnes-Hut energy, gradient and metric, GMRES; no |V|^2 storage"))
		return;
	s35_knot(x, e, N);
	memcpy(x0, x, sizeof(qaws_scalar) * 3 * N);
	memset(&d, 0, sizeof(d));
	d.dimension = 3; d.vertex_count = N; d.edges = e; d.edge_count = N;
	qaws_repulsive_curves_create(&d, &c);
	memset(&sd, 0, sizeof(sd));
	sd.hold_total_length = 1; sd.max_iterations = 25; sd.theta = (qaws_scalar)0.25;
	memset(&r, 0, sizeof(r));
	qaws_repulsive_curves_solve(c, &sd, x, &r);
	qaws_repulsive_curves_destroy(c);
	for (k = 0; k < 2; k++)
	{
		v.x0 = 24 + k * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1; v.xmax = v.ymax = 1;
		sprintf(label, k ? "after %u steps: E %.1f, length held" : "a lumpy (2, 5) torus knot, %u edges: E %.1f", k ? r.iterations : N,
			k ? (double)r.energy : (double)r.initial_energy);
		svg_panel(&s, &v, label);
		s35_draw(&s, &v, k ? x : x0, N, k ? blue : grey);
	}
	/* time per step against size */
	{
		static unsigned int const sizes[4] = { 125, 250, 500, 1000 };
		double tt[2][4], ymin = 1e9, ymax = -1e9, xy[8];
		unsigned int j;
		for (k = 0; k < 4; k++)
			for (j = 0; j < 2; j++)
			{
				/* (dense beyond 500 edges takes tens of seconds a step: not drawn) */
				tt[j][k] = j == 0 && sizes[k] > 500 ? 0 : log10(s35_step(sizes[k], j ? 0.25 : 0));
				if (j || sizes[k] <= 500) { ymin = fmin(ymin, tt[j][k]); ymax = fmax(ymax, tt[j][k]); }
			}
		v.x0 = 24 + 2 * (pw + 12); v.xmin = log10(100.0); v.xmax = log10(1300.0); v.ymin = floor(ymin) - 0.2; v.ymax = ymax + 0.3;
		svg_panel(&s, &v, "seconds per step against edges");
		for (j = 0; j < 2; j++)
		{
			unsigned int n = j ? 4 : 3;
			char const* col = j ? "#cf222e" : "#0969da";
			for (k = 0; k < n; k++) { xy[2 * k] = vx(&v, log10((double)sizes[k])); xy[2 * k + 1] = vy(&v, tt[j][k]); svg_circle(&s, xy[2 * k], xy[2 * k + 1], 3, col, col); }
			svg_polyline(&s, xy, (int)n, col, 2, 1, 0);
		}
		svg_text(&s, v.x0 + 14, v.y0 + 44, 12, "#0969da", "start", "dense (metric and KKT by LU)");
		svg_text(&s, v.x0 + 14, v.y0 + 62, 12, "#cf222e", "start", "hierarchical, theta 0.25");
		for (k = 0; k < 4; k++)
		{
			sprintf(label, "%u", sizes[k]);
			svg_text(&s, vx(&v, log10((double)sizes[k])), v.y0 + ph + 16, 11, "#57606a", "middle", label);
		}
		for (k = 0; k < 8; k++)
		{
			double y = -3.0 + k;
			if (y < v.ymin || y > v.ymax) continue;
			sprintf(label, "1e%d", (int)y);
			svg_line(&s, v.x0, vy(&v, y), v.x0 + pw, vy(&v, y), "#d0d7de", 0.6, 1);
			svg_text(&s, v.x0 + pw - 6, vy(&v, y) - 3, 11, "#57606a", "end", label);
		}
		sprintf(label, "at 500 edges: %.2f s dense, %.2f s hierarchical", pow(10, tt[0][2]), pow(10, tt[1][2]));
		svg_text(&s, v.x0 + pw / 2, v.y0 + ph + 34, 12, "#57606a", "middle", label);
	}
	svg_close(&s);
	free(x); free(x0); free(e);
	printf("figure 35: hierarchical solve, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
