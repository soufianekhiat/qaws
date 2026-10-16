/* Figure 32: fractional Sobolev descent of tangent-point energies
   (qaws_repulsive_curves_solve, R2). A lumpy trefoil with its total length
   held, energy-based (nothing stops a strand passing through another
   between the vertices: it ends unknotted) and collision-checked (every
   accepted motion validated: it stays a trefoil); open strands pinned at
   both ends, their edge lengths held, around a fixed obstacle loop (a
   pinned component); a Hopf link relaxed, collision-checked. Grey: the
   start. */

static void s32_draw(svg* s, viewport const* v, qaws_scalar const* x, unsigned int dim, unsigned int const* e, unsigned int ne,
	projection const* pr, char const* color, double width, double opacity)
{
	unsigned int i;
	for (i = 0; i < ne; i++)
	{
		double a[2], b[2];
		unsigned int p = e[2 * i], q = e[2 * i + 1];
		if (pr)
		{
			project(pr, x[3 * p], x[3 * p + 1], x[3 * p + 2], &a[0], &a[1]);
			project(pr, x[3 * q], x[3 * q + 1], x[3 * q + 2], &b[0], &b[1]);
		}
		else
		{
			a[0] = vx(v, x[2 * p]); a[1] = vy(v, x[2 * p + 1]); b[0] = vx(v, x[2 * q]); b[1] = vy(v, x[2 * q + 1]);
		}
		svg_line(s, a[0], a[1], b[0], b[1], color, width, opacity);
	}
}

/* the projection fitting two point sets (3D) into the panel */
static void s32_fit(projection* pr, viewport const* v, qaws_scalar const* a, qaws_scalar const* b, unsigned int n)
{
	double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 }, p[2];
	unsigned int i, k;
	projection u;
	u.cx = u.cy = 0; u.scale = 1; u.zscale = 1;
	for (k = 0; k < 2; k++)
		for (i = 0; i < n; i++)
		{
			qaws_scalar const* x = k ? b : a;
			project(&u, x[3 * i], x[3 * i + 1], x[3 * i + 2], &p[0], &p[1]);
			lo[0] = fmin(lo[0], p[0]); lo[1] = fmin(lo[1], p[1]); hi[0] = fmax(hi[0], p[0]); hi[1] = fmax(hi[1], p[1]);
		}
	pr->zscale = 1;
	pr->scale = 0.8 * fmin(v->w / (hi[0] - lo[0]), (v->h - 50) / (hi[1] - lo[1]));
	pr->cx = v->x0 + v->w / 2 - pr->scale * 0.5 * (lo[0] + hi[0]);
	pr->cy = v->y0 + v->h / 2 + 5 - pr->scale * 0.5 * (lo[1] + hi[1]);
}

static void demo_repulsive_solve(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 300, ph = 300, pi = 3.14159265358979323846;
	char label[200], sub[200];
	unsigned int k;
	if (!svg_open(&s, "showcase/b2d32_repulsive_solve.svg", (int)(48 + 4 * (pw + 12)), (int)(110 + ph + 30),
		"Fractional Sobolev descent of tangent-point energies",
		"qaws_repulsive_curves_solve: grey the start; energy-based steps may let strands pass through each other, collision-checked steps never do"))
		return;
	/* the trefoil, twice */
	for (k = 0; k < 2; k++)
	{
		enum { N = 48 };
		qaws_scalar x[3 * N], x0[3 * N];
		unsigned int e[2 * N], i;
		qaws_repulsive_curves_desc d;
		qaws_repulsive_curves* c = NULL;
		qaws_repulsive_curves_solve_desc sd;
		qaws_repulsive_report r;
		viewport v;
		projection pr;
		v.x0 = 24 + k * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1; v.xmax = v.ymax = 1;
		pr.cx = v.x0 + pw / 2; pr.cy = v.y0 + ph / 2 + 10; pr.scale = 38; pr.zscale = 1;
		for (i = 0; i < N; i++)
		{
			double t = 2 * pi * i / N;
			x[3 * i] = (qaws_scalar)(sin(t) + 2 * sin(2 * t) + 0.3 * cos(5 * t)); x[3 * i + 1] = (qaws_scalar)(cos(t) - 2 * cos(2 * t));
			x[3 * i + 2] = (qaws_scalar)(-0.6 * sin(3 * t));
			e[2 * i] = i; e[2 * i + 1] = (i + 1) % N;
		}
		memcpy(x0, x, sizeof(x));
		memset(&d, 0, sizeof(d));
		d.dimension = 3; d.vertex_count = N; d.edges = e; d.edge_count = N;
		memset(&sd, 0, sizeof(sd));
		sd.hold_total_length = 1; sd.collision_checked = (int)k; sd.max_iterations = 150;
		memset(&r, 0, sizeof(r));
		if (qaws_repulsive_curves_create(&c, &d) == QAWS_STATUS_OK)
			qaws_repulsive_curves_solve(&r, x, c, &sd);
		sprintf(label, k ? "48 vertices, collision-checked: still a trefoil" : "48 vertices, energy-based: passed through itself");
		s32_fit(&pr, &v, x0, x, N);
		svg_panel(&s, &v, label);
		s32_draw(&s, &v, x0, 3, e, N, &pr, "#8c959f", 1.2, 0.8);
		s32_draw(&s, &v, x, 3, e, N, &pr, k ? "#0969da" : "#cf222e", 2.4, 1);
		sprintf(sub, "E %.2f -> %.2f in %u steps (reason %d)", (double)r.initial_energy, (double)r.energy, r.iterations, (int)r.reason);
		svg_text(&s, v.x0 + 10, v.y0 + ph - 10, 12, "#57606a", "start", sub);
		qaws_repulsive_curves_destroy(c);
	}
	/* pinned strands round an obstacle (2D) */
	{
		enum { S = 5, M = 25, O = 24 };
		qaws_scalar x[2 * (S * M + O)], x0[2 * (S * M + O)];
		unsigned int e[2 * (S * (M - 1) + O)], pins[2 * S + O], i, j, ne = 0, np = 0;
		qaws_repulsive_curves_desc d;
		qaws_repulsive_curves* c = NULL;
		qaws_repulsive_curves_solve_desc sd;
		qaws_repulsive_report r;
		viewport v;
		v.x0 = 24 + 2 * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = -1.4; v.xmax = 1.4; v.ymin = -1.15; v.ymax = 1.65;
		for (j = 0; j < S; j++)
			for (i = 0; i < M; i++)
			{
				double t = (double)i / (M - 1), y0 = -0.8 + 0.4 * j;
				x[2 * (M * j + i)] = (qaws_scalar)(-1 + 2 * t);
				x[2 * (M * j + i) + 1] = (qaws_scalar)(y0 + 0.08 * sin(2 * pi * t * (j + 1)));
				if (i + 1 < M) { e[2 * ne] = M * j + i; e[2 * ne + 1] = M * j + i + 1; ne++; }
				if (i == 0 || i + 1 == M) pins[np++] = M * j + i;
			}
		for (i = 0; i < O; i++)
		{
			double t = 2 * pi * i / O;
			x[2 * (S * M + i)] = (qaws_scalar)(0.18 * cos(t)); x[2 * (S * M + i) + 1] = (qaws_scalar)(0.2 + 0.12 * sin(t));
			e[2 * ne] = S * M + i; e[2 * ne + 1] = S * M + (i + 1) % O; ne++;
			pins[np++] = S * M + i;
		}
		memcpy(x0, x, sizeof(x));
		memset(&d, 0, sizeof(d));
		d.dimension = 2; d.vertex_count = S * M + O; d.edges = e; d.edge_count = ne;
		memset(&sd, 0, sizeof(sd));
		sd.pinned = pins; sd.pinned_count = np; sd.hold_edge_lengths = 1; sd.collision_checked = 1; sd.max_iterations = 120;
		memset(&r, 0, sizeof(r));
		if (qaws_repulsive_curves_create(&c, &d) == QAWS_STATUS_OK)
			qaws_repulsive_curves_solve(&r, x, c, &sd);
		svg_panel(&s, &v, "strands pinned at both ends, round a fixed loop");
		s32_draw(&s, &v, x0, 2, e, ne, NULL, "#8c959f", 1.2, 0.8);
		s32_draw(&s, &v, x, 2, e, ne - O, NULL, "#0969da", 2.2, 1);
		s32_draw(&s, &v, x + 0, 2, e + 2 * (ne - O), O, NULL, "#1b1f24", 2.2, 1);
		for (i = 0; i < 2 * S; i++) svg_circle(&s, vx(&v, x[2 * pins[i]]), vy(&v, x[2 * pins[i] + 1]), 3, "#1b1f24", "none");
		sprintf(sub, "E %.1f -> %.1f in %u steps (reason %d)", (double)r.initial_energy, (double)r.energy, r.iterations, (int)r.reason);
		svg_text(&s, v.x0 + 10, v.y0 + ph - 10, 12, "#57606a", "start", sub);
		qaws_repulsive_curves_destroy(c);
	}
	/* a Hopf link */
	{
		enum { M = 32 };
		qaws_scalar x[3 * 2 * M], x0[3 * 2 * M];
		unsigned int e[2 * 2 * M], i;
		qaws_repulsive_curves_desc d;
		qaws_repulsive_curves* c = NULL;
		qaws_repulsive_curves_solve_desc sd;
		qaws_repulsive_report r;
		viewport v;
		projection pr;
		v.x0 = 24 + 3 * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1; v.xmax = v.ymax = 1;
		pr.cx = v.x0 + pw / 2 - 10; pr.cy = v.y0 + ph / 2 + 10; pr.scale = 70; pr.zscale = 1;
		for (i = 0; i < M; i++)
		{
			double t = 2 * pi * i / M;
			x[3 * i] = (qaws_scalar)cos(t); x[3 * i + 1] = (qaws_scalar)sin(t); x[3 * i + 2] = 0;
			x[3 * (M + i)] = (qaws_scalar)(0.9 + 0.8 * cos(t)); x[3 * (M + i) + 1] = 0; x[3 * (M + i) + 2] = (qaws_scalar)(0.45 * sin(t));
			e[2 * i] = i; e[2 * i + 1] = (i + 1) % M;
			e[2 * (M + i)] = M + i; e[2 * (M + i) + 1] = M + (i + 1) % M;
		}
		memcpy(x0, x, sizeof(x));
		memset(&d, 0, sizeof(d));
		d.dimension = 3; d.vertex_count = 2 * M; d.edges = e; d.edge_count = 2 * M;
		memset(&sd, 0, sizeof(sd));
		sd.hold_edge_lengths = 1; sd.collision_checked = 1; sd.max_iterations = 60;
		memset(&r, 0, sizeof(r));
		if (qaws_repulsive_curves_create(&c, &d) == QAWS_STATUS_OK)
			qaws_repulsive_curves_solve(&r, x, c, &sd);
		s32_fit(&pr, &v, x0, x, 2 * M);
		svg_panel(&s, &v, "a Hopf link, edge lengths held, collision-checked");
		s32_draw(&s, &v, x0, 3, e, 2 * M, &pr, "#8c959f", 1.2, 0.8);
		s32_draw(&s, &v, x, 3, e, M, &pr, "#0969da", 2.4, 1);
		s32_draw(&s, &v, x, 3, e + 2 * M, M, &pr, "#cf222e", 2.4, 1);
		sprintf(sub, "E %.2f -> %.2f in %u steps (reason %d)", (double)r.initial_energy, (double)r.energy, r.iterations, (int)r.reason);
		svg_text(&s, v.x0 + 10, v.y0 + ph - 10, 12, "#57606a", "start", sub);
		qaws_repulsive_curves_destroy(c);
	}
	svg_close(&s);
	printf("figure 32: fractional Sobolev descent, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
