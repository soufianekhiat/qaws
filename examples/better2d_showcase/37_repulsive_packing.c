/* Figure 37: curves packed into a region (qaws_repulsive_curves_solve,
   collision-checked). A container is a pinned loop: the curves inside are
   repelled by it and, the motion checked, can never cross it. Left: one
   long strand folded inside a disc, every edge's length held; it spreads
   until it fills the disc evenly. Right: three strands in a square. */

static void s37_draw2(svg* s, viewport const* v, qaws_scalar const* x, unsigned int const* e, unsigned int ne, char const* color, double width)
{
	unsigned int i;
	for (i = 0; i < ne; i++)
		svg_line(s, vx(v, x[2 * e[2 * i]]), vy(v, x[2 * e[2 * i] + 1]), vx(v, x[2 * e[2 * i + 1]]), vy(v, x[2 * e[2 * i + 1] + 1]), color, width, 1);
}

/* solve: n vertices (2D), edges; the first nfix vertices pinned (the container) */
static void s37_solve(qaws_scalar* x, unsigned int n, unsigned int const* e, unsigned int ne, unsigned int nfix, unsigned int steps, qaws_repulsive_report* r)
{
	qaws_repulsive_curves_desc cd;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_curves_solve_desc d;
	unsigned int* pins = (unsigned int*)malloc(sizeof(unsigned int) * nfix), i;
	for (i = 0; i < nfix; i++) pins[i] = i;
	memset(&cd, 0, sizeof(cd));
	cd.dimension = 2; cd.vertex_count = n; cd.edges = e; cd.edge_count = ne;
	qaws_repulsive_curves_create(&cd, &c);
	memset(&d, 0, sizeof(d));
	d.pinned = pins; d.pinned_count = nfix; d.hold_edge_lengths = 1; d.collision_checked = 1; d.max_iterations = steps;
	memset(r, 0, sizeof(*r));
	qaws_repulsive_curves_solve(c, &d, x, r);
	qaws_repulsive_curves_destroy(c);
	free(pins);
}

static void demo_repulsive_packing(void)
{
	enum { NC = 48, NS = 140, SQ = 40, NT = 50 };
	clock_t t0 = clock();
	svg s;
	double const pw = 330, ph = 330;
	char label[200];
	viewport v;
	unsigned int i, k;
	if (!svg_open(&s, "showcase/b2d37_repulsive_packing.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + ph + 30),
		"Repulsive curves packed into a region",
		"collision-checked: a pinned loop contains them; every edge's length held; they spread to fill it"))
		return;
	/* a disc: the container (NC edges) and one folded strand (NS edges) */
	{
		qaws_scalar x[2 * (NC + NS + 1)], x0[2 * (NC + NS + 1)];
		unsigned int e[2 * (NC + NS)], n = NC + NS + 1;
		qaws_repulsive_report r;
		for (i = 0; i < NC; i++)
		{
			double t = 2 * 3.14159265358979323846 * i / NC;
			x[2 * i] = (qaws_scalar)cos(t); x[2 * i + 1] = (qaws_scalar)sin(t);
			e[2 * i] = i; e[2 * i + 1] = (i + 1) % NC;
		}
		/* a serpentine in the middle of the disc: rows 0.06 apart */
		for (i = 0; i <= NS; i++)
		{
			unsigned int row = i / 20, col = i % 20;
			double xx = -0.38 + 0.04 * (row % 2 ? 19 - col : col), yy = -0.21 + 0.06 * row;
			x[2 * (NC + i)] = (qaws_scalar)xx; x[2 * (NC + i) + 1] = (qaws_scalar)yy;
			if (i < NS) { e[2 * (NC + i)] = NC + i; e[2 * (NC + i) + 1] = NC + i + 1; }
		}
		memcpy(x0, x, sizeof(x));
		s37_solve(x, n, e, NC + NS, NC, 120, &r);
		for (k = 0; k < 2; k++)
		{
			v.x0 = 24 + k * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1.1; v.xmax = v.ymax = 1.1;
			if (k) sprintf(label, "after %u steps: E %.0f, never crossing the disc", r.iterations, (double)r.energy);
			else sprintf(label, "a folded strand (%u edges) in a disc: E %.0f", NS, (double)r.initial_energy);
			svg_panel(&s, &v, label);
			s37_draw2(&s, &v, k ? x : x0, e, NC, "#57606a", 2.4);
			s37_draw2(&s, &v, k ? x : x0, e + 2 * NC, NS, "#0969da", 1.8);
		}
	}
	/* a square: the container (4 SQ edges) and three strands (NT edges each) */
	{
		static char const* const colors[3] = { "#cf222e", "#1a7f37", "#8250df" };
		unsigned int nb = 4 * SQ, n = nb + 3 * (NT + 1), ne = nb + 3 * NT, j;
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
		qaws_repulsive_report r;
		for (i = 0; i < nb; i++)
		{
			unsigned int side = i / SQ;
			double t = (double)(i % SQ) / SQ, px = side == 0 ? -1 + 2 * t : side == 1 ? 1 : side == 2 ? 1 - 2 * t : -1, py = side == 0 ? -1 : side == 1 ? -1 + 2 * t : side == 2 ? 1 : 1 - 2 * t;
			x[2 * i] = (qaws_scalar)px; x[2 * i + 1] = (qaws_scalar)py;
			e[2 * i] = i; e[2 * i + 1] = (i + 1) % nb;
		}
		for (k = 0; k < 3; k++)
			for (j = 0; j <= NT; j++)
			{
				unsigned int vtx = nb + k * (NT + 1) + j;
				double t = (double)j / NT;
				/* a small serpentine each: rows of 10 edges, 0.05 apart */
				unsigned int row = j / 10, col = j % 10;
				x[2 * vtx] = (qaws_scalar)(-0.25 + 0.05 * (row % 2 ? 10 - col : col) + 0.0 * t); x[2 * vtx + 1] = (qaws_scalar)(-0.45 + 0.33 * k + 0.05 * row);
				if (j < NT) { e[2 * (nb + k * NT + j)] = vtx; e[2 * (nb + k * NT + j) + 1] = vtx + 1; }
			}
		s37_solve(x, n, e, ne, nb, 120, &r);
		v.x0 = 24 + 2 * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1.25; v.xmax = v.ymax = 1.15;
		sprintf(label, "three strands in a square: E %.0f -> %.0f", (double)r.initial_energy, (double)r.energy);
		svg_panel(&s, &v, label);
		s37_draw2(&s, &v, x, e, nb, "#57606a", 2.4);
		for (k = 0; k < 3; k++) s37_draw2(&s, &v, x, e + 2 * (nb + k * NT), NT, colors[k], 2);
		free(x); free(e);
	}
	svg_close(&s);
	printf("figure 37: repulsive packing, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
