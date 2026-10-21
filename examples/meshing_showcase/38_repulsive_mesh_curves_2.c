/*
 * Figure 38: repulsive curves on a mesh, distances along it, on a terrain
 * of Gaussian hills and a pit (qaws_repulsive_mesh_curves.h).
 *   - a tree (six branches from a root, each forking in two), crumpled at
 *     the start, its edge lengths held: it spreads over the hills;
 *   - a flower-shaped loop, its edge lengths held: it opens into a smooth
 *     loop draped over them.
 * The start in thin dark lines, the result thick; the terrain coloured by
 * height.
 */

static double ms_rm2_height(double x, double y)
{
	static double const g[6][4] = {   /* centre x, y, height, width */
		{ -0.7, -0.5, 0.55, 0.38 }, { 0.65, -0.6, 0.45, 0.32 }, { 0.2, 0.75, 0.6, 0.42 },
		{ -0.8, 0.7, 0.35, 0.28 }, { 0.95, 0.45, 0.4, 0.3 }, { 0.05, -0.05, -0.3, 0.3 } };
	double z = 0;
	unsigned int k;
	for (k = 0; k < 6; k++)
	{
		double dx = x - g[k][0], dy = y - g[k][1];
		z += g[k][2] * exp(-(dx * dx + dy * dy) / (2 * g[k][3] * g[k][3]));
	}
	return z;
}

static qaws_mesh_3d* ms_rm2_terrain(unsigned int n, double half)
{
	unsigned int i, j, k = 0, nv = (n + 1) * (n + 1);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n);
	qaws_mesh_3d* m;
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			double* p = x + 3 * (j * (n + 1) + i);
			p[0] = -half + 2 * half * i / n; p[1] = -half + 2 * half * j / n; p[2] = ms_rm2_height(p[0], p[1]);
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int v = j * (n + 1) + i;
			/* (diagonals alternating, no grain) */
			if ((i + j) & 1) { t[k++] = v; t[k++] = v + 1; t[k++] = v + n + 2; t[k++] = v; t[k++] = v + n + 2; t[k++] = v + n + 1; }
			else { t[k++] = v; t[k++] = v + 1; t[k++] = v + n + 1; t[k++] = v + 1; t[k++] = v + n + 2; t[k++] = v + n + 1; }
		}
	m = ms_make(x, nv, t, NULL, 2 * n * n);
	free(x);
	free(t);
	return m;
}

/* edges drawn lifted a little above the terrain */
static void ms_rm2_draw_edges(er_canvas* c, qaws_scalar const* x, unsigned int const* e, unsigned int ne, double const* col, double width)
{
	unsigned int i, k;
	for (i = 0; i < ne; i++)
	{
		double a[3], b[3];
		for (k = 0; k < 3; k++) { a[k] = (double)x[3 * e[2 * i] + k]; b[k] = (double)x[3 * e[2 * i + 1] + k]; }
		a[2] += 0.02; b[2] += 0.02;
		er_line(c, a, b, col, width);
	}
}

/* solve from the points p (xy; z from the terrain), draw start and end in a panel */
static void ms_rm2_run(er_canvas* c, int px, qaws_mesh_3d const* m, double const* vrgb, qaws_scalar* p, unsigned int nv, unsigned int const* e,
	unsigned int ne, char const* title)
{
	static double const dark[3] = { 40, 40, 48 }, wire[3] = { 90, 96, 110 }, orange[3] = { 245, 245, 250 };
	qaws_repulsive_curves_desc cd;
	qaws_repulsive_curves* cu = NULL;
	qaws_mesh_location* L = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_repulsive_mesh_curves_desc md;
	qaws_repulsive_curves_solve_desc sd;
	qaws_repulsive_mesh_curves_solver* s = NULL;
	qaws_repulsive_report r;
	char sub[200];
	unsigned int i;
	for (i = 0; i < nv; i++) p[3 * i + 2] = (qaws_scalar)ms_rm2_height((double)p[3 * i], (double)p[3 * i + 1]);
	memset(&cd, 0, sizeof cd);
	cd.dimension = 3;
	cd.vertex_count = nv;
	cd.edges = e;
	cd.edge_count = ne;
	qaws_repulsive_curves_create(&cu, &cd);
	ms_rm_locate(L, m, p, nv);
	memset(&md, 0, sizeof md);
	md.mesh = m;
	memset(&sd, 0, sizeof sd);
	sd.hold_edge_lengths = 1;
	sd.max_iterations = 200;
	memset(&r, 0, sizeof r);
	qaws_repulsive_mesh_curves_solve(&r, L, cu, &md, &sd);
	qaws_repulsive_mesh_curves_solver_create(&s, cu, L, &md, NULL);
	qaws_repulsive_mesh_curves_solver_get_positions(x, s);
	sprintf(sub, "%u vertices, %u steps, energy %.0f to %.0f", nv, r.iterations, (double)r.initial_energy, (double)r.energy);
	ms_panel(c, px, 0, 900, 820, m, 0.5, 0.75, title, sub);
	ms_draw(c, m, NULL, vrgb, 1, wire, 0.2);
	ms_rm2_draw_edges(c, p, e, ne, dark, 1.6);
	ms_rm2_draw_edges(c, x, e, ne, dark, 6.0);
	ms_rm2_draw_edges(c, x, e, ne, orange, 3.4);
	printf("  %s: %s\n", title, sub);
	qaws_repulsive_mesh_curves_solver_destroy(s);
	qaws_repulsive_curves_destroy(cu);
	free(L);
	free(x);
}

static void demo_repulsive_mesh_curves_2(void)
{
	er_canvas* c = er_create(1800, 820, 3, 252, 252, 250);
	qaws_mesh_3d* m = ms_rm2_terrain(44, 1.6);
	unsigned int nmv = qaws_mesh_3d_get_vertex_count(m), i;
	double* vrgb = (double*)malloc(sizeof(double) * 3 * nmv);
	qaws_scalar* mp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nmv);
	if (!c)
		return;
	qaws_mesh_3d_get_positions(mp, m);
	for (i = 0; i < nmv; i++)
		ms_colormap(0.15 + 0.75 * ((double)mp[3 * i + 2] + 0.3) / 0.9, vrgb + 3 * i);
	{
		/* the tree: root 0; six branches of 8 vertices, each forking into two twigs of 5, zigzagging tight */
		enum { NB = 6, NA = 8, NT = 5, NV = 1 + NB * (NA + 2 * NT) };
		qaws_scalar p[3 * NV];
		unsigned int e[2 * (NV - 1)], v = 1, ne = 0, b, k, f;
		p[0] = 0; p[1] = 0.05f;
		for (b = 0; b < NB; b++)
		{
			double th = 2 * MS_PI * b / NB + 0.3, ux = cos(th), uy = sin(th), tipx = 0, tipy = 0;
			unsigned int prev = 0, tip;
			for (k = 1; k <= NA; k++, v++)
			{
				double rr = 0.035 * k, side = (k & 1 ? 1 : -1) * 0.03;
				p[3 * v] = (qaws_scalar)(p[0] + rr * ux - side * uy); p[3 * v + 1] = (qaws_scalar)(p[1] + rr * uy + side * ux);
				e[2 * ne] = prev; e[2 * ne + 1] = v; ne++;
				prev = v;
				tipx = (double)p[3 * v]; tipy = (double)p[3 * v + 1];
			}
			tip = prev;
			for (f = 0; f < 2; f++)
			{
				double ph = th + (f ? 0.6 : -0.6), wx = cos(ph), wy = sin(ph);
				prev = tip;
				for (k = 1; k <= NT; k++, v++)
				{
					double rr = 0.035 * k, side = (k & 1 ? 1 : -1) * 0.025;
					p[3 * v] = (qaws_scalar)(tipx + rr * wx - side * wy); p[3 * v + 1] = (qaws_scalar)(tipy + rr * wy + side * wx);
					e[2 * ne] = prev; e[2 * ne + 1] = v; ne++;
					prev = v;
				}
			}
		}
		ms_rm2_run(c, 0, m, vrgb, p, NV, e, ne, "a tree spreading over the hills");
	}
	{
		/* a flower: a closed loop of 50 vertices, seven petals about the pit */
		enum { N = 50 };
		qaws_scalar p[3 * N];
		unsigned int e[2 * N];
		for (i = 0; i < N; i++)
		{
			double t = 2 * MS_PI * i / N, rr = 0.32 + 0.16 * cos(7 * t);
			p[3 * i] = (qaws_scalar)(0.05 + rr * cos(t)); p[3 * i + 1] = (qaws_scalar)(-0.05 + rr * sin(t));
			e[2 * i] = i; e[2 * i + 1] = (i + 1) % N;
		}
		ms_rm2_run(c, 900, m, vrgb, p, N, e, N, "a flower opening over them");
	}
	er_save(c, "showcase/mesh38_repulsive_mesh_curves_2.png");
	printf("mesh38_repulsive_mesh_curves_2.png\n");
	free(vrgb);
	free(mp);
	qaws_mesh_3d_destroy(m);
	er_destroy(c);
}
