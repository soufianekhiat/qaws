/*
 * Figure 37: repulsive curves on a mesh, distances along it
 * (qaws_repulsive_mesh_curves.h).
 *   - a strip folded over, a loop on each sheet across a small gap, each
 *     loop's edges held: measured by chords the loops push each other
 *     apart through the gap (the top one toward the fold); measured along
 *     the surface they are far apart, round the fold, and slide away from
 *     it, both the same way;
 *   - a star-shaped loop on a sphere, its length held, smoothed into a
 *     round one, its vertices walked across the triangles.
 */

#define MS_RM_GAP 0.15

/* a strip 1 wide folded over: z = 0 from x = 0 to 2, a half circle, back at z = MS_RM_GAP */
static void ms_rm_fold_at(double u, double w, double* p)
{
	double rr = 0.5 * MS_RM_GAP, L = 4 + MS_PI * rr, s = u * L;
	p[1] = w;
	if (s <= 2) { p[0] = s; p[2] = 0; }
	else if (s <= 2 + MS_PI * rr) { double a = (s - 2) / rr; p[0] = 2 + rr * sin(a); p[2] = rr - rr * cos(a); }
	else { p[0] = 2 - (s - 2 - MS_PI * rr); p[2] = MS_RM_GAP; }
}

static qaws_mesh_3d* ms_rm_fold(unsigned int nu, unsigned int nw)
{
	unsigned int i, j, k = 0, nv = (nu + 1) * (nw + 1);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nw);
	qaws_mesh_3d* m;
	for (j = 0; j <= nw; j++)
		for (i = 0; i <= nu; i++)
			ms_rm_fold_at((double)i / nu, (double)j / nw, x + 3 * (j * (nu + 1) + i));
	for (j = 0; j < nw; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int v = j * (nu + 1) + i;
			t[k++] = v; t[k++] = v + 1; t[k++] = v + nu + 2;
			t[k++] = v; t[k++] = v + nu + 2; t[k++] = v + nu + 1;
		}
	m = ms_make(x, nv, t, NULL, 2 * nu * nw);
	free(x);
	free(t);
	return m;
}

/* a unit sphere: rings between single poles on z */
static qaws_mesh_3d* ms_rm_sphere(unsigned int nu, unsigned int nv)
{
	unsigned int n = nu * (nv - 1) + 2, i, j, k = 0, nf = 2 * nu * (nv - 1);
	double* x = (double*)malloc(sizeof(double) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_mesh_3d* m;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = MS_PI * j / nv, ph = 2 * MS_PI * i / nu;
			double* p = x + 3 * ((j - 1) * nu + i);
			p[0] = sin(th) * cos(ph); p[1] = sin(th) * sin(ph); p[2] = cos(th);
		}
	x[3 * (n - 2)] = 0; x[3 * (n - 2) + 1] = 0; x[3 * (n - 2) + 2] = 1;
	x[3 * (n - 1)] = 0; x[3 * (n - 1) + 1] = 0; x[3 * (n - 1) + 2] = -1;
	for (i = 0; i < nu; i++)
	{
		unsigned int a = i, b = (i + 1) % nu, o = (nv - 2) * nu;
		t[k++] = n - 2; t[k++] = a; t[k++] = b;
		t[k++] = n - 1; t[k++] = o + b; t[k++] = o + a;
	}
	for (j = 0; j + 2 < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = b;
			t[k++] = a; t[k++] = a + nu; t[k++] = b + nu;
		}
	m = ms_make(x, n, t, NULL, k / 3);
	free(x);
	free(t);
	return m;
}

static qaws_repulsive_curves* ms_rm_loops(unsigned int n, unsigned int loops)
{
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * n * loops), l, i;
	qaws_repulsive_curves_desc d;
	qaws_repulsive_curves* c = NULL;
	for (l = 0; l < loops; l++)
		for (i = 0; i < n; i++)
		{
			e[2 * (l * n + i)] = l * n + i;
			e[2 * (l * n + i) + 1] = l * n + (i + 1) % n;
		}
	memset(&d, 0, sizeof d);
	d.dimension = 3;
	d.vertex_count = n * loops;
	d.edges = e;
	d.edge_count = n * loops;
	qaws_repulsive_curves_create(&c, &d);
	free(e);
	return c;
}

static void ms_rm_locate(qaws_mesh_location* out, qaws_mesh_3d const* m, qaws_scalar const* p, unsigned int n)
{
	qaws_mesh_reference_desc rd;
	qaws_mesh_reference* r = NULL;
	memset(&rd, 0, sizeof rd);
	rd.interpolant = QAWS_INTERPOLANT_LINEAR;
	qaws_mesh_reference_create(&r, m, &rd);
	qaws_mesh_reference_project(out, NULL, r, p, NULL, n);
	qaws_mesh_reference_destroy(r);
}

/* closed loops of n vertices each (3 per vertex), lifted along lift (NULL: none) by the scale s */
static void ms_rm_draw_loops(er_canvas* c, qaws_scalar const* x, unsigned int n, unsigned int loops, double const* const* col, double width, int radial)
{
	unsigned int l, i, k;
	for (l = 0; l < loops; l++)
		for (i = 0; i < n; i++)
		{
			double a[3], b[3], sa = 1, sb = 1;
			qaws_scalar const *p = x + 3 * (l * n + i), *q = x + 3 * (l * n + (i + 1) % n);
			if (radial)
			{
				sa = 1.012 / sqrt((double)p[0] * p[0] + (double)p[1] * p[1] + (double)p[2] * p[2]);
				sb = 1.012 / sqrt((double)q[0] * q[0] + (double)q[1] * q[1] + (double)q[2] * q[2]);
			}
			for (k = 0; k < 3; k++) { a[k] = sa * (double)p[k]; b[k] = sb * (double)q[k]; }
			er_line(c, a, b, col[l], width);
		}
}

static void demo_repulsive_mesh_curves(void)
{
	static double const wire[3] = { 150, 156, 170 }, grey[3] = { 120, 120, 120 }, blue[3] = { 40, 90, 200 }, orange[3] = { 230, 110, 30 },
		col[3] = { 200, 206, 220 }, dark[3] = { 60, 66, 80 };
	double const* two[2] = { blue, orange };
	double const* greys[2] = { grey, grey };
	er_canvas* c = er_create(1800, 720, 3, 252, 252, 250);
	char sub[200];
	int pass;
	if (!c)
		return;
	for (pass = 0; pass < 2; pass++)
	{
		/* the fold: a loop on each sheet, offset, each loop's edges held */
		qaws_mesh_3d* m = ms_rm_fold(56, 14);
		qaws_repulsive_curves* cu = ms_rm_loops(14, 2);
		qaws_mesh_location L[28];
		qaws_scalar p[84], x[84];
		qaws_repulsive_mesh_curves_desc md;
		qaws_repulsive_curves_solve_desc sd;
		qaws_repulsive_mesh_curves_solver* s = NULL;
		qaws_repulsive_report r;
		unsigned int i;
		double cx[2] = { 0, 0 };
		for (i = 0; i < 28; i++)
		{
			double t = 2 * MS_PI * i / 14;
			int top = i >= 14;
			p[3 * i] = (qaws_scalar)((top ? 1.1 : 0.9) + 0.22 * cos(t)); p[3 * i + 1] = (qaws_scalar)(0.5 + 0.22 * sin(t)); p[3 * i + 2] = (qaws_scalar)(top ? MS_RM_GAP : 0);
		}
		ms_rm_locate(L, m, p, 28);
		memset(&md, 0, sizeof md);
		md.mesh = m;
		md.euclidean = !pass;
		memset(&sd, 0, sizeof sd);
		sd.hold_edge_lengths = 1;
		sd.max_iterations = 60;
		memset(&r, 0, sizeof r);
		qaws_repulsive_mesh_curves_solve(&r, L, cu, &md, &sd);
		qaws_repulsive_mesh_curves_solver_create(&s, cu, L, &md, NULL);
		qaws_repulsive_mesh_curves_solver_get_positions(x, s);
		for (i = 0; i < 28; i++) cx[i >= 14] += (double)x[3 * i] / 14;
		sprintf(sub, "%u steps, centres x 0.90, 1.10 to %.2f, %.2f", r.iterations, cx[0], cx[1]);
		ms_panel(c, 600 * pass, 0, 600, 720, m, 0.55, 0.32, pass ? "along the surface" : "by chords", sub);
		ms_draw(c, m, NULL, NULL, 0, wire, 0.6);
		ms_rm_draw_loops(c, p, 14, 2, greys, 1.2, 0);
		ms_rm_draw_loops(c, x, 14, 2, two, 3.4, 0);
		printf("  fold, %s: %s\n", pass ? "geodesic" : "euclidean", sub);
		qaws_repulsive_mesh_curves_solver_destroy(s);
		qaws_repulsive_curves_destroy(cu);
		qaws_mesh_3d_destroy(m);
	}
	{
		/* a star on a sphere, its length held */
		qaws_mesh_3d* m = ms_rm_sphere(48, 24);
		qaws_repulsive_curves* cu = ms_rm_loops(48, 1);
		qaws_mesh_location L[48], L0[48];
		qaws_scalar p[144], x[144];
		qaws_repulsive_mesh_curves_desc md;
		qaws_repulsive_curves_solve_desc sd;
		qaws_repulsive_mesh_curves_solver* s = NULL;
		qaws_repulsive_report r;
		double const* one[1] = { orange };
		double const* g1[1] = { grey };
		unsigned int i, crossed = 0;
		for (i = 0; i < 48; i++)
		{
			double t = 2 * MS_PI * i / 48, th = 0.45 + 0.12 * cos(5 * t);
			p[3 * i] = (qaws_scalar)cos(th); p[3 * i + 1] = (qaws_scalar)(sin(th) * cos(t)); p[3 * i + 2] = (qaws_scalar)(sin(th) * sin(t));
		}
		ms_rm_locate(L, m, p, 48);
		memcpy(L0, L, sizeof L);
		memset(&md, 0, sizeof md);
		md.mesh = m;
		memset(&sd, 0, sizeof sd);
		sd.hold_total_length = 1;
		sd.max_iterations = 60;
		memset(&r, 0, sizeof r);
		qaws_repulsive_mesh_curves_solve(&r, L, cu, &md, &sd);
		qaws_repulsive_mesh_curves_solver_create(&s, cu, L, &md, NULL);
		qaws_repulsive_mesh_curves_solver_get_positions(x, s);
		for (i = 0; i < 48; i++) crossed += L[i].patch != L0[i].patch;
		sprintf(sub, "%u steps, energy %.0f to %.0f, %u of 48 walked", r.iterations, (double)r.initial_energy, (double)r.energy, crossed);
		ms_panel(c, 1200, 0, 600, 720, m, -1.2, 0.3, "a star on a sphere", sub);
		ms_draw(c, m, col, NULL, 1, dark, 0.25);
		ms_rm_draw_loops(c, p, 48, 1, g1, 2.0, 1);
		ms_rm_draw_loops(c, x, 48, 1, one, 3.4, 1);
		printf("  sphere: %s\n", sub);
		qaws_repulsive_mesh_curves_solver_destroy(s);
		qaws_repulsive_curves_destroy(cu);
		qaws_mesh_3d_destroy(m);
	}
	er_save(c, "showcase/mesh37_repulsive_mesh_curves.png");
	printf("mesh37_repulsive_mesh_curves.png\n");
	er_destroy(c);
}
