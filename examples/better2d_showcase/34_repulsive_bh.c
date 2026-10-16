/* Figure 34: Barnes-Hut tangent-point energies and gradients
   (qaws_repulsive_curves_bh, qaws_repulsive_mesh_bh, R4). The relative
   error against the direct evaluator as theta falls, for a trefoil and a
   lumpy sphere; the time against size for trefoils, direct against
   theta 0.25. Uses s33_sphere (figure 33). */

static double s34_err(qaws_scalar const* a, qaws_scalar const* b, unsigned int n)
{
	double e = 0, m = 0;
	unsigned int i;
	for (i = 0; i < n; i++) { e = fmax(e, fabs((double)a[i] - b[i])); m = fmax(m, fabs((double)b[i])); }
	return e / m;
}

static void s34_trefoil(qaws_scalar* x, unsigned int* e, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * 3.14159265358979323846 * i / n;
		x[3 * i] = (qaws_scalar)(sin(t) + 2 * sin(2 * t) + 0.13 * cos(5 * t)); x[3 * i + 1] = (qaws_scalar)(cos(t) - 2 * cos(2 * t) + 0.07 * sin(4 * t)); x[3 * i + 2] = (qaws_scalar)(-0.6 * sin(3 * t));
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % n;
	}
}

/* energy and gradient, direct (theta < 0) or Barnes-Hut */
static double s34_eval(qaws_repulsive_curves* c, qaws_repulsive_mesh* m, qaws_scalar const* x, double theta, qaws_scalar* g, unsigned int n)
{
	qaws_scalar E = 0;
	memset(g, 0, sizeof(qaws_scalar) * n);
	if (theta < 0) { if (c) qaws_repulsive_curves_gradient(g, &E, c, x, 1); else qaws_repulsive_mesh_gradient(g, &E, m, x, 1); }
	else if (c) qaws_repulsive_curves_bh(g, &E, NULL, c, x, (qaws_scalar)theta, 1);
	else qaws_repulsive_mesh_bh(g, &E, NULL, m, x, (qaws_scalar)theta, 1);
	return (double)E;
}

static void demo_repulsive_bh(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 470, ph = 330, bottom = -16;
	enum { NT = 12 };
	static char const* const colors[4] = { "#0969da", "#54aeff", "#cf222e", "#ff8182" };
	static char const* const names[4] = { "trefoil (300 edges): energy", "trefoil: gradient", "lumpy sphere (512 faces): energy", "lumpy sphere: gradient" };
	double err[4][NT], xy[2 * NT];
	char label[200];
	viewport v;
	unsigned int i, k;
	if (!svg_open(&s, "showcase/b2d34_repulsive_bh.svg", (int)(48 + 2 * (pw + 12)), (int)(110 + ph + 60),
		"Repulsive energies by Barnes-Hut",
		"qaws_repulsive_curves_bh / _mesh_bh: far nodes taken whole; theta 0 is the direct sum"))
		return;
	/* errors against theta */
	for (k = 0; k < 2; k++)
	{
		qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 300), *g0 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 300), *g1 = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * 300);
		qaws_repulsive_curves* c = NULL;
		qaws_repulsive_mesh* m = NULL;
		unsigned int n, j;
		double E0;
		if (k == 0)
		{
			unsigned int e[600];
			qaws_repulsive_curves_desc d;
			s34_trefoil(x, e, 300);
			memset(&d, 0, sizeof(d));
			d.dimension = 3; d.vertex_count = 300; d.edges = e; d.edge_count = 300;
			qaws_repulsive_curves_create(&c, &d);
			n = 900;
		}
		else
		{
			double xd[3 * 260];
			unsigned int tris[3 * 520], nv, nf;
			qaws_repulsive_mesh_desc d;
			s33_sphere(3, xd, tris, &nv, &nf);
			for (j = 0; j < nv; j++)
			{
				double r = 1 + 0.15 * sin(3 * xd[3 * j]) * cos(2 * xd[3 * j + 1]);
				x[3 * j] = (qaws_scalar)(xd[3 * j] * r); x[3 * j + 1] = (qaws_scalar)(xd[3 * j + 1] * r); x[3 * j + 2] = (qaws_scalar)(xd[3 * j + 2] * r);
			}
			memset(&d, 0, sizeof(d));
			d.vertex_count = nv; d.triangles = tris; d.triangle_count = nf;
			qaws_repulsive_mesh_create(&m, &d);
			n = 3 * nv;
		}
		E0 = s34_eval(c, m, x, -1, g0, n);
		for (j = 0; j < NT; j++)
		{
			double E1 = s34_eval(c, m, x, 0.05 * (j + 1), g1, n);
			err[2 * k][j] = fmax(log10(fabs(E1 - E0) / E0), bottom);
			err[2 * k + 1][j] = fmax(log10(s34_err(g1, g0, n)), bottom);
		}
		qaws_repulsive_curves_destroy(c); qaws_repulsive_mesh_destroy(m);
		free(x); free(g0); free(g1);
	}
	v.x0 = 24; v.y0 = 84; v.w = pw; v.h = ph; v.xmin = 0; v.xmax = 0.65; v.ymin = bottom - 0.5; v.ymax = 0.5;
	svg_panel(&s, &v, "relative error against the direct sum");
	for (i = 0; i < 4; i++)
	{
		for (k = 0; k < NT; k++) { xy[2 * k] = vx(&v, 0.05 * (k + 1)); xy[2 * k + 1] = vy(&v, err[i][k]); }
		svg_polyline(&s, xy, NT, colors[i], 2, 1, i & 1);
		svg_text(&s, v.x0 + pw - 14, v.y0 + 0.42 * ph + 18 * i, 12, colors[i], "end", names[i]);
	}
	for (k = 0; k <= 16; k += 4)
	{
		sprintf(label, "1e-%u", k);
		svg_line(&s, v.x0, vy(&v, -(double)k), v.x0 + pw, vy(&v, -(double)k), "#d0d7de", 0.6, 1);
		svg_text(&s, v.x0 + pw - 6, vy(&v, -(double)k) - 3, 11, "#57606a", "end", label);
	}
	for (k = 1; k <= 6; k++)
	{
		sprintf(label, "%.1f", 0.1 * k);
		svg_text(&s, vx(&v, 0.1 * k), v.y0 + ph + 16, 11, "#57606a", "middle", label);
	}
	svg_text(&s, v.x0 + pw / 2, v.y0 + ph + 34, 12, "#57606a", "middle", "theta (a node is far when its size is at most theta times its distance)");
	/* time against size, trefoils */
	{
		static unsigned int const sizes[5] = { 125, 250, 500, 1000, 2000 };
		double tt[2][5], ymin = 1e9, ymax = -1e9;
		unsigned int j;
		for (k = 0; k < 5; k++)
		{
			unsigned int nn = sizes[k];
			qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn), *g = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn);
			unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * nn);
			qaws_repulsive_curves* c = NULL;
			qaws_repulsive_curves_desc d;
			s34_trefoil(x, e, nn);
			memset(&d, 0, sizeof(d));
			d.dimension = 3; d.vertex_count = nn; d.edges = e; d.edge_count = nn;
			qaws_repulsive_curves_create(&c, &d);
			for (j = 0; j < 2; j++)
			{
				unsigned int rep = 0;
				clock_t a = clock();
				do { s34_eval(c, NULL, x, j ? 0.25 : -1, g, 3 * nn); rep++; } while (clock() - a < CLOCKS_PER_SEC / 5);
				tt[j][k] = log10((double)(clock() - a) / CLOCKS_PER_SEC / rep);
				ymin = fmin(ymin, tt[j][k]); ymax = fmax(ymax, tt[j][k]);
			}
			qaws_repulsive_curves_destroy(c);
			free(x); free(g); free(e);
		}
		v.x0 = 24 + pw + 12; v.xmin = log10(100.0); v.xmax = log10(2500.0); v.ymin = floor(ymin) - 0.3; v.ymax = ymax + 0.3;
		svg_panel(&s, &v, "seconds per energy and gradient, trefoils");
		for (j = 0; j < 2; j++)
		{
			for (k = 0; k < 5; k++) { xy[2 * k] = vx(&v, log10((double)sizes[k])); xy[2 * k + 1] = vy(&v, tt[j][k]); svg_circle(&s, xy[2 * k], xy[2 * k + 1], 3, colors[2 * j], colors[2 * j]); }
			svg_polyline(&s, xy, 5, colors[2 * j], 2, 1, 0);
		}
		svg_text(&s, v.x0 + 14, v.y0 + 44, 12, colors[0], "start", "direct (all pairs)");
		svg_text(&s, v.x0 + 14, v.y0 + 62, 12, colors[2], "start", "Barnes-Hut, theta 0.25");
		for (k = 0; k < 5; k++)
		{
			sprintf(label, "%u", sizes[k]);
			svg_text(&s, vx(&v, log10((double)sizes[k])), v.y0 + ph + 16, 11, "#57606a", "middle", label);
		}
		for (k = 0; k < 7; k++)
		{
			double y = -5.0 + k;
			if (y < v.ymin || y > v.ymax) continue;
			sprintf(label, "1e%d", (int)y);
			svg_line(&s, v.x0, vy(&v, y), v.x0 + pw, vy(&v, y), "#d0d7de", 0.6, 1);
			svg_text(&s, v.x0 + pw - 6, vy(&v, y) - 3, 11, "#57606a", "end", label);
		}
		sprintf(label, "edges (log scale); at 2000: %.2f s direct, %.2f s Barnes-Hut", pow(10, tt[0][4]), pow(10, tt[1][4]));
		svg_text(&s, v.x0 + pw / 2, v.y0 + ph + 34, 12, "#57606a", "middle", label);
	}
	svg_close(&s);
	printf("figure 34: Barnes-Hut, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
