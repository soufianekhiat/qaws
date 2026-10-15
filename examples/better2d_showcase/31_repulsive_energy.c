/* Figure 31: tangent-point energies and their gradients
   (qaws_repulsive_curves_gradient, R1: the energy layer only, no solver).
   Arrows: minus the coordinate gradient at each vertex (scaled per panel),
   the direction the energy falls fastest in the plain L2 sense; it is not
   the Sobolev descent direction a solver will use. An arch beside
   a closed loop (self and mutual repulsion); a network with a junction; a
   trefoil polygon in 3D, projected. */

/* one panel: polylines (by their vertex ranges) and minus-gradient arrows */
static void s31_panel(svg* s, viewport const* v, qaws_scalar const* x, unsigned int dim, unsigned int const* edges, unsigned int ne,
	qaws_scalar const* g, unsigned int nv, projection const* pr, char const* label)
{
	double gmax = 0, span = v->xmax - v->xmin;
	unsigned int i;
	char text[200];
	svg_panel(s, v, label);
	for (i = 0; i < nv; i++)
		gmax = fmax(gmax, sqrt((double)g[dim * i] * g[dim * i] + (double)g[dim * i + 1] * g[dim * i + 1] + (dim == 3 ? (double)g[dim * i + 2] * g[dim * i + 2] : 0)));
	for (i = 0; i < ne; i++)
	{
		double a[2], b[2];
		unsigned int p = edges[2 * i], q = edges[2 * i + 1];
		if (pr)
		{
			project(pr, x[3 * p], x[3 * p + 1], x[3 * p + 2], &a[0], &a[1]);
			project(pr, x[3 * q], x[3 * q + 1], x[3 * q + 2], &b[0], &b[1]);
		}
		else
		{
			a[0] = vx(v, x[2 * p]); a[1] = vy(v, x[2 * p + 1]); b[0] = vx(v, x[2 * q]); b[1] = vy(v, x[2 * q + 1]);
		}
		svg_line(s, a[0], a[1], b[0], b[1], "#0969da", 2.2, 1);
	}
	for (i = 0; i < nv && gmax > 0; i++)
	{
		double k = 0.09 * span / gmax, a[2], b[2];
		if (pr)
		{
			project(pr, x[3 * i], x[3 * i + 1], x[3 * i + 2], &a[0], &a[1]);
			project(pr, x[3 * i] - k * g[3 * i], x[3 * i + 1] - k * g[3 * i + 1], x[3 * i + 2] - k * g[3 * i + 2], &b[0], &b[1]);
		}
		else
		{
			a[0] = vx(v, x[2 * i]); a[1] = vy(v, x[2 * i + 1]);
			b[0] = vx(v, x[2 * i] - k * g[2 * i]); b[1] = vy(v, x[2 * i + 1] - k * g[2 * i + 1]);
		}
		svg_line(s, a[0], a[1], b[0], b[1], "#cf222e", 1.4, 0.9);
		svg_circle(s, a[0], a[1], 1.8, "#0969da", "none");
	}
	sprintf(text, "arrows: minus the gradient (longest %.3g)", gmax);
	svg_text(s, v->x0 + 10, v->y0 + v->h - 10, 12, "#57606a", "start", text);
}

static void demo_repulsive_energy(void)
{
	clock_t t0 = clock();
	svg s;
	double const pw = 380, ph = 380, pi = 3.14159265358979323846;
	qaws_scalar x[3 * 120], g[3 * 120];
	unsigned int e[2 * 120], n, ne, i;
	char label[200];
	if (!svg_open(&s, "showcase/b2d31_repulsive_energy.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + ph + 20),
		"Tangent-point energies and their gradients",
		"qaws_repulsive_curves (Repulsive Curves, Eq. 17): exact derivatives of the discrete energy; the arrows are the plain L2 gradient, which crowds at bends and kinks: the Sobolev metric of a solver (R2) spreads it"))
		return;
	/* a U nearly closing on itself, and a loop beside it (2D) */
	{
		qaws_repulsive_curves_desc d;
		qaws_repulsive_curves* c = NULL;
		qaws_scalar E = 0;
		viewport v = { 24, 84, 0, 0, -1.4, 1.4, -1.4, 1.4 };
		v.w = pw; v.h = ph;
		n = 0; ne = 0;
		for (i = 0; i < 40; i++)
		{
			double t = pi * (0.15 + 0.7 * i / 39.0);
			x[2 * n] = (qaws_scalar)(0.9 * cos(t + pi) * 0.35 - 0.3); x[2 * n + 1] = (qaws_scalar)(-1.1 * sin(t + pi) - 0.1);
			if (i) { e[2 * ne] = n - 1; e[2 * ne + 1] = n; ne++; }
			n++;
		}
		for (i = 0; i < 30; i++)
		{
			double t = 2 * pi * i / 30;
			x[2 * n] = (qaws_scalar)(0.55 + 0.35 * cos(t)); x[2 * n + 1] = (qaws_scalar)(0.15 + 0.5 * sin(t));
			e[2 * ne] = n; e[2 * ne + 1] = 40 + (i + 1) % 30; ne++;
			n++;
		}
		memset(&d, 0, sizeof(d));
		d.dimension = 2; d.vertex_count = n; d.edges = e; d.edge_count = ne;
		memset(g, 0, sizeof(qaws_scalar) * 2 * n);
		if (qaws_repulsive_curves_create(&d, &c) == QAWS_STATUS_OK)
			qaws_repulsive_curves_gradient(c, x, 1, g, &E);
		sprintf(label, "an arch beside a loop (2D): E = %.4g", (double)E);
		s31_panel(&s, &v, x, 2, e, ne, g, n, NULL, label);
		qaws_repulsive_curves_destroy(c);
	}
	/* a network: three branches from one junction, and a polyline across */
	{
		qaws_repulsive_curves_desc d;
		qaws_repulsive_curves* c = NULL;
		qaws_scalar E = 0;
		viewport v = { 0, 84, 0, 0, -1.4, 1.4, -1.4, 1.4 };
		unsigned int b, k;
		v.x0 = 24 + (pw + 12); v.w = pw; v.h = ph;
		n = 1; ne = 0;
		x[0] = 0; x[1] = (qaws_scalar)0.2;
		for (b = 0; b < 3; b++)
		{
			double a = pi / 2 + 2 * pi * b / 3;
			for (k = 1; k <= 12; k++)
			{
				double r = 0.09 * k;
				x[2 * n] = (qaws_scalar)(r * cos(a + 0.25 * sin(0.4 * k))); x[2 * n + 1] = (qaws_scalar)(0.2 + r * sin(a + 0.25 * sin(0.4 * k)));
				e[2 * ne] = k == 1 ? 0 : n - 1; e[2 * ne + 1] = n; ne++;
				n++;
			}
		}
		for (k = 0; k < 24; k++)
		{
			x[2 * n] = (qaws_scalar)(-1.2 + 0.1 * k); x[2 * n + 1] = (qaws_scalar)(-0.75 + 0.12 * sin(0.5 * k));
			if (k) { e[2 * ne] = n - 1; e[2 * ne + 1] = n; ne++; }
			n++;
		}
		memset(&d, 0, sizeof(d));
		d.dimension = 2; d.vertex_count = n; d.edges = e; d.edge_count = ne;
		memset(g, 0, sizeof(qaws_scalar) * 2 * n);
		if (qaws_repulsive_curves_create(&d, &c) == QAWS_STATUS_OK)
			qaws_repulsive_curves_gradient(c, x, 1, g, &E);
		sprintf(label, "a junction of three branches and a polyline: E = %.4g", (double)E);
		s31_panel(&s, &v, x, 2, e, ne, g, n, NULL, label);
		qaws_repulsive_curves_destroy(c);
	}
	/* a trefoil polygon in 3D */
	{
		qaws_repulsive_curves_desc d;
		qaws_repulsive_curves* c = NULL;
		qaws_scalar E = 0;
		viewport v = { 0, 84, 0, 0, -1.4, 1.4, -1.4, 1.4 };
		projection pr;
		v.x0 = 24 + 2 * (pw + 12); v.w = pw; v.h = ph;
		pr.cx = v.x0 + pw / 2; pr.cy = v.y0 + ph / 2; pr.scale = 52; pr.zscale = 1;
		n = 60;
		for (i = 0; i < n; i++)
		{
			double t = 2 * pi * i / n;
			x[3 * i] = (qaws_scalar)(sin(t) + 2 * sin(2 * t)); x[3 * i + 1] = (qaws_scalar)(cos(t) - 2 * cos(2 * t)); x[3 * i + 2] = (qaws_scalar)(-sin(3 * t));
			e[2 * i] = i; e[2 * i + 1] = (i + 1) % n;
		}
		ne = n;
		memset(&d, 0, sizeof(d));
		d.dimension = 3; d.vertex_count = n; d.edges = e; d.edge_count = ne;
		memset(g, 0, sizeof(qaws_scalar) * 3 * n);
		if (qaws_repulsive_curves_create(&d, &c) == QAWS_STATUS_OK)
			qaws_repulsive_curves_gradient(c, x, 1, g, &E);
		sprintf(label, "a trefoil polygon (3D, projected): E = %.4g", (double)E);
		v.xmin = -4; v.xmax = 4;
		s31_panel(&s, &v, x, 3, e, ne, g, n, &pr, label);
		qaws_repulsive_curves_destroy(c);
	}
	svg_close(&s);
	printf("figure 31: tangent-point energies and gradients, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
