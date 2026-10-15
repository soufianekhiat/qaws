/* Figure 36: R5. Curves on a surface (qaws_repulsive_surface_curves_solve):
   five loops on a torus, their vertices' surface parameters the unknowns,
   the energy measured in 3D, the length held; they spread over the torus
   and cross its seams. Native control points (qaws_repulsive_controls_solve):
   two cubic B-splines, their end control points pinned, relaxed through
   their control points: the result still two B-splines. */

/* the torus (radii 1 and 0.4, axis z) as a light wireframe */
static void s36_torus(svg* s, projection const* pr)
{
	unsigned int i, j;
	for (i = 0; i < 36; i++)
		for (j = 0; j < 18; j++)
		{
			double a0 = 2 * 3.14159265358979323846 * i / 36, a1 = 2 * 3.14159265358979323846 * (i + 1) / 36;
			double b0 = 2 * 3.14159265358979323846 * j / 18, b1 = 2 * 3.14159265358979323846 * (j + 1) / 18;
			double p[2], q[2], r[2];
			project(pr, (1 + 0.4 * cos(b0)) * cos(a0), (1 + 0.4 * cos(b0)) * sin(a0), 0.4 * sin(b0), &p[0], &p[1]);
			project(pr, (1 + 0.4 * cos(b0)) * cos(a1), (1 + 0.4 * cos(b0)) * sin(a1), 0.4 * sin(b0), &q[0], &q[1]);
			project(pr, (1 + 0.4 * cos(b1)) * cos(a0), (1 + 0.4 * cos(b1)) * sin(a0), 0.4 * sin(b1), &r[0], &r[1]);
			svg_line(s, p[0], p[1], q[0], q[1], "#d0d7de", 0.6, 1);
			svg_line(s, p[0], p[1], r[0], r[1], "#d0d7de", 0.6, 1);
		}
}

static void s36_polyline3(svg* s, projection const* pr, qaws_scalar const* x, unsigned int const* e, unsigned int ne, char const* color, double width)
{
	unsigned int i;
	for (i = 0; i < ne; i++)
	{
		double a[2], b[2];
		project(pr, x[3 * e[2 * i]], x[3 * e[2 * i] + 1], x[3 * e[2 * i] + 2], &a[0], &a[1]);
		project(pr, x[3 * e[2 * i + 1]], x[3 * e[2 * i + 1] + 1], x[3 * e[2 * i + 1] + 2], &b[0], &b[1]);
		svg_line(s, a[0], a[1], b[0], b[1], color, width, 1);
	}
}

static void demo_repulsive_surface(void)
{
	enum { K = 5, N = 30 };
	clock_t t0 = clock();
	svg s;
	double const pw = 330, ph = 330;
	static char const* const colors[K] = { "#0969da", "#cf222e", "#1a7f37", "#8250df", "#bc4c00" };
	char label[200];
	viewport v;
	projection pr;
	unsigned int k, i;
	if (!svg_open(&s, "showcase/b2d36_repulsive_surface.svg", (int)(48 + 3 * (pw + 12)), (int)(110 + ph + 30),
		"Repulsive curves on a surface, and on control points",
		"qaws_repulsive_surface_curves_solve: the energy in 3D, the unknowns the surface parameters; qaws_repulsive_controls_solve"))
		return;
	/* loops on a torus */
	{
		qaws_curve* profile = NULL;
		qaws_surface* torus = NULL;
		qaws_arc_segment a;
		qaws_arc_desc ad;
		qaws_surface_revolution_desc rd;
		qaws_repulsive_surface_desc sd;
		qaws_repulsive_curves_desc cd;
		qaws_repulsive_curves* c = NULL;
		qaws_repulsive_curves_solve_desc d;
		qaws_repulsive_report r;
		qaws_scalar uv[2 * K * N], x[3 * K * N], x0[3 * K * N];
		unsigned int e[2 * K * N];
		qaws_range ru, rv;
		memset(&a, 0, sizeof(a));
		a.center[0] = 1; a.radius = (qaws_scalar)0.4; a.angle_end = (qaws_scalar)(2 * 3.14159265358979323846);
		memset(&ad, 0, sizeof(ad));
		ad.dimension = QAWS_DIMENSION_2D; ad.segments = &a; ad.segment_count = 1;
		qaws_curve_create_arc(&ad, &profile);
		memset(&rd, 0, sizeof(rd));
		rd.profile = profile; rd.axis_direction.z = 1;
		qaws_surface_create_revolution(&rd, &torus);
		ru = qaws_surface_get_u_range(torus); rv = qaws_surface_get_v_range(torus);
		/* five small loops crowded near one spot, around the tube */
		for (k = 0; k < K; k++)
			for (i = 0; i < N; i++)
			{
				double t = 2 * 3.14159265358979323846 * i / N;
				uv[2 * (k * N + i)] = (qaws_scalar)((double)ru.min_value + ((double)ru.max_value - ru.min_value) * (0.96 + 0.012 * k + 0.02 * sin(t)));
				uv[2 * (k * N + i) + 1] = (qaws_scalar)((double)rv.min_value + ((double)rv.max_value - rv.min_value) * ((double)i / N));
				e[2 * (k * N + i)] = k * N + i; e[2 * (k * N + i) + 1] = k * N + (i + 1) % N;
			}
		memset(&sd, 0, sizeof(sd));
		sd.surface = torus; sd.periodic_u = 1; sd.periodic_v = 1;
		qaws_repulsive_surface_positions(&sd, uv, K * N, x0);
		memset(&cd, 0, sizeof(cd));
		cd.dimension = 3; cd.vertex_count = K * N; cd.edges = e; cd.edge_count = K * N;
		qaws_repulsive_curves_create(&cd, &c);
		memset(&d, 0, sizeof(d));
		d.hold_edge_lengths = 1; d.max_iterations = 60;
		memset(&r, 0, sizeof(r));
		qaws_repulsive_surface_curves_solve(c, &sd, &d, uv, &r);
		qaws_repulsive_surface_positions(&sd, uv, K * N, x);
		for (i = 0; i < 2; i++)
		{
			v.x0 = 24 + i * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph; v.xmin = v.ymin = -1; v.xmax = v.ymax = 1;
			pr.cx = v.x0 + pw / 2; pr.cy = v.y0 + ph / 2 + 10; pr.scale = 82; pr.zscale = 1;
			sprintf(label, i ? "after %u steps: E %.1f, every edge's length held" : "five loops on a torus: E %.0f", i ? r.iterations : 0, i ? (double)r.energy : (double)r.initial_energy);
			if (!i) sprintf(label, "five loops on a torus, crowded: E %.0f", (double)r.initial_energy);
			svg_panel(&s, &v, label);
			s36_torus(&s, &pr);
			for (k = 0; k < K; k++) s36_polyline3(&s, &pr, i ? x : x0, e + 2 * k * N, N, colors[k], 2.2);
		}
		svg_text(&s, 24 + pw + 12 + 10, 84 + ph - 10, 12, "#57606a", "start", "the loops stay on the torus exactly; seams crossed");
		qaws_repulsive_curves_destroy(c);
		qaws_surface_destroy(torus);
		qaws_curve_destroy(profile);
	}
	/* two B-splines through their control points */
	{
		static qaws_scalar const knots[14] = { 0, 0, 0, 0, (qaws_scalar)(1.0 / 7), (qaws_scalar)(2.0 / 7), (qaws_scalar)(3.0 / 7), (qaws_scalar)(4.0 / 7),
			(qaws_scalar)(5.0 / 7), (qaws_scalar)(6.0 / 7), 1, 1, 1, 1 };
		qaws_scalar p[2][30];
		qaws_curve* cv[2] = { NULL, NULL };
		qaws_bspline_desc bd;
		qaws_repulsive_controls_desc cd;
		qaws_repulsive_controls* ctl = NULL;
		qaws_repulsive_curves_solve_desc d;
		qaws_repulsive_report r;
		qaws_scalar P[60], P0[60];
		unsigned int j;
		static unsigned int const ends[4] = { 0, 9, 10, 19 };
		for (k = 0; k < 2; k++)
			for (i = 0; i < 10; i++)
			{
				double t = (double)i / 9;
				p[k][3 * i] = (qaws_scalar)(-0.8 + 1.6 * t); p[k][3 * i + 1] = (qaws_scalar)((k ? 0.12 : -0.12) + 0.35 * sin(6 * t + k)); p[k][3 * i + 2] = (qaws_scalar)(0.3 * cos(5 * t + 2 * k));
			}
		for (k = 0; k < 2; k++)
		{
			memset(&bd, 0, sizeof(bd));
			bd.dimension = QAWS_DIMENSION_3D; bd.degree = 3; bd.control_points = p[k]; bd.control_point_count = 10; bd.knots = knots; bd.knot_count = 14;
			qaws_curve_create_bspline(&bd, &cv[k]);
		}
		memset(&cd, 0, sizeof(cd));
		cd.polyline.curves = (qaws_curve const* const*)cv; cd.polyline.curve_count = 2; cd.polyline.dimension = 3; cd.polyline.edge_length = (qaws_scalar)0.06;
		qaws_repulsive_controls_create(&cd, &ctl);
		qaws_repulsive_controls_read(ctl, P0);
		memcpy(P, P0, sizeof(P));
		memset(&d, 0, sizeof(d));
		d.hold_total_length = 1; d.max_iterations = 40; d.pinned = ends; d.pinned_count = 4;   /* (each end control pinned) */
		memset(&r, 0, sizeof(r));
		qaws_repulsive_controls_solve(ctl, &d, P, &r);
		v.x0 = 24 + 2 * (pw + 12); v.y0 = 84; v.w = pw; v.h = ph;
		pr.cx = v.x0 + pw / 2; pr.cy = v.y0 + ph / 2; pr.scale = 120; pr.zscale = 1;
		sprintf(label, "two B-splines by their control points: E %.1f -> %.1f", (double)r.initial_energy, (double)r.energy);
		svg_panel(&s, &v, label);
		for (j = 0; j < 2; j++)
		{
			qaws_scalar const* Q = j ? P : P0;
			for (k = 0; k < 2; k++)
			{
				qaws_curve* out = NULL;
				double xy[2 * 121];
				unsigned int m;
				qaws_repulsive_controls_to_curve(ctl, Q, k, &out);
				for (m = 0; out && m <= 120; m++)
				{
					qaws_eval_result_3d er;
					qaws_curve_evaluate_3d(out, (qaws_scalar)(m / 120.0), QAWS_EVAL_FLAG_POSITION, &er);
					project(&pr, er.position.x, er.position.y, er.position.z, &xy[2 * m], &xy[2 * m + 1]);
				}
				svg_polyline(&s, xy, 121, j ? colors[k] : "#afb8c1", j ? 2.2 : 1.4, 1, 0);
				for (m = 0; j && m < 10; m++)
				{
					double cx, cy;
					project(&pr, Q[30 * k + 3 * m], Q[30 * k + 3 * m + 1], Q[30 * k + 3 * m + 2], &cx, &cy);
					svg_circle(&s, cx, cy, 2.6, "#ffffff", colors[k]);
				}
				qaws_curve_destroy(out);
			}
		}
		svg_text(&s, v.x0 + 10, v.y0 + ph - 10, 12, "#57606a", "start", "grey: the start; dots: the control points found");
		qaws_repulsive_controls_destroy(ctl);
		qaws_curve_destroy(cv[0]); qaws_curve_destroy(cv[1]);
	}
	svg_close(&s);
	printf("figure 36: repulsive curves on a surface and control points, %.1f s\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
}
