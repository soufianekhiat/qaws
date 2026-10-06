/* ================================================================== */
/*  8. Curve networks: Gordon "affects" and loft "affected by"         */
/* ================================================================== */

static qaws_scalar const g_net_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_curve* net_curve(qaws_scalar const* p)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = 6;
	d.knots = g_net_knots;
	d.knot_count = 10;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

static double net_height(double x, double y)
{
	return 0.35 * sin(1.1 * x) * cos(0.9 * y) + 0.15 * x;
}

typedef struct affects_ctx
{
	qaws_surface const* surface;
	qaws_diff_views const* views;
} affects_ctx;

static double affects_value(void const* user, qaws_scalar u, qaws_scalar v)
{
	affects_ctx const* a = (affects_ctx const*)user;
	qaws_surface_jet p, t;
	qaws_surface_eval_tangent(NULL, a->surface, u, v, 0, 0, QAWS_SJET_P, a->views, &p, &t);
	return sqrt(t.d[0].x * t.d[0].x + t.d[0].y * t.d[0].y + t.d[0].z * t.d[0].z);
}

static void draw_curve3(svg* s, projection const* pr, qaws_curve const* c, char const* color, double width)
{
	double xy[2 * 80];
	qaws_range r = qaws_curve_get_parameter_range(c);
	int k;
	for (k = 0; k < 80; k++)
	{
		qaws_eval_result_3d e;
		qaws_curve_evaluate_3d(c, r.min_value + (r.max_value - r.min_value) * k / (qaws_scalar)79, QAWS_EVAL_FLAG_POSITION, &e);
		project(pr, e.position.x, e.position.y, e.position.z, &xy[2 * k], &xy[2 * k + 1]);
	}
	svg_polyline(s, xy, 80, color, width, 1, 0);
}

static void demo_networks(void)
{
	static qaws_scalar const params[3] = { 0, 0.5f, 1 };
	qaws_scalar net[6][18];
	qaws_curve* curves[6];
	qaws_curve const* uc[3];
	qaws_curve const* vc[3];
	qaws_surface* gordon = NULL;
	qaws_surface_gordon_desc gd;
	svg s;
	int i, n;
	unsigned int sel_curve = 4, sel_cp = 3;  /* middle v-curve, interior control point */
	viewport a = { 30, 80, 560, 500, 0, 0, 0, 0 };
	viewport b = { 610, 80, 560, 500, 0, 0, 0, 0 };
	projection pa = { 310, 230, 95, 1.6 };

	for (i = 0; i < 6; i++)
	{
		for (n = 0; n < 6; n++)
		{
			double sv = n * 0.6, fixed = (i % 3) * 1.5;
			qaws_scalar* p = &net[i][3 * n];
			if (i < 3) { p[0] = (qaws_scalar)sv; p[1] = (qaws_scalar)fixed; }
			else { p[0] = (qaws_scalar)fixed; p[1] = (qaws_scalar)sv; }
			p[2] = (qaws_scalar)net_height(p[0], p[1]);
		}
		curves[i] = net_curve(net[i]);
	}
	for (i = 0; i < 3; i++)
	{
		uc[i] = curves[i];
		vc[i] = curves[3 + i];
	}
	gd.u_curves = uc;
	gd.u_curve_count = 3;
	gd.v_params = params;
	gd.v_curves = vc;
	gd.v_curve_count = 3;
	gd.u_params = params;
	qaws_surface_create_gordon(&gd, &gordon);

	svg_open(&s, "showcase/8_networks.svg", 1200, 600, "Curve networks: Gordon affects / loft affected by",
		"Left: tangent of one control point of the middle v-curve -> |dS| over a Gordon surface. Right: adjoint of one loft point -> influence per section.");

	/* Gordon: affects map through child views (child 4 = middle v-curve). */
	svg_panel(&s, &a, "Gordon: where does the red control point act?");
	{
		qaws_scalar dir[6][18];
		qaws_field_view fv[6];
		qaws_diff_views child[6], views;
		affects_ctx actx;
		double sx, sy;
		memset(dir, 0, sizeof(dir));
		dir[sel_curve][3 * sel_cp + 2] = 1;  /* move up */
		for (i = 0; i < 6; i++)
			child[i] = one_field(&fv[i], QAWS_FIELD_CONTROL_POINTS, dir[i], 6, 3);
		for (i = 0; i < 6; i++)
			child[i].fields = &fv[i];
		views.fields = NULL;
		views.field_count = 0;
		views.children = child;
		views.child_count = 6;
		actx.surface = gordon;
		actx.views = &views;
		draw_quads(&s, &pa, gordon, 26, 26, affects_value, &actx, 0, 1);
		for (i = 0; i < 6; i++)
			draw_curve3(&s, &pa, curves[i], i == (int)sel_curve ? "#cf222e" : "#24292f", i == (int)sel_curve ? 3 : 1.6);
		project(&pa, net[sel_curve][3 * sel_cp], net[sel_curve][3 * sel_cp + 1], net[sel_curve][3 * sel_cp + 2], &sx, &sy);
		svg_line(&s, sx, sy, sx, sy - 45, "#cf222e", 2.5, 1);
		svg_circle(&s, sx, sy, 6, "#cf222e", "#ffffff");
		svg_colorbar(&s, a.x0 + 12, a.y0 + a.h - 24, 200, 8, "|dS| = 0", "1");
		svg_text(&s, a.x0 + 12, a.y0 + a.h - 34, 12, "#57606a", "start",
			"the edit spreads along the v-curve and fades with the Catmull-Rom blend in u");
	}

	/* Loft: affected-by for one surface point through the section curves. */
	svg_panel(&s, &b, "Loft: which sections move the marked point?");
	{
		qaws_scalar sec[5][18], bars[5][18];
		qaws_curve* sections[5];
		qaws_curve const* sp[5];
		qaws_surface* loft = NULL;
		qaws_surface_loft_desc ld;
		qaws_field_view fv[5];
		qaws_diff_views child[5], views;
		qaws_surface_jet ybar, p, t;
		qaws_scalar u0 = (qaws_scalar)0.5, v0 = (qaws_scalar)0.42;
		double infl[5], imax = 0, sx, sy;
		projection pb = { 890, 180, 60, 1.6 };
		char buf[64], col[32];

		for (i = 0; i < 5; i++)
		{
			for (n = 0; n < 6; n++)
			{
				qaws_scalar* q = &sec[i][3 * n];
				q[0] = (qaws_scalar)(n * 0.6);
				q[1] = (qaws_scalar)(i * 1.0);
				q[2] = (qaws_scalar)(0.5 * sin(1.3 * n + i) * (0.4 + 0.15 * i));
			}
			sections[i] = net_curve(sec[i]);
			sp[i] = sections[i];
		}
		ld.sections = sp;
		ld.section_count = 5;
		ld.v_parameters = NULL;
		qaws_surface_create_loft(&ld, &loft);

		memset(bars, 0, sizeof(bars));
		for (i = 0; i < 5; i++)
			child[i] = one_field(&fv[i], QAWS_FIELD_CONTROL_POINTS, bars[i], 6, 3);
		for (i = 0; i < 5; i++)
			child[i].fields = &fv[i];
		views.fields = NULL;
		views.field_count = 0;
		views.children = child;
		views.child_count = 5;
		memset(&ybar, 0, sizeof(ybar));
		ybar.d[0] = v3(0, 0, 1);  /* height of the marked point */
		qaws_surface_eval_adjoint(NULL, loft, u0, v0, QAWS_SJET_P, &ybar, &views, NULL, NULL);
		for (i = 0; i < 5; i++)
		{
			infl[i] = 0;
			for (n = 0; n < 6; n++)
				infl[i] += fabs(bars[i][3 * n + 2]);
			if (infl[i] > imax) imax = infl[i];
		}

		draw_quads(&s, &pb, loft, 26, 26, vase_height_value, NULL, 0, 1);
		for (i = 0; i < 5; i++)
		{
			heat(infl[i] / imax, col);
			draw_curve3(&s, &pb, sections[i], col, 2 + 4 * infl[i] / imax);
		}
		qaws_surface_eval_tangent(NULL, loft, u0, v0, 0, 0, QAWS_SJET_P, NULL, &p, &t);
		project(&pb, p.d[0].x, p.d[0].y, p.d[0].z, &sx, &sy);
		svg_circle(&s, sx, sy, 7, "#cf222e", "#ffffff");

		/* bar chart of influence per section (signed sum shown as magnitude) */
		fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"175\" fill=\"#ffffff\" stroke=\"#d0d7de\"/>\n", b.x0 + 10, b.y0 + b.h - 185, b.w - 20);
		svg_text(&s, b.x0 + 20, b.y0 + b.h - 168, 12, "#57606a", "start", "sum over each section of |d height / d control point z|");
		for (i = 0; i < 5; i++)
		{
			double bx = b.x0 + 40 + i * 100, bh = 95 * infl[i] / imax;
			heat(infl[i] / imax, col);
			fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"60\" height=\"%.1f\" fill=\"%s\"/>\n", bx, b.y0 + b.h - 40 - bh, bh, col);
			sprintf(buf, "section %d", i);
			svg_text(&s, bx + 30, b.y0 + b.h - 24, 11, "#24292f", "middle", buf);
			sprintf(buf, "%.2f", infl[i]);
			svg_text(&s, bx + 30, b.y0 + b.h - 46 - bh, 11, "#24292f", "middle", buf);
		}
		qaws_surface_destroy(loft);
		for (i = 0; i < 5; i++)
			qaws_curve_destroy(sections[i]);
	}

	svg_close(&s);
	qaws_surface_destroy(gordon);
	for (i = 0; i < 6; i++)
		qaws_curve_destroy(curves[i]);
	printf("8_networks: written\n");
}

