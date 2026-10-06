/* ================================================================== */
/*  7. Coons patch fairing by editing two boundaries                  */
/* ================================================================== */

#define CB_CP 6
#define CB_GRID 9
#define CB_ITERS 400

static qaws_scalar const g_cb_knots[10] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };

static qaws_curve* cb_curve(qaws_scalar const* p)
{
	qaws_bspline_desc d;
	qaws_curve* c = NULL;
	memset(&d, 0, sizeof(d));
	d.dimension = QAWS_DIMENSION_3D;
	d.degree = 3;
	d.control_points = p;
	d.control_point_count = CB_CP;
	d.knots = g_cb_knots;
	d.knot_count = 10;
	qaws_curve_create_bspline(&d, &c);
	return c;
}

typedef struct cb_patch
{
	qaws_curve* curves[4];
	qaws_surface* surface;
} cb_patch;

static void cb_build(qaws_scalar (*p)[CB_CP * 3], cb_patch* out)
{
	qaws_surface_coons_desc d;
	int i;
	for (i = 0; i < 4; i++)
		out->curves[i] = cb_curve(p[i]);
	memset(&d, 0, sizeof(d));
	d.c0 = out->curves[0];
	d.c1 = out->curves[1];
	d.d0 = out->curves[2];
	d.d1 = out->curves[3];
	out->surface = NULL;
	qaws_surface_create_coons(&d, &out->surface);
}

static void cb_destroy(cb_patch* p)
{
	int i;
	qaws_surface_destroy(p->surface);
	for (i = 0; i < 4; i++)
		qaws_curve_destroy(p->curves[i]);
}

static double cb_abs_mean(void const* user, qaws_scalar u, qaws_scalar v)
{
	qaws_surface_jet j;
	qaws_surface_geometry g;
	qaws_surface_eval_jet((qaws_surface const*)user, u, v, QAWS_SJET_ORDER2, &j);
	qaws_surface_geometry_eval(&j, NULL, NULL, &g, NULL, NULL, NULL);
	return fabs(g.mean);
}

/* Mean-curvature energy on a grid and its adjoint into boundary views. */
static double cb_energy(cb_patch const* patch, qaws_diff_views* views)
{
	double e = 0;
	int i, j;
	for (i = 0; i < CB_GRID; i++)
		for (j = 0; j < CB_GRID; j++)
		{
			qaws_scalar u = (qaws_scalar)((i + 0.5) / CB_GRID), v = (qaws_scalar)((j + 0.5) / CB_GRID);
			qaws_surface_jet p, ybar;
			qaws_surface_geometry g, gbar;
			qaws_surface_eval_jet(patch->surface, u, v, QAWS_SJET_ORDER2, &p);
			qaws_surface_geometry_eval(&p, NULL, NULL, &g, NULL, NULL, NULL);
			e += g.mean * g.mean / (CB_GRID * CB_GRID);
			if (views)
			{
				memset(&gbar, 0, sizeof(gbar));
				gbar.mean = (qaws_scalar)(2 * g.mean / (CB_GRID * CB_GRID));
				memset(&ybar, 0, sizeof(ybar));
				qaws_surface_geometry_adjoint(&p, &gbar, &ybar, NULL);
				qaws_surface_eval_adjoint(NULL, patch->surface, u, v, ybar.channels, &ybar, views, NULL, NULL);
			}
		}
	return e;
}

static void cb_draw_boundaries(svg* s, projection const* pr, cb_patch const* patch)
{
	int i, k;
	for (i = 0; i < 4; i++)
	{
		double xy[2 * 80];
		qaws_range r = qaws_curve_get_parameter_range(patch->curves[i]);
		for (k = 0; k < 80; k++)
		{
			qaws_eval_result_3d e;
			qaws_curve_evaluate_3d(patch->curves[i], r.min_value + (r.max_value - r.min_value) * k / (qaws_scalar)79,
				QAWS_EVAL_FLAG_POSITION, &e);
			project(pr, e.position.x, e.position.y, e.position.z, &xy[2 * k], &xy[2 * k + 1]);
		}
		svg_polyline(s, xy, 80, (i == 1 || i == 3) ? "#cf222e" : "#24292f", (i == 1 || i == 3) ? 3 : 2, 1, 0);
	}
}

static void demo_coons(void)
{
	qaws_scalar p[4][CB_CP * 3], initial[4][CB_CP * 3];
	qaws_scalar grad[4][CB_CP * 3];
	unsigned char interior[CB_CP] = { 0, 1, 1, 1, 1, 0 };
	double energy[CB_ITERS + 1];
	adam opt;
	int i, n, it;
	svg s;
	char buf[200];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < 4; i++)
		for (n = 0; n < CB_CP; n++)
		{
			qaws_scalar t = (qaws_scalar)n * (qaws_scalar)0.6;
			qaws_scalar* q = &p[i][3 * n];
			qaws_scalar amp = (i == 1 || i == 3) ? (qaws_scalar)0.7 : (qaws_scalar)0.12;
			qaws_scalar wave = (n == 0 || n == CB_CP - 1) ? (qaws_scalar)0 : (qaws_scalar)(amp * sin(2.3 * n + i));
			if (i == 0) { q[0] = t; q[1] = 0; }
			else if (i == 1) { q[0] = t; q[1] = 3; }
			else if (i == 2) { q[0] = 0; q[1] = t; }
			else { q[0] = 3; q[1] = t; }
			q[2] = wave;
		}
	memcpy(initial, p, sizeof(p));

	for (it = 0; it <= CB_ITERS; it++)
	{
		cb_patch patch;
		qaws_field_view fv[4];
		qaws_diff_views child[4], views;
		cb_build(p, &patch);
		memset(grad, 0, sizeof(grad));
		for (i = 0; i < 4; i++)
		{
			child[i] = one_field(&fv[i], QAWS_FIELD_CONTROL_POINTS, grad[i], CB_CP, 3);
			child[i].fields = &fv[i];
			fv[i].active = interior;          /* end points stay on the corners */
			fv[i].component_mask = 1u << 2;   /* only heights move */
		}
		views.fields = NULL;
		views.field_count = 0;
		views.children = child;
		views.child_count = 4;
		/* Only boundaries c1 and d1 (children 1 and 3) are editable. */
		child[0].field_count = 0;
		child[2].field_count = 0;
		energy[it] = cb_energy(&patch, &views);
		cb_destroy(&patch);
		if (it < CB_ITERS)
			adam_step(&opt, &p[0][0], &grad[0][0], 4 * CB_CP * 3, 0.01);
	}

	svg_open(&s, "showcase/7_coons.svg", 1200, 600, "Coons patch fairing by editing two boundaries",
		"Minimize mean(H^2) over the patch; adjoint: mean curvature -> geometry -> Coons rule -> child curves. Editable: red boundaries, interior heights only (masks).");
	{
		viewport a = { 30, 80, 560, 400, 0, 0, 0, 0 };
		viewport b = { 610, 80, 560, 400, 0, 0, 0, 0 };
		viewport lv = { 30, 495, 1140, 95, 0, 0, 0, 0 };
		projection pa = { 310, 250, 72, 1.6 }, pb = { 890, 250, 72, 1.6 };
		cb_patch before, after;
		cb_build(initial, &before);
		cb_build(p, &after);
		svg_panel(&s, &a, "before: |mean curvature|");
		draw_quads(&s, &pa, before.surface, 24, 24, cb_abs_mean, before.surface, 0, 0.8);
		cb_draw_boundaries(&s, &pa, &before);
		svg_panel(&s, &b, "after: two boundaries edited by the optimizer");
		draw_quads(&s, &pb, after.surface, 24, 24, cb_abs_mean, after.surface, 0, 0.8);
		cb_draw_boundaries(&s, &pb, &after);
		svg_colorbar(&s, b.x0 + b.w - 230, b.y0 + b.h - 26, 200, 8, "|H| = 0", "0.8");
		svg_loss_plot(&s, &lv, energy, CB_ITERS + 1, "#cf222e", "mean(H^2) (log)");
		cb_destroy(&before);
		cb_destroy(&after);
	}
	svg_close(&s);
	sprintf(buf, "%.4e -> %.4e", energy[0], energy[CB_ITERS]);
	printf("7_coons: mean(H^2) %s\n", buf);
}

