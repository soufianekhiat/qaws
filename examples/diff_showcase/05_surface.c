/* ================================================================== */
/*  5. Surface: masked height-field fit + curvature sensitivity       */
/* ================================================================== */

#define SN 6
#define SG 18
#define S_ITERS 400

static double target_height(double x, double y)
{
	double r2 = (x - 0.45) * (x - 0.45) + (y - 0.55) * (y - 0.55);
	return 0.55 * exp(-r2 / 0.08) - 0.25 * exp(-((x - 0.8) * (x - 0.8) + (y - 0.2) * (y - 0.2)) / 0.05);
}

static qaws_surface* bspline_surface(qaws_scalar const* cps)
{
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = SN;
	d.v_point_count = SN;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

/* Isometric projection of a world point into a viewport. */
static void iso(viewport const* v, double x, double y, double z, double* sx, double* sy)
{
	double px = (x - y) * 0.866, py = (x + y) * 0.5 - z * 1.4;
	*sx = v->x0 + v->w / 2 + px * v->w * 0.5;
	*sy = v->y0 + v->h * 0.34 + py * v->h * 0.42;
}

static void draw_surface(svg* s, viewport const* v, qaws_surface const* surf, double zlo, double zhi)
{
	int i, j;
	char col[32];
	for (i = 0; i < 24; i++)
		for (j = 0; j < 24; j++)
		{
			qaws_scalar u0 = (qaws_scalar)(i / 24.0), u1 = (qaws_scalar)((i + 1) / 24.0);
			qaws_scalar v0 = (qaws_scalar)(j / 24.0), v1 = (qaws_scalar)((j + 1) / 24.0);
			qaws_surface_eval_result r[4];
			double x[4], y[4], zc = 0;
			int k;
			qaws_surface_evaluate(surf, u0, v0, QAWS_SURFACE_EVAL_POSITION, &r[0]);
			qaws_surface_evaluate(surf, u1, v0, QAWS_SURFACE_EVAL_POSITION, &r[1]);
			qaws_surface_evaluate(surf, u1, v1, QAWS_SURFACE_EVAL_POSITION, &r[2]);
			qaws_surface_evaluate(surf, u0, v1, QAWS_SURFACE_EVAL_POSITION, &r[3]);
			for (k = 0; k < 4; k++)
			{
				iso(v, r[k].position.x, r[k].position.y, r[k].position.z, &x[k], &y[k]);
				zc += r[k].position.z / 4;
			}
			heat((zc - zlo) / (zhi - zlo), col);
			fprintf(s->f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" stroke=\"#ffffff\" stroke-width=\"0.4\" stroke-opacity=\"0.6\"/>\n",
				x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], col);
		}
}

static void demo_surface(void)
{
	qaws_scalar cps[SN * SN * 3], grad[SN * SN * 3];
	qaws_scalar us[SG * SG], vs[SG * SG];
	double tz[SG * SG], loss[S_ITERS + 1];
	viewport a = { 30, 80, 370, 420, 0, 0, 0, 0 };
	viewport b = { 415, 80, 370, 420, 0, 0, 0, 0 };
	viewport c = { 800, 80, 370, 420, 0, 0, 0, 0 };
	viewport lv = { 30, 515, 755, 75, 0, 0, 0, 0 };
	adam opt;
	int i, j, it;
	svg s;
	char buf[200];

	memset(&opt, 0, sizeof(opt));
	for (i = 0; i < SN; i++)
		for (j = 0; j < SN; j++)
		{
			qaws_scalar* p = &cps[(i * SN + j) * 3];
			p[0] = (qaws_scalar)i / (SN - 1);
			p[1] = (qaws_scalar)j / (SN - 1);
			p[2] = 0;
		}
	for (i = 0; i < SG; i++)
		for (j = 0; j < SG; j++)
		{
			us[i * SG + j] = (qaws_scalar)(i / (SG - 1.0));
			vs[i * SG + j] = (qaws_scalar)(j / (SG - 1.0));
		}

	svg_open(&s, "showcase/5_surface.svg", 1200, 600, "Surfaces: masked height-field fit and curvature sensitivity",
		"Bicubic B-spline (6x6 net). Fit: only z is active (component mask), x/y stay on the grid. Right: d(Gaussian curvature at the red point)/d(control point z).");
	svg_panel(&s, &a, "initial (flat) + target samples");
	{
		qaws_surface* s0 = bspline_surface(cps);
		draw_surface(&s, &a, s0, -0.3, 0.6);
		qaws_surface_destroy(s0);
	}

	for (it = 0; it <= S_ITERS; it++)
	{
		qaws_surface* sf = bspline_surface(cps);
		static qaws_surface_jet primal[SG * SG], tangent[SG * SG], ybar[SG * SG];
		qaws_field_view fv;
		qaws_diff_views views;
		double l = 0;

		qaws_surface_eval_batch_tangent(NULL, sf, us, vs, NULL, NULL, SG * SG, QAWS_SJET_P, NULL, primal, tangent);
		memset(ybar, 0, sizeof(ybar));
		for (i = 0; i < SG * SG; i++)
		{
			double dz;
			if (it == 0)
				tz[i] = target_height(primal[i].d[0].x, primal[i].d[0].y);
			dz = primal[i].d[0].z - tz[i];
			l += dz * dz / (SG * SG);
			ybar[i].d[0].z = (qaws_scalar)(2 * dz / (SG * SG));
		}
		loss[it] = l;

		memset(grad, 0, sizeof(grad));
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, SN * SN, 3);
		fv.component_mask = 1u << 2; /* z only */
		qaws_surface_eval_batch_adjoint(NULL, sf, us, vs, SG * SG, QAWS_SJET_P, ybar, &views, NULL, NULL);
		qaws_surface_destroy(sf);
		if (it < S_ITERS)
			adam_step(&opt, cps, grad, SN * SN * 3, 0.01);
	}

	for (i = 0; i < SG * SG; i += 3)
	{
		double x, y;
		iso(&a, us[i], vs[i], tz[i], &x, &y);
		svg_circle(&s, x, y, 1.6, "#24292f", "none");
	}

	svg_panel(&s, &b, "after fitting (z only)");
	{
		qaws_surface* s1 = bspline_surface(cps);
		draw_surface(&s, &b, s1, -0.3, 0.6);
		qaws_surface_destroy(s1);
	}

	/* Gaussian curvature sensitivity: K at one point -> geometry adjoint ->
	   surface adjoint on the control net. */
	svg_panel(&s, &c, "dK / d(z_ij) at the red point");
	{
		qaws_surface* s1 = bspline_surface(cps);
		qaws_scalar u0 = (qaws_scalar)0.30, v0 = (qaws_scalar)0.70;
		qaws_surface_jet p, ybar;
		qaws_surface_geometry g, gbar;
		qaws_field_view fv;
		qaws_diff_views views;
		double gmax = 0, x, y;
		char col[32];

		qaws_surface_eval_jet(s1, u0, v0, QAWS_SJET_ORDER2, &p);
		qaws_surface_geometry_eval(&p, NULL, NULL, &g, NULL, NULL, NULL);
		memset(&gbar, 0, sizeof(gbar));
		gbar.gaussian = 1;
		memset(&ybar, 0, sizeof(ybar));
		qaws_surface_geometry_adjoint(&p, &gbar, &ybar, NULL);
		memset(grad, 0, sizeof(grad));
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, grad, SN * SN, 3);
		qaws_surface_eval_adjoint(NULL, s1, u0, v0, ybar.channels, &ybar, &views, NULL, NULL);
		for (i = 0; i < SN * SN; i++)
			if (fabs(grad[i * 3 + 2]) > gmax) gmax = fabs(grad[i * 3 + 2]);

		/* Top view of the control net: one cell per control point, colored
		   by dK/dz (diverging: blue lowers K, red raises it). */
		{
			double gx0 = c.x0 + 45, gy0 = c.y0 + 45, cell = (c.w - 90) / SN;
			for (i = 0; i < SN; i++)
				for (j = 0; j < SN; j++)
				{
					double sv = grad[(i * SN + j) * 3 + 2] / gmax;
					x = gx0 + i * cell;
					y = gy0 + (SN - 1 - j) * cell;
					heat(0.5 + 0.5 * sv, col);
					fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\" stroke=\"#ffffff\" stroke-width=\"2\"/>\n",
						x, y, cell, cell, sv == 0 ? "#eaeef2" : col);
					sprintf(buf, "%+.2f", sv);
					svg_text(&s, x + cell / 2, y + cell / 2 + 4, 11, "#1b1f24", "middle", sv == 0 ? "0" : buf);
				}
			x = gx0 + p.d[0].x * (SN - 1) * cell + cell / 2;
			y = gy0 + (1 - p.d[0].y) * (SN - 1) * cell + cell / 2;
			svg_circle(&s, x, y, 7, "#cf222e", "#ffffff");
			svg_text(&s, gx0, gy0 - 8, 11, "#57606a", "start", "top view of the 6x6 control net (normalized dK/dz)");
		}
		sprintf(buf, "K = %.3f, mean = %.3f", (double)g.gaussian, (double)g.mean);
		svg_text(&s, c.x0 + 12, c.y0 + c.h - 34, 12, "#24292f", "start", buf);
		svg_colorbar(&s, c.x0 + 12, c.y0 + c.h - 24, 200, 8, "lowers K", "raises K");
		qaws_surface_destroy(s1);
	}
	svg_loss_plot(&s, &lv, loss, S_ITERS + 1, "#cf222e", "height-fit MSE (log)");
	svg_close(&s);
	printf("5_surface: fit MSE %.4e -> %.4e\n", loss[0], loss[S_ITERS]);
}

