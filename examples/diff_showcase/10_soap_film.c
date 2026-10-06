/* ================================================================== */
/*  10. Soap film: area minimization with gradient descent vs Newton  */
/*      steps built from Hessian-vector products                      */
/* ================================================================== */

#define SF_N 10
#define SF_CP (SF_N * SF_N)
#define SF_GD_ITERS 60
#define SF_NEWTON_ITERS 12

static qaws_surface* film_surface(qaws_scalar const* cps)
{
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = SF_N;
	d.v_point_count = SF_N;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

static void film_init(qaws_scalar* cps, unsigned char* active)
{
	int i, j;
	for (i = 0; i < SF_N; i++)
		for (j = 0; j < SF_N; j++)
		{
			double x = i / (double)(SF_N - 1), y = j / (double)(SF_N - 1);
			qaws_scalar* p = &cps[(i * SF_N + j) * 3];
			int boundary = i == 0 || j == 0 || i == SF_N - 1 || j == SF_N - 1;
			p[0] = (qaws_scalar)x;
			p[1] = (qaws_scalar)y;
			/* wavy frame; interior starts as a tilted bump far from minimal */
			p[2] = boundary ? (qaws_scalar)(0.3 * sin(2 * PI * x) * cos(PI * y) + 0.2 * sin(3 * PI * y) * (x - 0.5))
			                : (qaws_scalar)(0.45 * sin(PI * x) * sin(PI * y));
			active[i * SF_N + j] = (unsigned char)!boundary;
		}
}

static double film_area(qaws_scalar const* cps)
{
	qaws_surface* s = film_surface(cps);
	qaws_scalar a = 0;
	qaws_surface_functional_eval(NULL, s, QAWS_FUNCTIONAL_AREA, 0, NULL, &a, NULL, NULL);
	qaws_surface_destroy(s);
	return a;
}

/* Masked view: interior control points, z component only. */
static qaws_diff_views film_view(qaws_field_view* fv, qaws_scalar* data, unsigned char const* active)
{
	qaws_diff_views v = one_field(fv, QAWS_FIELD_CONTROL_POINTS, data, SF_CP, 3);
	fv->active = active;
	fv->component_mask = 1u << 2;
	return v;
}

static void film_gradient(qaws_scalar const* cps, unsigned char const* active, qaws_scalar* g)
{
	qaws_surface* s = film_surface(cps);
	qaws_field_view fv;
	qaws_diff_views v;
	memset(g, 0, sizeof(qaws_scalar) * SF_CP * 3);
	v = film_view(&fv, g, active);
	qaws_surface_functional_gradient(NULL, s, QAWS_FUNCTIONAL_AREA, 0, &v, NULL);
	qaws_surface_destroy(s);
}

static void film_hvp(qaws_surface const* s, unsigned char const* active, qaws_scalar* dir, qaws_scalar* out)
{
	qaws_field_view fd, fo;
	qaws_diff_views vd = film_view(&fd, dir, active), vo;
	memset(out, 0, sizeof(qaws_scalar) * SF_CP * 3);
	vo = film_view(&fo, out, active);
	qaws_surface_functional_hvp(NULL, s, QAWS_FUNCTIONAL_AREA, 0, &vd, &vo);
}

static double vdot(qaws_scalar const* a, qaws_scalar const* b, int n)
{
	double s = 0;
	int i;
	for (i = 0; i < n; i++)
		s += (double)a[i] * b[i];
	return s;
}

static double film_mean_abs(void const* user, qaws_scalar u, qaws_scalar v)
{
	qaws_surface_jet j;
	qaws_surface_geometry g;
	qaws_surface_eval_jet((qaws_surface const*)user, u, v, QAWS_SJET_ORDER2, &j);
	qaws_surface_geometry_eval(&j, NULL, NULL, &g, NULL, NULL, NULL);
	return fabs(g.mean);
}

static void demo_soap_film(void)
{
	enum { N3 = SF_CP * 3 };
	qaws_scalar init[N3], gd[N3], nt[N3], g[N3];
	unsigned char active[SF_CP];
	double e_gd[SF_GD_ITERS + 1], e_nt[SF_GD_ITERS + 1];
	int it, i, hvp_calls = 0;
	svg s;
	char buf[220];

	film_init(init, active);
	memcpy(gd, init, sizeof(init));
	memcpy(nt, init, sizeof(init));

	/* Plain gradient descent with a fixed step. */
	for (it = 0; it <= SF_GD_ITERS; it++)
	{
		e_gd[it] = film_area(gd);
		if (it == SF_GD_ITERS)
			break;
		film_gradient(gd, active, g);
		for (i = 0; i < N3; i++)
			gd[i] -= (qaws_scalar)(1.5 * g[i]);
	}

	/* Newton-CG: solve H p = -g with conjugate gradients on HVPs, then a
	   backtracking line search on the area. */
	for (it = 0; it <= SF_NEWTON_ITERS; it++)
	{
		qaws_scalar p[N3], r[N3], dd[N3], hd[N3];
		qaws_surface* sf;
		double rr, e0;
		int k;
		e_nt[it] = e0 = film_area(nt);
		if (it == SF_NEWTON_ITERS)
			break;
		film_gradient(nt, active, g);
		sf = film_surface(nt);
		memset(p, 0, sizeof(p));
		for (i = 0; i < N3; i++)
			r[i] = -g[i];
		memcpy(dd, r, sizeof(r));
		rr = vdot(r, r, N3);
		for (k = 0; k < 25 && rr > 1e-24; k++)
		{
			double dhd, alpha, rr_new;
			film_hvp(sf, active, dd, hd);
			hvp_calls++;
			dhd = vdot(dd, hd, N3);
			if (dhd <= 0)
			{
				if (k == 0)
					memcpy(p, r, sizeof(p));
				break;
			}
			alpha = rr / dhd;
			for (i = 0; i < N3; i++)
			{
				p[i] += (qaws_scalar)(alpha * dd[i]);
				r[i] -= (qaws_scalar)(alpha * hd[i]);
			}
			rr_new = vdot(r, r, N3);
			for (i = 0; i < N3; i++)
				dd[i] = (qaws_scalar)(r[i] + rr_new / rr * dd[i]);
			rr = rr_new;
		}
		qaws_surface_destroy(sf);
		{
			double step = 1;
			qaws_scalar trial[N3];
			for (k = 0; k < 20; k++)
			{
				for (i = 0; i < N3; i++)
					trial[i] = (qaws_scalar)(nt[i] + step * p[i]);
				if (film_area(trial) < e0)
					break;
				step *= 0.5;
			}
			if (k < 20)
				memcpy(nt, trial, sizeof(trial));
		}
	}

	svg_open(&s, "showcase/10_soap_film.svg", 1200, 600, "Soap film: area minimization with Hessian-vector products",
		"Bicubic B-spline, boundary fixed, interior heights free (masks). Newton steps solve H p = -g by conjugate gradients on the library's direct HVPs.");
	{
		viewport a = { 30, 80, 370, 500, 0, 0, 0, 0 };
		viewport b = { 415, 80, 370, 500, 0, 0, 0, 0 };
		viewport lv = { 800, 80, 370, 300, 0, 0, 0, 0 };
		projection pa = { 215, 300, 200, 1.0 }, pb = { 600, 300, 200, 1.0 };
		qaws_surface* s0 = film_surface(init);
		qaws_surface* s1 = film_surface(nt);
		double emin = e_nt[SF_NEWTON_ITERS], xy0[2 * (SF_GD_ITERS + 1)], xy1[2 * (SF_NEWTON_ITERS + 1)];
		double lo = 1e300, hi = -1e300;

		svg_panel(&s, &a, "initial interior: |mean curvature|");
		draw_quads(&s, &pa, s0, 26, 26, film_mean_abs, s0, 0, 1.5);
		svg_panel(&s, &b, "after Newton: mean curvature driven toward 0");
		draw_quads(&s, &pb, s1, 26, 26, film_mean_abs, s1, 0, 1.5);
		svg_colorbar(&s, b.x0 + 12, b.y0 + b.h - 24, 200, 8, "|H| = 0", "1.5");

		/* Convergence: log10(E - E_min) per iteration. */
		svg_panel(&s, &lv, "log10(area - final area) per iteration");
		for (i = 0; i <= SF_GD_ITERS; i++)
		{
			double l = log10(e_gd[i] - emin + 1e-16);
			if (l < lo) lo = l;
			if (l > hi) hi = l;
		}
		for (i = 0; i <= SF_NEWTON_ITERS; i++)
		{
			double l = log10(e_nt[i] - emin + 1e-16);
			if (l < lo) lo = l;
			if (l > hi) hi = l;
		}
		if (lo < -14) lo = -14;
		for (i = 0; i <= SF_GD_ITERS; i++)
		{
			double l = log10(e_gd[i] - emin + 1e-16);
			if (l < lo) l = lo;
			xy0[2 * i] = lv.x0 + 14 + (lv.w - 28) * i / (double)SF_GD_ITERS;
			xy0[2 * i + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - l) / (hi - lo);
		}
		for (i = 0; i <= SF_NEWTON_ITERS; i++)
		{
			double l = log10(e_nt[i] - emin + 1e-16);
			if (l < lo) l = lo;
			xy1[2 * i] = lv.x0 + 14 + (lv.w - 28) * i / (double)SF_GD_ITERS;
			xy1[2 * i + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - l) / (hi - lo);
		}
		svg_polyline(&s, xy0, SF_GD_ITERS + 1, "#8c959f", 2, 1, 0);
		svg_polyline(&s, xy1, SF_NEWTON_ITERS + 1, "#0969da", 2.6, 1, 0);
		for (i = 0; i <= SF_NEWTON_ITERS; i++)
			svg_circle(&s, xy1[2 * i], xy1[2 * i + 1], 3, "#0969da", "#ffffff");
		svg_text(&s, lv.x0 + lv.w - 14, xy0[2 * SF_GD_ITERS + 1] - 8, 12, "#57606a", "end", "gradient descent");
		svg_text(&s, lv.x0 + 120, lv.y0 + 60, 12, "#0969da", "start", "Newton-CG on HVPs");

		sprintf(buf, "initial area %.6f", e_gd[0]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 26, 13, "#24292f", "start", buf);
		sprintf(buf, "gradient descent, %d steps: %.8f", SF_GD_ITERS, e_gd[SF_GD_ITERS]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 48, 13, "#57606a", "start", buf);
		sprintf(buf, "Newton-CG, %d steps (%d HVPs): %.8f", SF_NEWTON_ITERS, hvp_calls, e_nt[SF_NEWTON_ITERS]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 70, 13, "#0969da", "start", buf);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 98, 12, "#57606a", "start", "API: qaws_surface_functional_gradient / _hvp");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 116, 12, "#57606a", "start", "masks: interior elements, z component");
		qaws_surface_destroy(s0);
		qaws_surface_destroy(s1);
	}
	svg_close(&s);
	printf("10_soap_film: area %.8f -> GD %.8f, Newton %.8f (%d HVPs)\n", e_gd[0], e_gd[SF_GD_ITERS], e_nt[SF_NEWTON_ITERS], hvp_calls);
}


