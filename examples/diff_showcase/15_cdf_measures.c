/* ================================================================== */
/*  15. Inverse-CDF sampling: arc length, curvature, density          */
/* ================================================================== */

#define CDF_N 40

/* rho = 0.25 + 3 exp(-|x - c|^2 / (2 s^2)): a bright spot in the plane */
static qaws_scalar cdf_spot(qaws_vec3 p, void* user, qaws_vec3* g, qaws_scalar* H)
{
	double cx = 4.2, cy = 1.0, s2 = 0.55 * 0.55, dx = p.x - cx, dy = p.y - cy;
	double e = 3.0 * exp(-(dx * dx + dy * dy) / (2 * s2));
	(void)user;
	g->x = (qaws_scalar)(-e * dx / s2);
	g->y = (qaws_scalar)(-e * dy / s2);
	g->z = 0;
	H[0] = (qaws_scalar)(e * (dx * dx / s2 - 1) / s2);
	H[1] = (qaws_scalar)(e * dx * dy / (s2 * s2));
	H[2] = 0;
	H[3] = (qaws_scalar)(e * (dy * dy / s2 - 1) / s2);
	H[4] = 0;
	H[5] = 0;
	return (qaws_scalar)(0.25 + e);
}

static void demo_cdf_measures(void)
{
	static qaws_scalar const cps[2 * 10] = { 0.0, 0.3, 0.6, 1.5, 0.9, 0.2, 1.3, 1.6, 2.6, 1.7, 3.4, 0.1, 4.6, 0.4, 5.2, 1.7, 6.4, 1.2, 7.2, 0.4 };
	qaws_sample_measure_desc curv = { QAWS_MEASURE_CURVATURE, (qaws_scalar)0.4, NULL, NULL };
	qaws_sample_measure_desc dens = { QAWS_MEASURE_DENSITY, 0, cdf_spot, NULL };
	qaws_sample_measure_desc const* ms[3] = { NULL, &curv, &dens };
	char const* names[3] = { "arc length: rho = 1", "curvature: rho = sqrt(0.4^2 + kappa^2)", "density field: rho(C), a bright spot" };
	qaws_cdf_target tg[CDF_N];
	qaws_curve* c = bspline_2d(cps, 10);
	double xy[2 * 400];
	char buf[160];
	int m, i;
	svg s;
	for (i = 0; i < CDF_N; i++)
	{
		tg[i].distance = 0;
		tg[i].fraction = (qaws_scalar)((double)i / (CDF_N - 1));
	}
	svg_open(&s, "showcase/15_cdf_measures.svg", 1200, 560, "Inverse-CDF sampling of a curve under three measures",
		"40 samples at equal measure, t solving M(t) = i/39 M_total; arrows: first order sample tangents for the red control point moving up.");
	for (m = 0; m < 3; m++)
	{
		viewport v = { 30 + m * 390, 80, 370, 440, -0.4, 7.6, -0.6, 2.4 };
		qaws_cdf_sample val[CDF_N], t1[CDF_N];
		qaws_scalar dir[20], total = 0;
		qaws_field_view fv;
		qaws_diff_views views;
		unsigned int k_star = 5;
		svg_panel(&s, &v, names[m]);
		if (m == 2)
		{
			/* the density as a heat background, clipped to the panel */
			int gx, gy;
			fprintf(s.f, "<clipPath id=\"cdfclip\"><rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\"/></clipPath><g clip-path=\"url(#cdfclip)\">\n",
				v.x0, v.y0 + 28, v.w, v.h - 28);
			for (gy = 0; gy < 40; gy++)
				for (gx = 0; gx < 40; gx++)
				{
					qaws_vec3 p, g;
					qaws_scalar H[6];
					double x0 = v.xmin + (v.xmax - v.xmin) * gx / 40.0, y0 = v.ymin + (v.ymax - v.ymin) * gy / 40.0, r;
					char col[32];
					p = v3((qaws_scalar)(x0 + 0.1), (qaws_scalar)(y0 + 0.037), 0);
					r = cdf_spot(p, NULL, &g, H);
					if (r < 0.3)
						continue;
					heat((r - 0.25) / 3.0, col);
					fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\" opacity=\"0.35\"/>\n",
						vx(&v, x0), vy(&v, y0 + (v.ymax - v.ymin) / 40.0), v.w / 40.0 + 0.5, v.h / 40.0 + 0.5, col);
				}
			fprintf(s.f, "</g>\n");
		}
		curve_polyline(c, &v, xy, 400);
		svg_polyline(&s, xy, 400, "#57606a", 2, 1, 0);
		memset(dir, 0, sizeof(dir));
		dir[2 * k_star + 1] = 1;
		views = one_field(&fv, QAWS_FIELD_CONTROL_POINTS, dir, 10, 2);
		qaws_curve_cdf_sample_tangent(NULL, c, ms[m], tg, NULL, NULL, CDF_N, 0, &views, val, t1, NULL, &total);
		for (i = 0; i < CDF_N; i++)
		{
			double x0 = vx(&v, val[i].position.x), y0 = vy(&v, val[i].position.y), k1 = 0.5 * v.w / (v.xmax - v.xmin);
			svg_line(&s, x0, y0, x0 + k1 * t1[i].position.x, y0 - k1 * t1[i].position.y, "#0969da", 1.4, 0.85);
			svg_circle(&s, x0, y0, 3.2, "#cf222e", "#ffffff");
		}
		svg_circle(&s, vx(&v, cps[2 * k_star]), vy(&v, cps[2 * k_star + 1]), 5, "#cf222e", "#24292f");
		sprintf(buf, "total measure %.3f", (double)total);
		svg_text(&s, v.x0 + 12, v.y0 + v.h - 12, 12, "#57606a", "start", buf);
	}
	svg_close(&s);
	qaws_curve_destroy(c);
	printf("15_cdf_measures: written\n");
}

