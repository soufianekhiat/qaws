/* ================================================================== */
/*  10. Exact ray casting of a NURBS surface                          */
/* ================================================================== */

#define RC_W 120
#define RC_H 90

static void demo_ray_cast(void)
{
	/* a rational biquadratic wave, 4 x 4 control points, 2 x 2 patches */
	qaws_vec3 cps[16];
	qaws_scalar ws[16], kn[7] = { 0, 0, 0, 1, 2, 2, 2 };
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	qaws_exact_surface* e = NULL;
	qaws_exact_desc desc;
	static double depth[RC_W * RC_H], uu[RC_W * RC_H], vv[RC_W * RC_H];
	static int nhit[RC_W * RC_H];
	unsigned int i, j, total = 0, multi = 0, failed = 0;
	double tmin = HUGE_VAL, tmax = -HUGE_VAL, cell = 4.0;
	/* view: rays along dir, the image plane spanned by ex, ey (all dyadic) */
	double dir[3] = { -0.75, -1, -0.375 }, ex[3] = { 1, -0.75, 0 }, ey[3] = { -0.1875, -0.25, 1 }, o[3] = { 1.5, 1.5, 0 };
	double amin = HUGE_VAL, amax = -HUGE_VAL, bmin = HUGE_VAL, bmax = -HUGE_VAL, half, ac, bc;
	char buf[300];
	svg sv;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		{
			cps[i * 4 + j].x = (qaws_scalar)i;
			cps[i * 4 + j].y = (qaws_scalar)j;
			cps[i * 4 + j].z = (qaws_scalar)(((i + j) % 2) ? 1.25 : -0.75) * (qaws_scalar)(i == 0 || i == 3 ? 0.5 : 1);
			ws[i * 4 + j] = (qaws_scalar)((i == 1 && j == 2) ? 4 : 1);
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cps;
	d.u_point_count = 4;
	d.v_point_count = 4;
	d.weights = ws;
	d.u_knots = kn;
	d.u_knot_count = 7;
	d.v_knots = kn;
	d.v_knot_count = 7;
	qaws_surface_create_nurbs(&d, &s);
	qaws_exact_desc_default(&desc);
	qaws_exact_surface_prepare(&desc, s, &e, NULL);
	/* the window: the control net projected on (ex, ey), with a margin */
	for (i = 0; i < 16; i++)
	{
		double r[3] = { cps[i].x - o[0], cps[i].y - o[1], cps[i].z - o[2] };
		double pa = (r[0] * ex[0] + r[1] * ex[1] + r[2] * ex[2]) / (ex[0] * ex[0] + ex[1] * ex[1] + ex[2] * ex[2]);
		double pb = (r[0] * ey[0] + r[1] * ey[1] + r[2] * ey[2]) / (ey[0] * ey[0] + ey[1] * ey[1] + ey[2] * ey[2]);
		if (pa < amin) amin = pa;
		if (pa > amax) amax = pa;
		if (pb < bmin) bmin = pb;
		if (pb > bmax) bmax = pb;
	}
	ac = 0.5 * (amin + amax);
	bc = 0.5 * (bmin + bmax);
	half = 0.55 * ((amax - amin) / RC_W > (bmax - bmin) / RC_H ? (amax - amin) : (bmax - bmin) * RC_W / RC_H);
	for (j = 0; j < RC_H; j++)
		for (i = 0; i < RC_W; i++)
		{
			double a = ldexp(nearbyint(ldexp(ac - half + 2 * half * (i + 0.5) / RC_W, 12)), -12);
			double b = ldexp(nearbyint(ldexp(bc + half * RC_H / RC_W - 2 * half * RC_H / RC_W * (j + 0.5) / RC_H, 12)), -12);
			double p0[3], p1[3];
			qaws_exact_surface_hit hits[8];
			unsigned int n = 0, c, k = j * RC_W + i;
			for (c = 0; c < 3; c++)
			{
				p0[c] = o[c] + a * ex[c] + b * ey[c] - 4 * dir[c];
				p1[c] = p0[c] + dir[c];
			}
			nhit[k] = 0;
			if (qaws_exact_surface_line_hits(e, p0, p1, hits, 8, &n) != QAWS_STATUS_OK)
			{
				nhit[k] = -1;
				failed++;
				continue;
			}
			nhit[k] = (int)n;
			if (n)
			{
				depth[k] = hits[0].t_lo;
				uu[k] = 0.5 * (hits[0].u_lo + hits[0].u_hi);
				vv[k] = 0.5 * (hits[0].v_lo + hits[0].v_hi);
				if (depth[k] < tmin) tmin = depth[k];
				if (depth[k] > tmax) tmax = depth[k];
				total += n;
				multi += n > 1;
			}
		}
	svg_open(&sv, "showcase/exact10_ray_cast.svg", 1040, 500, "Exact ray casting of a NURBS surface",
		"Every pixel is a certified line / surface intersection: Bezout elimination per patch, roots isolated, u, v and t enclosed.");
	fprintf(sv.f, "<text x=\"20\" y=\"96\" font-size=\"14\" font-weight=\"600\" fill=\"#24292f\">first hit: depth shading, exact (u, v) checker</text>\n");
	fprintf(sv.f, "<text x=\"530\" y=\"96\" font-size=\"14\" font-weight=\"600\" fill=\"#24292f\">intersections per ray: grey 1, blue 2, purple 3, orange 4</text>\n");
	for (j = 0; j < RC_H; j++)
		for (i = 0; i < RC_W; i++)
		{
			unsigned int k = j * RC_W + i;
			char col[32];
			double x = 20 + i * cell, y = 106 + j * cell;
			if (nhit[k] > 0)
			{
				double sh = (depth[k] - tmin) / (tmax - tmin + 1e-300);
				int chk = ((int)floor(uu[k] * 4) + (int)floor(vv[k] * 4)) & 1;
				int g = (int)(235 - 150 * sh) - (chk ? 30 : 0);
				sprintf(col, "rgb(%d,%d,%d)", g, g + (chk ? 0 : 8), g + 20 > 255 ? 255 : g + 20);
				fprintf(sv.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n", x, y, cell + 0.3, cell + 0.3, col);
			}
			else if (nhit[k] < 0)
				fprintf(sv.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#cf222e\"/>\n", x, y, cell + 0.3, cell + 0.3);
			if (nhit[k] != 0)
			{
				static char const* const ncol[5] = { "#cf222e", "#d0d7de", "#0969da", "#8250df", "#bc4c00" };
				fprintf(sv.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n", 530 + i * cell, y, cell + 0.3, cell + 0.3,
					ncol[nhit[k] < 0 ? 0 : (nhit[k] > 4 ? 4 : nhit[k])]);
			}
		}
	sprintf(buf, "%d x %d rays: %u hits, %u rays crossing the surface more than once, %u uncertified (tangent rays)", RC_W, RC_H, total, multi, failed);
	svg_text(&sv, 20, 484, 13, "#57606a", "start", buf);
	svg_close(&sv);
	printf("exact10_ray_cast: %u hits, %u multi, %u uncertified\n", total, multi, failed);
	qaws_exact_surface_destroy(e);
	qaws_surface_destroy(s);
}
