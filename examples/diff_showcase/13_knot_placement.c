/* ================================================================== */
/*  13. Knot placement: fitting a surface with a crease by moving its  */
/*      knots together with its control points                        */
/* ================================================================== */

#define KP_NU 8
#define KP_NV 6
#define KP_UK (KP_NU + 4)
#define KP_VK (KP_NV + 4)
#define KP_SU 48
#define KP_SV 36
#define KP_ITERS 500

static double kp_target(double u, double v)
{
	/* a crease along u = 0.68 and a gentle wave across v */
	return 0.35 * tanh(18.0 * (u - 0.68)) + 0.12 * sin(2.0 * PI * v) * (1.0 - u);
}

static qaws_surface* kp_surface(qaws_scalar const* cps, qaws_scalar const* uk, qaws_scalar const* vk)
{
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	memset(&d, 0, sizeof(d));
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = (qaws_vec3 const*)cps;
	d.u_point_count = KP_NU;
	d.v_point_count = KP_NV;
	d.u_knots = uk;
	d.u_knot_count = KP_UK;
	d.v_knots = vk;
	d.v_knot_count = KP_VK;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

static void kp_samples(qaws_scalar* us, qaws_scalar* vs)
{
	int i, j;
	for (i = 0; i < KP_SU; i++)
		for (j = 0; j < KP_SV; j++)
		{
			us[i * KP_SV + j] = (qaws_scalar)((i + 0.5) / KP_SU);
			vs[i * KP_SV + j] = (qaws_scalar)((j + 0.5) / KP_SV);
		}
}

/* Mean squared height error and its gradient (z of control points, and
   interior knots when move_knots). */
static double kp_loss(qaws_scalar const* cps, qaws_scalar const* uk, qaws_scalar const* vk, int move_knots,
	qaws_scalar* g_cp, qaws_scalar* g_uk, qaws_scalar* g_vk)
{
	static qaws_scalar us[KP_SU * KP_SV], vs[KP_SU * KP_SV];
	static qaws_surface_jet jets[KP_SU * KP_SV], bars[KP_SU * KP_SV];
	static int init = 0;
	unsigned char uact[KP_UK], vact[KP_VK];
	qaws_surface* s = kp_surface(cps, uk, vk);
	qaws_field_view fv[3];
	qaws_diff_views views;
	double loss = 0;
	int k, n = KP_SU * KP_SV;

	if (!init)
	{
		kp_samples(us, vs);
		init = 1;
	}
	for (k = 0; k < n; k++)
	{
		qaws_surface_jet j;
		double e;
		qaws_surface_eval_jet(s, us[k], vs[k], QAWS_SJET_P, &j);
		jets[k] = j;
		e = j.d[0].z - kp_target(us[k], vs[k]);
		loss += e * e;
		memset(&bars[k], 0, sizeof(bars[k]));
		bars[k].d[0].z = (qaws_scalar)(2 * e / n);
		bars[k].channels = QAWS_SJET_P;
	}
	if (g_cp)
	{
		memset(g_cp, 0, sizeof(qaws_scalar) * KP_NU * KP_NV * 3);
		memset(g_uk, 0, sizeof(qaws_scalar) * KP_UK);
		memset(g_vk, 0, sizeof(qaws_scalar) * KP_VK);
		for (k = 0; k < KP_UK; k++) uact[k] = (unsigned char)(k > 3 && k < KP_UK - 4);
		for (k = 0; k < KP_VK; k++) vact[k] = (unsigned char)(k > 3 && k < KP_VK - 4);
		fv[0] = qaws_field_view_make(QAWS_FIELD_CONTROL_POINTS, g_cp, KP_NU * KP_NV, 3);
		fv[0].component_mask = 1u << 2;
		fv[1] = qaws_field_view_make(QAWS_FIELD_U_KNOTS, g_uk, KP_UK, 1);
		fv[1].active = uact;
		fv[2] = qaws_field_view_make(QAWS_FIELD_V_KNOTS, g_vk, KP_VK, 1);
		fv[2].active = vact;
		views.fields = fv;
		views.field_count = move_knots ? 3u : 1u;
		views.children = NULL;
		views.child_count = 0;
		qaws_surface_eval_batch_adjoint(NULL, s, us, vs, (unsigned int)n, QAWS_SJET_P, bars, &views, NULL, NULL);
	}
	qaws_surface_destroy(s);
	return loss / n;
}

/* Interior knots stay ordered with a minimum gap. */
static void kp_project(qaws_scalar* k, int count)
{
	int i, lo = 4, hi = count - 5;
	for (i = lo; i <= hi; i++)
	{
		double min = (i == lo ? 0.0 : k[i - 1]) + 0.025;
		if (k[i] < min) k[i] = (qaws_scalar)min;
	}
	for (i = hi; i >= lo; i--)
	{
		double max = (i == hi ? 1.0 : k[i + 1]) - 0.025;
		if (k[i] > max) k[i] = (qaws_scalar)max;
	}
}

static void kp_fit(qaws_scalar* cps, qaws_scalar* uk, qaws_scalar* vk, int move_knots, double* loss)
{
	static qaws_scalar g_cp[KP_NU * KP_NV * 3], g_uk[KP_UK], g_vk[KP_VK];
	qaws_scalar z[KP_NU * KP_NV], gz[KP_NU * KP_NV];
	adam a_cp, a_k;
	int it, i;
	memset(&a_cp, 0, sizeof(a_cp));
	memset(&a_k, 0, sizeof(a_k));
	for (it = 0; it <= KP_ITERS; it++)
	{
		loss[it] = kp_loss(cps, uk, vk, move_knots, g_cp, g_uk, g_vk);
		if (it == KP_ITERS)
			break;
		for (i = 0; i < KP_NU * KP_NV; i++)
		{
			z[i] = cps[i * 3 + 2];
			gz[i] = g_cp[i * 3 + 2];
		}
		adam_step(&a_cp, z, gz, KP_NU * KP_NV, 0.02 * (1.0 - 0.7 * it / (double)KP_ITERS));
		for (i = 0; i < KP_NU * KP_NV; i++)
			cps[i * 3 + 2] = z[i];
		if (move_knots)
		{
			qaws_scalar kk[KP_UK + KP_VK], gk[KP_UK + KP_VK];
			memcpy(kk, uk, sizeof(qaws_scalar) * KP_UK);
			memcpy(kk + KP_UK, vk, sizeof(qaws_scalar) * KP_VK);
			memcpy(gk, g_uk, sizeof(qaws_scalar) * KP_UK);
			memcpy(gk + KP_UK, g_vk, sizeof(qaws_scalar) * KP_VK);
			adam_step(&a_k, kk, gk, KP_UK + KP_VK, 0.004 * (1.0 - 0.7 * it / (double)KP_ITERS));
			memcpy(uk, kk, sizeof(qaws_scalar) * KP_UK);
			memcpy(vk, kk + KP_UK, sizeof(qaws_scalar) * KP_VK);
			kp_project(uk, KP_UK);
			kp_project(vk, KP_VK);
		}
	}
}

static void kp_error_map(svg* s, viewport const* v, qaws_scalar const* cps, qaws_scalar const* uk, qaws_scalar const* vk,
	double vmax, char* buf)
{
	qaws_surface* sf = kp_surface(cps, uk, vk);
	int i, j;
	double cw = v->w / 60.0, ch = (v->h - 30) / 45.0;
	for (i = 0; i < 60; i++)
		for (j = 0; j < 45; j++)
		{
			qaws_surface_jet jt;
			double u = (i + 0.5) / 60.0, w = (j + 0.5) / 45.0, e;
			char col[32];
			qaws_surface_eval_jet(sf, (qaws_scalar)u, (qaws_scalar)w, QAWS_SJET_P, &jt);
			e = fabs(jt.d[0].z - kp_target(u, w));
			heat(e / vmax, col);
			fprintf(s->f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\"/>\n",
				v->x0 + i * cw, v->y0 + 30 + (44 - j) * ch, cw + 0.4, ch + 0.4, col);
		}
	for (i = 4; i < KP_UK - 4; i++)
		svg_line(s, v->x0 + uk[i] * v->w, v->y0 + 30, v->x0 + uk[i] * v->w, v->y0 + v->h, "#ffffff", 2, 0.95);
	for (i = 4; i < KP_VK - 4; i++)
		svg_line(s, v->x0, v->y0 + v->h - vk[i] * (v->h - 30), v->x0 + v->w, v->y0 + v->h - vk[i] * (v->h - 30), "#ffffff", 2, 0.95);
	svg_line(s, v->x0 + 0.68 * v->w, v->y0 + 30, v->x0 + 0.68 * v->w, v->y0 + v->h, "#1b1f24", 1.2, 0.6);
	qaws_surface_destroy(sf);
	(void)buf;
}

static void demo_knot_placement(void)
{
	static qaws_scalar cps0[KP_NU * KP_NV * 3], cps1[KP_NU * KP_NV * 3];
	qaws_scalar uk0[KP_UK], vk0[KP_VK], uk1[KP_UK], vk1[KP_VK];
	static double l0[KP_ITERS + 1], l1[KP_ITERS + 1];
	int i, j;
	svg s;
	char buf[256];

	for (i = 0; i < KP_NU; i++)
		for (j = 0; j < KP_NV; j++)
		{
			qaws_scalar* p = &cps0[(i * KP_NV + j) * 3];
			p[0] = (qaws_scalar)(i / (double)(KP_NU - 1));
			p[1] = (qaws_scalar)(j / (double)(KP_NV - 1));
			p[2] = 0;
		}
	for (i = 0; i < KP_UK; i++)
		uk0[i] = (qaws_scalar)(i < 4 ? 0.0 : (i >= KP_UK - 4 ? 1.0 : (i - 3) / (double)(KP_NU - 3)));
	for (i = 0; i < KP_VK; i++)
		vk0[i] = (qaws_scalar)(i < 4 ? 0.0 : (i >= KP_VK - 4 ? 1.0 : (i - 3) / (double)(KP_NV - 3)));
	memcpy(cps1, cps0, sizeof(cps0));
	memcpy(uk1, uk0, sizeof(uk0));
	memcpy(vk1, vk0, sizeof(vk0));

	kp_fit(cps0, uk0, vk0, 0, l0);
	kp_fit(cps1, uk1, vk1, 1, l1);

	svg_open(&s, "showcase/13_knot_placement.svg", 1200, 600, "Knot placement: fitting a crease by moving knots",
		"Bicubic B-spline (8 x 6) fit to a height field with a crease at u = 0.68 (dark line). Interior U and V knots "
		"(white lines) move with the control points through knot adjoints.");
	{
		viewport a = { 24, 80, 370, 470, 0, 1, 0, 1 }, b = { 410, 80, 370, 470, 0, 1, 0, 1 };
		viewport lv = { 800, 80, 376, 300, 0, 0, 0, 0 };
		double vmax = 0.06, xy0[2 * (KP_ITERS + 1)], xy1[2 * (KP_ITERS + 1)], lo = 1e300, hi = -1e300;
		int it;
		svg_panel(&s, &a, "control points only, uniform knots: |error|");
		kp_error_map(&s, &a, cps0, uk0, vk0, vmax, buf);
		svg_panel(&s, &b, "control points + knots: |error|");
		kp_error_map(&s, &b, cps1, uk1, vk1, vmax, buf);
		svg_colorbar(&s, b.x0 + 12, b.y0 + b.h + 8, 200, 8, "0", "0.06");

		for (it = 0; it <= KP_ITERS; it++)
		{
			double a0 = log10(sqrt(l0[it])), a1 = log10(sqrt(l1[it]));
			if (a0 < lo) lo = a0;
			if (a1 < lo) lo = a1;
			if (a0 > hi) hi = a0;
			if (a1 > hi) hi = a1;
		}
		svg_panel(&s, &lv, "log10 RMS error per iteration");
		for (it = 0; it <= KP_ITERS; it++)
		{
			xy0[2 * it] = xy1[2 * it] = lv.x0 + 14 + (lv.w - 28) * it / (double)KP_ITERS;
			xy0[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(sqrt(l0[it]))) / (hi - lo);
			xy1[2 * it + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(sqrt(l1[it]))) / (hi - lo);
		}
		svg_polyline(&s, xy0, KP_ITERS + 1, "#8c959f", 2, 1, 0);
		svg_polyline(&s, xy1, KP_ITERS + 1, "#0969da", 2.6, 1, 0);
		svg_text(&s, lv.x0 + lv.w - 14, xy0[2 * KP_ITERS + 1] - 8, 11, "#57606a", "end", "control points only");
		svg_text(&s, lv.x0 + lv.w - 14, xy1[2 * KP_ITERS + 1] + 16, 11, "#0969da", "end", "control points + knots");

		sprintf(buf, "RMS: uniform knots %.4f   moved knots %.4f", sqrt(l0[KP_ITERS]), sqrt(l1[KP_ITERS]));
		svg_text(&s, lv.x0, lv.y0 + lv.h + 26, 13, "#24292f", "start", buf);
		sprintf(buf, "U knots: %.3f %.3f %.3f %.3f", uk1[4], uk1[5], uk1[6], uk1[7]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 48, 13, "#0969da", "start", buf);
		sprintf(buf, "V knots: %.3f %.3f", vk1[4], vk1[5]);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 70, 13, "#0969da", "start", buf);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 98, 12, "#57606a", "start", "qaws_surface_eval_batch_adjoint with views:");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 116, 12, "#57606a", "start", "control points (z mask), U_KNOTS, V_KNOTS");
		svg_text(&s, lv.x0, lv.y0 + lv.h + 134, 12, "#57606a", "start", "(interior knots active, ordered by projection)");
	}
	svg_close(&s);
	printf("13_knot_placement: RMS uniform %.5f, moved knots %.5f, u knots %.3f %.3f %.3f %.3f\n",
		sqrt(l0[KP_ITERS]), sqrt(l1[KP_ITERS]), uk1[4], uk1[5], uk1[6], uk1[7]);
}

