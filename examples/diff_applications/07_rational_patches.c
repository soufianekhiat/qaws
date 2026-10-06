/* ================================================================== */
/*  7. Rational patches: freeing the NURBS weights after ADMM         */
/* ================================================================== */

#define RP_ITERS 300

typedef struct rational_run
{
	double rms_admm, rms_ctl, rms_rat, gap_rat, wmin, wmax;
	double hist_ctl[RP_ITERS / 10 + 1], hist_rat[RP_ITERS / 10 + 1];
} rational_run;

static void rational_experiment(double amp, rational_run* out, patch_fit* rat_out, qaws_scalar (*w_out)[MP_CP], mesh* m_out)
{
	static patch_fit admm[MP_PATCHES], ctl[MP_PATCHES];
	static qaws_scalar wc[MP_PATCHES][MP_CP];
	double rms[MP_OUTER], gap[MP_OUTER];
	int nodes, f, k;
	g_blob_amp = amp;
	mesh_blob(m_out);
	mp_setup(m_out, admm, &nodes);
	mp_solve(m_out, admm, nodes, 0.03, rms, gap);
	out->rms_admm = rms[MP_OUTER - 1];
	for (f = 0; f < MP_PATCHES; f++)
	{
		size_t n = (size_t)admm[f].n;
		ctl[f] = admm[f];
		rat_out[f] = admm[f];
		/* own copies of the foot points */
		ctl[f].uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		rat_out[f].uv = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n);
		memcpy(ctl[f].uv, admm[f].uv, sizeof(qaws_scalar) * 2 * n);
		memcpy(rat_out[f].uv, admm[f].uv, sizeof(qaws_scalar) * 2 * n);
	}
	out->rms_ctl = mp_refine_rational(m_out, ctl, nodes, wc, RP_ITERS, 0, out->hist_ctl);
	out->rms_rat = mp_refine_rational(m_out, rat_out, nodes, w_out, RP_ITERS, 1, out->hist_rat);
	out->gap_rat = mp_seam_gap(rat_out, nodes);
	out->wmin = 1e30;
	out->wmax = 0;
	for (f = 0; f < MP_PATCHES; f++)
		for (k = 0; k < MP_CP; k++)
		{
			if (w_out[f][k] < out->wmin) out->wmin = w_out[f][k];
			if (w_out[f][k] > out->wmax) out->wmax = w_out[f][k];
		}
	for (f = 0; f < MP_PATCHES; f++)
	{
		free(admm[f].pts); free(admm[f].uv);
		free(ctl[f].uv);
	}
}

static void app_rational_patches(void)
{
	static patch_fit rat_blob[MP_PATCHES], rat_sphere[MP_PATCHES];
	static qaws_scalar w_blob[MP_PATCHES][MP_CP], w_sphere[MP_PATCHES][MP_CP];
	rational_run rb, rs;
	mesh mb, ms;
	svg s;
	char buf[256];
	int q;
	rational_experiment(1.0, &rb, rat_blob, w_blob, &mb);
	rational_experiment(0.0, &rs, rat_sphere, w_sphere, &ms);

	svg_open(&s, "showcase/app7_rational_patches.svg", 1240, 700, "Rational patches: freeing the NURBS weights",
		"After ADMM: refine control points only (weights 1) or control points + log-weights through the rational "
		"surface adjoints; seam copies stay averaged, so the patch set stays closed.");
	for (q = 0; q < 2; q++)
	{
		rational_run const* r = q ? &rs : &rb;
		qaws_scalar (*w)[MP_CP] = q ? w_sphere : w_blob;
		viewport lv = { 20 + q * 610, 90, 590, 250, 0, 0, 0, 0 };
		double lo = 1e300, hi = -1e300, x0[2 * (RP_ITERS / 10 + 1)], x1[2 * (RP_ITERS / 10 + 1)];
		int n = RP_ITERS / 10 + 1, k, f, i, j;
		for (k = 0; k < n; k++)
		{
			double a = log10(r->hist_ctl[k]), b = log10(r->hist_rat[k]);
			if (a < lo) lo = a; if (b < lo) lo = b;
			if (a > hi) hi = a; if (b > hi) hi = b;
		}
		sprintf(buf, "%s: log10 RMS distance per step", q ? "sphere mesh" : "blob mesh");
		svg_panel(&s, &lv, buf);
		for (k = 0; k < n; k++)
		{
			x0[2 * k] = x1[2 * k] = lv.x0 + 14 + (lv.w - 28) * k / (double)(n - 1);
			x0[2 * k + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(r->hist_ctl[k])) / (hi - lo + 1e-12);
			x1[2 * k + 1] = lv.y0 + 34 + (lv.h - 54) * (hi - log10(r->hist_rat[k])) / (hi - lo + 1e-12);
		}
		svg_polyline(&s, x0, n, "#8c959f", 2.2, 1, 0);
		svg_polyline(&s, x1, n, "#0969da", 2.6, 1, 0);
		svg_text(&s, lv.x0 + lv.w - 14, x0[2 * n - 1] - 8, 11, "#57606a", "end", "control points only");
		svg_text(&s, lv.x0 + lv.w - 14, x1[2 * n - 1] + 16, 11, "#0969da", "end", "control points + weights");
		sprintf(buf, "RMS: ADMM %.2e, control points %.2e, + weights %.2e (%.0f%% lower)", r->rms_admm, r->rms_ctl, r->rms_rat,
			100.0 * (1 - r->rms_rat / r->rms_ctl));
		svg_text(&s, lv.x0, lv.y0 + lv.h + 24, 13, "#24292f", "start", buf);
		sprintf(buf, "weights %.3f .. %.3f, max seam gap %.1e", r->wmin, r->wmax, r->gap_rat);
		svg_text(&s, lv.x0, lv.y0 + lv.h + 44, 13, "#57606a", "start", buf);
		/* weight grids per patch */
		for (f = 0; f < MP_PATCHES; f++)
		{
			double gx = lv.x0 + 8 + f * 97, gy = lv.y0 + lv.h + 70, cs = 12;
			double wl = q ? 0.96 : 0.75, wh = q ? 1.04 : 1.25;
			for (i = 0; i < MP_N; i++)
				for (j = 0; j < MP_N; j++)
				{
					char col[32];
					heat((w[f][i * MP_N + j] - wl) / (wh - wl), col);
					fprintf(s.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n",
						gx + j * cs, gy + i * cs, cs - 1, cs - 1, col);
				}
			sprintf(buf, "patch %d", f + 1);
			svg_text(&s, gx, gy + MP_N * cs + 14, 11, "#57606a", "start", buf);
		}
		sprintf(buf, "%.2f", q ? 0.96 : 0.75);
		{
			char hb[16];
			sprintf(hb, "%.2f", q ? 1.04 : 1.25);
			svg_colorbar(&s, lv.x0 + 8, lv.y0 + lv.h + 175, 220, 8, buf, hb);
		}
	}
	svg_text(&s, 20, 660, 13, "#57606a", "start",
		"gradient: 2 (S - X) at the foot points -> qaws_surface_eval_batch_adjoint with CONTROL_POINTS and WEIGHTS views "
		"(rational quotient rule), weights parameterized by their log to stay positive");
	svg_close(&s);
	printf("7_rational: blob admm %.5f ctl %.5f rat %.5f (w %.3f..%.3f, gap %.1e); sphere admm %.5f ctl %.5f rat %.5f (w %.3f..%.3f)\n",
		rb.rms_admm, rb.rms_ctl, rb.rms_rat, rb.wmin, rb.wmax, rb.gap_rat, rs.rms_admm, rs.rms_ctl, rs.rms_rat, rs.wmin, rs.wmax);
}

