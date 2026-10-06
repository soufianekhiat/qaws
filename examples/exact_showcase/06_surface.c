/* ================================================================== */
/*  6. Exact NURBS surface normals: up, down or exactly vertical      */
/* ================================================================== */

#define ES_N 120

static void demo_surface(void)
{
	static char const* const facing_colors[3] = { "#f2a900", "#24292f", "#0969da" };   /* down, vertical, up */
	/* x(u) and y(v) each hold still over one span (three equal control
	   coordinates) and x runs back over another: n_z = x'(u) y'(v) is then
	   exactly zero on whole bands and negative on the overhang. Separable
	   weights a_i b_j keep that factorization for the rational surface. */
	static double const fx[8] = { 0, 1, 2, 2, 2, 1.5, 1, 1.5 };
	static double const gy[8] = { 0, 1, 2, 2, 2, 3, 4, 5 };
	static double const a[8] = { 1, 2, 1, 3, 1, 2, 1, 4 }, b[8] = { 1, 3, 1, 2, 2, 1, 3, 1 };
	static signed char exact_map[ES_N * ES_N], f64_map[ES_N * ES_N], wrong_map[ES_N * ES_N];
	static char const* const wrong_colors[3] = { "#f6f8fa", "#f6f8fa", "#cf222e" };
	qaws_vec3 cps[64];
	qaws_scalar ws[64], knots[11] = { 0, 0, 0, 1, 2, 3, 4, 5, 6, 6, 6 };
	qaws_surface_nurbs_desc d;
	qaws_surface* surf = NULL;
	qaws_exact_surface* es = NULL;
	qaws_exact_desc desc;
	qaws_exact_report rep;
	projection pr = { 335, 170, 64, 1.0 };
	unsigned int i, j, pu = 0, pv = 0;
	int n_wrong = 0, n_zero = 0, n_f64_zero = 0;
	char buf[300];
	svg sv;
	for (i = 0; i < 8; i++)
		for (j = 0; j < 8; j++)
		{
			cps[i * 8 + j].x = (qaws_scalar)fx[i];
			cps[i * 8 + j].y = (qaws_scalar)gy[j];
			cps[i * 8 + j].z = (qaws_scalar)ldexp(nearbyint(ldexp(0.9 * sin(0.9 * i) * cos(0.7 * j) + 0.15 * i, 20)), -20);
			ws[i * 8 + j] = (qaws_scalar)(a[i] * b[j]);
		}
	memset(&d, 0, sizeof(d));
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cps;
	d.u_point_count = 8;
	d.v_point_count = 8;
	d.weights = ws;
	d.u_knots = knots;
	d.u_knot_count = 11;
	d.v_knots = knots;
	d.v_knot_count = 11;
	qaws_surface_create_nurbs(&d, &surf);
	qaws_exact_desc_default(&desc);
	qaws_exact_surface_prepare(&desc, surf, &es, &rep);
	qaws_exact_surface_patch_count(es, &pu, &pv);

	/* classify cell centers (u, v on the 2^-24 lattice) */
	for (j = 0; j < ES_N; j++)
		for (i = 0; i < ES_N; i++)
		{
			double u = ldexp(nearbyint(ldexp(6.0 * (i + 0.5) / ES_N, 20)), -20), v = ldexp(nearbyint(ldexp(6.0 * (j + 0.5) / ES_N, 20)), -20);
			double pos[3], nrm[3], nz;
			qaws_surface_eval_result r;
			int k = (int)(j * ES_N + i);
			qaws_exact_surface_evaluate(es, u, v, 0, pos, nrm);
			exact_map[k] = (signed char)(nrm[2] > 0 ? 1 : (nrm[2] < 0 ? -1 : 0));
			qaws_surface_evaluate(surf, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV, &r);
			nz = (double)r.du.x * r.dv.y - (double)r.du.y * r.dv.x;
			f64_map[k] = (signed char)(nz > 0 ? 1 : (nz < 0 ? -1 : 0));
			wrong_map[k] = (signed char)(f64_map[k] != exact_map[k] ? 1 : 0);
			n_wrong += wrong_map[k];
			n_zero += exact_map[k] == 0;
			n_f64_zero += f64_map[k] == 0;
		}

	svg_open(&sv, "showcase/exact6_surface.svg", 1280, 560, "Exact NURBS surface normals: up, down or exactly vertical",
		"A rational biquadratic NURBS surface (36 patches) with vertical bands and an overhang. Blue: n_z &gt; 0, orange: n_z &lt; 0, black: n_z = 0 (vertical).");

	/* 3D view: cells back to front, colored by the exact class */
	{
		int n = ES_N / 3, *order = (int*)malloc(sizeof(int) * (size_t)(n * n)), k, m;
		double* depth = (double*)malloc(sizeof(double) * (size_t)(n * n));
		svg_text(&sv, 20, 100, 14, "#24292f", "start", "the surface, colored by the exact sign of n_z");
		for (k = 0; k < n * n; k++)
		{
			double pos[3];
			qaws_exact_surface_evaluate(es, 6.0 * (k / n + 0.5) / n, 6.0 * (k % n + 0.5) / n, 0, pos, NULL);
			depth[k] = pos[0] + pos[1] - 0.3 * pos[2];
			order[k] = k;
		}
		for (k = 1; k < n * n; k++)
		{
			int key = order[k];
			for (m = k - 1; m >= 0 && depth[order[m]] > depth[key]; m--)
				order[m + 1] = order[m];
			order[m + 1] = key;
		}
		for (k = 0; k < n * n; k++)
		{
			int ci = order[k] / n, cj = order[k] % n, c, cls;
			double x[4], y[4], nrm[3], pos[3];
			static int const du[4] = { 0, 1, 1, 0 }, dv[4] = { 0, 0, 1, 1 };
			for (c = 0; c < 4; c++)
			{
				qaws_exact_surface_evaluate(es, 6.0 * (ci + du[c]) / n, 6.0 * (cj + dv[c]) / n, 0, pos, NULL);
				project(&pr, pos[0], pos[1], pos[2], &x[c], &y[c]);
			}
			qaws_exact_surface_evaluate(es, 6.0 * (ci + 0.5) / n, 6.0 * (cj + 0.5) / n, 0, pos, nrm);
			cls = nrm[2] > 0 ? 2 : (nrm[2] < 0 ? 0 : 1);
			fprintf(sv.f, "<polygon points=\"%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f\" fill=\"%s\" fill-opacity=\"%s\" stroke=\"#ffffff\" stroke-width=\"0.4\"/>\n",
				x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3], facing_colors[cls], cls == 1 ? "0.85" : "0.55");
		}
		free(order);
		free(depth);
		sprintf(buf, "%u x %u exact patches, integers up to %u bits", pu, pv, rep.storage_bits);
		svg_text(&sv, 20, 530, 12, "#57606a", "start", buf);
	}
	or_map(&sv, 500, 110, 3.0, ES_N, exact_map, facing_colors, "exact sign of n_z over (u, v)");
	or_map(&sv, 880, 110, 3.0, ES_N, f64_map, facing_colors, "f64 evaluator: sign of du x dv . z");
	sprintf(buf, "exact: %d of %d cells vertical (n_z = 0 exactly)", n_zero, ES_N * ES_N);
	svg_text(&sv, 500, 490, 12, "#57606a", "start", buf);
	sprintf(buf, "f64: %d cells wrong (%d called vertical)", n_wrong, n_f64_zero);
	svg_text(&sv, 880, 490, 12, "#cf222e", "start", buf);
	(void)wrong_colors;
	svg_close(&sv);
	printf("exact6_surface: %u x %u patches, %u bits; vertical %d, f64 wrong %d (f64 zero %d)\n", pu, pv, rep.storage_bits, n_zero, n_wrong,
		n_f64_zero);
	qaws_exact_surface_destroy(es);
	qaws_surface_destroy(surf);
}
