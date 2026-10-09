/* Figure 2: the heightfield of figure 1 certified: every contour / gradient
   crossing as an exact enclosure from one batched call, against the float
   batch and the pairwise exact call over every contour / gradient pair. */

static void demo_exact_heightfield(void)
{
	static double f[(HF_GRID + 1) * (HF_GRID + 1)];
	hf_scene sc;
	qaws_exact_curve** ec;
	qaws_exact_batch_hit* eh;
	qaws_curve_batch_hit_2d* fh = NULL;
	qaws_curve_batch_stats fst;
	qaws_exact_batch_stats est;
	qaws_exact_batch_desc d;
	qaws_exact_desc ed;
	qaws_exact_pair pbuf[256];
	unsigned int i, j, ne = 0, nf = 0, np = 0, cap = 1 << 15, quantized = 0, pairwise_failed = 0;
	double t0, t_exact, t_pair, t_float, t_dummy, widest = 0;
	qaws_status s;
	svg sv;
	viewport v, p;
	char sub[320], buf[200];
	double* xy = (double*)malloc(2 * 400 * sizeof(double));

	for (j = 0; j <= HF_GRID; j++)
		for (i = 0; i <= HF_GRID; i++)
			f[j * (HF_GRID + 1) + i] = hf_height(-1 + 2.0 * i / HF_GRID, -1 + 2.0 * j / HF_GRID, NULL, NULL);
	hf_build(&sc, 24, f);

	/* the float batch, for reference */
	hf_time(&sc, 0, &t_float, &t_dummy, &nf, &np, &fh, &fst);

	/* exact curves on the 2^-24 lattice */
	qaws_exact_desc_default(&ed);
	ed.space_exp2 = -24;
	ec = (qaws_exact_curve**)calloc(sc.count, sizeof(qaws_exact_curve*));
	for (i = 0; i < sc.count; i++)
	{
		qaws_exact_report rep;
		qaws_exact_curve_prepare(&ed, sc.curves[i], &ec[i], &rep);
		quantized += (rep.flags & QAWS_EXACT_FLAG_INPUT_QUANTIZED) != 0;
	}
	eh = (qaws_exact_batch_hit*)malloc(cap * sizeof(qaws_exact_batch_hit));
	memset(&d, 0, sizeof(d));
	d.curves = (qaws_exact_curve const* const*)ec;
	d.curve_count = sc.count;
	d.families = sc.family;
	t0 = hf_now();
	s = qaws_exact_curve_batch_hits(&d, eh, cap, &ne, &est);
	t_exact = hf_now() - t0;
	for (i = 0; i < ne && i < cap; i++)
	{
		double w = eh[i].pair.a_hi - eh[i].pair.a_lo;
		if (w > widest)
			widest = w;
	}

	/* pairwise exact over every contour / gradient pair */
	np = 0;
	t0 = hf_now();
	for (i = 0; i < sc.count; i++)
		for (j = i + 1; j < sc.count; j++)
			if (sc.family[i] != sc.family[j])
			{
				unsigned int k = 0;
				if (qaws_exact_curve_curve_hits(ec[i], ec[j], pbuf, 256, &k) != QAWS_STATUS_OK)
					pairwise_failed++;
				np += k;
			}
	t_pair = hf_now() - t0;
	printf("    exact heightfield: %u curves (%u quantized onto the lattice), %u spans, %u candidate span pairs; certified %u crossings in %.3f s (status %d, %u uncertified pairs, widest enclosure %.1e); float batch %u in %.4f s; pairwise exact %u in %.3f s (%u failed pairs)\n",
		sc.count, quantized, est.span_count, est.candidate_count, ne, t_exact, (int)s, est.uncertified_count, widest, nf, t_float, np, t_pair, pairwise_failed);

	sprintf(sub, "the %u curves of figure 1 as exact rational splines; %u crossings certified in one call (%.0f ms; pairwise exact %.0f ms)",
		sc.count, ne, t_exact * 1000, t_pair * 1000);
	if (!svg_open(&sv, "showcase/batch2_exact_heightfield.svg", 1400, 760, "Certified batched intersections", sub))
		return;
	v.x0 = 30; v.y0 = 80; v.w = 660; v.h = 660; v.xmin = -1; v.xmax = 1; v.ymin = -1; v.ymax = 1;
	svg_panel(&sv, &v, NULL);
	for (i = 0; i < sc.count; i++)
	{
		curve_polyline(sc.curves[i], &v, xy, 400);
		svg_polyline(&sv, xy, 400, sc.family[i] == 0 ? "#57606a" : "#54aeff", 1.0, 0.8, 0);
	}
	/* certified crossings: the point at the enclosure centre on the contour */
	for (i = 0; i < ne && i < cap; i++)
	{
		qaws_eval_result_2d r;
		qaws_exact_batch_hit const* h = &eh[i];
		qaws_curve_evaluate_2d(sc.curves[h->curve_a], (qaws_scalar)((h->pair.a_lo + h->pair.a_hi) / 2), QAWS_EVAL_FLAG_POSITION, &r);
		svg_circle(&sv, vx(&v, r.position.x), vy(&v, r.position.y), 2.6, h->pair.kind == QAWS_EXACT_HIT_POINT ? "#8250df" : "#1a7f37", "#ffffff");
	}
	svg_text(&sv, v.x0 + 10, v.y0 + 20, 13, "#24292f", "start", "green = certified crossing (isolating enclosure on both curves)");

	/* timings */
	p.x0 = 760; p.y0 = 80; p.w = 610; p.h = 300;
	svg_panel(&sv, &p, "time for every contour / gradient crossing");
	{
		char const* names[3] = { "float batch", "exact batch (certified)", "exact, pair by pair" };
		double ts[3];
		char const* cols[3] = { "#0969da", "#1a7f37", "#8250df" };
		double tmax;
		ts[0] = t_float; ts[1] = t_exact; ts[2] = t_pair;
		tmax = ts[2] > ts[1] ? ts[2] : ts[1];
		for (i = 0; i < 3; i++)
		{
			double y = p.y0 + 60 + i * 72, w = (p.w - 260) * ts[i] / tmax;
			fprintf(sv.f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"34\" fill=\"%s\"/>\n", p.x0 + 200, y, w < 2 ? 2 : w, cols[i]);
			svg_text(&sv, p.x0 + 190, y + 22, 13, "#24292f", "end", names[i]);
			sprintf(buf, "%.1f ms", ts[i] * 1000);
			svg_text(&sv, p.x0 + 208 + (w < 2 ? 2 : w), y + 22, 13, "#24292f", "start", buf);
		}
	}
	p.y0 = 410; p.h = 290;
	svg_panel(&sv, &p, "what the certified batch did");
	{
		char lines[7][160];
		sprintf(lines[0], "%u exact curves, %u rational Bezier spans", sc.count, est.span_count);
		sprintf(lines[1], "%u grid cells over sound control boxes of all spans", est.cell_count);
		sprintf(lines[2], "%u span pairs with overlapping boxes (of %.0f possible)", est.candidate_count, 0.5 * est.span_count * est.span_count);
		sprintf(lines[3], "%u crossings certified, %u curve pairs uncertified", ne, est.uncertified_count);
		sprintf(lines[4], "widest parameter enclosure %.1e", widest);
		sprintf(lines[5], "float batch on the original curves: %u crossings", nf);
		sprintf(lines[6], "inputs quantized onto the 2^-24 lattice: %u curves", quantized);
		for (i = 0; i < 7; i++)
			svg_text(&sv, p.x0 + 20, p.y0 + 50 + i * 30, 14, "#24292f", "start", lines[i]);
	}
	svg_close(&sv);
	printf("  -> showcase/batch2_exact_heightfield.svg\n");
	for (i = 0; i < sc.count; i++)
		qaws_exact_curve_destroy(ec[i]);
	free(ec);
	free(eh);
	free(fh);
	free(xy);
	hf_free(&sc);
}
