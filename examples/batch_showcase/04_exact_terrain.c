/* Figure 4: the contours of figure 3 certified: the terrain and the planes
   as exact surfaces, every contour branch from one batched call, its
   certified points on the float contour polylines. */

static void demo_exact_terrain(void)
{
	qaws_surface* sf[1 + TR_LEVELS];
	qaws_exact_surface* es[1 + TR_LEVELS];
	unsigned int fam[1 + TR_LEVELS], i, j, nc = 0, np = 0, enp = 0, enb = 0, agree = 0;
	static qaws_surface_batch_curve cv[1024];
	static qaws_ssi_point pt[1 << 17];
	static qaws_exact_ssi_point ept[1 << 16];
	static qaws_exact_ssi_batch_branch ebr[1024];
	qaws_surface_batch_desc sd;
	qaws_exact_ssi_batch_desc ed;
	qaws_exact_batch_stats est;
	qaws_exact_desc xd;
	qaws_status st;
	double t0, t_float, t_exact, zlo = 1e30, zhi = -1e30, widest = 0;
	svg s;
	projection pr;
	char sub[400], buf[200];

	sf[0] = tr_terrain();
	fam[0] = 0;
	for (i = 0; i <= 40; i++)
		for (j = 0; j <= 40; j++)
		{
			double z = tr_value(sf[0], (qaws_scalar)(i / 40.0), (qaws_scalar)(j / 40.0));
			if (z < zlo) zlo = z;
			if (z > zhi) zhi = z;
		}
	for (i = 0; i < TR_LEVELS; i++)
	{
		/* dyadic levels: the planes are exact as given */
		double z = ldexp(nearbyint(ldexp(zlo + (zhi - zlo) * (i + 0.5) / TR_LEVELS, 12)), -12);
		sf[1 + i] = tr_plane(z);
		fam[1 + i] = 1;
	}
	memset(&sd, 0, sizeof(sd));
	sd.surfaces = (qaws_surface const* const*)sf;
	sd.surface_count = 1 + TR_LEVELS;
	sd.families = fam;
	t0 = hf_now();
	qaws_surface_batch_find_intersections(&sd, cv, 1024, &nc, pt, 1 << 17, &np, NULL);
	t_float = hf_now() - t0;

	qaws_exact_desc_default(&xd);
	xd.space_exp2 = -20;
	for (i = 0; i <= TR_LEVELS; i++)
		qaws_exact_surface_prepare(&xd, sf[i], &es[i], NULL);
	memset(&ed, 0, sizeof(ed));
	ed.surfaces = (qaws_exact_surface const* const*)es;
	ed.surface_count = 1 + TR_LEVELS;
	ed.families = fam;
	t0 = hf_now();
	st = qaws_exact_surface_batch_hits(&ed, ept, 1 << 16, &enp, ebr, 1024, &enb, &est);
	t_exact = hf_now() - t0;
	for (i = 0; i < enp; i++)
		if (ept[i].u1_hi - ept[i].u1_lo > widest)
			widest = ept[i].u1_hi - ept[i].u1_lo;
	printf("    exact terrain: %u patches, %u candidate patch pairs; %u certified branches, %u points in %.2f s (status %d, %u uncertified, widest u enclosure %.1e); float %u curves in %.3f s\n",
		est.patch_count, est.candidate_count, enb, enp, t_exact, (int)st, est.uncertified_count, widest, nc, t_float);

	/* per level: as many float curves as certified branches? */
	for (i = 0; i < TR_LEVELS; i++)
	{
		unsigned int fc = 0, ec = 0, k;
		for (k = 0; k < nc && k < 1024; k++)
			fc += cv[k].surface_b == 1 + i;
		for (k = 0; k < enb; k++)
			ec += ebr[k].surface_b == 1 + i;
		agree += fc == ec;
	}
	printf("    float curve and certified branch counts agree on %u of %u levels\n", agree, TR_LEVELS);
	sprintf(sub, "terrain x %u planes as exact surfaces: %u certified contour branches, %u certified points in one call (%.2f s; %u candidate patch pairs over %u patches); float batch %.0f ms, same branch count on %u of %u levels",
		TR_LEVELS, enb, enp, t_exact, est.candidate_count, est.patch_count, t_float * 1000, agree, TR_LEVELS);
	if (!svg_open(&s, "showcase/batch4_exact_terrain.svg", 1400, 760, "Certified batched surface intersections", sub))
		return;
	pr.cx = 700; pr.cy = 380; pr.scale = 300; pr.zscale = 1.0;
	draw_quads(&s, &pr, sf[0], 50, 50, tr_value, sf[0], zlo, zhi);
	for (i = 0; i < nc && i < 1024; i++)
	{
		double* xy = (double*)malloc(2 * (cv[i].count + 1) * sizeof(double));
		unsigned int k;
		for (k = 0; k < cv[i].count; k++)
		{
			qaws_ssi_point const* p = &pt[cv[i].first + k];
			project(&pr, p->position.x, p->position.y, p->position.z + 0.004, &xy[2 * k], &xy[2 * k + 1]);
		}
		if (cv[i].closed)
		{
			xy[2 * k] = xy[0];
			xy[2 * k + 1] = xy[1];
			k++;
		}
		svg_polyline(&s, xy, (int)k, "#1b1f24", 1.0, 0.7, 0);
		free(xy);
	}
	/* certified points: the terrain at the centre of each (u1, v1) enclosure */
	for (i = 0; i < enb; i++)
	{
		unsigned int k;
		for (k = 0; k < ebr[i].branch.count; k++)
		{
			qaws_exact_ssi_point const* p = &ept[ebr[i].branch.first + k];
			qaws_surface_eval_result r;
			double x, y;
			qaws_surface_evaluate(sf[0], (qaws_scalar)((p->u1_lo + p->u1_hi) / 2), (qaws_scalar)((p->v1_lo + p->v1_hi) / 2), QAWS_SURFACE_EVAL_POSITION, &r);
			project(&pr, r.position.x, r.position.y, r.position.z + 0.004, &x, &y);
			svg_circle(&s, x, y, 2.4, ebr[i].branch.closed ? "#1a7f37" : "#8250df", "#ffffff");
		}
	}
	sprintf(buf, "black: float contour polylines; dots: certified points (green on closed branches, purple on open ones); widest enclosure %.0e", widest);
	svg_text(&s, 700, 735, 14, "#24292f", "middle", buf);
	svg_close(&s);
	printf("  -> showcase/batch4_exact_terrain.svg\n");
	for (i = 0; i <= TR_LEVELS; i++)
	{
		qaws_exact_surface_destroy(es[i]);
		qaws_surface_destroy(sf[i]);
	}
}
