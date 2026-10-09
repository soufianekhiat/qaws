/* Figure 5: points snapped to the nearest contour line of the heightfield
   of figure 1, every point in one call on a prepared contour set. */

#define CL_POINTS 800

static void demo_closest(void)
{
	static double f[(HF_GRID + 1) * (HF_GRID + 1)];
	static qaws_scalar pts[2 * CL_POINTS];
	static qaws_closest_point cp[CL_POINTS];
	hf_scene sc;
	qaws_curve const** contours;
	qaws_curve_set* set = NULL;
	qaws_curve_batch_desc d;
	qaws_curve_batch_stats st;
	unsigned int i, j, n = 0, reps = 0;
	unsigned long long state = 0x6A09E667F3BCC909ull;
	double t0, t_query, t_set;
	svg s;
	viewport v;
	char sub[320];
	double* xy = (double*)malloc(2 * 400 * sizeof(double));

	for (j = 0; j <= HF_GRID; j++)
		for (i = 0; i <= HF_GRID; i++)
			f[j * (HF_GRID + 1) + i] = hf_height(-1 + 2.0 * i / HF_GRID, -1 + 2.0 * j / HF_GRID, NULL, NULL);
	hf_build(&sc, 24, f);
	contours = (qaws_curve const**)malloc(sc.count * sizeof(qaws_curve*));
	for (i = 0; i < sc.count; i++)
		if (sc.family[i] == 0)
			contours[n++] = sc.curves[i];
	for (i = 0; i < CL_POINTS; i++)
	{
		double r[2];
		unsigned int k;
		for (k = 0; k < 2; k++)
		{
			state ^= state << 13; state ^= state >> 7; state ^= state << 17;
			r[k] = (double)(state >> 11) / 9007199254740992.0;
		}
		pts[2 * i] = (qaws_scalar)(-0.98 + 1.96 * r[0]);
		pts[2 * i + 1] = (qaws_scalar)(-0.98 + 1.96 * r[1]);
	}
	/* the contours prepared once, then the queries */
	memset(&d, 0, sizeof(d));
	d.curves = contours;
	d.curve_count = n;
	t0 = hf_now();
	qaws_curve_set_create(&d, &set);
	t_set = hf_now() - t0;
	t0 = hf_now();
	do
	{
		qaws_curve_set_find_closest(set, pts, CL_POINTS, 0, cp, &st);
		reps++;
	} while (hf_now() - t0 < 0.05);
	t_query = (hf_now() - t0) / reps;
	printf("    closest: %u points snapped to %u contours (%u segments) in %.4f s (set prepared in %.4f s, %u segments refined)\n", CL_POINTS, n,
		qaws_curve_set_get_segment_count(set), t_query, t_set, st.candidate_count);

	sprintf(sub, "%u points snapped to the nearest of %u contour lines in one call on a prepared set (%.1f ms, %u candidate segments refined of %u)",
		CL_POINTS, n, t_query * 1000, st.candidate_count, qaws_curve_set_get_segment_count(set));
	if (!svg_open(&s, "showcase/batch5_closest.svg", 1400, 760, "Batched closest points", sub))
		return;
	v.x0 = 370; v.y0 = 80; v.w = 660; v.h = 660; v.xmin = -1; v.xmax = 1; v.ymin = -1; v.ymax = 1;
	svg_panel(&s, &v, NULL);
	for (i = 0; i < n; i++)
	{
		curve_polyline(contours[i], &v, xy, 400);
		svg_polyline(&s, xy, 400, "#24292f", 1.1, 0.8, 0);
	}
	for (i = 0; i < CL_POINTS; i++)
	{
		char col[16];
		double t = cp[i].distance / 0.08;
		if (cp[i].curve == QAWS_CURVE_BATCH_NONE)
			continue;
		heat(t > 1 ? 1 : t, col);
		svg_line(&s, vx(&v, pts[2 * i]), vy(&v, pts[2 * i + 1]), vx(&v, cp[i].position.x), vy(&v, cp[i].position.y), col, 1.1, 0.9);
		svg_circle(&s, vx(&v, pts[2 * i]), vy(&v, pts[2 * i + 1]), 1.8, col, col);
	}
	svg_text(&s, v.x0 + 10, v.y0 + 20, 13, "#24292f", "start", "each point joined to its nearest contour point, colored by distance");
	svg_close(&s);
	printf("  -> showcase/batch5_closest.svg\n");
	qaws_curve_set_destroy(set);
	free((void*)contours);
	free(xy);
	hf_free(&sc);
}
