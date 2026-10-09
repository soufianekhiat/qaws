/* Figure 6: points around the terrain of figure 3 snapped to its nearest
   point, every point in one call on a prepared terrain set. */

#define SCL_POINTS 400

static void demo_surface_closest(void)
{
	qaws_surface* terrain = tr_terrain();
	qaws_surface_set* set = NULL;
	qaws_surface_batch_desc d;
	qaws_surface_batch_stats st;
	static qaws_scalar pts[3 * SCL_POINTS];
	static qaws_surface_batch_closest cp[SCL_POINTS];
	unsigned int i, reps = 0;
	unsigned long long state = 0xBB67AE8584CAA73Bull;
	double t0, t_query, zlo = 1e30, zhi = -1e30;
	svg s;
	projection pr;
	char sub[320];
	for (i = 0; i <= 1600; i++)
	{
		double z = tr_value(terrain, (qaws_scalar)((i % 41) / 40.0), (qaws_scalar)((i / 41) / 40.0));
		if (z < zlo) zlo = z;
		if (z > zhi) zhi = z;
	}
	for (i = 0; i < SCL_POINTS; i++)
	{
		double r[3];
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			state ^= state << 13; state ^= state >> 7; state ^= state << 17;
			r[k] = (double)(state >> 11) / 9007199254740992.0;
		}
		pts[3 * i] = (qaws_scalar)(-0.95 + 1.9 * r[0]);
		pts[3 * i + 1] = (qaws_scalar)(-0.95 + 1.9 * r[1]);
		pts[3 * i + 2] = (qaws_scalar)(zlo - 0.15 + (zhi - zlo + 0.25) * r[2]);
	}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_surface const* const*)&terrain;
	d.surface_count = 1;
	qaws_surface_set_create(&d, &set);
	t0 = hf_now();
	do
	{
		qaws_surface_set_find_closest(set, pts, SCL_POINTS, 0, cp, &st);
		reps++;
	} while (hf_now() - t0 < 0.05);
	t_query = (hf_now() - t0) / reps;
	printf("    surface closest: %u points snapped to the terrain (%u patches) in %.4f s (%u patches refined)\n", SCL_POINTS,
		qaws_surface_set_get_patch_count(set), t_query, st.candidate_count);
	sprintf(sub, "%u points above and below the terrain joined to their nearest terrain point, one call on a prepared set (%.1f ms; %u patch refinements over %u patches)",
		SCL_POINTS, t_query * 1000, st.candidate_count, qaws_surface_set_get_patch_count(set));
	if (!svg_open(&s, "showcase/batch6_surface_closest.svg", 1400, 760, "Batched closest points on a surface", sub))
		return;
	pr.cx = 700; pr.cy = 380; pr.scale = 300; pr.zscale = 1.0;
	draw_quads(&s, &pr, terrain, 50, 50, tr_value, terrain, zlo, zhi);
	for (i = 0; i < SCL_POINTS; i++)
	{
		double x0, y0, x1, y1;
		int above = pts[3 * i + 2] > cp[i].position.z;
		project(&pr, pts[3 * i], pts[3 * i + 1], pts[3 * i + 2], &x0, &y0);
		project(&pr, cp[i].position.x, cp[i].position.y, cp[i].position.z, &x1, &y1);
		svg_line(&s, x0, y0, x1, y1, above ? "#0969da" : "#8250df", 0.9, above ? 0.75 : 0.45);
		svg_circle(&s, x0, y0, 1.8, above ? "#0969da" : "#8250df", "none");
		svg_circle(&s, x1, y1, 1.6, "#cf222e", "#ffffff");
	}
	svg_text(&s, 700, 735, 14, "#24292f", "middle", "blue: points above the terrain, purple: below; red: their nearest terrain points");
	svg_close(&s);
	printf("  -> showcase/batch6_surface_closest.svg\n");
	qaws_surface_set_destroy(set);
	qaws_surface_destroy(terrain);
}
