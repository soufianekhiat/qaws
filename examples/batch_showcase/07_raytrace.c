/* Figure 7: the terrain of figure 3 ray traced by the batch: one call for
   the camera rays, one for the shadow rays toward the light. */

#define RT_W 240
#define RT_H 150

static void rt_normalize(double* v)
{
	double l = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] /= l; v[1] /= l; v[2] /= l;
}

static void demo_raytrace(void)
{
	qaws_surface* terrain = tr_terrain();
	qaws_surface_set* set = NULL;
	qaws_surface_batch_desc d;
	qaws_surface_batch_stats st1, st2;
	static qaws_scalar org[3 * RT_W * RT_H], dir[3 * RT_W * RT_H], sorg[3 * RT_W * RT_H], sdir[3 * RT_W * RT_H];
	static qaws_surface_ray_hit hit[RT_W * RT_H], shadow[RT_W * RT_H];
	double eye[3] = { 2.3, -2.7, 2.1 }, at[3] = { 0, 0, -0.1 }, up[3] = { 0, 0, 1 }, fw[3], rt[3], uv[3], light[3] = { -0.85, -0.45, 0.32 };
	double t0, t_primary, t_shadow, zlo = 1e30, zhi = -1e30;
	unsigned int i, j, k, nhit = 0, nshadow = 0;
	svg s;
	char sub[320];
	for (i = 0; i <= 1600; i++)
	{
		double z = tr_value(terrain, (qaws_scalar)((i % 41) / 40.0), (qaws_scalar)((i / 41) / 40.0));
		if (z < zlo) zlo = z;
		if (z > zhi) zhi = z;
	}
	rt_normalize(light);
	for (k = 0; k < 3; k++)
		fw[k] = at[k] - eye[k];
	rt_normalize(fw);
	rt[0] = fw[1] * up[2] - fw[2] * up[1]; rt[1] = fw[2] * up[0] - fw[0] * up[2]; rt[2] = fw[0] * up[1] - fw[1] * up[0];
	rt_normalize(rt);
	uv[0] = rt[1] * fw[2] - rt[2] * fw[1]; uv[1] = rt[2] * fw[0] - rt[0] * fw[2]; uv[2] = rt[0] * fw[1] - rt[1] * fw[0];
	for (j = 0; j < RT_H; j++)
		for (i = 0; i < RT_W; i++)
		{
			double sx = (i + 0.5) / RT_W * 2 - 1, sy = 1 - (j + 0.5) / RT_H * 2, a = 0.4;
			unsigned int p = j * RT_W + i;
			for (k = 0; k < 3; k++)
			{
				org[3 * p + k] = (qaws_scalar)eye[k];
				dir[3 * p + k] = (qaws_scalar)(fw[k] + a * sx * rt[k] + a * sy * (double)RT_H / RT_W * uv[k]);
			}
		}
	memset(&d, 0, sizeof(d));
	d.surfaces = (qaws_surface const* const*)&terrain;
	d.surface_count = 1;
	qaws_surface_set_create(&d, &set);
	t0 = hf_now();
	qaws_surface_set_raycast(set, org, dir, RT_W * RT_H, 0, hit, &st1);
	t_primary = hf_now() - t0;
	/* shadow rays from every hit, slightly lifted, toward the light */
	for (i = 0; i < RT_W * RT_H; i++)
	{
		if (hit[i].surface == QAWS_CURVE_BATCH_NONE)
		{
			sdir[3 * i] = sdir[3 * i + 1] = sdir[3 * i + 2] = 0;
			continue;
		}
		nhit++;
		sorg[3 * i] = (qaws_scalar)(hit[i].position.x + 1e-4 * light[0]);
		sorg[3 * i + 1] = (qaws_scalar)(hit[i].position.y + 1e-4 * light[1]);
		sorg[3 * i + 2] = (qaws_scalar)(hit[i].position.z + 1e-4 * light[2]);
		for (k = 0; k < 3; k++)
			sdir[3 * i + k] = (qaws_scalar)light[k];
	}
	t0 = hf_now();
	qaws_surface_set_raycast(set, sorg, sdir, RT_W * RT_H, 0, shadow, &st2);
	t_shadow = hf_now() - t0;
	for (i = 0; i < RT_W * RT_H; i++)
		nshadow += hit[i].surface != QAWS_CURVE_BATCH_NONE && shadow[i].surface != QAWS_CURVE_BATCH_NONE;
	printf("    ray trace: %u x %u camera rays in %.3f s (%u hits), %u shadow rays in %.3f s (%u in shadow)\n", RT_W, RT_H, t_primary, nhit, nhit, t_shadow, nshadow);

	sprintf(sub, "the B-spline terrain ray traced by the batch: %u camera rays (%.0f ms) and %u shadow rays (%.0f ms), shading from the exact surface normal",
		RT_W * RT_H, t_primary * 1000, nhit, t_shadow * 1000);
	if (!svg_open(&s, "showcase/batch7_raytrace.svg", 1400, 760, "Batched ray casting", sub))
		return;
	{
		double px = 1320.0 / RT_W, x0 = 40, y0 = 80;
		for (j = 0; j < RT_H; j++)
			for (i = 0; i < RT_W; i++)
			{
				unsigned int p = j * RT_W + i;
				char col[32];
				if (hit[p].surface == QAWS_CURVE_BATCH_NONE)
				{
					/* sky */
					double t = (double)j / RT_H;
					sprintf(col, "#%02x%02x%02x", (int)(205 + 30 * t), (int)(225 + 20 * t), 250);
				}
				else
				{
					qaws_surface_eval_result r;
					double n[3], lam, hc;
					int rr = 0, gg = 0, bb = 0;
					qaws_surface_evaluate(terrain, hit[p].u, hit[p].v, QAWS_SURFACE_EVAL_NORMAL, &r);
					n[0] = r.normal.x; n[1] = r.normal.y; n[2] = r.normal.z;
					if (n[2] < 0) { n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2]; }
					lam = n[0] * light[0] + n[1] * light[1] + n[2] * light[2];
					if (lam < 0) lam = 0;
					if (shadow[p].surface != QAWS_CURVE_BATCH_NONE)
						lam *= 0.35;
					hc = (hit[p].position.z - zlo) / (zhi - zlo);
					heat(hc < 0 ? 0 : (hc > 1 ? 1 : hc), col);
					sscanf(col, "rgb(%d,%d,%d)", &rr, &gg, &bb);
					rr = (int)(rr * (0.42 + 0.58 * lam));
					gg = (int)(gg * (0.42 + 0.58 * lam));
					bb = (int)(bb * (0.42 + 0.58 * lam));
					sprintf(col, "#%02x%02x%02x", rr, gg, bb);
				}
				fprintf(s.f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\" shape-rendering=\"crispEdges\"/>\n", x0 + i * px, y0 + j * px, px + 0.6,
					px + 0.6, col);
			}
	}
	svg_close(&s);
	printf("  -> showcase/batch7_raytrace.svg\n");
	qaws_surface_set_destroy(set);
	qaws_surface_destroy(terrain);
}
