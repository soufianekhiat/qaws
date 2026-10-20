/*
 * Figure 17: sizing by local feature size.
 *   - a block with a fin 0.02 thick, remeshed by LOCAL at 0.05 (angles
 *     bounded to 30..110) with uniform sizes and with the feature-size field
 *     at 3 cells across (gradation 0.3);
 *   - a clip (a C, its arms 0.1 thick, the slot between them 0.05 wide),
 *     the same two ways;
 *   - each vertex coloured by the mean length of its edges, on a log scale
 *     from 0.005 (blue) to 0.05 (red); printed: the vertex counts, the
 *     widest triangle on the fin's narrow sides, the field's time.
 */

/* a profile in (x, z), counter-clockwise, extruded over y in [-w, w]; its caps as convex polygons of the profile's points (sizes
   in cs, indices in ci) */
static qaws_mesh_3d* ms_extrude(double const* pr, unsigned int np, double w, unsigned int const* cs, unsigned int nc, unsigned int const* ci)
{
	double x[3 * 2 * 16];
	unsigned int f[4 * 16 + 2 * 32], fs[16 + 2 * 8], nf = 0, n = 0, k, c, i, at;
	for (k = 0; k < 2 * np; k++)
	{
		x[3 * k] = pr[2 * (k % np)];
		x[3 * k + 1] = k < np ? -w : w;
		x[3 * k + 2] = pr[2 * (k % np) + 1];
	}
	for (k = 0; k < np; k++)
	{
		f[n++] = k; f[n++] = np + k; f[n++] = np + (k + 1) % np; f[n++] = (k + 1) % np;
		fs[nf++] = 4;
	}
	for (c = 0, at = 0; c < nc; at += cs[c++])
	{
		for (i = 0; i < cs[c]; i++)
			f[n++] = ci[at + i];
		fs[nf++] = cs[c];
		for (i = 0; i < cs[c]; i++)
			f[n++] = np + ci[at + cs[c] - 1 - i];
		fs[nf++] = cs[c];
	}
	return ms_make(x, 2 * np, f, fs, nf);
}

/* each vertex coloured by its edges' mean length, log scale on [lo, hi] */
static double* ms_length_colours(qaws_mesh_3d const* m, double lo, double hi)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), ne = qaws_mesh_3d_get_edge_count(m), i;
	unsigned int *e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1)), *cnt = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
	double *x = ms_positions(m), *sum = (double*)calloc(nv + 1, sizeof(double)), *rgb = (double*)malloc(sizeof(double) * 3 * (nv + 1));
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		double const *a = x + 3 * e[2 * i], *b = x + 3 * e[2 * i + 1];
		double l = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
		sum[e[2 * i]] += l;
		sum[e[2 * i + 1]] += l;
		cnt[e[2 * i]]++;
		cnt[e[2 * i + 1]]++;
	}
	for (i = 0; i < nv; i++)
		ms_colormap(cnt[i] ? log(sum[i] / cnt[i] / lo) / log(hi / lo) : 0, rgb + 3 * i);
	free(e); free(cnt); free(x); free(sum);
	return rgb;
}

/* the widest z extent of the triangles on the fin's narrow sides (|z| < 0.01, y = +-0.15 or x = 0.15 past the block) */
static double ms_fin_span(qaws_mesh_3d const* m)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, k;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	double *x = ms_positions(m), worst = 0;
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		double c[3] = { 0, 0, 0 }, zl = HUGE_VAL, zh = -HUGE_VAL;
		for (k = 0; k < 3; k++)
		{
			double const* v = x + 3 * t[3 * i + k];
			c[0] += v[0] / 3; c[1] += v[1] / 3; c[2] += v[2] / 3;
			zl = fmin(zl, v[2]);
			zh = fmax(zh, v[2]);
		}
		if (c[0] > -0.2 && fabs(c[2]) < 0.0101 && (fabs(fabs(c[1]) - 0.15) < 1e-4 || c[0] > 0.15 - 1e-4))
			worst = fmax(worst, zh - zl);
	}
	free(t); free(x);
	return worst;
}

static void demo_feature_size(void)
{
	static double const fin[16] = { -0.5, -0.2, -0.2, -0.2, -0.2, -0.01, 0.15, -0.01, 0.15, 0.01, -0.2, 0.01, -0.2, 0.2, -0.5, 0.2 };
	static unsigned int const fin_cs[2] = { 6, 4 }, fin_ci[10] = { 0, 1, 2, 5, 6, 7, 2, 3, 4, 5 };
	static double const clip[20] = { -0.3, -0.125, -0.2, -0.125, 0.3, -0.125, 0.3, -0.025, -0.2, -0.025, -0.2, 0.025, 0.3, 0.025, 0.3, 0.125, -0.2, 0.125,
		-0.3, 0.125 };
	static unsigned int const clip_cs[3] = { 6, 4, 4 }, clip_ci[14] = { 0, 1, 4, 5, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8 };
	static double const wire[3] = { 40, 45, 60 };
	static double const view[2][2] = { { 0.75, 0.55 }, { 1.05, 0.3 } };
	er_canvas* c = er_create(2200, 640, 3, 252, 252, 250);
	qaws_mesh_3d* shape[2];
	unsigned int s;
	if (!c)
		return;
	shape[0] = ms_extrude(fin, 8, 0.15, fin_cs, 2, fin_ci);
	shape[1] = ms_extrude(clip, 10, 0.15, clip_cs, 3, clip_ci);
	for (s = 0; s < 2; s++)
	{
		qaws_mesh_reference* r = NULL;
		qaws_size_field* sf = NULL;
		qaws_feature_size_desc fd;
		qaws_mesh_report rep;
		qaws_mesh_3d* out[2];
		clock_t t0;
		double build;
		unsigned int k;
		if (!shape[s])
			continue;
		qaws_mesh_3d_detect_features(shape[s], (qaws_scalar)(MS_PI / 6), 0, 0);
		qaws_mesh_reference_create(&r, shape[s], NULL);
		memset(&fd, 0, sizeof fd);
		fd.cells_across = 3;
		fd.min_length = (qaws_scalar)0.002;
		fd.max_length = (qaws_scalar)0.05;
		fd.gradation = (qaws_scalar)0.3;
		t0 = clock();
		if (r)
			qaws_size_field_create_feature_size(&sf, r, &fd);
		build = (double)(clock() - t0) / CLOCKS_PER_SEC;
		for (k = 0; k < 2; k++)
		{
			char title[100], sub[200];
			double* rgb;
			out[k] = r ? ms_remesh_bounded(shape[s], r, 0.05, k ? sf : NULL, 30, 110, 0, &rep) : NULL;
			if (!out[k])
				continue;
			sprintf(title, "%s, %s", s ? "a clip" : "a fin", k ? "feature size" : "uniform 0.05");
			if (s == 0)
				sprintf(sub, "%u vertices, thin sides %.4f / 0.02", qaws_mesh_3d_get_vertex_count(out[k]), ms_fin_span(out[k]));
			else
				sprintf(sub, "%u vertices, angles %.0f..%.0f", qaws_mesh_3d_get_vertex_count(out[k]), (double)rep.quality.min_angle,
					(double)rep.quality.max_angle);
			ms_panel(c, 1100 * (int)s + 550 * (int)k, 0, 550, 640, out[k], view[s][0], view[s][1], title, sub);
			if (k)
			{
				sprintf(sub, "3 cells across, field in %.2f s", build);
				er_text(c, 1100 * (int)s + 550 + 14, 68, sub, 2, 95, 95, 95);
			}
			rgb = ms_length_colours(out[k], 0.005, 0.05);
			ms_draw(c, out[k], NULL, rgb, 0, wire, 0.5);
			printf("  %s: %u vertices, angles %.1f..%.1f, mean edge %.4f%s\n", title, qaws_mesh_3d_get_vertex_count(out[k]), (double)rep.quality.min_angle,
				(double)rep.quality.max_angle, (double)rep.quality.mean_edge, k ? "" : " (uniform)");
			if (s == 0)
				printf("    the fin's narrow sides: the widest triangle spans %.4f of 0.02\n", ms_fin_span(out[k]));
			free(rgb);
		}
		if (sf)
			printf("  %s: the field built in %.2f s\n", s ? "the clip" : "the fin", build);
		qaws_mesh_3d_destroy(out[0]);
		qaws_mesh_3d_destroy(out[1]);
		qaws_size_field_destroy(sf);
		qaws_mesh_reference_destroy(r);
		qaws_mesh_3d_destroy(shape[s]);
	}
	er_text(c, 14, 610, "edge length: 0.005 blue .. 0.05 red (log)", 2, 95, 95, 95);
	er_save(c, "showcase/mesh17_feature_size.png");
	er_destroy(c);
}
