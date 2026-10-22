/*
 * Figure 41: improvement of tetrahedral meshes (qaws_volume_mesh_improve),
 * cut open, cells coloured by their smallest dihedral angle (blue below 10
 * degrees to red at 70, as figure 6).
 *   - a cube's Kuhn grid, its inside vertices shaken into slivers;
 *   - the same improved: smoothing, then edge and face removals where it
 *     stalls, the border untouched;
 *   - the blob's tetrahedra (figure 6's operator) improved further;
 *   - the grid's smallest dihedral angles per cell, before and after.
 */

#define MS_VI_BINS 28

static int ms_vi_keep(double const* p)
{
	return p[1] > 0.5;
}

/* the unit cube's Kuhn grid of n^3 cubes, the inside vertices shaken by up to amp cells */
static qaws_volume_mesh* ms_vi_kuhn(unsigned int n, double amp)
{
	static unsigned int const perm[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
	unsigned int m1 = n + 1, nv = m1 * m1 * m1, nt = 6 * n * n * n, i, j, k, p, q = 0, seed = 777u;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nt);
	qaws_volume_mesh_desc d;
	qaws_volume_mesh* m = NULL;
	for (k = 0; k < m1; k++)
		for (j = 0; j < m1; j++)
			for (i = 0; i < m1; i++)
			{
				qaws_scalar* p3 = x + 3 * (i + m1 * (j + m1 * k));
				int inside = i > 0 && j > 0 && k > 0 && i < n && j < n && k < n;
				unsigned int c;
				p3[0] = (qaws_scalar)((double)i / n);
				p3[1] = (qaws_scalar)((double)j / n);
				p3[2] = (qaws_scalar)((double)k / n);
				for (c = 0; inside && c < 3; c++)
				{
					seed = seed * 1664525u + 1013904223u;
					p3[c] += (qaws_scalar)(amp / n * (((seed >> 8) & 0xFFFF) / 32768.0 - 1));
				}
			}
	for (k = 0; k < n; k++)
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
				for (p = 0; p < 6; p++)
				{
					unsigned int c[3] = { i, j, k }, v[4], s;
					v[0] = i + m1 * (j + m1 * k);
					for (s = 0; s < 3; s++)
					{
						c[perm[p][s]]++;
						v[s + 1] = c[0] + m1 * (c[1] + m1 * c[2]);
					}
					if (p == 1 || p == 2 || p == 5)   /* (the odd permutations' paths turn the other way) */
					{
						s = v[1]; v[1] = v[2]; v[2] = s;
					}
					memcpy(t + 4 * q++, v, sizeof v);
				}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.cell_vertices = t;
	d.cell_count = nt;
	qaws_volume_mesh_create(&m, &d);
	free(x);
	free(t);
	return m;
}

/* the tetrahedra's smallest dihedral angles binned over [0, 70] degrees (the last bin: 70 and more) */
static void ms_vi_bin(qaws_volume_mesh const* v, unsigned int* hist)
{
	static unsigned char const F[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k, j;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	memset(hist, 0, sizeof(unsigned int) * MS_VI_BINS);
	qaws_volume_mesh_get_positions(x, v);
	qaws_volume_mesh_get_cells(cv, NULL, v);
	for (i = 0; i < nc; i++)
	{
		double nrm[4][3], mind = 180;
		int b;
		for (k = 0; k < 4; k++)
		{
			qaws_scalar const *a = x + 3 * cv[4 * i + F[k][0]], *p = x + 3 * cv[4 * i + F[k][1]], *d = x + 3 * cv[4 * i + F[k][2]];
			double u[3] = { (double)p[0] - a[0], (double)p[1] - a[1], (double)p[2] - a[2] }, w[3] = { (double)d[0] - a[0], (double)d[1] - a[1], (double)d[2] - a[2] };
			nrm[k][0] = u[1] * w[2] - u[2] * w[1];
			nrm[k][1] = u[2] * w[0] - u[0] * w[2];
			nrm[k][2] = u[0] * w[1] - u[1] * w[0];
		}
		for (k = 0; k < 4; k++)
			for (j = k + 1; j < 4; j++)
			{
				double cs = -(nrm[k][0] * nrm[j][0] + nrm[k][1] * nrm[j][1] + nrm[k][2] * nrm[j][2]) /
					sqrt((nrm[k][0] * nrm[k][0] + nrm[k][1] * nrm[k][1] + nrm[k][2] * nrm[k][2]) * (nrm[j][0] * nrm[j][0] + nrm[j][1] * nrm[j][1] + nrm[j][2] * nrm[j][2]));
				mind = fmin(mind, acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180 / MS_PI);
			}
		b = (int)(mind / 70 * (MS_VI_BINS - 1));
		hist[b < 0 ? 0 : (b >= MS_VI_BINS ? MS_VI_BINS - 1 : b)]++;
	}
	free(x);
	free(cv);
}

/* two histograms side by side (before grey, after blue), counts on a square-root scale */
static void ms_vi_histogram(er_canvas* c, int px, int py, int pw, int ph, unsigned int const* h0, unsigned int const* h1, char const* title, char const* sub)
{
	static double const grey[3] = { 170, 170, 170 }, blue[3] = { 70, 110, 180 }, axis[3] = { 60, 60, 60 };
	double box[12] = { -0.03, 0, -0.05, 1.03, 0, -0.05, -0.03, 0, 0.62, 1.03, 0, 0.62 }, top = 1, p0[3], p1[3];
	unsigned int i, s;
	char lab[32];
	for (i = 0; i < MS_VI_BINS; i++) top = fmax(top, fmax(h0[i], h1[i]));
	er_panel(c, px, py + 70, pw, ph - 110, box, 4, 0, 0, 0.05);
	er_text(c, px + 14, py + 12, title, 3, 30, 30, 30);
	er_text(c, px + 14, py + 46, sub, 2, 95, 95, 95);
	for (s = 0; s < 2; s++)
		for (i = 0; i < MS_VI_BINS; i++)
		{
			unsigned int const* h = s ? h1 : h0;
			if (!h[i])
				continue;
			p0[0] = p1[0] = (i + 0.3 + 0.4 * s) / MS_VI_BINS;
			p0[1] = p1[1] = 0;
			p0[2] = 0;
			p1[2] = 0.58 * sqrt(h[i] / top);
			er_line(c, p0, p1, s ? blue : grey, 4);
		}
	p0[0] = 0; p0[1] = 0; p0[2] = 0;
	p1[0] = 1; p1[1] = 0; p1[2] = 0;
	er_line(c, p0, p1, axis, 1.4);
	for (i = 0; i <= 7; i++)
	{
		double p[3] = { i / 7.0 * (MS_VI_BINS - 1) / MS_VI_BINS, 0, 0 }, sc[3];
		er_project(c, p, sc);
		sprintf(lab, "%u", 10 * i);
		er_text(c, (int)(sc[0] / c->ss) - 8, (int)(sc[1] / c->ss) + 8, lab, 2, 80, 80, 80);
	}
}

static void demo_volume_improve(void)
{
	static double const wire[3] = { 40, 45, 60 }, tan[3] = { 232, 196, 120 }, tint[3] = { 178, 192, 212 };
	er_canvas* c = er_create(1800, 700, 3, 252, 252, 250);
	qaws_volume_mesh *g = ms_vi_kuhn(6, 0.32), *gi = NULL, *b = NULL, *bi = NULL;
	qaws_volume_improve_report r;
	qaws_volume_quality q;
	unsigned int h0[MS_VI_BINS], h1[MS_VI_BINS];
	char sub[200];
	if (!c || !g)
		return;
	qaws_volume_mesh_measure_quality(&q, g);
	sprintf(sub, "%u cells, dihedral %.1f..%.1f", q.tetrahedra, (double)q.min_dihedral, (double)q.max_dihedral);
	ms_vpanel(c, 0, 0, 450, 700, g, 0.6, 0.5, "slivers", sub);
	ms_cutaway(c, g, ms_vi_keep, tan, wire);
	printf("  shaken grid: %s\n", sub);
	qaws_volume_mesh_improve(&gi, &r, g, NULL);
	if (gi)
	{
		sprintf(sub, "dihedral %.1f..%.1f, %u flips", (double)r.quality.min_dihedral, (double)r.quality.max_dihedral,
			r.edge_removals + r.face_removals);
		ms_vpanel(c, 450, 0, 450, 700, gi, 0.6, 0.5, "improved", sub);
		ms_cutaway(c, gi, ms_vi_keep, tan, wire);
		printf("  improved: %s\n", sub);
		ms_vi_bin(g, h0);
		ms_vi_bin(gi, h1);
		ms_vi_histogram(c, 1350, 0, 450, 700, h0, h1, "smallest dihedral", "per cell, degrees: before, after");
	}
	{
		qaws_mesh_3d* blob = ms_blob(2);
		qaws_tetrahedralize_desc td;
		memset(&td, 0, sizeof td);
		if (blob)
			qaws_volume_mesh_tetrahedralize(&b, NULL, blob, &td);
		if (b)
		{
			qaws_volume_mesh_measure_quality(&q, b);
			qaws_volume_mesh_improve(&bi, &r, b, NULL);
		}
		if (bi)
		{
			sprintf(sub, "dihedral %.1f -> %.1f, %u flips", (double)q.min_dihedral, (double)r.quality.min_dihedral, r.edge_removals + r.face_removals);
			ms_vpanel(c, 900, 0, 450, 700, bi, 0.4, 0.6, "a blob's cells", sub);
			ms_cutaway(c, bi, ms_keep_half, tint, wire);
			printf("  blob: %s, %u moves\n", sub, r.moves);
		}
		qaws_mesh_3d_destroy(blob);
	}
	qaws_volume_mesh_destroy(g);
	qaws_volume_mesh_destroy(gi);
	qaws_volume_mesh_destroy(b);
	qaws_volume_mesh_destroy(bi);
	er_save(c, "showcase/mesh41_volume_improve.png");
	printf("mesh41_volume_improve.png\n");
	er_destroy(c);
}
