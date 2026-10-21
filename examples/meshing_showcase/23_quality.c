/*
 * Figure 23: quality reports (qaws_*_measure_elements,
 * qaws_quality_summarize). Above, three meshes coloured by a measure per
 * element, the ten worst ringed; below, each measure's histogram from the
 * summary, the solver's threshold in red and the count past it.
 *   - a 2D mesh of a jittered grid, by shape (mean ratio; threshold 0.5);
 *   - the blob of figure 2, by smallest angle (threshold 30 degrees);
 *   - the blob's tetrahedra cut open, by smallest dihedral (threshold 15).
 */

#define MS_Q_BINS 24

/* a histogram of the summary's bins over [lo, hi], the threshold marked */
static void ms_q_histogram(er_canvas* c, int px, int py, int pw, int ph, unsigned int const* hist, double lo, double hi, double threshold,
	char const* title, char const* sub)
{
	static double const bar[3] = { 90, 120, 170 }, axis[3] = { 60, 60, 60 }, red[3] = { 214, 39, 40 };
	double box[12] = { -0.03, 0, -0.05, 1.03, 0, -0.05, -0.03, 0, 0.62, 1.03, 0, 0.62 }, top = 1, p0[3], p1[3];
	unsigned int i;
	char lab[32];
	for (i = 0; i < MS_Q_BINS; i++)
		top = fmax(top, hist[i]);
	er_panel(c, px, py + 70, pw, ph - 110, box, 4, 0, 0, 0.05);
	er_text(c, px + 14, py + 12, title, 3, 30, 30, 30);
	er_text(c, px + 14, py + 46, sub, 2, 95, 95, 95);
	for (i = 0; i < MS_Q_BINS; i++)
		if (hist[i])
		{
			p0[0] = p1[0] = (i + 0.5) / MS_Q_BINS;
			p0[1] = p1[1] = 0;
			p0[2] = 0;
			p1[2] = 0.58 * hist[i] / top;
			er_line(c, p0, p1, bar, 9);
		}
	p0[0] = 0; p0[1] = 0; p0[2] = 0;
	p1[0] = 1; p1[1] = 0; p1[2] = 0;
	er_line(c, p0, p1, axis, 1.4);
	p0[0] = p1[0] = (threshold - lo) / (hi - lo);
	p1[2] = 0.6;
	er_line(c, p0, p1, red, 1.8);
	for (i = 0; i <= 4; i++)
	{
		double p[3] = { i / 4.0, 0, 0 }, s[3];
		er_project(c, p, s);
		sprintf(lab, "%g", lo + (hi - lo) * i / 4);
		er_text(c, (int)(s[0] / c->ss) - 10, (int)(s[1] / c->ss) + 8, lab, 2, 80, 80, 80);
	}
}

/* each face's centroid (3 per face) */
static double* ms_q_centroids(qaws_mesh_3d const* m)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), f, q, i, k;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	double *x = ms_positions(m), *cen = (double*)calloc(3 * (size_t)nf + 1, sizeof(double));
	qaws_mesh_3d_get_faces(fv, fs, m);
	for (f = 0, q = 0; f < nf; q += fs[f++])
		for (i = 0; i < fs[f]; i++)
			for (k = 0; k < 3; k++)
				cen[3 * f + k] += x[3 * fv[q + i] + k] / fs[f];
	free(fv); free(fs); free(x);
	return cen;
}

/* faces coloured by values mapped from [lo, hi] (bad to good), the worst ringed */
static void ms_q_faces(er_canvas* c, qaws_mesh_3d const* m, qaws_scalar const* v, double lo, double hi, unsigned int const* worst, unsigned int nw)
{
	static double const ring[3] = { 20, 20, 20 }, wire[3] = { 70, 70, 80 };
	unsigned int nf = qaws_mesh_3d_get_face_count(m), f;
	double *rgb = (double*)malloc(sizeof(double) * 3 * (nf + 1)), *cen = ms_q_centroids(m);
	for (f = 0; f < nf; f++)
		ms_quality_colour((v[f] - lo) / (hi - lo), rgb + 3 * f);
	ms_draw_faces(c, m, rgb);
	ms_draw(c, m, NULL, NULL, 0, wire, 0.6);
	for (f = 0; f < nw; f++)
		er_point(c, cen + 3 * worst[f], ring, 5);
	free(rgb);
	free(cen);
}

/* the kept tetrahedra's exposed faces: the border tinted, the cut coloured by the cell's value in [lo, hi]; the worst kept ones ringed */
static void ms_q_cutaway(er_canvas* c, qaws_volume_mesh const* v, qaws_scalar const* val, double lo, double hi, unsigned int const* worst,
	unsigned int nw)
{
	static unsigned char const F[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	static double const tint[3] = { 178, 192, 212 }, wire[3] = { 40, 45, 60 }, ring[3] = { 20, 20, 20 };
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k, j;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	unsigned char* kept = (unsigned char*)malloc(nc);
	double* cen = (double*)calloc(3 * (size_t)nc, sizeof(double));
	qaws_volume_mesh_get_positions(ps, v);
	for (i = 0; i < 3 * nv; i++)
		x[i] = (double)ps[i];
	qaws_volume_mesh_get_cells(cv, NULL, v);
	qaws_volume_mesh_get_cell_adjacency(adj, v);
	for (i = 0; i < nc; i++)
	{
		for (k = 0; k < 4; k++)
			for (j = 0; j < 3; j++)
				cen[3 * i + j] += x[3 * cv[4 * i + k] + j] / 4;
		kept[i] = (unsigned char)ms_keep_half(cen + 3 * i);
	}
	for (i = 0; i < nc; i++)
	{
		double col[3];
		if (!kept[i])
			continue;
		ms_quality_colour((val[i] - lo) / (hi - lo), col);
		for (k = 0; k < 4; k++)
		{
			unsigned int o = adj[4 * i + k];
			double const *a = x + 3 * cv[4 * i + F[k][0]], *b = x + 3 * cv[4 * i + F[k][1]], *d = x + 3 * cv[4 * i + F[k][2]];
			if (o != QAWS_MESH_NONE && kept[o])
				continue;
			er_triangle(c, a, b, d, NULL, NULL, NULL, o == QAWS_MESH_NONE ? tint : col, o == QAWS_MESH_NONE ? tint : col, o == QAWS_MESH_NONE ? tint : col);
			er_line(c, a, b, wire, 0.8);
			er_line(c, b, d, wire, 0.8);
			er_line(c, d, a, wire, 0.8);
		}
	}
	for (i = 0; i < nw; i++)
		if (kept[worst[i]])
			er_point(c, cen + 3 * worst[i], ring, 5);
	free(ps); free(x); free(cv); free(adj); free(kept); free(cen);
}

static void demo_quality(void)
{
	er_canvas* c = er_create(2200, 1150, 3, 252, 252, 250);
	unsigned int hist[MS_Q_BINS], worst[10], i, j, n = 14, s = 11;
	qaws_quality_summary sum;
	qaws_quality_summary_desc d;
	char sub[200];
	if (!c)
		return;
	memset(&d, 0, sizeof d);
	d.check = 1;
	d.bins = MS_Q_BINS;
	d.worst_count = 10;
	/* a 2D mesh: a jittered grid of triangles, drawn at z = 0 through its corner store */
	{
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * (n + 1) * (n + 1));
		unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n), t = 0;
		qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 2 * n * n);
		qaws_mesh_2d_desc md;
		qaws_mesh_2d* m = NULL;
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				int edge = i == 0 || j == 0 || i == n || j == n;
				x[2 * (j * (n + 1) + i)] = (qaws_scalar)(i + (edge ? 0 : 0.7 * (ms_uniform(&s) - 0.5)));
				x[2 * (j * (n + 1) + i) + 1] = (qaws_scalar)(j + (edge ? 0 : 0.7 * (ms_uniform(&s) - 0.5)));
			}
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int a = j * (n + 1) + i, b = a + 1, cc = a + n + 2, e = a + n + 1;
				f[t++] = a; f[t++] = b; f[t++] = cc;
				f[t++] = a; f[t++] = cc; f[t++] = e;
			}
		memset(&md, 0, sizeof md);
		md.positions = x;
		md.vertex_count = (n + 1) * (n + 1);
		md.face_vertices = f;
		md.face_count = 2 * n * n;
		qaws_mesh_2d_create(&m, &md);
		if (m)
		{
			qaws_mesh_3d const* m3 = qaws_mesh_2d_get_mesh(m);
			qaws_mesh_2d_measure_elements(v, m, QAWS_MEASURE_SHAPE);
			d.measure = QAWS_MEASURE_SHAPE;
			d.threshold = (qaws_scalar)0.5;
			d.range[0] = 0;
			d.range[1] = 1;
			qaws_quality_summarize(&sum, hist, worst, v, 2 * n * n, &d);
			sprintf(sub, "%u triangles, shape %.2f..%.2f, median %.2f", sum.count, (double)sum.min, (double)sum.max, (double)sum.median);
			ms_panel(c, 0, 0, 733, 700, m3, 0, 1.5707963, "a 2D mesh, by shape", sub);
			ms_q_faces(c, m3, v, 0, 1, worst, sum.worst_found);
			sprintf(sub, "%u under 0.5 (mean ratio)", sum.failing);
			ms_q_histogram(c, 0, 700, 733, 450, hist, 0, 1, 0.5, "shape", sub);
			printf("  2D: %s; %s\n", "shape", sub);
			qaws_mesh_2d_destroy(m);
		}
		free(x); free(f); free(v);
	}
	/* the blob by smallest angle, and its tetrahedra by smallest dihedral */
	{
		qaws_mesh_3d* b = ms_blob(3);
		qaws_scalar* v = NULL;
		qaws_volume_mesh* tv;
		if (!b)
			return;
		v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (qaws_mesh_3d_get_face_count(b) + 1));
		qaws_mesh_3d_measure_elements(v, b, QAWS_MEASURE_MIN_ANGLE);
		d.measure = QAWS_MEASURE_MIN_ANGLE;
		d.threshold = 30;
		d.range[0] = 0;
		d.range[1] = 60;
		qaws_quality_summarize(&sum, hist, worst, v, qaws_mesh_3d_get_face_count(b), &d);
		sprintf(sub, "%u triangles, smallest angle %.1f..%.1f", sum.count, (double)sum.min, (double)sum.max);
		ms_panel(c, 733, 0, 733, 700, b, 0.6, 0.45, "a surface, by smallest angle", sub);
		ms_q_faces(c, b, v, 0, 60, worst, sum.worst_found);
		sprintf(sub, "%u under 30 degrees, median %.1f", sum.failing, (double)sum.median);
		ms_q_histogram(c, 733, 700, 733, 450, hist, 0, 60, 30, "smallest angle", sub);
		printf("  surface: %s\n", sub);
		free(v);
		tv = ms_tets(b, sub, "blob");
		if (tv)
		{
			unsigned int nc = qaws_volume_mesh_get_cell_count(tv);
			v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (nc + 1));
			qaws_volume_mesh_measure_elements(v, tv, QAWS_MEASURE_MIN_ANGLE);
			d.threshold = 15;
			d.range[0] = 0;
			d.range[1] = 90;
			qaws_quality_summarize(&sum, hist, worst, v, nc, &d);
			sprintf(sub, "%u tetrahedra, dihedral %.1f..%.1f", sum.count, (double)sum.min, (double)sum.max);
			ms_vpanel(c, 1466, 0, 733, 700, tv, 0.6, 0.45, "a volume, by smallest dihedral", sub);
			ms_q_cutaway(c, tv, v, 0, 70, worst, sum.worst_found);
			sprintf(sub, "%u under 15 degrees, median %.1f", sum.failing, (double)sum.median);
			ms_q_histogram(c, 1466, 700, 733, 450, hist, 0, 90, 15, "smallest dihedral", sub);
			printf("  volume: %s\n", sub);
			free(v);
			qaws_volume_mesh_destroy(tv);
		}
		qaws_mesh_3d_destroy(b);
	}
	er_save(c, "showcase/mesh23_quality.png");
	printf("mesh23_quality.png\n");
	er_destroy(c);
}
