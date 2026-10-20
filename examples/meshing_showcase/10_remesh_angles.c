/*
 * Figure 10: remeshing with bounded angles.
 *   - the blob, the adaptive bumpy sphere and the can of figure 3, remeshed
 *     by LOCAL's defaults (top) and with the angle-bounded pass at Wang et
 *     al.'s [35, 86] degrees (middle), each triangle coloured by its smallest
 *     angle (blue at 60, through white to dark red at 30, red under 30);
 *   - the corner angles of both, as histograms (grey: default, orange:
 *     bounded), the bounds drawn;
 *   - printed: the same per panel, and the blob with a tolerance of 2e-3
 *     (its vertices, its sampled distance to the reference).
 */

static qaws_mesh_3d* ms_remesh_bounded(qaws_mesh_3d const* m, qaws_mesh_reference const* r, double L, qaws_size_field const* sf, double lo,
	double hi, double tol, qaws_mesh_report* rep)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	d.size = sf;
	d.min_angle = (qaws_scalar)(lo * MS_PI / 180);
	d.max_angle = (qaws_scalar)(hi * MS_PI / 180);
	d.tolerance = (qaws_scalar)tol;
	if (qaws_mesh_remesh(&out, rep, m, &d) != QAWS_STATUS_OK)
		printf("  remesh failed\n");
	return out;
}

/* each corner's angle in degrees (3 per triangle, caller frees) */
static double* ms_corner_angles(qaws_mesh_3d const* m, unsigned int* n)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, k;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt ? nt : 1));
	double *x = ms_positions(m), *a = (double*)malloc(sizeof(double) * 3 * (nt ? nt : 1));
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < 3 * nt; i++)
	{
		double const *P = x + 3 * t[i], *A = x + 3 * t[3 * (i / 3) + (i + 1) % 3], *B = x + 3 * t[3 * (i / 3) + (i + 2) % 3];
		double u[3], v[3], uu = 0, vv = 0, uv = 0, c;
		for (k = 0; k < 3; k++)
		{
			u[k] = A[k] - P[k];
			v[k] = B[k] - P[k];
			uu += u[k] * u[k]; vv += v[k] * v[k]; uv += u[k] * v[k];
		}
		c = uv / sqrt(uu * vv);
		a[i] = acos(c < -1 ? -1 : c > 1 ? 1 : c) * 180 / MS_PI;
	}
	free(t);
	free(x);
	*n = 3 * nt;
	return a;
}

/* the share of the vertices off the border and features of valence 5 to 7 */
static double ms_v567(qaws_mesh_3d const* m)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nfe = qaws_mesh_3d_get_feature_edge_count(m), i, n = 0, k = 0;
	unsigned int* fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
	unsigned char* on = (unsigned char*)calloc(nv + 1, 1);
	qaws_mesh_3d_get_feature_edges(fe, m);
	for (i = 0; i < 2 * nfe; i++)
		on[fe[i]] = 1;
	for (i = 0; i < nv; i++)
	{
		unsigned int ring[64], c = 0;
		if (on[i] || qaws_mesh_3d_is_boundary_vertex(m, i))
			continue;
		qaws_mesh_3d_get_vertex_ring(ring, &c, m, 64, i);
		n++;
		k += c >= 5 && c <= 7;
	}
	free(fe);
	free(on);
	return n ? (double)k / n : 1;
}

/* the largest distance to the reference of the faces' centres, edge midpoints and quarter points */
static double ms_sampled_deviation(qaws_mesh_3d const* m, qaws_mesh_reference const* r)
{
	unsigned int nf = qaws_mesh_3d_get_triangle_count(m), i, j, h, k, n = 0;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nf ? nf : 1));
	double* x = ms_positions(m);
	qaws_scalar *s = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 30 * (nf ? nf : 1)), *d = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 10 * (nf ? nf : 1));
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * 10 * (nf ? nf : 1));
	double dmax = 0;
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nf; i++)
	{
		for (k = 0; k < 3; k++)
			s[3 * n + k] = (qaws_scalar)((x[3 * t[3 * i] + k] + x[3 * t[3 * i + 1] + k] + x[3 * t[3 * i + 2] + k]) / 3);
		n++;
		for (j = 0; j < 3; j++)
			for (h = 1; h <= 3; h++, n++)
				for (k = 0; k < 3; k++)
					s[3 * n + k] = (qaws_scalar)((1 - 0.25 * h) * x[3 * t[3 * i + j] + k] + 0.25 * h * x[3 * t[3 * i + (j + 1) % 3] + k]);
	}
	qaws_mesh_reference_project(l, d, r, s, NULL, n);
	for (i = 0; i < n; i++)
		dmax = fmax(dmax, (double)d[i]);
	free(t); free(x); free(s); free(d); free(l);
	return dmax;
}

/* the mesh, each triangle coloured by its smallest angle, its edges drawn */
static void ms_draw_by_angle(er_canvas* c, qaws_mesh_3d const* m)
{
	static double const wire[3] = { 40, 45, 60 }, red[3] = { 214, 39, 40 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), na, i;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt ? nt : 1));
	double *x = ms_positions(m), *n = (double*)malloc(sizeof(double) * 3 * (nv ? nv : 1)), *a = ms_corner_angles(m, &na);
	qaws_scalar* ns = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv ? nv : 1));
	qaws_mesh_3d_get_triangles(t, m);
	qaws_mesh_3d_get_vertex_normals(ns, m);
	for (i = 0; i < 3 * nv; i++)
		n[i] = (double)ns[i];
	for (i = 0; i < nt; i++)
	{
		double lo = fmin(a[3 * i], fmin(a[3 * i + 1], a[3 * i + 2])), col[3];
		if (lo < 30)
			memcpy(col, red, sizeof col);
		else
			ms_colormap((60 - lo) / 30, col);
		er_triangle(c, x + 3 * t[3 * i], x + 3 * t[3 * i + 1], x + 3 * t[3 * i + 2], n + 3 * t[3 * i], n + 3 * t[3 * i + 1], n + 3 * t[3 * i + 2], col, col, col);
	}
	ms_draw(c, m, NULL, NULL, 0, wire, 0.8);
	free(t); free(x); free(n); free(a); free(ns);
}

/* both meshes' corner angles as histograms of 2-degree bins over [20, 130], the bounds as lines */
static void ms_angle_histogram(er_canvas* c, int px, int py, int pw, int ph, qaws_mesh_3d const* a, qaws_mesh_3d const* b, char const* title)
{
	static double const grey[3] = { 120, 130, 150 }, orange[3] = { 232, 120, 30 }, axis[3] = { 60, 60, 60 }, bound[3] = { 214, 39, 40 };
	double box[12] = { 18, 0, -4, 132, 0, -4, 18, 0, 64, 132, 0, 64 }, h[2][55], top = 1;
	unsigned int i, j, k, na, nb;
	double *aa = ms_corner_angles(a, &na), *ab = ms_corner_angles(b, &nb);
	char lab[16];
	memset(h, 0, sizeof h);
	for (i = 0; i < na; i++)
		if (aa[i] >= 20 && aa[i] < 130)
			h[0][(int)((aa[i] - 20) / 2)]++;
	for (i = 0; i < nb; i++)
		if (ab[i] >= 20 && ab[i] < 130)
			h[1][(int)((ab[i] - 20) / 2)]++;
	for (j = 0; j < 2; j++)
		for (i = 0; i < 55; i++)
			top = fmax(top, h[j][i]);
	er_panel(c, px, py + 70, pw, ph - 110, box, 4, 0, 0, 0.04);
	er_text(c, px + 14, py + 12, title, 3, 30, 30, 30);
	for (j = 0; j < 2; j++)
		for (i = 0; i < 55; i++)
			if (h[j][i] > 0)
			{
				double p0[3], p1[3];
				p0[0] = p1[0] = 20 + 2 * i + (j ? 1.35 : 0.65);
				p0[1] = p1[1] = 0;
				p0[2] = 0;
				p1[2] = 60 * h[j][i] / top;
				er_line(c, p0, p1, j ? orange : grey, 3.2);
			}
	{
		double p0[3] = { 20, 0, 0 }, p1[3] = { 130, 0, 0 };
		er_line(c, p0, p1, axis, 1.4);
	}
	for (k = 0; k < 2; k++)
	{
		double x0 = k ? 86 : 35, p0[3], p1[3];
		p0[0] = p1[0] = x0;
		p0[1] = p1[1] = 0;
		p0[2] = 0;
		p1[2] = 62;
		er_line(c, p0, p1, bound, 1.6);
	}
	for (k = 30; k <= 120; k += 30)
	{
		double p[3], s[3];
		p[0] = k; p[1] = 0; p[2] = 0;
		er_project(c, p, s);
		sprintf(lab, "%u", k);
		er_text(c, (int)(s[0] / c->ss) - 9, (int)(s[1] / c->ss) + 8, lab, 2, 80, 80, 80);
	}
	{
		unsigned int oa = 0, ob = 0;
		for (i = 0; i < na; i++)
			oa += aa[i] < 35 || aa[i] > 86;
		for (i = 0; i < nb; i++)
			ob += ab[i] < 35 || ab[i] > 86;
		sprintf(lab, "%u out", oa);
		er_text(c, px + 14, py + 46, "default", 2, 120, 130, 150);
		er_text(c, px + 110, py + 46, lab, 2, 120, 130, 150);
		sprintf(lab, "%u out", ob);
		er_text(c, px + 230, py + 46, "[35, 86]", 2, 232, 120, 30);
		er_text(c, px + 340, py + 46, lab, 2, 232, 120, 30);
	}
	free(aa);
	free(ab);
}

static void ms_angle_sub(char* sub, qaws_mesh_3d const* m, qaws_mesh_report const* rep)
{
	sprintf(sub, "angles %.1f..%.1f, V567 %.1f%%, off by %.1e", (double)rep->quality.min_angle, (double)rep->quality.max_angle, 100 * ms_v567(m),
		(double)rep->max_deviation);
}

static void demo_remesh_angles(void)
{
	static char const* const name[3] = { "the blob", "adaptive bumpy sphere", "the can" };
	static double const view[3][2] = { { 0.6, 0.5 }, { 0.9, 0.35 }, { 0.35, 0.5 } }, L[3] = { 0.11, 0, 0.12 };
	er_canvas* c = er_create(1950, 1720, 3, 252, 252, 250);
	qaws_mesh_3d *in[3], *def[3] = { NULL, NULL, NULL }, *bnd[3] = { NULL, NULL, NULL }, *tol = NULL;
	qaws_mesh_reference* r[3] = { NULL, NULL, NULL };
	qaws_size_field* sf = NULL;
	qaws_curvature_size_desc cd;
	qaws_mesh_report rep;
	char title[120], sub[200];
	unsigned int i;
	if (!c)
		return;
	in[0] = ms_blob(1);
	in[1] = ms_bumpy();
	in[2] = ms_can(8, 1, 1.6);
	qaws_mesh_3d_detect_features(in[2], (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	for (i = 0; i < 3; i++)
		qaws_mesh_reference_create(&r[i], in[i], NULL);
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)(12 * MS_PI / 180);
	cd.min_length = (qaws_scalar)0.02;
	cd.max_length = (qaws_scalar)0.3;
	cd.gradation = (qaws_scalar)0.35;
	if (r[1])
		qaws_size_field_create_curvature(&sf, r[1], &cd);
	for (i = 0; i < 3; i++)
	{
		if (!r[i] || (i == 1 && !sf))
			continue;
		def[i] = ms_remesh_bounded(in[i], r[i], L[i], i == 1 ? sf : NULL, 0, 0, 0, &rep);
		if (def[i])
		{
			ms_angle_sub(sub, def[i], &rep);
			sprintf(title, "%s, default", name[i]);
			ms_panel(c, 650 * (int)i, 0, 650, 640, def[i], view[i][0], view[i][1], title, sub);
			ms_draw_by_angle(c, def[i]);
			printf("  %s, default: %u triangles, %s\n", name[i], qaws_mesh_3d_get_face_count(def[i]), sub);
		}
		bnd[i] = ms_remesh_bounded(in[i], r[i], L[i], i == 1 ? sf : NULL, 35, 86, 0, &rep);
		if (bnd[i])
		{
			ms_angle_sub(sub, bnd[i], &rep);
			sprintf(title, "%s, [35, 86]", name[i]);
			ms_panel(c, 650 * (int)i, 640, 650, 640, bnd[i], view[i][0], view[i][1], title, sub);
			ms_draw_by_angle(c, bnd[i]);
			printf("  %s, angles in [35, 86]: %u triangles, %s\n", name[i], qaws_mesh_3d_get_face_count(bnd[i]), sub);
		}
		if (def[i] && bnd[i])
		{
			sprintf(title, "angles: %s", name[i]);
			ms_angle_histogram(c, 650 * (int)i, 1280, 650, 440, def[i], bnd[i], title);
		}
	}
	/* the blob within 2e-3 of its reference */
	if (r[0] && def[0])
	{
		tol = ms_remesh_bounded(in[0], r[0], L[0], NULL, 0, 0, 2e-3, &rep);
		if (tol)
			printf("  the blob, tolerance 2e-3: %u vertices (default %u), sampled distance %.3e (default %.3e), angles %.1f..%.1f\n",
				qaws_mesh_3d_get_vertex_count(tol), qaws_mesh_3d_get_vertex_count(def[0]), ms_sampled_deviation(tol, r[0]), ms_sampled_deviation(def[0], r[0]),
				(double)rep.quality.min_angle, (double)rep.quality.max_angle);
	}
	er_save(c, "showcase/mesh10_remesh_angles.png");
	printf("mesh10_remesh_angles.png\n");
	qaws_size_field_destroy(sf);
	for (i = 0; i < 3; i++)
	{
		qaws_mesh_reference_destroy(r[i]);
		qaws_mesh_3d_destroy(in[i]);
		qaws_mesh_3d_destroy(def[i]);
		qaws_mesh_3d_destroy(bnd[i]);
	}
	qaws_mesh_3d_destroy(tol);
	er_destroy(c);
}
