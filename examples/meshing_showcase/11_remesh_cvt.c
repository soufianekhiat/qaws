/*
 * Figure 11: remeshing by a centroidal Voronoi tessellation (CVT).
 *   - the blob by LOCAL and by CVT at the same target, the vertices
 *     coloured by valence (5 blue, 6 white, 7 red): CVT's lattice is
 *     regular but for a few defects;
 *   - the can by CVT with the dual of its Delaunay mesh on the reference
 *     (circumcentres projected; the solver's own restricted Voronoi cells
 *     are not exported), the caps and the side tinted: no cell crosses a rim;
 *   - the cube with creases by CVT, the crease seeds (locked) marked.
 */

#include <time.h>

/* valence colours: 4 and below dark blue, 5 blue, 6 white, 7 red, 8 and above dark red */
static void ms_valence_rgb(unsigned int n, double* rgb)
{
	static double const k[5][3] = { { 30, 40, 120 }, { 59, 110, 220 }, { 240, 240, 240 }, { 214, 39, 40 }, { 110, 0, 20 } };
	unsigned int i = n <= 4 ? 0 : n >= 8 ? 4 : n - 4;
	rgb[0] = k[i][0]; rgb[1] = k[i][1]; rgb[2] = k[i][2];
}

/* interior vertices of valence 6 (share), the edges' coefficient of variation */
static void ms_regularity(qaws_mesh_3d const* m, double* v6, double* cv)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), ne = qaws_mesh_3d_get_edge_count(m), i, in = 0, r6 = 0;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
	double* x = ms_positions(m), s = 0, s2 = 0;
	qaws_mesh_3d_get_edges(m, e);
	for (i = 0; i < ne; i++)
	{
		double const *a = x + 3 * e[2 * i], *b = x + 3 * e[2 * i + 1];
		double l = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
		s += l;
		s2 += l * l;
	}
	for (i = 0; i < nv; i++)
	{
		unsigned int ring[64], n = 0;
		if (qaws_mesh_3d_is_boundary_vertex(m, i))
			continue;
		qaws_mesh_3d_get_vertex_ring(m, i, ring, 64, &n);
		in++;
		r6 += n == 6;
	}
	*v6 = in ? 100.0 * r6 / in : 0;
	*cv = ne ? sqrt(fmax(0, s2 / ne - (s / ne) * (s / ne))) / (s / ne) : 0;
	free(e);
	free(x);
}

static qaws_mesh_3d* ms_remesh_by(qaws_mesh_3d const* m, qaws_mesh_reference const* r, qaws_remesh_method method, double L, qaws_mesh_report* rep,
	double* seconds)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	clock_t t0 = clock();
	memset(&d, 0, sizeof d);
	d.method = method;
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	if (qaws_mesh_remesh(m, &d, &out, rep) != QAWS_STATUS_OK)
		printf("  remesh failed\n");
	*seconds = (double)(clock() - t0) / CLOCKS_PER_SEC;
	return out;
}

/* a remeshed blob in a panel, its vertices coloured by valence; the numbers printed and as the subtitle */
static void ms_valence_panel(er_canvas* c, int x0, qaws_mesh_3d const* m, qaws_mesh_report const* rep, double sec, char const* title)
{
	static double const tint[3] = { 178, 192, 212 }, wire[3] = { 40, 45, 60 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i;
	double v6, cv, *p = ms_positions(m), rgb[3];
	char sub[200];
	ms_regularity(m, &v6, &cv);
	sprintf(sub, "valence 6: %.1f%%, angles %.0f..%.0f", v6, (double)rep->quality.min_angle, (double)rep->quality.max_angle);
	ms_panel(c, x0, 0, 550, 640, m, 0.6, 0.5, title, sub);
	sprintf(sub, "edge CV %.3f, %.2f s", cv, sec);
	er_text(c, x0 + 14, 68, sub, 2, 95, 95, 95);
	ms_draw(c, m, tint, NULL, 1, wire, 0.8);
	for (i = 0; i < nv; i++)
	{
		unsigned int ring[64], n = 0;
		qaws_mesh_3d_get_vertex_ring(m, i, ring, 64, &n);
		ms_valence_rgb(n, rgb);
		er_point(c, p + 3 * i, rgb, n == 6 ? 2.2 : 3.6);
	}
	printf("  %s: %u vertices, valence 6 at %.1f%%, irregular %u, angles %.1f..%.1f, edges mean %.4f CV %.3f, deviation %.2e, %.2f s\n", title, nv, v6,
		rep->quality.irregular, (double)rep->quality.min_angle, (double)rep->quality.max_angle, (double)rep->quality.mean_edge, cv, (double)rep->max_deviation,
		sec);
	free(p);
}

/* the 3D circumcentre of triangle (a, b, c) */
static void ms_circumcentre(double const* a, double const* b, double const* c, double* o)
{
	double u[3], v[3], n[3], t1[3], t2[3], uu, vv, nn;
	unsigned int k;
	for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; }
	n[0] = u[1] * v[2] - u[2] * v[1]; n[1] = u[2] * v[0] - u[0] * v[2]; n[2] = u[0] * v[1] - u[1] * v[0];
	uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
	vv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
	/* (|v|^2 n x u + |u|^2 v x n) / (2 |n|^2) */
	t1[0] = n[1] * u[2] - n[2] * u[1]; t1[1] = n[2] * u[0] - n[0] * u[2]; t1[2] = n[0] * u[1] - n[1] * u[0];
	t2[0] = v[1] * n[2] - v[2] * n[1]; t2[1] = v[2] * n[0] - v[0] * n[2]; t2[2] = v[0] * n[1] - v[1] * n[0];
	for (k = 0; k < 3; k++)
		o[k] = a[k] + (nn > 0 ? (vv * t1[k] + uu * t2[k]) / (2 * nn) : 0);
}

static void demo_remesh_cvt(void)
{
	static double const wire[3] = { 40, 45, 60 }, red[3] = { 214, 39, 40 }, dark[3] = { 25, 25, 30 }, white[3] = { 245, 245, 245 };
	static double const tints[4][3] = { { 226, 200, 160 }, { 170, 200, 170 }, { 178, 192, 212 }, { 205, 185, 220 } };
	er_canvas* c = er_create(2200, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *can = ms_can(8, 1, 1.6), *loc = NULL, *cvt = NULL, *rcan = NULL, *cube = NULL, *rcube = NULL;
	qaws_mesh_reference *rb = NULL, *rc = NULL;
	qaws_mesh_report rep;
	double sec;
	char sub[200];
	unsigned int i;
	if (!c || !blob || !can)
		return;
	qaws_mesh_reference_create(blob, NULL, &rb);
	loc = ms_remesh_by(blob, rb, QAWS_REMESH_LOCAL, 0.11, &rep, &sec);
	if (loc)
		ms_valence_panel(c, 0, loc, &rep, sec, "the blob by LOCAL");
	cvt = ms_remesh_by(blob, rb, QAWS_REMESH_CVT, 0.11, &rep, &sec);
	if (cvt)
		ms_valence_panel(c, 550, cvt, &rep, sec, "the blob by CVT");
	/* the dual of the can's Delaunay mesh: per edge, its faces' circumcentres joined (to the edge's middle on a crease), on the
	   reference */
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	qaws_mesh_reference_create(can, NULL, &rc);
	rcan = rc ? ms_remesh_by(can, rc, QAWS_REMESH_CVT, 0.12, &rep, &sec) : NULL;
	if (rcan)
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(rcan), nt = qaws_mesh_3d_get_triangle_count(rcan), ne = qaws_mesh_3d_get_edge_count(rcan), k, q;
		unsigned int *t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nt), *e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
		double *p = ms_positions(rcan), *cc = (double*)malloc(sizeof(double) * 3 * (nt + ne)), v6, cv;
		qaws_scalar *cs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nt + ne)), *cp = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nt + ne));
		qaws_mesh_location *foot = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nv), *hint = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (nt + ne)),
			*cl = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (nt + ne));
		unsigned char* feat = (unsigned char*)calloc(ne + 1, 1);
		qaws_mesh_3d_get_triangles(rcan, t);
		qaws_mesh_3d_get_edges(rcan, e);
		qaws_mesh_3d_get_footpoints(rcan, foot);
		{
			unsigned int nfe = qaws_mesh_3d_get_feature_edge_count(rcan), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
			qaws_mesh_3d_get_feature_edges(rcan, fe);
			for (q = 0; q < nfe; q++)
				for (i = 0; i < ne; i++)
					if ((e[2 * i] == fe[2 * q] && e[2 * i + 1] == fe[2 * q + 1]) || (e[2 * i] == fe[2 * q + 1] && e[2 * i + 1] == fe[2 * q]))
						feat[i] = 1;
			free(fe);
		}
		ms_regularity(rcan, &v6, &cv);
		sprintf(sub, "valence 6: %.1f%%, angles %.0f..%.0f", v6, (double)rep.quality.min_angle, (double)rep.quality.max_angle);
		ms_panel(c, 1100, 0, 550, 640, rcan, 0.35, 0.5, "the can's Delaunay dual", sub);
		sprintf(sub, "caps and side tinted, edge CV %.3f, %.2f s", cv, sec);
		er_text(c, 1100 + 14, 68, sub, 2, 95, 95, 95);
		/* the caps and the side tinted by height, flat */
		for (i = 0; i < nt; i++)
		{
			double const *a = p + 3 * t[3 * i], *b = p + 3 * t[3 * i + 1], *d = p + 3 * t[3 * i + 2];
			double z = (a[2] + b[2] + d[2]) / 3;
			double const* col = tints[z > 0.79 ? 0 : z < -0.79 ? 1 : 3];
			er_triangle(c, a, b, d, NULL, NULL, NULL, col, col, col);
		}
		/* circumcentres (faces) and edge middles, projected onto the reference */
		for (i = 0; i < nt; i++)
		{
			ms_circumcentre(p + 3 * t[3 * i], p + 3 * t[3 * i + 1], p + 3 * t[3 * i + 2], cc + 3 * i);
			hint[i] = foot[t[3 * i]];
		}
		for (i = 0; i < ne; i++)
		{
			for (k = 0; k < 3; k++)
				cc[3 * (nt + i) + k] = 0.5 * (p[3 * e[2 * i] + k] + p[3 * e[2 * i + 1] + k]);
			hint[nt + i] = foot[e[2 * i]];
		}
		for (i = 0; i < 3 * (nt + ne); i++)
			cs[i] = (qaws_scalar)cc[i];
		qaws_mesh_reference_project(rc, cs, nt + ne, hint, cl, NULL);
		qaws_mesh_reference_evaluate(rc, cl, nt + ne, cp, NULL, NULL);
		for (i = 0; i < 3 * (nt + ne); i++)
			cc[i] = (double)cp[i];
		/* each edge's dual: its two faces' centres, or a face's centre to the edge's middle on a crease */
		for (i = 0; i < ne; i++)
		{
			unsigned int f[2], nf = 0;
			for (q = 0; q < nt && nf < 2; q++)
				for (k = 0; k < 3; k++)
				{
					unsigned int u = t[3 * q + k], v = t[3 * q + (k + 1) % 3];
					if ((u == e[2 * i] && v == e[2 * i + 1]) || (v == e[2 * i] && u == e[2 * i + 1]))
						f[nf++] = q;
				}
			if (nf == 2 && !feat[i])
				er_line(c, cc + 3 * f[0], cc + 3 * f[1], dark, 1.6);
			else
				for (k = 0; k < nf; k++)
					er_line(c, cc + 3 * f[k], cc + 3 * (nt + i), dark, 1.6);
		}
		for (i = 0; i < nv; i++)
			er_point(c, p + 3 * i, red, 2.0);
		printf("  the can by CVT: %u vertices, valence 6 at %.1f%%, angles %.1f..%.1f, edge CV %.3f, %.2f s\n", nv, v6, (double)rep.quality.min_angle,
			(double)rep.quality.max_angle, cv, sec);
		free(t); free(e); free(p); free(cc); free(cs); free(cp); free(foot); free(hint); free(cl); free(feat);
	}
	/* the cube with its creases, the crease seeds marked */
	{
		static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
		static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
		static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
		cube = ms_make(x, 8, f, fs, 6);
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)(30 * MS_PI / 180), 0, 0);
		rcube = ms_remesh_by(cube, NULL, QAWS_REMESH_CVT, 0.25, &rep, &sec);
		if (rcube)
		{
			unsigned int nv = qaws_mesh_3d_get_vertex_count(rcube), nfe = qaws_mesh_3d_get_feature_edge_count(rcube);
			unsigned int* fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
			unsigned char* on = (unsigned char*)calloc(nv, 1);
			double* p = ms_positions(rcube), v6, cv;
			unsigned int seeds = 0;
			ms_regularity(rcube, &v6, &cv);
			qaws_mesh_3d_get_feature_edges(rcube, fe);
			for (i = 0; i < nfe; i++)
				on[fe[2 * i]] = on[fe[2 * i + 1]] = 1;
			for (i = 0; i < nv; i++)
				seeds += on[i];
			sprintf(sub, "%u crease seeds (locked), %u edges", seeds, nfe);
			ms_panel(c, 1650, 0, 550, 640, rcube, 0.62, 0.55, "the cube with creases", sub);
			sprintf(sub, "valence 6 %.1f%%, angles %.0f..%.0f, CV %.3f", v6, (double)rep.quality.min_angle, (double)rep.quality.max_angle, cv);
			er_text(c, 1650 + 14, 68, sub, 2, 95, 95, 95);
			ms_draw(c, rcube, tints[2], NULL, 0, wire, 1.0);
			for (i = 0; i < nfe; i++)
				er_line(c, p + 3 * fe[2 * i], p + 3 * fe[2 * i + 1], red, 2.4);
			for (i = 0; i < nv; i++)
				if (on[i])
					er_point(c, p + 3 * i, qaws_mesh_3d_is_feature_corner(rcube, i) ? dark : white, qaws_mesh_3d_is_feature_corner(rcube, i) ? 6.0 : 4.0);
			printf("  the cube by CVT: %u vertices, %u crease edges, angles %.1f..%.1f, edge CV %.3f, %.2f s\n", nv, nfe, (double)rep.quality.min_angle,
				(double)rep.quality.max_angle, cv, sec);
			free(fe); free(on); free(p);
		}
	}
	er_save(c, "showcase/mesh11_remesh_cvt.png");
	printf("mesh11_remesh_cvt.png\n");
	qaws_mesh_reference_destroy(rb);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(loc);
	qaws_mesh_3d_destroy(cvt);
	qaws_mesh_3d_destroy(rcan);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(rcube);
	er_destroy(c);
}
