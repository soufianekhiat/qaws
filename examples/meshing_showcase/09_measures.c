/*
 * Figure 9: measures in the forms remeshing papers report, the deviation
 * bounded from both sides.
 *   - figure 3's uniform blob (LOCAL at 0.11): faces by Q, the vertices off
 *     their ideal valence dotted (blue below it, red above);
 *   - figure 4's can in quads (INSTANT, pure, at 0.14), turned toward its
 *     worst quad (ringed): quads by their smallest scaled Jacobian (red
 *     under 0.5), irregular vertices dotted (the library's rule), the
 *     creases' corners and edges regular;
 *   - the blob's faces by their centres' distance to the reference, the
 *     worst point of each direction marked, both certified intervals;
 *   - the baseline of figures 3 and 4's meshes (the full table on stdout).
 */

#include <time.h>

/* per face: a triangle's Q, a quad's smallest scaled Jacobian, -1 otherwise */
static double* ms_face_quality(qaws_mesh_3d const* m)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), f, q, i, k;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	double *x = ms_positions(m), *out = (double*)malloc(sizeof(double) * (nf + 1));
	qaws_mesh_3d_get_faces(m, fv, fs);
	for (f = 0, q = 0; f < nf; q += fs[f++])
	{
		double const* p[4];
		out[f] = -1;
		for (i = 0; i < fs[f] && i < 4; i++)
			p[i] = x + 3 * fv[q + i];
		if (fs[f] == 3)
		{
			double l[3], s = 0, h = 0, u[3], v[3], n[3];
			for (i = 0; i < 3; i++)
			{
				double const *a = p[i], *b = p[(i + 1) % 3];
				l[i] = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
				s += 0.5 * l[i];
				h = fmax(h, l[i]);
			}
			for (k = 0; k < 3; k++)
			{
				u[k] = p[1][k] - p[0][k];
				v[k] = p[2][k] - p[0][k];
			}
			n[0] = u[1] * v[2] - u[2] * v[1];
			n[1] = u[2] * v[0] - u[0] * v[2];
			n[2] = u[0] * v[1] - u[1] * v[0];
			out[f] = s * h > 0 ? sqrt(3.0) * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / (s * h) : 0;
		}
		else if (fs[f] == 4)
		{
			double n[3] = { 0, 0, 0 }, ln, j = 1e300;
			for (i = 0; i < 4; i++)
			{
				double const *a = p[i], *b = p[(i + 1) % 4];
				n[0] += a[1] * b[2] - a[2] * b[1];
				n[1] += a[2] * b[0] - a[0] * b[2];
				n[2] += a[0] * b[1] - a[1] * b[0];
			}
			ln = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
			for (i = 0; i < 4; i++)
			{
				double en[3], ep[3], c[3], d;
				for (k = 0; k < 3; k++)
				{
					en[k] = p[(i + 1) % 4][k] - p[i][k];
					ep[k] = p[(i + 3) % 4][k] - p[i][k];
				}
				c[0] = en[1] * ep[2] - en[2] * ep[1];
				c[1] = en[2] * ep[0] - en[0] * ep[2];
				c[2] = en[0] * ep[1] - en[1] * ep[0];
				d = sqrt((en[0] * en[0] + en[1] * en[1] + en[2] * en[2]) * (ep[0] * ep[0] + ep[1] * ep[1] + ep[2] * ep[2])) * ln;
				j = fmin(j, d > 0 ? (c[0] * n[0] + c[1] * n[1] + c[2] * n[2]) / d : 0);
			}
			out[f] = j;
		}
	}
	free(fv); free(fs); free(x);
	return out;
}

/* a quality in [0, 1]: light when good, through amber to orange; below 0 (failed) red */
static void ms_quality_colour(double t, double* rgb)
{
	static double const good[3] = { 236, 238, 232 }, mid[3] = { 247, 200, 110 }, poor[3] = { 232, 120, 40 }, bad[3] = { 200, 20, 30 };
	unsigned int k;
	t = t > 1 ? 1 : t;
	for (k = 0; k < 3; k++)
		rgb[k] = t < 0 ? bad[k] : t > 0.5 ? mid[k] + (t - 0.5) / 0.5 * (good[k] - mid[k]) : poor[k] + t / 0.5 * (mid[k] - poor[k]);
}

/* the faces drawn flat, a colour each (3 per face), fanned from their first vertex */
static void ms_draw_faces(er_canvas* c, qaws_mesh_3d const* m, double const* rgb)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), f, q, i;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	double* x = ms_positions(m);
	qaws_mesh_3d_get_faces(m, fv, fs);
	for (f = 0, q = 0; f < nf; q += fs[f++])
		for (i = 1; i + 1 < fs[f]; i++)
			er_triangle(c, x + 3 * fv[q], x + 3 * fv[q + i], x + 3 * fv[q + i + 1], NULL, NULL, NULL, rgb + 3 * f, rgb + 3 * f, rgb + 3 * f);
	free(fv); free(fs); free(x);
}

/* the vertices off their ideal valence: blue below it, red above */
static void ms_dot_irregular(er_canvas* c, qaws_mesh_3d const* m, unsigned int const* ideal)
{
	static double const red[3] = { 214, 39, 40 }, blue[3] = { 31, 90, 200 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), v, ring[64], n;
	double* x = ms_positions(m);
	for (v = 0; v < nv; v++)
	{
		n = 0;
		qaws_mesh_3d_get_vertex_ring(m, v, ring, 64, &n);
		if (ideal[v] && n != ideal[v])
			er_point(c, x + 3 * v, n < ideal[v] ? blue : red, 4.5);
	}
	free(x);
}

typedef struct ms_measured
{
	char const* name;
	qaws_mesh_3d* mesh;
	qaws_mesh_reference const* ref;     /* the surface it was made on */
	qaws_shape_quality s;
	qaws_mesh_deviation d;
	double seconds;
} ms_measured;

static void ms_measure(ms_measured* r)
{
	clock_t c0 = clock();
	memset(&r->s, 0, sizeof r->s);
	memset(&r->d, 0, sizeof r->d);
	if (!r->mesh || !r->ref)
		return;
	qaws_mesh_3d_measure_shape(r->mesh, &r->s, NULL);
	qaws_mesh_3d_measure_deviation(r->mesh, r->ref, NULL, &r->d);
	r->seconds = (double)(clock() - c0) / CLOCKS_PER_SEC;
}

static void demo_measures(void)
{
	static double const wire[3] = { 40, 45, 60 }, red[3] = { 214, 39, 40 }, blue[3] = { 31, 90, 200 }, dark[3] = { 20, 20, 20 };
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *bumpy = ms_bumpy(), *can = ms_can(8, 1, 1.6), *torus = ms_torus(24, 10, 2, 0.75), *cube = NULL;
	qaws_mesh_reference *rb = NULL, *ry = NULL, *rc = NULL, *rq = NULL, *rt = NULL;
	qaws_size_field* sf = NULL;
	qaws_curvature_size_desc cd;
	qaws_mesh_report rep;
	ms_measured row[8];
	char sub[200], line[200];
	unsigned int i, k;
	if (!c || !blob || !bumpy || !can || !torus)
		return;
	{
		static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
		static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
		static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
		cube = ms_make(x, 8, f, fs, 6);
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)(30 * MS_PI / 180), 0, 0);
	}
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	qaws_mesh_reference_create(blob, NULL, &rb);
	qaws_mesh_reference_create(bumpy, NULL, &ry);
	qaws_mesh_reference_create(can, NULL, &rc);
	qaws_mesh_reference_create(cube, NULL, &rq);
	qaws_mesh_reference_create(torus, NULL, &rt);
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)(12 * MS_PI / 180);
	cd.min_length = (qaws_scalar)0.02;
	cd.max_length = (qaws_scalar)0.3;
	cd.gradation = (qaws_scalar)0.35;
	if (ry)
		qaws_size_field_create_curvature(ry, &cd, &sf);
	/* figures 3 and 4's meshes, each with the surface it was made on */
	memset(row, 0, sizeof row);
	row[0].name = "uniform blob";
	row[0].mesh = rb ? ms_remesh(blob, rb, 0.11, NULL, &rep) : NULL;
	row[0].ref = rb;
	row[1].name = "adaptive bumpy";
	row[1].mesh = sf ? ms_remesh(bumpy, ry, 0, sf, &rep) : NULL;
	row[1].ref = ry;
	row[2].name = "can";
	row[2].mesh = rc ? ms_remesh(can, rc, 0.12, NULL, &rep) : NULL;
	row[2].ref = rc;
	row[3].name = "cube";
	row[3].mesh = ms_remesh(cube, NULL, 0.22, NULL, &rep);
	row[3].ref = rq;
	row[4].name = "blob quad-dom";
	row[4].mesh = ms_quads(blob, 0.16, 0, &rep);
	row[4].ref = rb;
	row[5].name = "blob pure";
	row[5].mesh = ms_quads(blob, 0.11, 1, &rep);
	row[5].ref = rb;
	row[6].name = "can pure";
	row[6].mesh = ms_quads(can, 0.14, 1, &rep);
	row[6].ref = rc;
	row[7].name = "torus pure";
	row[7].mesh = ms_quads(torus, 0.3, 1, &rep);
	row[7].ref = rt;
	for (i = 0; i < 8; i++)
		ms_measure(&row[i]);
	printf("  %-15s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %4s  %-17s %-17s %8s %6s %6s %6s %6s %7s\n", "mesh", "#V", "#F", "Qmin", "Qavg", "amin",
		"amean", "amax", "%<30", "%>90", "V567", "irr", "H to %diag", "H from %diag", "RMS", "AD", "ED", "SJmin", "SJavg", "time s");
	for (i = 0; i < 8; i++)
	{
		ms_measured const* r = &row[i];
		double dg = r->d.diagonal > 0 ? 100 / r->d.diagonal : 0;
		if (!r->mesh)
			continue;
		printf("  %-15s %6u %6u %6.3f %6.3f %6.1f %6.1f %6.1f %6.2f %6.2f %6.3f %4u  %.4f-%.4f%s %.4f-%.4f%s %8.2e %6.2f %6.3f %6.3f %6.3f %7.3f\n", r->name,
			qaws_mesh_3d_get_vertex_count(r->mesh), qaws_mesh_3d_get_face_count(r->mesh), (double)r->s.min_quality, (double)r->s.mean_quality,
			(double)r->s.min_angle, (double)r->s.mean_min_angle, (double)r->s.max_angle, 100 * (double)r->s.below_30, 100 * (double)r->s.above_90,
			(double)r->s.regular_share, r->s.irregular, dg * r->d.to_reference[0], dg * r->d.to_reference[1], r->d.converged ? " " : "*",
			dg * r->d.from_reference[0], dg * r->d.from_reference[1], r->d.converged ? " " : "*", (double)r->d.rms, (double)r->s.angle_deviation,
			(double)r->s.edge_deviation, (double)r->s.min_jacobian, (double)r->s.mean_jacobian, r->seconds);
	}
	/* (a) the blob by Q, its irregular vertices */
	if (row[0].mesh)
	{
		qaws_mesh_3d* m = row[0].mesh;
		unsigned int nf = qaws_mesh_3d_get_face_count(m), *ideal = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_vertex_count(m));
		double *q = ms_face_quality(m), *col = (double*)malloc(sizeof(double) * 3 * nf);
		qaws_shape_quality s;
		for (i = 0; i < nf; i++)
			ms_quality_colour((q[i] - 0.6) / 0.4, col + 3 * i);
		qaws_mesh_3d_measure_shape(m, &s, ideal);
		sprintf(sub, "Q min %.2f mean %.3f, %u irregular", (double)s.min_quality, (double)s.mean_quality, s.irregular);
		ms_panel(c, 0, 0, 550, 640, m, 0.6, 0.5, "triangles by Q", sub);
		ms_draw_faces(c, m, col);
		ms_draw(c, m, NULL, NULL, 0, wire, 0.7);
		ms_dot_irregular(c, m, ideal);
		free(q); free(col); free(ideal);
	}
	/* (b) the can's quads by their smallest scaled Jacobian, turned toward the worst one (ringed) */
	if (row[6].mesh)
	{
		qaws_mesh_3d* m = row[6].mesh;
		unsigned int nf = qaws_mesh_3d_get_face_count(m), nv = qaws_mesh_3d_get_vertex_count(m), nfe = qaws_mesh_3d_get_feature_edge_count(m);
		unsigned int *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1)), *ideal = (unsigned int*)malloc(sizeof(unsigned int) * nv);
		unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
		unsigned int w = 0, wo = 0, o = 0;
		double *q = ms_face_quality(m), *col = (double*)malloc(sizeof(double) * 3 * nf), *p = ms_positions(m), cen[3] = { 0, 0, 0 }, dir[3], l;
		qaws_shape_quality s;
		for (i = 0; i < nf; i++)
			ms_quality_colour(q[i] < 0.5 ? -1 : (q[i] - 0.5) / 0.5, col + 3 * i);
		qaws_mesh_3d_get_faces(m, fv, fs);
		for (i = 0; i < nf; o += fs[i++])
			if (q[i] < q[w])
			{
				w = i;
				wo = o;
			}
		/* from the can's centre toward the worst quad's: the view's direction, tilted 40 degrees at most */
		for (i = 0; i < fs[w]; i++)
			for (k = 0; k < 3; k++)
				cen[k] += p[3 * fv[wo + i] + k] / fs[w];
		memcpy(dir, cen, sizeof dir);
		for (i = 0; i < nv; i++)
			for (k = 0; k < 3; k++)
				dir[k] -= p[3 * i + k] / nv;
		l = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
		qaws_mesh_3d_measure_shape(m, &s, ideal);
		sprintf(sub, "SJ min %.2f mean %.3f, %u irregular", (double)s.min_jacobian, (double)s.mean_jacobian, s.irregular);
		ms_panel(c, 550, 0, 550, 640, m, atan2(-dir[0], -dir[1]), asin(fmax(-0.7, fmin(0.7, l > 0 ? dir[2] / l : 0.5))), "quads by SJ, the worst ringed", sub);
		ms_draw_faces(c, m, col);
		ms_draw(c, m, NULL, NULL, 0, wire, 0.8);
		qaws_mesh_3d_get_feature_edges(m, fe);
		for (i = 0; i < nfe; i++)
			er_line(c, p + 3 * fe[2 * i], p + 3 * fe[2 * i + 1], dark, 2.0);
		ms_dot_irregular(c, m, ideal);
		{
			/* the ring: a circle in the view's plane, moved toward the eye (the view is orthographic) so no face hides it */
			double yaw = atan2(-dir[0], -dir[1]), pitch = asin(fmax(-0.7, fmin(0.7, l > 0 ? dir[2] / l : 0.5))), a[3], b[3];
			double right[3], up[3], eye[3];
			right[0] = cos(yaw); right[1] = -sin(yaw); right[2] = 0;
			up[0] = sin(pitch) * sin(yaw); up[1] = sin(pitch) * cos(yaw); up[2] = cos(pitch);
			eye[0] = -cos(pitch) * sin(yaw); eye[1] = -cos(pitch) * cos(yaw); eye[2] = sin(pitch);
			for (i = 0; i < 32; i++)
			{
				double t0 = 2 * MS_PI * i / 32, t1 = 2 * MS_PI * (i + 1) / 32;
				for (k = 0; k < 3; k++)
				{
					a[k] = cen[k] + 0.5 * eye[k] + 0.14 * (cos(t0) * right[k] + sin(t0) * up[k]);
					b[k] = cen[k] + 0.5 * eye[k] + 0.14 * (cos(t1) * right[k] + sin(t1) * up[k]);
				}
				er_line(c, a, b, dark, 3.0);
			}
		}
		printf("  can pure: worst quad %u, SJ %.3f, its centre %.3f %.3f %.3f off the can's\n", w, q[w], dir[0], dir[1], dir[2]);
		free(q); free(col); free(ideal); free(fe); free(p); free(fv); free(fs);
	}
	/* (c) the blob's deviation: faces by their centres' distance, the worst point each way */
	if (row[0].mesh && rb)
	{
		qaws_mesh_3d* m = row[0].mesh;
		qaws_mesh_deviation const* d = &row[0].d;
		unsigned int nf = qaws_mesh_3d_get_face_count(m), *t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
		double *x = ms_positions(m), *col = (double*)malloc(sizeof(double) * 3 * nf), wt[3], wf[3];
		qaws_scalar *cen = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nf), *dist = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
		qaws_mesh_location* loc = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * nf);
		qaws_mesh_3d_get_triangles(m, t);
		for (i = 0; i < nf; i++)
			for (k = 0; k < 3; k++)
				cen[3 * i + k] = (qaws_scalar)((x[3 * t[3 * i] + k] + x[3 * t[3 * i + 1] + k] + x[3 * t[3 * i + 2] + k]) / 3);
		qaws_mesh_reference_project(rb, cen, nf, NULL, loc, dist);
		for (i = 0; i < nf; i++)
			ms_colormap(0.5 + 0.5 * dist[i] / d->to_reference[1], col + 3 * i);
		for (k = 0; k < 3; k++)
		{
			wt[k] = (double)d->to_reference_at[k];
			wf[k] = (double)d->from_reference_at[k];
		}
		sprintf(sub, "to %.5f-%.5f, from %.5f-%.5f", (double)d->to_reference[0], (double)d->to_reference[1], (double)d->from_reference[0],
			(double)d->from_reference[1]);
		ms_panel(c, 1100, 0, 550, 640, m, 0.6, 0.5, "deviation, both ways", sub);
		ms_draw_faces(c, m, col);
		ms_draw(c, m, NULL, NULL, 0, wire, 0.5);
		er_point(c, wt, red, 9);
		er_point(c, wf, blue, 9);
		free(t); free(x); free(col); free(cen); free(dist); free(loc);
	}
	/* the baseline table */
	er_text(c, 1664, 12, "the baseline (figures 3 and 4)", 3, 30, 30, 30);
	er_text(c, 1664, 46, "Hausdorff both ways in % of the diagonal, certified", 2, 95, 95, 95);
	er_text(c, 1664, 92, "triangles      faces Qmin Qavg  amin amean amax  >90%  V567 irr  H %diag", 2, 30, 30, 30);
	for (i = 0; i < 4; i++)
	{
		ms_measured const* r = &row[i];
		double dg = r->d.diagonal > 0 ? 100 / r->d.diagonal : 0;
		if (!r->mesh)
			continue;
		sprintf(line, "%-14s %5u %4.2f %4.2f %5.1f %5.1f %5.1f %4.1f %5.3f %3u  %.3f-%.3f", r->name, qaws_mesh_3d_get_face_count(r->mesh),
			(double)r->s.min_quality, (double)r->s.mean_quality, (double)r->s.min_angle, (double)r->s.mean_min_angle, (double)r->s.max_angle,
			100 * (double)r->s.above_90, (double)r->s.regular_share, r->s.irregular, dg * r->d.hausdorff[0], dg * r->d.hausdorff[1]);
		er_text(c, 1664, 124 + 30 * i, line, 2, 60, 60, 70);
	}
	er_text(c, 1664, 266, "quads          faces SJmin SJavg  AD    ED   irr  H %diag", 2, 30, 30, 30);
	for (i = 4; i < 8; i++)
	{
		ms_measured const* r = &row[i];
		double dg = r->d.diagonal > 0 ? 100 / r->d.diagonal : 0;
		if (!r->mesh)
			continue;
		sprintf(line, "%-14s %5u %5.2f %5.3f %5.1f %5.3f %4u  %.3f-%.3f", r->name, qaws_mesh_3d_get_face_count(r->mesh), (double)r->s.min_jacobian,
			(double)r->s.mean_jacobian, (double)r->s.angle_deviation, (double)r->s.edge_deviation, r->s.irregular, dg * r->d.hausdorff[0],
			dg * r->d.hausdorff[1]);
		er_text(c, 1664, 298 + 30 * (i - 4), line, 2, 60, 60, 70);
	}
	er_text(c, 1664, 440, "irr: vertices off their ideal valence, from the", 2, 95, 95, 95);
	er_text(c, 1664, 466, "corners' angles: box corners 3, crease edges 4,", 2, 95, 95, 95);
	er_text(c, 1664, 492, "a straight border 3; smooth vertices 6 or 4", 2, 95, 95, 95);
	er_text(c, 1664, 530, "dots: blue below the ideal, red above", 2, 95, 95, 95);
	er_text(c, 1664, 556, "worst points: red on the mesh, blue on the reference", 2, 95, 95, 95);
	er_text(c, 1664, 582, "ring: the can's worst quad (red: SJ under 0.5)", 2, 95, 95, 95);
	er_save(c, "showcase/mesh9_measures.png");
	printf("mesh9_measures.png\n");
	for (i = 0; i < 8; i++)
		qaws_mesh_3d_destroy(row[i].mesh);
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(rb);
	qaws_mesh_reference_destroy(ry);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_reference_destroy(rq);
	qaws_mesh_reference_destroy(rt);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(bumpy);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(torus);
	qaws_mesh_3d_destroy(cube);
	er_destroy(c);
}
