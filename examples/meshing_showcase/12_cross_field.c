/*
 * Figure 12: smoothest cross fields (Knoppel et al. 2013) on a sparse LDL^T.
 *   - the blob and the can (dense remeshes, the can's rims held): Instant
 *     Meshes' LOCAL field above, the globally SMOOTHEST below; singular faces
 *     dotted (red +1, blue -1);
 *   - the can all quads by INSTANT, which now keeps the global field, from
 *     above and below: the irregular vertices dotted (red 3, blue 5 and
 *     more) and counted per cap;
 *   - a torus: the smoothest field (pairs of singularities on a fat torus)
 *     and the one drawn to curvature (none, along the parallels).
 */

#include "mesh/qaws_cross_field.h"
#include <time.h>

/* crosses at the vertices (half-size h), the singular faces dotted; returns the counts of +1 and -1 */
static void ms_cross_field(er_canvas* c, qaws_mesh_3d const* m, qaws_scalar const* dir, int const* idx, double h, unsigned int* plus,
	unsigned int* minus)
{
	static double const red[3] = { 214, 39, 40 }, blue[3] = { 31, 90, 200 }, dark[3] = { 60, 60, 70 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_face_count(m), i, k;
	double* p = ms_positions(m);
	qaws_scalar* nrm = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	qaws_mesh_3d_get_vertex_normals(m, nrm);
	qaws_mesh_3d_get_triangles(m, fv);
	for (i = 0; i < nv; i++)
	{
		double const* x = p + 3 * i;
		double d0[3], d1[3], a[3], b[3];
		d0[0] = dir[3 * i]; d0[1] = dir[3 * i + 1]; d0[2] = dir[3 * i + 2];
		d1[0] = nrm[3 * i + 1] * d0[2] - nrm[3 * i + 2] * d0[1];
		d1[1] = nrm[3 * i + 2] * d0[0] - nrm[3 * i] * d0[2];
		d1[2] = nrm[3 * i] * d0[1] - nrm[3 * i + 1] * d0[0];
		for (k = 0; k < 3; k++) { a[k] = x[k] - h * d0[k]; b[k] = x[k] + h * d0[k]; }
		er_line(c, a, b, dark, 1.2);
		for (k = 0; k < 3; k++) { a[k] = x[k] - h * d1[k]; b[k] = x[k] + h * d1[k]; }
		er_line(c, a, b, dark, 1.2);
	}
	*plus = *minus = 0;
	for (i = 0; i < nf; i++)
		if (idx[i])
		{
			double cen[3];
			for (k = 0; k < 3; k++)
				cen[k] = (p[3 * fv[3 * i] + k] + p[3 * fv[3 * i + 1] + k] + p[3 * fv[3 * i + 2] + k]) / 3;
			er_point(c, cen, idx[i] > 0 ? red : blue, 7);
			*plus += idx[i] > 0;
			*minus += idx[i] < 0;
		}
	free(p); free(nrm); free(fv);
}

/* one field panel: the mesh shaded, its field (LOCAL or SMOOTHEST, maybe drawn to curvature) as crosses */
static void ms_field_panel(er_canvas* c, int x, int y, qaws_mesh_3d const* m, qaws_direction_field_method method, int aligned, double h,
	double yaw, double pitch, char const* title, double const* tint)
{
	qaws_direction_field_desc d;
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nf = qaws_mesh_3d_get_face_count(m), plus = 0, minus = 0;
	qaws_scalar* dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	int* idx = (int*)malloc(sizeof(int) * nf);
	char sub[200];
	clock_t t0;
	double ms;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.method = method;
	d.curvature_aligned = aligned;
	t0 = clock();
	st = qaws_mesh_3d_direction_field(m, &d, dir, idx);
	ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
	ms_panel(c, x, y, 550, 600, m, yaw, pitch, title, NULL);
	ms_draw(c, m, tint, NULL, 1, NULL, 0);
	if (st == QAWS_STATUS_OK)
		ms_cross_field(c, m, dir, idx, h, &plus, &minus);
	sprintf(sub, "%u vertices: %u of +1, %u of -1; %.0f ms", nv, plus, minus, ms);
	er_text(c, x + 14, y + 46, sub, 2, 95, 95, 95);
	printf("  %s: %s\n", title, sub);
	free(dir);
	free(idx);
}

/* the quads' irregular vertices dotted (red: 3 faces, blue: 5 or more; c NULL: only counted), those on the can's top and
   bottom caps (rims included) counted */
static unsigned int ms_irregular(er_canvas* c, qaws_mesh_3d const* q, unsigned int* top, unsigned int* bottom)
{
	static double const red[3] = { 214, 39, 40 }, blue[3] = { 31, 90, 200 };
	unsigned int nv = qaws_mesh_3d_get_vertex_count(q), nf = qaws_mesh_3d_get_face_count(q), nc = qaws_mesh_3d_get_corner_count(q), i, n = 0;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	unsigned int* val = (unsigned int*)calloc(nv, sizeof(unsigned int));
	double* p = ms_positions(q);
	qaws_mesh_3d_get_faces(q, fv, fs);
	for (i = 0; i < nc; i++)
		val[fv[i]]++;
	*top = *bottom = 0;
	for (i = 0; i < nv; i++)
		if (val[i] != 4)
		{
			if (c)
				er_point(c, p + 3 * i, val[i] < 4 ? red : blue, 7);
			*top += p[3 * i + 2] > 0.8 - 1e-6;
			*bottom += p[3 * i + 2] < -0.8 + 1e-6;
			n++;
		}
	free(fv); free(fs); free(val); free(p);
	return n;
}

static void demo_cross_field(void)
{
	static double const tint[3] = { 178, 192, 212 }, warm[3] = { 226, 200, 160 }, green[3] = { 170, 205, 160 }, wire[3] = { 40, 45, 60 };
	er_canvas* c = er_create(2200, 1200, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *can = ms_can(8, 1, 1.6), *torus = ms_torus(48, 18, 2, 0.75), *db = NULL, *dc = NULL, *q = NULL;
	qaws_remesh_desc rd;
	qaws_quadrangulate_desc qd;
	qaws_mesh_report rep;
	char sub[200];
	if (!c || !blob || !can || !torus)
		return;
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	memset(&rd, 0, sizeof rd);
	rd.target_length = (qaws_scalar)0.09;
	qaws_mesh_remesh(blob, &rd, &db, NULL);
	rd.target_length = (qaws_scalar)0.07;
	qaws_mesh_remesh(can, &rd, &dc, NULL);
	/* columns 1 and 2: LOCAL above, SMOOTHEST below */
	if (db)
	{
		ms_field_panel(c, 0, 0, db, QAWS_DIRECTION_FIELD_LOCAL, 0, 0.032, 0.6, 0.5, "the blob, LOCAL", tint);
		ms_field_panel(c, 0, 600, db, QAWS_DIRECTION_FIELD_SMOOTHEST, 0, 0.032, 0.6, 0.5, "the blob, SMOOTHEST", tint);
	}
	if (dc)
	{
		ms_field_panel(c, 550, 0, dc, QAWS_DIRECTION_FIELD_LOCAL, 0, 0.026, 0.35, 0.75, "the can, LOCAL", warm);
		ms_field_panel(c, 550, 600, dc, QAWS_DIRECTION_FIELD_SMOOTHEST, 0, 0.026, 0.35, 0.75, "the can, SMOOTHEST", warm);
		/* the factor's fill */
		{
			mf_hier h;
			cf_system* sys = NULL;
			double fl = 0, fa = 0;
			clock_t t0 = clock();
			if (mf_build(&h, &dc->hm, NULL) == QAWS_STATUS_OK)
			{
				if (cf_create(&dc->hm, h.lv[0].nrm, &h.lv[0], 4, &sys) == QAWS_STATUS_OK)
				{
					cf_fill(sys, &fl, &fa);
					cf_free(sys);
				}
				mf_free(&h);
			}
			printf("  can's factor: nnz(L) %.0f, nnz(A) %.0f (lower), fill %.1f; built and factored in %.0f ms\n", fl, fa, fa > 0 ? fl / fa : 0,
				1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC);
		}
	}
	/* column 3: the can all quads by INSTANT, its caps from above and below */
	memset(&qd, 0, sizeof qd);
	qd.target_length = (qaws_scalar)0.14;
	qd.pure = 1;
	{
		clock_t t0 = clock();
		qaws_mesh_quadrangulate(can, &qd, &q, &rep);
		if (q)
		{
			unsigned int nfe = qaws_mesh_3d_get_feature_edge_count(q), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1)), i, n, top, bottom;
			double* p = ms_positions(q);
			static double const red[3] = { 214, 39, 40 };
			char sub2[200];
			double sec = (double)(clock() - t0) / CLOCKS_PER_SEC;
			qaws_mesh_3d_get_feature_edges(q, fe);
			n = ms_irregular(NULL, q, &top, &bottom);
			/* the count before this item was measured on the LOCAL field alone: no setting brings that field back */
			sprintf(sub, "%u irregular (22 measured before); top cap %u", rep.quality.irregular, top);
			sprintf(sub2, "%u quads; bottom cap %u irregular; %.2f s", rep.quality.quads, bottom, sec);
			ms_panel(c, 1100, 0, 550, 600, q, 0.35, 1.0, "INSTANT quads, the top", sub);
			ms_draw(c, q, warm, NULL, 1, wire, 1.0);
			for (i = 0; i < nfe; i++)
				er_line(c, p + 3 * fe[2 * i], p + 3 * fe[2 * i + 1], red, 2.2);
			ms_irregular(c, q, &top, &bottom);
			ms_panel(c, 1100, 600, 550, 600, q, 0.35, -1.0, "INSTANT quads, the bottom", sub2);
			ms_draw(c, q, warm, NULL, 1, wire, 1.0);
			for (i = 0; i < nfe; i++)
				er_line(c, p + 3 * fe[2 * i], p + 3 * fe[2 * i + 1], red, 2.2);
			ms_irregular(c, q, &top, &bottom);
			printf("  can all quads: %s; %s; %u vertices not of 4 faces (side %u), %u crease edges\n", sub, sub2, n, n - top - bottom, nfe);
			free(fe);
			free(p);
		}
	}
	/* column 4: the torus, smoothest and drawn to curvature */
	ms_field_panel(c, 1650, 0, torus, QAWS_DIRECTION_FIELD_SMOOTHEST, 0, 0.09, 0.5, 0.75, "a torus, SMOOTHEST", green);
	ms_field_panel(c, 1650, 600, torus, QAWS_DIRECTION_FIELD_SMOOTHEST, 1, 0.09, 0.5, 0.75, "drawn to curvature", green);
	er_save(c, "showcase/mesh12_cross_field.png");
	printf("mesh12_cross_field.png\n");
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(torus);
	qaws_mesh_3d_destroy(db);
	qaws_mesh_3d_destroy(dc);
	qaws_mesh_3d_destroy(q);
	er_destroy(c);
}
