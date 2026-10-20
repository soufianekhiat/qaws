/*
 * Figure 4: quads (Instant Meshes on the reference).
 *   - the blob's cross field: a cross per vertex of a dense remesh, the
 *     singularities marked (red: index +1, blue: -1);
 *   - the blob quadrangulated: quad-dominant, then all quads;
 *   - the can: quads along its creases;
 *   - a torus: all quads, no singularity needed.
 */

static qaws_mesh_3d* ms_quads(qaws_mesh_3d const* m, double L, int pure, qaws_mesh_report* rep)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* out = NULL;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.method = QAWS_QUADS_INSTANT;
	d.target_length = (qaws_scalar)L;
	d.pure = pure;
	st = qaws_mesh_quadrangulate(&out, rep, m, &d);
	if (st != QAWS_STATUS_OK)
		printf("  quadrangulate failed: %s\n", qaws_status_to_string(st));
	return out;
}

static void ms_quad_sub(char* sub, qaws_mesh_3d const* q, qaws_mesh_report const* rep)
{
	sprintf(sub, "%u quads, %u others, %u singular", rep->quality.quads, rep->quality.triangles + rep->quality.polygons, rep->quality.irregular);
	(void)q;
}

static void demo_quads(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, warm[3] = { 226, 200, 160 }, green[3] = { 170, 205, 160 };
	static double const red[3] = { 214, 39, 40 }, blue[3] = { 31, 90, 200 }, dark[3] = { 60, 60, 70 };
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *can = ms_can(8, 1, 1.6), *torus = ms_torus(24, 10, 2, 0.75), *dense = NULL, *q1 = NULL, *q2 = NULL, *q3 = NULL, *q4 = NULL;
	qaws_mesh_report rep;
	char sub[200];
	unsigned int i;
	if (!c || !blob || !can || !torus)
		return;
	/* the cross field on a dense remesh of the blob */
	{
		qaws_remesh_desc rd;
		qaws_scalar* dir;
		int* idx;
		double* p;
		unsigned int nv, nf, plus = 0, minus = 0, *fv;
		memset(&rd, 0, sizeof rd);
		rd.target_length = (qaws_scalar)0.09;
		qaws_mesh_remesh(&dense, NULL, blob, &rd);
		nv = qaws_mesh_3d_get_vertex_count(dense);
		nf = qaws_mesh_3d_get_face_count(dense);
		dir = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		idx = (int*)malloc(sizeof(int) * nf);
		fv = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
		qaws_mesh_3d_direction_field(dir, idx, dense, NULL);
		for (i = 0; i < nf; i++)
		{
			plus += idx[i] == 1;
			minus += idx[i] == -1;
		}
		sprintf(sub, "on %u vertices: %u singularities +1, %u -1", nv, plus, minus);
		ms_panel(c, 0, 0, 550, 640, dense, 0.6, 0.5, "a cross field", sub);
		ms_draw(c, dense, tint, NULL, 1, NULL, 0);
		p = ms_positions(dense);
		{
			qaws_scalar* nrm = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			qaws_mesh_3d_get_vertex_normals(nrm, dense);
			for (i = 0; i < nv; i++)
			{
				double const* x = p + 3 * i;
				double d0[3], d1[3], a[3], b[3], h = 0.032;
				unsigned int k;
				d0[0] = dir[3 * i]; d0[1] = dir[3 * i + 1]; d0[2] = dir[3 * i + 2];
				d1[0] = nrm[3 * i + 1] * d0[2] - nrm[3 * i + 2] * d0[1];
				d1[1] = nrm[3 * i + 2] * d0[0] - nrm[3 * i] * d0[2];
				d1[2] = nrm[3 * i] * d0[1] - nrm[3 * i + 1] * d0[0];
				for (k = 0; k < 3; k++) { a[k] = x[k] - h * d0[k]; b[k] = x[k] + h * d0[k]; }
				er_line(c, a, b, dark, 1.2);
				for (k = 0; k < 3; k++) { a[k] = x[k] - h * d1[k]; b[k] = x[k] + h * d1[k]; }
				er_line(c, a, b, dark, 1.2);
			}
			free(nrm);
		}
		qaws_mesh_3d_get_triangles(fv, dense);
		for (i = 0; i < nf; i++)
			if (idx[i])
			{
				double cen[3];
				unsigned int k;
				for (k = 0; k < 3; k++)
					cen[k] = (p[3 * fv[3 * i] + k] + p[3 * fv[3 * i + 1] + k] + p[3 * fv[3 * i + 2] + k]) / 3;
				er_point(c, cen, idx[i] > 0 ? red : blue, 7);
			}
		free(dir); free(idx); free(fv); free(p);
	}
	q1 = ms_quads(blob, 0.16, 0, &rep);
	if (q1)
	{
		ms_quad_sub(sub, q1, &rep);
		ms_panel(c, 550, 0, 550, 640, q1, 0.6, 0.5, "Instant quads", sub);
		ms_draw(c, q1, tint, NULL, 1, wire, 1.1);
		printf("  blob quads: %u quads, %u triangles, %u polygons, singular %u, deviation %.2e\n", rep.quality.quads, rep.quality.triangles,
			rep.quality.polygons, rep.quality.irregular, (double)rep.max_deviation);
	}
	q2 = ms_quads(blob, 0.11, 1, &rep);
	if (q2)
	{
		ms_quad_sub(sub, q2, &rep);
		ms_panel(c, 1100, 0, 550, 640, q2, 0.6, 0.5, "all quads", sub);
		ms_draw(c, q2, tint, NULL, 1, wire, 1.0);
		printf("  blob pure: %u quads, %u others, singular %u\n", rep.quality.quads, rep.quality.triangles + rep.quality.polygons, rep.quality.irregular);
	}
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	q3 = ms_quads(can, 0.14, 1, &rep);
	if (q3)
	{
		unsigned int nfe = qaws_mesh_3d_get_feature_edge_count(q3), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
		double* p = ms_positions(q3);
		ms_quad_sub(sub, q3, &rep);
		ms_panel(c, 1650, 0, 550, 640, q3, 0.35, 0.5, "along the creases", sub);
		ms_draw(c, q3, warm, NULL, 1, wire, 1.0);
		qaws_mesh_3d_get_feature_edges(fe, q3);
		for (i = 0; i < nfe; i++)
			er_line(c, p + 3 * fe[2 * i], p + 3 * fe[2 * i + 1], red, 2.4);
		printf("  can quads: %u quads, %u others, singular %u, creases %u\n", rep.quality.quads, rep.quality.triangles + rep.quality.polygons,
			rep.quality.irregular, nfe);
		free(fe);
		free(p);
	}
	q4 = ms_quads(torus, 0.3, 1, &rep);
	if (q4)
	{
		ms_quad_sub(sub, q4, &rep);
		ms_panel(c, 2200, 0, 550, 640, q4, 0.5, 0.75, "a torus", sub);
		ms_draw(c, q4, green, NULL, 1, wire, 1.0);
		printf("  torus quads: %u quads, %u others, singular %u\n", rep.quality.quads, rep.quality.triangles + rep.quality.polygons, rep.quality.irregular);
	}
	er_save(c, "showcase/mesh4_quads.png");
	printf("mesh4_quads.png\n");
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(torus);
	qaws_mesh_3d_destroy(dense);
	qaws_mesh_3d_destroy(q1);
	qaws_mesh_3d_destroy(q2);
	qaws_mesh_3d_destroy(q3);
	qaws_mesh_3d_destroy(q4);
	er_destroy(c);
}
