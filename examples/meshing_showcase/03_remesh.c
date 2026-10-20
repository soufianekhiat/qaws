/*
 * Figure 3: remeshing on the reference.
 *   - the coarse blob, and remeshed uniformly on its PN reference: smooth,
 *     the vertices on the surface the coarse mesh implies;
 *   - a bumpy sphere remeshed adaptively by a curvature size field (chord
 *     tolerance, gradation): dense on the bump, coarse where flat;
 *   - the coarse can with its creases, remeshed: rims kept sharp;
 *   - a cube with a colour per face (corner attributes with seams on its
 *     edges) remeshed: the colours transferred, the seams kept.
 */

static qaws_mesh_3d* ms_remesh(qaws_mesh_3d const* m, qaws_mesh_reference const* r, double L, qaws_size_field const* sf, qaws_mesh_report* rep)
{
	qaws_remesh_desc d;
	qaws_mesh_3d* out = NULL;
	memset(&d, 0, sizeof d);
	d.method = QAWS_REMESH_LOCAL;
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	d.size = sf;
	d.iterations = 12;
	if (qaws_mesh_remesh(&out, rep, m, &d) != QAWS_STATUS_OK)
		printf("  remesh failed\n");
	return out;
}

/* a sphere with a sharp bump and a dent: curvature from flat to high */
static qaws_mesh_3d* ms_bumpy(void)
{
	qaws_mesh_3d* m = ms_blob(2);
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i;
	double* p = ms_positions(m);
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	for (i = 0; i < n; i++)
	{
		double* q = p + 3 * i, l = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]), b, d, s;
		q[0] /= l; q[1] /= l; q[2] /= l;
		b = (q[0] - 0.3) * (q[0] - 0.3) + (q[1] + 0.2) * (q[1] + 0.2) + (q[2] - 0.93) * (q[2] - 0.93);
		d = (q[0] + 0.75) * (q[0] + 0.75) + (q[1] + 0.55) * (q[1] + 0.55) + (q[2] - 0.35) * (q[2] - 0.35);
		s = 1 + 0.65 * exp(-9 * b) - 0.22 * exp(-10 * d);
		ps[3 * i] = (qaws_scalar)(q[0] * s);
		ps[3 * i + 1] = (qaws_scalar)(q[1] * s);
		ps[3 * i + 2] = (qaws_scalar)(q[2] * s);
	}
	qaws_mesh_3d_set_positions(m, ps);
	free(p);
	free(ps);
	return m;
}

static void demo_remesh(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, warm[3] = { 226, 200, 160 }, red[3] = { 214, 39, 40 };
	er_canvas* c = er_create(2750, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *bumpy = ms_bumpy(), *can = ms_can(8, 1, 1.6), *uni = NULL, *ada = NULL, *rcan = NULL, *cube = NULL, *rcube = NULL;
	qaws_mesh_reference *rb = NULL, *rc = NULL, *ry = NULL;
	qaws_size_field* sf = NULL;
	qaws_mesh_report rep;
	qaws_curvature_size_desc cd;
	char sub[200];
	unsigned int i;
	if (!c || !blob || !can)
		return;
	qaws_mesh_reference_create(&rb, blob, NULL);
	ms_panel(c, 0, 0, 550, 640, blob, 0.6, 0.5, "the coarse blob", "80 triangles");
	ms_draw(c, blob, tint, NULL, 0, wire, 1.6);
	uni = ms_remesh(blob, rb, 0.11, NULL, &rep);
	if (uni)
	{
		sprintf(sub, "%u triangles, min angle %.0f, off by %.1e", qaws_mesh_3d_get_face_count(uni), (double)rep.quality.min_angle, (double)rep.max_deviation);
		ms_panel(c, 550, 0, 550, 640, uni, 0.6, 0.5, "remeshed on PN", sub);
		ms_draw(c, uni, tint, NULL, 1, wire, 0.9);
		printf("  uniform: %u triangles, %u splits %u collapses %u flips %u moves, angles %.1f..%.1f, edges %.3f..%.3f mean %.3f, irregular %u, deviation %.2e\n",
			qaws_mesh_3d_get_face_count(uni), rep.splits, rep.collapses, rep.flips, rep.moves, (double)rep.quality.min_angle, (double)rep.quality.max_angle,
			(double)rep.quality.min_edge, (double)rep.quality.max_edge, (double)rep.quality.mean_edge, rep.quality.irregular, (double)rep.max_deviation);
	}
	memset(&cd, 0, sizeof cd);
	cd.angle_tolerance = (qaws_scalar)(12 * MS_PI / 180);
	cd.min_length = (qaws_scalar)0.02;
	cd.max_length = (qaws_scalar)0.3;
	cd.gradation = (qaws_scalar)0.35;
	qaws_mesh_reference_create(&ry, bumpy, NULL);
	if (ry)
		qaws_size_field_create_curvature(&sf, ry, &cd);
	ada = sf ? ms_remesh(bumpy, ry, 0, sf, &rep) : NULL;
	if (ada)
	{
		sprintf(sub, "normals turn < 12 deg: %u triangles", qaws_mesh_3d_get_face_count(ada));
		ms_panel(c, 1100, 0, 550, 640, ada, 0.9, 0.35, "adaptive by curvature", sub);
		ms_draw(c, ada, tint, NULL, 1, wire, 0.9);
		printf("  adaptive: %u triangles, angles %.1f..%.1f, edges %.3f..%.3f, deviation %.2e\n", qaws_mesh_3d_get_face_count(ada),
			(double)rep.quality.min_angle, (double)rep.quality.max_angle, (double)rep.quality.min_edge, (double)rep.quality.max_edge, (double)rep.max_deviation);
	}
	/* the can, its creases kept */
	qaws_mesh_3d_detect_features(can, (qaws_scalar)(50 * MS_PI / 180), 0, 0);
	qaws_mesh_reference_create(&rc, can, NULL);
	rcan = rc ? ms_remesh(can, rc, 0.12, NULL, &rep) : NULL;
	if (rcan)
	{
		unsigned int nfe = qaws_mesh_3d_get_feature_edge_count(rcan), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfe + 1));
		double* p = ms_positions(rcan);
		sprintf(sub, "%u triangles, %u crease edges", qaws_mesh_3d_get_face_count(rcan), nfe);
		ms_panel(c, 1650, 0, 550, 640, rcan, 0.35, 0.5, "the can remeshed", sub);
		ms_draw(c, rcan, warm, NULL, 1, wire, 0.9);
		qaws_mesh_3d_get_feature_edges(fe, rcan);
		for (i = 0; i < nfe; i++)
			er_line(c, p + 3 * fe[2 * i], p + 3 * fe[2 * i + 1], red, 2.4);
		printf("  can: %u triangles, %u creases, angles %.1f..%.1f, deviation %.2e\n", qaws_mesh_3d_get_face_count(rcan), nfe,
			(double)rep.quality.min_angle, (double)rep.quality.max_angle, (double)rep.max_deviation);
		free(fe);
		free(p);
	}
	/* a cube coloured per face (corners), its seams on the creases */
	{
		static double const x[24] = { -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1 };
		static unsigned int const f[24] = { 0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7 };
		static unsigned int const fs[6] = { 4, 4, 4, 4, 4, 4 };
		static double const pal[6][3] = { { 230, 97, 1 }, { 253, 184, 99 }, { 178, 171, 210 }, { 94, 60, 153 }, { 120, 198, 121 }, { 49, 130, 189 } };
		qaws_scalar cols[72];
		qaws_attribute_desc ad;
		unsigned int a, nc;
		cube = ms_make(x, 8, f, fs, 6);
		for (i = 0; i < 24; i++)
		{
			cols[3 * i] = (qaws_scalar)(pal[i / 4][0] / 255);
			cols[3 * i + 1] = (qaws_scalar)(pal[i / 4][1] / 255);
			cols[3 * i + 2] = (qaws_scalar)(pal[i / 4][2] / 255);
		}
		memset(&ad, 0, sizeof ad);
		ad.name = "colour";
		ad.element = QAWS_MESH_CORNER;
		ad.components = 3;
		ad.rule = QAWS_ATTRIBUTE_LINEAR;
		qaws_mesh_3d_add_attribute(&a, cube, cols, &ad);
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)(30 * MS_PI / 180), 0, 0);
		rcube = ms_remesh(cube, NULL, 0.22, NULL, &rep);
		if (rcube && qaws_mesh_3d_find_attribute(&a, rcube, "colour") == QAWS_STATUS_OK)
		{
			/* corners coloured: drawn per face from its corners */
			unsigned int nf = qaws_mesh_3d_get_face_count(rcube), *fv, *fz, k, q;
			qaws_scalar* cc;
			double* p = ms_positions(rcube);
			nc = qaws_mesh_3d_get_corner_count(rcube);
			fv = (unsigned int*)malloc(sizeof(unsigned int) * nc);
			fz = (unsigned int*)malloc(sizeof(unsigned int) * nf);
			cc = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nc);
			qaws_mesh_3d_get_faces(fv, fz, rcube);
			qaws_mesh_3d_get_attribute(cc, rcube, a);
			sprintf(sub, "%u triangles: colour seams kept", nf);
			ms_panel(c, 2200, 0, 550, 640, rcube, 0.62, 0.55, "a corner attribute", sub);
			for (k = 0, q = 0; k < nf; k++, q += fz[k - 1])
			{
				double c0[3], c1[3], c2[3];
				unsigned int j;
				for (j = 0; j < 3; j++)
				{
					c0[j] = 255 * cc[3 * q + j];
					c1[j] = 255 * cc[3 * (q + 1) + j];
					c2[j] = 255 * cc[3 * (q + 2) + j];
				}
				er_triangle(c, p + 3 * fv[q], p + 3 * fv[q + 1], p + 3 * fv[q + 2], NULL, NULL, NULL, c0, c1, c2);
			}
			ms_draw(c, rcube, NULL, NULL, 0, wire, 0.9);
			free(fv); free(fz); free(cc); free(p);
		}
	}
	er_save(c, "showcase/mesh3_remesh.png");
	printf("mesh3_remesh.png\n");
	qaws_size_field_destroy(sf);
	qaws_mesh_reference_destroy(rb);
	qaws_mesh_reference_destroy(rc);
	qaws_mesh_reference_destroy(ry);
	qaws_mesh_3d_destroy(bumpy);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(can);
	qaws_mesh_3d_destroy(uni);
	qaws_mesh_3d_destroy(ada);
	qaws_mesh_3d_destroy(rcan);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(rcube);
	er_destroy(c);
}
