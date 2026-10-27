/*
 * Figure 78: the generation methods beside REMESH.
 *   STRUCTURED, iso-parametric grids: the wavy patch in quads, the same to a
 *   chord tolerance (the grid's counts doubled until every cell's centre is
 *   within it), the NURBS torus (no irregular vertex), the box shell (each
 *   face a grid, the counts matched across its edges).
 *   PARAMETRIC, below: each face's grid remeshed in (u, v) in the surface's
 *   metric (its first fundamental form over the size squared), lifted: the
 *   torus, the box shell (faces meeting on the grids' border nodes), the
 *   patch to a chord tolerance (sizes from its principal curvatures);
 *   RESTRICTED_DELAUNAY, last: the same patch made Delaunay and refined
 *   (Chew): every angle at least 30 degrees.
 */

static qaws_mesh_3d* ms_gm_make(qaws_surface const* s, qaws_shell_atlas const* a, qaws_generate_method method, qaws_elements el, double L,
	double chord, qaws_mesh_report* rep)
{
	qaws_generate_desc d;
	qaws_mesh_3d* out = NULL;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.method = method;
	d.elements = el;
	d.target_length = (qaws_scalar)L;
	d.chord_tolerance = (qaws_scalar)chord;
	memset(rep, 0, sizeof *rep);
	st = s ? qaws_mesh_generate_surface(&out, rep, s, &d) : qaws_mesh_generate_shell(&out, rep, a, &d);
	if (st != QAWS_STATUS_OK)
		printf("  generate failed: %s\n", qaws_status_to_string(st));
	return out;
}

static void demo_generate_methods(void)
{
	static double const wire[3] = { 40, 45, 60 }, tint[3] = { 178, 192, 212 }, gold[3] = { 232, 196, 120 }, green[3] = { 170, 205, 160 };
	static double const pal[6][3] = { { 230, 97, 1 }, { 253, 184, 99 }, { 178, 171, 210 }, { 94, 60, 153 }, { 120, 198, 121 }, { 49, 130, 189 } };
	er_canvas* c = er_create(2200, 1280, 3, 252, 252, 250);
	qaws_surface *wavy = ms_wavy(), *torus = ms_nurbs_torus(2, 0.7), *box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d* g;
	qaws_mesh_report rep;
	char sub[200];
	unsigned int i;
	if (!c || !wavy || !torus)
		return;
	ms_box_faces(box, 1, 0.7, 0.5);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	if ((g = ms_gm_make(wavy, NULL, QAWS_GENERATE_STRUCTURED, QAWS_ELEMENTS_QUADS, 0.1, 0, &rep)) != NULL)
	{
		sprintf(sub, "%u quads, %u irregular", rep.quality.quads, rep.quality.irregular);
		ms_panel(c, 0, 0, 550, 640, g, 0.6, 0.65, "STRUCTURED patch", sub);
		ms_draw(c, g, tint, NULL, 1, wire, 0.9);
		printf("  structured patch: %u quads, edges %.3f..%.3f\n", rep.quality.quads, (double)rep.quality.min_edge, (double)rep.quality.max_edge);
	}
	qaws_mesh_3d_destroy(g);
	if ((g = ms_gm_make(wavy, NULL, QAWS_GENERATE_STRUCTURED, QAWS_ELEMENTS_TRIANGLES, 0, 0.01, &rep)) != NULL)
	{
		sprintf(sub, "chord 0.01: %u triangles", qaws_mesh_3d_get_face_count(g));
		ms_panel(c, 550, 0, 550, 640, g, 0.6, 0.65, "to a chord tolerance", sub);
		ms_draw(c, g, gold, NULL, 1, wire, 0.7);
		printf("  structured patch by chord: %u triangles, edges %.3f..%.3f\n", qaws_mesh_3d_get_face_count(g), (double)rep.quality.min_edge,
			(double)rep.quality.max_edge);
	}
	qaws_mesh_3d_destroy(g);
	if ((g = ms_gm_make(torus, NULL, QAWS_GENERATE_STRUCTURED, QAWS_ELEMENTS_QUADS, 0.25, 0, &rep)) != NULL)
	{
		sprintf(sub, "%u quads, %u irregular", rep.quality.quads, rep.quality.irregular);
		ms_panel(c, 1100, 0, 550, 640, g, 0.5, 0.75, "the NURBS torus", sub);
		ms_draw(c, g, green, NULL, 1, wire, 0.9);
		printf("  structured torus: %u quads, %u irregular\n", rep.quality.quads, rep.quality.irregular);
	}
	qaws_mesh_3d_destroy(g);
	if ((g = atlas ? ms_gm_make(NULL, atlas, QAWS_GENERATE_STRUCTURED, QAWS_ELEMENTS_QUADS, 0.15, 0, &rep) : NULL) != NULL)
	{
		unsigned int nf = qaws_mesh_3d_get_face_count(g), f, a, q = 0;
		qaws_scalar* fid = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
		unsigned int* fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
		unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(g));
		double* p = ms_positions(g);
		qaws_mesh_3d_find_attribute(&a, g, "face");
		qaws_mesh_3d_get_attribute(fid, g, a);
		qaws_mesh_3d_get_faces(fv, fs, g);
		sprintf(sub, "6 faces: %u quads, %u irregular", rep.quality.quads, rep.quality.irregular);
		ms_panel(c, 1650, 0, 550, 640, g, 0.62, 0.55, "the box shell", sub);
		for (f = 0; f < nf; f++)
		{
			unsigned int k, id = (unsigned int)fid[f] % 6;
			for (k = 1; k + 1 < fs[f]; k++)
				er_triangle(c, p + 3 * fv[q], p + 3 * fv[q + k], p + 3 * fv[q + k + 1], NULL, NULL, NULL, pal[id], pal[id], pal[id]);
			q += fs[f];
		}
		ms_draw(c, g, NULL, NULL, 0, wire, 1.0);
		printf("  structured box: %u quads, %u irregular\n", rep.quality.quads, rep.quality.irregular);
		free(fid); free(fs); free(fv); free(p);
	}
	qaws_mesh_3d_destroy(g);
	if ((g = ms_gm_make(torus, NULL, QAWS_GENERATE_PARAMETRIC, QAWS_ELEMENTS_TRIANGLES, 0.25, 0, &rep)) != NULL)
	{
		sprintf(sub, "%u triangles, angles %.0f..%.0f", qaws_mesh_3d_get_face_count(g), (double)rep.quality.min_angle, (double)rep.quality.max_angle);
		ms_panel(c, 0, 640, 550, 640, g, 0.5, 0.75, "PARAMETRIC torus", sub);
		ms_draw(c, g, green, NULL, 1, wire, 0.9);
		printf("  parametric torus: %u triangles, angles %.1f..%.1f, edges %.3f..%.3f\n", qaws_mesh_3d_get_face_count(g), (double)rep.quality.min_angle,
			(double)rep.quality.max_angle, (double)rep.quality.min_edge, (double)rep.quality.max_edge);
	}
	qaws_mesh_3d_destroy(g);
	if ((g = atlas ? ms_gm_make(NULL, atlas, QAWS_GENERATE_PARAMETRIC, QAWS_ELEMENTS_TRIANGLES, 0.15, 0, &rep) : NULL) != NULL)
	{
		sprintf(sub, "6 faces: %u triangles, angles %.0f..%.0f", qaws_mesh_3d_get_face_count(g), (double)rep.quality.min_angle, (double)rep.quality.max_angle);
		ms_panel(c, 550, 640, 550, 640, g, 0.62, 0.55, "the box shell", sub);
		ms_draw(c, g, tint, NULL, 0, wire, 1.0);
		printf("  parametric box: %u triangles, angles %.1f..%.1f\n", qaws_mesh_3d_get_face_count(g), (double)rep.quality.min_angle, (double)rep.quality.max_angle);
	}
	qaws_mesh_3d_destroy(g);
	{
		qaws_generate_desc d;
		g = NULL;
		memset(&d, 0, sizeof d);
		d.method = QAWS_GENERATE_PARAMETRIC;
		d.chord_tolerance = (qaws_scalar)0.003;
		d.min_length = (qaws_scalar)0.02;
		d.max_length = (qaws_scalar)0.25;
		if (qaws_mesh_generate_surface(&g, &rep, wavy, &d) == QAWS_STATUS_OK && g)
		{
			sprintf(sub, "chord 0.003: %u triangles", qaws_mesh_3d_get_face_count(g));
			ms_panel(c, 1100, 640, 550, 640, g, 0.6, 0.65, "adaptive by curvature", sub);
			ms_draw(c, g, gold, NULL, 1, wire, 0.8);
			printf("  parametric patch by chord: %u triangles, edges %.3f..%.3f, angles %.1f..%.1f\n", qaws_mesh_3d_get_face_count(g),
				(double)rep.quality.min_edge, (double)rep.quality.max_edge, (double)rep.quality.min_angle, (double)rep.quality.max_angle);
		}
		qaws_mesh_3d_destroy(g);
	}
	{
		qaws_generate_desc d;
		g = NULL;
		memset(&d, 0, sizeof d);
		d.method = QAWS_GENERATE_RESTRICTED_DELAUNAY;
		d.chord_tolerance = (qaws_scalar)0.003;
		d.min_length = (qaws_scalar)0.02;
		d.max_length = (qaws_scalar)0.25;
		if (qaws_mesh_generate_surface(&g, &rep, wavy, &d) == QAWS_STATUS_OK && g)
		{
			sprintf(sub, "%u triangles, angles %.1f..%.1f", qaws_mesh_3d_get_face_count(g), (double)rep.quality.min_angle, (double)rep.quality.max_angle);
			ms_panel(c, 1650, 640, 550, 640, g, 0.6, 0.65, "RESTRICTED_DELAUNAY", sub);
			ms_draw(c, g, tint, NULL, 1, wire, 0.8);
			printf("  delaunay patch: %u triangles, angles %.1f..%.1f\n", qaws_mesh_3d_get_face_count(g), (double)rep.quality.min_angle,
				(double)rep.quality.max_angle);
		}
		qaws_mesh_3d_destroy(g);
	}
	er_save(c, "showcase/mesh78_generate_methods.png");
	printf("mesh78_generate_methods.png\n");
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
	qaws_surface_destroy(wavy);
	qaws_surface_destroy(torus);
	er_destroy(c);
}
