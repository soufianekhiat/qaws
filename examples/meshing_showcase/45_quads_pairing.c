/*
 * Figure 45: quads from the input's own faces.
 *   - figure 2's blob, triangles;
 *   - PAIRING: its triangles paired across their sides into quads (the ones
 *     no pair takes in red);
 *   - PAIRING, pure: the pairs split round their centres, all quads;
 *   - SPLIT: every triangle into three quads on the reference.
 */

/* a mesh's faces, quads in rgb, other faces in red, flat shaded, outlined */
static void ms_qp_draw(er_canvas* c, qaws_mesh_3d const* m, double const* rgb, double const* wire)
{
	static double const red[3] = { 214, 60, 60 };
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), f, q = 0, k;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	double* x = ms_positions(m);
	qaws_mesh_3d_get_faces(fv, fs, m);
	for (f = 0; f < nf; q += fs[f++])
	{
		double const* col = fs[f] == 4 ? rgb : red;
		for (k = 1; k + 1 < fs[f]; k++)
			er_triangle(c, x + 3 * fv[q], x + 3 * fv[q + k], x + 3 * fv[q + k + 1], NULL, NULL, NULL, col, col, col);
		for (k = 0; k < fs[f]; k++)
			er_line(c, x + 3 * fv[q + k], x + 3 * fv[q + (k + 1) % fs[f]], wire, 0.6);
	}
	free(fv);
	free(fs);
	free(x);
}

static void demo_quads_pairing(void)
{
	static double const wire[3] = { 50, 55, 70 }, skin[3] = { 200, 205, 215 }, blue[3] = { 140, 170, 220 }, green[3] = { 150, 200, 160 }, tan[3] = { 225, 200, 150 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(2), *a = NULL, *b = NULL, *s = NULL;
	qaws_quadrangulate_desc d;
	qaws_mesh_report r;
	char sub[200];
	double const yaw = 0.5, pitch = 0.45;
	if (!c || !blob)
		return;
	sprintf(sub, "%u triangles", qaws_mesh_3d_get_face_count(blob));
	ms_panel(c, 0, 0, 450, 640, blob, yaw, pitch, "triangles", sub);
	ms_draw(c, blob, skin, NULL, 0, wire, 0.6);
	memset(&d, 0, sizeof d);
	d.method = QAWS_QUADS_PAIRING;
	if (qaws_mesh_quadrangulate(&a, &r, blob, &d) == QAWS_STATUS_OK)
	{
		sprintf(sub, "%u quads, %u triangles left", r.quality.quads, r.quality.triangles);
		ms_panel(c, 450, 0, 450, 640, a, yaw, pitch, "paired", sub);
		ms_qp_draw(c, a, blue, wire);
		printf("  paired: %s\n", sub);
	}
	d.pure = 1;
	if (qaws_mesh_quadrangulate(&b, &r, blob, &d) == QAWS_STATUS_OK)
	{
		sprintf(sub, "%u quads, angles %.0f..%.0f", r.quality.quads, (double)r.quality.min_angle, (double)r.quality.max_angle);
		ms_panel(c, 900, 0, 450, 640, b, yaw, pitch, "paired, pure", sub);
		ms_qp_draw(c, b, green, wire);
		printf("  pure: %s\n", sub);
	}
	d.method = QAWS_QUADS_SPLIT;
	d.pure = 0;
	if (qaws_mesh_quadrangulate(&s, &r, blob, &d) == QAWS_STATUS_OK)
	{
		sprintf(sub, "%u quads, angles %.0f..%.0f", r.quality.quads, (double)r.quality.min_angle, (double)r.quality.max_angle);
		ms_panel(c, 1350, 0, 450, 640, s, yaw, pitch, "split", sub);
		ms_qp_draw(c, s, tan, wire);
		printf("  split: %s\n", sub);
	}
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(s);
	er_save(c, "showcase/mesh45_quads_pairing.png");
	printf("mesh45_quads_pairing.png\n");
	er_destroy(c);
}
