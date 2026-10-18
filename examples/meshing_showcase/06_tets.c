/*
 * Figure 6: tetrahedra (qaws_volume_mesh_tetrahedralize), cut open.
 *   - a NURBS torus's surface mesh filled: half its cells shown;
 *   - a box of six faces (a shell's triangles): a corner taken out;
 *   - a blob remeshed on its PN reference: half its cells shown.
 * The cut's faces are coloured by their cell's smallest dihedral angle
 * (blue below 10 degrees to red at 70), the border's in plain colour.
 */

static int ms_keep_half(double const* p)
{
	return p[1] > 0.0;
}

static int ms_keep_box(double const* p)
{
	return !(p[0] > 0 && p[1] < 0 && p[2] > 0);
}

/* the faces of the kept cells that no other kept cell shares, drawn */
static void ms_cutaway(er_canvas* c, qaws_volume_mesh const* v, int (*keep)(double const*), double const* rgb, double const* wire)
{
	static unsigned char const F[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k, j;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	unsigned char* kept = (unsigned char*)malloc(nc);
	qaws_volume_mesh_get_positions(v, ps);
	for (i = 0; i < 3 * nv; i++)
		x[i] = (double)ps[i];
	qaws_volume_mesh_get_cells(v, cv, NULL);
	qaws_volume_mesh_get_cell_adjacency(v, adj);
	for (i = 0; i < nc; i++)
	{
		double cen[3] = { 0, 0, 0 };
		for (k = 0; k < 4; k++)
			for (j = 0; j < 3; j++)
				cen[j] += x[3 * cv[4 * i + k] + j] / 4;
		kept[i] = (unsigned char)keep(cen);
	}
	for (i = 0; i < nc; i++)
	{
		double nrm[4][3], mind = 180, col[3];
		if (!kept[i])
			continue;
		for (k = 0; k < 4; k++)
		{
			double const *a = x + 3 * cv[4 * i + F[k][0]], *b = x + 3 * cv[4 * i + F[k][1]], *d = x + 3 * cv[4 * i + F[k][2]];
			double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, w[3] = { d[0] - a[0], d[1] - a[1], d[2] - a[2] };
			nrm[k][0] = u[1] * w[2] - u[2] * w[1];
			nrm[k][1] = u[2] * w[0] - u[0] * w[2];
			nrm[k][2] = u[0] * w[1] - u[1] * w[0];
		}
		for (k = 0; k < 4; k++)
			for (j = k + 1; j < 4; j++)
			{
				double cs = -(nrm[k][0] * nrm[j][0] + nrm[k][1] * nrm[j][1] + nrm[k][2] * nrm[j][2])
					/ sqrt((nrm[k][0] * nrm[k][0] + nrm[k][1] * nrm[k][1] + nrm[k][2] * nrm[k][2]) * (nrm[j][0] * nrm[j][0] + nrm[j][1] * nrm[j][1] + nrm[j][2] * nrm[j][2]));
				mind = fmin(mind, acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180 / MS_PI);
			}
		ms_colormap((mind - 10) / 60, col);
		for (k = 0; k < 4; k++)
		{
			unsigned int o = adj[4 * i + k];
			double const *a = x + 3 * cv[4 * i + F[k][0]], *b = x + 3 * cv[4 * i + F[k][1]], *d = x + 3 * cv[4 * i + F[k][2]];
			if (o != QAWS_MESH_NONE && kept[o])
				continue;
			er_triangle(c, a, b, d, NULL, NULL, NULL, o == QAWS_MESH_NONE ? rgb : col, o == QAWS_MESH_NONE ? rgb : col, o == QAWS_MESH_NONE ? rgb : col);
			er_line(c, a, b, wire, 0.8);
			er_line(c, b, d, wire, 0.8);
			er_line(c, d, a, wire, 0.8);
		}
	}
	free(ps); free(x); free(cv); free(adj); free(kept);
}

static qaws_volume_mesh* ms_tets(qaws_mesh_3d const* b, char* sub, char const* what)
{
	qaws_tetrahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_volume_quality q;
	qaws_status st;
	memset(&d, 0, sizeof d);
	st = qaws_volume_mesh_tetrahedralize(b, &d, &v, NULL);
	if (st != QAWS_STATUS_OK)
	{
		printf("  tetrahedralize failed: %s\n", qaws_status_to_string(st));
		return NULL;
	}
	qaws_volume_mesh_measure_quality(v, &q);
	sprintf(sub, "%u tetrahedra, dihedral %.0f..%.0f", q.tetrahedra, (double)q.min_dihedral, (double)q.max_dihedral);
	printf("  %s: %u tetrahedra, %u vertices, dihedral %.1f..%.1f, radius-edge %.2f, volume %.6f, %s\n", what, q.tetrahedra,
		qaws_volume_mesh_get_vertex_count(v), (double)q.min_dihedral, (double)q.max_dihedral, (double)q.max_radius_edge, (double)q.volume,
		qaws_status_to_string(qaws_volume_mesh_validate(v)));
	return v;
}

/* a panel framed on a volume mesh's vertices */
static void ms_vpanel(er_canvas* c, int x, int y, int w, int h, qaws_volume_mesh const* v, double yaw, double pitch, char const* title, char const* sub)
{
	unsigned int n = qaws_volume_mesh_get_vertex_count(v), i;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	double* p = (double*)malloc(sizeof(double) * 3 * n);
	qaws_volume_mesh_get_positions(v, ps);
	for (i = 0; i < 3 * n; i++)
		p[i] = (double)ps[i];
	er_panel(c, x, y + 72, w, h - 72, p, n, yaw, pitch, 0.06);
	er_text(c, x + 14, y + 12, title, 3, 30, 30, 30);
	er_text(c, x + 14, y + 46, sub, 2, 95, 95, 95);
	free(ps);
	free(p);
}

static void demo_tets(void)
{
	static double const wire[3] = { 40, 45, 60 }, green[3] = { 170, 205, 160 }, tan[3] = { 232, 196, 120 }, tint[3] = { 178, 192, 212 };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	qaws_surface *torus = ms_nurbs_torus(2, 0.7), *box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d *ts = NULL, *bs = NULL, *blob = ms_blob(1), *rblob = NULL;
	qaws_mesh_reference* rb = NULL;
	qaws_volume_mesh *v1 = NULL, *v2 = NULL, *v3 = NULL;
	char sub[200];
	if (!c || !torus || !blob)
		return;
	ts = ms_generate(torus, NULL, QAWS_ELEMENTS_TRIANGLES, 0.25, 0, NULL);
	v1 = ts ? ms_tets(ts, sub, "torus") : NULL;
	if (v1)
	{
		ms_vpanel(c, 0, 0, 733, 700, v1, 0.5, 0.75, "a NURBS torus, filled", sub);
		ms_cutaway(c, v1, ms_keep_half, green, wire);
	}
	ms_box_faces(box, 1, 0.7, 0.5);
	qaws_shell_atlas_create((qaws_surface const* const*)box, 6, 0, &atlas);
	bs = atlas ? ms_generate(NULL, atlas, QAWS_ELEMENTS_TRIANGLES, 0.15, 0, NULL) : NULL;
	v2 = bs ? ms_tets(bs, sub, "box") : NULL;
	if (v2)
	{
		ms_vpanel(c, 733, 0, 733, 700, v2, 0.62, 0.55, "a box, a corner out", sub);
		ms_cutaway(c, v2, ms_keep_box, tan, wire);
	}
	qaws_mesh_reference_create(blob, NULL, &rb);
	rblob = rb ? ms_remesh(blob, rb, 0.12, NULL, NULL) : NULL;
	v3 = rblob ? ms_tets(rblob, sub, "blob") : NULL;
	if (v3)
	{
		ms_vpanel(c, 1466, 0, 733, 700, v3, 0.4, 0.6, "a blob, remeshed, filled", sub);
		ms_cutaway(c, v3, ms_keep_half, tint, wire);
	}
	er_save(c, "showcase/mesh6_tets.png");
	printf("mesh6_tets.png\n");
	qaws_volume_mesh_destroy(v1);
	qaws_volume_mesh_destroy(v2);
	qaws_volume_mesh_destroy(v3);
	qaws_mesh_3d_destroy(ts);
	qaws_mesh_3d_destroy(bs);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(rblob);
	qaws_mesh_reference_destroy(rb);
	qaws_shell_atlas_destroy(atlas);
	qaws_surface_destroy(torus);
	{
		unsigned int i;
		for (i = 0; i < 6; i++)
			qaws_surface_destroy(box[i]);
	}
	er_destroy(c);
}
