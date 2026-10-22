/*
 * Figure 43: tetrahedra keeping the border (qaws_tetrahedralize_desc
 * keep_border).
 *   - Schonhardt's polyhedron: a twisted prism cut along its reflex
 *     diagonals, which no tetrahedra of its own vertices fill;
 *   - the recovery alone: points added on the border (red), its triangles
 *     split; the cells drawn shrunk apart;
 *   - kept: the points moved inside (red), the border's 8 triangles as
 *     given;
 *   - a sphere's rings filled, its 224 triangles kept, cut open.
 */

static qaws_mesh_3d* ms_tb_prism(void)
{
	double x[18];
	static unsigned int const t[24] = { 0, 2, 1, 3, 4, 5, 0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4, 2, 0, 3, 2, 3, 5 };
	unsigned int i;
	for (i = 0; i < 3; i++)
	{
		double th = MS_PI / 2 + 2 * MS_PI * i / 3;
		x[3 * i] = cos(th);
		x[3 * i + 1] = sin(th);
		x[3 * i + 2] = 0;
		x[9 + 3 * i] = cos(th + MS_PI / 6);
		x[9 + 3 * i + 1] = sin(th + MS_PI / 6);
		x[9 + 3 * i + 2] = 1.2;
	}
	return ms_make(x, 6, t, NULL, 8);
}

/* the cells shrunk toward their centroids, each a colour; vertices not at one of the surface s's as red dots */
static void ms_tb_exploded(er_canvas* c, qaws_volume_mesh const* v, qaws_mesh_3d const* s, double const* wire)
{
	static unsigned char const F[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	static double const red[3] = { 214, 39, 40 };
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k, j;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	qaws_volume_mesh_get_positions(ps, v);
	qaws_volume_mesh_get_cells(cv, NULL, v);
	for (i = 0; i < nc; i++)
	{
		double p[4][3], cen[3] = { 0, 0, 0 }, col[3];
		for (k = 0; k < 4; k++)
			for (j = 0; j < 3; j++) cen[j] += (double)ps[3 * cv[4 * i + k] + j] / 4;
		for (k = 0; k < 4; k++)
			for (j = 0; j < 3; j++) p[k][j] = cen[j] + 0.55 * ((double)ps[3 * cv[4 * i + k] + j] - cen[j]);
		ms_colormap((double)i / (nc > 1 ? nc - 1 : 1), col);
		for (k = 0; k < 4; k++)
		{
			er_triangle(c, p[F[k][0]], p[F[k][1]], p[F[k][2]], NULL, NULL, NULL, col, col, col);
			er_line(c, p[F[k][0]], p[F[k][1]], wire, 0.7);
			er_line(c, p[F[k][1]], p[F[k][2]], wire, 0.7);
			er_line(c, p[F[k][2]], p[F[k][0]], wire, 0.7);
		}
	}
	for (i = 0; i < nv; i++)
	{
		double p[3] = { (double)ps[3 * i], (double)ps[3 * i + 1], (double)ps[3 * i + 2] };
		unsigned int n = qaws_mesh_3d_get_vertex_count(s), found = 0;
		for (k = 0; k < n; k++)
		{
			qaws_vec3 q = qaws_mesh_3d_get_vertex(s, k);
			found |= fabs((double)q.x - p[0]) + fabs((double)q.y - p[1]) + fabs((double)q.z - p[2]) < 1e-9;
		}
		if (!found)
			er_point(c, p, red, 5);
	}
	free(ps);
	free(cv);
}

static void demo_tetrahedra_border(void)
{
	static double const wire[3] = { 60, 66, 80 }, skin[3] = { 200, 205, 215 }, tint[3] = { 178, 192, 212 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	qaws_mesh_3d *p = ms_tb_prism(), *s = ms_rm_sphere(16, 8);
	qaws_tetrahedralize_desc d;
	qaws_mesh_report r;
	qaws_volume_mesh *a = NULL, *b = NULL, *v = NULL;
	char sub[200];
	double const yaw = 0.5, pitch = 0.8;
	if (!c || !p || !s)
		return;
	ms_panel(c, 0, 0, 450, 640, p, yaw, pitch, "Schonhardt", "8 triangles, 3 reflex edges");
	ms_draw(c, p, skin, NULL, 0, wire, 1.2);
	memset(&d, 0, sizeof d);
	d.border_only = 1;
	qaws_volume_mesh_tetrahedralize(&a, &r, p, &d);
	if (a)
	{
		sprintf(sub, "%u cells, %u points on the border", qaws_volume_mesh_get_cell_count(a), r.splits);
		ms_panel(c, 450, 0, 450, 640, p, yaw, pitch, "recovered", sub);
		ms_tb_exploded(c, a, p, wire);
		printf("  recovered: %s\n", sub);
	}
	d.keep_border = 1;
	qaws_volume_mesh_tetrahedralize(&b, &r, p, &d);
	if (b)
	{
		sprintf(sub, "%u cells, %u points moved inside", qaws_volume_mesh_get_cell_count(b), r.collapses);
		ms_panel(c, 900, 0, 450, 640, p, yaw, pitch, "border kept", sub);
		ms_tb_exploded(c, b, p, wire);
		printf("  kept: %s\n", sub);
	}
	d.border_only = 0;
	qaws_volume_mesh_tetrahedralize(&v, &r, s, &d);
	if (v)
	{
		qaws_volume_quality q;
		qaws_volume_mesh_measure_quality(&q, v);
		sprintf(sub, "224 triangles kept, %u points inside", r.collapses);
		ms_vpanel(c, 1350, 0, 450, 640, v, 0.4, 0.6, "a sphere's rings", sub);
		ms_cutaway(c, v, ms_keep_half, tint, wire);
		printf("  sphere: %s, dihedral %.1f..%.1f\n", sub, (double)q.min_dihedral, (double)q.max_dihedral);
	}
	qaws_volume_mesh_destroy(a);
	qaws_volume_mesh_destroy(b);
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(p);
	qaws_mesh_3d_destroy(s);
	er_save(c, "showcase/mesh43_tetrahedra_border.png");
	printf("mesh43_tetrahedra_border.png\n");
	er_destroy(c);
}
