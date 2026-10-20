/*
 * Figure 7: hexahedra and frame fields.
 *   - a box of six faces in hexahedra (GRID), a corner taken out: the
 *     staircase of cubes and the layer padding it to the surface, sharp
 *     edges kept;
 *   - a NURBS torus in hexahedra, half shown;
 *   - the frame field of the torus's tetrahedra: the singular edges (red),
 *     where hexahedra following the field would meet three or five round an
 *     edge, the frames in the torus's middle plane (their three axes).
 * The cut faces are coloured by their cell's smallest scaled Jacobian (blue
 * at 0.2 to red at 1).
 */

/* the faces of the kept cells (hexahedra) that no other kept cell shares, drawn */
static void ms_hex_cutaway(er_canvas* c, qaws_volume_mesh const* v, int (*keep)(double const*), double const* rgb, double const* wire)
{
	static unsigned char const F[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	static unsigned char const nb[8][3] = { { 1, 3, 4 }, { 2, 0, 5 }, { 3, 1, 6 }, { 0, 2, 7 }, { 7, 5, 0 }, { 4, 6, 1 }, { 5, 7, 2 }, { 6, 4, 3 } };
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k, j;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 8 * nc);
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nc);
	unsigned char* kept = (unsigned char*)malloc(nc);
	qaws_volume_mesh_get_positions(ps, v);
	for (i = 0; i < 3 * nv; i++)
		x[i] = (double)ps[i];
	qaws_volume_mesh_get_cells(cv, NULL, v);
	qaws_volume_mesh_get_cell_adjacency(adj, v);
	for (i = 0; i < nc; i++)
	{
		double cen[3] = { 0, 0, 0 };
		for (k = 0; k < 8; k++)
			for (j = 0; j < 3; j++)
				cen[j] += x[3 * cv[8 * i + k] + j] / 8;
		kept[i] = (unsigned char)keep(cen);
	}
	for (i = 0; i < nc; i++)
	{
		double sj = 1, col[3];
		if (!kept[i])
			continue;
		for (k = 0; k < 8; k++)
		{
			double const* p = x + 3 * cv[8 * i + k];
			double a[3], b[3], d[3], l;
			for (j = 0; j < 3; j++)
			{
				a[j] = x[3 * cv[8 * i + nb[k][0]] + j] - p[j];
				b[j] = x[3 * cv[8 * i + nb[k][1]] + j] - p[j];
				d[j] = x[3 * cv[8 * i + nb[k][2]] + j] - p[j];
			}
			l = sqrt((a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) * (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
			sj = fmin(sj, (a[0] * (b[1] * d[2] - b[2] * d[1]) + a[1] * (b[2] * d[0] - b[0] * d[2]) + a[2] * (b[0] * d[1] - b[1] * d[0])) / l);
		}
		ms_colormap((sj - 0.2) / 0.8, col);
		for (k = 0; k < 6; k++)
		{
			unsigned int o = adj[6 * i + k];
			double const *p0 = x + 3 * cv[8 * i + F[k][0]], *p1 = x + 3 * cv[8 * i + F[k][1]], *p2 = x + 3 * cv[8 * i + F[k][2]], *p3 = x + 3 * cv[8 * i + F[k][3]];
			double const* fc = o == QAWS_MESH_NONE ? rgb : col;
			if (o != QAWS_MESH_NONE && kept[o])
				continue;
			er_triangle(c, p0, p1, p2, NULL, NULL, NULL, fc, fc, fc);
			er_triangle(c, p0, p2, p3, NULL, NULL, NULL, fc, fc, fc);
			er_line(c, p0, p1, wire, 0.8);
			er_line(c, p1, p2, wire, 0.8);
			er_line(c, p2, p3, wire, 0.8);
			er_line(c, p3, p0, wire, 0.8);
		}
	}
	free(ps); free(x); free(cv); free(adj); free(kept);
}

static qaws_volume_mesh* ms_hexes(qaws_mesh_3d const* b, double h, char* sub, char const* what)
{
	qaws_hexahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_volume_quality q;
	qaws_status st;
	memset(&d, 0, sizeof d);
	d.cell_size = (qaws_scalar)h;
	st = qaws_volume_mesh_hexahedralize(&v, NULL, b, &d);
	if (st != QAWS_STATUS_OK)
	{
		printf("  hexahedralize failed: %s\n", qaws_status_to_string(st));
		return NULL;
	}
	qaws_volume_mesh_measure_quality(&q, v);
	sprintf(sub, "%u hexahedra, scaled Jacobian from %.2f", q.hexahedra, (double)q.min_scaled_jacobian);
	printf("  %s: %u hexahedra, %u vertices, scaled Jacobian %.3f, volume %.6f, %s\n", what, q.hexahedra, qaws_volume_mesh_get_vertex_count(v),
		(double)q.min_scaled_jacobian, (double)q.volume, qaws_status_to_string(qaws_volume_mesh_validate(v)));
	return v;
}

static void demo_hexes(void)
{
	static double const wire[3] = { 40, 45, 60 }, tan[3] = { 232, 196, 120 }, green[3] = { 170, 205, 160 }, grey[3] = { 175, 180, 190 };
	static double const red[3] = { 200, 30, 40 }, axis[3][3] = { { 220, 90, 60 }, { 70, 160, 80 }, { 60, 100, 210 } };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	qaws_surface *torus = ms_nurbs_torus(2, 0.7), *box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d *ts = NULL, *bs = NULL;
	qaws_volume_mesh *h1 = NULL, *h2 = NULL, *tt = NULL;
	char sub[200];
	if (!c || !torus)
		return;
	ms_box_faces(box, 1, 0.7, 0.5);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	bs = atlas ? ms_generate(NULL, atlas, QAWS_ELEMENTS_TRIANGLES, 0.1, 0, NULL) : NULL;
	h1 = bs ? ms_hexes(bs, 0.15, sub, "box") : NULL;
	if (h1)
	{
		ms_vpanel(c, 0, 0, 733, 700, h1, 0.62, 0.55, "a box in hexahedra", sub);
		ms_hex_cutaway(c, h1, ms_keep_box, tan, wire);
	}
	ts = ms_generate(torus, NULL, QAWS_ELEMENTS_TRIANGLES, 0.2, 0, NULL);
	h2 = ts ? ms_hexes(ts, 0.22, sub, "torus") : NULL;
	if (h2)
	{
		ms_vpanel(c, 733, 0, 733, 700, h2, 0.5, 0.75, "a torus in hexahedra", sub);
		ms_hex_cutaway(c, h2, ms_keep_half, green, wire);
	}
	/* the frame field of the torus's tetrahedra */
	if (ts)
	{
		qaws_tetrahedralize_desc td;
		memset(&td, 0, sizeof td);
		qaws_volume_mesh_tetrahedralize(&tt, NULL, ts, &td);
	}
	if (tt)
	{
		unsigned int nv = qaws_volume_mesh_get_vertex_count(tt), nc = qaws_volume_mesh_get_cell_count(tt), i, k, ns = 0;
		static unsigned char const E[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } };
		qaws_scalar *ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *fr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * nc);
		double* x = (double*)malloc(sizeof(double) * 3 * nv);
		unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
		unsigned char* sing = (unsigned char*)malloc(nc);
		qaws_mesh_3d* tb = NULL;
		qaws_volume_mesh_get_positions(ps, tt);
		for (i = 0; i < 3 * nv; i++)
			x[i] = (double)ps[i];
		qaws_volume_mesh_get_cells(cv, NULL, tt);
		qaws_volume_mesh_frame_field(fr, sing, tt, NULL);
		for (i = 0; i < nc; i++)
			for (k = 0; k < 6; k++)
				ns += (sing[i] >> k) & 1;
		sprintf(sub, "%u tetrahedra, %u singular edge sides", nc, ns);
		ms_vpanel(c, 1466, 0, 733, 700, tt, 0.5, 0.75, "its frame field", sub);
		/* the surface as lines only (nothing hides the inside), the singular edges, the frames in the middle plane */
		qaws_volume_mesh_get_boundary(&tb, tt);
		if (tb)
			ms_draw(c, tb, NULL, NULL, 0, grey, 0.6);
		for (i = 0; i < nc; i++)
		{
			double cen[3] = { 0, 0, 0 };
			for (k = 0; k < 6; k++)
				if ((sing[i] >> k) & 1)
					er_line(c, x + 3 * cv[4 * i + E[k][0]], x + 3 * cv[4 * i + E[k][1]], red, 2.5);
			for (k = 0; k < 4; k++)
			{
				cen[0] += x[3 * cv[4 * i + k]] / 4;
				cen[1] += x[3 * cv[4 * i + k] + 1] / 4;
				cen[2] += x[3 * cv[4 * i + k] + 2] / 4;
			}
			if (fabs(cen[2]) < 0.06 && (i % 3) == 0)
				for (k = 0; k < 3; k++)
				{
					double a[3], b[3];
					unsigned int j;
					for (j = 0; j < 3; j++)
					{
						a[j] = cen[j] - 0.07 * (double)fr[9 * i + 3 * k + j];
						b[j] = cen[j] + 0.07 * (double)fr[9 * i + 3 * k + j];
					}
					er_line(c, a, b, axis[k], 1.2);
				}
		}
		printf("  torus frame field: %u tetrahedra, %u singular edge sides\n", nc, ns);
		qaws_mesh_3d_destroy(tb);
		free(ps); free(fr); free(x); free(cv); free(sing);
	}
	er_save(c, "showcase/mesh7_hexes.png");
	printf("mesh7_hexes.png\n");
	qaws_volume_mesh_destroy(h1);
	qaws_volume_mesh_destroy(h2);
	qaws_volume_mesh_destroy(tt);
	qaws_mesh_3d_destroy(ts);
	qaws_mesh_3d_destroy(bs);
	qaws_shell_atlas_destroy(atlas);
	qaws_surface_destroy(torus);
	{
		unsigned int i;
		for (i = 0; i < 6; i++)
			qaws_surface_destroy(box[i]);
	}
	er_destroy(c);
}
