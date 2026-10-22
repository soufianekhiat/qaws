/*
 * Figure 40: subdivision of meshes (qaws_mesh_3d_subdivide).
 *   - the cage: a ring of eight blocks with a post and a stub, quads;
 *   - Catmull-Clark, three levels, smooth (genus one kept);
 *   - the cage's edges semi-sharp (sharpness 1.5): rounded, not melted;
 *   - Loop on the cage cut into triangles, the base's rim an infinitely
 *     sharp crease (the bottom stays flat), a colour per vertex carried.
 */

/* the voxels' outer faces as quads on the lattice x + 4 y + 16 z, unused lattice points dropped */
static qaws_mesh_3d* ms_sd_cage(int triangles)
{
	static int const off[6][4][3] = { { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 1 }, { 0, 1, 0 } }, { { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 }, { 1, 0, 1 } },
		{ { 0, 0, 0 }, { 1, 0, 0 }, { 1, 0, 1 }, { 0, 0, 1 } }, { { 0, 1, 0 }, { 0, 1, 1 }, { 1, 1, 1 }, { 1, 1, 0 } },
		{ { 0, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 1, 0, 0 } }, { { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } } };
	static int const dir[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
	unsigned char occ[3][3][3];
	unsigned int f[6 * 27 * 4], fs[6 * 27 * 2], map[64], nf = 0, nq = 0, nv = 0, i, j, k, d, q;
	double x[64 * 3];
	memset(occ, 0, sizeof occ);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			occ[i][j][0] = (unsigned char)!(i == 1 && j == 1);
	occ[0][0][1] = occ[0][0][2] = occ[2][2][1] = 1;
	for (i = 0; i < 64; i++) map[i] = 0xFFFFFFFFu;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				for (d = 0; occ[i][j][k] && d < 6; d++)
				{
					int a = (int)i + dir[d][0], b = (int)j + dir[d][1], c = (int)k + dir[d][2];
					unsigned int v[4];
					if (a >= 0 && a < 3 && b >= 0 && b < 3 && c >= 0 && c < 3 && occ[a][b][c])
						continue;
					for (q = 0; q < 4; q++)
					{
						unsigned int l = (i + off[d][q][0]) + 4 * (j + off[d][q][1]) + 16 * (k + off[d][q][2]);
						if (map[l] == 0xFFFFFFFFu)
						{
							x[3 * nv] = (double)(l % 4) - 1.5;
							x[3 * nv + 1] = (double)(l / 4 % 4) - 1.5;
							x[3 * nv + 2] = (double)(l / 16);
							map[l] = nv++;
						}
						v[q] = map[l];
					}
					if (triangles)
					{
						unsigned int t[6] = { v[0], v[1], v[2], v[0], v[2], v[3] };
						memcpy(f + 3 * nf, t, sizeof t);
						nf += 2;
					}
					else
					{
						memcpy(f + 4 * nq, v, sizeof v);
						fs[nq++] = 4;
					}
				}
	return triangles ? ms_make(x, nv, f, NULL, nf) : ms_make(x, nv, f, fs, nq);
}

static void demo_subdivide(void)
{
	static double const cage[3] = { 200, 205, 215 }, smooth[3] = { 120, 150, 210 }, sharp[3] = { 130, 190, 140 }, wire[3] = { 60, 66, 80 },
		light[3] = { 150, 156, 170 };
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	char sub[200];
	qaws_mesh_subdivide_desc d;
	qaws_mesh_3d *q = ms_sd_cage(0), *o = NULL;
	double const yaw = 2.6, pitch = 0.85;
	if (!c || !q)
		return;
	ms_panel(c, 0, 0, 450, 640, q, yaw, pitch, "the cage", "a ring, a post, a stub: quads");
	ms_draw(c, q, cage, NULL, 0, wire, 0.9);
	memset(&d, 0, sizeof d);
	d.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
	d.levels = 3;
	qaws_mesh_3d_subdivide(&o, q, &d);
	if (o)
	{
		sprintf(sub, "3 levels: %u quads", qaws_mesh_3d_get_face_count(o));
		ms_panel(c, 450, 0, 450, 640, o, yaw, pitch, "Catmull-Clark", sub);
		ms_draw(c, o, smooth, NULL, 1, light, 0.3);
		printf("  Catmull-Clark: %s\n", sub);
		qaws_mesh_3d_destroy(o);
		o = NULL;
	}
	{
		/* every cage edge at a right angle a crease of sharpness 1.5 */
		unsigned int ne, i;
		qaws_scalar* sh;
		qaws_mesh_3d_detect_features(q, (qaws_scalar)0.5, 0, 0);
		ne = qaws_mesh_3d_get_feature_edge_count(q);
		sh = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ne ? ne : 1));
		for (i = 0; i < ne; i++) sh[i] = (qaws_scalar)1.5;
		d.sharpness = sh;
		qaws_mesh_3d_subdivide(&o, q, &d);
		d.sharpness = NULL;
		free(sh);
		if (o)
		{
			sprintf(sub, "%u edges sharpness 1.5", ne);
			ms_panel(c, 900, 0, 450, 640, o, yaw, pitch, "semi-sharp creases", sub);
			ms_draw(c, o, sharp, NULL, 1, light, 0.3);
			printf("  semi-sharp: %s\n", sub);
			qaws_mesh_3d_destroy(o);
			o = NULL;
		}
	}
	{
		/* Loop on triangles, the base's rim a crease, a colour per vertex */
		qaws_mesh_3d* t = ms_sd_cage(1);
		unsigned int nv = qaws_mesh_3d_get_vertex_count(t), ne, i, n = 0, a;
		qaws_scalar* col = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		unsigned int* fe;
		qaws_attribute_desc ad;
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(t, i);
			double rgb[3];
			ms_colormap(0.5 + 0.5 * sin(1.3 * (double)p.x + 0.9 * (double)p.y) * cos(0.8 * (double)p.z), rgb);
			col[3 * i] = (qaws_scalar)rgb[0];
			col[3 * i + 1] = (qaws_scalar)rgb[1];
			col[3 * i + 2] = (qaws_scalar)rgb[2];
		}
		memset(&ad, 0, sizeof ad);
		ad.name = "colour";
		ad.element = QAWS_MESH_VERTEX;
		ad.components = 3;
		ad.rule = QAWS_ATTRIBUTE_LINEAR;
		qaws_mesh_3d_add_attribute(&a, t, col, &ad);
		free(col);
		qaws_mesh_3d_detect_features(t, (qaws_scalar)0.5, 0, 0);
		ne = qaws_mesh_3d_get_feature_edge_count(t);
		fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne ? ne : 1));
		qaws_mesh_3d_get_feature_edges(fe, t);
		for (i = 0; i < ne; i++)
			if (qaws_mesh_3d_get_vertex(t, fe[2 * i]).z == 0 && qaws_mesh_3d_get_vertex(t, fe[2 * i + 1]).z == 0)
			{
				fe[2 * n] = fe[2 * i];
				fe[2 * n + 1] = fe[2 * i + 1];
				n++;
			}
		qaws_mesh_3d_set_feature_edges(t, fe, n);
		free(fe);
		d.scheme = QAWS_MESH_SUBDIVIDE_LOOP;
		d.levels = 3;
		qaws_mesh_3d_subdivide(&o, t, &d);
		if (o)
		{
			unsigned int m = qaws_mesh_3d_get_vertex_count(o), id;
			qaws_scalar* cs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * m);
			double* rgb = (double*)malloc(sizeof(double) * 3 * m);
			qaws_mesh_3d_find_attribute(&id, o, "colour");
			qaws_mesh_3d_get_attribute(cs, o, id);
			for (i = 0; i < 3 * m; i++) rgb[i] = (double)cs[i];
			sprintf(sub, "%u triangles, the rim sharp", qaws_mesh_3d_get_face_count(o));
			ms_panel(c, 1350, 0, 450, 640, o, yaw, pitch, "Loop, a colour carried", sub);
			ms_draw(c, o, NULL, rgb, 1, NULL, 0);
			printf("  Loop: %s\n", sub);
			free(cs);
			free(rgb);
			qaws_mesh_3d_destroy(o);
		}
		qaws_mesh_3d_destroy(t);
	}
	qaws_mesh_3d_destroy(q);
	er_save(c, "showcase/mesh40_subdivide.png");
	printf("mesh40_subdivide.png\n");
	er_destroy(c);
}
