/*
 * Figure 65: attributes on volume meshes, cut open.
 *   - a ball's sRGB colours carried into its tetrahedra
 *     (qaws_volume_mesh_tetrahedralize): each inside vertex the colour at
 *     its nearest point of the surface;
 *   - a box's side labels (one per face) on its tetrahedra's border facets
 *     (QAWS_MESH_FACET: boundary conditions), the inner facets 0 (grey);
 *   - two materials, an inclusion in a block (QAWS_MESH_CELL by the cells'
 *     centroids), before and after qaws_volume_mesh_improve: the facets
 *     between materials stay as walls, each material's volume exact.
 * Printed: cells, the materials' volumes, the worst dihedral angles.
 */

typedef struct ms_va_view
{
	int mode;                       /* 0 a vertex colour "color", 1 a facet label "side", 2 a cell material "material" */
	double const* vcol;             /* per vertex, 0..255 */
	qaws_scalar const* lab;         /* per facet or per cell */
} ms_va_view;

static void ms_va_label(double v, double* rgb)
{
	static double const k[7][3] = { { 200, 200, 200 }, { 230, 90, 80 }, { 240, 170, 60 }, { 90, 180, 100 }, { 70, 150, 220 }, { 150, 110, 210 }, { 230, 120, 180 } };
	unsigned int i = v < 0 ? 0 : v > 6 ? 6 : (unsigned int)v, j;
	for (j = 0; j < 3; j++) rgb[j] = k[i][j];
}

/* the tetrahedra whose centroid keep accepts, their exposed facets coloured by the view */
static void ms_va_cut(er_canvas* c, qaws_volume_mesh const* v, int (*keep)(double const*), ms_va_view const* view)
{
	static unsigned char const F[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	static double const wire[3] = { 60, 60, 70 };
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k, j;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	unsigned char* kept = (unsigned char*)malloc(nc);
	qaws_volume_mesh_get_positions(ps, v);
	for (i = 0; i < 3 * nv; i++) x[i] = (double)ps[i];
	qaws_volume_mesh_get_cells(cv, NULL, v);
	qaws_volume_mesh_get_cell_adjacency(adj, v);
	for (i = 0; i < nc; i++)
	{
		double cen[3] = { 0, 0, 0 };
		for (k = 0; k < 4; k++)
			for (j = 0; j < 3; j++) cen[j] += x[3 * cv[4 * i + k] + j] / 4;
		kept[i] = (unsigned char)keep(cen);
	}
	for (i = 0; i < nc; i++)
		for (k = 0; kept[i] && k < 4; k++)
		{
			unsigned int o = adj[4 * i + k], a = cv[4 * i + F[k][0]], b = cv[4 * i + F[k][1]], d = cv[4 * i + F[k][2]];
			double col[3];
			if (o != QAWS_MESH_NONE && kept[o])
				continue;
			if (view->mode == 0)
				er_triangle(c, x + 3 * a, x + 3 * b, x + 3 * d, NULL, NULL, NULL, view->vcol + 3 * a, view->vcol + 3 * b, view->vcol + 3 * d);
			else
			{
				ms_va_label((double)view->lab[view->mode == 1 ? 4 * i + k : i], col);
				er_triangle(c, x + 3 * a, x + 3 * b, x + 3 * d, NULL, NULL, NULL, col, col, col);
			}
			er_line(c, x + 3 * a, x + 3 * b, wire, 0.6);
			er_line(c, x + 3 * b, x + 3 * d, wire, 0.6);
			er_line(c, x + 3 * d, x + 3 * a, wire, 0.6);
		}
	free(ps); free(x); free(cv); free(adj); free(kept);
}

static int ms_va_front(double const* p)
{
	return p[1] > 0;
}

static int ms_va_half(double const* p)
{
	return p[1] > 0.5;
}

static int ms_va_corner(double const* p)
{
	return !(p[0] > 0.5 && p[1] < 0.5 && p[2] > 0.5);
}

static void ms_va_panel(er_canvas* c, int px, qaws_volume_mesh const* v, double yaw, double pitch, char const* title, char const* sub)
{
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), i;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	qaws_volume_mesh_get_positions(ps, v);
	for (i = 0; i < 3 * nv; i++) x[i] = (double)ps[i];
	er_panel(c, px, 72, 600, 640 - 72, x, nv, yaw, pitch, 0.06);
	er_text(c, px + 14, 12, title, 3, 30, 30, 30);
	er_text(c, px + 14, 46, sub, 2, 95, 95, 95);
	free(ps); free(x);
}

static qaws_scalar* ms_va_get(qaws_volume_mesh const* v, char const* name, unsigned int count)
{
	unsigned int id;
	qaws_scalar* o = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (count + 1));
	if (qaws_volume_mesh_find_attribute(&id, v, name) != QAWS_STATUS_OK)
	{
		free(o);
		return NULL;
	}
	qaws_volume_mesh_get_attribute(o, v, id);
	return o;
}

/* the materials' volumes (1, 2), the cells' worst dihedral angle (degrees) */
static void ms_va_measure(qaws_volume_mesh const* v, double* vol, double* worst)
{
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k;
	qaws_scalar *ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv), *mat = ms_va_get(v, "material", nc);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	qaws_volume_quality q;
	qaws_volume_mesh_get_positions(ps, v);
	qaws_volume_mesh_get_cells(cv, NULL, v);
	vol[0] = vol[1] = 0;
	for (i = 0; mat && i < nc; i++)
	{
		double a[3], b[3], d[3];
		qaws_scalar const* p0 = ps + 3 * cv[4 * i];
		for (k = 0; k < 3; k++)
		{
			a[k] = (double)ps[3 * cv[4 * i + 1] + k] - (double)p0[k];
			b[k] = (double)ps[3 * cv[4 * i + 2] + k] - (double)p0[k];
			d[k] = (double)ps[3 * cv[4 * i + 3] + k] - (double)p0[k];
		}
		vol[(int)mat[i] - 1] += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) + a[2] * (b[0] * d[1] - b[1] * d[0])) / 6;
	}
	qaws_volume_mesh_measure_quality(&q, v);
	*worst = (double)q.min_dihedral;
	free(ps); free(mat); free(cv);
}

/* the unit box's surface, n by n squares per side (two triangles each, outward); per face its side (1 + 2 axis + far) */
static qaws_mesh_3d* ms_va_grid_box(unsigned int n, double* side)
{
	unsigned int nmax = 6 * (n + 1) * (n + 1), nv = 0, nt = 0, s, i, j, k;
	double* x = (double*)malloc(sizeof(double) * 3 * nmax);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * 6 * n * n);
	unsigned int* id = (unsigned int*)malloc(sizeof(unsigned int) * (n + 1) * (n + 1));
	qaws_mesh_3d* m;
	for (s = 0; s < 6; s++)
	{
		unsigned int ax = s / 2, far = s % 2, u = (ax + 1) % 3, w = (ax + 2) % 3;
		for (j = 0; j <= n; j++)
			for (i = 0; i <= n; i++)
			{
				double p[3];
				p[ax] = far; p[u] = (double)i / n; p[w] = (double)j / n;
				for (k = 0; k < nv; k++)
					if (x[3 * k] == p[0] && x[3 * k + 1] == p[1] && x[3 * k + 2] == p[2])
						break;
				if (k == nv)
				{
					memcpy(x + 3 * nv, p, sizeof p);
					nv++;
				}
				id[j * (n + 1) + i] = k;
			}
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
			{
				unsigned int a = id[j * (n + 1) + i], b = id[j * (n + 1) + i + 1], c = id[(j + 1) * (n + 1) + i + 1], d = id[(j + 1) * (n + 1) + i];
				/* (u, w, ax) right-handed: counter-clockwise in (u, w) faces +ax, the far side's outward */
				if (far)
				{
					t[3 * nt] = a; t[3 * nt + 1] = b; t[3 * nt + 2] = c; side[nt++] = 1 + 2 * ax + far;
					t[3 * nt] = a; t[3 * nt + 1] = c; t[3 * nt + 2] = d; side[nt++] = 1 + 2 * ax + far;
				}
				else
				{
					t[3 * nt] = a; t[3 * nt + 1] = c; t[3 * nt + 2] = b; side[nt++] = 1 + 2 * ax + far;
					t[3 * nt] = a; t[3 * nt + 1] = d; t[3 * nt + 2] = c; side[nt++] = 1 + 2 * ax + far;
				}
			}
	}
	m = ms_make(x, nv, t, NULL, nt);
	free(x); free(t); free(id);
	return m;
}

static void demo_volume_attributes(void)
{
	er_canvas* c = er_create(2400, 640, 3, 252, 252, 250);
	qaws_tetrahedralize_desc td;
	char sub[200];
	if (!c)
		return;
	memset(&td, 0, sizeof td);
	{
		/* a coloured ball filled */
		qaws_mesh_3d* ball = ms_at_ball(1, 0, 0, 0);
		qaws_volume_mesh* v = NULL;
		td.target_length = (qaws_scalar)0.16;
		if (ball && qaws_volume_mesh_tetrahedralize(&v, NULL, ball, &td) == QAWS_STATUS_OK && v)
		{
			unsigned int nv = qaws_volume_mesh_get_vertex_count(v), i;
			qaws_scalar* col = ms_va_get(v, "color", 3 * nv);
			double* vc = (double*)malloc(sizeof(double) * 3 * (nv + 1));
			ms_va_view view;
			for (i = 0; col && i < 3 * nv; i++) vc[i] = 255 * fmin(1, fmax(0, (double)col[i]));
			view.mode = 0; view.vcol = vc; view.lab = NULL;
			sprintf(sub, "%u cells, inside the nearest colour", qaws_volume_mesh_get_cell_count(v));
			ms_va_panel(c, 0, v, -0.5, 0.5, "a ball's colours, inside", sub);
			if (col)
				ms_va_cut(c, v, ms_va_front, &view);
			printf("  ball: %s\n", sub);
			free(col); free(vc);
		}
		qaws_mesh_3d_destroy(ball);
		qaws_volume_mesh_destroy(v);
	}
	{
		/* a box's side labels on its border facets */
		double side[6 * 2 * 64];
		qaws_mesh_3d* box = ms_va_grid_box(8, side);
		qaws_volume_mesh* v = NULL;
		ms_at_add(box, "side", QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_NEAREST, side, qaws_mesh_3d_get_face_count(box));
		td.target_length = (qaws_scalar)0.12;
		if (qaws_volume_mesh_tetrahedralize(&v, NULL, box, &td) == QAWS_STATUS_OK && v)
		{
			qaws_scalar* lab = ms_va_get(v, "side", qaws_volume_mesh_get_facet_count(v));
			ms_va_view view;
			view.mode = 1; view.vcol = NULL; view.lab = lab;
			sprintf(sub, "%u cells, inner facets 0 (grey)", qaws_volume_mesh_get_cell_count(v));
			ms_va_panel(c, 600, v, -0.6, 0.55, "side labels on facets", sub);
			if (lab)
				ms_va_cut(c, v, ms_va_corner, &view);
			printf("  box: %s\n", sub);
			free(lab);
		}
		qaws_mesh_3d_destroy(box);
		qaws_volume_mesh_destroy(v);
	}
	{
		/* an inclusion: material 2 inside a ball of radius 0.33 at the block's centre, before and after improvement */
		double side[6 * 2 * 100];
		qaws_mesh_3d* box = ms_va_grid_box(10, side);
		qaws_volume_mesh *v = NULL, *o = NULL;
		td.target_length = (qaws_scalar)0.1;
		if (qaws_volume_mesh_tetrahedralize(&v, NULL, box, &td) == QAWS_STATUS_OK && v)
		{
			unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k;
			qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
			double* mat = (double*)malloc(sizeof(double) * nc);
			qaws_volume_improve_desc id;
			double vb[2], va[2], wb, wa;
			qaws_volume_mesh_get_positions(ps, v);
			qaws_volume_mesh_get_cells(cv, NULL, v);
			for (i = 0; i < nc; i++)
			{
				double cen[3] = { 0, 0, 0 };
				for (k = 0; k < 4; k++)
				{
					cen[0] += (double)ps[3 * cv[4 * i + k]] / 4 - 0.125;
					cen[1] += (double)ps[3 * cv[4 * i + k] + 1] / 4 - 0.125;
					cen[2] += (double)ps[3 * cv[4 * i + k] + 2] / 4 - 0.125;
				}
				mat[i] = cen[0] * cen[0] + cen[1] * cen[1] + cen[2] * cen[2] < 0.33 * 0.33 ? 2 : 1;
			}
			{
				qaws_attribute_desc d;
				qaws_scalar* s = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc);
				for (i = 0; i < nc; i++) s[i] = (qaws_scalar)mat[i];
				d.name = "material"; d.element = QAWS_MESH_CELL; d.components = 1; d.rule = QAWS_ATTRIBUTE_NEAREST;
				qaws_volume_mesh_add_attribute(NULL, v, s, &d);
				free(s);
			}
			memset(&id, 0, sizeof id);
			id.target_dihedral = 50;
			qaws_volume_mesh_improve(&o, NULL, v, &id);
			ms_va_measure(v, vb, &wb);
			{
				qaws_scalar* lab = ms_va_get(v, "material", nc);
				ms_va_view view;
				view.mode = 2; view.vcol = NULL; view.lab = lab;
				sprintf(sub, "%.6f + %.6f, worst %.1f deg", vb[0], vb[1], wb);
				ms_va_panel(c, 1200, v, -0.5, 0.5, "two materials, as made", sub);
				ms_va_cut(c, v, ms_va_half, &view);
				printf("  materials: %s\n", sub);
				free(lab);
			}
			if (o)
			{
				qaws_scalar* lab = ms_va_get(o, "material", qaws_volume_mesh_get_cell_count(o));
				ms_va_view view;
				ms_va_measure(o, va, &wa);
				view.mode = 2; view.vcol = NULL; view.lab = lab;
				sprintf(sub, "%.6f + %.6f, worst %.1f deg", va[0], va[1], wa);
				ms_va_panel(c, 1800, o, -0.5, 0.5, "improved, interface kept", sub);
				if (lab)
					ms_va_cut(c, o, ms_va_half, &view);
				printf("  improved: %s\n", sub);
				free(lab);
			}
			free(ps); free(cv); free(mat);
		}
		qaws_mesh_3d_destroy(box);
		qaws_volume_mesh_destroy(v);
		qaws_volume_mesh_destroy(o);
	}
	er_save(c, "showcase/mesh65_volume_attributes.png");
	printf("mesh65_volume_attributes.png\n");
	er_destroy(c);
}
