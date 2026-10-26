/*
 * Figure 66: volume attributes through hexahedra, boundary layers and the
 * volume connector, cut open (vertex colours, sRGB).
 *   - a coloured ball as hexahedra (qaws_volume_mesh_hexahedralize, GRID):
 *     each vertex the colour at its nearest surface point;
 *   - boundary layers grown inside it (qaws_volume_mesh_grow_layers): each
 *     stack of prisms its surface vertex's colour, layer after layer;
 *   - two blocks of hexahedra, one blue and one orange, joined
 *     (qaws_volume_mesh_connect): the bridge's layers go from one colour to
 *     the other (OKLAB), its cells' materials meet at the middle.
 * Printed: cells per panel.
 */

/* the cells whose centroid keep accepts, their exposed facets in vertex colours (any cell type) */
static void ms_vb_cut(er_canvas* c, qaws_volume_mesh const* v, int (*keep)(double const*), double const* vcol)
{
	static unsigned char const tet[4][4] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
	static unsigned char const hex[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	static unsigned char const pri[5][4] = { { 0, 2, 1 }, { 3, 4, 5 }, { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 2, 0, 3, 5 } };
	static unsigned char const pyr[5][4] = { { 0, 3, 2, 1 }, { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } };
	static unsigned char const pri_n[5] = { 3, 3, 4, 4, 4 }, pyr_n[5] = { 4, 3, 3, 3, 3 };
	static double const wire[3] = { 60, 60, 70 };
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), nk = qaws_volume_mesh_get_corner_count(v);
	unsigned int nf = qaws_volume_mesh_get_facet_count(v), i, k, j, ck = 0, cf = 0;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * nk);
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	unsigned int* first = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
	qaws_cell_type* ty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * nc);
	unsigned char* kept = (unsigned char*)malloc(nc);
	unsigned int* cellf = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	qaws_volume_mesh_get_positions(ps, v);
	for (i = 0; i < 3 * nv; i++) x[i] = (double)ps[i];
	qaws_volume_mesh_get_cells(cv, ty, v);
	qaws_volume_mesh_get_cell_adjacency(adj, v);
	for (i = 0; i < nc; i++)
	{
		unsigned int n = ty[i] == QAWS_CELL_TETRAHEDRON ? 4 : ty[i] == QAWS_CELL_HEXAHEDRON ? 8 : ty[i] == QAWS_CELL_PRISM ? 6 : 5;
		unsigned int m = ty[i] == QAWS_CELL_TETRAHEDRON ? 4 : ty[i] == QAWS_CELL_HEXAHEDRON ? 6 : 5;
		double cen[3] = { 0, 0, 0 };
		first[i] = ck;
		for (k = 0; k < n; k++)
			for (j = 0; j < 3; j++) cen[j] += x[3 * cv[ck + k] + j] / n;
		kept[i] = (unsigned char)keep(cen);
		for (k = 0; k < m; k++) cellf[cf + k] = i;
		ck += n;
		cf += m;
	}
	for (i = 0, cf = 0; i < nc; i++)
	{
		unsigned int m = ty[i] == QAWS_CELL_TETRAHEDRON ? 4 : ty[i] == QAWS_CELL_HEXAHEDRON ? 6 : 5;
		for (k = 0; k < m; k++, cf++)
		{
			unsigned int o = adj[cf], n, q;
			unsigned char const* l = ty[i] == QAWS_CELL_TETRAHEDRON ? tet[k] : ty[i] == QAWS_CELL_HEXAHEDRON ? hex[k] : ty[i] == QAWS_CELL_PRISM ? pri[k] : pyr[k];
			n = ty[i] == QAWS_CELL_TETRAHEDRON ? 3 : ty[i] == QAWS_CELL_HEXAHEDRON ? 4 : ty[i] == QAWS_CELL_PRISM ? pri_n[k] : pyr_n[k];
			if (!kept[i] || (o != QAWS_MESH_NONE && kept[cellf[o]]))
				continue;
			for (q = 1; q + 1 < n; q++)
			{
				unsigned int a = cv[first[i] + l[0]], b = cv[first[i] + l[q]], d = cv[first[i] + l[q + 1]];
				er_triangle(c, x + 3 * a, x + 3 * b, x + 3 * d, NULL, NULL, NULL, vcol + 3 * a, vcol + 3 * b, vcol + 3 * d);
			}
			for (q = 0; q < n; q++)
				er_line(c, x + 3 * cv[first[i] + l[q]], x + 3 * cv[first[i] + l[(q + 1) % n]], wire, 0.6);
		}
	}
	free(ps); free(x); free(cv); free(adj); free(first); free(ty); free(kept); free(cellf);
}

/* a volume's "color" per vertex as 0..255 (NULL: none) */
static double* ms_vb_colours(qaws_volume_mesh const* v)
{
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), i;
	qaws_scalar* col = ms_va_get(v, "color", 3 * nv);
	double* o;
	if (!col)
		return NULL;
	o = (double*)malloc(sizeof(double) * 3 * (nv + 1));
	for (i = 0; i < 3 * nv; i++) o[i] = 255 * fmin(1, fmax(0, (double)col[i]));
	free(col);
	return o;
}

static int ms_vb_low(double const* p)
{
	return p[1] < 1.0;
}

static void demo_volume_attributes_more(void)
{
	er_canvas* c = er_create(1800, 640, 3, 252, 252, 250);
	char sub[200];
	if (!c)
		return;
	{
		/* a coloured ball as hexahedra */
		qaws_mesh_3d* ball = ms_at_ball(1, 0, 0, 0);
		qaws_volume_mesh* v = NULL;
		qaws_hexahedralize_desc hd;
		memset(&hd, 0, sizeof hd);
		hd.cell_size = (qaws_scalar)0.2;
		if (ball && qaws_volume_mesh_hexahedralize(&v, NULL, ball, &hd) == QAWS_STATUS_OK && v)
		{
			double* vc = ms_vb_colours(v);
			sprintf(sub, "%u hexahedra, the nearest surface colour", qaws_volume_mesh_get_cell_count(v));
			ms_va_panel(c, 0, v, -0.5, 0.5, "a ball's colours in hexahedra", sub);
			if (vc)
				ms_vb_cut(c, v, ms_va_front, vc);
			printf("  hexahedra: %s\n", sub);
			free(vc);
		}
		qaws_mesh_3d_destroy(ball);
		qaws_volume_mesh_destroy(v);
	}
	{
		/* boundary layers grown inside it */
		qaws_mesh_3d* ball = ms_at_ball(1, 0, 0, 0);
		qaws_volume_mesh* v = NULL;
		qaws_boundary_layers_desc bd;
		memset(&bd, 0, sizeof bd);
		bd.layer_count = 5;
		bd.first_height = (qaws_scalar)0.03;
		bd.growth = (qaws_scalar)1.4;
		if (ball && qaws_volume_mesh_grow_layers(&v, NULL, ball, &bd) == QAWS_STATUS_OK && v)
		{
			double* vc = ms_vb_colours(v);
			sprintf(sub, "%u prisms, each stack one colour", qaws_volume_mesh_get_cell_count(v));
			ms_va_panel(c, 600, v, -0.5, 0.5, "boundary layers, coloured", sub);
			if (vc)
				ms_vb_cut(c, v, ms_va_front, vc);
			printf("  layers: %s\n", sub);
			free(vc);
		}
		qaws_mesh_3d_destroy(ball);
		qaws_volume_mesh_destroy(v);
	}
	{
		/* two blocks of 4 x 4 hexahedra a unit and a half apart, blue and orange, joined */
		qaws_volume_mesh *a = NULL, *b = NULL, *m = NULL;
		unsigned int s, n = 4;
		for (s = 0; s < 2; s++)
		{
			unsigned int w = n + 1, i, j, k, q = 0;
			qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * w * w * 2);
			unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 8 * n * n);
			qaws_cell_type* ty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * n * n);
			double* col = (double*)malloc(sizeof(double) * 3 * w * w * 2);
			qaws_volume_mesh_desc d;
			qaws_volume_mesh* blk = NULL;
			for (k = 0; k < 2; k++)
				for (j = 0; j < w; j++)
					for (i = 0; i < w; i++)
					{
						unsigned int v = (k * w + j) * w + i;
						x[3 * v] = (qaws_scalar)(0.5 * i);
						x[3 * v + 1] = (qaws_scalar)(0.5 * j);
						x[3 * v + 2] = (qaws_scalar)(s ? 2.5 + 0.5 * k : 0.5 * k);
						col[3 * v] = s ? 0.98 : 0.15;
						col[3 * v + 1] = s ? 0.55 : 0.45;
						col[3 * v + 2] = s ? 0.15 : 0.95;
					}
			for (j = 0; j < n; j++)
				for (i = 0; i < n; i++, q++)
				{
					unsigned int b0 = j * w + i, t0 = w * w + b0;
					unsigned int const h8[8] = { b0, b0 + 1, b0 + w + 1, b0 + w, t0, t0 + 1, t0 + w + 1, t0 + w };
					memcpy(cv + 8 * q, h8, sizeof h8);
					ty[q] = QAWS_CELL_HEXAHEDRON;
				}
			memset(&d, 0, sizeof d);
			d.positions = x; d.vertex_count = 2 * w * w; d.cell_vertices = cv; d.cell_types = ty; d.cell_count = n * n;
			qaws_volume_mesh_create(&blk, &d);
			if (blk)
			{
				qaws_attribute_desc ad;
				qaws_scalar* cs = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * w * w * 2);
				for (i = 0; i < 3 * w * w * 2; i++) cs[i] = (qaws_scalar)col[i];
				ad.name = "color"; ad.element = QAWS_MESH_VERTEX; ad.components = 3; ad.rule = QAWS_ATTRIBUTE_OKLAB;
				qaws_volume_mesh_add_attribute(NULL, blk, cs, &ad);
				free(cs);
			}
			if (s)
				b = blk;
			else
				a = blk;
			free(x); free(cv); free(ty); free(col);
		}
		if (a && b)
		{
			unsigned int top[16], bot[16], i;
			qaws_connect_desc cd;
			for (i = 0; i < n * n; i++)
			{
				top[i] = 6 * i + 1;
				bot[i] = 6 * i;
			}
			memset(&cd, 0, sizeof cd);
			cd.rows = 6;
			if (qaws_volume_mesh_connect(&m, NULL, a, top, b, bot, n * n, n * n, &cd) == QAWS_STATUS_OK && m)
			{
				double* vc = ms_vb_colours(m);
				sprintf(sub, "%u cells, 6 bridge layers blended in Oklab", qaws_volume_mesh_get_cell_count(m));
				ms_va_panel(c, 1200, m, -0.6, 0.35, "two blocks joined", sub);
				if (vc)
					ms_vb_cut(c, m, ms_vb_low, vc);
				printf("  joined: %s\n", sub);
				free(vc);
			}
		}
		qaws_volume_mesh_destroy(a); qaws_volume_mesh_destroy(b); qaws_volume_mesh_destroy(m);
	}
	er_save(c, "showcase/mesh66_volume_attributes_more.png");
	printf("mesh66_volume_attributes_more.png\n");
	er_destroy(c);
}
