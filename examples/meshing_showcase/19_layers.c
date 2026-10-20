/*
 * Figure 19: boundary layers (qaws_volume_mesh_grow_layers), cut open.
 *   - a blob remeshed on its PN reference, eight layers of prisms grown inward;
 *   - a box meshed in quads, layers of hexahedra, a corner taken out;
 *   - the blob's layers and tetrahedra filling what they leave inside (the
 *     last layer's offset surface tetrahedralized): a mesh for a flow solver.
 * The layers' cut faces go from dark (by the surface) to light (the last
 * layer), the tetrahedra's are pale green, the surface plain.
 */

/* the faces of a cell type, facet by facet (VTK's, as qaws_volume_mesh.h), 3 or 4 vertices each (a -1 ends a triangle) */
static int const ms_facets[4][6][4] = {
	{ { 1, 2, 3, -1 }, { 0, 3, 2, -1 }, { 0, 1, 3, -1 }, { 0, 2, 1, -1 } },
	{ { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } },
	{ { 0, 2, 1, -1 }, { 3, 4, 5, -1 }, { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 2, 0, 3, 5 } },
	{ { 0, 3, 2, 1 }, { 0, 1, 4, -1 }, { 1, 2, 4, -1 }, { 2, 3, 4, -1 }, { 3, 0, 4, -1 } } };
static unsigned int const ms_ncorner[4] = { 4, 8, 6, 5 }, ms_nfacet[4] = { 4, 6, 5, 5 };

/* the kept cells' faces no other kept cell shares; layer cells (layer_cells of them a layer, layers of them first) coloured
   by their layer, the rest by `fill`, the border's faces `surface` */
static void ms_layer_cutaway(er_canvas* c, qaws_volume_mesh const* v, int (*keep)(double const*), unsigned int layer_cells, unsigned int layers,
	double const* surface, double const* fill, double const* wire)
{
	unsigned int nv = qaws_volume_mesh_get_vertex_count(v), nc = qaws_volume_mesh_get_cell_count(v), i, k, j, *off, *f0;
	qaws_scalar* ps = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	double* x = (double*)malloc(sizeof(double) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * qaws_volume_mesh_get_corner_count(v));
	unsigned int* adj = (unsigned int*)malloc(sizeof(unsigned int) * qaws_volume_mesh_get_facet_count(v));
	qaws_cell_type* ty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * nc);
	unsigned char* kept = (unsigned char*)malloc(nc);
	off = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
	f0 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
	qaws_volume_mesh_get_positions(ps, v);
	for (i = 0; i < 3 * nv; i++)
		x[i] = (double)ps[i];
	qaws_volume_mesh_get_cells(cv, ty, v);
	qaws_volume_mesh_get_cell_adjacency(adj, v);
	off[0] = f0[0] = 0;
	for (i = 0; i < nc; i++)
	{
		double cen[3] = { 0, 0, 0 };
		off[i + 1] = off[i] + ms_ncorner[ty[i]];
		f0[i + 1] = f0[i] + ms_nfacet[ty[i]];
		for (k = 0; k < ms_ncorner[ty[i]]; k++)
			for (j = 0; j < 3; j++)
				cen[j] += x[3 * cv[off[i] + k] + j] / ms_ncorner[ty[i]];
		kept[i] = (unsigned char)keep(cen);
	}
	for (i = 0; i < nc; i++)
	{
		double col[3];
		if (!kept[i])
			continue;
		if (i < layer_cells * layers)
		{
			double t = layers > 1 ? (double)(i / layer_cells) / (layers - 1) : 0;
			col[0] = 70 + 150 * t; col[1] = 95 + 125 * t; col[2] = 150 + 90 * t;
		}
		else
			memcpy(col, fill, sizeof col);
		for (k = 0; k < ms_nfacet[ty[i]]; k++)
		{
			unsigned int o = adj[f0[i] + k], n = ms_facets[ty[i]][k][3] < 0 ? 3 : 4, q;
			double const* p[4];
			double const* rgb = o == QAWS_MESH_NONE ? surface : col;
			if (o != QAWS_MESH_NONE && kept[o])
				continue;
			for (q = 0; q < n; q++)
				p[q] = x + 3 * cv[off[i] + ms_facets[ty[i]][k][q]];
			er_triangle(c, p[0], p[1], p[2], NULL, NULL, NULL, rgb, rgb, rgb);
			if (n == 4)
				er_triangle(c, p[0], p[2], p[3], NULL, NULL, NULL, rgb, rgb, rgb);
			for (q = 0; q < n; q++)
				er_line(c, p[q], p[(q + 1) % n], wire, 0.7);
		}
	}
	free(ps); free(x); free(cv); free(adj); free(ty); free(kept); free(off); free(f0);
}

static int ms_keep_front(double const* p)
{
	return p[1] > 0.0;
}

static int ms_keep_corner(double const* p)
{
	return !(p[0] > 0 && p[1] < 0 && p[2] > 0);
}

static qaws_volume_mesh* ms_layers(qaws_mesh_3d const* s, unsigned int n, double h, double g, char* sub, char const* what)
{
	qaws_boundary_layers_desc d;
	qaws_mesh_report rep;
	qaws_volume_mesh* v = NULL;
	qaws_volume_quality q;
	qaws_status st;
	memset(&d, 0, sizeof d);
	memset(&rep, 0, sizeof rep);
	d.layer_count = n; d.first_height = (qaws_scalar)h; d.growth = (qaws_scalar)g;
	st = qaws_volume_mesh_grow_layers(&v, &rep, s, &d);
	if (st != QAWS_STATUS_OK)
	{
		printf("  %s: grow_layers failed: %s\n", what, qaws_status_to_string(st));
		return NULL;
	}
	qaws_volume_mesh_measure_quality(&q, v);
	sprintf(sub, "%u layers, %u %s, %u thinned", n, q.prisms ? q.prisms : q.hexahedra, q.prisms ? "prisms" : "hexahedra", rep.moves);
	printf("  %s: %u prisms, %u hexahedra, volume %.6f, %u stacks thinned, %u halving rounds, %s\n", what, q.prisms, q.hexahedra, (double)q.volume, rep.moves,
		rep.iterations, qaws_status_to_string(qaws_volume_mesh_validate(v)));
	return v;
}

/* the last layer's copy of the surface (same faces, the last layer's positions), tetrahedralized, and merged with the layers */
static qaws_volume_mesh* ms_hybrid(qaws_mesh_3d const* s, qaws_volume_mesh const* layers, unsigned int n, char* sub)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(s), nf = qaws_mesh_3d_get_face_count(s), nlc = qaws_volume_mesh_get_cell_count(layers);
	unsigned int nlv = qaws_volume_mesh_get_vertex_count(layers), i, *fv, *fs, nt = 0, ntv = 0;
	qaws_scalar* lx = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nlv);
	qaws_mesh_3d_desc md;
	qaws_mesh_3d* inner = NULL;
	qaws_volume_mesh *tets = NULL, *all = NULL;
	qaws_tetrahedralize_desc td;
	fv = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(s));
	fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	qaws_volume_mesh_get_positions(lx, layers);
	qaws_mesh_3d_get_faces(fv, fs, s);
	memset(&md, 0, sizeof md);
	md.positions = lx + 3 * (size_t)n * nv; md.vertex_count = nv; md.face_vertices = fv; md.face_sizes = fs; md.face_count = nf;
	qaws_mesh_3d_create(&inner, &md);
	memset(&td, 0, sizeof td);
	if (inner && qaws_volume_mesh_tetrahedralize(&tets, NULL, inner, &td) == QAWS_STATUS_OK)
	{
		/* the tetrahedra's border vertices are the last layer's (the same positions): matched, the others appended */
		unsigned int tv = qaws_volume_mesh_get_vertex_count(tets), tc = qaws_volume_mesh_get_cell_count(tets), j, k;
		qaws_scalar* tx = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * tv);
		unsigned int *map = (unsigned int*)malloc(sizeof(unsigned int) * tv), *tcv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * tc);
		qaws_scalar* ax = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nlv + tv));
		unsigned int* acv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_volume_mesh_get_corner_count(layers) + 4 * tc));
		qaws_cell_type* aty = (qaws_cell_type*)malloc(sizeof(qaws_cell_type) * (nlc + tc));
		qaws_volume_mesh_desc vd;
		qaws_volume_mesh_get_positions(tx, tets);
		qaws_volume_mesh_get_cells(tcv, NULL, tets);
		memcpy(ax, lx, sizeof(qaws_scalar) * 3 * nlv);
		qaws_volume_mesh_get_cells(acv, aty, layers);
		ntv = nlv;
		for (j = 0; j < tv; j++)
		{
			map[j] = 0xFFFFFFFFu;
			for (k = 0; k < nv && map[j] == 0xFFFFFFFFu; k++)
				if (fabs((double)tx[3 * j] - lx[3 * ((size_t)n * nv + k)]) + fabs((double)tx[3 * j + 1] - lx[3 * ((size_t)n * nv + k) + 1])
					+ fabs((double)tx[3 * j + 2] - lx[3 * ((size_t)n * nv + k) + 2]) < 1e-9)   /* (the lattice may round them) */
					map[j] = n * nv + k;
			if (map[j] == 0xFFFFFFFFu)
			{
				memcpy(ax + 3 * (size_t)ntv, tx + 3 * (size_t)j, sizeof(qaws_scalar) * 3);
				map[j] = ntv++;
			}
		}
		nt = qaws_volume_mesh_get_corner_count(layers);
		for (j = 0; j < 4 * tc; j++)
			acv[nt + j] = map[tcv[j]];
		for (j = 0; j < tc; j++)
			aty[nlc + j] = QAWS_CELL_TETRAHEDRON;
		memset(&vd, 0, sizeof vd);
		vd.positions = ax; vd.vertex_count = ntv; vd.cell_vertices = acv; vd.cell_types = aty; vd.cell_count = nlc + tc;
		if (qaws_volume_mesh_create(&all, &vd) == QAWS_STATUS_OK)
		{
			qaws_volume_quality q;
			qaws_volume_mesh_measure_quality(&q, all);
			sprintf(sub, "%u prisms, %u tetrahedra, one mesh", q.prisms, q.tetrahedra);
			printf("  hybrid: %u prisms, %u tetrahedra, volume %.6f, %s\n", q.prisms, q.tetrahedra, (double)q.volume, qaws_status_to_string(qaws_volume_mesh_validate(all)));
		}
		free(tx); free(map); free(tcv); free(ax); free(acv); free(aty);
	}
	(void)i;
	qaws_volume_mesh_destroy(tets);
	qaws_mesh_3d_destroy(inner);
	free(lx); free(fv); free(fs);
	return all;
}

static void demo_layers(void)
{
	static double const wire[3] = { 35, 40, 55 }, sand[3] = { 226, 200, 150 }, leaf[3] = { 196, 226, 186 };
	er_canvas* c = er_create(2200, 700, 3, 252, 252, 250);
	qaws_mesh_3d *blob = ms_blob(1), *rblob = NULL, *bq = NULL;
	qaws_mesh_reference* rb = NULL;
	qaws_surface* box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_volume_mesh *v1 = NULL, *v2 = NULL, *v3 = NULL;
	char sub[200];
	if (!c || !blob)
		return;
	qaws_mesh_reference_create(&rb, blob, NULL);
	rblob = rb ? ms_remesh(blob, rb, 0.1, NULL, NULL) : NULL;
	v1 = rblob ? ms_layers(rblob, 8, 0.006, 1.3, sub, "blob") : NULL;
	if (v1)
	{
		ms_vpanel(c, 0, 0, 733, 700, v1, 0.5, 0.7, "a blob, prisms grown inward", sub);
		ms_layer_cutaway(c, v1, ms_keep_front, qaws_mesh_3d_get_face_count(rblob), 8, sand, leaf, wire);
	}
	ms_box_faces(box, 1, 0.7, 0.5);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	bq = atlas ? ms_generate(NULL, atlas, QAWS_ELEMENTS_QUADS, 0.12, 0, NULL) : NULL;
	v2 = bq ? ms_layers(bq, 6, 0.008, 1.35, sub, "box") : NULL;
	if (v2)
	{
		ms_vpanel(c, 733, 0, 733, 700, v2, 0.62, 0.55, "a box in quads, hexahedra layers", sub);
		ms_layer_cutaway(c, v2, ms_keep_corner, qaws_mesh_3d_get_face_count(bq), 6, sand, leaf, wire);
	}
	v3 = v1 ? ms_hybrid(rblob, v1, 8, sub) : NULL;
	if (v3)
	{
		ms_vpanel(c, 1466, 0, 733, 700, v3, 0.5, 0.7, "layers and tetrahedra inside", sub);
		ms_layer_cutaway(c, v3, ms_keep_front, qaws_mesh_3d_get_face_count(rblob), 8, sand, leaf, wire);
	}
	er_save(c, "showcase/mesh19_layers.png");
	printf("mesh19_layers.png\n");
	qaws_volume_mesh_destroy(v1);
	qaws_volume_mesh_destroy(v2);
	qaws_volume_mesh_destroy(v3);
	qaws_mesh_3d_destroy(blob);
	qaws_mesh_3d_destroy(rblob);
	qaws_mesh_3d_destroy(bq);
	qaws_mesh_reference_destroy(rb);
	qaws_shell_atlas_destroy(atlas);
	{
		unsigned int i;
		for (i = 0; i < 6; i++)
			qaws_surface_destroy(box[i]);
	}
	er_destroy(c);
}
