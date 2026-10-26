/*
 * Test 182: attributes on volume meshes (qaws_volume_mesh_add_attribute and
 * the rest), through cloning, the border, improvement and tetrahedralization
 *
 *   - two tetrahedra sharing a facet: a value per vertex, corner, cell and
 *     facet read back, set, found, removed; a clone keeps them; the border's
 *     faces carry each facet's own value, its cell's, its corners', and the
 *     vertices' (exact);
 *   - an edge or a face element, a NORMAL of 2 components refused;
 *   - a cube's tetrahedra, two materials (cells by their centroid's side of
 *     x = 0.5), corners carrying their cell's material: improvement (a high
 *     target, flips and smoothing made) keeps each material's volume
 *     exactly, every corner its cell's material, the vertex values as given;
 *   - tetrahedralization of a cube carrying x + 2 y + 3 z per vertex and a
 *     side label per face: every border vertex the field exactly, every
 *     border facet its side's label, the inner facets 0.
 */

#include "test_common.h"
#include "qaws_volume_mesh.h"
#include "qaws_volume_improve.h"
#include "qaws_tetrahedralize.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static unsigned int va_add(qaws_volume_mesh* m, char const* name, qaws_mesh_element e, unsigned int n, qaws_attribute_rule rule, qaws_scalar const* v)
{
	qaws_attribute_desc d;
	unsigned int id = QAWS_MESH_NONE;
	d.name = name; d.element = e; d.components = n; d.rule = rule;
	qaws_volume_mesh_add_attribute(&id, m, v, &d);
	return id;
}

static qaws_scalar* va_get(qaws_volume_mesh const* m, char const* name, unsigned int count)
{
	unsigned int id;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (count + 1));
	if (qaws_volume_mesh_find_attribute(&id, m, name) != QAWS_STATUS_OK || qaws_volume_mesh_get_attribute(v, m, id) != QAWS_STATUS_OK)
	{
		free(v);
		return NULL;
	}
	return v;
}

static qaws_scalar* va_surface_get(qaws_mesh_3d const* m, char const* name, unsigned int count)
{
	unsigned int id;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (count + 1));
	if (qaws_mesh_3d_find_attribute(&id, m, name) != QAWS_STATUS_OK || qaws_mesh_3d_get_attribute(v, m, id) != QAWS_STATUS_OK)
	{
		free(v);
		return NULL;
	}
	return v;
}

static void test_api(void)
{
	/* two tetrahedra on the facet 0 1 2, apexes 3 above and 4 below */
	static unsigned int const cv[8] = { 0, 1, 2, 3, 0, 2, 1, 4 };
	qaws_scalar x[15] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.3f, 0.3f, 1, 0.3f, 0.3f, -1 };
	qaws_scalar vv[5], cc[8], ce[2] = { 7, 9 }, fa[8], back[8];
	qaws_volume_mesh_desc d;
	qaws_volume_mesh *m = NULL, *k = NULL;
	qaws_mesh_3d* b = NULL;
	qaws_attribute_desc ad;
	unsigned int i, id, nb, nc, rf, rm, bad = 0, n_found;
	qaws_status s1, s2, s3;
	char msg[400];
	for (i = 0; i < 5; i++) vv[i] = (qaws_scalar)(10 + i);
	for (i = 0; i < 8; i++)
	{
		cc[i] = (qaws_scalar)(100 + i);
		fa[i] = (qaws_scalar)(1000 + i);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 5; d.cell_vertices = cv; d.cell_count = 2;
	qaws_volume_mesh_create(&m, &d);
	va_add(m, "v", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, vv);
	va_add(m, "c", QAWS_MESH_CORNER, 1, QAWS_ATTRIBUTE_LINEAR, cc);
	id = va_add(m, "cell", QAWS_MESH_CELL, 1, QAWS_ATTRIBUTE_NEAREST, NULL);
	qaws_volume_mesh_set_attribute(m, id, ce);
	va_add(m, "f", QAWS_MESH_FACET, 1, QAWS_ATTRIBUTE_NEAREST, fa);
	va_add(m, "gone", QAWS_MESH_CELL, 1, QAWS_ATTRIBUTE_NEAREST, NULL);
	qaws_volume_mesh_find_attribute(&id, m, "gone");
	qaws_volume_mesh_remove_attribute(m, id);
	n_found = qaws_volume_mesh_get_attribute_count(m);
	qaws_volume_mesh_find_attribute(&id, m, "cell");
	qaws_volume_mesh_get_attribute(back, m, id);
	bad += back[0] != 7 || back[1] != 9;
	/* the clone, then its border: 6 faces */
	qaws_volume_mesh_clone(&k, m);
	if (k && qaws_volume_mesh_get_boundary(&b, k) == QAWS_STATUS_OK && b)
	{
		qaws_scalar *bv = va_surface_get(b, "v", qaws_mesh_3d_get_vertex_count(b)), *bvert = va_surface_get(b, "vertex", qaws_mesh_3d_get_vertex_count(b));
		qaws_scalar *bc = va_surface_get(b, "c", qaws_mesh_3d_get_corner_count(b)), *bcell = va_surface_get(b, "cell", qaws_mesh_3d_get_face_count(b));
		qaws_scalar *bf = va_surface_get(b, "f", qaws_mesh_3d_get_face_count(b));
		unsigned int f[18];
		nb = qaws_mesh_3d_get_face_count(b);
		qaws_mesh_3d_get_triangles(f, b);
		for (i = 0; bv && bvert && i < qaws_mesh_3d_get_vertex_count(b); i++)
			bad += bv[i] != 10 + bvert[i];
		/* faces 0 1 2 the upper cell's (facets 1, 2, 3 of cell 0 are the border), 3 4 5 the lower one's */
		for (i = 0; bcell && bf && bc && i < nb; i++)
		{
			unsigned int cell = (unsigned int)bcell[i], j, cvv;
			bad += (bcell[i] == 0 ? 7 : 9) != (bcell[i] == 0 ? ce[0] : ce[1]);
			bad += bf[i] < 1000 + 4 * cell || bf[i] > 1003 + 4 * cell;
			/* each corner's value is that of its vertex's corner in the cell */
			for (j = 0; j < 3; j++)
			{
				unsigned int vtx = (unsigned int)bvert[f[3 * i + j]], q;
				cvv = 0;
				for (q = 0; q < 4; q++)
					if (cv[4 * cell + q] == vtx)
						cvv = 100 + 4 * cell + q;
				bad += bc[3 * i + j] != (qaws_scalar)cvv;
			}
		}
		free(bv); free(bvert); free(bc); free(bcell); free(bf);
	}
	else
		bad++;
	nc = qaws_volume_mesh_get_corner_count(m);
	ad.name = "e"; ad.element = QAWS_MESH_EDGE; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	s1 = qaws_volume_mesh_add_attribute(NULL, m, NULL, &ad);
	ad.element = QAWS_MESH_FACE;
	s2 = qaws_volume_mesh_add_attribute(NULL, m, NULL, &ad);
	ad.element = QAWS_MESH_CELL; ad.components = 2; ad.rule = QAWS_ATTRIBUTE_NORMAL;
	s3 = qaws_volume_mesh_add_attribute(NULL, m, NULL, &ad);
	rf = qaws_volume_mesh_get_attribute_count(k);
	rm = (unsigned int)(s1 == QAWS_STATUS_INVALID_ARGUMENT && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT);
	sprintf(msg, "two tetrahedra: %u attributes after a removal (4), the clone's %u, %u corners, border values off in %u places; edge, face, a NORMAL of 2 refused: %s",
		n_found, rf, nc, bad, rm ? "yes" : "no");
	printf("  %s\n", msg);
	TEST_ASSERT(n_found == 4 && rf == 4 && nc == 8 && bad == 0 && rm, msg);
	qaws_volume_mesh_destroy(m); qaws_volume_mesh_destroy(k); qaws_mesh_3d_destroy(b);
}

/* a box's surface (triangles, outward) */
static qaws_mesh_3d* va_box(double h)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_scalar x[24];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? h : 0);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? h : 0);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? h : 0);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

/* tetrahedron volume */
static double va_tet(qaws_scalar const* x, unsigned int const* v)
{
	double a[3], b[3], c[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		a[k] = (double)x[3 * v[1] + k] - (double)x[3 * v[0] + k];
		b[k] = (double)x[3 * v[2] + k] - (double)x[3 * v[0] + k];
		c[k] = (double)x[3 * v[3] + k] - (double)x[3 * v[0] + k];
	}
	return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
}

/* the volume per material (1, 2) and how many corners differ from their cell's material */
static void va_materials(qaws_volume_mesh const* m, double* vol, unsigned int* off)
{
	unsigned int nv = qaws_volume_mesh_get_vertex_count(m), nc = qaws_volume_mesh_get_cell_count(m), c;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
	qaws_scalar *mat = va_get(m, "material", nc), *cor = va_get(m, "corner", 4 * nc);
	vol[0] = vol[1] = 0;
	*off = 0;
	qaws_volume_mesh_get_positions(x, m);
	qaws_volume_mesh_get_cells(cv, NULL, m);
	for (c = 0; mat && cor && c < nc; c++)
	{
		unsigned int k;
		vol[(int)mat[c] - 1] += va_tet(x, cv + 4 * c);
		for (k = 0; k < 4; k++) *off += cor[4 * c + k] != mat[c];
	}
	if (!mat || !cor)
		*off = 1000000;
	free(x); free(cv); free(mat); free(cor);
}

static void test_improve(void)
{
	qaws_mesh_3d* box = va_box(1);
	qaws_volume_mesh *t = NULL, *o = NULL;
	qaws_tetrahedralize_desc td;
	qaws_volume_improve_desc id;
	qaws_volume_improve_report r;
	double before[2], after[2];
	unsigned int off0 = 0, off1 = 0, vbad = 0;
	char msg[400];
	memset(&td, 0, sizeof td);
	td.target_length = (qaws_scalar)0.3;
	memset(&r, 0, sizeof r);
	if (qaws_volume_mesh_tetrahedralize(&t, NULL, box, &td) == QAWS_STATUS_OK && t)
	{
		unsigned int nv = qaws_volume_mesh_get_vertex_count(t), nc = qaws_volume_mesh_get_cell_count(t), c, k;
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		unsigned int* cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
		qaws_scalar *mat = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nc), *cor = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 4 * nc), *vt = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
		qaws_volume_mesh_get_positions(x, t);
		qaws_volume_mesh_get_cells(cv, NULL, t);
		for (c = 0; c < nc; c++)
		{
			double cx = 0;
			for (k = 0; k < 4; k++) cx += (double)x[3 * cv[4 * c + k]] / 4;
			mat[c] = (qaws_scalar)(cx < 0.5 ? 1 : 2);
			for (k = 0; k < 4; k++) cor[4 * c + k] = mat[c];
		}
		for (k = 0; k < nv; k++) vt[k] = (qaws_scalar)k;
		va_add(t, "material", QAWS_MESH_CELL, 1, QAWS_ATTRIBUTE_NEAREST, mat);
		va_add(t, "corner", QAWS_MESH_CORNER, 1, QAWS_ATTRIBUTE_LINEAR, cor);
		va_add(t, "id", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, vt);
		va_materials(t, before, &off0);
		memset(&id, 0, sizeof id);
		id.target_dihedral = 60;
		if (qaws_volume_mesh_improve(&o, &r, t, &id) == QAWS_STATUS_OK && o)
		{
			qaws_scalar* vo = va_get(o, "id", qaws_volume_mesh_get_vertex_count(o));
			va_materials(o, after, &off1);
			for (k = 0; vo && k < qaws_volume_mesh_get_vertex_count(o); k++) vbad += vo[k] != (qaws_scalar)k;
			if (!vo)
				vbad = 1;
			free(vo);
		}
		free(x); free(cv); free(mat); free(cor); free(vt);
	}
	sprintf(msg, "improved (%u edge removals, %u face removals, %u moves): material volumes %.15f %.15f -> %.15f %.15f, corners off their material %u -> %u, vertex values off %u",
		r.edge_removals, r.face_removals, r.moves, before[0], before[1], after[0], after[1], off0, off1, vbad);
	printf("  %s\n", msg);
	TEST_ASSERT(o && r.edge_removals + r.face_removals > 0 && r.moves > 0 && fabs(before[0] - after[0]) < 1e-12 && fabs(before[1] - after[1]) < 1e-12 &&
		fabs(before[0] + before[1] - 1) < 1e-12 && off0 == 0 && off1 == 0 && vbad == 0, msg);
	qaws_mesh_3d_destroy(box); qaws_volume_mesh_destroy(t); qaws_volume_mesh_destroy(o);
}

static void test_tetrahedralize(void)
{
	qaws_mesh_3d* box = va_box(1);
	qaws_volume_mesh* t = NULL;
	qaws_tetrahedralize_desc td;
	qaws_attribute_desc ad;
	unsigned int nborder = 0, fbad = 0, inner_bad = 0, nfb = 0, i;
	double ferr = 0;
	char msg[400];
	{
		/* the field per vertex, each face's side (1 + 2 axis + (on the far side)) */
		qaws_scalar fv[8], side[12];
		unsigned int f[36];
		qaws_mesh_3d_get_triangles(f, box);
		for (i = 0; i < 8; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(box, i);
			fv[i] = p.x + 2 * p.y + 3 * p.z;
		}
		for (i = 0; i < 12; i++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(box, f[3 * i]), b = qaws_mesh_3d_get_vertex(box, f[3 * i + 1]), c = qaws_mesh_3d_get_vertex(box, f[3 * i + 2]);
			unsigned int ax = a.x == b.x && b.x == c.x ? 0 : a.y == b.y && b.y == c.y ? 1 : 2;
			double coord = ax == 0 ? (double)a.x : ax == 1 ? (double)a.y : (double)a.z;
			side[i] = (qaws_scalar)(1 + 2 * ax + (coord > 0.5));
		}
		ad.name = "field"; ad.element = QAWS_MESH_VERTEX; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_LINEAR;
		qaws_mesh_3d_add_attribute(NULL, box, fv, &ad);
		ad.name = "side"; ad.element = QAWS_MESH_FACE; ad.rule = QAWS_ATTRIBUTE_NEAREST;
		qaws_mesh_3d_add_attribute(NULL, box, side, &ad);
	}
	memset(&td, 0, sizeof td);
	td.target_length = (qaws_scalar)0.35;
	if (qaws_volume_mesh_tetrahedralize(&t, NULL, box, &td) == QAWS_STATUS_OK && t)
	{
		unsigned int nv = qaws_volume_mesh_get_vertex_count(t), nc = qaws_volume_mesh_get_cell_count(t), nf = qaws_volume_mesh_get_facet_count(t), c, k;
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		unsigned int *cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc), *adj = (unsigned int*)malloc(sizeof(unsigned int) * nf);
		qaws_scalar *field = va_get(t, "field", nv), *side = va_get(t, "side", nf);
		static unsigned char const tf[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
		qaws_volume_mesh_get_positions(x, t);
		qaws_volume_mesh_get_cells(cv, NULL, t);
		qaws_volume_mesh_get_cell_adjacency(adj, t);
		for (i = 0; field && i < nv; i++)
		{
			double px = (double)x[3 * i], py = (double)x[3 * i + 1], pz = (double)x[3 * i + 2];
			int on = fabs(px) < 1e-12 || fabs(px - 1) < 1e-12 || fabs(py) < 1e-12 || fabs(py - 1) < 1e-12 || fabs(pz) < 1e-12 || fabs(pz - 1) < 1e-12;
			if (!on)
				continue;
			nborder++;
			ferr = fmax(ferr, fabs((double)field[i] - (px + 2 * py + 3 * pz)));
		}
		for (c = 0; side && c < nc; c++)
			for (k = 0; k < 4; k++)
			{
				unsigned int f = 4 * c + k, j, ax = 3;
				double cen[3] = { 0, 0, 0 }, want;
				if (adj[f] != QAWS_MESH_NONE)
				{
					inner_bad += side[f] != 0;
					continue;
				}
				nfb++;
				for (j = 0; j < 3; j++)
				{
					unsigned int v = cv[4 * c + tf[k][j]];
					cen[0] += (double)x[3 * v] / 3; cen[1] += (double)x[3 * v + 1] / 3; cen[2] += (double)x[3 * v + 2] / 3;
				}
				for (j = 0; j < 3; j++)
					if (fabs(cen[j]) < 1e-9 || fabs(cen[j] - 1) < 1e-9)
						ax = j;
				want = ax < 3 ? 1 + 2 * ax + (cen[ax] > 0.5) : -1;
				fbad += (double)side[f] != want;
			}
		if (!field || !side)
			fbad = 1000000;
		free(x); free(cv); free(adj); free(field); free(side);
	}
	sprintf(msg, "tetrahedralized: %u border vertices, the field off by %.1e; %u border facets, %u off their side's label, %u inner facets nonzero",
		nborder, ferr, nfb, fbad, inner_bad);
	printf("  %s\n", msg);
	TEST_ASSERT(t && nborder >= 8 && ferr < 1e-12 && nfb > 0 && fbad == 0 && inner_bad == 0, msg);
	qaws_mesh_3d_destroy(box); qaws_volume_mesh_destroy(t);
}

int test_182_volume_attributes_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 182: attributes on volume meshes\n");
	test_api();
	test_improve();
	test_tetrahedralize();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
