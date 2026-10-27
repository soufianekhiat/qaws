/*
 * Test 186: attribute promotion (qaws_mesh_promote.h) on a unit cube of 12 triangles, each against a brute force
 *
 *   - a vertex height to the faces by AVERAGE: each face the mean of its corners' heights;
 *   - a face number to the vertices by MAXIMUM with an index attribute: each vertex the greatest number of the faces about
 *     it, its index that face;
 *   - a corner value to the vertices by AVERAGE: each vertex the mean of its corners';
 *   - the height to the edges by SUM: each edge its two ends' sum;
 *   - pieces: the face numbers averaged over each side (a piece attribute on the faces, the same element): each side's two
 *     faces 2 s + 0.5;
 *   - the whole mesh (reduce): values 1 2 2 3 3 3 7 7 give MODE 3 (from vertex 3), MEDIAN 3 (the higher middle), RMS
 *     sqrt(134 / 8), MAXIMUM 7 (vertex 6);
 *   - a NORMAL averaged to the faces stays unit;
 *   - refused: the same element without a piece, no such attribute, a new name another attribute holds.
 */

#include "test_common.h"
#include "qaws_mesh_promote.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static unsigned int const pr_tri[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };

static qaws_mesh_3d* pr_cube(double* x)
{
	qaws_scalar xs[24];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (i == 1 || i == 2 || i == 5 || i == 6) ? 1 : 0;
		x[3 * i + 1] = (i == 2 || i == 3 || i == 6 || i == 7) ? 1 : 0;
		x[3 * i + 2] = i >= 4 ? 1 : 0;
	}
	for (i = 0; i < 24; i++)
		xs[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = xs; d.vertex_count = 8; d.face_vertices = pr_tri; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static void pr_add(qaws_mesh_3d* m, char const* name, qaws_mesh_element el, unsigned int n, qaws_attribute_rule rule, qaws_scalar const* v)
{
	qaws_attribute_desc d;
	unsigned int id;
	d.name = name; d.element = el; d.components = n; d.rule = rule;
	qaws_mesh_3d_add_attribute(&id, m, v, &d);
}

static qaws_scalar* pr_get(qaws_mesh_3d const* m, char const* name, unsigned int count)
{
	unsigned int id;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (count + 1));
	if (qaws_mesh_3d_find_attribute(&id, m, name) != QAWS_STATUS_OK)
	{
		free(v);
		return NULL;
	}
	qaws_mesh_3d_get_attribute(v, m, id);
	return v;
}

static qaws_status pr_promote(qaws_mesh_3d* m, char const* name, qaws_mesh_element el, qaws_promote_method method, char const* new_name,
	char const* piece, char const* index_name)
{
	qaws_attribute_promote_desc d;
	memset(&d, 0, sizeof d);
	d.name = name; d.element = el; d.method = method; d.new_name = new_name; d.piece = piece; d.index_name = index_name;
	return qaws_mesh_3d_promote_attribute(m, &d);
}

static void test_relations(void)
{
	double x[24];
	qaws_mesh_3d* m = pr_cube(x);
	qaws_scalar h[8], fid[12], cv[36], side[12], *fa, *vm, *vw, *va, *es, *ps;
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), *ev = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne), i, k, j;
	double w1 = 0, w2 = 0, w3 = 0, w4 = 0, w5 = 0;
	qaws_status s[6];
	char msg[300];
	for (i = 0; i < 8; i++)
		h[i] = (qaws_scalar)x[3 * i + 2];
	for (i = 0; i < 12; i++)
	{
		fid[i] = (qaws_scalar)i;
		side[i] = (qaws_scalar)(i / 2);
	}
	for (i = 0; i < 36; i++)
		cv[i] = (qaws_scalar)(i * 0.5);
	pr_add(m, "h", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, h);
	pr_add(m, "id", QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_LINEAR, fid);
	pr_add(m, "c", QAWS_MESH_CORNER, 1, QAWS_ATTRIBUTE_LINEAR, cv);
	pr_add(m, "side", QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_NEAREST, side);
	s[0] = pr_promote(m, "h", QAWS_MESH_FACE, QAWS_PROMOTE_AVERAGE, "h_face", NULL, NULL);
	s[1] = pr_promote(m, "id", QAWS_MESH_VERTEX, QAWS_PROMOTE_MAXIMUM, "id_max", NULL, "id_who");
	s[2] = pr_promote(m, "c", QAWS_MESH_VERTEX, QAWS_PROMOTE_AVERAGE, "c_vertex", NULL, NULL);
	s[3] = pr_promote(m, "h", QAWS_MESH_EDGE, QAWS_PROMOTE_SUM, "h_edge", NULL, NULL);
	s[4] = pr_promote(m, "id", QAWS_MESH_FACE, QAWS_PROMOTE_AVERAGE, "id_side", "side", NULL);
	fa = pr_get(m, "h_face", 12); vm = pr_get(m, "id_max", 8); vw = pr_get(m, "id_who", 8); va = pr_get(m, "c_vertex", 8); es = pr_get(m, "h_edge", ne);
	ps = pr_get(m, "id_side", 12);
	qaws_mesh_3d_get_edges(ev, m);
	for (i = 0; fa && i < 12; i++)
		w1 = fmax(w1, fabs((double)fa[i] - (x[3 * pr_tri[3 * i] + 2] + x[3 * pr_tri[3 * i + 1] + 2] + x[3 * pr_tri[3 * i + 2] + 2]) / 3));
	for (i = 0; vm && vw && va && i < 8; i++)
	{
		double best = -1, sum = 0;
		unsigned int cnt = 0;
		for (k = 0; k < 12; k++)
			for (j = 0; j < 3; j++)
				if (pr_tri[3 * k + j] == i)
				{
					best = fmax(best, k);
					sum += 0.5 * (3 * k + j);
					cnt++;
				}
		w2 = fmax(w2, fabs((double)vm[i] - best) + fabs((double)vw[i] - best));
		w3 = fmax(w3, fabs((double)va[i] - sum / cnt));
	}
	for (i = 0; es && i < ne; i++)
		w4 = fmax(w4, fabs((double)es[i] - (x[3 * ev[2 * i] + 2] + x[3 * ev[2 * i + 1] + 2])));
	for (i = 0; ps && i < 12; i++)
		w5 = fmax(w5, fabs((double)ps[i] - (2 * (i / 2) + 0.5)));
	sprintf(msg, "relations: vertex to face %.1e, face to vertex max and index %.1e, corner to vertex %.1e, vertex to edge sum %.1e, pieces %.1e",
		w1, w2, w3, w4, w5);
	printf("  %s\n", msg);
	TEST_ASSERT(s[0] == QAWS_STATUS_OK && s[1] == QAWS_STATUS_OK && s[2] == QAWS_STATUS_OK && s[3] == QAWS_STATUS_OK && s[4] == QAWS_STATUS_OK &&
		fa && vm && vw && va && es && ps && w1 < 1e-12 && w2 == 0 && w3 < 1e-12 && w4 == 0 && w5 < 1e-12, msg);
	free(fa); free(vm); free(vw); free(va); free(es); free(ps); free(ev);
	qaws_mesh_3d_destroy(m);
}

static void test_reduce_and_rules(void)
{
	double x[24];
	qaws_mesh_3d* m = pr_cube(x);
	qaws_scalar v[8] = { 1, 2, 2, 3, 3, 3, 7, 7 }, nrm[24], mode, med, rms, mx, *fn, h[8];
	int imode = -2, imx = -2;
	unsigned int i, id;
	double wu = 0;
	qaws_status r1, r2, r3, r4;
	char msg[300];
	pr_add(m, "v", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, v);
	qaws_mesh_3d_reduce_attribute(&mode, &imode, m, "v", QAWS_PROMOTE_MODE);
	qaws_mesh_3d_reduce_attribute(&med, NULL, m, "v", QAWS_PROMOTE_MEDIAN);
	qaws_mesh_3d_reduce_attribute(&rms, NULL, m, "v", QAWS_PROMOTE_RMS);
	qaws_mesh_3d_reduce_attribute(&mx, &imx, m, "v", QAWS_PROMOTE_MAXIMUM);
	/* (vertex normals: out of the cube's centre) */
	for (i = 0; i < 8; i++)
	{
		double l = sqrt(3 * 0.25);
		nrm[3 * i] = (qaws_scalar)((x[3 * i] - 0.5) / l);
		nrm[3 * i + 1] = (qaws_scalar)((x[3 * i + 1] - 0.5) / l);
		nrm[3 * i + 2] = (qaws_scalar)((x[3 * i + 2] - 0.5) / l);
		h[i] = 0;
	}
	pr_add(m, "n", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_NORMAL, nrm);
	pr_promote(m, "n", QAWS_MESH_FACE, QAWS_PROMOTE_AVERAGE, "n_face", NULL, NULL);
	fn = pr_get(m, "n_face", 36);
	for (i = 0; fn && i < 12; i++)
		wu = fmax(wu, fabs(sqrt((double)fn[3 * i] * fn[3 * i] + (double)fn[3 * i + 1] * fn[3 * i + 1] + (double)fn[3 * i + 2] * fn[3 * i + 2]) - 1));
	pr_add(m, "z", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, h);
	r1 = pr_promote(m, "v", QAWS_MESH_VERTEX, QAWS_PROMOTE_SUM, "vv", NULL, NULL);
	r2 = pr_promote(m, "nothing", QAWS_MESH_FACE, QAWS_PROMOTE_SUM, NULL, NULL, NULL);
	r3 = pr_promote(m, "v", QAWS_MESH_FACE, QAWS_PROMOTE_SUM, "z", NULL, NULL);
	r4 = qaws_mesh_3d_find_attribute(&id, m, "v");
	sprintf(msg, "reduce: mode %.0f (vertex %d), median %.0f, rms %.6f (%.6f), max %.0f (vertex %d); a NORMAL to faces unit within %.1e; refused %d %d %d, the original kept %d",
		(double)mode, imode, (double)med, (double)rms, sqrt(134.0 / 8), (double)mx, imx, wu, (int)r1, (int)r2, (int)r3, r4 == QAWS_STATUS_OK);
	printf("  %s\n", msg);
	TEST_ASSERT(mode == 3 && imode == 3 && med == 3 && fabs((double)rms - sqrt(134.0 / 8)) < 1e-12 && mx == 7 && imx == 6 && fn && wu < 1e-12 &&
		r1 == QAWS_STATUS_INVALID_ARGUMENT && r2 == QAWS_STATUS_INVALID_ARGUMENT && r3 == QAWS_STATUS_INVALID_ARGUMENT && r4 == QAWS_STATUS_OK, msg);
	free(fn);
	qaws_mesh_3d_destroy(m);
}

int test_186_mesh_promote_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 186: attribute promotion\n");
	test_relations();
	test_reduce_and_rules();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
