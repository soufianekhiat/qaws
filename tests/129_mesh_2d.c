/*
 * Test 129: planar meshes in the corner store (qaws_mesh_2d_get_mesh)
 *
 *   - an annulus triangulated: its view a valid mesh at z = 0 with the
 *     same triangles, one component with two boundary loops, its boundary
 *     edges features; quadrangulated: the view's faces quads;
 *   - a square from arrays, edited: an edge split, a face split, a flip, a
 *     collapse, each checked and the getters kept in step (triangles,
 *     boundary edges, the area); an attribute through the view; remeshing
 *     refused on it.
 */

#include "test_common.h"
#include <math.h>
#include <string.h>

static double t2_area(qaws_mesh_2d const* m)
{
	unsigned int i, nt = qaws_mesh_2d_get_triangle_count(m);
	double a = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int t[3];
		qaws_vec2 p, q, r;
		qaws_mesh_2d_get_triangle(m, i, t);
		p = qaws_mesh_2d_get_vertex(m, t[0]);
		q = qaws_mesh_2d_get_vertex(m, t[1]);
		r = qaws_mesh_2d_get_vertex(m, t[2]);
		a += 0.5 * (((double)q.x - p.x) * ((double)r.y - p.y) - ((double)q.y - p.y) * ((double)r.x - p.x));
	}
	return a;
}

static void test_annulus(void)
{
	qaws_curve *o = NULL, *h = NULL, *hr = NULL;
	qaws_curve const* cs[2];
	qaws_path_2d p[2];
	qaws_triangulate_desc d;
	qaws_mesh_2d* m = NULL;
	unsigned int i;
	char msg[300];
	qaws_curve_create_ellipse_2d(0, 0, 2, 2, 0, &o);
	qaws_curve_create_ellipse_2d((qaws_scalar)0.3, 0, 1, 1, 0, &h);
	qaws_curve_reverse(h, &hr);
	cs[0] = o;
	cs[1] = hr;
	for (i = 0; i < 2; i++)
	{
		p[i].curves = &cs[i];
		p[i].curve_count = 1;
		p[i].closed = 1;
	}
	memset(&d, 0, sizeof d);
	d.paths = p;
	d.path_count = 2;
	d.fill_rule = QAWS_FILL_NON_ZERO;
	d.tolerance = (qaws_scalar)2e-3;
	d.min_angle = 25;
	qaws_triangulate_region_2d(&d, &m);
	if (!m)
		TEST_ASSERT(0, "annulus: triangulated");
	else
	{
		qaws_mesh_3d const* v = qaws_mesh_2d_get_mesh(m);
		qaws_mesh_topology top;
		unsigned int nv = qaws_mesh_3d_get_vertex_count(v), nf = qaws_mesh_3d_get_face_count(v), same = 1, flat = 1;
		unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nf + 1));
		qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
		qaws_mesh_3d_get_topology(v, &top);
		qaws_mesh_3d_get_faces(v, fv, NULL);
		qaws_mesh_3d_get_positions(v, x);
		for (i = 0; i < qaws_mesh_2d_get_triangle_count(m) && i < nf; i++)
		{
			unsigned int t[3];
			qaws_mesh_2d_get_triangle(m, i, t);
			same &= t[0] == fv[3 * i] && t[1] == fv[3 * i + 1] && t[2] == fv[3 * i + 2];
		}
		for (i = 0; i < nv; i++)
			flat &= x[3 * i + 2] == 0 && (i >= qaws_mesh_2d_get_vertex_count(m)
				|| (x[3 * i] == qaws_mesh_2d_get_vertex(m, i).x && x[3 * i + 1] == qaws_mesh_2d_get_vertex(m, i).y));
		sprintf(msg, "annulus view: validate %s, %u faces for %u triangles (same %u), at z = 0 %u, %u component(s), %u boundary loops, %u feature edges for %u boundary edges",
			qaws_status_to_string(qaws_mesh_3d_validate(v)), nf, qaws_mesh_2d_get_triangle_count(m), same, flat, top.components, top.boundary_loops,
			qaws_mesh_3d_get_feature_edge_count(v), qaws_mesh_2d_get_boundary_edge_count(m));
		TEST_ASSERT(qaws_mesh_3d_validate(v) == QAWS_STATUS_OK && nf == qaws_mesh_2d_get_triangle_count(m) && same && flat && top.components == 1
			&& top.boundary_loops == 2 && qaws_mesh_3d_get_feature_edge_count(v) == qaws_mesh_2d_get_boundary_edge_count(m), msg);
		free(fv);
		free(x);
		/* quadrangulated: the view follows */
		if (qaws_mesh_2d_quadrangulate(m, 1, 2) == QAWS_STATUS_OK)
		{
			qaws_mesh_quality q;
			v = qaws_mesh_2d_get_mesh(m);
			qaws_mesh_3d_measure_quality(v, &q);
			sprintf(msg, "annulus quadrangulated: the view %u quads (the mesh %u), %u triangles", q.quads, qaws_mesh_2d_get_quad_count(m), q.triangles);
			TEST_ASSERT(q.quads == qaws_mesh_2d_get_quad_count(m) && q.quads > 0 && q.triangles == 0 && qaws_mesh_3d_validate(v) == QAWS_STATUS_OK, msg);
		}
		else
			TEST_ASSERT(0, "annulus: quadrangulated");
	}
	qaws_mesh_2d_destroy(m);
	qaws_curve_destroy(o);
	qaws_curve_destroy(h);
	qaws_curve_destroy(hr);
}

static void test_edits(void)
{
	static qaws_scalar const p[8] = { 0, 0, 1, 0, 1, 1, 0, 1 };
	static unsigned int const f[6] = { 0, 1, 2, 0, 2, 3 };
	qaws_mesh_2d_desc d;
	qaws_mesh_2d* m = NULL;
	qaws_status s1, s2, s3, s4, s5, s6;
	unsigned int nv, nt1, nt2, nt3, nb, a = QAWS_MESH_NONE;
	char msg[400];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = 4;
	d.face_vertices = f;
	d.face_count = 2;
	s1 = qaws_mesh_2d_create(&d, &m);
	if (s1 != QAWS_STATUS_OK || !m)
	{
		TEST_ASSERT(0, "square from arrays: created");
		return;
	}
	nb = qaws_mesh_2d_get_boundary_edge_count(m);
	/* the diagonal split at its midpoint, a face split at its centroid, a flip, a collapse */
	s2 = qaws_mesh_2d_split_edge(m, 0, 2, NULL, &nv);
	nt1 = qaws_mesh_2d_get_triangle_count(m);
	s3 = qaws_mesh_2d_split_face(m, 0, NULL, NULL);
	nt2 = qaws_mesh_2d_get_triangle_count(m);
	s4 = qaws_mesh_2d_flip_edge(m, 0, 1);                    /* a border edge: refused */
	s5 = qaws_mesh_3d_can_collapse_edge(qaws_mesh_2d_get_mesh(m), 4, 0) ? qaws_mesh_2d_collapse_edge(m, 5, 4, NULL) : QAWS_STATUS_OK;
	nt3 = qaws_mesh_2d_get_triangle_count(m);
	{
		qaws_attribute_desc ad;
		qaws_scalar* vals;
		unsigned int i, nvv = qaws_mesh_3d_get_vertex_count(qaws_mesh_2d_get_mesh(m));
		vals = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nvv);
		for (i = 0; i < nvv; i++)
			vals[i] = (qaws_scalar)i;
		memset(&ad, 0, sizeof ad);
		ad.name = "id";
		ad.element = QAWS_MESH_VERTEX;
		ad.components = 1;
		s6 = qaws_mesh_2d_add_attribute(m, &ad, vals, &a);
		free(vals);
	}
	sprintf(msg, "square edited: created %s, %u boundary edges; split edge %s (vertex %u, %u triangles), split face %s (%u), a border flip %s, a collapse %s (%u), area %.12f, boundary edges %u, validate %s, attribute %s, remesh %s",
		qaws_status_to_string(s1), nb, qaws_status_to_string(s2), nv, nt1, qaws_status_to_string(s3), nt2, qaws_status_to_string(s4), qaws_status_to_string(s5), nt3,
		t2_area(m), qaws_mesh_2d_get_boundary_edge_count(m), qaws_status_to_string(qaws_mesh_3d_validate(qaws_mesh_2d_get_mesh(m))), qaws_status_to_string(s6),
		qaws_status_to_string(qaws_mesh_2d_remesh(m, (qaws_scalar)0.1, NULL, NULL, 2)));
	TEST_ASSERT(nb == 4 && s2 == QAWS_STATUS_OK && nv == 4 && nt1 == 4 && s3 == QAWS_STATUS_OK && nt2 == 6 && s4 == QAWS_STATUS_INVALID_TOPOLOGY
		&& s5 == QAWS_STATUS_OK && fabs(t2_area(m) - 1) < 1e-12 && qaws_mesh_2d_get_boundary_edge_count(m) == 4
		&& qaws_mesh_3d_validate(qaws_mesh_2d_get_mesh(m)) == QAWS_STATUS_OK && s6 == QAWS_STATUS_OK && a == 0
		&& qaws_mesh_2d_remesh(m, (qaws_scalar)0.1, NULL, NULL, 2) == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	qaws_mesh_2d_destroy(m);
}

int test_129_mesh_2d_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 129: planar meshes in the corner store\n");
	test_annulus();
	test_edits();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
