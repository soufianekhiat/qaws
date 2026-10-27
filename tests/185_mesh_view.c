/*
 * Test 185: meshes from host buffers (qaws_mesh_view.h), a box of 12 triangles
 *
 *   - an interleaved f32 vertex buffer (position, then a uv) with 16-bit indices: the mesh made from it has the box's
 *     positions and triangles, its volume 1; planar f64 positions with 32-bit indices padded to 4 per triangle: the same;
 *   - refusals: an index past the vertices, a position not finite, 2 components;
 *   - a frame's positions (the box moved by (1, 2, 3)) into the mesh and into a BVH made from the view: the closest point
 *     to (5, 0.5, 0.5) at distance 4 from the box at rest, sqrt(17.5) after the refit; a refit of another count refused, the BVH as
 *     it was;
 *   - a borrowed vertex motion between two frames of the host's buffer (the box's top rising from 1 to 3) against a
 *     ceiling at 2: contact at 0.5, exact; the host's second frame rewritten (the top to 5): the same motion now meets the
 *     ceiling at 0.25, read from the buffer as it is.
 */

#include "test_common.h"
#include "qaws_mesh_view.h"
#include "qaws_mesh_sweep.h"
#include <math.h>
#include <string.h>

typedef struct mv_vertex
{
	float pos[3];
	float uv[2];
} mv_vertex;

static unsigned int const mv_tri[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };

/* the unit box's corners, z of the top four given */
static void mv_corners(double* x, double top)
{
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (i == 1 || i == 2 || i == 5 || i == 6) ? 1 : 0;
		x[3 * i + 1] = (i == 2 || i == 3 || i == 6 || i == 7) ? 1 : 0;
		x[3 * i + 2] = i >= 4 ? top : 0;
	}
}

static double mv_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = m ? qaws_mesh_3d_get_triangle_count(m) : 0, i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		v += ((double)a.x * ((double)b.y * e.z - (double)b.z * e.y) + (double)a.y * ((double)b.z * e.x - (double)b.x * e.z) +
			(double)a.z * ((double)b.x * e.y - (double)b.y * e.x)) / 6;
	}
	return v;
}

static void test_views(void)
{
	mv_vertex vb[8];
	uint16_t ib[36];
	double x[24], planes[3][8];
	uint32_t ib4[48];
	qaws_mesh_view v, w;
	qaws_mesh_3d *m = NULL, *m2 = NULL;
	unsigned int i, k, same = 1;
	qaws_status bad1, bad2, bad3;
	char msg[300];
	mv_corners(x, 1);
	for (i = 0; i < 8; i++)
	{
		for (k = 0; k < 3; k++)
		{
			vb[i].pos[k] = (float)x[3 * i + k];
			planes[k][i] = x[3 * i + k];
		}
		vb[i].uv[0] = vb[i].uv[1] = 0.5f;
	}
	for (i = 0; i < 36; i++)
		ib[i] = (uint16_t)mv_tri[i];
	for (i = 0; i < 12; i++)
	{
		for (k = 0; k < 3; k++)
			ib4[4 * i + k] = mv_tri[3 * i + k];
		ib4[4 * i + 3] = 0xDEADu;
	}
	memset(&v, 0, sizeof v);
	v.positions.data = vb;
	v.positions.type = QAWS_ELEMENT_F32;
	v.positions.components = 3;
	v.positions.count = 8;
	v.positions.stride = sizeof(mv_vertex);
	v.indices = ib;
	v.index_type = QAWS_INDEX_U16;
	v.triangle_count = 12;
	qaws_mesh_3d_create_from_view(&m, &v);
	memset(&w, 0, sizeof w);
	for (k = 0; k < 3; k++)
		w.positions.planes[k] = planes[k];
	w.positions.type = QAWS_ELEMENT_F64;
	w.positions.components = 3;
	w.positions.count = 8;
	w.indices = ib4;
	w.index_type = QAWS_INDEX_U32;
	w.index_stride = 4 * sizeof(uint32_t);
	w.triangle_count = 12;
	qaws_mesh_3d_create_from_view(&m2, &w);
	for (i = 0; m && m2 && i < 12; i++)
	{
		unsigned int a[3], b[3];
		qaws_mesh_3d_get_triangle(a, m, i);
		qaws_mesh_3d_get_triangle(b, m2, i);
		for (k = 0; k < 3; k++)
			same &= a[k] == mv_tri[3 * i + k] && b[k] == mv_tri[3 * i + k];
	}
	ib[5] = 8;
	bad1 = qaws_mesh_view_validate(&v);
	ib[5] = (uint16_t)mv_tri[5];
	planes[1][3] = NAN;
	bad2 = qaws_mesh_view_validate(&w);
	planes[1][3] = 1;
	v.positions.components = 2;
	bad3 = qaws_mesh_view_validate(&v);
	sprintf(msg, "an interleaved f32 buffer with 16-bit indices: volume %.6f, planar f64 with padded 32-bit: %.6f, triangles the same %u; refusals %d %d %d",
		mv_volume(m), mv_volume(m2), same, (int)bad1, (int)bad2, (int)bad3);
	printf("  %s\n", msg);
	TEST_ASSERT(m && m2 && fabs(mv_volume(m) - 1) < 1e-12 && fabs(mv_volume(m2) - 1) < 1e-12 && same && bad1 == QAWS_STATUS_INVALID_ARGUMENT &&
		bad2 == QAWS_STATUS_INVALID_ARGUMENT && bad3 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	qaws_mesh_3d_destroy(m); qaws_mesh_3d_destroy(m2);
}

static void test_frames(void)
{
	mv_vertex vb[8];
	uint16_t ib[36];
	double x[24];
	qaws_mesh_view v;
	qaws_mesh_3d* m = NULL;
	qaws_mesh_bvh* b = NULL;
	qaws_mesh_hit h0, h1, h2;
	qaws_scalar q[3] = { 5, (qaws_scalar)0.5, (qaws_scalar)0.5 };
	qaws_array_view short_view;
	qaws_vec3 c0;
	qaws_status s1, s2;
	unsigned int i, k;
	char msg[300];
	mv_corners(x, 1);
	for (i = 0; i < 8; i++)
		for (k = 0; k < 3; k++)
			vb[i].pos[k] = (float)x[3 * i + k];
	for (i = 0; i < 36; i++)
		ib[i] = (uint16_t)mv_tri[i];
	memset(&v, 0, sizeof v);
	v.positions.data = vb;
	v.positions.type = QAWS_ELEMENT_F32;
	v.positions.components = 3;
	v.positions.count = 8;
	v.positions.stride = sizeof(mv_vertex);
	v.indices = ib;
	v.index_type = QAWS_INDEX_U16;
	v.triangle_count = 12;
	qaws_mesh_3d_create_from_view(&m, &v);
	qaws_mesh_bvh_create_view(&b, &v, NULL);
	qaws_mesh_bvh_closest(&h0, b, q, 1, 0, NULL);
	/* the next frame: the host moves its buffer */
	for (i = 0; i < 8; i++)
	{
		vb[i].pos[0] += 1;
		vb[i].pos[1] += 2;
		vb[i].pos[2] += 3;
	}
	s1 = qaws_mesh_3d_set_positions_from_view(m, &v.positions);
	qaws_mesh_bvh_refit_view(b, &v.positions);
	qaws_mesh_bvh_closest(&h1, b, q, 1, 0, NULL);
	short_view = v.positions;
	short_view.count = 7;
	s2 = qaws_mesh_bvh_refit_view(b, &short_view);
	qaws_mesh_bvh_closest(&h2, b, q, 1, 0, NULL);
	c0 = qaws_mesh_3d_get_vertex(m, 0);
	sprintf(msg, "frames: the mesh's corner 0 at (%.1f %.1f %.1f) (1 2 3); closest from (5, 0.5, 0.5) at %.6f (4), refit %.6f (sqrt 17.5), refused %d and still %.6f",
		(double)c0.x, (double)c0.y, (double)c0.z, (double)h0.distance, (double)h1.distance, (int)s2, (double)h2.distance);
	printf("  %s\n", msg);
	TEST_ASSERT(m && b && s1 == QAWS_STATUS_OK && c0.x == 1 && c0.y == 2 && c0.z == 3 && fabs((double)h0.distance - 4) < 1e-12 &&
		fabs((double)h1.distance - sqrt(17.5)) < 1e-12 && s2 == QAWS_STATUS_INVALID_ARGUMENT && h2.distance == h1.distance, msg);
	qaws_mesh_bvh_destroy(b);
	qaws_mesh_3d_destroy(m);
}

static void test_borrowed_motion(void)
{
	mv_vertex f0[8], f1[8];
	double x[24], top;
	qaws_array_view p0, p1;
	qaws_motion* mo = NULL;
	qaws_mesh_3d *box = NULL, *ceil = NULL;
	qaws_mesh_3d_desc d;
	qaws_scalar cx[24];
	qaws_sweep_desc sd;
	qaws_contact_report r1, r2;
	unsigned int i, k;
	char msg[300];
	mv_corners(x, 1);
	for (i = 0; i < 8; i++)
		for (k = 0; k < 3; k++)
		{
			f0[i].pos[k] = (float)x[3 * i + k];
			f1[i].pos[k] = (float)x[3 * i + k] * (k == 2 ? 3.0f : 1.0f);
		}
	/* a ceiling over the box's middle, from z 2 to 2.5 */
	for (i = 0; i < 8; i++)
		for (k = 0; k < 3; k++)
			cx[3 * i + k] = (qaws_scalar)(k < 2 ? 0.25 + 0.5 * x[3 * i + k] : (x[3 * i + 2] > 0 ? 2.5 : 2));
	memset(&d, 0, sizeof d);
	d.positions = cx; d.vertex_count = 8; d.face_vertices = mv_tri; d.face_count = 12;
	qaws_mesh_3d_create(&ceil, &d);
	for (i = 0; i < 24; i++)
		cx[i] = (qaws_scalar)x[i];
	qaws_mesh_3d_create(&box, &d);
	memset(&p0, 0, sizeof p0);
	p0.data = f0;
	p0.type = QAWS_ELEMENT_F32;
	p0.components = 3;
	p0.count = 8;
	p0.stride = sizeof(mv_vertex);
	p1 = p0;
	p1.data = f1;
	qaws_motion_create_vertices_borrowed(&mo, &p0, &p1);
	memset(&sd, 0, sizeof sd);
	sd.tolerance = (qaws_scalar)1e-4;
	qaws_mesh_3d_first_contact(&r1, box, mo, ceil, NULL, &sd);
	/* the host rewrites its second frame: the top to 5 */
	top = 5;
	for (i = 4; i < 8; i++)
		f1[i].pos[2] = (float)top;
	qaws_mesh_3d_first_contact(&r2, box, mo, ceil, NULL, &sd);
	sprintf(msg, "a borrowed vertex motion: contact at %.17g (0.5, exact %d); its frame rewritten: %.17g (0.25, exact %d)", (double)r1.time, r1.exact, (double)r2.time, r2.exact);
	printf("  %s\n", msg);
	TEST_ASSERT(mo && r1.contact && r1.exact && (double)r1.time == 0.5 && r2.contact && r2.exact && (double)r2.time == 0.25, msg);
	qaws_motion_destroy(mo);
	qaws_mesh_3d_destroy(box); qaws_mesh_3d_destroy(ceil);
}

int test_185_mesh_view_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 185: meshes from host buffers\n");
	test_views();
	test_frames();
	test_borrowed_motion();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
