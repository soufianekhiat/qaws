/*
 * Test 46: Subdivision surfaces (Catmull-Clark & Loop)
 *
 * Tests subdivision surface creation from control meshes,
 * evaluation, and OBJ export.
 */

#include "test_common.h"
#include "qaws_surface_subdiv.h"
#include "qaws_surface.h"

/* Write a subdivision mesh directly to OBJ (bypasses parametric eval which
   doesn't handle closed surfaces well). */
static void obj_subdiv_mesh(obj_writer* obj, qaws_surface const* surface)
{
	qaws_vec3 const* verts;
	qaws_vec3 const* normals;
	unsigned int const* tris;
	unsigned int vert_count, tri_count;
	unsigned int first_v, first_n;
	unsigned int i;

	if (qaws_surface_subdiv_get_mesh(surface, &verts, &normals, &tris,
		&vert_count, &tri_count) != QAWS_STATUS_OK)
		return;

	first_v = obj->vertex_count + 1;
	first_n = obj->normal_count + 1;

	for (i = 0; i < vert_count; i++)
	{
		obj_vertex(obj, verts[i]);
		obj_normal(obj, normals[i]);
	}

	for (i = 0; i < tri_count; i++)
	{
		unsigned int i0 = first_v + tris[i * 3 + 0];
		unsigned int i1 = first_v + tris[i * 3 + 1];
		unsigned int i2 = first_v + tris[i * 3 + 2];
		unsigned int n0 = first_n + tris[i * 3 + 0];
		unsigned int n1 = first_n + tris[i * 3 + 1];
		unsigned int n2 = first_n + tris[i * 3 + 2];
		fprintf(obj->fp, "f %u//%u %u//%u %u//%u\n",
			i0, n0, i1, n1, i2, n2);
	}
}

static void test_catmull_clark_cube(void)
{
	/* Create a Catmull-Clark subdivision surface from a cube.
	   8 vertices, 6 quad faces. After subdivision, this should
	   produce a rounded shape approaching a sphere. */
	qaws_surface_subdiv_desc desc;
	qaws_surface* surface = NULL;
	qaws_status status;
	qaws_surface_eval_result r;

	printf("  test_catmull_clark_cube...\n");

	/* Cube vertices */
	qaws_vec3 verts[8] = {
		{-1, -1, -1}, { 1, -1, -1}, { 1,  1, -1}, {-1,  1, -1},
		{-1, -1,  1}, { 1, -1,  1}, { 1,  1,  1}, {-1,  1,  1}
	};

	/* 6 quad faces (CCW from outside) */
	unsigned int indices[] = {
		0, 3, 2, 1,  /* bottom (-z) */
		4, 5, 6, 7,  /* top (+z) */
		0, 1, 5, 4,  /* front (-y) */
		2, 3, 7, 6,  /* back (+y) */
		0, 4, 7, 3,  /* left (-x) */
		1, 2, 6, 5   /* right (+x) */
	};
	unsigned int face_sizes[] = {4, 4, 4, 4, 4, 4};

	desc.vertices = verts;
	desc.vertex_count = 8;
	desc.face_indices = indices;
	desc.face_sizes = face_sizes;
	desc.face_count = 6;
	desc.scheme = QAWS_SUBDIV_CATMULL_CLARK;
	desc.subdivision_level = 3;

	status = qaws_surface_create_subdiv(&desc, &surface);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(surface != NULL, "CC subdivision surface created");
	TEST_ASSERT(qaws_surface_get_kind(surface) == QAWS_SURFACE_KIND_SUBDIV,
		"surface kind is SUBDIV");

	/* Evaluate at center */
	memset(&r, 0, sizeof(r));
	status = qaws_surface_evaluate(surface, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
	TEST_ASSERT_STATUS(status);
	printf("    center: (%.4f, %.4f, %.4f)\n",
		(double)r.position.x, (double)r.position.y, (double)r.position.z);
	printf("    normal: (%.4f, %.4f, %.4f)\n",
		(double)r.normal.x, (double)r.normal.y, (double)r.normal.z);

	/* Check corners and edges */
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, (qaws_scalar)0, (qaws_scalar)0,
		QAWS_SURFACE_EVAL_POSITION, &r);
	printf("    corner(0,0): (%.4f, %.4f, %.4f)\n",
		(double)r.position.x, (double)r.position.y, (double)r.position.z);

	/* Verify mesh data is available */
	{
		qaws_vec3 const* mv;
		qaws_vec3 const* mn;
		unsigned int const* mt;
		unsigned int mvc, mtc;
		status = qaws_surface_subdiv_get_mesh(surface, &mv, &mn, &mt, &mvc, &mtc);
		TEST_ASSERT_STATUS(status);
		printf("    mesh: %u vertices, %u triangles\n", mvc, mtc);
		TEST_ASSERT(mvc > 0, "mesh has vertices");
		TEST_ASSERT(mtc > 0, "mesh has triangles");
	}

	/* OBJ export using direct mesh */
	{
		obj_writer obj;
		svg_ensure_output_dir();
		if (obj_open(&obj, OBJ_OUTPUT_DIR "/46_subdiv_cube.obj",
			OBJ_OUTPUT_DIR "/46_subdiv_cube.mtl"))
		{
			obj_material(&obj, "subdiv_cube", 0.9, 0.5, 0.3);
			obj_group(&obj, "cc_cube");
			obj_use_material(&obj, "subdiv_cube");
			obj_subdiv_mesh(&obj, surface);
			obj_close(&obj);
			printf("    wrote 46_subdiv_cube.obj\n");
		}
	}

	if (surface) qaws_surface_destroy(surface);
}

static void test_loop_tetrahedron(void)
{
	/* Create a Loop subdivision surface from a tetrahedron.
	   4 vertices, 4 triangle faces. After subdivision, this should
	   produce a rounded shape. */
	qaws_surface_subdiv_desc desc;
	qaws_surface* surface = NULL;
	qaws_status status;
	qaws_surface_eval_result r;

	printf("  test_loop_tetrahedron...\n");

	qaws_vec3 verts[4] = {
		{ 1,  1,  1},
		{ 1, -1, -1},
		{-1,  1, -1},
		{-1, -1,  1}
	};

	unsigned int indices[] = {
		0, 1, 2,
		0, 2, 3,
		0, 3, 1,
		1, 3, 2
	};
	unsigned int face_sizes[] = {3, 3, 3, 3};

	desc.vertices = verts;
	desc.vertex_count = 4;
	desc.face_indices = indices;
	desc.face_sizes = face_sizes;
	desc.face_count = 4;
	desc.scheme = QAWS_SUBDIV_LOOP;
	desc.subdivision_level = 3;

	status = qaws_surface_create_subdiv(&desc, &surface);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(surface != NULL, "Loop subdivision surface created");

	/* Evaluate */
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_NORMAL, &r);
	printf("    center: (%.4f, %.4f, %.4f)\n",
		(double)r.position.x, (double)r.position.y, (double)r.position.z);

	/* OBJ export using direct mesh */
	{
		obj_writer obj;
		svg_ensure_output_dir();
		if (obj_open(&obj, OBJ_OUTPUT_DIR "/46_subdiv_tetra.obj",
			OBJ_OUTPUT_DIR "/46_subdiv_tetra.mtl"))
		{
			obj_material(&obj, "subdiv_tetra", 0.3, 0.9, 0.5);
			obj_group(&obj, "loop_tetra");
			obj_use_material(&obj, "subdiv_tetra");
			obj_subdiv_mesh(&obj, surface);
			obj_close(&obj);
			printf("    wrote 46_subdiv_tetra.obj\n");
		}
	}

	if (surface) qaws_surface_destroy(surface);
}

static void test_subdiv_plane(void)
{
	/* Subdivide a single quad (plane) - should remain flat. */
	qaws_surface_subdiv_desc desc;
	qaws_surface* surface = NULL;
	qaws_status status;
	qaws_surface_eval_result r;

	printf("  test_subdiv_plane...\n");

	qaws_vec3 verts[4] = {
		{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}
	};
	unsigned int indices[] = {0, 1, 2, 3};
	unsigned int face_sizes[] = {4};

	desc.vertices = verts;
	desc.vertex_count = 4;
	desc.face_indices = indices;
	desc.face_sizes = face_sizes;
	desc.face_count = 1;
	desc.scheme = QAWS_SUBDIV_CATMULL_CLARK;
	desc.subdivision_level = 2;

	status = qaws_surface_create_subdiv(&desc, &surface);
	TEST_ASSERT_STATUS(status);

	/* Evaluate at center - z should be ~0 */
	memset(&r, 0, sizeof(r));
	qaws_surface_evaluate(surface, (qaws_scalar)0.5, (qaws_scalar)0.5,
		QAWS_SURFACE_EVAL_POSITION, &r);
	printf("    center z: %.6f (should be ~0)\n", (double)r.position.z);
	{
		qaws_scalar abs_z = r.position.z < 0 ? -r.position.z : r.position.z;
		TEST_ASSERT(abs_z < (qaws_scalar)0.01, "plane stays flat");
	}

	if (surface) qaws_surface_destroy(surface);
}

int test_46_subdiv_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 46: Subdivision surfaces\n");
	test_catmull_clark_cube();
	test_loop_tetrahedron();
	test_subdiv_plane();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
