/*
 * Test 47: B-rep (Boundary Representation) shell
 *
 * Tests B-rep topology creation, adjacency queries, validation,
 * volume and area computation.
 */

#include "test_common.h"
#include "qaws_brep.h"
#include "qaws_surface_bilinear.h"
#include "qaws_surface.h"

static void test_brep_box(void)
{
	/* Build a box B-rep from 6 bilinear faces.
	   Should form a valid closed manifold shell. */
	qaws_brep_shell* shell = NULL;
	qaws_status status;
	qaws_brep_vertex_id vids[8];
	qaws_brep_edge_id eids[12];
	qaws_brep_face_id fids[6];
	qaws_surface* faces[6] = {NULL};
	unsigned int val_flags = 0;
	int chi = 0;
	unsigned int i;

	printf("  test_brep_box...\n");

	/* Box vertices (unit cube offset from origin to avoid degenerate faces):
	   Cube from (0.5, 0.5, 0.5) to (1.5, 1.5, 1.5), volume = 1.0 */
	qaws_vec3 positions[8] = {
		{(qaws_scalar)0.5,(qaws_scalar)0.5,(qaws_scalar)0.5},
		{(qaws_scalar)1.5,(qaws_scalar)0.5,(qaws_scalar)0.5},
		{(qaws_scalar)1.5,(qaws_scalar)1.5,(qaws_scalar)0.5},
		{(qaws_scalar)0.5,(qaws_scalar)1.5,(qaws_scalar)0.5},
		{(qaws_scalar)0.5,(qaws_scalar)0.5,(qaws_scalar)1.5},
		{(qaws_scalar)1.5,(qaws_scalar)0.5,(qaws_scalar)1.5},
		{(qaws_scalar)1.5,(qaws_scalar)1.5,(qaws_scalar)1.5},
		{(qaws_scalar)0.5,(qaws_scalar)1.5,(qaws_scalar)1.5}
	};

	status = qaws_brep_create(&shell);
	TEST_ASSERT_STATUS(status);
	TEST_ASSERT(shell != NULL, "shell created");

	/* Add vertices */
	for (i = 0; i < 8; i++)
	{
		status = qaws_brep_add_vertex(shell, positions[i], &vids[i]);
		TEST_ASSERT_STATUS(status);
	}
	TEST_ASSERT(qaws_brep_get_vertex_count(shell) == 8, "8 vertices");

	/* Add edges (12 edges of a cube) */
	{
		qaws_brep_edge_desc edesc;
		/* Bottom face edges */
		edesc.curve_3d = NULL;
		edesc.v_start = 0; edesc.v_end = 1; qaws_brep_add_edge(shell, &edesc, &eids[0]);
		edesc.v_start = 1; edesc.v_end = 2; qaws_brep_add_edge(shell, &edesc, &eids[1]);
		edesc.v_start = 2; edesc.v_end = 3; qaws_brep_add_edge(shell, &edesc, &eids[2]);
		edesc.v_start = 3; edesc.v_end = 0; qaws_brep_add_edge(shell, &edesc, &eids[3]);
		/* Top face edges */
		edesc.v_start = 4; edesc.v_end = 5; qaws_brep_add_edge(shell, &edesc, &eids[4]);
		edesc.v_start = 5; edesc.v_end = 6; qaws_brep_add_edge(shell, &edesc, &eids[5]);
		edesc.v_start = 6; edesc.v_end = 7; qaws_brep_add_edge(shell, &edesc, &eids[6]);
		edesc.v_start = 7; edesc.v_end = 4; qaws_brep_add_edge(shell, &edesc, &eids[7]);
		/* Vertical edges */
		edesc.v_start = 0; edesc.v_end = 4; qaws_brep_add_edge(shell, &edesc, &eids[8]);
		edesc.v_start = 1; edesc.v_end = 5; qaws_brep_add_edge(shell, &edesc, &eids[9]);
		edesc.v_start = 2; edesc.v_end = 6; qaws_brep_add_edge(shell, &edesc, &eids[10]);
		edesc.v_start = 3; edesc.v_end = 7; qaws_brep_add_edge(shell, &edesc, &eids[11]);
	}
	TEST_ASSERT(qaws_brep_get_edge_count(shell) == 12, "12 edges");

	/* Create bilinear surfaces for each face */
	{
		qaws_surface_bilinear_desc bdesc;
		/* Bottom: z=0 */
		bdesc.p00 = positions[0]; bdesc.p10 = positions[1];
		bdesc.p01 = positions[3]; bdesc.p11 = positions[2];
		qaws_surface_create_bilinear(&bdesc, &faces[0]);

		/* Top: z=1 */
		bdesc.p00 = positions[4]; bdesc.p10 = positions[5];
		bdesc.p01 = positions[7]; bdesc.p11 = positions[6];
		qaws_surface_create_bilinear(&bdesc, &faces[1]);

		/* Front: y=0 */
		bdesc.p00 = positions[0]; bdesc.p10 = positions[1];
		bdesc.p01 = positions[4]; bdesc.p11 = positions[5];
		qaws_surface_create_bilinear(&bdesc, &faces[2]);

		/* Back: y=1 */
		bdesc.p00 = positions[3]; bdesc.p10 = positions[2];
		bdesc.p01 = positions[7]; bdesc.p11 = positions[6];
		qaws_surface_create_bilinear(&bdesc, &faces[3]);

		/* Left: x=0 */
		bdesc.p00 = positions[0]; bdesc.p10 = positions[3];
		bdesc.p01 = positions[4]; bdesc.p11 = positions[7];
		qaws_surface_create_bilinear(&bdesc, &faces[4]);

		/* Right: x=1 */
		bdesc.p00 = positions[1]; bdesc.p10 = positions[2];
		bdesc.p01 = positions[5]; bdesc.p11 = positions[6];
		qaws_surface_create_bilinear(&bdesc, &faces[5]);
	}

	/* Add faces with edge references */
	{
		qaws_brep_face_desc fdesc;
		qaws_brep_edge_id fedges[4];
		int foris[4];

		/* Bottom face (z=min): du×dv points +z (inward), need flip */
		fdesc.surface = faces[0]; fdesc.is_reversed = 1;
		fedges[0] = eids[0]; fedges[1] = eids[1]; fedges[2] = eids[2]; fedges[3] = eids[3];
		foris[0] = 1; foris[1] = 1; foris[2] = 1; foris[3] = 1;
		qaws_brep_add_face(shell, &fdesc, fedges, foris, 4, &fids[0]);

		/* Top face (z=max): du×dv points +z (outward), OK */
		fdesc.surface = faces[1]; fdesc.is_reversed = 0;
		fedges[0] = eids[4]; fedges[1] = eids[5]; fedges[2] = eids[6]; fedges[3] = eids[7];
		foris[0] = 1; foris[1] = 1; foris[2] = 1; foris[3] = 1;
		qaws_brep_add_face(shell, &fdesc, fedges, foris, 4, &fids[1]);

		/* Front face (y=min): du×dv points -y (outward), OK */
		fdesc.surface = faces[2]; fdesc.is_reversed = 0;
		fedges[0] = eids[0]; fedges[1] = eids[9]; fedges[2] = eids[4]; fedges[3] = eids[8];
		foris[0] = -1; foris[1] = 1; foris[2] = -1; foris[3] = -1;
		qaws_brep_add_face(shell, &fdesc, fedges, foris, 4, &fids[2]);

		/* Back face (y=max): du×dv points -y (inward), need flip */
		fdesc.surface = faces[3]; fdesc.is_reversed = 1;
		fedges[0] = eids[2]; fedges[1] = eids[11]; fedges[2] = eids[6]; fedges[3] = eids[10];
		foris[0] = -1; foris[1] = 1; foris[2] = -1; foris[3] = -1;
		qaws_brep_add_face(shell, &fdesc, fedges, foris, 4, &fids[3]);

		/* Left face (x=min): du×dv points +x (inward), need flip */
		fdesc.surface = faces[4]; fdesc.is_reversed = 1;
		fedges[0] = eids[3]; fedges[1] = eids[8]; fedges[2] = eids[7]; fedges[3] = eids[11];
		foris[0] = -1; foris[1] = 1; foris[2] = -1; foris[3] = -1;
		qaws_brep_add_face(shell, &fdesc, fedges, foris, 4, &fids[4]);

		/* Right face (x=max): du×dv points +x (outward), OK */
		fdesc.surface = faces[5]; fdesc.is_reversed = 0;
		fedges[0] = eids[1]; fedges[1] = eids[10]; fedges[2] = eids[5]; fedges[3] = eids[9];
		foris[0] = -1; foris[1] = 1; foris[2] = -1; foris[3] = -1;
		qaws_brep_add_face(shell, &fdesc, fedges, foris, 4, &fids[5]);
	}
	TEST_ASSERT(qaws_brep_get_face_count(shell) == 6, "6 faces");

	/* Euler characteristic: V - E + F = 8 - 12 + 6 = 2 (sphere topology) */
	status = qaws_brep_euler_characteristic(shell, &chi);
	TEST_ASSERT_STATUS(status);
	printf("    Euler characteristic: %d (expected 2)\n", chi);
	TEST_ASSERT(chi == 2, "Euler V-E+F = 2");

	/* Volume should be 1.0 for a unit cube */
	{
		qaws_scalar volume = 0;
		status = qaws_brep_compute_volume(shell, &volume);
		TEST_ASSERT_STATUS(status);
		printf("    Volume: %.4f (expected 1.0)\n", (double)volume);
		{
			qaws_scalar diff = volume - (qaws_scalar)1.0;
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.1, "volume ~ 1.0");
		}
	}

	/* Surface area should be 6.0 for a unit cube */
	{
		qaws_scalar area = 0;
		status = qaws_brep_compute_surface_area(shell, &area);
		TEST_ASSERT_STATUS(status);
		printf("    Surface area: %.4f (expected 6.0)\n", (double)area);
		{
			qaws_scalar diff = area - (qaws_scalar)6.0;
			if (diff < 0) diff = -diff;
			TEST_ASSERT(diff < (qaws_scalar)0.5, "area ~ 6.0");
		}
	}

	/* Validation */
	status = qaws_brep_validate(shell, &val_flags);
	TEST_ASSERT_STATUS(status);
	printf("    Validation flags: 0x%x\n", val_flags);

	/* Topology queries */
	{
		qaws_brep_vertex_id vs, ve;
		qaws_brep_face_id fl, fr;

		status = qaws_brep_get_edge_vertices(shell, eids[0], &vs, &ve);
		TEST_ASSERT_STATUS(status);
		printf("    Edge 0: v%u -> v%u\n", vs, ve);

		status = qaws_brep_get_edge_faces(shell, eids[0], &fl, &fr);
		TEST_ASSERT_STATUS(status);
		printf("    Edge 0 faces: f%u, f%u\n", fl, fr);
	}

	/* Cleanup */
	qaws_brep_destroy(shell);
	for (i = 0; i < 6; i++)
		if (faces[i]) qaws_surface_destroy(faces[i]);
}

static void test_brep_adjacency(void)
{
	/* Test adjacency queries on a simple shell */
	qaws_brep_shell* shell = NULL;
	qaws_status status;
	qaws_brep_vertex_id vid;
	qaws_brep_edge_id eid;
	qaws_brep_face_id fid;
	qaws_vec3 pos;

	printf("  test_brep_adjacency...\n");

	status = qaws_brep_create(&shell);
	TEST_ASSERT_STATUS(status);

	/* Add 3 vertices forming a triangle */
	pos.x = 0; pos.y = 0; pos.z = 0;
	qaws_brep_add_vertex(shell, pos, &vid);
	pos.x = 1;
	qaws_brep_add_vertex(shell, pos, &vid);
	pos.x = 0; pos.y = 1;
	qaws_brep_add_vertex(shell, pos, &vid);

	/* Add 3 edges */
	{
		qaws_brep_edge_desc edesc;
		edesc.curve_3d = NULL;
		edesc.v_start = 0; edesc.v_end = 1;
		qaws_brep_add_edge(shell, &edesc, &eid);
		edesc.v_start = 1; edesc.v_end = 2;
		qaws_brep_add_edge(shell, &edesc, &eid);
		edesc.v_start = 2; edesc.v_end = 0;
		qaws_brep_add_edge(shell, &edesc, &eid);
	}

	/* Add 1 face */
	{
		qaws_brep_face_desc fdesc;
		qaws_brep_edge_id fedges[3] = {0, 1, 2};
		int foris[3] = {1, 1, 1};
		fdesc.surface = NULL;
		fdesc.is_reversed = 0;
		qaws_brep_add_face(shell, &fdesc, fedges, foris, 3, &fid);
	}

	/* Query vertex edges */
	{
		qaws_brep_edge_id vedges[8];
		unsigned int count = 0;
		status = qaws_brep_get_vertex_edges(shell, 0, vedges, 8, &count);
		TEST_ASSERT_STATUS(status);
		printf("    Vertex 0 has %u edges\n", count);
		TEST_ASSERT(count == 2, "vertex 0 has 2 edges");
	}

	/* Query face edges */
	{
		qaws_brep_edge_id fedges[8];
		unsigned int count = 0;
		status = qaws_brep_get_face_edges(shell, 0, fedges, 8, &count);
		TEST_ASSERT_STATUS(status);
		printf("    Face 0 has %u edges\n", count);
		TEST_ASSERT(count == 3, "face 0 has 3 edges");
	}

	/* Query vertex position */
	{
		qaws_vec3 vpos;
		status = qaws_brep_get_vertex_position(shell, 1, &vpos);
		TEST_ASSERT_STATUS(status);
		TEST_ASSERT(approx_eq(vpos.x, (qaws_scalar)1.0), "vertex 1 x=1");
	}

	qaws_brep_destroy(shell);
}

int test_47_brep_main(void)
{
	g_pass = 0;
	g_fail = 0;

	printf("Test 47: B-rep shell\n");
	test_brep_box();
	test_brep_adjacency();

	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
