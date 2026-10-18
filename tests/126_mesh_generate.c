/*
 * Test 126: meshes of exact geometry (qaws_mesh_generate_surface/_shell)
 *
 *   - a bicubic B-spline patch in triangles: valid, one border loop whose
 *     vertices are on the domain's sides, every vertex the surface at its
 *     footpoint (face 0, u, v), edges near the target; in quads: all quads;
 *     adaptive by a chord tolerance: edges graded within the bounds;
 *   - a rational NURBS torus in quads: closed, genus 1, every vertex on the
 *     torus;
 *   - a box of six bilinear faces as a shell, in quads: closed, genus 0,
 *     every face on the box face its "face" attribute names, the eight
 *     corners vertices, no vertex near an edge but off it;
 *   - refusals: the reserved methods, no description, a negative length.
 */

#include "test_common.h"
#include "qaws_mesh_generate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static qaws_surface* tg_patch(void)
{
	qaws_vec3 cp[49];
	qaws_surface_bspline_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 7; i++)
		for (j = 0; j < 7; j++)
		{
			cp[i * 7 + j].x = (qaws_scalar)(i / 3.0 - 1);
			cp[i * 7 + j].y = (qaws_scalar)(j / 3.0 - 1);
			cp[i * 7 + j].z = (qaws_scalar)(0.55 * sin(1.7 * i) * cos(1.3 * j));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = cp;
	d.u_point_count = 7;
	d.v_point_count = 7;
	qaws_surface_create_bspline(&d, &s);
	return s;
}

/* the torus of radii R, r: the unit circle's 9 rational points times the profile circle's */
static qaws_surface* tg_torus(double R, double r)
{
	static double const cx[9] = { 1, 1, 0, -1, -1, -1, 0, 1, 1 }, cy[9] = { 0, 1, 1, 1, 0, -1, -1, -1, 0 };
	static qaws_scalar const kn[12] = { 0, 0, 0, 0.25f, 0.25f, 0.5f, 0.5f, 0.75f, 0.75f, 1, 1, 1 };
	qaws_vec3 cp[81];
	qaws_scalar w[81];
	qaws_surface_nurbs_desc d;
	qaws_surface* s = NULL;
	unsigned int i, j;
	for (i = 0; i < 9; i++)
		for (j = 0; j < 9; j++)
		{
			double rad = R + r * cx[j];
			cp[i * 9 + j].x = (qaws_scalar)(rad * cx[i]);
			cp[i * 9 + j].y = (qaws_scalar)(rad * cy[i]);
			cp[i * 9 + j].z = (qaws_scalar)(r * cy[j]);
			w[i * 9 + j] = (qaws_scalar)((i & 1 ? sqrt(0.5) : 1) * (j & 1 ? sqrt(0.5) : 1));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 2;
	d.v_degree = 2;
	d.control_points = cp;
	d.u_point_count = 9;
	d.v_point_count = 9;
	d.weights = w;
	d.u_knots = kn;
	d.u_knot_count = 12;
	d.v_knots = kn;
	d.v_knot_count = 12;
	qaws_surface_create_nurbs(&d, &s);
	return s;
}

/* the box [-h0, h0] x [-h1, h1] x [-h2, h2]: face 2k is at -h_k, face 2k + 1 at +h_k ... in the order bottom, top, front, back, left, right */
static int const tg_axis[6] = { 2, 2, 1, 1, 0, 0 };
static double const tg_side[6] = { -1, 1, -1, 1, -1, 1 };

static void tg_box(qaws_surface** f, double const* h)
{
	static int const q[6][4][3] = {
		{ { -1, -1, -1 }, { -1, 1, -1 }, { 1, -1, -1 }, { 1, 1, -1 } }, { { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { 1, -1, -1 }, { -1, -1, 1 }, { 1, -1, 1 } }, { { -1, 1, -1 }, { -1, 1, 1 }, { 1, 1, -1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { -1, -1, 1 }, { -1, 1, -1 }, { -1, 1, 1 } }, { { 1, -1, -1 }, { 1, 1, -1 }, { 1, -1, 1 }, { 1, 1, 1 } } };
	unsigned int i, k;
	for (i = 0; i < 6; i++)
	{
		qaws_surface_bilinear_desc d;
		qaws_vec3* p[4];
		p[0] = &d.p00; p[1] = &d.p10; p[2] = &d.p01; p[3] = &d.p11;
		for (k = 0; k < 4; k++)
		{
			p[k]->x = (qaws_scalar)(h[0] * q[i][k][0]);
			p[k]->y = (qaws_scalar)(h[1] * q[i][k][1]);
			p[k]->z = (qaws_scalar)(h[2] * q[i][k][2]);
		}
		f[i] = NULL;
		qaws_surface_create_bilinear(&d, &f[i]);
	}
}

/* a generated mesh's arrays */
typedef struct tg_mesh
{
	unsigned int nv, nf;
	double* x;
	unsigned int *fv, *fs;
	qaws_mesh_location* foot;
	qaws_mesh_topology top;
	qaws_mesh_quality q;
} tg_mesh;

static void tg_read(qaws_mesh_3d const* m, tg_mesh* t)
{
	qaws_scalar* p;
	unsigned int i;
	memset(t, 0, sizeof *t);
	t->nv = qaws_mesh_3d_get_vertex_count(m);
	t->nf = qaws_mesh_3d_get_face_count(m);
	t->x = (double*)malloc(sizeof(double) * 3 * (t->nv + 1));
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (t->nv + 1));
	t->fv = (unsigned int*)malloc(sizeof(unsigned int) * (qaws_mesh_3d_get_corner_count(m) + 1));
	t->fs = (unsigned int*)malloc(sizeof(unsigned int) * (t->nf + 1));
	t->foot = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (t->nv + 1));
	qaws_mesh_3d_get_positions(m, p);
	for (i = 0; i < 3 * t->nv; i++)
		t->x[i] = (double)p[i];
	free(p);
	qaws_mesh_3d_get_faces(m, t->fv, t->fs);
	if (qaws_mesh_3d_get_footpoints(m, t->foot) != QAWS_STATUS_OK)
		for (i = 0; i < t->nv; i++)
			t->foot[i].patch = QAWS_MESH_NONE;
	qaws_mesh_3d_get_topology(m, &t->top);
	qaws_mesh_3d_measure_quality(m, &t->q);
}

static void tg_free(tg_mesh* t)
{
	free(t->x); free(t->fv); free(t->fs); free(t->foot);
}

/* the largest distance from a vertex to its face's surface at its footpoint (HUGE_VAL when one has none) */
static double tg_foot_error(tg_mesh const* t, qaws_surface const* one, qaws_shell_atlas const* shell)
{
	double worst = 0;
	unsigned int i, k;
	for (i = 0; i < t->nv; i++)
	{
		qaws_surface const* s = one ? (t->foot[i].patch == 0 ? one : NULL) : qaws_shell_atlas_get_face(shell, t->foot[i].patch);
		qaws_surface_eval_result e;
		double d = 0, y[3];
		if (t->foot[i].patch == QAWS_MESH_NONE || !s)
			return HUGE_VAL;
		memset(&e, 0, sizeof e);
		qaws_surface_evaluate(s, t->foot[i].a, t->foot[i].b, QAWS_SURFACE_EVAL_POSITION, &e);
		y[0] = (double)e.position.x;
		y[1] = (double)e.position.y;
		y[2] = (double)e.position.z;
		for (k = 0; k < 3; k++)
			d += (y[k] - t->x[3 * i + k]) * (y[k] - t->x[3 * i + k]);
		worst = fmax(worst, sqrt(d));
	}
	return worst;
}

static void test_patch(void)
{
	qaws_surface* s = tg_patch();
	qaws_range ur = qaws_surface_get_u_range(s), vr = qaws_surface_get_v_range(s);
	qaws_generate_desc d;
	qaws_mesh_report rep;
	qaws_mesh_3d *tri = NULL, *quad = NULL, *ada = NULL;
	qaws_status s1, s2, s3;
	tg_mesh t;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.target_length = (qaws_scalar)0.15;
	s1 = qaws_mesh_generate_surface(s, &d, &tri, &rep);
	if (s1 != QAWS_STATUS_OK || !tri)
		TEST_ASSERT(0, "patch: triangles generated");
	else
	{
		unsigned int* across = (unsigned int*)malloc(sizeof(unsigned int) * qaws_mesh_3d_get_corner_count(tri));
		unsigned int f, c0, k, nborder = 0, off = 0;
		double err;
		tg_read(tri, &t);
		err = tg_foot_error(&t, s, NULL);
		/* the border's vertices on the domain's sides */
		qaws_mesh_3d_get_face_adjacency(tri, across);
		for (f = 0, c0 = 0; f < t.nf; c0 += t.fs[f], f++)
			for (k = 0; k < t.fs[f]; k++)
				if (across[c0 + k] == QAWS_MESH_NONE)
				{
					qaws_mesh_location const* l = &t.foot[t.fv[c0 + k]];
					double u = (double)l->a, v = (double)l->b, tol = 1e-9;
					nborder++;
					off += !(fabs(u - (double)ur.min_value) < tol || fabs(u - (double)ur.max_value) < tol || fabs(v - (double)vr.min_value) < tol
						|| fabs(v - (double)vr.max_value) < tol);
				}
		sprintf(msg, "patch triangles: %u faces (%u triangles), %u border loop(s), %u border edges with %u ends off the sides, footpoints to %.2e, "
			"mean edge %.3f for 0.15", t.nf, t.q.triangles, t.top.boundary_loops, nborder, off, err, (double)t.q.mean_edge);
		TEST_ASSERT(qaws_mesh_3d_validate(tri) == QAWS_STATUS_OK && t.q.triangles == t.nf && t.top.boundary_loops == 1 && t.top.components == 1
			&& nborder > 20 && off == 0 && err < 1e-9 && fabs((double)t.q.mean_edge - 0.15) < 0.4 * 0.15, msg);
		free(across);
		tg_free(&t);
	}
	d.elements = QAWS_ELEMENTS_QUADS;
	s2 = qaws_mesh_generate_surface(s, &d, &quad, NULL);
	if (s2 != QAWS_STATUS_OK || !quad)
		TEST_ASSERT(0, "patch: quads generated");
	else
	{
		double err;
		tg_read(quad, &t);
		err = tg_foot_error(&t, s, NULL);
		sprintf(msg, "patch quads: %u faces (%u quads), %u border loop(s), footpoints to %.2e", t.nf, t.q.quads, t.top.boundary_loops, err);
		TEST_ASSERT(qaws_mesh_3d_validate(quad) == QAWS_STATUS_OK && t.q.quads == t.nf && t.nf > 50 && t.top.boundary_loops == 1 && err < 1e-9, msg);
		tg_free(&t);
	}
	memset(&d, 0, sizeof d);
	d.chord_tolerance = (qaws_scalar)0.002;
	d.min_length = (qaws_scalar)0.03;
	d.max_length = (qaws_scalar)0.3;
	d.gradation = (qaws_scalar)0.4;
	s3 = qaws_mesh_generate_surface(s, &d, &ada, &rep);
	if (s3 != QAWS_STATUS_OK || !ada)
		TEST_ASSERT(0, "patch: adaptive generated");
	else
	{
		double err;
		tg_read(ada, &t);
		err = tg_foot_error(&t, s, NULL);
		sprintf(msg, "patch adaptive: %u triangles, edges %.3f..%.3f within 0.03..0.3, footpoints to %.2e", t.nf, (double)t.q.min_edge,
			(double)t.q.max_edge, err);
		TEST_ASSERT(qaws_mesh_3d_validate(ada) == QAWS_STATUS_OK && (double)t.q.min_edge > 0.5 * 0.03 && (double)t.q.max_edge < 1.5 * 0.3
			&& t.q.max_edge > 2 * t.q.min_edge && err < 1e-9, msg);
		tg_free(&t);
	}
	qaws_mesh_3d_destroy(tri);
	qaws_mesh_3d_destroy(quad);
	qaws_mesh_3d_destroy(ada);
	qaws_surface_destroy(s);
}

static void test_torus(void)
{
	qaws_surface* s = tg_torus(2, 0.7);
	qaws_generate_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_status st;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.elements = QAWS_ELEMENTS_QUADS;
	d.target_length = (qaws_scalar)0.25;
	st = qaws_mesh_generate_surface(s, &d, &m, NULL);
	if (st != QAWS_STATUS_OK || !m)
		TEST_ASSERT(0, "torus: generated");
	else
	{
		tg_mesh t;
		double off = 0, err;
		unsigned int i;
		tg_read(m, &t);
		err = tg_foot_error(&t, s, NULL);
		for (i = 0; i < t.nv; i++)
		{
			double const* p = t.x + 3 * i;
			double rho = sqrt(p[0] * p[0] + p[1] * p[1]);
			off = fmax(off, fabs(sqrt((rho - 2) * (rho - 2) + p[2] * p[2]) - 0.7));
		}
		sprintf(msg, "torus quads: %u faces (%u quads), closed %d genus %u, off the torus %.2e, footpoints to %.2e", t.nf, t.q.quads, t.top.closed,
			t.top.genus, off, err);
		TEST_ASSERT(qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.quads == t.nf && t.nf > 300 && t.top.closed && t.top.genus == 1 && off < 1e-9
			&& err < 1e-9, msg);
		tg_free(&t);
	}
	qaws_mesh_3d_destroy(m);
	qaws_surface_destroy(s);
}

static void test_box(void)
{
	static double const h[3] = { 1, 0.7, 0.5 };
	qaws_surface* box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_generate_desc d;
	qaws_mesh_report rep;
	qaws_mesh_3d* m = NULL;
	qaws_status st;
	char msg[400];
	unsigned int i;
	tg_box(box, h);
	st = qaws_shell_atlas_create((qaws_surface const* const*)box, 6, 0, &atlas);
	memset(&d, 0, sizeof d);
	d.elements = QAWS_ELEMENTS_QUADS;
	d.target_length = (qaws_scalar)0.15;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_generate_shell(atlas, &d, &m, &rep);
	if (st != QAWS_STATUS_OK || !m)
		TEST_ASSERT(0, "box: generated");
	else
	{
		tg_mesh t;
		unsigned int a, f, c0, k, offface = 0, corners = 0, near = 0, fa = QAWS_MESH_NONE;
		qaws_scalar* fid = NULL;
		double err;
		tg_read(m, &t);
		err = tg_foot_error(&t, NULL, atlas);
		/* every face on the box face its attribute names */
		if (qaws_mesh_3d_find_attribute(m, "face", &fa) == QAWS_STATUS_OK)
		{
			fid = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (t.nf + 1));
			qaws_mesh_3d_get_attribute(m, fa, fid);
			for (f = 0, c0 = 0; f < t.nf; c0 += t.fs[f], f++)
			{
				unsigned int F = (unsigned int)fid[f];
				for (k = 0; k < t.fs[f]; k++)
					offface += F >= 6 || fabs(t.x[3 * t.fv[c0 + k] + tg_axis[F]] - tg_side[F] * h[tg_axis[F]]) > 1e-9;
			}
			free(fid);
		}
		/* the corners; vertices near an edge (two coordinates near their extremes) but off it */
		for (a = 0; a < 8; a++)
			for (i = 0; i < t.nv; i++)
			{
				double const* p = t.x + 3 * i;
				if (fabs(p[0] - (a & 1 ? h[0] : -h[0])) + fabs(p[1] - (a & 2 ? h[1] : -h[1])) + fabs(p[2] - (a & 4 ? h[2] : -h[2])) < 1e-12)
				{
					corners++;
					break;
				}
			}
		for (i = 0; i < t.nv; i++)
		{
			unsigned int on = 0, close = 0;
			for (k = 0; k < 3; k++)
			{
				double gap = h[k] - fabs(t.x[3 * i + k]);
				on += gap < 1e-9;
				close += gap < 0.3 * 0.15;
			}
			near += close >= 2 && on < 2;
		}
		sprintf(msg, "box quads: %u faces (%u quads), closed %d genus %u, %u corners of 8, %u vertices off their face, %u near an edge but off it, "
			"footpoints to %.2e", t.nf, t.q.quads, t.top.closed, t.top.genus, corners, offface, near, err);
		TEST_ASSERT(qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.quads == t.nf && t.nf > 200 && t.top.closed && t.top.genus == 0
			&& fa != QAWS_MESH_NONE && offface == 0 && corners == 8 && near == 0 && err < 1e-9, msg);
		tg_free(&t);
	}
	qaws_mesh_3d_destroy(m);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
}

static void test_refusals(void)
{
	qaws_surface* s = tg_patch();
	qaws_generate_desc d;
	qaws_mesh_3d* out = NULL;
	qaws_status s1, s2, s3;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.method = QAWS_GENERATE_PARAMETRIC;
	s1 = qaws_mesh_generate_surface(s, &d, &out, NULL);
	s2 = qaws_mesh_generate_surface(s, NULL, &out, NULL);
	d.method = QAWS_GENERATE_REMESH;
	d.target_length = (qaws_scalar)-1;
	s3 = qaws_mesh_generate_surface(s, &d, &out, NULL);
	sprintf(msg, "refusals: PARAMETRIC %s, no description %s, a negative length %s", qaws_status_to_string(s1), qaws_status_to_string(s2),
		qaws_status_to_string(s3));
	TEST_ASSERT(s1 == QAWS_STATUS_UNSUPPORTED_OPERATION && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT && !out, msg);
	qaws_surface_destroy(s);
}

int test_126_mesh_generate_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 126: meshes of exact geometry\n");
	test_patch();
	test_torus();
	test_box();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
