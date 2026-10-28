/*
 * Test 126: meshes of exact geometry (qaws_mesh_generate_surface/_shell)
 *
 *   - a bicubic B-spline patch in triangles: valid, one border loop whose
 *     vertices are on the domain's sides, every vertex the surface at its
 *     footpoint (face 0, u, v), edges near the target; in quads: all quads;
 *     adaptive by a chord tolerance: edges graded within the bounds;
 *   - a rational NURBS torus in quads: closed, genus 1, every vertex on the
 *     torus, a grid along its parallels (no irregular vertex);
 *   - a box of six bilinear faces as a shell, in quads: closed, genus 0,
 *     every face on the box face its "face" attribute names, the eight
 *     corners vertices, no vertex near an edge but off it, grids on the
 *     faces (irregular vertices only at the corners);
 *   - STRUCTURED: the torus in quads, a grid (closed, genus 1, no irregular
 *     vertex, on the torus exactly); with a chord tolerance, every cell's
 *     centre within it of its corners' mean; the box shell in quads, its
 *     faces' grids meeting (closed, genus 0, the 8 corners its only irregular
 *     vertices); the patch in triangles; the patch to a chord of 0.002, its
 *     lines placed by need: few triangles, every one's centre within the
 *     chord of the patch;
 *   - PARAMETRIC (triangles meshed in (u, v) in the surface's metric): the
 *     torus and the box shell closed, every vertex on its face exactly, the
 *     edges near the target in space, angles fair; the patch to a chord
 *     tolerance, finer where it bends;
 *   - RESTRICTED_DELAUNAY: the torus, the box shell and the patch to a chord
 *     tolerance: closed (but the patch), every vertex on its face, no edge
 *     off the features facing angles summing past pi (Delaunay), every
 *     angle at least 30 degrees but between two feature edges (Chew);
 *   - trimmed faces: a square plate with a round hole and the disk filling
 *     it (another plane's parametrization, the circle their shared edge), by
 *     PARAMETRIC and RESTRICTED_DELAUNAY: valid, the plate's outer square
 *     the only border (the circle shared), every vertex on the plane, the
 *     circle's vertices on the circle; STRUCTURED on a plane trimmed to an
 *     annular sector (four trim pieces: its grid the Coons map of them): all
 *     quads, no irregular vertex, the inner arc's vertices on its circle; on
 *     the plate with a hole (not four-sided) refused;
 *   - trimmed faces of a surface closed in u: a cylinder trimmed by its rims
 *     (loops wrapping round, a band) with a hole across its seam, filled by
 *     a disk of the same surface, by PARAMETRIC and RESTRICTED_DELAUNAY:
 *     valid, the rims the only borders, every vertex on the cylinder at
 *     its footpoint; STRUCTURED on a band between two rims (the top one
 *     starting a period and a quarter on) and on a four-sided patch across
 *     the seam (a piece given a period on): all quads, one piece, the rims
 *     (the patch's loop) the only borders, every vertex on the cylinder at
 *     its footpoint, no edge across the cylinder;
 *   - footpoints through a Boolean: the union of two generated boxes keeps
 *     (face, u, v) footpoints, each on its operand's face at the vertex,
 *     known at the cuts' vertices too;
 *   - refusals: RESTRICTED_DELAUNAY and PARAMETRIC in quads, no
 *     description, a negative length.
 *
 * By default the torus at 0.4 and the box at 0.25 (written at 0.25 and
 * 0.15, run so with QAWS_TEST_HEAVY=1).
 */

#include "test_common.h"
#include "qaws_mesh_generate.h"
#include "qaws_mesh_boolean.h"
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
	qaws_surface_create_bspline(&s, &d);
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
	qaws_surface_create_nurbs(&s, &d);
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
		qaws_surface_create_bilinear(&f[i], &d);
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
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < 3 * t->nv; i++)
		t->x[i] = (double)p[i];
	free(p);
	qaws_mesh_3d_get_faces(t->fv, t->fs, m);
	if (qaws_mesh_3d_get_footpoints(t->foot, m) != QAWS_STATUS_OK)
		for (i = 0; i < t->nv; i++)
			t->foot[i].patch = QAWS_MESH_NONE;
	qaws_mesh_3d_get_topology(&t->top, m);
	qaws_mesh_3d_measure_quality(&t->q, m);
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
		qaws_surface_evaluate(&e, s, t->foot[i].a, t->foot[i].b, QAWS_SURFACE_EVAL_POSITION);
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
	s1 = qaws_mesh_generate_surface(&tri, &rep, s, &d);
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
		qaws_mesh_3d_get_face_adjacency(across, tri);
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
	s2 = qaws_mesh_generate_surface(&quad, NULL, s, &d);
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
	s3 = qaws_mesh_generate_surface(&ada, &rep, s, &d);
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
	double const L = test_heavy() ? 0.25 : 0.4;   /* (by default 0.4: 2.5 times fewer quads; the face count follows) */
	memset(&d, 0, sizeof d);
	d.elements = QAWS_ELEMENTS_QUADS;
	d.target_length = (qaws_scalar)L;
	st = qaws_mesh_generate_surface(&m, NULL, s, &d);
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
		sprintf(msg, "torus quads: %u faces (%u quads), %u irregular vertices, closed %d genus %u, off the torus %.2e, footpoints to %.2e", t.nf,
			t.q.quads, t.q.irregular, t.top.closed, t.top.genus, off, err);
		TEST_ASSERT(qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.quads == t.nf && t.nf > 300 * pow(0.25 / L, 2) && t.q.irregular == 0 && t.top.closed
			&& t.top.genus == 1 && off < TEST_BOUND(1e-9, 32) && err < TEST_BOUND(1e-9, 32), msg);
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
	double const L = test_heavy() ? 0.15 : 0.25;   /* (by default 0.25: 2.8 times fewer quads; the face count and the edge band follow) */
	unsigned int i;
	tg_box(box, h);
	st = qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	memset(&d, 0, sizeof d);
	d.elements = QAWS_ELEMENTS_QUADS;
	d.target_length = (qaws_scalar)L;
	if (st == QAWS_STATUS_OK)
		st = qaws_mesh_generate_shell(&m, &rep, atlas, &d);
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
		if (qaws_mesh_3d_find_attribute(&fa, m, "face") == QAWS_STATUS_OK)
		{
			fid = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (t.nf + 1));
			qaws_mesh_3d_get_attribute(fid, m, fa);
			for (f = 0, c0 = 0; f < t.nf; c0 += t.fs[f], f++)
			{
				unsigned int F = (unsigned int)fid[f];
				for (k = 0; k < t.fs[f]; k++)
					offface += F >= 6 || fabs(t.x[3 * t.fv[c0 + k] + tg_axis[F]] - tg_side[F] * h[tg_axis[F]]) > TEST_BOUND(1e-9, 32);
			}
			free(fid);
		}
		/* the corners; vertices near an edge (two coordinates near their extremes) but off it */
		for (a = 0; a < 8; a++)
			for (i = 0; i < t.nv; i++)
			{
				double const* p = t.x + 3 * i;
				if (fabs(p[0] - (a & 1 ? h[0] : -h[0])) + fabs(p[1] - (a & 2 ? h[1] : -h[1])) + fabs(p[2] - (a & 4 ? h[2] : -h[2])) < TEST_BOUND(1e-12, 32))
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
				on += gap < TEST_BOUND(1e-9, 32);
				close += gap < 0.3 * L;
			}
			near += close >= 2 && on < 2;
		}
		sprintf(msg, "box quads: %u faces (%u quads), %u irregular vertices, closed %d genus %u, %u corners of 8, %u vertices off their face, "
			"%u near an edge but off it, footpoints to %.2e", t.nf, t.q.quads, t.q.irregular, t.top.closed, t.top.genus, corners, offface, near, err);
		TEST_ASSERT(qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.quads == t.nf && t.nf > 200 * pow(0.15 / L, 2) && t.q.irregular == 8 && t.top.closed && t.top.genus == 0
			&& fa != QAWS_MESH_NONE && offface == 0 && corners == 8 && near == 0 && err < TEST_BOUND(1e-9, 32), msg);
		tg_free(&t);
	}
	qaws_mesh_3d_destroy(m);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
}

/* STRUCTURED on one surface or a shell: the mesh read, checked by the caller */
static qaws_mesh_3d* tg_structured(qaws_surface const* s, qaws_shell_atlas const* a, qaws_elements el, double L, double chord, tg_mesh* t)
{
	qaws_generate_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.method = QAWS_GENERATE_STRUCTURED;
	d.elements = el;
	d.target_length = (qaws_scalar)L;
	d.chord_tolerance = (qaws_scalar)chord;
	if ((s ? qaws_mesh_generate_surface(&m, NULL, s, &d) : qaws_mesh_generate_shell(&m, NULL, a, &d)) == QAWS_STATUS_OK && m)
		tg_read(m, t);
	else
		memset(t, 0, sizeof *t);
	return m;
}

static void test_structured(void)
{
	static double const h[3] = { 1, 0.7, 0.5 };
	qaws_surface *torus = tg_torus(2, 0.7), *patch = tg_patch(), *box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d* m;
	tg_mesh t;
	char msg[400];
	unsigned int i;
	double off = 0, err, chord = 0, pu, pv;
	/* the torus, a grid */
	m = tg_structured(torus, NULL, QAWS_ELEMENTS_QUADS, 0.4, 0, &t);
	err = m ? tg_foot_error(&t, torus, NULL) : HUGE_VAL;
	for (i = 0; m && i < t.nv; i++)
	{
		double const* q = t.x + 3 * i;
		double rho = sqrt(q[0] * q[0] + q[1] * q[1]);
		off = fmax(off, fabs(sqrt((rho - 2) * (rho - 2) + q[2] * q[2]) - 0.7));
	}
	sprintf(msg, "STRUCTURED torus: %u faces (%u quads), %u irregular, closed %d genus %u, off %.2e, footpoints %.2e, edges %.3g..%.3g", t.nf, t.q.quads,
		t.q.irregular, t.top.closed, t.top.genus, off, err, (double)t.q.min_edge, (double)t.q.max_edge);
	printf("  %s\n", msg);
	TEST_ASSERT(m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.quads == t.nf && t.q.irregular == 0 && t.top.closed && t.top.genus == 1
		&& off < TEST_BOUND(1e-9, 32) && err < TEST_BOUND(1e-9, 32) && t.q.max_edge < 0.42, msg);
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	/* to a chord tolerance: each quad's centre (by its corners' footpoints) on the torus near its corners' mean */
	m = tg_structured(torus, NULL, QAWS_ELEMENTS_QUADS, 0, 0.01, &t);
	pu = (double)qaws_surface_get_u_range(torus).max_value - (double)qaws_surface_get_u_range(torus).min_value;
	pv = (double)qaws_surface_get_v_range(torus).max_value - (double)qaws_surface_get_v_range(torus).min_value;
	for (i = 0; m && i < t.nf; i++)
	{
		unsigned int const* q = t.fv + 4 * i;
		double mid[3] = { 0, 0, 0 }, u = 0, v = 0, y[3], dd = 0;
		qaws_surface_eval_result e;
		unsigned int k, z;
		if (t.fs[i] != 4)
			continue;
		for (z = 0; z < 4; z++)
		{
			double du = (double)t.foot[q[z]].a - (double)t.foot[q[0]].a, dv = (double)t.foot[q[z]].b - (double)t.foot[q[0]].b;
			/* (corners across the seam: their parameters brought next to the first's) */
			du -= pu * floor(du / pu + 0.5);
			dv -= pv * floor(dv / pv + 0.5);
			u += 0.25 * ((double)t.foot[q[0]].a + du);
			v += 0.25 * ((double)t.foot[q[0]].b + dv);
			for (k = 0; k < 3; k++)
				mid[k] += 0.25 * t.x[3 * q[z] + k];
		}
		memset(&e, 0, sizeof e);
		qaws_surface_evaluate(&e, torus, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION);
		y[0] = (double)e.position.x; y[1] = (double)e.position.y; y[2] = (double)e.position.z;
		for (k = 0; k < 3; k++)
			dd += (y[k] - mid[k]) * (y[k] - mid[k]);
		chord = fmax(chord, sqrt(dd));
	}
	sprintf(msg, "STRUCTURED torus to a chord of 0.01: %u faces, closed %d, the cells' centres off their corners' mean by %.4f at most", t.nf,
		t.top.closed, chord);
	printf("  %s\n", msg);
	TEST_ASSERT(m && t.top.closed && t.top.genus == 1 && t.q.quads == t.nf && chord <= 0.01 && chord > 0.001, msg);
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	/* the box shell */
	tg_box(box, h);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	m = atlas ? tg_structured(NULL, atlas, QAWS_ELEMENTS_QUADS, 0.25, 0, &t) : NULL;
	err = m ? tg_foot_error(&t, NULL, atlas) : HUGE_VAL;
	sprintf(msg, "STRUCTURED box shell: %u faces (%u quads), %u irregular, closed %d genus %u, footpoints %.2e, edges %.3g..%.3g", t.nf, t.q.quads,
		t.q.irregular, t.top.closed, t.top.genus, err, (double)t.q.min_edge, (double)t.q.max_edge);
	printf("  %s\n", msg);
	TEST_ASSERT(m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.quads == t.nf && t.q.irregular == 8 && t.top.closed && t.top.genus == 0
		&& err < TEST_BOUND(1e-9, 32) && t.q.max_edge <= 0.25 + 1e-9, msg);
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	/* the patch to a chord of 0.002: lines by need */
	m = tg_structured(patch, NULL, QAWS_ELEMENTS_TRIANGLES, 0, 0.002, &t);
	{
		double worst = 0;
		unsigned int c0 = 0, z;
		for (i = 0; m && i < t.nf; c0 += t.fs[i], i++)
		{
			double u = 0, v = 0, c[3] = { 0, 0, 0 }, dd = 0;
			qaws_surface_eval_result e;
			for (z = 0; z < 3; z++)
			{
				unsigned int q = t.fv[c0 + z];
				u += (double)t.foot[q].a / 3;
				v += (double)t.foot[q].b / 3;
				c[0] += t.x[3 * q] / 3; c[1] += t.x[3 * q + 1] / 3; c[2] += t.x[3 * q + 2] / 3;
			}
			memset(&e, 0, sizeof e);
			qaws_surface_evaluate(&e, patch, (qaws_scalar)u, (qaws_scalar)v, QAWS_SURFACE_EVAL_POSITION);
			dd = ((double)e.position.x - c[0]) * ((double)e.position.x - c[0]) + ((double)e.position.y - c[1]) * ((double)e.position.y - c[1])
				+ ((double)e.position.z - c[2]) * ((double)e.position.z - c[2]);
			worst = fmax(worst, sqrt(dd));
		}
		sprintf(msg, "STRUCTURED patch to a chord of 0.002: %u triangles, the centres off the patch by %.4f at most", t.nf, worst);
		printf("  %s\n", msg);
		TEST_ASSERT(m && t.q.triangles == t.nf && t.nf < 8000 && worst <= 0.002, msg);
	}
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	/* the patch in triangles */
	m = tg_structured(patch, NULL, QAWS_ELEMENTS_TRIANGLES, 0.3, 0, &t);
	err = m ? tg_foot_error(&t, patch, NULL) : HUGE_VAL;
	sprintf(msg, "STRUCTURED patch in triangles: %u triangles, valid %d, one border %d, footpoints %.2e, angles %.3g..%.3g", t.nf,
		m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, t.top.boundary_loops, err, (double)t.q.min_angle, (double)t.q.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.triangles == t.nf && t.top.boundary_loops == 1 && err < TEST_BOUND(1e-9, 32), msg);
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
	qaws_surface_destroy(torus);
	qaws_surface_destroy(patch);
}

static qaws_mesh_3d* tg_parametric(qaws_surface const* s, qaws_shell_atlas const* a, double L, double chord, tg_mesh* t)
{
	qaws_generate_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.method = QAWS_GENERATE_PARAMETRIC;
	d.target_length = (qaws_scalar)L;
	d.chord_tolerance = (qaws_scalar)chord;
	d.min_length = (qaws_scalar)(chord > 0 ? 0.02 : 0);
	d.max_length = (qaws_scalar)(chord > 0 ? 0.3 : 0);
	if ((s ? qaws_mesh_generate_surface(&m, NULL, s, &d) : qaws_mesh_generate_shell(&m, NULL, a, &d)) == QAWS_STATUS_OK && m)
		tg_read(m, t);
	else
		memset(t, 0, sizeof *t);
	return m;
}

static void test_parametric(void)
{
	static double const h[3] = { 1, 0.7, 0.5 };
	qaws_surface *torus = tg_torus(2, 0.7), *patch = tg_patch(), *box[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d* m;
	tg_mesh t;
	char msg[400];
	unsigned int i;
	double err;
	m = tg_parametric(torus, NULL, 0.3, 0, &t);
	err = m ? tg_foot_error(&t, torus, NULL) : HUGE_VAL;
	sprintf(msg, "PARAMETRIC torus to 0.3: %u triangles, valid %d, closed %d genus %u, footpoints %.2e, edges %.3g..%.3g mean %.3g, angles %.3g..%.3g",
		t.nf, m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, t.top.closed, t.top.genus, err, (double)t.q.min_edge, (double)t.q.max_edge,
		(double)t.q.mean_edge, (double)t.q.min_angle, (double)t.q.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.q.triangles == t.nf && t.top.closed && t.top.genus == 1 && err < TEST_BOUND(1e-9, 32)
		&& fabs(t.q.mean_edge - 0.3) < 0.08 && t.q.min_angle > 20, msg);
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	tg_box(box, h);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	m = atlas ? tg_parametric(NULL, atlas, 0.25, 0, &t) : NULL;
	err = m ? tg_foot_error(&t, NULL, atlas) : HUGE_VAL;
	sprintf(msg, "PARAMETRIC box shell to 0.25: %u triangles, valid %d, closed %d genus %u, footpoints %.2e, edges %.3g..%.3g mean %.3g, angles %.3g..%.3g",
		t.nf, m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, t.top.closed, t.top.genus, err, (double)t.q.min_edge, (double)t.q.max_edge,
		(double)t.q.mean_edge, (double)t.q.min_angle, (double)t.q.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && t.top.closed && t.top.genus == 0 && err < TEST_BOUND(1e-9, 32)
		&& fabs(t.q.mean_edge - 0.25) < 0.07 && t.q.min_angle > 20, msg);
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	m = tg_parametric(patch, NULL, 0, 0.005, &t);
	err = m ? tg_foot_error(&t, patch, NULL) : HUGE_VAL;
	sprintf(msg, "PARAMETRIC patch to a chord of 0.005: %u triangles, valid %d, footpoints %.2e, edges %.3g..%.3g, angles %.3g..%.3g", t.nf,
		m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, err, (double)t.q.min_edge, (double)t.q.max_edge, (double)t.q.min_angle, (double)t.q.max_angle);
	printf("  %s\n", msg);
	TEST_ASSERT(m && qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && err < TEST_BOUND(1e-9, 32) && t.q.max_edge > 2 * t.q.min_edge, msg);
	if (m)
		tg_free(&t);
	qaws_mesh_3d_destroy(m);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
	qaws_surface_destroy(torus);
	qaws_surface_destroy(patch);
}

/* RESTRICTED_DELAUNAY's checks: the edges off the features whose facing angles sum past pi, the smallest angle not
   between two feature edges */
static void tg_delaunay_check(qaws_mesh_3d const* m, unsigned int* bad_edges, double* min_angle)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), ne = qaws_mesh_3d_get_feature_edge_count(m), i, j, k, *t, *fe;
	unsigned long long* key;
	double* x;
	qaws_scalar* p;
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m);
	t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt + 1));
	fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne + 1));
	key = (unsigned long long*)malloc(sizeof(unsigned long long) * 3 * (nt + 1));
	x = (double*)malloc(sizeof(double) * 3 * (nv + 1));
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv + 1));
	qaws_mesh_3d_get_triangles(t, m);
	qaws_mesh_3d_get_feature_edges(fe, m);
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < 3 * nv; i++)
		x[i] = (double)p[i];
	*bad_edges = 0;
	*min_angle = 180;
	/* (each triangle's corners' angles; each edge's facing angle, keyed by its sorted ends) */
	for (i = 0; i < nt; i++)
		for (k = 0; k < 3; k++)
		{
			unsigned int a = t[3 * i + k], b = t[3 * i + (k + 1) % 3], c = t[3 * i + (k + 2) % 3], lo = a < b ? a : b, hi = a < b ? b : a, fa = 0, fb = 0;
			double u[3], w[3], ang;
			for (j = 0; j < 3; j++)
			{
				u[j] = x[3 * a + j] - x[3 * c + j];
				w[j] = x[3 * b + j] - x[3 * c + j];
			}
			ang = acos(fmax(-1, fmin(1, (u[0] * w[0] + u[1] * w[1] + u[2] * w[2]) / sqrt((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2])))));
			/* (the angle at c is between edges c-a and c-b: kept off the bound when both are features) */
			for (j = 0; j < ne; j++)
			{
				fa |= (fe[2 * j] == c && fe[2 * j + 1] == a) || (fe[2 * j] == a && fe[2 * j + 1] == c);
				fb |= (fe[2 * j] == c && fe[2 * j + 1] == b) || (fe[2 * j] == b && fe[2 * j + 1] == c);
			}
			if (!(fa && fb))
				*min_angle = fmin(*min_angle, ang * 180 / 3.14159265358979);
			key[3 * i + k] = ((unsigned long long)lo << 32 | hi);
			(void)ang;
		}
	for (i = 0; i < 3 * nt; i++)
		for (j = i + 1; j < 3 * nt; j++)
			if (key[i] == key[j])
			{
				unsigned int a = (unsigned int)(key[i] >> 32), b = (unsigned int)(key[i] & 0xffffffffu), q, feature = 0;
				double s = 0;
				unsigned int sides[2], z;
				sides[0] = i;
				sides[1] = j;
				for (q = 0; q < ne; q++)
					feature |= (fe[2 * q] == a && fe[2 * q + 1] == b) || (fe[2 * q] == b && fe[2 * q + 1] == a);
				if (feature)
					break;
				for (z = 0; z < 2; z++)
				{
					unsigned int tri = sides[z] / 3, kk = sides[z] % 3, c = t[3 * tri + (kk + 2) % 3];
					double u[3], w[3];
					for (q = 0; q < 3; q++)
					{
						u[q] = x[3 * a + q] - x[3 * c + q];
						w[q] = x[3 * b + q] - x[3 * c + q];
					}
					s += acos(fmax(-1, fmin(1, (u[0] * w[0] + u[1] * w[1] + u[2] * w[2]) / sqrt((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2])))));
				}
				*bad_edges += s > 3.14159265358979 + 1e-9;
				break;
			}
	free(t); free(fe); free(key); free(x); free(p);
}

static void test_restricted_delaunay(void)
{
	static double const h[3] = { 1, 0.7, 0.5 };
	qaws_surface *torus = tg_torus(2, 0.7), *patch = tg_patch(), *box[6];
	qaws_shell_atlas* atlas = NULL;
	unsigned int i, bad;
	double amin, err;
	char msg[400];
	tg_box(box, h);
	qaws_shell_atlas_create(&atlas, (qaws_surface const* const*)box, 6, 0);
	for (i = 0; i < 3; i++)
	{
		qaws_generate_desc d;
		qaws_mesh_3d* m = NULL;
		tg_mesh t;
		char const* name = i == 0 ? "torus" : i == 1 ? "box shell" : "patch to a chord of 0.005";
		memset(&d, 0, sizeof d);
		d.method = QAWS_GENERATE_RESTRICTED_DELAUNAY;
		d.target_length = (qaws_scalar)(i == 0 ? 0.3 : i == 1 ? 0.25 : 0);
		d.chord_tolerance = (qaws_scalar)(i == 2 ? 0.005 : 0);
		d.min_length = (qaws_scalar)(i == 2 ? 0.02 : 0);
		d.max_length = (qaws_scalar)(i == 2 ? 0.3 : 0);
		if ((i == 1 ? qaws_mesh_generate_shell(&m, NULL, atlas, &d) : qaws_mesh_generate_surface(&m, NULL, i ? patch : torus, &d)) != QAWS_STATUS_OK || !m)
		{
			sprintf(msg, "RESTRICTED_DELAUNAY %s: generated", name);
			TEST_ASSERT(0, msg);
			continue;
		}
		tg_read(m, &t);
		err = tg_foot_error(&t, i == 1 ? NULL : i ? patch : torus, i == 1 ? atlas : NULL);
		tg_delaunay_check(m, &bad, &amin);
		sprintf(msg, "RESTRICTED_DELAUNAY %s: %u triangles, valid %d, closed %d, footpoints %.2e, %u edges not Delaunay, smallest angle %.2f",
			name, t.nf, qaws_mesh_3d_validate(m) == QAWS_STATUS_OK, t.top.closed, err, bad, amin);
		printf("  %s\n", msg);
		TEST_ASSERT(qaws_mesh_3d_validate(m) == QAWS_STATUS_OK && (i == 2 || t.top.closed) && err < TEST_BOUND(1e-9, 32) && bad == 0 && amin >= 30 - 1e-6, msg);
		tg_free(&t);
		qaws_mesh_3d_destroy(m);
	}
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(box[i]);
	qaws_surface_destroy(torus);
	qaws_surface_destroy(patch);
}

static qaws_curve* tg_line(double x0, double y0, double x1, double y1)
{
	qaws_bezier_desc bd;
	qaws_scalar cp[4];
	qaws_curve* c = NULL;
	cp[0] = (qaws_scalar)x0; cp[1] = (qaws_scalar)y0; cp[2] = (qaws_scalar)x1; cp[3] = (qaws_scalar)y1;
	memset(&bd, 0, sizeof bd);
	bd.dimension = QAWS_DIMENSION_2D; bd.degree = 1; bd.control_points = cp; bd.control_point_count = 2;
	qaws_curve_create_bezier(&c, &bd);
	return c;
}

static qaws_curve* tg_arc(double cx, double cy, double r, double a0, double a1)
{
	qaws_arc_segment g;
	qaws_arc_desc ad;
	qaws_curve* c = NULL;
	memset(&g, 0, sizeof g);
	g.center[0] = (qaws_scalar)cx; g.center[1] = (qaws_scalar)cy; g.radius = (qaws_scalar)r;
	g.angle_start = (qaws_scalar)a0; g.angle_end = (qaws_scalar)a1;
	memset(&ad, 0, sizeof ad);
	ad.dimension = QAWS_DIMENSION_2D; ad.segments = &g; ad.segment_count = 1;
	qaws_curve_create_arc(&c, &ad);
	return c;
}

static void test_trimmed(void)
{
	/* the plate: z = 0 as (x, y) = (-2 + 4u, -2 + 4v), trimmed to [0.125, 0.875]^2 less the circle of 0.2 round (0.5, 0.5)
	   (radius 0.8 in space); the disk: its own plane (x, y) = (1 - 2v, -1 + 2u), the circle of 0.4 round (0.5, 0.5) */
	static double const A[4][3] = { { -2, -2, 0 }, { 2, -2, 0 }, { -2, 2, 0 }, { 2, 2, 0 } }, B[4][3] = { { 1, -1, 0 }, { 1, 1, 0 }, { -1, -1, 0 }, { -1, 1, 0 } };
	qaws_surface* s[2];
	qaws_curve *oc[4], *hole, *disk;
	qaws_curve const *ol[4], *hl[1], *dl[1];
	qaws_path_2d pa[2], pb[1];
	qaws_shell_face_desc fd[2];
	qaws_shell_atlas* a = NULL;
	unsigned int i, k, m;
	char msg[400];
	for (k = 0; k < 2; k++)
	{
		qaws_surface_bilinear_desc d;
		double const (*P)[3] = k ? B : A;
		memset(&d, 0, sizeof d);
		d.p00.x = (qaws_scalar)P[0][0]; d.p00.y = (qaws_scalar)P[0][1];
		d.p10.x = (qaws_scalar)P[1][0]; d.p10.y = (qaws_scalar)P[1][1];
		d.p01.x = (qaws_scalar)P[2][0]; d.p01.y = (qaws_scalar)P[2][1];
		d.p11.x = (qaws_scalar)P[3][0]; d.p11.y = (qaws_scalar)P[3][1];
		s[k] = NULL;
		qaws_surface_create_bilinear(&s[k], &d);
	}
	oc[0] = tg_line(0.125, 0.125, 0.875, 0.125);
	oc[1] = tg_line(0.875, 0.125, 0.875, 0.875);
	oc[2] = tg_line(0.875, 0.875, 0.125, 0.875);
	oc[3] = tg_line(0.125, 0.875, 0.125, 0.125);
	hole = tg_arc(0.5, 0.5, 0.2, 0, 2 * 3.14159265358979323846);
	disk = tg_arc(0.5, 0.5, 0.4, 1.0, 1.0 + 2 * 3.14159265358979323846);
	for (i = 0; i < 4; i++)
		ol[i] = oc[i];
	hl[0] = hole;
	dl[0] = disk;
	pa[0].curves = ol; pa[0].curve_count = 4; pa[0].closed = 1;
	pa[1].curves = hl; pa[1].curve_count = 1; pa[1].closed = 1;
	pb[0].curves = dl; pb[0].curve_count = 1; pb[0].closed = 1;
	fd[0].surface = s[0]; fd[0].trims = pa; fd[0].trim_count = 2;
	fd[1].surface = s[1]; fd[1].trims = pb; fd[1].trim_count = 1;
	qaws_shell_atlas_create_trimmed(&a, fd, 2, 0);
	for (m = 0; a && m < 2; m++)
	{
		qaws_generate_desc d;
		qaws_mesh_3d* out = NULL;
		qaws_status st;
		tg_mesh t;
		double offz = 0, oncircle = HUGE_VAL, err;
		unsigned int near = 0;
		memset(&d, 0, sizeof d);
		d.method = m ? QAWS_GENERATE_RESTRICTED_DELAUNAY : QAWS_GENERATE_PARAMETRIC;
		d.target_length = (qaws_scalar)0.25;
		st = qaws_mesh_generate_shell(&out, NULL, a, &d);
		if (st != QAWS_STATUS_OK || !out)
		{
			sprintf(msg, "trimmed plate and disk by %s: %s", m ? "RESTRICTED_DELAUNAY" : "PARAMETRIC", qaws_status_to_string(st));
			TEST_ASSERT(0, msg);
			continue;
		}
		tg_read(out, &t);
		err = tg_foot_error(&t, NULL, a);
		oncircle = 0;
		for (i = 0; i < t.nv; i++)
		{
			double r = sqrt(t.x[3 * i] * t.x[3 * i] + t.x[3 * i + 1] * t.x[3 * i + 1]);
			offz = fmax(offz, fabs(t.x[3 * i + 2]));
			/* (the circle's vertices: within a tenth of the size of radius 0.8, then on it) */
			if (fabs(r - 0.8) < 0.025)
			{
				near++;
				oncircle = fmax(oncircle, fabs(r - 0.8));
			}
		}
		sprintf(msg, "trimmed plate and disk by %s: %u triangles, valid %d, border loops %u, genus %u, off the plane %.1e, footpoints %.1e, "
			"%u vertices on the circle (to %.1e), angles %.3g..%.3g", m ? "RESTRICTED_DELAUNAY" : "PARAMETRIC", t.nf, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK,
			t.top.boundary_loops, t.top.genus, offz, err, near, oncircle, (double)t.q.min_angle, (double)t.q.max_angle);
		printf("  %s\n", msg);
		TEST_ASSERT(qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.top.boundary_loops == 1 && t.top.genus == 0 && t.top.components == 1 && offz < 1e-12
			&& err < TEST_BOUND(1e-9, 32) && near >= 16 && oncircle < 1e-6 && t.q.min_angle > (m ? 29.99 : 20), msg);
		tg_free(&t);
		qaws_mesh_3d_destroy(out);
	}
	if (a)
	{
		/* STRUCTURED: the plate (not four-sided) refused */
		qaws_generate_desc d;
		qaws_mesh_3d* out = NULL;
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.method = QAWS_GENERATE_STRUCTURED;
		d.target_length = (qaws_scalar)0.25;
		st = qaws_mesh_generate_shell(&out, NULL, a, &d);
		sprintf(msg, "STRUCTURED on the plate with a hole: %s", qaws_status_to_string(st));
		TEST_ASSERT(st == QAWS_STATUS_UNSUPPORTED_OPERATION && !out, msg);
	}
	{
		/* an annular sector of the plate's plane: radii 0.1 and 0.3 round (0.5, 0.5) (0.4 and 1.2 in space), a quarter turn */
		qaws_curve *sc[4];
		qaws_curve const* sl[4];
		qaws_path_2d ps[1];
		qaws_shell_face_desc fs1[1];
		qaws_shell_atlas* as = NULL;
		qaws_generate_desc d;
		qaws_mesh_3d* out = NULL;
		qaws_status st;
		double hp = 0.5 * 3.14159265358979323846;
		sc[0] = tg_line(0.6, 0.5, 0.8, 0.5);
		sc[1] = tg_arc(0.5, 0.5, 0.3, 0, hp);
		sc[2] = tg_line(0.5, 0.8, 0.5, 0.6);
		sc[3] = tg_arc(0.5, 0.5, 0.1, hp, 0);
		for (i = 0; i < 4; i++)
			sl[i] = sc[i];
		ps[0].curves = sl; ps[0].curve_count = 4; ps[0].closed = 1;
		fs1[0].surface = s[0]; fs1[0].trims = ps; fs1[0].trim_count = 1;
		qaws_shell_atlas_create_trimmed(&as, fs1, 1, 0);
		memset(&d, 0, sizeof d);
		d.method = QAWS_GENERATE_STRUCTURED;
		d.elements = QAWS_ELEMENTS_QUADS;
		d.target_length = (qaws_scalar)0.1;
		st = as ? qaws_mesh_generate_shell(&out, NULL, as, &d) : QAWS_STATUS_INVALID_ARGUMENT;
		if (st == QAWS_STATUS_OK && out)
		{
			tg_mesh t;
			double inner = 0, err;
			unsigned int on = 0;
			tg_read(out, &t);
			err = tg_foot_error(&t, NULL, as);
			for (i = 0; i < t.nv; i++)
			{
				double r = sqrt(t.x[3 * i] * t.x[3 * i] + t.x[3 * i + 1] * t.x[3 * i + 1]);
				if (fabs(r - 0.4) < 0.01)
				{
					on++;
					inner = fmax(inner, fabs(r - 0.4));
				}
			}
			sprintf(msg, "STRUCTURED annular sector: %u faces (%u quads), %u irregular, border loops %u, footpoints %.1e, %u on the inner arc (to %.1e)",
				t.nf, t.q.quads, t.q.irregular, t.top.boundary_loops, err, on, inner);
			printf("  %s\n", msg);
			TEST_ASSERT(qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.q.quads == t.nf && t.q.irregular == 0 && t.top.boundary_loops == 1
				&& err < TEST_BOUND(1e-9, 32) && on >= 3 && inner < 1e-9, msg);
			tg_free(&t);
		}
		else
		{
			sprintf(msg, "STRUCTURED annular sector: %s", qaws_status_to_string(st));
			TEST_ASSERT(0, msg);
		}
		qaws_mesh_3d_destroy(out);
		qaws_shell_atlas_destroy(as);
		for (i = 0; i < 4; i++)
			qaws_curve_destroy(sc[i]);
	}
	qaws_shell_atlas_destroy(a);
	for (i = 0; i < 4; i++)
		qaws_curve_destroy(oc[i]);
	qaws_curve_destroy(hole);
	qaws_curve_destroy(disk);
	qaws_surface_destroy(s[0]);
	qaws_surface_destroy(s[1]);
}

static void test_periodic_trimmed(void)
{
	/* the cylinder of radius 1/2 and height pi round z (u round it, period 1), trimmed by its rims and a hole of radius 0.1
	   on the seam; the disk filling it, its circle starting in another copy */
	qaws_curve *prof = tg_line(0.5, 0, 0.5, -3.14159265358979323846), *bot = tg_line(0, 0, 1, 0), *top = tg_line(1, 1, 0, 1);
	qaws_curve *hole = tg_arc(0, 0.5, 0.1, 0, 2 * 3.14159265358979323846), *disk = tg_arc(1, 0.5, 0.1, 2, 2 + 2 * 3.14159265358979323846);
	qaws_curve const *l1[1], *l2[1], *l3[1], *l4[1];
	qaws_surface* s = NULL;
	qaws_surface_revolution_desc rd;
	qaws_path_2d pa[3], pb[1];
	qaws_shell_face_desc fd[2];
	qaws_shell_atlas* a = NULL;
	unsigned int m, i;
	char msg[400];
	memset(&rd, 0, sizeof rd);
	rd.profile = prof;
	rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&s, &rd);
	l1[0] = bot; l2[0] = top; l3[0] = hole; l4[0] = disk;
	pa[0].curves = l1; pa[0].curve_count = 1; pa[0].closed = 1;
	pa[1].curves = l2; pa[1].curve_count = 1; pa[1].closed = 1;
	pa[2].curves = l3; pa[2].curve_count = 1; pa[2].closed = 1;
	pb[0].curves = l4; pb[0].curve_count = 1; pb[0].closed = 1;
	fd[0].surface = s; fd[0].trims = pa; fd[0].trim_count = 3;
	fd[1].surface = s; fd[1].trims = pb; fd[1].trim_count = 1;
	qaws_shell_atlas_create_trimmed(&a, fd, 2, 0);
	for (m = 0; m < 2; m++)
	{
		qaws_generate_desc d;
		qaws_mesh_3d* out = NULL;
		qaws_status st;
		memset(&d, 0, sizeof d);
		d.method = m ? QAWS_GENERATE_RESTRICTED_DELAUNAY : QAWS_GENERATE_PARAMETRIC;
		d.target_length = (qaws_scalar)0.25;
		st = a ? qaws_mesh_generate_shell(&out, NULL, a, &d) : QAWS_STATUS_INVALID_ARGUMENT;
		if (st == QAWS_STATUS_OK && out)
		{
			tg_mesh t;
			double off = 0, err;
			tg_read(out, &t);
			err = tg_foot_error(&t, NULL, a);
			for (i = 0; i < t.nv; i++)
				off = fmax(off, fabs(sqrt(t.x[3 * i] * t.x[3 * i] + t.x[3 * i + 1] * t.x[3 * i + 1]) - 0.5));
			sprintf(msg, "trimmed cylinder (a band, a hole on the seam) by %s: %u triangles, valid %d, border loops %u, genus %u, components %u, off the cylinder %.1e, "
				"footpoints %.1e, angles %.3g..%.3g", m ? "RESTRICTED_DELAUNAY" : "PARAMETRIC", t.nf, qaws_mesh_3d_validate(out) == QAWS_STATUS_OK, t.top.boundary_loops,
				t.top.genus, t.top.components, off, err, (double)t.q.min_angle, (double)t.q.max_angle);
			printf("  %s\n", msg);
			TEST_ASSERT(qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.top.boundary_loops == 2 && t.top.genus == 0 && t.top.components == 1
				&& off < TEST_BOUND(1e-9, 32) && err < TEST_BOUND(1e-9, 32) && t.q.min_angle > (m ? 29.99 : 20), msg);
			tg_free(&t);
		}
		else
		{
			sprintf(msg, "trimmed cylinder by %s: %s", m ? "RESTRICTED_DELAUNAY" : "PARAMETRIC", qaws_status_to_string(st));
			TEST_ASSERT(0, msg);
		}
		qaws_mesh_3d_destroy(out);
	}
	qaws_shell_atlas_destroy(a);
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof); qaws_curve_destroy(bot); qaws_curve_destroy(top); qaws_curve_destroy(hole); qaws_curve_destroy(disk);
}

/* STRUCTURED on trimmed faces of the cylinder of test_periodic_trimmed: a band between two rims (the top one starting a
   period and a quarter on: the cut joining their starts slants) and a four-sided patch across the seam (one piece given
   a period on) */
static void test_periodic_structured(void)
{
	qaws_curve* prof = tg_line(0.5, 0, 0.5, -3.14159265358979323846);
	qaws_curve* c[6];
	qaws_curve const *l1[1], *l2[1], *l4[4];
	qaws_surface* s = NULL;
	qaws_surface_revolution_desc rd;
	qaws_path_2d pa[2];
	qaws_shell_face_desc fd[1];
	unsigned int m, i;
	char msg[400];
	memset(&rd, 0, sizeof rd);
	rd.profile = prof;
	rd.axis_direction.z = 1;
	qaws_surface_create_revolution(&s, &rd);
	c[0] = tg_line(0, 0.1, 1, 0.1);
	c[1] = tg_line(2.25, 0.9, 1.25, 0.9);
	c[2] = tg_line(-0.2, 0.3, 0.2, 0.3);
	c[3] = tg_line(0.2, 0.3, 0.2, 0.7);
	c[4] = tg_line(1.2, 0.7, 0.8, 0.7);
	c[5] = tg_line(0.8, 0.7, 0.8, 0.3);
	l1[0] = c[0]; l2[0] = c[1];
	for (i = 0; i < 4; i++)
		l4[i] = c[2 + i];
	for (m = 0; m < 2; m++)
	{
		qaws_shell_atlas* a = NULL;
		qaws_generate_desc d;
		qaws_mesh_3d* out = NULL;
		qaws_status st;
		if (m == 0)
		{
			pa[0].curves = l1; pa[0].curve_count = 1; pa[0].closed = 1;
			pa[1].curves = l2; pa[1].curve_count = 1; pa[1].closed = 1;
			fd[0].trim_count = 2;
		}
		else
		{
			pa[0].curves = l4; pa[0].curve_count = 4; pa[0].closed = 1;
			fd[0].trim_count = 1;
		}
		fd[0].surface = s; fd[0].trims = pa;
		qaws_shell_atlas_create_trimmed(&a, fd, 1, 0);
		memset(&d, 0, sizeof d);
		d.method = QAWS_GENERATE_STRUCTURED;
		d.elements = QAWS_ELEMENTS_QUADS;
		d.target_length = (qaws_scalar)0.2;
		st = a ? qaws_mesh_generate_shell(&out, NULL, a, &d) : QAWS_STATUS_INVALID_ARGUMENT;
		if (st == QAWS_STATUS_OK && out)
		{
			tg_mesh t;
			double off = 0, longest = 0, err, z0 = 1e300, z1 = -1e300;
			unsigned int f, k;
			tg_read(out, &t);
			err = tg_foot_error(&t, NULL, a);
			for (i = 0; i < t.nv; i++)
			{
				off = fmax(off, fabs(sqrt(t.x[3 * i] * t.x[3 * i] + t.x[3 * i + 1] * t.x[3 * i + 1]) - 0.5));
				z0 = fmin(z0, t.x[3 * i + 2]);
				z1 = fmax(z1, t.x[3 * i + 2]);
			}
			for (f = 0, k = 0; f < t.nf; k += t.fs[f++])
				for (i = 0; i < t.fs[f]; i++)
				{
					double const *p = t.x + 3 * t.fv[k + i], *q = t.x + 3 * t.fv[k + (i + 1) % t.fs[f]];
					longest = fmax(longest, sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2])));
				}
			sprintf(msg, "STRUCTURED on the cylinder, %s: %u faces (%u quads), %u irregular, border loops %u, genus %u, components %u, off it %.1e, footpoints %.1e, "
				"longest edge %.3f, height %.4f..%.4f", m ? "a patch across the seam" : "a band", t.nf, t.q.quads, t.q.irregular, t.top.boundary_loops, t.top.genus,
				t.top.components, off, err, longest, z0, z1);
			printf("  %s\n", msg);
			TEST_ASSERT(qaws_mesh_3d_validate(out) == QAWS_STATUS_OK && t.q.quads == t.nf && t.top.boundary_loops == (m ? 1u : 2u) && t.top.genus == 0
				&& t.top.components == 1 && off < TEST_BOUND(1e-9, 32) && err < TEST_BOUND(1e-9, 32) && longest < 0.3
				&& fabs(fabs(z1 - z0) - (m ? 0.4 : 0.8) * 3.14159265358979323846) < 1e-9 && (m || t.q.irregular == 0), msg);
			tg_free(&t);
		}
		else
		{
			sprintf(msg, "STRUCTURED on the cylinder, %s: %s", m ? "a patch across the seam" : "a band", qaws_status_to_string(st));
			TEST_ASSERT(0, msg);
		}
		qaws_mesh_3d_destroy(out);
		qaws_shell_atlas_destroy(a);
	}
	qaws_surface_destroy(s);
	qaws_curve_destroy(prof);
	for (i = 0; i < 6; i++)
		qaws_curve_destroy(c[i]);
}

static void test_boolean_footpoints(void)
{
	static double const ha[3] = { 1, 0.7, 0.5 }, hb[3] = { 0.6, 1.2, 0.3 };
	qaws_surface *fa[6], *fb[6];
	qaws_shell_atlas *sa = NULL, *sb = NULL;
	qaws_mesh_3d *ma = NULL, *mb = NULL, *o = NULL;
	qaws_generate_desc d;
	unsigned int i, k, j, known = 0, cut = 0;
	double worst = 0;
	char msg[300];
	tg_box(fa, ha);
	tg_box(fb, hb);
	qaws_shell_atlas_create(&sa, (qaws_surface const* const*)fa, 6, 0);
	qaws_shell_atlas_create(&sb, (qaws_surface const* const*)fb, 6, 0);
	memset(&d, 0, sizeof d);
	d.target_length = (qaws_scalar)0.25;
	if (sa && sb && qaws_mesh_generate_shell(&ma, NULL, sa, &d) == QAWS_STATUS_OK && qaws_mesh_generate_shell(&mb, NULL, sb, &d) == QAWS_STATUS_OK
		&& qaws_mesh_3d_boolean(&o, NULL, ma, mb, QAWS_BOOLEAN_UNION, NULL) == QAWS_STATUS_OK && o)
	{
		tg_mesh t, ta, tb;
		tg_read(o, &t);
		tg_read(ma, &ta);
		tg_read(mb, &tb);
		/* each known footpoint on A's or B's face (the patches numbered alike): the nearer */
		for (i = 0; i < t.nv; i++)
		{
			double best = HUGE_VAL;
			int input = 0;
			if (t.foot[i].patch == QAWS_MESH_NONE)
				continue;
			known++;
			for (k = 0; k < 2; k++)
			{
				qaws_surface const* s = qaws_shell_atlas_get_face(k ? sb : sa, t.foot[i].patch);
				qaws_surface_eval_result e;
				double dd;
				if (!s)
					continue;
				memset(&e, 0, sizeof e);
				qaws_surface_evaluate(&e, s, t.foot[i].a, t.foot[i].b, QAWS_SURFACE_EVAL_POSITION);
				dd = sqrt(((double)e.position.x - t.x[3 * i]) * ((double)e.position.x - t.x[3 * i]) + ((double)e.position.y - t.x[3 * i + 1]) * ((double)e.position.y - t.x[3 * i + 1])
					+ ((double)e.position.z - t.x[3 * i + 2]) * ((double)e.position.z - t.x[3 * i + 2]));
				best = fmin(best, dd);
			}
			worst = fmax(worst, best);
			/* (a cut's vertex: at no input vertex) */
			for (k = 0; k < 2 && !input; k++)
			{
				tg_mesh const* q = k ? &tb : &ta;
				for (j = 0; j < q->nv && !input; j++)
					input = fabs(q->x[3 * j] - t.x[3 * i]) + fabs(q->x[3 * j + 1] - t.x[3 * i + 1]) + fabs(q->x[3 * j + 2] - t.x[3 * i + 2]) < 1e-12;
			}
			cut += !input;
		}
		sprintf(msg, "Boolean footpoints: %u of %u vertices known (%u at the cuts), on their faces to %.1e", known, t.nv, cut, worst);
		printf("  %s\n", msg);
		TEST_ASSERT(known * 10 >= t.nv * 9 && cut > 0 && worst < TEST_BOUND(1e-9, 32), msg);
		tg_free(&t);
		tg_free(&ta);
		tg_free(&tb);
	}
	else
		TEST_ASSERT(0, "Boolean footpoints: the boxes and their union");
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(ma);
	qaws_mesh_3d_destroy(mb);
	qaws_shell_atlas_destroy(sa);
	qaws_shell_atlas_destroy(sb);
	for (i = 0; i < 6; i++)
	{
		qaws_surface_destroy(fa[i]);
		qaws_surface_destroy(fb[i]);
	}
}

static void test_refusals(void)
{
	qaws_surface* s = tg_patch();
	qaws_generate_desc d;
	qaws_mesh_3d* out = NULL;
	qaws_status s1, s2, s3, s4;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.method = QAWS_GENERATE_RESTRICTED_DELAUNAY;
	d.elements = QAWS_ELEMENTS_QUADS;
	s1 = qaws_mesh_generate_surface(&out, NULL, s, &d);
	d.method = QAWS_GENERATE_PARAMETRIC;
	s4 = qaws_mesh_generate_surface(&out, NULL, s, &d);
	d.elements = QAWS_ELEMENTS_TRIANGLES;
	s2 = qaws_mesh_generate_surface(&out, NULL, s, NULL);
	d.method = QAWS_GENERATE_REMESH;
	d.target_length = (qaws_scalar)-1;
	s3 = qaws_mesh_generate_surface(&out, NULL, s, &d);
	sprintf(msg, "refusals: RESTRICTED_DELAUNAY in quads %s, PARAMETRIC in quads %s, no description %s, a negative length %s", qaws_status_to_string(s1),
		qaws_status_to_string(s4), qaws_status_to_string(s2),
		qaws_status_to_string(s3));
	TEST_ASSERT(s1 == QAWS_STATUS_UNSUPPORTED_OPERATION && s4 == QAWS_STATUS_UNSUPPORTED_OPERATION && s2 == QAWS_STATUS_INVALID_ARGUMENT && s3 == QAWS_STATUS_INVALID_ARGUMENT && !out, msg);
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
	test_structured();
	test_parametric();
	test_restricted_delaunay();
	test_trimmed();
	test_periodic_trimmed();
	test_periodic_structured();
	test_boolean_footpoints();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
