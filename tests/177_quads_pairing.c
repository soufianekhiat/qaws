/*
 * Test 177: quads from the input's own faces (QAWS_QUADS_PAIRING and
 * QAWS_QUADS_SPLIT), small meshes
 *
 *   - a grid of squares cut into triangles (alternating diagonals): PAIRING
 *     pairs them all back, every quad a square;
 *   - a sphere's rings: PAIRING leaves only the poles' odd triangles; pure
 *     makes it all quads;
 *   - a box of triangles with its edges as features: no pair across an edge
 *     (every quad flat);
 *   - SPLIT: an octahedron into 24 quads, a cube into 24, the new points on
 *     the reference, the cube's faces flat; the input's vertices kept;
 *   - LAYOUT still refused.
 */

#include "test_common.h"
#include "qaws_quadrangulate.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define QP_PI 3.14159265358979323846

static qaws_mesh_3d* qp_make(qaws_scalar const* x, unsigned int nv, unsigned int const* fv, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = nv;
	d.face_vertices = fv;
	d.face_sizes = fs;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static qaws_mesh_3d* qp_run(qaws_mesh_3d const* m, qaws_quad_method method, int pure, qaws_mesh_report* r)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_3d* o = NULL;
	memset(&d, 0, sizeof d);
	d.method = method;
	d.pure = pure;
	if (qaws_mesh_quadrangulate(&o, r, m, &d) != QAWS_STATUS_OK)
		return NULL;
	return o;
}

/* the largest distance of a quad's corner from its corners' mean plane (Newell normal), over the quads */
static double qp_nonplanar(qaws_mesh_3d const* m)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), f, q = 0, k;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * nc), *fs = (unsigned int*)malloc(sizeof(unsigned int) * nf);
	double worst = 0;
	qaws_mesh_3d_get_faces(fv, fs, m);
	for (f = 0; f < nf; q += fs[f++])
	{
		double n[3] = { 0, 0, 0 }, c[3] = { 0, 0, 0 }, l;
		if (fs[f] != 4)
			continue;
		for (k = 0; k < 4; k++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, fv[q + k]), b = qaws_mesh_3d_get_vertex(m, fv[q + (k + 1) % 4]);
			n[0] += ((double)a.y - b.y) * ((double)a.z + b.z);
			n[1] += ((double)a.z - b.z) * ((double)a.x + b.x);
			n[2] += ((double)a.x - b.x) * ((double)a.y + b.y);
			c[0] += (double)a.x / 4;
			c[1] += (double)a.y / 4;
			c[2] += (double)a.z / 4;
		}
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (k = 0; l > 0 && k < 4; k++)
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, fv[q + k]);
			worst = fmax(worst, fabs(((double)a.x - c[0]) * n[0] + ((double)a.y - c[1]) * n[1] + ((double)a.z - c[2]) * n[2]) / l);
		}
	}
	free(fv);
	free(fs);
	return worst;
}

static void test_grid(void)
{
	unsigned int n = 6, i, j, k = 0, t[6 * 36];
	qaws_scalar x[3 * 49];
	qaws_mesh_3d *m, *o;
	qaws_mesh_report r;
	char msg[300];
	for (j = 0; j <= n; j++)
		for (i = 0; i <= n; i++)
		{
			x[3 * (j * 7 + i)] = (qaws_scalar)i;
			x[3 * (j * 7 + i) + 1] = (qaws_scalar)j;
			x[3 * (j * 7 + i) + 2] = 0;
		}
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++)
		{
			unsigned int a = j * 7 + i, b = a + 1, c = a + 8, d = a + 7;
			if ((i + j) % 2)
			{
				t[k++] = a; t[k++] = b; t[k++] = c;
				t[k++] = a; t[k++] = c; t[k++] = d;
			}
			else
			{
				t[k++] = a; t[k++] = b; t[k++] = d;
				t[k++] = b; t[k++] = c; t[k++] = d;
			}
		}
	m = qp_make(x, 49, t, NULL, 72);
	memset(&r, 0, sizeof r);
	o = qp_run(m, QAWS_QUADS_PAIRING, 0, &r);
	sprintf(msg, "a grid of 72 triangles paired: %u quads, %u triangles, angles %.1f..%.1f, vertices %u (49)", r.quality.quads, r.quality.triangles,
		(double)r.quality.min_angle, (double)r.quality.max_angle, o ? qaws_mesh_3d_get_vertex_count(o) : 0);
	printf("  %s\n", msg);
	TEST_ASSERT(o && r.quality.quads == 36 && r.quality.triangles == 0 && fabs((double)r.quality.min_angle - 90) < 1e-9 && fabs((double)r.quality.max_angle - 90) < 1e-9 &&
		qaws_mesh_3d_get_vertex_count(o) == 49, msg);
	qaws_mesh_3d_destroy(o);
	qaws_mesh_3d_destroy(m);
}

/* a unit sphere: rings between single poles */
static qaws_mesh_3d* qp_sphere(unsigned int nu, unsigned int nv)
{
	unsigned int nn = nu * (nv - 1) + 2, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nn);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * (nv - 1));
	qaws_mesh_3d* m;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = QP_PI * j / nv, ph = 2 * QP_PI * i / nu;
			qaws_scalar* p = x + 3 * ((j - 1) * nu + i);
			p[0] = (qaws_scalar)(sin(th) * cos(ph)); p[1] = (qaws_scalar)(sin(th) * sin(ph)); p[2] = (qaws_scalar)cos(th);
		}
	x[3 * (nn - 2)] = 0; x[3 * (nn - 2) + 1] = 0; x[3 * (nn - 2) + 2] = 1;
	x[3 * (nn - 1)] = 0; x[3 * (nn - 1) + 1] = 0; x[3 * (nn - 1) + 2] = -1;
	for (i = 0; i < nu; i++)
	{
		unsigned int a = i, b = (i + 1) % nu, o = (nv - 2) * nu;
		t[k++] = nn - 2; t[k++] = a; t[k++] = b;
		t[k++] = nn - 1; t[k++] = o + b; t[k++] = o + a;
	}
	for (j = 0; j + 2 < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = b;
			t[k++] = a; t[k++] = a + nu; t[k++] = b + nu;
		}
	m = qp_make(x, nn, t, NULL, k / 3);
	free(x);
	free(t);
	return m;
}

static void test_sphere(void)
{
	qaws_mesh_3d *s = qp_sphere(16, 8), *a, *b;
	qaws_mesh_report ra, rb;
	char msg[300];
	memset(&ra, 0, sizeof ra);
	memset(&rb, 0, sizeof rb);
	a = qp_run(s, QAWS_QUADS_PAIRING, 0, &ra);
	b = qp_run(s, QAWS_QUADS_PAIRING, 1, &rb);
	sprintf(msg, "a sphere's 224 triangles paired: %u quads, %u triangles left (the poles' 32 can pair at most to 16); pure: %u quads, %u others",
		ra.quality.quads, ra.quality.triangles, rb.quality.quads, rb.quality.triangles + rb.quality.polygons);
	printf("  %s\n", msg);
	TEST_ASSERT(a && b && ra.quality.quads == 96 + 16 && ra.quality.triangles == 0 && rb.quality.triangles + rb.quality.polygons == 0 && rb.quality.quads == 4 * ra.quality.quads, msg);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(s);
}

static qaws_mesh_3d* qp_box(int triangles)
{
	static unsigned int const q[24] = { 0, 2, 3, 1, 4, 5, 7, 6, 0, 1, 5, 4, 2, 6, 7, 3, 0, 4, 6, 2, 1, 3, 7, 5 };
	qaws_scalar x[24];
	unsigned int t[36], fs[6] = { 4, 4, 4, 4, 4, 4 }, i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)(i & 1 ? 1 : -1);
		x[3 * i + 1] = (qaws_scalar)(i & 2 ? 1 : -1);
		x[3 * i + 2] = (qaws_scalar)(i & 4 ? 1 : -1);
	}
	if (!triangles)
		return qp_make(x, 8, q, fs, 6);
	for (i = 0; i < 6; i++)
	{
		/* each square cut along its diagonal from its second corner */
		t[6 * i] = q[4 * i + 1]; t[6 * i + 1] = q[4 * i + 2]; t[6 * i + 2] = q[4 * i + 3];
		t[6 * i + 3] = q[4 * i + 1]; t[6 * i + 4] = q[4 * i + 3]; t[6 * i + 5] = q[4 * i];
	}
	return qp_make(x, 8, t, NULL, 12);
}

static void test_features_split(void)
{
	qaws_mesh_3d *box = qp_box(1), *cube = qp_box(0), *a, *b, *c;
	qaws_mesh_report ra, rb, rc;
	static qaws_scalar const ox[18] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1 };
	static unsigned int const ot[24] = { 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5 };
	qaws_mesh_3d* oct = qp_make(ox, 6, ot, NULL, 8);
	unsigned int i, kept = 0;
	char msg[400];
	memset(&ra, 0, sizeof ra);
	memset(&rb, 0, sizeof rb);
	memset(&rc, 0, sizeof rc);
	qaws_mesh_3d_detect_features(box, (qaws_scalar)0.5, 0, 0);
	qaws_mesh_3d_detect_features(cube, (qaws_scalar)0.5, 0, 0);
	a = qp_run(box, QAWS_QUADS_PAIRING, 0, &ra);
	b = qp_run(cube, QAWS_QUADS_SPLIT, 0, &rb);
	c = qp_run(oct, QAWS_QUADS_SPLIT, 0, &rc);
	for (i = 0; b && i < 8; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(b, i), q = qaws_mesh_3d_get_vertex(cube, i);
		kept += p.x == q.x && p.y == q.y && p.z == q.z;
	}
	sprintf(msg, "a box of triangles, its edges features, paired: %u quads, %u triangles, off plane %.1e; the cube split: %u quads, %u vertices, off plane %.1e, the 8 corners kept %u, %u feature edges; an octahedron split: %u quads, max deviation %.1e",
		ra.quality.quads, ra.quality.triangles, a ? qp_nonplanar(a) : -1, rb.quality.quads, b ? qaws_mesh_3d_get_vertex_count(b) : 0, b ? qp_nonplanar(b) : -1, kept,
		b ? qaws_mesh_3d_get_feature_edge_count(b) : 0, rc.quality.quads, (double)rc.max_deviation);
	printf("  %s\n", msg);
	TEST_ASSERT(a && ra.quality.quads == 6 && ra.quality.triangles == 0 && qp_nonplanar(a) < 1e-12 && b && rb.quality.quads == 24 && qaws_mesh_3d_get_vertex_count(b) == 26 &&
		qp_nonplanar(b) < 1e-9 && kept == 8 && qaws_mesh_3d_get_feature_edge_count(b) == 24 && c && rc.quality.quads == 24 && (double)rc.max_deviation < 0.3, msg);
	qaws_mesh_3d_destroy(a);
	qaws_mesh_3d_destroy(b);
	qaws_mesh_3d_destroy(c);
	qaws_mesh_3d_destroy(box);
	qaws_mesh_3d_destroy(cube);
	qaws_mesh_3d_destroy(oct);
}

int test_177_quads_pairing_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 177: quads from the input's own faces\n");
	test_grid();
	test_sphere();
	test_features_split();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
