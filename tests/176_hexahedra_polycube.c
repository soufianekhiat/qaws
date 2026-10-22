/*
 * Test 176: hexahedra of a polycube map (QAWS_HEXES_POLYCUBE), small meshes
 *
 *   - a box turned in space: its polycube is itself, the hexahedra a grid
 *     along it (count, near cubes; the volume to 1e-4: the border's vertices
 *     lie on its faces, the hexahedra's faces inside a little off plane);
 *   - an L of three blocks: its two arms' faces on their own planes, the
 *     hexahedra a grid following the L;
 *   - a bent bar: its inner and outer faces labelled by their nearest axes,
 *     so its polycube turns a corner at 45 degrees (a polycube's limit: the
 *     hexahedra there are poorer), the grid following it;
 *   - a sphere: deformed into a cube (its field's singular edges are why
 *     FRAME_FIELD refuses it), filled, valid, its volume near the ball's.
 */

#include "test_common.h"
#include "qaws_hexahedralize.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define PC_PI 3.14159265358979323846

/* the border of the lattice cells n[0] x n[1] x n[2] that keep() keeps, as quads, mapped by f from the unit cube */
static qaws_mesh_3d* pc_block(unsigned int const* n, int (*keep)(unsigned int, unsigned int, unsigned int), void (*f)(double const* p, double* x))
{
	static unsigned char const F[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	static int const FN[6][3] = { { 0, 0, -1 }, { 0, 0, 1 }, { 0, -1, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { -1, 0, 0 } };
	static int const CO[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
	unsigned int nv = (n[0] + 1) * (n[1] + 1) * (n[2] + 1), *map = (unsigned int*)malloc(sizeof(unsigned int) * nv), i, j, k, q, used = 0, nf = 0;
	unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * 6 * n[0] * n[1] * n[2]);
	unsigned int* fs = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n[0] * n[1] * n[2]);
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < nv; i++) map[i] = QAWS_MESH_NONE;
	for (k = 0; k < n[2]; k++)
		for (j = 0; j < n[1]; j++)
			for (i = 0; i < n[0]; i++)
				for (q = 0; q < 6 && keep(i, j, k); q++)
				{
					long a = (long)i + FN[q][0], b = (long)j + FN[q][1], c = (long)k + FN[q][2];
					unsigned int e;
					if (a >= 0 && b >= 0 && c >= 0 && a < (long)n[0] && b < (long)n[1] && c < (long)n[2] && keep((unsigned int)a, (unsigned int)b, (unsigned int)c))
						continue;
					for (e = 0; e < 4; e++)
					{
						unsigned int ci = i + CO[F[q][e]][0], cj = j + CO[F[q][e]][1], ck = k + CO[F[q][e]][2];
						unsigned int v = ci + (n[0] + 1) * (cj + (n[1] + 1) * ck);
						if (map[v] == QAWS_MESH_NONE)
						{
							double p[3], y[3];
							p[0] = (double)ci / n[0];
							p[1] = (double)cj / n[1];
							p[2] = (double)ck / n[2];
							f(p, y);
							x[3 * used] = (qaws_scalar)y[0];
							x[3 * used + 1] = (qaws_scalar)y[1];
							x[3 * used + 2] = (qaws_scalar)y[2];
							map[v] = used++;
						}
						fv[4 * nf + e] = map[v];
					}
					fs[nf++] = 4;
				}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = used;
	d.face_vertices = fv;
	d.face_sizes = fs;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(map); free(fv); free(fs); free(x);
	return m;
}

static int pc_all(unsigned int i, unsigned int j, unsigned int k)
{
	(void)i; (void)j; (void)k;
	return 1;
}

/* the L: of 2 x 2 x 1 blocks, all but the corner (1, 1) */
static int pc_l(unsigned int i, unsigned int j, unsigned int k)
{
	(void)k;
	return !(i >= 4 && j >= 4);
}

/* 2 x 1 x 1, turned by 0.5 about z and 0.3 about x */
static void pc_turned(double const* p, double* x)
{
	double a = p[0] * 2 - 1, b = p[1] - 0.5, c = p[2] - 0.5, y0 = a * cos(0.5) - b * sin(0.5), y1 = a * sin(0.5) + b * cos(0.5);
	x[0] = y0;
	x[1] = y1 * cos(0.3) - c * sin(0.3);
	x[2] = y1 * sin(0.3) + c * cos(0.3);
}

/* the L's 2 x 2 x 0.5 frame, turned by 0.4 about z */
static void pc_lmap(double const* p, double* x)
{
	double a = 2 * p[0], b = 2 * p[1];
	x[0] = a * cos(0.4) - b * sin(0.4);
	x[1] = a * sin(0.4) + b * cos(0.4);
	x[2] = 0.5 * p[2];
}

static void pc_bent(double const* p, double* x)
{
	double r = 1 + p[0], t = PC_PI / 2 * p[1];
	x[0] = r * cos(t);
	x[1] = r * sin(t);
	x[2] = p[2];
}

/* a unit sphere: rings between single poles */
static qaws_mesh_3d* pc_sphere(unsigned int nu, unsigned int nv)
{
	unsigned int n = nu * (nv - 1) + 2, i, j, k = 0;
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * (nv - 1));
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (j = 1; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double th = PC_PI * j / nv, ph = 2 * PC_PI * i / nu;
			qaws_scalar* p = x + 3 * ((j - 1) * nu + i);
			p[0] = (qaws_scalar)(sin(th) * cos(ph)); p[1] = (qaws_scalar)(sin(th) * sin(ph)); p[2] = (qaws_scalar)cos(th);
		}
	x[3 * (n - 2)] = 0; x[3 * (n - 2) + 1] = 0; x[3 * (n - 2) + 2] = 1;
	x[3 * (n - 1)] = 0; x[3 * (n - 1) + 1] = 0; x[3 * (n - 1) + 2] = -1;
	for (i = 0; i < nu; i++)
	{
		unsigned int a = i, b = (i + 1) % nu, o = (nv - 2) * nu;
		t[k++] = n - 2; t[k++] = a; t[k++] = b;
		t[k++] = n - 1; t[k++] = o + b; t[k++] = o + a;
	}
	for (j = 0; j + 2 < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu;
			t[k++] = a; t[k++] = b + nu; t[k++] = b;
			t[k++] = a; t[k++] = a + nu; t[k++] = b + nu;
		}
	memset(&d, 0, sizeof d);
	d.positions = x;
	d.vertex_count = n;
	d.face_vertices = t;
	d.face_count = k / 3;
	qaws_mesh_3d_create(&m, &d);
	free(x);
	free(t);
	return m;
}

typedef struct pc_out
{
	qaws_status st;
	qaws_volume_quality q;
	int valid;
} pc_out;

static pc_out pc_run(qaws_mesh_3d const* s, qaws_hexahedralize_method method, double h)
{
	qaws_hexahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	pc_out o;
	memset(&o, 0, sizeof o);
	memset(&d, 0, sizeof d);
	d.method = method;
	d.cell_size = (qaws_scalar)h;
	o.st = qaws_volume_mesh_hexahedralize(&v, NULL, s, &d);
	if (v)
	{
		qaws_volume_mesh_measure_quality(&o.q, v);
		o.valid = qaws_volume_mesh_validate(v) == QAWS_STATUS_OK;
		qaws_volume_mesh_destroy(v);
	}
	return o;
}

static void test_boxes(void)
{
	static unsigned int const nb[3] = { 4, 2, 2 }, nl[3] = { 8, 8, 2 }, nbent[3] = { 2, 6, 2 };
	qaws_mesh_3d *box = pc_block(nb, pc_all, pc_turned), *l = pc_block(nl, pc_l, pc_lmap), *bent = pc_block(nbent, pc_all, pc_bent);
	pc_out a = pc_run(box, QAWS_HEXES_POLYCUBE, 0.25), b = pc_run(l, QAWS_HEXES_POLYCUBE, 0.25), c = pc_run(bent, QAWS_HEXES_POLYCUBE, 0.25);
	char msg[400];
	sprintf(msg, "a turned box: %u hexahedra (128), scaled Jacobian from %.4f, volume %.10f (2), valid %d", a.q.hexahedra, (double)a.q.min_scaled_jacobian,
		(double)a.q.volume, a.valid);
	printf("  %s\n", msg);
	TEST_ASSERT(a.st == QAWS_STATUS_OK && a.q.hexahedra == 128 && a.q.min_scaled_jacobian > 0.98 && fabs((double)a.q.volume - 2) < 2e-4 && a.valid, msg);
	sprintf(msg, "an L of three blocks: %u hexahedra (96), scaled Jacobian from %.4f, volume %.10f (1.5), valid %d", b.q.hexahedra, (double)b.q.min_scaled_jacobian,
		(double)b.q.volume, b.valid);
	printf("  %s\n", msg);
	TEST_ASSERT(b.st == QAWS_STATUS_OK && b.q.hexahedra == 96 && b.q.min_scaled_jacobian > 0.95 && fabs((double)b.q.volume - 1.5) < 3e-4 && b.valid, msg);
	sprintf(msg, "a bent bar, its polycube cornered: %u hexahedra, scaled Jacobian from %.4f, volume %.5f (3 pi / 4 = %.5f), valid %d", c.q.hexahedra,
		(double)c.q.min_scaled_jacobian, (double)c.q.volume, 0.75 * PC_PI, c.valid);
	printf("  %s\n", msg);
	TEST_ASSERT(c.st == QAWS_STATUS_OK && c.q.min_scaled_jacobian > 0.2 && fabs((double)c.q.volume - 0.75 * PC_PI) < 0.05 * 0.75 * PC_PI && c.valid, msg);
	qaws_mesh_3d_destroy(box);
	qaws_mesh_3d_destroy(l);
	qaws_mesh_3d_destroy(bent);
}

static void test_sphere(void)
{
	qaws_mesh_3d* s = pc_sphere(16, 8);
	pc_out a = pc_run(s, QAWS_HEXES_POLYCUBE, 0.3), f = pc_run(s, QAWS_HEXES_FRAME_FIELD, 0.3);
	char msg[300];
	sprintf(msg, "a sphere into a cube: %s, %u hexahedra, scaled Jacobian from %.3f, volume %.4f (the ball 4.1888), valid %d; FRAME_FIELD %s",
		qaws_status_to_string(a.st), a.q.hexahedra, (double)a.q.min_scaled_jacobian, (double)a.q.volume, a.valid, qaws_status_to_string(f.st));
	printf("  %s\n", msg);
	TEST_ASSERT(a.st == QAWS_STATUS_OK && a.q.hexahedra > 50 && a.q.min_scaled_jacobian > 0 && fabs((double)a.q.volume - 4.18879) < 0.15 * 4.18879 && a.valid &&
		f.st == QAWS_STATUS_UNSUPPORTED_OPERATION, msg);
	qaws_mesh_3d_destroy(s);
}

int test_176_hexahedra_polycube_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 176: hexahedra of a polycube map\n");
	test_boxes();
	test_sphere();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
