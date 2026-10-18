/*
 * Test 128: hexahedra (GRID) and frame fields
 *
 *   - a box generated as a shell, in hexahedra: valid, all hexahedra, its
 *     volume exactly the box's, the scaled Jacobian bounded, the border
 *     closed and on the box's faces, the box's eight corners kept;
 *   - a NURBS torus's surface in hexahedra: valid, the border closed of
 *     genus 1, the volume near the torus's;
 *   - the frame field of the box's tetrahedra: no singular edge, every
 *     frame the world's axes; of the torus's: orthonormal frames, the
 *     border's cells holding an axis along its normal;
 *   - FRAME_FIELD: a turned box in perfect cubes along it, its volume exact;
 *     a bent bar; a round torus (singular edges) refused;
 *   - refusals: an open surface, the reserved methods, a field on hexahedra.
 */

#include "test_common.h"
#include "qaws_hexahedralize.h"
#include "qaws_tetrahedralize.h"
#include "qaws_mesh_generate.h"
#include "qaws_mesh_field.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static double const th_h[3] = { 1, 0.7, 0.5 };

static qaws_mesh_3d* th_box(qaws_shell_atlas** out_atlas, qaws_surface** f)
{
	static int const q[6][4][3] = {
		{ { -1, -1, -1 }, { -1, 1, -1 }, { 1, -1, -1 }, { 1, 1, -1 } }, { { -1, -1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { 1, -1, -1 }, { -1, -1, 1 }, { 1, -1, 1 } }, { { -1, 1, -1 }, { -1, 1, 1 }, { 1, 1, -1 }, { 1, 1, 1 } },
		{ { -1, -1, -1 }, { -1, -1, 1 }, { -1, 1, -1 }, { -1, 1, 1 } }, { { 1, -1, -1 }, { 1, 1, -1 }, { 1, -1, 1 }, { 1, 1, 1 } } };
	qaws_generate_desc g;
	qaws_mesh_3d* s = NULL;
	unsigned int i, k;
	for (i = 0; i < 6; i++)
	{
		qaws_surface_bilinear_desc bd;
		qaws_vec3* p[4];
		p[0] = &bd.p00; p[1] = &bd.p10; p[2] = &bd.p01; p[3] = &bd.p11;
		for (k = 0; k < 4; k++)
		{
			p[k]->x = (qaws_scalar)(th_h[0] * q[i][k][0]);
			p[k]->y = (qaws_scalar)(th_h[1] * q[i][k][1]);
			p[k]->z = (qaws_scalar)(th_h[2] * q[i][k][2]);
		}
		f[i] = NULL;
		qaws_surface_create_bilinear(&bd, &f[i]);
	}
	*out_atlas = NULL;
	if (qaws_shell_atlas_create((qaws_surface const* const*)f, 6, 0, out_atlas) != QAWS_STATUS_OK)
		return NULL;
	memset(&g, 0, sizeof g);
	g.target_length = (qaws_scalar)0.1;
	qaws_mesh_generate_shell(*out_atlas, &g, &s, NULL);
	return s;
}

static qaws_surface* th_torus(double R, double r)
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

static void test_box(void)
{
	qaws_surface* f[6];
	qaws_shell_atlas* atlas = NULL;
	qaws_mesh_3d* s = th_box(&atlas, f);
	qaws_hexahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_status st;
	char msg[400];
	unsigned int i;
	memset(&d, 0, sizeof d);
	d.cell_size = (qaws_scalar)0.15;
	st = s ? qaws_volume_mesh_hexahedralize(s, &d, &v, NULL) : QAWS_STATUS_INVALID_ARGUMENT;
	if (st != QAWS_STATUS_OK || !v)
		TEST_ASSERT(0, "box: hexahedralized");
	else
	{
		qaws_volume_quality q;
		qaws_mesh_3d* b = NULL;
		qaws_mesh_topology top;
		unsigned int nb, off = 0, corners = 0, a;
		qaws_scalar* x;
		qaws_volume_mesh_measure_quality(v, &q);
		qaws_volume_mesh_get_boundary(v, &b);
		qaws_mesh_3d_get_topology(b, &top);
		nb = qaws_mesh_3d_get_vertex_count(b);
		x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nb);
		qaws_mesh_3d_get_positions(b, x);
		for (i = 0; i < nb; i++)
		{
			unsigned int on = 0, k;
			for (k = 0; k < 3; k++)
				on += fabs(fabs((double)x[3 * i + k]) - th_h[k]) < 1e-9;
			off += on == 0;
		}
		for (a = 0; a < 8; a++)
			for (i = 0; i < nb; i++)
				if (fabs((double)x[3 * i] - (a & 1 ? 1 : -1) * th_h[0]) + fabs((double)x[3 * i + 1] - (a & 2 ? 1 : -1) * th_h[1])
					+ fabs((double)x[3 * i + 2] - (a & 4 ? 1 : -1) * th_h[2]) < 1e-9)
				{
					corners++;
					break;
				}
		sprintf(msg, "box: %u hexahedra (of %u cells), validate %s, volume %.12f, scaled Jacobian from %.3f, border closed %d genus %u, %u border vertices off the faces, %u corners of 8",
			q.hexahedra, qaws_volume_mesh_get_cell_count(v), qaws_status_to_string(qaws_volume_mesh_validate(v)), (double)q.volume, (double)q.min_scaled_jacobian,
			top.closed, top.genus, off, corners);
		TEST_ASSERT(qaws_volume_mesh_validate(v) == QAWS_STATUS_OK && q.hexahedra == qaws_volume_mesh_get_cell_count(v) && q.hexahedra > 300
			&& fabs((double)q.volume - 2.8) < 1e-9 && (double)q.min_scaled_jacobian > 0.3 && top.closed && top.genus == 0 && off == 0 && corners == 8, msg);
		free(x);
		qaws_mesh_3d_destroy(b);
	}
	/* its tetrahedra's frame field: the world's axes, no singular edge */
	{
		qaws_tetrahedralize_desc td;
		qaws_volume_mesh* t = NULL;
		memset(&td, 0, sizeof td);
		if (s && qaws_volume_mesh_tetrahedralize(s, &td, &t, NULL) == QAWS_STATUS_OK)
		{
			unsigned int nc = qaws_volume_mesh_get_cell_count(t), ns = 0;
			qaws_scalar* fr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * nc);
			unsigned char* sing = (unsigned char*)malloc(nc);
			double off = 0;
			st = qaws_volume_mesh_frame_field(t, NULL, fr, sing);
			for (i = 0; i < nc; i++)
			{
				unsigned int k;
				ns += sing[i] != 0;
				for (k = 0; k < 3; k++)
					off = fmax(off, 1 - fmax(fabs((double)fr[9 * i + 3 * k]), fmax(fabs((double)fr[9 * i + 3 * k + 1]), fabs((double)fr[9 * i + 3 * k + 2]))));
			}
			sprintf(msg, "box frame field: %s, %u cells by singular edges, axes off the world's by %.2e", qaws_status_to_string(st), ns, off);
			TEST_ASSERT(st == QAWS_STATUS_OK && ns == 0 && off < 1e-9, msg);
			free(fr);
			free(sing);
		}
		else
			TEST_ASSERT(0, "box: tetrahedralized for its field");
		qaws_volume_mesh_destroy(t);
	}
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(s);
	qaws_shell_atlas_destroy(atlas);
	for (i = 0; i < 6; i++)
		qaws_surface_destroy(f[i]);
}

static void test_torus(void)
{
	qaws_surface* t = th_torus(2, 0.7);
	qaws_generate_desc g;
	qaws_mesh_3d* s = NULL;
	qaws_hexahedralize_desc d;
	qaws_volume_mesh* v = NULL;
	qaws_status st;
	char msg[400];
	memset(&g, 0, sizeof g);
	g.target_length = (qaws_scalar)0.2;
	st = qaws_mesh_generate_surface(t, &g, &s, NULL);
	memset(&d, 0, sizeof d);
	d.cell_size = (qaws_scalar)0.22;
	if (st == QAWS_STATUS_OK)
		st = qaws_volume_mesh_hexahedralize(s, &d, &v, NULL);
	if (st != QAWS_STATUS_OK || !v)
		TEST_ASSERT(0, "torus: hexahedralized");
	else
	{
		qaws_volume_quality q;
		qaws_mesh_3d* b = NULL;
		qaws_mesh_topology top;
		double exact = 2 * M_PI * M_PI * 2 * 0.7 * 0.7;
		qaws_volume_mesh_measure_quality(v, &q);
		qaws_volume_mesh_get_boundary(v, &b);
		qaws_mesh_3d_get_topology(b, &top);
		sprintf(msg, "torus: %u hexahedra, validate %s, volume %.4f for %.4f, scaled Jacobian from %.3f, border closed %d genus %u", q.hexahedra,
			qaws_status_to_string(qaws_volume_mesh_validate(v)), (double)q.volume, exact, (double)q.min_scaled_jacobian, top.closed, top.genus);
		TEST_ASSERT(qaws_volume_mesh_validate(v) == QAWS_STATUS_OK && fabs((double)q.volume - exact) < 0.05 * exact && (double)q.min_scaled_jacobian > 0.2
			&& top.closed && top.genus == 1, msg);
		qaws_mesh_3d_destroy(b);
	}
	/* its tetrahedra's field: orthonormal, the border's normals held */
	{
		qaws_tetrahedralize_desc td;
		qaws_volume_mesh* tt = NULL;
		memset(&td, 0, sizeof td);
		if (s && qaws_volume_mesh_tetrahedralize(s, &td, &tt, NULL) == QAWS_STATUS_OK)
		{
			unsigned int nc = qaws_volume_mesh_get_cell_count(tt), nv = qaws_volume_mesh_get_vertex_count(tt), i, k, j;
			qaws_scalar *fr = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 9 * nc), *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
			unsigned int *cv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc), *adj = (unsigned int*)malloc(sizeof(unsigned int) * 4 * nc);
			double orth = 0, along = 0;
			unsigned int nborder = 0;
			static unsigned char const F[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
			st = qaws_volume_mesh_frame_field(tt, NULL, fr, NULL);
			qaws_volume_mesh_get_positions(tt, x);
			qaws_volume_mesh_get_cells(tt, cv, NULL);
			qaws_volume_mesh_get_cell_adjacency(tt, adj);
			for (i = 0; i < nc; i++)
			{
				qaws_scalar const* f = fr + 9 * i;
				for (k = 0; k < 3; k++)
					for (j = 0; j < 3; j++)
						orth = fmax(orth, fabs((double)(f[3 * k] * f[3 * j] + f[3 * k + 1] * f[3 * j + 1] + f[3 * k + 2] * f[3 * j + 2]) - (k == j)));
				for (k = 0; k < 4; k++)
				{
					qaws_scalar const *a = x + 3 * cv[4 * i + F[k][0]], *b = x + 3 * cv[4 * i + F[k][1]], *c = x + 3 * cv[4 * i + F[k][2]];
					double u[3], w[3], n[3], l, best = 0;
					if (adj[4 * i + k] != QAWS_MESH_NONE)
						continue;
					for (j = 0; j < 3; j++)
					{
						u[j] = (double)(b[j] - a[j]);
						w[j] = (double)(c[j] - a[j]);
					}
					n[0] = u[1] * w[2] - u[2] * w[1];
					n[1] = u[2] * w[0] - u[0] * w[2];
					n[2] = u[0] * w[1] - u[1] * w[0];
					l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
					for (j = 0; j < 3; j++)
						best = fmax(best, fabs((double)(f[3 * j] * n[0] + f[3 * j + 1] * n[1] + f[3 * j + 2] * n[2])) / l);
					along = fmax(along, 1 - best);
					nborder++;
					break;   /* (the first border facet holds an axis) */
				}
			}
			sprintf(msg, "torus frame field: %s, orthonormal to %.1e, %u border cells, an axis along the normal to %.1e", qaws_status_to_string(st), orth, nborder, along);
			TEST_ASSERT(st == QAWS_STATUS_OK && orth < 1e-9 && nborder > 100 && along < 1e-9, msg);
			free(fr); free(x); free(cv); free(adj);
		}
		else
			TEST_ASSERT(0, "torus: tetrahedralized for its field");
		qaws_volume_mesh_destroy(tt);
	}
	qaws_volume_mesh_destroy(v);
	qaws_mesh_3d_destroy(s);
	qaws_surface_destroy(t);
}

/* the border of the lattice block n[0] x n[1] x n[2] as quads, mapped by f */
static qaws_mesh_3d* th_lattice(unsigned int const* n, void (*f)(double const* p, double* x))
{
	static unsigned char const F[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
	static int const FN[6][3] = { { 0, 0, -1 }, { 0, 0, 1 }, { 0, -1, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { -1, 0, 0 } };
	static int const CO[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
	unsigned int nv = (n[0] + 1) * (n[1] + 1) * (n[2] + 1), *map = (unsigned int*)malloc(sizeof(unsigned int) * nv), i, j, k, q, used = 0, nf = 0;
	unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * 4 * 2 * (n[0] * n[1] + n[1] * n[2] + n[0] * n[2]));
	unsigned int* fs = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (n[0] * n[1] + n[1] * n[2] + n[0] * n[2]));
	qaws_scalar* x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < nv; i++)
		map[i] = QAWS_MESH_NONE;
	for (k = 0; k < n[2]; k++)
		for (j = 0; j < n[1]; j++)
			for (i = 0; i < n[0]; i++)
				for (q = 0; q < 6; q++)
				{
					long a = (long)i + FN[q][0], b = (long)j + FN[q][1], c = (long)k + FN[q][2];
					unsigned int e;
					if (a >= 0 && b >= 0 && c >= 0 && a < (long)n[0] && b < (long)n[1] && c < (long)n[2])
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
	qaws_mesh_3d_create(&d, &m);
	free(map); free(fv); free(fs); free(x);
	return m;
}

static void th_turned_box(double const* p, double* x)
{
	double a = p[0] * 2 - 1, b = (p[1] * 2 - 1) * 0.7, c = (p[2] * 2 - 1) * 0.5, y0 = a * cos(0.5) - b * sin(0.5), y1 = a * sin(0.5) + b * cos(0.5);
	x[0] = y0;
	x[1] = y1 * cos(0.3) - c * sin(0.3);
	x[2] = y1 * sin(0.3) + c * cos(0.3);
}

static void th_bent_bar(double const* p, double* x)
{
	double r = 1 + p[0], t = M_PI / 2 * p[1];
	x[0] = r * cos(t);
	x[1] = r * sin(t);
	x[2] = p[2];
}

static void test_frame_field(void)
{
	static unsigned int const nbox[3] = { 16, 12, 8 }, nbent[3] = { 8, 20, 8 };
	qaws_mesh_3d *box = th_lattice(nbox, th_turned_box), *bent = th_lattice(nbent, th_bent_bar), *ts = NULL;
	qaws_surface* t = th_torus(2, 0.7);
	qaws_generate_desc g;
	qaws_hexahedralize_desc d;
	qaws_volume_mesh *v1 = NULL, *v2 = NULL, *v3 = NULL;
	qaws_volume_quality q1, q2;
	qaws_status s1, s2, s3;
	char msg[400];
	memset(&d, 0, sizeof d);
	d.method = QAWS_HEXES_FRAME_FIELD;
	d.cell_size = (qaws_scalar)0.15;
	s1 = qaws_volume_mesh_hexahedralize(box, &d, &v1, NULL);
	s2 = qaws_volume_mesh_hexahedralize(bent, &d, &v2, NULL);
	memset(&g, 0, sizeof g);
	g.target_length = (qaws_scalar)0.25;
	qaws_mesh_generate_surface(t, &g, &ts, NULL);
	d.cell_size = (qaws_scalar)0.25;
	s3 = qaws_volume_mesh_hexahedralize(ts, &d, &v3, NULL);
	memset(&q1, 0, sizeof q1);
	memset(&q2, 0, sizeof q2);
	if (v1)
		qaws_volume_mesh_measure_quality(v1, &q1);
	if (v2)
		qaws_volume_mesh_measure_quality(v2, &q2);
	sprintf(msg, "turned box along its field: %s, %u hexahedra, validate %s, scaled Jacobian from %.4f, volume %.12f, edges %.3f..%.3f",
		qaws_status_to_string(s1), q1.hexahedra, v1 ? qaws_status_to_string(qaws_volume_mesh_validate(v1)) : "-", (double)q1.min_scaled_jacobian,
		(double)q1.volume, (double)q1.min_edge, (double)q1.max_edge);
	TEST_ASSERT(s1 == QAWS_STATUS_OK && v1 && qaws_volume_mesh_validate(v1) == QAWS_STATUS_OK && q1.hexahedra == qaws_volume_mesh_get_cell_count(v1)
		&& q1.hexahedra > 500 && (double)q1.min_scaled_jacobian > 0.99 && fabs((double)q1.volume - 2.8) < 1e-9 && (double)q1.min_edge > 0.12
		&& (double)q1.max_edge < 0.18, msg);
	sprintf(msg, "bent bar along its field: %s, %u hexahedra, validate %s, scaled Jacobian from %.4f, volume %.5f for 3 pi / 4", qaws_status_to_string(s2),
		q2.hexahedra, v2 ? qaws_status_to_string(qaws_volume_mesh_validate(v2)) : "-", (double)q2.min_scaled_jacobian, (double)q2.volume);
	TEST_ASSERT(s2 == QAWS_STATUS_OK && v2 && qaws_volume_mesh_validate(v2) == QAWS_STATUS_OK && (double)q2.min_scaled_jacobian > 0.9
		&& fabs((double)q2.volume - 0.75 * M_PI) < 0.01 * 0.75 * M_PI, msg);
	sprintf(msg, "round torus along its field (singular edges): %s", qaws_status_to_string(s3));
	TEST_ASSERT(s3 == QAWS_STATUS_UNSUPPORTED_OPERATION && !v3, msg);
	qaws_volume_mesh_destroy(v1);
	qaws_volume_mesh_destroy(v2);
	qaws_volume_mesh_destroy(v3);
	qaws_mesh_3d_destroy(box);
	qaws_mesh_3d_destroy(bent);
	qaws_mesh_3d_destroy(ts);
	qaws_surface_destroy(t);
}

static void test_refusals(void)
{
	static qaws_scalar const p[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 }, cube[24] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1 };
	static unsigned int const f[4] = { 0, 1, 2, 3 }, hex[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	unsigned int fs = 4;
	qaws_mesh_3d_desc md;
	qaws_mesh_3d* open = NULL;
	qaws_hexahedralize_desc d;
	qaws_volume_mesh_desc vd;
	qaws_volume_mesh *v = NULL, *h = NULL;
	qaws_cell_type ht = QAWS_CELL_HEXAHEDRON;
	qaws_scalar fr[9];
	qaws_status s1, s2, s3;
	char msg[200];
	memset(&md, 0, sizeof md);
	md.positions = p;
	md.vertex_count = 4;
	md.face_vertices = f;
	md.face_sizes = &fs;
	md.face_count = 1;
	qaws_mesh_3d_create(&md, &open);
	memset(&d, 0, sizeof d);
	s1 = qaws_volume_mesh_hexahedralize(open, &d, &v, NULL);
	d.method = QAWS_HEXES_POLYCUBE;
	s2 = qaws_volume_mesh_hexahedralize(open, &d, &v, NULL);
	memset(&vd, 0, sizeof vd);
	vd.positions = cube;
	vd.vertex_count = 8;
	vd.cell_vertices = hex;
	vd.cell_types = &ht;
	vd.cell_count = 1;
	qaws_volume_mesh_create(&vd, &h);
	s3 = qaws_volume_mesh_frame_field(h, NULL, fr, NULL);
	sprintf(msg, "refusals: an open surface %s, POLYCUBE %s, a field on hexahedra %s", qaws_status_to_string(s1), qaws_status_to_string(s2),
		qaws_status_to_string(s3));
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_TOPOLOGY && s2 == QAWS_STATUS_UNSUPPORTED_OPERATION && s3 == QAWS_STATUS_UNSUPPORTED_OPERATION && !v, msg);
	qaws_mesh_3d_destroy(open);
	qaws_volume_mesh_destroy(h);
}

int test_128_hexahedralize_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 128: hexahedra and frame fields\n");
	test_box();
	test_torus();
	test_frame_field();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
