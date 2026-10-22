/*
 * Test 136: quads by integer-grid maps (GLOBAL)
 *
 *   - the map itself on the can (heavy): every triangle positive in uv,
 *     the lattice's isolines the quads' edges;
 *   - a torus (a grid: 0 irregular), a creased box (each face a grid, its
 *     8 corners of valence 3, all regular by the feature-aware measure),
 *     the showcase's blob (irregular as many as its cones, 8), the can
 *     with its rims creased (irregular no more than INSTANT's), a wavy open
 *     patch (its border on isolines: the border vertices of valence 3, its
 *     corners 2);
 *   - tori of radii 1 and 0.7 at 0.22, 1 and 0.5 at 0.3: their periods
 *     quantized together, every straight edge loop once round the hole or
 *     once round the tube (the flow's own periods made them spirals);
 *   - every result valid, all quads, footpoints on the reference, scaled
 *     Jacobian above 0, made by the map (report iterations 1); against
 *     INSTANT (pure) on the same inputs: irregular and AD no higher. By
 *     default the box, the can and the patch, against INSTANT's numbers
 *     measured below; the torus, the blob and INSTANT run alongside, heavy.
 *
 * Measured (GLOBAL against INSTANT, pure, both finished alike): torus at
 * 0.3 640 quads both, 0 irregular, SJ min 0.994 (0.985), AD 3.13 (3.60);
 * box at 0.1 376 quads both, SJ 1, AD 0.00 (0.04); blob at 0.11 1397
 * quads (INSTANT 1276, from twice the target split), 8 irregular both, SJ
 * min 0.815 (0.772), AD 2.74 (3.27), ED 0.139 (0.175); can at 0.14 766
 * quads both, 8 irregular, SJ min 0.853 (0.801), AD 2.23 (2.23), ED 0.129
 * (0.138); wavy patch at 0.1 498 quads both, 2 irregular, AD 2.97 (3.22).
 *   - a sweep of the five shapes at 0.75, 1.25, 1.6 and 2.2 times these
 *     targets (three of them by default, all heavy): every one by the map,
 *     valid, all quads, SJ above 0.
 *
 * Over the sweep the map made 14 of 20 with the flow's integers alone (the
 * rest INSTANT's: folds left on the can at 1.25 and 1.6, the box at 2.2
 * from corner steps read in the wrong plane, a cone wound to 90 degrees on
 * the patch at 0.75 and 1.25, a transition off at a box corner at 0.75);
 * now 20 of 20, the can and the patch by the quantized integers: can at
 * 1.25 512 quads, SJ min 0.860 (INSTANT 0.178), AD 2.41 (2.78); patch at
 * 1.25 2 irregular (3), AD 3.59 (4.08); box at 2.2 AD 0.00.
 */

#include "test_common.h"
#include "qaws_quadrangulate.h"
#include "qaws_mesh_generate.h"
#include "qaws_mesh_measure.h"
#include "mesh/qaws_quad_global.h"
#include <time.h>

typedef struct tg_res
{
	qaws_shape_quality s;
	unsigned int quads, faces, border3, border, made;
	double foot, ms;
	int valid;
} tg_res;

static qaws_mesh_3d* tg_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
	unsigned int i;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(p);
	return m;
}

/* the showcase's blob: an icosahedron split once (edges split, then new-old edges flipped), pushed in and out smoothly */
static qaws_mesh_3d* tg_blob(void)
{
	static unsigned int const f[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
		3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double g = (1 + sqrt(5.0)) / 2, x[36];
	double const v[36] = { -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, 0, 0, -1, g, 0, 1, g, 0, -1, -g, 0, 1, -g, g, 0, -1, g, 0, 1, -g, 0, -1, -g, 0, 1 };
	qaws_mesh_3d* m;
	unsigned int i, n, ne, nv0, *e, *src, w;
	qaws_scalar* p;
	for (i = 0; i < 36; i++)
		x[i] = v[i] / sqrt(1 + g * g);
	m = tg_make(x, 12, f, 20);
	ne = qaws_mesh_3d_get_edge_count(m);
	nv0 = qaws_mesh_3d_get_vertex_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	src = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
		if (qaws_mesh_3d_split_edge(&w, m, NULL, NULL, e[2 * i], e[2 * i + 1]) == QAWS_STATUS_OK)
		{
			src[2 * (w - nv0)] = e[2 * i];
			src[2 * (w - nv0) + 1] = e[2 * i + 1];
		}
	free(e);
	ne = qaws_mesh_3d_get_edge_count(m);
	e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
	{
		unsigned int a = e[2 * i], b = e[2 * i + 1], nw, old;
		if ((a >= nv0) == (b >= nv0))
			continue;
		nw = a >= nv0 ? a : b;
		old = a >= nv0 ? b : a;
		if (old != src[2 * (nw - nv0)] && old != src[2 * (nw - nv0) + 1] && qaws_mesh_3d_can_flip_edge(m, a, b))
			qaws_mesh_3d_flip_edge(m, a, b);
	}
	free(e);
	free(src);
	n = qaws_mesh_3d_get_vertex_count(m);
	p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * n);
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < n; i++)
	{
		double q[3], l, s;
		q[0] = (double)p[3 * i]; q[1] = (double)p[3 * i + 1]; q[2] = (double)p[3 * i + 2];
		l = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		q[0] /= l; q[1] /= l; q[2] /= l;
		s = 1 + 0.32 * sin(2.2 * q[0] + 0.4) * cos(1.7 * q[1]) + 0.2 * q[2] * q[2];
		p[3 * i] = (qaws_scalar)(q[0] * 1.25 * s); p[3 * i + 1] = (qaws_scalar)(q[1] * s); p[3 * i + 2] = (qaws_scalar)(q[2] * 0.9 * s);
	}
	qaws_mesh_3d_set_positions(m, p);
	free(p);
	return m;
}

/* the showcase's can: 8 sides of radius 1, height 1.6, each cap a fan; its rims creases */
static qaws_mesh_3d* tg_can(void)
{
	unsigned int n = 8, i, k = 0, fv[3 * 32];
	double x[3 * 18];
	qaws_mesh_3d* m;
	for (i = 0; i < n; i++)
	{
		double a = 2 * M_PI * i / n;
		x[3 * i] = cos(a); x[3 * i + 1] = sin(a); x[3 * i + 2] = -0.8;
		x[3 * (n + i)] = cos(a); x[3 * (n + i) + 1] = sin(a); x[3 * (n + i) + 2] = 0.8;
	}
	x[6 * n] = 0; x[6 * n + 1] = 0; x[6 * n + 2] = -0.8;
	x[6 * n + 3] = 0; x[6 * n + 4] = 0; x[6 * n + 5] = 0.8;
	for (i = 0; i < n; i++)
	{
		unsigned int a = i, b = (i + 1) % n;
		fv[k++] = a; fv[k++] = b; fv[k++] = n + b;
		fv[k++] = a; fv[k++] = n + b; fv[k++] = n + a;
		fv[k++] = 2 * n; fv[k++] = b; fv[k++] = a;
		fv[k++] = 2 * n + 1; fv[k++] = n + a; fv[k++] = n + b;
	}
	m = tg_make(x, 2 * n + 2, fv, 4 * n);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(50 * M_PI / 180), 0, 0);
	return m;
}

/* test 125's torus: nu x nv quads of radii R and r split in triangles, every other column shifted half a step */
static qaws_mesh_3d* tg_torus(unsigned int nu, unsigned int nv, double R, double r)
{
	double* x = (double*)malloc(sizeof(double) * 3 * nu * nv);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * nu * nv), i, j, k = 0;
	qaws_mesh_3d* m;
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			double u = 2 * M_PI * i / nu, v = 2 * M_PI * (j + 0.5 * (i & 1)) / nv;
			x[3 * (j * nu + i)] = (R + r * cos(v)) * cos(u);
			x[3 * (j * nu + i) + 1] = (R + r * cos(v)) * sin(u);
			x[3 * (j * nu + i) + 2] = r * sin(v);
		}
	for (j = 0; j < nv; j++)
		for (i = 0; i < nu; i++)
		{
			unsigned int a = j * nu + i, b = j * nu + (i + 1) % nu, c = ((j + 1) % nv) * nu + (i + 1) % nu, d = ((j + 1) % nv) * nu + i;
			f[k++] = a; f[k++] = b; f[k++] = c;
			f[k++] = a; f[k++] = c; f[k++] = d;
		}
	m = tg_make(x, nu * nv, f, 2 * nu * nv);
	free(x);
	free(f);
	return m;
}

/* a box of sides 1, 0.8, 0.6, each face an n x n grid of triangle pairs, its edges creases */
static qaws_mesh_3d* tg_box(unsigned int n)
{
	unsigned int N = n + 1, *id = (unsigned int*)malloc(sizeof(unsigned int) * N * N * N), nv = 0, nf = 0, a, s, u, v, k;
	double* x = (double*)malloc(sizeof(double) * 3 * N * N * N);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 36 * n * n);
	static double const size[3] = { 1, 0.8, 0.6 };
	qaws_mesh_3d* m;
	for (k = 0; k < N * N * N; k++)
	{
		unsigned int p[3];
		p[0] = k % N;
		p[1] = (k / N) % N;
		p[2] = k / (N * N);
		id[k] = 0xFFFFFFFFu;
		if (p[0] % n && p[1] % n && p[2] % n)
			continue;
		id[k] = nv;
		for (a = 0; a < 3; a++)
			x[3 * nv + a] = size[a] * ((double)p[a] / n - 0.5);
		nv++;
	}
	for (a = 0; a < 3; a++)
		for (s = 0; s < 2; s++)
			for (u = 0; u < n; u++)
				for (v = 0; v < n; v++)
				{
					unsigned int q[4], c, b = (a + 1) % 3, d = (a + 2) % 3;
					for (c = 0; c < 4; c++)
					{
						unsigned int p[3];
						p[a] = s * n;
						p[b] = u + (c == 1 || c == 2);
						p[d] = v + (c >= 2);
						q[c] = id[p[0] + N * (p[1] + N * p[2])];
					}
					if (!s)
					{
						c = q[1];
						q[1] = q[3];
						q[3] = c;
					}
					f[3 * nf] = q[0]; f[3 * nf + 1] = q[1]; f[3 * nf + 2] = q[2]; nf++;
					f[3 * nf] = q[0]; f[3 * nf + 1] = q[2]; f[3 * nf + 2] = q[3]; nf++;
				}
	m = tg_make(x, nv, f, nf);
	qaws_mesh_3d_detect_features(m, (qaws_scalar)(50 * M_PI / 180), 0, 0);
	free(id);
	free(x);
	free(f);
	return m;
}

/* the showcase's wavy patch (figure 5) sampled on an n x n grid, each cell split along alternating diagonals */
static qaws_mesh_3d* tg_wavy(unsigned int n)
{
	qaws_vec3 cp[49], *s = (qaws_vec3*)malloc(sizeof(qaws_vec3) * n * n);
	qaws_surface_bspline_desc d;
	qaws_surface* srf = NULL;
	qaws_mesh_3d* out = NULL;
	double* x = (double*)malloc(sizeof(double) * 3 * n * n);
	unsigned int* f = (unsigned int*)malloc(sizeof(unsigned int) * 6 * n * n), i, j, k = 0, got = 0;
	for (i = 0; i < 7; i++)
		for (j = 0; j < 7; j++)
		{
			cp[i * 7 + j].x = (qaws_scalar)(i / 3.0 - 1);
			cp[i * 7 + j].y = (qaws_scalar)(j / 3.0 - 1);
			cp[i * 7 + j].z = (qaws_scalar)(0.55 * sin(1.7 * i) * cos(1.3 * j) + (i == 3 && j == 3 ? 0.6 : 0));
		}
	memset(&d, 0, sizeof d);
	d.u_degree = 3;
	d.v_degree = 3;
	d.control_points = cp;
	d.u_point_count = 7;
	d.v_point_count = 7;
	if (qaws_surface_create_bspline(&srf, &d) == QAWS_STATUS_OK && qaws_surface_sample(s, &got, srf, n, n, n * n) == QAWS_STATUS_OK && got == n * n)
	{
		for (i = 0; i < n * n; i++)
		{
			x[3 * i] = (double)s[i].x;
			x[3 * i + 1] = (double)s[i].y;
			x[3 * i + 2] = (double)s[i].z;
		}
		for (i = 0; i + 1 < n; i++)
			for (j = 0; j + 1 < n; j++)
			{
				unsigned int a = i * n + j, b = (i + 1) * n + j, c = (i + 1) * n + j + 1, e = i * n + j + 1;
				if ((i + j) & 1)
				{
					f[k++] = a; f[k++] = b; f[k++] = c;
					f[k++] = a; f[k++] = c; f[k++] = e;
				}
				else
				{
					f[k++] = a; f[k++] = b; f[k++] = e;
					f[k++] = b; f[k++] = c; f[k++] = e;
				}
			}
		out = tg_make(x, n * n, f, k / 3);
	}
	qaws_surface_destroy(srf);
	free(s);
	free(x);
	free(f);
	return out;
}

/* the largest distance from a vertex to the reference at its footpoint */
static double tg_foot(qaws_mesh_3d const* q, qaws_mesh_reference const* r)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(q), i, k;
	qaws_mesh_location* l = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (nv + 1));
	qaws_scalar *x = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (3 * nv + 3)), *y = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (3 * nv + 3));
	double worst = 0;
	qaws_mesh_3d_get_positions(x, q);
	if (qaws_mesh_3d_get_footpoints(l, q) != QAWS_STATUS_OK || qaws_mesh_reference_evaluate(y, NULL, NULL, r, l, nv) != QAWS_STATUS_OK)
		worst = HUGE_VAL;
	for (i = 0; i < nv && worst < HUGE_VAL; i++)
		for (k = 0; k < 3; k++)
			worst = fmax(worst, fabs((double)x[3 * i + k] - (double)y[3 * i + k]));
	free(l); free(x); free(y);
	return worst;
}

/* quads by method at L (pure), measured; NULL on failure */
static qaws_mesh_3d* tg_quads(qaws_mesh_3d const* m, double L, qaws_quad_method method, tg_res* s)
{
	qaws_quadrangulate_desc d;
	qaws_mesh_reference* r = NULL;
	qaws_mesh_3d* out = NULL;
	qaws_mesh_report rep;
	qaws_mesh_quality q;
	clock_t t0;
	memset(s, 0, sizeof *s);
	memset(&d, 0, sizeof d);
	if (qaws_mesh_reference_create(&r, m, NULL) != QAWS_STATUS_OK)
		return NULL;
	d.method = method;
	d.reference = r;
	d.target_length = (qaws_scalar)L;
	d.pure = 1;
	t0 = clock();
	if (qaws_mesh_quadrangulate(&out, &rep, m, &d) == QAWS_STATUS_OK)
	{
		unsigned int nv = qaws_mesh_3d_get_vertex_count(out), ne = qaws_mesh_3d_get_edge_count(out), nc = qaws_mesh_3d_get_corner_count(out);
		unsigned int nf = qaws_mesh_3d_get_face_count(out), i, f, c0, *e, *deg, *fv, *fs, *across;
		unsigned char* bd;
		s->ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
		qaws_mesh_3d_measure_shape(&s->s, NULL, out);
		qaws_mesh_3d_measure_quality(&q, out);
		s->quads = q.quads;
		s->faces = qaws_mesh_3d_get_face_count(out);
		s->made = rep.iterations;
		s->valid = qaws_mesh_3d_validate(out) == QAWS_STATUS_OK;
		s->foot = tg_foot(out, r);
		/* border vertices and those of valence 3 among them (a straight border of a grid) */
		e = (unsigned int*)malloc(sizeof(unsigned int) * (2 * ne + 2));
		deg = (unsigned int*)calloc(nv + 1, sizeof(unsigned int));
		bd = (unsigned char*)calloc(nv + 1, 1);
		qaws_mesh_3d_get_edges(e, out);
		for (i = 0; i < ne; i++)
		{
			deg[e[2 * i]]++;
			deg[e[2 * i + 1]]++;
		}
		fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
		fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
		across = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
		qaws_mesh_3d_get_faces(fv, fs, out);
		qaws_mesh_3d_get_face_adjacency(across, out);
		for (f = 0, c0 = 0; f < nf; c0 += fs[f], f++)
			for (i = 0; i < fs[f]; i++)
				if (across[c0 + i] == QAWS_MESH_NONE)
					bd[fv[c0 + i]] = bd[fv[c0 + (i + 1) % fs[f]]] = 1;
		free(fv); free(fs); free(across);
		for (i = 0; i < nv; i++)
			if (bd[i])
			{
				s->border++;
				s->border3 += deg[i] == 3;
			}
		free(e); free(deg); free(bd);
	}
	qaws_mesh_reference_destroy(r);
	return out;
}

static void tg_print(char const* name, char const* method, tg_res const* s)
{
	printf("  %s %s: %u faces (%u quads), irregular %u, SJ min %.3f mean %.3f, AD %.2f deg, ED %.3f, %s, %.0f ms\n", name, method, s->faces,
		s->quads, s->s.irregular, (double)s->s.min_jacobian, (double)s->s.mean_jacobian, (double)s->s.angle_deviation, (double)s->s.edge_deviation,
		s->made ? "by the map" : "INSTANT's", s->ms);
}

/* ---- the map ------------------------------------------------------------ */

/* a quad torus's straight edge loops (through valence 4 vertices, the edge across from the one it came by): each one's
   turns round the hole and round the tube (radius R); those not once round exactly one of them (a spiral) counted */
static void tg_windings(qaws_mesh_3d const* q, double R, unsigned int* loops, unsigned int* spirals)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(q), nc = qaws_mesh_3d_get_corner_count(q), f, i, s, j;
	unsigned int *fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
	unsigned int *c0 = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1)), *pos = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
	unsigned char* seen = (unsigned char*)calloc(nc + 1, 1);
	*loops = *spirals = 0;
	qaws_mesh_3d_get_faces(fv, fs, q);
	for (f = 0, i = 0; f < nf; f++)
		for (j = 0; j < fs[f]; j++, i++)
		{
			c0[i] = i - j;   /* (each corner's face start and place in it) */
			pos[i] = j | (fs[f] << 8);
		}
	for (s = 0; s < nc; s++)
	{
		unsigned int e = s, n = 0;
		double wt = 0, wp = 0;
		if (seen[s])
			continue;
		do
		{
			unsigned int k = pos[e] & 255, m = pos[e] >> 8, a = fv[e], b = fv[c0[e] + (k + 1) % m], c = fv[c0[e] + (k + 2) % m], g;
			qaws_vec3 A = qaws_mesh_3d_get_vertex(q, a), B = qaws_mesh_3d_get_vertex(q, b);
			double dt = atan2(B.y, B.x) - atan2(A.y, A.x), dp = atan2(B.z, sqrt(B.x * B.x + B.y * B.y) - R) - atan2(A.z, sqrt(A.x * A.x + A.y * A.y) - R);
			wt += dt - 2 * M_PI * floor(dt / (2 * M_PI) + 0.5);
			wp += dp - 2 * M_PI * floor(dp / (2 * M_PI) + 0.5);
			seen[e] = 1;
			/* the face across b -> c (it has c -> b), then the corner after b there: the straight way on; b -> a seen too */
			for (g = 0; g < nc; g++)
				if (fv[g] == c && fv[c0[g] + ((pos[g] & 255) + 1) % (pos[g] >> 8)] == b)
					break;
			for (i = 0; i < nc; i++)
				if (fv[i] == b && fv[c0[i] + ((pos[i] & 255) + 1) % (pos[i] >> 8)] == a)
					seen[i] = 1;
			if (g == nc)
				break;
			{
				unsigned int d = fv[c0[g] + ((pos[g] & 255) + 2) % (pos[g] >> 8)];
				for (i = 0; i < nc; i++)
					if (fv[i] == b && fv[c0[i] + ((pos[i] & 255) + 1) % (pos[i] >> 8)] == d)
						break;
			}
			e = i;
		} while (e < nc && e != s && ++n < nc);
		{
			int it = (int)floor(wt / (2 * M_PI) + 0.5), ip = (int)floor(wp / (2 * M_PI) + 0.5);
			(*loops)++;
			*spirals += e != s || !((abs(it) == 1 && ip == 0) || (it == 0 && abs(ip) == 1));
		}
	}
	free(fv); free(fs); free(c0); free(pos); free(seen);
}

/* tori of uneven radii whose flow's periods made spirals (a seam the quads cross shifted): the periods quantized
   together, every edge loop closes once round the hole or once round the tube */
static void test_periods(void)
{
	static double const cases[2][3] = { { 0.7, 0.22, 0 }, { 0.5, 0.3, 0 } };
	unsigned int k;
	for (k = 0; k < 2; k++)
	{
		qaws_mesh_3d* m = tg_torus(96, 48, 1, cases[k][0]);
		tg_res s;
		qaws_mesh_3d* q = tg_quads(m, cases[k][1], QAWS_QUADS_GLOBAL, &s);
		unsigned int loops = 0, spirals = 1;
		char msg[300];
		if (q)
			tg_windings(q, 1, &loops, &spirals);
		sprintf(msg, "a torus of radii 1 and %.1f at %.2f: %u quads, irregular %u, %u edge loops, %u spiral", cases[k][0], cases[k][1], s.quads,
			s.s.irregular, loops, spirals);
		printf("  %s\n", msg);
		TEST_ASSERT(q && s.made == 1 && s.s.irregular == 0 && loops > 10 && spirals == 0, msg);
		qaws_mesh_3d_destroy(q);
		qaws_mesh_3d_destroy(m);
	}
}

static void test_map(void)
{
	qaws_mesh_3d* blob = tg_can();
	qaws_quadrangulate_desc d;
	qg_map map;
	unsigned int i, neg = 0;
	double amin = HUGE_VAL;
	char msg[200];
	memset(&d, 0, sizeof d);
	d.method = QAWS_QUADS_GLOBAL;
	d.target_length = (qaws_scalar)0.14;
	d.pure = 1;
	TEST_ASSERT(qg_parametrization(blob, &d, &map) == QAWS_STATUS_OK && map.nt > 0, "can: the map made");
	for (i = 0; i < map.nt; i++)
	{
		double const* u = map.uv + 6 * (size_t)i;
		double a = 0.5 * ((u[2] - u[0]) * (u[5] - u[1]) - (u[3] - u[1]) * (u[4] - u[0]));
		amin = fmin(amin, a);
		neg += a <= 0;
	}
	sprintf(msg, "can's map: %u triangles, none folded (%u), smallest uv area %.2e", map.nt, neg, amin);
	TEST_ASSERT(map.nt > 0 && neg == 0, msg);
	qg_map_free(&map);
	qaws_mesh_3d_destroy(blob);
}

/* ---- through qaws_mesh_quadrangulate ------------------------------------ */

static void test_shapes(void)
{
	static char const* const name[5] = { "torus", "box", "blob", "can", "wavy patch" };
	static double const L[5] = { 0.3, 0.1, 0.11, 0.14, 0.1 };
	qaws_mesh_3d* in[5] = { NULL, NULL, NULL, NULL, NULL };
	/* INSTANT's numbers at these targets (measured with it run alongside, heavy): irregular, AD */
	static unsigned int const q_irregular[5] = { 0, 0, 8, 8, 2 };
	static double const q_ad[5] = { 3.60, 0.04, 3.27, 2.23, 3.22 };
	tg_res g[5], q[5];
	unsigned int i;
	int heavy = test_heavy();
	char msg[600];
	memset(g, 0, sizeof g);
	memset(q, 0, sizeof q);
	for (i = 0; i < 5; i++)
	{
		qaws_mesh_3d *a, *b;
		if (!heavy && (i == 0 || i == 2))
			continue;
		in[i] = i == 0 ? tg_torus(30, 12, 2, 0.75) : i == 1 ? tg_box(4) : i == 2 ? tg_blob() : i == 3 ? tg_can() : tg_wavy(41);
		a = in[i] ? tg_quads(in[i], L[i], QAWS_QUADS_GLOBAL, &g[i]) : NULL;
		b = in[i] && heavy ? tg_quads(in[i], L[i], QAWS_QUADS_INSTANT, &q[i]) : NULL;
		if (!heavy)
		{
			q[i].s.irregular = q_irregular[i];
			q[i].s.angle_deviation = (qaws_scalar)q_ad[i];
		}
		if (!a || (heavy && !b))
			TEST_ASSERT(0, "quadrangulated");
		else
		{
			tg_print(name[i], "GLOBAL", &g[i]);
			if (heavy)
				tg_print(name[i], "INSTANT", &q[i]);
		}
		qaws_mesh_3d_destroy(a);
		qaws_mesh_3d_destroy(b);
		qaws_mesh_3d_destroy(in[i]);
		if (!a || (heavy && !b))
			return;
		/* every result: made by the map, valid, all quads, on the reference, no quad turned; no worse than INSTANT */
		sprintf(msg, "%s at %.2f: by the map %u, valid %d, %u of %u quads, footpoints to %.1e, SJ min %.3f; irregular %u (INSTANT %u), "
			"AD %.2f (INSTANT %.2f)", name[i], L[i], g[i].made, g[i].valid, g[i].quads, g[i].faces, g[i].foot, (double)g[i].s.min_jacobian,
			g[i].s.irregular, q[i].s.irregular, (double)g[i].s.angle_deviation, (double)q[i].s.angle_deviation);
		TEST_ASSERT(g[i].made == 1 && g[i].valid && g[i].quads == g[i].faces && g[i].foot < TEST_BOUND(1e-9, 32) && g[i].s.min_jacobian > 0
			&& g[i].s.irregular <= q[i].s.irregular && g[i].s.angle_deviation <= 1.02 * q[i].s.angle_deviation + 0.01, msg);
	}
	/* the shapes' own counts: the torus a grid, the box's faces grids (its corners regular by the feature-aware measure), the
	   blob's irregular vertices its 8 cones, the patch's border along isolines (border vertices of valence 3 but its 4
	   corners) */
	sprintf(msg, "torus %u irregular (0), box %u (0), blob %u (8), can %u (8), wavy patch border: %u of %u at valence 3 (all but 4)",
		g[0].s.irregular, g[1].s.irregular, g[2].s.irregular, g[3].s.irregular, g[4].border3, g[4].border);
	TEST_ASSERT((!heavy || (g[0].s.irregular == 0 && g[2].s.irregular == 8)) && g[1].s.irregular == 0 && g[3].s.irregular <= 8
		&& g[4].border3 + 4 == g[4].border, msg);
}

/* ---- a sweep of targets ------------------------------------------------- */

/* the shapes at 0.75, 1.25, 1.6 and 2.2 times their targets (heavy; by default two the flow's integers alone gave no
   map for: the box at 2.2, its corner steps read in the wrong plane; the can at 1.25, folds by its caps' cones, made by
   the quantized integers; heavy adds the patch at 1.25, a cone wound a full turn short): each made by the map, valid, all
   quads, no quad turned */
static void test_sweep(void)
{
	static char const* const name[5] = { "torus", "box", "blob", "can", "wavy patch" };
	static double const L[5] = { 0.3, 0.1, 0.11, 0.14, 0.1 }, all[4] = { 0.75, 1.25, 1.6, 2.2 };
	static unsigned char const light[5][4] = { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 0, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 0, 0 } };
	unsigned int i, j, runs = 0, made = 0, good = 0;
	int heavy = test_heavy();
	char msg[300], fail[400] = "";
	for (i = 0; i < 5; i++)
	{
		qaws_mesh_3d* in;
		for (j = 0; j < 4 && !heavy && !light[i][j]; j++) {}
		if (j == 4)
			continue;
		in = i == 0 ? tg_torus(30, 12, 2, 0.75) : i == 1 ? tg_box(4) : i == 2 ? tg_blob() : i == 3 ? tg_can() : tg_wavy(41);
		for (j = 0; j < 4 && in; j++)
		{
			tg_res s;
			qaws_mesh_3d* q;
			if (!heavy && !light[i][j])
				continue;
			q = tg_quads(in, L[i] * all[j], QAWS_QUADS_GLOBAL, &s);
			runs++;
			made += q && s.made == 1;
			good += q && s.valid && s.quads == s.faces && s.faces && s.s.min_jacobian > 0;
			if (!q || s.made != 1)
				sprintf(fail + strlen(fail), " %s x%.2f", name[i], all[j]);
			if (q)
				printf("  sweep %s x%.2f: %u quads, irregular %u, SJ min %.3f, AD %.2f, %s\n", name[i], all[j], s.faces, s.s.irregular,
					(double)s.s.min_jacobian, (double)s.s.angle_deviation, s.made ? "by the map" : "INSTANT's");
			qaws_mesh_3d_destroy(q);
		}
		qaws_mesh_3d_destroy(in);
	}
	sprintf(msg, "sweep: %u of %u by the map (missed:%s), %u valid, all quads, SJ > 0", made, runs, fail, good);
	TEST_ASSERT(runs && made == runs && good == runs, msg);
}

int test_136_quad_global_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 136: quads by integer-grid maps\n");
	if (test_heavy())
		test_map();   /* (the can by the map is checked by default through test_shapes) */
	test_shapes();
	test_sweep();
	test_periods();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
