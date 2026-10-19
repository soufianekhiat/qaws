/*
 * Octahedral frame fields on tetrahedral meshes (qaws_mesh_field.h), after
 * Ray, Sokolov and Levy 2016: each frame as its function on the sphere
 * (band-4 spherical harmonics, 9 numbers), the smoothest such field over
 * the cells' adjacency found by one sparse linear solve with the border's
 * normals held (a border cell's function confined to the frames turning
 * about its normal: 2 numbers), each cell's frame then the one nearest its
 * function (gradient ascent on the rotations, from a neighbour's frame,
 * breadth first from the border); then Gauss-Seidel sweeps over the
 * facets' adjacency, neighbours matched by the signed permutation of their
 * axes nearest the cell's. Singular edges where the matchings round an edge
 * compose to a turn.
 *
 * shortcut: the projection and sweeps are local (Ray et al. follow with a
 * non-linear smoothing of the whole field; Palmer, Bommes and Solomon 2020
 * project algebraically), singular curves follow the tetrahedra's edges.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* matrices are 3 x 3, column-major: m[3 j + i] is row i of column j (column j the axis j) */

static double ff_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void ff_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double ff_unit(double* v)
{
	double l = sqrt(ff_dot(v, v));
	if (l > 0)
	{
		v[0] /= l; v[1] /= l; v[2] /= l;
	}
	return l;
}

/* the signed permutation P (P[i] = j + 1 or -(j + 1)) taking b's axes onto a's: axis i of a ~ sign * axis j of b, a
   rotation (determinant 1) nearest b^T a */
static void ff_match(double const* a, double const* b, int* P)
{
	static int const perm[6][3] = { { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 }, { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 } };
	double D[3][3], best = -1;
	int p, i, bp = 0;
	for (i = 0; i < 3; i++)
	{
		int j;
		for (j = 0; j < 3; j++)
			D[i][j] = ff_dot(a + 3 * i, b + 3 * j);
	}
	for (p = 0; p < 6; p++)
	{
		double s = fabs(D[0][perm[p][0]]) + fabs(D[1][perm[p][1]]) + fabs(D[2][perm[p][2]]);
		if (s > best)
		{
			best = s;
			bp = p;
		}
	}
	for (i = 0; i < 3; i++)
		P[i] = D[i][perm[bp][i]] < 0 ? -(perm[bp][i] + 1) : perm[bp][i] + 1;
	/* a proper rotation: an odd permutation with even signs (or the reverse) flips the weakest axis */
	{
		int odd = bp == 1 || bp == 2 || bp == 5, neg = (P[0] < 0) + (P[1] < 0) + (P[2] < 0), w = 0;
		if ((odd + neg) & 1)
		{
			for (i = 1; i < 3; i++)
				if (fabs(D[i][perm[bp][i]]) < fabs(D[w][perm[bp][w]]))
					w = i;
			P[w] = -P[w];
		}
	}
}

/* b's axes matched to a's, added to m */
static void ff_add_matched(double const* a, double const* b, double* m)
{
	int P[3], i, k;
	ff_match(a, b, P);
	for (i = 0; i < 3; i++)
	{
		int j = (P[i] < 0 ? -P[i] : P[i]) - 1;
		double s = P[i] < 0 ? -1 : 1;
		for (k = 0; k < 3; k++)
			m[3 * i + k] += s * b[3 * j + k];
	}
}

/* the rotation nearest m (polar decomposition, Newton: R <- (R + R^-T) / 2) */
static void ff_polar(double* m)
{
	int it;
	for (it = 0; it < 12; it++)
	{
		double c[9], det, inv[9];
		int k;
		/* the cofactors: R^-T = cof(R) / det(R) */
		ff_cross(m + 3, m + 6, c);
		ff_cross(m + 6, m, c + 3);
		ff_cross(m, m + 3, c + 6);
		det = ff_dot(m, c);
		if (!(fabs(det) > 1e-300))
			break;
		for (k = 0; k < 9; k++)
		{
			inv[k] = c[k] / det;
			m[k] = 0.5 * (m[k] + inv[k]);
		}
	}
	ff_unit(m);
	ff_unit(m + 3);
	ff_cross(m, m + 3, m + 6);
	ff_unit(m + 6);
	ff_cross(m + 6, m, m + 3);
}

/* a frame with axis 0 along n (and axis 1 along n2 made orthogonal, when given) */
static void ff_from_normals(double const* n, double const* n2, double* m)
{
	double t[3];
	memcpy(m, n, sizeof(double) * 3);
	ff_unit(m);
	if (n2)
		memcpy(t, n2, sizeof t);
	else
	{
		t[0] = fabs(m[0]) < 0.6 ? 1 : 0;
		t[1] = fabs(m[0]) < 0.6 ? 0 : 1;
		t[2] = 0;
	}
	{
		double d = ff_dot(t, m);
		t[0] -= d * m[0]; t[1] -= d * m[1]; t[2] -= d * m[2];
	}
	ff_unit(t);
	memcpy(m + 3, t, sizeof t);
	ff_cross(m, m + 3, m + 6);
}

/* the averaged frame bent back to the cell's constraints: the axis nearest each normal set onto it */
static void ff_constrain(double* m, double const* n, unsigned int nn)
{
	double keep[9];
	int i0 = 0, i1 = -1, i;
	if (!nn)
		return;
	memcpy(keep, m, sizeof keep);
	for (i = 1; i < 3; i++)
		if (fabs(ff_dot(m + 3 * i, n)) > fabs(ff_dot(m + 3 * i0, n)))
			i0 = i;
	if (nn > 1)
		for (i = 0; i < 3; i++)
			if (i != i0 && (i1 < 0 || fabs(ff_dot(m + 3 * i, n + 3)) > fabs(ff_dot(m + 3 * i1, n + 3))))
				i1 = i;
	{
		double a0[3], a1[3], a2[3], s = ff_dot(m + 3 * i0, n) < 0 ? -1 : 1;
		int k, i2;
		for (k = 0; k < 3; k++)
			a0[k] = s * n[k];
		if (nn > 1)
		{
			double s1 = ff_dot(m + 3 * i1, n + 3) < 0 ? -1 : 1, d;
			for (k = 0; k < 3; k++)
				a1[k] = s1 * n[3 + k];
			d = ff_dot(a1, a0);
			for (k = 0; k < 3; k++)
				a1[k] -= d * a0[k];
		}
		else
		{
			i1 = (i0 + 1) % 3;
			{
				double d = ff_dot(keep + 3 * i1, a0);
				for (k = 0; k < 3; k++)
					a1[k] = keep[3 * i1 + k] - d * a0[k];
			}
		}
		if (ff_unit(a1) == 0)
		{
			ff_from_normals(a0, NULL, m);
			return;
		}
		i2 = 3 - i0 - i1;
		ff_cross(a0, a1, a2);
		/* the third axis keeps the handedness: (i0, i1, i2) an even order of (0, 1, 2) or not */
		if (!((i1 - i0 + 3) % 3 == 1))
		{
			a2[0] = -a2[0]; a2[1] = -a2[1]; a2[2] = -a2[2];
		}
		memcpy(m + 3 * i0, a0, sizeof a0);
		memcpy(m + 3 * i1, a1, sizeof a1);
		memcpy(m + 3 * i2, a2, sizeof a2);
	}
}

/* ======================================================================== */
/*  Frames as functions: band-4 spherical harmonics                         */
/* ======================================================================== */

/* A frame (axes a_i) is the function sum_i (a_i . x)^4 on the sphere: a
   constant and a band-4 harmonic, so 9 coefficients, s_k = sum_i Y_k(a_i)
   (Funk-Hecke, up to one constant), a linear space where averaging is easy. */

/* the 15 quartic monomials x^a y^b z^c */
static unsigned char const ff_pw[15][3] = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 }, { 3, 1, 0 }, { 3, 0, 1 }, { 1, 3, 0 }, { 0, 3, 1 }, { 1, 0, 3 },
	{ 0, 1, 3 }, { 2, 2, 0 }, { 2, 0, 2 }, { 0, 2, 2 }, { 2, 1, 1 }, { 1, 2, 1 }, { 1, 1, 2 } };

static double ff_dfact(int n)
{
	double r = 1;
	for (; n > 1; n -= 2)
		r *= n;
	return r;
}

/* the integral of x^a y^b z^c over the unit sphere */
static double ff_sphere(int a, int b, int c)
{
	if ((a | b | c) & 1)
		return 0;
	return 4 * 3.14159265358979323846 * ff_dfact(a - 1) * ff_dfact(b - 1) * ff_dfact(c - 1) / ff_dfact(a + b + c + 1);
}

/* the real harmonics of band 4 (35z^4 - 30z^2 r^2 + 3r^4 and its kin), orthonormal on the sphere */
static void ff_basis(double Y[9][15])
{
	static signed char const h[9][15] = {
		{ 3, 3, 8, 0, 0, 0, 0, 0, 0, 6, -24, -24, 0, 0, 0 },     /* 35z^4 - 30z^2 r^2 + 3r^4 */
		{ 0, 0, 0, 0, -3, 0, 0, 4, 0, 0, 0, 0, 0, -3, 0 },       /* xz (7z^2 - 3r^2) */
		{ 0, 0, 0, 0, 0, 0, -3, 0, 4, 0, 0, 0, -3, 0, 0 },       /* yz (7z^2 - 3r^2) */
		{ -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 0 },       /* (x^2 - y^2)(7z^2 - r^2) */
		{ 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 6 },       /* xy (7z^2 - r^2) */
		{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0 },        /* xz (x^2 - 3y^2) */
		{ 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 3, 0, 0 },        /* yz (3x^2 - y^2) */
		{ 1, 1, 0, 0, 0, 0, 0, 0, 0, -6, 0, 0, 0, 0, 0 },        /* x^4 - 6x^2 y^2 + y^4 */
		{ 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };      /* xy (x^2 - y^2) */
	int k, i, j;
	for (k = 0; k < 9; k++)
	{
		double n = 0;
		for (i = 0; i < 15; i++)
			for (j = 0; j < 15; j++)
				n += h[k][i] * h[k][j] * ff_sphere(ff_pw[i][0] + ff_pw[j][0], ff_pw[i][1] + ff_pw[j][1], ff_pw[i][2] + ff_pw[j][2]);
		for (i = 0; i < 15; i++)
			Y[k][i] = h[k][i] / sqrt(n);
	}
}

/* a quartic (15 coefficients) at a, and its gradient (g may be NULL) */
static double ff_quartic(double const* p, double const* a, double* g)
{
	double v = 0, pw[3][5];
	int i, k;
	for (k = 0; k < 3; k++)
	{
		pw[k][0] = 1;
		for (i = 1; i < 5; i++)
			pw[k][i] = pw[k][i - 1] * a[k];
	}
	if (g)
		g[0] = g[1] = g[2] = 0;
	for (i = 0; i < 15; i++)
	{
		int ex = ff_pw[i][0], ey = ff_pw[i][1], ez = ff_pw[i][2];
		v += p[i] * pw[0][ex] * pw[1][ey] * pw[2][ez];
		if (g)
		{
			if (ex)
				g[0] += p[i] * ex * pw[0][ex - 1] * pw[1][ey] * pw[2][ez];
			if (ey)
				g[1] += p[i] * ey * pw[0][ex] * pw[1][ey - 1] * pw[2][ez];
			if (ez)
				g[2] += p[i] * ez * pw[0][ex] * pw[1][ey] * pw[2][ez - 1];
		}
	}
	return v;
}

/* a frame's coefficients */
static void ff_frame_sh(double Y[9][15], double const* R, double* s)
{
	int k, a;
	for (k = 0; k < 9; k++)
	{
		s[k] = 0;
		for (a = 0; a < 3; a++)
			s[k] += ff_quartic(Y[k], R + 3 * a, NULL);
	}
}

/* v turned about axis w by |w| (Rodrigues) */
static void ff_turn(double const* w, double* v)
{
	double t = sqrt(ff_dot(w, w)), k[3], c, s, kv[3], d;
	int i;
	if (!(t > 0))
		return;
	for (i = 0; i < 3; i++)
		k[i] = w[i] / t;
	c = cos(t);
	s = sin(t);
	ff_cross(k, v, kv);
	d = ff_dot(k, v);
	for (i = 0; i < 3; i++)
		v[i] = v[i] * c + kv[i] * s + k[i] * d * (1 - c);
}

/* the frame nearest coefficients q from R: gradient ascent of sum_i F_q(a_i) on the rotations (each axis turned by the
   torque sum_i a_i x grad F_q(a_i)), steps grown while they gain, halved when they do not */
static void ff_project(double Y[9][15], double const* q, double* R)
{
	double F[15], step = 0.5, val = 0;
	int i, k, it;
	for (i = 0; i < 15; i++)
	{
		F[i] = 0;
		for (k = 0; k < 9; k++)
			F[i] += q[k] * Y[k][i];
	}
	for (k = 0; k < 3; k++)
		val += ff_quartic(F, R + 3 * k, NULL);
	for (it = 0; it < 200 && step > 1e-9; it++)
	{
		double w[3] = { 0, 0, 0 }, trial[9], tv = 0, wl;
		for (k = 0; k < 3; k++)
		{
			double g[3], t[3];
			ff_quartic(F, R + 3 * k, g);
			ff_cross(R + 3 * k, g, t);
			w[0] += t[0];
			w[1] += t[1];
			w[2] += t[2];
		}
		wl = sqrt(ff_dot(w, w));
		if (!(wl > 1e-14))
			break;
		for (k = 0; k < 3; k++)
			w[k] *= step / wl;
		memcpy(trial, R, sizeof trial);
		for (k = 0; k < 3; k++)
		{
			ff_turn(w, trial + 3 * k);
			tv += ff_quartic(F, trial + 3 * k, NULL);
		}
		if (tv > val)
		{
			memcpy(R, trial, sizeof trial);
			val = tv;
			step = fmin(1, step * 1.5);
		}
		else
			step *= 0.5;
	}
}

/* the frames' linear system: per cell its coefficients q = J z + b, z free (9 inside; 2 where one normal holds an axis,
   q = A + B c + C s turning about it; none where two normals fix it, q = A) */
typedef struct ff_sys
{
	unsigned char* kind;
	unsigned int* off;
	double *A, *B, *C;
	unsigned int nz;
} ff_sys;

static void ff_q(ff_sys const* S, unsigned int c, double const* z, int with_b, double* q)
{
	unsigned int k;
	double const *A = S->A + 9 * (size_t)c, *B = S->B + 9 * (size_t)c, *C = S->C + 9 * (size_t)c;
	for (k = 0; k < 9; k++)
		q[k] = (with_b && S->kind[c] ? A[k] : 0)
			+ (!z ? 0 : S->kind[c] == 0 ? z[S->off[c] + k] : S->kind[c] == 1 ? B[k] * z[S->off[c]] + C[k] * z[S->off[c] + 1] : 0);
}

static void ff_qt(ff_sys const* S, unsigned int c, double const* d, double sign, double* out)
{
	unsigned int k;
	if (S->kind[c] == 0)
		for (k = 0; k < 9; k++)
			out[S->off[c] + k] += sign * d[k];
	else if (S->kind[c] == 1)
		for (k = 0; k < 9; k++)
		{
			out[S->off[c]] += sign * S->B[9 * (size_t)c + k] * d[k];
			out[S->off[c] + 1] += sign * S->C[9 * (size_t)c + k] * d[k];
		}
}

/* out = H v (with_b 0), or out = -sum D^T (b_c - b_o) (with_b 1, v unused): the smoothness sum over facets of |q_c - q_o|^2 */
static void ff_apply(ff_sys const* S, qaws_volume_mesh const* m, double const* v, int with_b, double* out)
{
	unsigned int c, f;
	memset(out, 0, sizeof(double) * S->nz);
	for (c = 0; c < m->nc; c++)
		for (f = 0; f < 4; f++)
		{
			unsigned int o = m->fo[m->f0[c] + f], u, k;
			double qc[9], qo[9], d[9];
			if (o == HM_NONE || o / 4 < c)
				continue;
			u = o / 4;
			ff_q(S, c, v, with_b, qc);
			ff_q(S, u, v, with_b, qo);
			for (k = 0; k < 9; k++)
				d[k] = qc[k] - qo[k];
			/* (with_b: v is NULL, so d is b_c - b_o) */
			ff_qt(S, c, d, with_b ? -1 : 1, out);
			ff_qt(S, u, d, with_b ? 1 : -1, out);
		}
}

/* the smoothest coefficients: preconditioned conjugate gradients on the normal equations; 0 on allocation failure */
static int ff_solve(ff_sys const* S, qaws_volume_mesh const* m, double* z)
{
	double *r = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)S->nz + 1))), *p = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), ((size_t)S->nz + 1)));
	double *Hp = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)S->nz + 1))), *dg = (double*)qaws_internal_zalloc((size_t)S->nz + 1, sizeof(double));
	double *zz = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), ((size_t)S->nz + 1))), rz, r0;
	unsigned int c, f, i, it;
	if (!r || !p || !Hp || !dg || !zz)
	{
		qaws_internal_free(r); qaws_internal_free(p); qaws_internal_free(Hp); qaws_internal_free(dg); qaws_internal_free(zz);
		return 0;
	}
	/* the diagonal: a facet adds J^T J's diagonal on each side */
	for (c = 0; c < m->nc; c++)
		for (f = 0; f < 4; f++)
			if (m->fo[m->f0[c] + f] != HM_NONE)
			{
				unsigned int k;
				if (S->kind[c] == 0)
					for (k = 0; k < 9; k++)
						dg[S->off[c] + k] += 1;
				else if (S->kind[c] == 1)
				{
					dg[S->off[c]] += ff_dot(S->B + 9 * (size_t)c, S->B + 9 * (size_t)c) + ff_dot(S->B + 9 * (size_t)c + 3, S->B + 9 * (size_t)c + 3)
						+ ff_dot(S->B + 9 * (size_t)c + 6, S->B + 9 * (size_t)c + 6);
					dg[S->off[c] + 1] += ff_dot(S->C + 9 * (size_t)c, S->C + 9 * (size_t)c) + ff_dot(S->C + 9 * (size_t)c + 3, S->C + 9 * (size_t)c + 3)
						+ ff_dot(S->C + 9 * (size_t)c + 6, S->C + 9 * (size_t)c + 6);
				}
			}
	for (i = 0; i < S->nz; i++)
		dg[i] = dg[i] > 0 ? 1 / dg[i] : 0;
	/* r = b - H z */
	ff_apply(S, m, NULL, 1, r);
	ff_apply(S, m, z, 0, Hp);
	for (i = 0; i < S->nz; i++)
	{
		r[i] -= Hp[i];
		zz[i] = dg[i] * r[i];
		p[i] = zz[i];
	}
	rz = 0;
	for (i = 0; i < S->nz; i++)
		rz += r[i] * zz[i];
	r0 = rz;
	for (it = 0; it < 4000 && rz > 1e-24 * r0 && rz > 0; it++)
	{
		double pHp = 0, alpha, rz2 = 0;
		ff_apply(S, m, p, 0, Hp);
		for (i = 0; i < S->nz; i++)
			pHp += p[i] * Hp[i];
		if (!(pHp > 0))
			break;
		alpha = rz / pHp;
		for (i = 0; i < S->nz; i++)
		{
			z[i] += alpha * p[i];
			r[i] -= alpha * Hp[i];
			zz[i] = dg[i] * r[i];
			rz2 += r[i] * zz[i];
		}
		for (i = 0; i < S->nz; i++)
			p[i] = zz[i] + rz2 / rz * p[i];
		rz = rz2;
	}
	qaws_internal_free(r); qaws_internal_free(p); qaws_internal_free(Hp); qaws_internal_free(dg); qaws_internal_free(zz);
	return 1;
}

qaws_status qaws_volume_mesh_frame_field(struct qaws_volume_mesh const* mesh, qaws_frame_field_desc const* desc, qaws_scalar* out_frames,
	unsigned char* out_singular)
{
	double* R;
	unsigned int c;
	qaws_status st;
	if (!mesh || (mesh->nc && !out_frames))
		return QAWS_STATUS_INVALID_ARGUMENT;
	R = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (9 * (size_t)mesh->nc + 1)));
	if (!R)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = ff_frames(mesh, desc && desc->sweeps ? desc->sweeps : 60, R, out_singular);
	for (c = 0; st == QAWS_STATUS_OK && c < 9 * mesh->nc; c++)
		out_frames[c] = (qaws_scalar)R[c];
	qaws_internal_free(R);
	return st;
}

int ff_comb(qaws_volume_mesh const* m, double* R)
{
	unsigned int *queue = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->nc + 1))), nq = 0, head, c, f;
	unsigned char* set = (unsigned char*)qaws_internal_zalloc((size_t)m->nc + 1, 1);
	int ok = queue && set;
	for (c = 0; ok && c < m->nc; c++)
	{
		if (set[c])
			continue;
		set[c] = 1;
		queue[nq++] = c;
		for (head = nq - 1; head < nq; head++)
			for (f = 0; f < 4; f++)
			{
				unsigned int o = m->fo[m->f0[queue[head]] + f], u, i, k;
				double const* a = R + 9 * (size_t)queue[head];
				double* b;
				double old[9];
				int P[3];
				if (o == HM_NONE || set[o / 4])
					continue;
				u = o / 4;
				b = R + 9 * (size_t)u;
				/* u's axes reordered and signed onto the cell it was reached from */
				ff_match(a, b, P);
				memcpy(old, b, sizeof old);
				for (i = 0; i < 3; i++)
				{
					int j = (P[i] < 0 ? -P[i] : P[i]) - 1;
					for (k = 0; k < 3; k++)
						b[3 * i + k] = (P[i] < 0 ? -1 : 1) * old[3 * j + k];
				}
				set[u] = 1;
				queue[nq++] = u;
			}
	}
	/* every facet agrees now, or a singular edge or a hole turns the frames */
	for (c = 0; ok && c < m->nc; c++)
		for (f = 0; f < 4 && ok; f++)
		{
			unsigned int o = m->fo[m->f0[c] + f];
			int P[3];
			if (o == HM_NONE)
				continue;
			ff_match(R + 9 * (size_t)c, R + 9 * (size_t)(o / 4), P);
			ok = P[0] == 1 && P[1] == 2 && P[2] == 3;
		}
	qaws_internal_free(queue);
	qaws_internal_free(set);
	return ok;
}

qaws_status ff_frames(qaws_volume_mesh const* m, unsigned int sweeps, double* R, unsigned char* out_singular)
{
	double* nrm = NULL;
	unsigned char* ncon = NULL;
	unsigned int c, s, i, k;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (c = 0; c < m->nc; c++)
		if (m->type[c] != QAWS_CELL_TETRAHEDRON)
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
	nrm = (double*)qaws_internal_scratch_zalloc(6 * (size_t)m->nc + 1, sizeof(double));
	ncon = (unsigned char*)qaws_internal_zalloc((size_t)m->nc + 1, 1);
	if (!nrm || !ncon)
		goto done;
	/* the border's normals in its cells: one, or two meeting near square */
	for (c = 0; c < m->nc; c++)
	{
		unsigned int f;
		for (f = 0; f < 4; f++)
		{
			unsigned int n;
			unsigned char const* loc;
			double u[3], w[3], q[3];
			double const *a, *b, *d;
			if (m->fo[m->f0[c] + f] != HM_NONE)
				continue;
			loc = vm_facet(QAWS_CELL_TETRAHEDRON, f, &n);
			a = m->x + 3 * (size_t)m->cv[m->c0[c] + loc[0]];
			b = m->x + 3 * (size_t)m->cv[m->c0[c] + loc[1]];
			d = m->x + 3 * (size_t)m->cv[m->c0[c] + loc[2]];
			for (k = 0; k < 3; k++)
			{
				u[k] = b[k] - a[k];
				w[k] = d[k] - a[k];
			}
			ff_cross(u, w, q);
			if (ff_unit(q) == 0)
				continue;
			if (ncon[c] == 0)
			{
				memcpy(nrm + 6 * (size_t)c, q, sizeof q);
				ncon[c] = 1;
			}
			else if (ncon[c] == 1 && fabs(ff_dot(q, nrm + 6 * (size_t)c)) < 0.5)
			{
				memcpy(nrm + 6 * (size_t)c + 3, q, sizeof q);
				ncon[c] = 2;
			}
		}
	}
	/* the start: the smoothest field of frame functions with the border's normals held (one linear solve over every
	   cell), then each cell's frame the one nearest its function, sought from a neighbour's (breadth first from the
	   border) */
	{
		double Y[9][15], *z = NULL;
		ff_sys S;
		unsigned int *queue = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->nc + 1))), nq = 0, head;
		unsigned char* set = (unsigned char*)qaws_internal_zalloc((size_t)m->nc + 1, 1);
		int ok;
		memset(&S, 0, sizeof S);
		ff_basis(Y);
		S.kind = ncon;
		S.off = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->nc + 1)));
		S.A = (double*)qaws_internal_zalloc(9 * (size_t)m->nc + 1, sizeof(double));
		S.B = (double*)qaws_internal_zalloc(9 * (size_t)m->nc + 1, sizeof(double));
		S.C = (double*)qaws_internal_zalloc(9 * (size_t)m->nc + 1, sizeof(double));
		ok = queue && set && S.off && S.A && S.B && S.C;
		for (c = 0; ok && c < m->nc; c++)
		{
			double R0[9], Rt[9], s0[9], s1[9], s2[9];
			S.off[c] = S.nz;
			S.nz += ncon[c] == 0 ? 9 : ncon[c] == 1 ? 2 : 0;
			if (!ncon[c])
				continue;
			ff_from_normals(nrm + 6 * (size_t)c, ncon[c] > 1 ? nrm + 6 * (size_t)c + 3 : NULL, R0);
			if (ncon[c] > 1)
			{
				ff_frame_sh(Y, R0, S.A + 9 * (size_t)c);
				continue;
			}
			/* about the normal: s(t) = A + B cos 4t + C sin 4t, from t = 0, pi / 8, pi / 4 */
			ff_frame_sh(Y, R0, s0);
			for (i = 0; i < 2; i++)
			{
				double t = (i + 1) * 3.14159265358979323846 / 8;
				memcpy(Rt, R0, sizeof Rt);
				for (k = 0; k < 3; k++)
				{
					Rt[3 + k] = R0[3 + k] * cos(t) + R0[6 + k] * sin(t);
					Rt[6 + k] = R0[6 + k] * cos(t) - R0[3 + k] * sin(t);
				}
				ff_frame_sh(Y, Rt, i ? s2 : s1);
			}
			for (k = 0; k < 9; k++)
			{
				S.A[9 * (size_t)c + k] = 0.5 * (s0[k] + s2[k]);
				S.B[9 * (size_t)c + k] = 0.5 * (s0[k] - s2[k]);
				S.C[9 * (size_t)c + k] = s1[k] - S.A[9 * (size_t)c + k];
			}
		}
		z = ok ? (double*)qaws_internal_zalloc((size_t)S.nz + 1, sizeof(double)) : NULL;
		ok = ok && z && ff_solve(&S, m, z);
		/* the border's frames: fixed, or turned about the normal by the angle the solve chose */
		for (c = 0; ok && c < m->nc; c++)
			if (ncon[c])
			{
				double* f = R + 9 * (size_t)c;
				ff_from_normals(nrm + 6 * (size_t)c, ncon[c] > 1 ? nrm + 6 * (size_t)c + 3 : NULL, f);
				if (ncon[c] == 1)
				{
					double t = atan2(z[S.off[c] + 1], z[S.off[c]]) / 4, u[3], v[3];
					memcpy(u, f + 3, sizeof u);
					memcpy(v, f + 6, sizeof v);
					for (k = 0; k < 3; k++)
					{
						f[3 + k] = u[k] * cos(t) + v[k] * sin(t);
						f[6 + k] = v[k] * cos(t) - u[k] * sin(t);
					}
				}
				set[c] = 1;
				queue[nq++] = c;
			}
		/* inside: breadth first, each cell's search from the cell it was reached from (the world's axes for a part the border
		   does not reach) */
		for (head = 0; ok && nq < m->nc; head++)
		{
			unsigned int f, from;
			if (head == nq)
			{
				for (c = 0; set[c]; c++) {}
				memset(R + 9 * (size_t)c, 0, sizeof(double) * 9);
				R[9 * (size_t)c] = R[9 * (size_t)c + 4] = R[9 * (size_t)c + 8] = 1;
				{
					double q[9];
					ff_q(&S, c, z, 1, q);
					ff_project(Y, q, R + 9 * (size_t)c);
				}
				set[c] = 1;
				queue[nq++] = c;
			}
			from = queue[head];
			for (f = 0; f < 4; f++)
			{
				unsigned int o = m->fo[m->f0[from] + f], u;
				double q[9];
				if (o == HM_NONE || set[o / 4])
					continue;
				u = o / 4;
				memcpy(R + 9 * (size_t)u, R + 9 * (size_t)from, sizeof(double) * 9);
				ff_q(&S, u, z, 1, q);
				ff_project(Y, q, R + 9 * (size_t)u);
				set[u] = 1;
				queue[nq++] = u;
			}
		}
		qaws_internal_free(z); qaws_internal_free(queue); qaws_internal_free(set); qaws_internal_free(S.off); qaws_internal_free(S.A); qaws_internal_free(S.B); qaws_internal_free(S.C);
		if (!ok)
			goto done;
	}
	/* smoothed: each cell its neighbours' frames matched to its own, averaged, made a rotation, bent to its normals */
	for (s = 0; s < sweeps; s++)
		for (i = 0; i < m->nc; i++)
		{
			unsigned int cc = s & 1 ? m->nc - 1 - i : i, f, n = 0;
			double acc[9];
			memset(acc, 0, sizeof acc);
			for (f = 0; f < 4; f++)
			{
				unsigned int o = m->fo[m->f0[cc] + f];
				if (o == HM_NONE)
					continue;
				ff_add_matched(R + 9 * (size_t)cc, R + 9 * (size_t)(o / 4), acc);   /* (tetrahedra: facet o of cell o / 4) */
				n++;
			}
			if (!n)
				continue;
			ff_polar(acc);
			ff_constrain(acc, nrm + 6 * (size_t)cc, ncon[cc]);
			memcpy(R + 9 * (size_t)cc, acc, sizeof acc);
		}
	/* singular edges: round each inner edge, the matchings composed */
	if (out_singular)
	{
		memset(out_singular, 0, m->nc);
		for (c = 0; c < m->nc; c++)
		{
			unsigned int const* v = m->cv + m->c0[c];
			unsigned int a, b;
			for (a = 0; a < 4; a++)
				for (b = a + 1; b < 4; b++)
				{
					/* the cells round edge {v[a], v[b]}: across the facet opposite one of the other two vertices, and on */
					unsigned int ea = v[a], eb = v[b], cur = c, prev = HM_NONE, guard, lowest = c, ring[64], nr = 0;
					int T[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 }, inner = 1;
					for (guard = 0; guard < 64; guard++)
					{
						unsigned int f, nx = HM_NONE;
						ring[nr++] = cur;
						for (f = 0; f < 4; f++)
						{
							unsigned int w = m->cv[m->c0[cur] + f], o;
							if (w == ea || w == eb)
								continue;
							o = m->fo[m->f0[cur] + f];
							if (o == HM_NONE)
							{
								inner = 0;
								break;
							}
							if (o / 4 != prev && nx == HM_NONE)
								nx = o / 4;
						}
						if (!inner || nx == HM_NONE)
						{
							inner = 0;
							break;
						}
						/* composed: R_cur ~ R_nx Q, so round the cycle R_c ~ R_c (Q_last ... Q_first); T = Q T */
						{
							int P[3], U[9], r, q, z;
							ff_match(R + 9 * (size_t)cur, R + 9 * (size_t)nx, P);
							for (r = 0; r < 3; r++)
								for (q = 0; q < 3; q++)
								{
									int sum = 0;
									/* Q's column z holds its sign at row |P[z]| - 1 */
									for (z = 0; z < 3; z++)
										if ((P[z] < 0 ? -P[z] : P[z]) - 1 == r)
											sum += (P[z] < 0 ? -1 : 1) * T[3 * q + z];
									U[3 * q + r] = sum;
								}
							memcpy(T, U, sizeof T);
						}
						lowest = nx < lowest ? nx : lowest;
						prev = cur;
						cur = nx;
						if (cur == c)
							break;
					}
					if (!inner || cur != c || lowest != c)
						continue;   /* (on the border, or seen from its lowest cell) */
					if (!(T[0] == 1 && T[4] == 1 && T[8] == 1))
						for (k = 0; k < nr; k++)
						{
							/* the bit of the edge in each cell round it: edges 0-1 0-2 0-3 1-2 1-3 2-3 */
							unsigned int const* w = m->cv + m->c0[ring[k]];
							unsigned int x0, x1, e = 0;
							for (x0 = 0; x0 < 4; x0++)
								for (x1 = x0 + 1; x1 < 4; x1++, e++)
									if ((w[x0] == ea && w[x1] == eb) || (w[x0] == eb && w[x1] == ea))
										out_singular[ring[k]] |= (unsigned char)(1u << e);
						}
				}
		}
	}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(nrm);
	qaws_internal_free(ncon);
	return st;
}
