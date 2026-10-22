/*
 * Curved volume elements (qaws_volume_curve.h): quadratic or cubic nodes
 * numbered from the cells' edges and faces, the border's projected onto a
 * mesh reference, the displacement spread inside, elements relaxed toward
 * straight until their Jacobians, bounded in the Bernstein basis, pass.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_volume_curve.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define VCU_ROUNDS 12u   /* halvings before an element is made straight */
#define VCU_MAXN 64u     /* nodes per element at most (the tricubic hexahedron) */

static unsigned char const vcu_tet_edge[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
static unsigned char const vcu_tet_face[4][3] = { { 0, 1, 3 }, { 1, 2, 3 }, { 0, 2, 3 }, { 0, 1, 2 } };
static unsigned char const vcu_hex_edge[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 },
	{ 2, 6 }, { 3, 7 } };
static unsigned char const vcu_hex_face[6][4] = { { 0, 3, 7, 4 }, { 1, 2, 6, 5 }, { 0, 1, 5, 4 }, { 3, 2, 6, 7 }, { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };

/* an element's nodes of degree p: per slot its lattice point (a hexahedron's (i, j, k) in [0, p]^3, a tetrahedron's
   barycentric counts summing to p) and what it lies on (0 corner, 1 edge, 2 face, 3 inside; the edge's or face's local
   index; along an edge its step from the edge's first corner, on a face its (a, b) in the face's frame: a from corner 0
   toward corner 1, b from corner 0 toward the last) */
typedef struct vcu_el
{
	unsigned int p, nn, hex;
	unsigned char lat[VCU_MAXN][4], kind[VCU_MAXN], ent[VCU_MAXN], a[VCU_MAXN], b[VCU_MAXN];
	unsigned char at[VCU_MAXN];   /* a lattice point's slot, by i + (p + 1) (j + (p + 1) k) */
} vcu_el;

static unsigned int vcu_lat_index(vcu_el const* e, unsigned char const* l)
{
	return l[0] + (e->p + 1) * (l[1] + (e->p + 1) * l[2]);
}

static void vcu_hex_corner(unsigned int q, unsigned int p, unsigned char* l)
{
	l[0] = (unsigned char)((q == 1 || q == 2 || q == 5 || q == 6) ? p : 0);
	l[1] = (unsigned char)((q == 2 || q == 3 || q == 6 || q == 7) ? p : 0);
	l[2] = (unsigned char)(q >= 4 ? p : 0);
	l[3] = 0;
}

static void vcu_tet_corner(unsigned int q, unsigned int p, unsigned char* l)
{
	unsigned int k;
	for (k = 0; k < 4; k++) l[k] = (unsigned char)(k == q ? p : 0);
}

/* the lattice point a / p of the way from corner c0 toward c1 and b / p toward c3 (a hexahedron's face or edge), or the
   barycentric mix (a tetrahedron's) */
static void vcu_mix(unsigned char* l, unsigned char const* l0, unsigned char const* l1, unsigned char const* l3, unsigned int a, unsigned int b, unsigned int p)
{
	unsigned int k;
	for (k = 0; k < 4; k++) l[k] = (unsigned char)(((int)l0[k] * (int)p + ((int)l1[k] - (int)l0[k]) * (int)a + ((int)l3[k] - (int)l0[k]) * (int)b) / (int)p);
}

static void vcu_slot(vcu_el* e, unsigned char const* l, unsigned int kind, unsigned int ent, unsigned int a, unsigned int b)
{
	unsigned int n = e->nn++;
	memcpy(e->lat[n], l, 4);
	e->kind[n] = (unsigned char)kind;
	e->ent[n] = (unsigned char)ent;
	e->a[n] = (unsigned char)a;
	e->b[n] = (unsigned char)b;
	e->at[vcu_lat_index(e, l)] = (unsigned char)n;
}

/* the layout of a hexahedron or a tetrahedron of degree p (2 or 3): corners, edges' nodes, faces', the inside's */
static void vcu_layout(vcu_el* e, int hex, unsigned int p)
{
	unsigned char c[8][4], l[4];
	unsigned int q, m, a, b, i, j, k;
	memset(e, 0, sizeof *e);
	e->p = p;
	e->hex = (unsigned int)hex;
	for (q = 0; q < (hex ? 8u : 4u); q++)
	{
		if (hex)
			vcu_hex_corner(q, p, c[q]);
		else
			vcu_tet_corner(q, p, c[q]);
		vcu_slot(e, c[q], 0, q, 0, 0);
	}
	for (q = 0; q < (hex ? 12u : 6u); q++)
		for (m = 1; m < p; m++)
		{
			unsigned char const* c0 = c[hex ? vcu_hex_edge[q][0] : vcu_tet_edge[q][0]];
			vcu_mix(l, c0, c[hex ? vcu_hex_edge[q][1] : vcu_tet_edge[q][1]], c0, m, 0, p);
			vcu_slot(e, l, 1, q, m, 0);
		}
	if (hex)
	{
		for (q = 0; q < 6; q++)
			for (b = 1; b < p; b++)
				for (a = 1; a < p; a++)
				{
					vcu_mix(l, c[vcu_hex_face[q][0]], c[vcu_hex_face[q][1]], c[vcu_hex_face[q][3]], a, b, p);
					vcu_slot(e, l, 2, q, a, b);
				}
		for (k = 1; k < p; k++)
			for (j = 1; j < p; j++)
				for (i = 1; i < p; i++)
				{
					l[0] = (unsigned char)i;
					l[1] = (unsigned char)j;
					l[2] = (unsigned char)k;
					l[3] = 0;
					vcu_slot(e, l, 3, 0, 0, 0);
				}
	}
	else if (p == 3)
		for (q = 0; q < 4; q++)
		{
			/* (a cubic tetrahedron's face: its centre, one node) */
			for (k = 0; k < 4; k++) l[k] = (unsigned char)(c[vcu_tet_face[q][0]][k] / 3 + c[vcu_tet_face[q][1]][k] / 3 + c[vcu_tet_face[q][2]][k] / 3);
			vcu_slot(e, l, 2, q, 1, 1);
		}
}

static double vcu_det(double const* a, double const* b, double const* c)
{
	return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

/* the 1D Lagrange basis on p + 1 equal steps of [-1, 1] at s: values and derivatives */
static void vcu_lag(unsigned int p, double s, double* v, double* dv)
{
	unsigned int m, r, u;
	for (m = 0; m <= p; m++)
	{
		double tm = -1 + 2.0 * m / p, val = 1, der = 0;
		for (r = 0; r <= p; r++)
			if (r != m)
			{
				double tr = -1 + 2.0 * r / p, term = 1 / (tm - tr);
				val *= (s - tr) / (tm - tr);
				for (u = 0; u <= p; u++)
					if (u != m && u != r)
					{
						double tu = -1 + 2.0 * u / p;
						term *= (s - tu) / (tm - tu);
					}
				der += term;
			}
		v[m] = val;
		dv[m] = der;
	}
}

/* the 1D factor of a simplex Lagrange basis: prod over m < a of (p x - m) / (m + 1), and its derivative */
static void vcu_simplex_factor(unsigned int p, unsigned int a, double x, double* v, double* dv)
{
	unsigned int m, r;
	double val = 1, der = 0;
	for (m = 0; m < a; m++) val *= (p * x - m) / (m + 1);
	for (r = 0; r < a; r++)
	{
		double term = (double)p / (r + 1);
		for (m = 0; m < a; m++)
			if (m != r)
				term *= (p * x - m) / (m + 1);
		der += term;
	}
	*v = val;
	*dv = der;
}

/* an element's Jacobian determinant: a hexahedron's at (s, t, u) in [-1, 1]^3, a tetrahedron's at barycentric pt */
static double vcu_el_det(vcu_el const* e, double const (*x)[3], double const* pt)
{
	double col[3][3];
	unsigned int n, d, k;
	memset(col, 0, sizeof col);
	if (e->hex)
	{
		double a[4], da[4], b[4], db[4], c[4], dc[4];
		vcu_lag(e->p, pt[0], a, da);
		vcu_lag(e->p, pt[1], b, db);
		vcu_lag(e->p, pt[2], c, dc);
		for (n = 0; n < e->nn; n++)
		{
			unsigned int i = e->lat[n][0], j = e->lat[n][1], kk = e->lat[n][2];
			double w0 = da[i] * b[j] * c[kk], w1 = a[i] * db[j] * c[kk], w2 = a[i] * b[j] * dc[kk];
			for (d = 0; d < 3; d++)
			{
				col[0][d] += w0 * x[n][d];
				col[1][d] += w1 * x[n][d];
				col[2][d] += w2 * x[n][d];
			}
		}
	}
	else
		for (n = 0; n < e->nn; n++)
		{
			double f[4], df[4], g[4];
			for (k = 0; k < 4; k++) vcu_simplex_factor(e->p, e->lat[n][k], pt[k], &f[k], &df[k]);
			/* (d phi / d lambda_k, then the parameters' derivatives: lambda_k+1 against lambda_0) */
			for (k = 0; k < 4; k++)
			{
				unsigned int r;
				g[k] = df[k];
				for (r = 0; r < 4; r++)
					if (r != k)
						g[k] *= f[r];
			}
			for (k = 0; k < 3; k++)
				for (d = 0; d < 3; d++) col[k][d] += x[n][d] * (g[k + 1] - g[0]);
		}
	return vcu_det(col[0], col[1], col[2]);
}

/* the Bernstein machinery for degree p elements: a tetrahedron's determinant of degree q = 3 (p - 1) collocated at its
   lattice points of degree q, the inverse of the collocation in tinv; a hexahedron's of degree 3 p - 1 per axis, the 1D
   collocation's inverse in hinv */
typedef struct vcu_bern
{
	unsigned int qt, nq, qh;
	unsigned char (*tl)[4];
	double* tinv;
	double hinv[10][10];
} vcu_bern;

/* Gauss-Jordan: n x n A (by rows) inverted into I; 0 when singular */
static int vcu_invert(double* A, double* I, unsigned int n)
{
	unsigned int r, q, i;
	for (r = 0; r < n; r++)
		for (q = 0; q < n; q++) I[r * n + q] = r == q;
	for (q = 0; q < n; q++)
	{
		unsigned int piv = q;
		double f;
		for (r = q + 1; r < n; r++)
			if (fabs(A[r * n + q]) > fabs(A[piv * n + q]))
				piv = r;
		if (!(fabs(A[piv * n + q]) > 0))
			return 0;
		for (i = 0; i < n; i++)
		{
			f = A[q * n + i]; A[q * n + i] = A[piv * n + i]; A[piv * n + i] = f;
			f = I[q * n + i]; I[q * n + i] = I[piv * n + i]; I[piv * n + i] = f;
		}
		f = A[q * n + q];
		for (i = 0; i < n; i++) { A[q * n + i] /= f; I[q * n + i] /= f; }
		for (r = 0; r < n; r++)
			if (r != q && A[r * n + q] != 0)
			{
				f = A[r * n + q];
				for (i = 0; i < n; i++) { A[r * n + i] -= f * A[q * n + i]; I[r * n + i] -= f * I[q * n + i]; }
			}
	}
	return 1;
}

static double vcu_fact(unsigned int n)
{
	double f = 1;
	while (n > 1) f *= n--;
	return f;
}

static void vcu_bern_free(vcu_bern* B)
{
	qaws_internal_free(B->tl);
	qaws_internal_free(B->tinv);
	memset(B, 0, sizeof *B);
}

/* 0 on allocation failure */
static int vcu_bern_init(vcu_bern* B, unsigned int p)
{
	unsigned int q = 3 * (p - 1), n = 0, i, j, k, r, c;
	double *A = NULL, H[100], HA[100];
	memset(B, 0, sizeof *B);
	B->qt = q;
	B->nq = (q + 1) * (q + 2) * (q + 3) / 6;
	B->qh = 3 * p - 1;
	B->tl = (unsigned char(*)[4])qaws_internal_malloc(4 * (size_t)B->nq);
	B->tinv = (double*)qaws_internal_malloc(sizeof(double) * (size_t)B->nq * B->nq);
	A = (double*)qaws_internal_scratch_malloc(sizeof(double) * (size_t)B->nq * B->nq);
	if (!B->tl || !B->tinv || !A)
	{
		qaws_internal_free(A);
		vcu_bern_free(B);
		return 0;
	}
	for (i = 0; i <= q; i++)
		for (j = 0; i + j <= q; j++)
			for (k = 0; i + j + k <= q; k++)
			{
				B->tl[n][0] = (unsigned char)i;
				B->tl[n][1] = (unsigned char)j;
				B->tl[n][2] = (unsigned char)k;
				B->tl[n][3] = (unsigned char)(q - i - j - k);
				n++;
			}
	for (r = 0; r < B->nq; r++)
		for (c = 0; c < B->nq; c++)
		{
			double v = vcu_fact(q);
			for (k = 0; k < 4; k++) v *= pow(B->tl[r][k] / (double)q, (double)B->tl[c][k]) / vcu_fact(B->tl[c][k]);
			A[r * B->nq + c] = v;
		}
	k = (unsigned int)vcu_invert(A, B->tinv, B->nq);
	qaws_internal_free(A);
	for (r = 0; r <= B->qh; r++)
		for (c = 0; c <= B->qh; c++)
		{
			double x = r / (double)B->qh;
			H[r * (B->qh + 1) + c] = vcu_fact(B->qh) / (vcu_fact(c) * vcu_fact(B->qh - c)) * pow(x, (double)c) * pow(1 - x, (double)(B->qh - c));
		}
	k &= (unsigned int)vcu_invert(H, HA, B->qh + 1);
	for (r = 0; r <= B->qh; r++)
		for (c = 0; c <= B->qh; c++) B->hinv[r][c] = HA[r * (B->qh + 1) + c];
	(void)k;
	return 1;
}

/* A tetrahedron's determinant (degree 3 (p - 1) in the barycentrics) in the Bernstein basis of that degree: its values at
   the lattice points solved for the coefficients; their least bounds the determinant from below over the whole element
   (Johnen, Remacle, Geuzaine 2013) */
static double vcu_tet_bound(vcu_bern const* B, vcu_el const* e, double const (*x)[3])
{
	double v[84], lo = HUGE_VAL;
	unsigned int r, c, k;
	for (r = 0; r < B->nq; r++)
	{
		double lam[4];
		for (k = 0; k < 4; k++) lam[k] = B->tl[r][k] / (double)B->qt;
		v[r] = vcu_el_det(e, x, lam);
	}
	for (r = 0; r < B->nq; r++)
	{
		double s = 0;
		for (c = 0; c < B->nq; c++) s += B->tinv[r * B->nq + c] * v[c];
		lo = fmin(lo, s);
	}
	return lo;
}

/* A hexahedron's determinant, of degree 3 p - 1 in each of s, t, u: in the tensor Bernstein basis of that degree its
   coefficients from its values at the grid points (one 1D solve along each axis). The straight hexahedron's in the same
   basis: where all its coefficients s_i are positive, the curved one's c_i >= r s_i for r the least c_i / s_i, so the
   ratio of determinants is >= r everywhere. */
static double vcu_hex_bound(vcu_bern const* B, vcu_el const* e, double const (*xc)[3], double const (*xs)[3], double* c, double* s, double const* lo,
	double const* hi)
{
	unsigned int m = B->qh + 1, n3 = m * m * m, i, j, k, p, q, ax;
	double r = HUGE_VAL, cmin = HUGE_VAL;
	int straight = 1;
	for (k = 0; k < m; k++)
		for (j = 0; j < m; j++)
			for (i = 0; i < m; i++)
			{
				double pt[3];
				pt[0] = lo[0] + (hi[0] - lo[0]) * i / B->qh;
				pt[1] = lo[1] + (hi[1] - lo[1]) * j / B->qh;
				pt[2] = lo[2] + (hi[2] - lo[2]) * k / B->qh;
				c[i + m * (j + m * k)] = vcu_el_det(e, xc, pt);
				s[i + m * (j + m * k)] = vcu_el_det(e, xs, pt);
			}
	for (ax = 0; ax < 3; ax++)
	{
		unsigned int st = ax == 0 ? 1 : ax == 1 ? m : m * m, line;
		for (line = 0; line < m * m; line++)
		{
			/* (the line's first index: the other two axes' digits) */
			unsigned int b0 = ax == 0 ? m * line : ax == 1 ? line % m + m * m * (line / m) : line;
			double fc[10], fs[10];
			for (p = 0; p < m; p++) { fc[p] = c[b0 + st * p]; fs[p] = s[b0 + st * p]; }
			for (q = 0; q < m; q++)
			{
				double vc = 0, vs = 0;
				for (p = 0; p < m; p++) { vc += B->hinv[q][p] * fc[p]; vs += B->hinv[q][p] * fs[p]; }
				c[b0 + st * q] = vc;
				s[b0 + st * q] = vs;
			}
		}
	}
	for (i = 0; i < n3; i++)
	{
		cmin = fmin(cmin, c[i]);
		if (s[i] > 0)
			r = fmin(r, c[i] / s[i]);
		else
			straight = 0;
	}
	return straight ? r : (cmin > 0 ? 1 : -1);
}

/* the bound over box lo hi, its eight halves' where it is under target (depth more times at most): the bound tightens
   as the boxes shrink (Johnen, Remacle, Geuzaine 2013) */
static double vcu_hex_bound_split(vcu_bern const* B, vcu_el const* e, double const (*xc)[3], double const (*xs)[3], double* work, double const* lo,
	double const* hi, unsigned int depth, double target)
{
	double r = vcu_hex_bound(B, e, xc, xs, work, work + 1000, lo, hi), mid[3], l2[3], h2[3];
	unsigned int o, k;
	if (r >= target || depth == 0)
		return r;
	for (k = 0; k < 3; k++) mid[k] = 0.5 * (lo[k] + hi[k]);
	r = HUGE_VAL;
	for (o = 0; o < 8; o++)
	{
		for (k = 0; k < 3; k++)
		{
			l2[k] = (o >> k) & 1 ? mid[k] : lo[k];
			h2[k] = (o >> k) & 1 ? hi[k] : mid[k];
		}
		r = fmin(r, vcu_hex_bound_split(B, e, xc, xs, work, l2, h2, depth - 1, target));
	}
	return r;
}

/* an element's smallest Jacobian over its straight self's, bounded below over the whole element (nodes en, positions X straight + D) */
static double vcu_ratio(vcu_bern const* B, vcu_el const* e, unsigned int const* en, double const* X, double const* D, double* work, double target)
{
	double xc[VCU_MAXN][3], xs[VCU_MAXN][3];
	unsigned int n, d;
	for (n = 0; n < e->nn; n++)
		for (d = 0; d < 3; d++)
		{
			xs[n][d] = X[3 * (size_t)en[n] + d];
			xc[n][d] = xs[n][d] + D[3 * (size_t)en[n] + d];
		}
	if (!e->hex)
	{
		/* (the straight tetrahedron's determinant is constant: the bound over it) */
		static double const mid[4] = { 0.25, 0.25, 0.25, 0.25 };
		double ds = vcu_el_det(e, (double const(*)[3])xs, mid), dc = vcu_tet_bound(B, e, (double const(*)[3])xc);
		return ds > 0 ? dc / ds : (dc > 0 ? 1 : -1);
	}
	{
		static double const lo[3] = { -1, -1, -1 }, hi[3] = { 1, 1, 1 };
		return vcu_hex_bound_split(B, e, (double const(*)[3])xc, (double const(*)[3])xs, work, lo, hi, 3, target);
	}
}

/* a hexahedron's inside nodes' displacements from its boundary's: Gordon and Hall's transfinite blend, the boolean sum
   of the linear blends along s, t and u (each term the boundary lattice nodes with the chosen axes moved to 0 or p,
   weighted by 1 - x or x; sign alternating with how many axes move), so they follow its curved faces (en its nodes) */
static void vcu_gordon_hall(vcu_el const* e, unsigned int const* en, double* D)
{
	unsigned int s, d, p = e->p;
	for (s = 0; s < e->nn; s++)
	{
		double x[3], y[3] = { 0, 0, 0 };
		unsigned int mask, sides, ax;
		if (e->kind[s] != 3)
			continue;
		for (ax = 0; ax < 3; ax++) x[ax] = e->lat[s][ax] / (double)p;
		for (mask = 1; mask < 8; mask++)
			for (sides = 0; sides < 8; sides++)
			{
				unsigned char l[4];
				double w = 1;
				int odd = 0;
				unsigned int g;
				if (sides & ~mask)
					continue;   /* (the sides of the axes that move only) */
				for (ax = 0; ax < 4; ax++) l[ax] = ax < 3 ? e->lat[s][ax] : 0;
				for (ax = 0; ax < 3; ax++)
					if (mask & (1u << ax))
					{
						int hi = (sides >> ax) & 1;
						w *= hi ? x[ax] : 1 - x[ax];
						l[ax] = (unsigned char)(hi ? p : 0);
						odd = !odd;
					}
				g = en[e->at[vcu_lat_index(e, l)]];
				for (d = 0; d < 3; d++) y[d] += (odd ? w : -w) * D[3 * (size_t)g + d];
			}
		for (d = 0; d < 3; d++) D[3 * (size_t)en[s] + d] = y[d];
	}
}

/* a face node's displacement blended from its face's edges' (slot s of an element with nodes en): a hexahedron's face by
   Coons (the edges at the node's (u, v), less the corners' bilinear blend), a cubic tetrahedron's centre by its six edge
   nodes' sum over 4 less its corners' over 6 (both exact for straight faces, frame-free: the cells sharing a face agree) */
static void vcu_coons(double* out, vcu_el const* e, unsigned int const* en, unsigned int s, double const* D)
{
	unsigned int p = e->p, k, d, g;
	unsigned char l[4];
	for (d = 0; d < 3; d++) out[d] = 0;
	if (!e->hex)
	{
		unsigned char const* f3 = vcu_tet_face[e->ent[s]];
		for (k = 0; k < 3; k++)
		{
			unsigned char c0[4], c1[4];
			unsigned int q2;
			vcu_tet_corner(f3[k], p, c0);
			vcu_tet_corner(f3[(k + 1) % 3], p, c1);
			for (q2 = 1; q2 < p; q2++)
			{
				vcu_mix(l, c0, c1, c0, q2, 0, p);
				g = en[e->at[vcu_lat_index(e, l)]];
				for (d = 0; d < 3; d++) out[d] += D[3 * (size_t)g + d] / 4;
			}
			for (d = 0; d < 3; d++) out[d] -= D[3 * (size_t)en[f3[k]] + d] / 6;
		}
		return;
	}
	{
		unsigned char const* f4 = vcu_hex_face[e->ent[s]];
		unsigned char cl[4][4];
		unsigned int a = e->a[s], b = e->b[s];
		double u = a / (double)p, v = b / (double)p, wq[4];
		for (k = 0; k < 4; k++) vcu_hex_corner(f4[k], p, cl[k]);
		/* bottom (c0 c1) at u, top (c3 c2) at u, left (c0 c3) at v, right (c1 c2) at v */
		vcu_mix(l, cl[0], cl[1], cl[0], a, 0, p);
		g = en[e->at[vcu_lat_index(e, l)]];
		for (d = 0; d < 3; d++) out[d] += (1 - v) * D[3 * (size_t)g + d];
		vcu_mix(l, cl[3], cl[2], cl[3], a, 0, p);
		g = en[e->at[vcu_lat_index(e, l)]];
		for (d = 0; d < 3; d++) out[d] += v * D[3 * (size_t)g + d];
		vcu_mix(l, cl[0], cl[3], cl[0], b, 0, p);
		g = en[e->at[vcu_lat_index(e, l)]];
		for (d = 0; d < 3; d++) out[d] += (1 - u) * D[3 * (size_t)g + d];
		vcu_mix(l, cl[1], cl[2], cl[1], b, 0, p);
		g = en[e->at[vcu_lat_index(e, l)]];
		for (d = 0; d < 3; d++) out[d] += u * D[3 * (size_t)g + d];
		wq[0] = (1 - u) * (1 - v);
		wq[1] = u * (1 - v);
		wq[2] = u * v;
		wq[3] = (1 - u) * v;
		for (k = 0; k < 4; k++)
			for (d = 0; d < 3; d++) out[d] -= wq[k] * D[3 * (size_t)en[f4[k]] + d];
	}
}

static int vcu_kinds(qaws_volume_mesh const* mesh, unsigned int* nt, unsigned int* nh)
{
	unsigned int c;
	*nt = *nh = 0;
	for (c = 0; c < mesh->nc; c++)
	{
		if (mesh->type[c] == QAWS_CELL_TETRAHEDRON) (*nt)++;
		else if (mesh->type[c] == QAWS_CELL_HEXAHEDRON) (*nh)++;
		else return 0;
	}
	return 1;
}

typedef struct vcu_key
{
	unsigned int k[4], c, f;   /* sorted vertices (an edge's last two HM_NONE, a triangle's last); its cell and local entity */
} vcu_key;

static int vcu_key_cmp(void const* p, void const* q)
{
	vcu_key const *a = (vcu_key const*)p, *b = (vcu_key const*)q;
	unsigned int i;
	for (i = 0; i < 4; i++)
		if (a->k[i] != b->k[i])
			return a->k[i] < b->k[i] ? -1 : 1;
	return 0;
}

static void vcu_sort4(unsigned int* s, unsigned int n)
{
	unsigned int i, j;
	for (i = 1; i < n; i++)
		for (j = i; j > 0 && s[j - 1] > s[j]; j--)
		{
			unsigned int t = s[j];
			s[j] = s[j - 1];
			s[j - 1] = t;
		}
}

/* the nodes of degree p numbered: per cell its element nodes (en, from eo), and the counts; keys (scratch, for the
   caller) sorted. Edges' nodes p - 1 each, from the smaller vertex; faces' (a hexahedron's (p - 1)^2 in the frame from
   its smallest vertex toward its smaller neighbour, a cubic tetrahedron's one), then each cell's inside. en NULL: only
   the counts */
static qaws_status vcu_number(unsigned int* out_nn, unsigned int* out_nen, unsigned int* en, unsigned int* eo, qaws_volume_mesh const* mesh,
	vcu_el const* lt, vcu_el const* lh)
{
	unsigned int p = lh->p, c, k, i, nk = 0, ne = 0, nfn = 0, nin = 0, nen = 0, fpq = (p - 1) * (p - 1), fpt = (p - 1) * (p - 2) / 2;
	unsigned int ipc = (p - 1) * (p - 1) * (p - 1);
	vcu_key* key = (vcu_key*)qaws_internal_scratch_malloc(sizeof(vcu_key) * (10 * (size_t)mesh->nc + 18 * (size_t)mesh->nc + 1));
	if (!key)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (c = 0; c < mesh->nc; c++)
	{
		unsigned int const* v = mesh->cv + mesh->c0[c];
		int tet = mesh->type[c] == QAWS_CELL_TETRAHEDRON;
		if (eo)
			eo[c] = nen;
		nen += tet ? lt->nn : lh->nn;
		for (k = 0; k < (tet ? 6u : 12u); k++)
		{
			key[nk].k[0] = v[tet ? vcu_tet_edge[k][0] : vcu_hex_edge[k][0]];
			key[nk].k[1] = v[tet ? vcu_tet_edge[k][1] : vcu_hex_edge[k][1]];
			key[nk].k[2] = key[nk].k[3] = HM_NONE;
			vcu_sort4(key[nk].k, 2);
			key[nk].c = c;
			key[nk++].f = k;
		}
		for (k = 0; (tet ? fpt : fpq) && k < (tet ? 4u : 6u); k++)
		{
			for (i = 0; i < (tet ? 3u : 4u); i++) key[nk].k[i] = v[tet ? vcu_tet_face[k][i] : vcu_hex_face[k][i]];
			if (tet)
				key[nk].k[3] = HM_NONE;
			vcu_sort4(key[nk].k, tet ? 3 : 4);
			key[nk].c = c;
			key[nk++].f = k;
		}
	}
	if (eo)
		eo[mesh->nc] = nen;
	qsort(key, nk, sizeof(vcu_key), vcu_key_cmp);
	/* the edges first, then the faces (two passes over the sorted keys) */
	{
		unsigned int id = 0, base = mesh->nv, pass;
		for (pass = 0; pass < 2; pass++)
		for (k = 0; k < nk; k++)
		{
			int fresh = k == 0 || vcu_key_cmp(&key[k], &key[k - 1]);
			unsigned int const* v = mesh->cv + mesh->c0[key[k].c];
			int tet = mesh->type[key[k].c] == QAWS_CELL_TETRAHEDRON;
			vcu_el const* e = tet ? lt : lh;
			unsigned int s;
			if ((key[k].k[2] == HM_NONE) != (pass == 0))
				continue;
			if (key[k].k[2] == HM_NONE)
			{
				/* an edge: its p - 1 nodes from its smaller vertex */
				unsigned int a0 = tet ? vcu_tet_edge[key[k].f][0] : vcu_hex_edge[key[k].f][0];
				if (fresh)
				{
					id = base;
					base += p - 1;
					ne++;
				}
				for (s = 0; en && s < e->nn; s++)
					if (e->kind[s] == 1 && e->ent[s] == key[k].f)
						en[eo[key[k].c] + s] = id + (v[a0] == key[k].k[0] ? e->a[s] : p - e->a[s]) - 1;
			}
			else
			{
				unsigned int cnt = key[k].k[3] == HM_NONE ? fpt : fpq;
				if (fresh)
				{
					id = base;
					base += cnt;
					nfn += cnt;
				}
				for (s = 0; en && s < e->nn; s++)
				{
					unsigned int g[4], o = 0, q, P0[2], PU[2], PV[2], w[2], sa, sb;
					static unsigned int const pc[4][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
					double det, ss, tt;
					if (e->kind[s] != 2 || e->ent[s] != key[k].f)
						continue;
					if (tet)
					{
						en[eo[key[k].c] + s] = id;
						continue;
					}
					/* the face's frame from its smallest vertex, toward its smaller neighbour: the node's (a, b) in it */
					for (q = 0; q < 4; q++)
					{
						g[q] = v[vcu_hex_face[key[k].f][q]];
						if (g[q] < g[o])
							o = q;
					}
					memcpy(P0, pc[o], sizeof P0);
					if (g[(o + 1) % 4] < g[(o + 3) % 4])
					{
						memcpy(PU, pc[(o + 1) % 4], sizeof PU);
						memcpy(PV, pc[(o + 3) % 4], sizeof PV);
					}
					else
					{
						memcpy(PU, pc[(o + 3) % 4], sizeof PU);
						memcpy(PV, pc[(o + 1) % 4], sizeof PV);
					}
					w[0] = e->a[s];
					w[1] = e->b[s];
					/* (w - p P0) = ss p (PU - P0) + tt p (PV - P0), in steps of 1 / p) */
					det = ((double)PU[0] - P0[0]) * ((double)PV[1] - P0[1]) - ((double)PU[1] - P0[1]) * ((double)PV[0] - P0[0]);
					ss = (((double)w[0] - (double)p * P0[0]) * ((double)PV[1] - P0[1]) - ((double)w[1] - (double)p * P0[1]) * ((double)PV[0] - P0[0])) / det;
					tt = ((((double)PU[0] - P0[0]) * ((double)w[1] - (double)p * P0[1])) - (((double)PU[1] - P0[1]) * ((double)w[0] - (double)p * P0[0]))) / det;
					sa = (unsigned int)(ss + 0.5);
					sb = (unsigned int)(tt + 0.5);
					en[eo[key[k].c] + s] = id + (sb - 1) * (p - 1) + (sa - 1);
				}
			}
		}
		for (c = 0; c < mesh->nc; c++)
		{
			unsigned int cnt = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? 0 : ipc, s, j = 0;
			vcu_el const* e = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? lt : lh;
			for (s = 0; en && s < e->nn; s++)
			{
				if (e->kind[s] == 0)
					en[eo[c] + s] = mesh->cv[mesh->c0[c] + e->ent[s]];
				else if (e->kind[s] == 3)
					en[eo[c] + s] = base + j++;
			}
			base += cnt;
			nin += cnt;
		}
		*out_nn = base;
	}
	*out_nen = nen;
	(void)ne;
	(void)nfn;
	(void)nin;
	qaws_internal_free(key);
	return QAWS_STATUS_OK;
}

static qaws_status vcu_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_volume_mesh const* mesh, unsigned int p)
{
	unsigned int nt, nh;
	vcu_el lt, lh;
	if (!out_node_count || !out_element_node_count || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!vcu_kinds(mesh, &nt, &nh))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	vcu_layout(&lt, 0, p);
	vcu_layout(&lh, 1, p);
	return vcu_number(out_node_count, out_element_node_count, NULL, NULL, mesh, &lt, &lh);
}

qaws_status qaws_volume_mesh_quadratic_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_volume_mesh const* mesh)
{
	return vcu_counts(out_node_count, out_element_node_count, mesh, 2);
}

qaws_status qaws_volume_mesh_cubic_node_counts(unsigned int* out_node_count, unsigned int* out_element_node_count, qaws_volume_mesh const* mesh)
{
	return vcu_counts(out_node_count, out_element_node_count, mesh, 3);
}

/* the straight point of lattice point l of a cell (corners at v) */
static void vcu_straight(double* out, qaws_volume_mesh const* mesh, unsigned int const* v, vcu_el const* e, unsigned char const* l)
{
	unsigned int d, q;
	for (d = 0; d < 3; d++) out[d] = 0;
	if (e->hex)
	{
		double s = l[0] / (double)e->p, t = l[1] / (double)e->p, u = l[2] / (double)e->p;
		for (q = 0; q < 8; q++)
		{
			unsigned char cl[4];
			double w;
			vcu_hex_corner(q, 1, cl);
			w = (cl[0] ? s : 1 - s) * (cl[1] ? t : 1 - t) * (cl[2] ? u : 1 - u);
			for (d = 0; d < 3; d++) out[d] += w * mesh->x[3 * (size_t)v[q] + d];
		}
	}
	else
		for (q = 0; q < 4; q++)
			for (d = 0; d < 3; d++) out[d] += l[q] / (double)e->p * mesh->x[3 * (size_t)v[q] + d];
}

/* the cell's local face whose corners are the facet's (loc, n of them) */
static unsigned int vcu_local_face(int tet, unsigned char const* loc, unsigned int n)
{
	unsigned int sv[4], fv[4], j, d;
	for (j = 0; j < n; j++) sv[j] = loc[j];
	vcu_sort4(sv, n);
	for (j = 0; j < (tet ? 4u : 6u); j++)
	{
		for (d = 0; d < n; d++) fv[d] = tet ? vcu_tet_face[j][d] : vcu_hex_face[j][d];
		vcu_sort4(fv, n);
		if (!memcmp(fv, sv, sizeof(unsigned int) * n))
			return j;
	}
	return HM_NONE;
}

static qaws_status vcu_nodes(qaws_scalar* out_positions, unsigned int* out_element_nodes, qaws_volume_curve_report* report,
	qaws_volume_mesh const* mesh, qaws_mesh_reference const* reference, qaws_volume_curve_desc const* desc, unsigned int p)
{
	unsigned int nt, nh, c, k, i, d, nn, nen, nb = 0, *eo = NULL, *en, *bord = NULL, *bc = NULL, *bs = NULL;
	unsigned char *isb = NULL, *flag = NULL;
	double *X = NULL, *D = NULL, *acc = NULL, *work = NULL, *arc = NULL, minj = desc && desc->min_jacobian > 0 ? (double)desc->min_jacobian : 0.1;
	unsigned int sweeps = desc && desc->smoothing ? desc->smoothing : 20, *cnt = NULL;
	qaws_volume_curve_report rep;
	vcu_el lt, lh;
	vcu_bern B;
	qaws_status st;
	memset(&rep, 0, sizeof rep);
	memset(&B, 0, sizeof B);
	if (!out_positions || !out_element_nodes || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!vcu_kinds(mesh, &nt, &nh))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	vcu_layout(&lt, 0, p);
	vcu_layout(&lh, 1, p);
	en = out_element_nodes;
	eo = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)mesh->nc + 1));
	if (!eo)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = vcu_number(&nn, &nen, en, eo, mesh, &lt, &lh);
	if (st != QAWS_STATUS_OK)
		goto done;
	isb = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nn + 1, 1);
	flag = (unsigned char*)qaws_internal_scratch_zalloc((size_t)mesh->nc + 1, 1);
	X = (double*)qaws_internal_scratch_zalloc(3 * (size_t)nn + 1, sizeof(double));
	D = (double*)qaws_internal_scratch_zalloc(3 * (size_t)nn + 1, sizeof(double));
	acc = (double*)qaws_internal_scratch_zalloc(3 * (size_t)nn + 1, sizeof(double));
	cnt = (unsigned int*)qaws_internal_scratch_zalloc((size_t)nn + 1, sizeof(unsigned int));
	bord = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nn + 1));
	bc = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nn + 1));   /* a border face node's cell and slot */
	bs = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nn + 1));
	work = (double*)qaws_internal_scratch_malloc(sizeof(double) * 2000);
	if (!isb || !flag || !X || !D || !acc || !cnt || !bord || !bc || !bs || !work || !vcu_bern_init(&B, p))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}

	/* the straight points */
	for (c = 0; c < mesh->nc; c++)
	{
		vcu_el const* e = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? &lt : &lh;
		for (k = 0; k < e->nn; k++) vcu_straight(X + 3 * (size_t)en[eo[c] + k], mesh, mesh->cv + mesh->c0[c], e, e->lat[k]);
	}

	/* the border's nodes: on facets with no cell across (1 an edge's, 2 a face's) */
	for (c = 0; c < mesh->nc; c++)
	{
		unsigned int nfa = mesh->f0[c + 1] - mesh->f0[c], f;
		int tet = mesh->type[c] == QAWS_CELL_TETRAHEDRON;
		vcu_el const* e = tet ? &lt : &lh;
		for (f = 0; f < nfa; f++)
		{
			unsigned int n, j, lf, s;
			unsigned char const* loc;
			if (mesh->fo[mesh->f0[c] + f] != HM_NONE)
				continue;
			loc = vm_facet(mesh->type[c], f, &n);
			lf = vcu_local_face(tet, loc, n);
			for (s = 0; s < e->nn; s++)
			{
				unsigned int g = en[eo[c] + s], on = 0;
				if (e->kind[s] == 1)
				{
					unsigned int a0 = tet ? vcu_tet_edge[e->ent[s]][0] : vcu_hex_edge[e->ent[s]][0], a1 = tet ? vcu_tet_edge[e->ent[s]][1] : vcu_hex_edge[e->ent[s]][1], h0 = 0, h1 = 0;
					for (j = 0; j < n; j++)
					{
						h0 |= loc[j] == a0;
						h1 |= loc[j] == a1;
					}
					on = h0 && h1;
				}
				else if (e->kind[s] == 2 && e->ent[s] == lf)
					on = 2;
				if (!on || isb[g])
					continue;
				isb[g] = (unsigned char)on;
				bc[nb] = c;
				bs[nb] = s;
				bord[nb++] = g;
			}
		}
	}
	rep.border_nodes = reference ? nb : 0;

	/* edge nodes, then face nodes (from their Coons points) onto the reference */
	if (reference && nb)
	{
		qaws_scalar* pts = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 4 * (size_t)nb);
		qaws_mesh_location* loc = (qaws_mesh_location*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_location) * (size_t)nb);
		unsigned int pass;
		if (!pts || !loc)
		{
			qaws_internal_free(pts);
			qaws_internal_free(loc);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		{
			/* the border's edge nodes at equal lengths along their curved edges (projecting the straight points bunches
			   them toward an arc's middle, and high degrees overshoot on uneven nodes): each edge's segment sampled, the
			   samples projected, the node's start at its fraction of the projected polyline's length */
			unsigned int K = 16, ne1 = 0, j, i2, m = 0;
			qaws_scalar* sp;
			qaws_mesh_location* sl;
			for (j = 0; j < nb; j++) ne1 += isb[bord[j]] == 1;
			sp = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 3 * ((size_t)ne1 * (K + 1) + 1));
			sl = (qaws_mesh_location*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_location) * ((size_t)ne1 * (K + 1) + 1));
			arc = (double*)qaws_internal_scratch_malloc(sizeof(double) * 3 * ((size_t)nb + 1));
			if (!sp || !sl || !arc)
				st = QAWS_STATUS_ALLOCATION_FAILURE;
			for (j = 0; st == QAWS_STATUS_OK && j < nb; j++)
				if (isb[bord[j]] == 1)
				{
					vcu_el const* e = mesh->type[bc[j]] == QAWS_CELL_TETRAHEDRON ? &lt : &lh;
					unsigned int ed = e->ent[bs[j]], a0 = e->hex ? vcu_hex_edge[ed][0] : vcu_tet_edge[ed][0], a1 = e->hex ? vcu_hex_edge[ed][1] : vcu_tet_edge[ed][1];
					double const *A = X + 3 * (size_t)en[eo[bc[j]] + a0], *Bp = X + 3 * (size_t)en[eo[bc[j]] + a1];
					for (i2 = 0; i2 <= K; i2++, m++)
						for (d = 0; d < 3; d++) sp[3 * (size_t)m + d] = (qaws_scalar)(A[d] + (Bp[d] - A[d]) * i2 / K);
				}
			if (st == QAWS_STATUS_OK && m)
				st = qaws_mesh_reference_project(sl, NULL, reference, sp, NULL, m);
			if (st == QAWS_STATUS_OK && m)
				st = qaws_mesh_reference_evaluate(sp, NULL, NULL, reference, sl, m);
			for (j = 0, m = 0; st == QAWS_STATUS_OK && j < nb; j++)
				if (isb[bord[j]] == 1)
				{
					vcu_el const* e = mesh->type[bc[j]] == QAWS_CELL_TETRAHEDRON ? &lt : &lh;
					double len[17], want, f;
					qaws_scalar const* q = sp + 3 * (size_t)m;
					len[0] = 0;
					for (i2 = 1; i2 <= K; i2++)
					{
						double l2 = 0;
						for (d = 0; d < 3; d++) l2 += ((double)q[3 * i2 + d] - q[3 * (i2 - 1) + d]) * ((double)q[3 * i2 + d] - q[3 * (i2 - 1) + d]);
						len[i2] = len[i2 - 1] + sqrt(l2);
					}
					want = len[K] * e->a[bs[j]] / e->p;
					for (i2 = 1; i2 < K && len[i2] < want; i2++)
						;
					f = len[i2] > len[i2 - 1] ? (want - len[i2 - 1]) / (len[i2] - len[i2 - 1]) : 0;
					for (d = 0; d < 3; d++) arc[3 * (size_t)j + d] = (1 - f) * q[3 * (i2 - 1) + d] + f * q[3 * i2 + d];
					m += K + 1;
				}
			qaws_internal_free(sp);
			qaws_internal_free(sl);
		}
		for (pass = 1; pass <= 2 && st == QAWS_STATUS_OK; pass++)
		{
			unsigned int m = 0, j;
			for (j = 0; j < nb; j++)
			{
				unsigned int s = bord[j];
				double pnt[3];
				if (isb[s] != pass)
					continue;
				for (d = 0; d < 3; d++) pnt[d] = pass == 1 ? arc[3 * (size_t)j + d] : X[3 * (size_t)s + d];
				if (pass == 2)
				{
					/* the straight point plus the face's Coons blend of its curved edges' displacements */
					double dc[3];
					vcu_coons(dc, mesh->type[bc[j]] == QAWS_CELL_TETRAHEDRON ? &lt : &lh, en + eo[bc[j]], bs[j], D);
					for (d = 0; d < 3; d++) pnt[d] += dc[d];
				}
				for (d = 0; d < 3; d++) pts[3 * (size_t)m + d] = (qaws_scalar)pnt[d];
				m++;
			}
			if (!m)
				continue;
			st = qaws_mesh_reference_project(loc, NULL, reference, pts, NULL, m);
			if (st == QAWS_STATUS_OK)
				st = qaws_mesh_reference_evaluate(pts, NULL, NULL, reference, loc, m);
			for (j = 0, m = 0; st == QAWS_STATUS_OK && j < nb; j++)
			{
				unsigned int s = bord[j];
				if (isb[s] != pass)
					continue;
				for (d = 0; d < 3; d++) D[3 * (size_t)s + d] = (double)pts[3 * (size_t)m + d] - X[3 * (size_t)s + d];
				m++;
			}
		}
		/* the straight nodes' distances: the sagitta */
		if (st == QAWS_STATUS_OK)
		{
			unsigned int j;
			qaws_scalar* dist = pts + 3 * (size_t)nb;
			for (j = 0; j < nb; j++)
				for (d = 0; d < 3; d++) pts[3 * (size_t)j + d] = (qaws_scalar)X[3 * (size_t)bord[j] + d];
			st = qaws_mesh_reference_project(loc, dist, reference, pts, NULL, nb);
			for (j = 0; st == QAWS_STATUS_OK && j < nb; j++)
				rep.max_sagitta = (qaws_scalar)fmax((double)rep.max_sagitta, (double)dist[j]);
		}
		qaws_internal_free(pts);
		qaws_internal_free(loc);
		if (st != QAWS_STATUS_OK)
			goto done;

		/* the displacement spread to the inside nodes: Jacobi averaging over the elements' other nodes (the vertices still: it fades inward) */
		for (k = 0; k < sweeps; k++)
		{
			memset(acc, 0, sizeof(double) * 3 * (size_t)nn);
			memset(cnt, 0, sizeof(unsigned int) * (size_t)nn);
			for (c = 0; c < mesh->nc; c++)
			{
				unsigned int n0 = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? 4u : 8u, n1 = eo[c + 1] - eo[c], a, b;
				for (a = n0; a < n1; a++)
				{
					unsigned int s = en[eo[c] + a];
					if (isb[s])
						continue;
					for (b = 0; b < n1; b++)
						if (b != a)
						{
							unsigned int t = en[eo[c] + b];
							for (d = 0; d < 3; d++) acc[3 * (size_t)s + d] += D[3 * (size_t)t + d];
							cnt[s]++;
						}
				}
			}
			for (i = mesh->nv; i < nn; i++)
				if (!isb[i] && cnt[i])
					for (d = 0; d < 3; d++) D[3 * (size_t)i + d] = acc[3 * (size_t)i + d] / cnt[i];
		}

		/* the faces inside the mesh from their edges (Coons), then a hexahedron's own inside nodes from its boundary
		   (Gordon and Hall): transfinite, so they follow the border's curving rather than lag behind it */
		for (c = 0; c < mesh->nc; c++)
		{
			vcu_el const* e = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? &lt : &lh;
			unsigned int s;
			for (s = 0; s < e->nn; s++)
				if (e->kind[s] == 2 && !isb[en[eo[c] + s]])
				{
					double dc[3];
					vcu_coons(dc, e, en + eo[c], s, D);
					for (d = 0; d < 3; d++) D[3 * (size_t)en[eo[c] + s] + d] = dc[d];
				}
		}
		for (c = 0; c < mesh->nc; c++)
			if (mesh->type[c] == QAWS_CELL_HEXAHEDRON)
				vcu_gordon_hall(&lh, en + eo[c], D);

		/* the Jacobians: failing elements' displacements halved, then straight */
		for (k = 0; k <= VCU_ROUNDS; k++)
		{
			unsigned int bad = 0;
			for (c = 0; c < mesh->nc; c++)
			{
				unsigned int a;
				vcu_el const* e = mesh->type[c] == QAWS_CELL_TETRAHEDRON ? &lt : &lh;
				if (vcu_ratio(&B, e, en + eo[c], X, D, work, minj) >= minj)
					continue;
				bad++;
				flag[c] = 1;
				for (a = e->hex ? 8u : 4u; a < e->nn; a++)
					for (d = 0; d < 3; d++) D[3 * (size_t)en[eo[c] + a] + d] *= k < VCU_ROUNDS ? 0.5 : 0;
			}
			if (!bad)
				break;
		}
		for (c = 0; c < mesh->nc; c++) rep.relaxed += flag[c];
		{
			unsigned int j;
			qaws_scalar* pts2 = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * 4 * (size_t)nb);
			qaws_mesh_location* loc2 = (qaws_mesh_location*)qaws_internal_scratch_malloc(sizeof(qaws_mesh_location) * (size_t)nb);
			if (!pts2 || !loc2)
			{
				qaws_internal_free(pts2);
				qaws_internal_free(loc2);
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			for (j = 0; j < nb; j++)
				for (d = 0; d < 3; d++) pts2[3 * (size_t)j + d] = (qaws_scalar)(X[3 * (size_t)bord[j] + d] + D[3 * (size_t)bord[j] + d]);
			st = qaws_mesh_reference_project(loc2, pts2 + 3 * (size_t)nb, reference, pts2, NULL, nb);
			for (j = 0; st == QAWS_STATUS_OK && j < nb; j++)
				rep.max_distance = (qaws_scalar)fmax((double)rep.max_distance, (double)pts2[3 * (size_t)nb + j]);
			qaws_internal_free(pts2);
			qaws_internal_free(loc2);
			if (st != QAWS_STATUS_OK)
				goto done;
		}
	}
	rep.min_jacobian = (qaws_scalar)HUGE_VAL;
	for (c = 0; c < mesh->nc; c++)
		rep.min_jacobian = (qaws_scalar)fmin((double)rep.min_jacobian, vcu_ratio(&B, mesh->type[c] == QAWS_CELL_TETRAHEDRON ? &lt : &lh, en + eo[c], X, D, work, minj));
	for (i = 0; i < 3 * nn; i++) out_positions[i] = (qaws_scalar)(X[i] + D[i]);
done:
	if (report)
		*report = rep;
	vcu_bern_free(&B);
	qaws_internal_free(eo);
	qaws_internal_free(isb);
	qaws_internal_free(flag);
	qaws_internal_free(X);
	qaws_internal_free(D);
	qaws_internal_free(acc);
	qaws_internal_free(cnt);
	qaws_internal_free(bord);
	qaws_internal_free(bc);
	qaws_internal_free(bs);
	qaws_internal_free(work);
	qaws_internal_free(arc);
	return st;
}

qaws_status qaws_volume_mesh_quadratic_nodes(qaws_scalar* out_positions, unsigned int* out_element_nodes, qaws_volume_curve_report* report,
	qaws_volume_mesh const* mesh, qaws_mesh_reference const* reference, qaws_volume_curve_desc const* desc)
{
	return vcu_nodes(out_positions, out_element_nodes, report, mesh, reference, desc, 2);
}

qaws_status qaws_volume_mesh_cubic_nodes(qaws_scalar* out_positions, unsigned int* out_element_nodes, qaws_volume_curve_report* report,
	qaws_volume_mesh const* mesh, qaws_mesh_reference const* reference, qaws_volume_curve_desc const* desc)
{
	return vcu_nodes(out_positions, out_element_nodes, report, mesh, reference, desc, 3);
}
