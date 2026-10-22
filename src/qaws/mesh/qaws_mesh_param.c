/*
 * Parameterization (qaws_mesh_param.h): the mesh cut along its seams into
 * wedges (a vertex's corners between seams or the border), LSCM as the
 * conformal energy's sparse system, ARAP's local-global rounds on the
 * wedges' cotangent Laplacian (qaws_mesh_smooth.c's sparse helpers).
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_param.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PM_NONE 0xFFFFFFFFu

static int pm_cmp_u64(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

static int pm_cmp_double(void const* p, void const* q)
{
	double a = *(double const*)p, b = *(double const*)q;
	return a < b ? -1 : a > b;
}

static unsigned int pm_find(unsigned int* par, unsigned int i)
{
	while (par[i] != i)
	{
		par[i] = par[par[i]];
		i = par[i];
	}
	return i;
}

static void pm_join(unsigned int* par, unsigned int a, unsigned int b)
{
	a = pm_find(par, a);
	b = pm_find(par, b);
	if (a != b)
		par[a < b ? b : a] = a < b ? a : b;
}

static unsigned long long pm_key(unsigned int a, unsigned int b)
{
	return a < b ? ((unsigned long long)a << 32 | b) : ((unsigned long long)b << 32 | a);
}

/* a sparse symmetric matrix from triplets (row, column, value; merged), rows sorted */
typedef struct pm_trip
{
	unsigned int r, c;
	double v;
} pm_trip;

static int pm_cmp_trip(void const* p, void const* q)
{
	pm_trip const *a = (pm_trip const*)p, *b = (pm_trip const*)q;
	return a->r != b->r ? (a->r < b->r ? -1 : 1) : a->c != b->c ? (a->c < b->c ? -1 : 1) : 0;
}

static int pm_csr(sm_csr* A, pm_trip* T, unsigned int nt, unsigned int n)
{
	unsigned int i, m = 0;
	memset(A, 0, sizeof *A);
	qsort(T, nt, sizeof *T, pm_cmp_trip);
	A->n = n;
	A->rp = (unsigned int*)qaws_internal_zalloc((size_t)n + 1, sizeof(unsigned int));
	A->ci = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	A->v = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nt + 1));
	if (!A->rp || !A->ci || !A->v)
	{
		sm_csr_free(A);
		return 0;
	}
	for (i = 0; i < nt; i++)
	{
		if (m && T[i].r == T[i - 1].r && T[i].c == T[i - 1].c)
		{
			A->v[m - 1] += T[i].v;
			continue;
		}
		A->ci[m] = T[i].c;
		A->v[m++] = T[i].v;
		A->rp[T[i].r + 1]++;
	}
	for (i = 0; i < n; i++)
		A->rp[i + 1] += A->rp[i];
	return 1;
}

/* triangle f laid out isometrically in the plane: X (3 corners x 2) */
static void pm_layout(double const* x, unsigned int const* t, unsigned int f, double (*X)[2])
{
	double const *a = x + 3 * (size_t)t[3 * (size_t)f], *b = x + 3 * (size_t)t[3 * (size_t)f + 1], *c = x + 3 * (size_t)t[3 * (size_t)f + 2];
	double e[3], g[3], l, px, n2;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		e[k] = b[k] - a[k];
		g[k] = c[k] - a[k];
	}
	l = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
	px = l > 0 ? (e[0] * g[0] + e[1] * g[1] + e[2] * g[2]) / l : 0;
	n2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2] - px * px;
	X[0][0] = X[0][1] = 0;
	X[1][0] = l;
	X[1][1] = 0;
	X[2][0] = px;
	X[2][1] = n2 > 0 ? sqrt(n2) : 0;
}

/* the cotangent at corner e of a triangle laid out in the plane (P: 3 corners x 2) */
static double pm_cot2(double const* P, unsigned int e)
{
	double u[2] = { P[2 * ((e + 1) % 3)] - P[2 * e], P[2 * ((e + 1) % 3) + 1] - P[2 * e + 1] };
	double w[2] = { P[2 * ((e + 2) % 3)] - P[2 * e], P[2 * ((e + 2) % 3) + 1] - P[2 * e + 1] };
	double cr = fabs(u[0] * w[1] - u[1] * w[0]);
	return cr > 0 ? (u[0] * w[0] + u[1] * w[1]) / cr : 0;
}

static double pm_angle(double const* o, double const* p, double const* q, unsigned int dim)
{
	double u[3] = { 0, 0, 0 }, w[3] = { 0, 0, 0 }, c[3];
	unsigned int k;
	for (k = 0; k < dim; k++)
	{
		u[k] = p[k] - o[k];
		w[k] = q[k] - o[k];
	}
	c[0] = u[1] * w[2] - u[2] * w[1];
	c[1] = u[2] * w[0] - u[0] * w[2];
	c[2] = u[0] * w[1] - u[1] * w[0];
	return atan2(sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]), u[0] * w[0] + u[1] * w[1] + u[2] * w[2]);
}

/* BFF (Sawhney and Crane 2017) on the wedges (W, positions xw, K their cotangent Laplacian, corner i's wedge wedge[i], bord[i]:
   the half-edge from corner i on the border; a pin per chart): the scale factors u from the curvatures (Yamabe: K u = the
   angle defects inside; disk: the border's turning set to 2 pi times its share of the border's length, else u = 0 there,
   the lengths kept), the border's new turning k - (K u)_B and lengths e^((u_i + u_j) / 2) l, closed by the least change of
   the lengths (weighted by them), its curve integrated; the inside its harmonic extension. Into uv (2 per wedge).
   INVALID_TOPOLOGY for a chart without exactly one border loop. */
static qaws_status pm_bff(double* uv, double const* xw, unsigned int W, unsigned int const* wedge, unsigned int nt,
	unsigned char const* bord, unsigned int const* chart, unsigned int C, sm_csr const* K, unsigned int const* pin, int disk,
	unsigned int const* cv, qaws_scalar const* turn)
{
	unsigned int n = 3 * nt, i, j, k, *nx = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)W + 1));
	unsigned int *wch = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)W + 1)), *loops = (unsigned int*)qaws_internal_zalloc((size_t)C + 1, sizeof(unsigned int));
	double *ang = (double*)qaws_internal_zalloc((size_t)W + 1, sizeof(double)), *u = (double*)qaws_internal_zalloc((size_t)W + 1, sizeof(double));
	double *rhs = (double*)qaws_internal_zalloc(2 * (size_t)W + 1, sizeof(double)), *ku = (double*)qaws_internal_zalloc((size_t)W + 1, sizeof(double));
	double *blen = (double*)qaws_internal_zalloc((size_t)W + 1, sizeof(double)), *ltot = (double*)qaws_internal_zalloc((size_t)C + 1, sizeof(double));
	unsigned char *isb = (unsigned char*)qaws_internal_zalloc((size_t)W + 1, 1), *move = (unsigned char*)qaws_internal_malloc((size_t)W + 1);
	unsigned char* seen = (unsigned char*)qaws_internal_zalloc((size_t)W + 1, 1);
	double *bt = (double*)qaws_internal_zalloc((size_t)W + 1, sizeof(double)), *bsum = (double*)qaws_internal_zalloc((size_t)C + 1, sizeof(double));
	sm_system S;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&S, 0, sizeof S);
	if (!nx || !wch || !loops || !ang || !u || !rhs || !ku || !blen || !ltot || !isb || !move || !seen || !bt || !bsum)
		goto done;
	disk = disk || turn;
	for (i = 0; i < W; i++) nx[i] = PM_NONE;
	/* the corners' angles summed per wedge; the border's loops (next wedge along each, interior on its left) */
	for (i = 0; i < n; i++)
	{
		unsigned int a = 3 * (i / 3) + (i % 3 + 1) % 3, b = 3 * (i / 3) + (i % 3 + 2) % 3;
		ang[wedge[i]] += pm_angle(xw + 3 * (size_t)wedge[i], xw + 3 * (size_t)wedge[a], xw + 3 * (size_t)wedge[b], 3);
		wch[wedge[i]] = chart[i / 3];
		if (bord[i])
		{
			unsigned int w0 = wedge[i], w1 = wedge[a];
			double l = sqrt((xw[3 * (size_t)w1] - xw[3 * (size_t)w0]) * (xw[3 * (size_t)w1] - xw[3 * (size_t)w0]) + (xw[3 * (size_t)w1 + 1] - xw[3 * (size_t)w0 + 1]) * (xw[3 * (size_t)w1 + 1] - xw[3 * (size_t)w0 + 1]) +
				(xw[3 * (size_t)w1 + 2] - xw[3 * (size_t)w0 + 2]) * (xw[3 * (size_t)w1 + 2] - xw[3 * (size_t)w0 + 2]));
			if (nx[w0] != PM_NONE)
			{
				st = QAWS_STATUS_INVALID_TOPOLOGY;   /* (two border edges leave a wedge: not a disk) */
				goto done;
			}
			nx[w0] = w1;
			isb[w0] = isb[w1] = 1;
			if (turn)
			{
				bt[w0] = (double)turn[cv[i]];   /* (the prescribed turning, once per border wedge) */
				bsum[chart[i / 3]] += bt[w0];
			}
			blen[w0] += 0.5 * l;
			blen[w1] += 0.5 * l;
			ltot[chart[i / 3]] += l;
		}
	}
	for (i = 0; i < W; i++)
		if (isb[i] && !seen[i])
		{
			unsigned int w = i, guard = 0;
			if (nx[i] == PM_NONE)
			{
				st = QAWS_STATUS_INVALID_TOPOLOGY;
				goto done;
			}
			do
			{
				seen[w] = 1;
				w = nx[w];
			} while (w != PM_NONE && w != i && ++guard <= W);
			if (w != i || ++loops[wch[i]] > 1)
			{
				st = QAWS_STATUS_INVALID_TOPOLOGY;   /* (a chart with holes, or an open loop) */
				goto done;
			}
		}
	/* the scale factors: K u = the angle defects inside; on the border u = 0 (lengths) or K u = k - the target turning (disk:
	   the prescribed turning, what it leaves of 2 pi spread by the border's length) */
	for (i = 0; i < W; i++)
	{
		move[i] = !isb[i] || disk;
		rhs[i] = isb[i] ? (disk ? (3.14159265358979323846 - ang[i]) - bt[i] - (6.283185307179586 - bsum[wch[i]]) * blen[i] / ltot[wch[i]] : 0) : 6.283185307179586 - ang[i];
	}
	for (i = 0; disk && i < C; i++) move[pin[i]] = 0;
	st = sm_system_factor(&S, K, xw, move);
	/* (a disc: the turning is a share of the border's length after scaling, which u sets: rounds until it settles) */
	for (k = 0; st == QAWS_STATUS_OK && k < (disk ? 80u : 1u); k++)
	{
		st = sm_system_solve(&S, u, rhs, 1);
		if (st != QAWS_STATUS_OK || !disk)
			break;
		memset(blen, 0, sizeof(double) * W);
		memset(ltot, 0, sizeof(double) * C);
		for (i = 0; i < W; i++)
			if (isb[i])
			{
				unsigned int v = nx[i];
				double l = exp(0.5 * (u[i] + u[v])) * sqrt((xw[3 * (size_t)v] - xw[3 * (size_t)i]) * (xw[3 * (size_t)v] - xw[3 * (size_t)i]) +
					(xw[3 * (size_t)v + 1] - xw[3 * (size_t)i + 1]) * (xw[3 * (size_t)v + 1] - xw[3 * (size_t)i + 1]) + (xw[3 * (size_t)v + 2] - xw[3 * (size_t)i + 2]) * (xw[3 * (size_t)v + 2] - xw[3 * (size_t)i + 2]));
				blen[i] += 0.5 * l;
				blen[v] += 0.5 * l;
				ltot[wch[i]] += l;
			}
		for (i = 0; i < W; i++)
			if (isb[i])
				rhs[i] = 0.5 * rhs[i] + 0.5 * ((3.14159265358979323846 - ang[i]) - bt[i] - (6.283185307179586 - bsum[wch[i]]) * blen[i] / ltot[wch[i]]);
		for (i = 0; i < C; i++) u[pin[i]] = 0;
	}
	sm_system_free(&S);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < W; i++)
		for (j = K->rp[i]; j < K->rp[i + 1]; j++) ku[i] += K->v[j] * u[K->ci[j]];
	/* each chart's border: turning k - (K u), lengths scaled, closed, integrated */
	memset(seen, 0, W);
	for (i = 0; i < W; i++)
		if (isb[i] && !seen[i])
		{
			double phi = 0, M[4] = { 0, 0, 0, 0 }, r[2] = { 0, 0 }, det, lam[2], px = 0, py = 0;
			unsigned int w = i;
			/* (pass 1: the closure's 2 x 2 system; pass 2: the curve) */
			for (k = 0; k < 2; k++)
			{
				w = i;
				phi = 0;
				px = py = 0;
				do
				{
					unsigned int v = nx[w];
					double l = sqrt((xw[3 * (size_t)v] - xw[3 * (size_t)w]) * (xw[3 * (size_t)v] - xw[3 * (size_t)w]) + (xw[3 * (size_t)v + 1] - xw[3 * (size_t)w + 1]) * (xw[3 * (size_t)v + 1] - xw[3 * (size_t)w + 1]) +
						(xw[3 * (size_t)v + 2] - xw[3 * (size_t)w + 2]) * (xw[3 * (size_t)v + 2] - xw[3 * (size_t)w + 2]));
					double lt = exp(0.5 * (u[w] + u[v])) * l, c, s;
					if (w != i)
						phi += (3.14159265358979323846 - ang[w]) - ku[w];
					c = cos(phi);
					s = sin(phi);
					if (k == 0)
					{
						M[0] += l * c * c;
						M[1] += l * c * s;
						M[3] += l * s * s;
						r[0] += lt * c;
						r[1] += lt * s;
					}
					else
					{
						double ls = lt - l * (c * lam[0] + s * lam[1]);
						uv[2 * (size_t)w] = px;
						uv[2 * (size_t)w + 1] = py;
						seen[w] = 1;
						px += ls * c;
						py += ls * s;
					}
					w = v;
				} while (w != i);
				if (k == 0)
				{
					M[2] = M[1];
					det = M[0] * M[3] - M[1] * M[2];
					lam[0] = det != 0 ? (M[3] * r[0] - M[1] * r[1]) / det : 0;
					lam[1] = det != 0 ? (M[0] * r[1] - M[2] * r[0]) / det : 0;
				}
			}
		}
	/* the inside: the curve's harmonic extension */
	for (i = 0; i < W; i++) move[i] = !isb[i];
	memset(rhs, 0, sizeof(double) * 2 * W);
	st = sm_system_factor(&S, K, xw, move);
	if (st == QAWS_STATUS_OK)
		st = sm_system_solve(&S, uv, rhs, 2);
done:
	sm_system_free(&S);
	qaws_internal_free(bt); qaws_internal_free(bsum);
	qaws_internal_free(nx); qaws_internal_free(wch); qaws_internal_free(loops); qaws_internal_free(ang); qaws_internal_free(u); qaws_internal_free(rhs);
	qaws_internal_free(ku); qaws_internal_free(blen); qaws_internal_free(ltot); qaws_internal_free(isb); qaws_internal_free(move); qaws_internal_free(seen);
	return st;
}

/* The charts packed by their shapes (as xatlas does): each rasterized at the pixel size of a 96-pixel square of their whole
   area over 0.7, grown by the gap and a pixel, placed largest first at its lowest (then leftmost) free spot in a strip,
   turned a quarter when that sits lower; the strip's width the best of a few around the square. uv turned in place for
   the charts turned; xy: each chart's box corner (its rotated box's) in uv units; *out_side the square's side. 0 on
   allocation failure */
static int pm_pack_bitmap(double* uv, unsigned int const* wedge, unsigned int const* chart, unsigned int n, unsigned int W, unsigned int C,
	double gap, double* xy, double* out_side)
{
	unsigned int c, i, k, r, *wch = NULL, *ord = NULL, *pix = NULL, *p0 = NULL, np = 0, capp = 0, gw, best_gw = 0, trial;
	double *bb = NULL, area = 0, ps, best_side = HUGE_VAL;
	unsigned char *occ = NULL, *rot = NULL;
	int ok = 0;
	wch = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)W + 1));
	ord = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)C + 1));
	bb = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 5, (size_t)C + 1));   /* min x, min y, max x, max y, area */
	p0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)C + 2));   /* per chart and turn its pixels' start */
	rot = (unsigned char*)qaws_internal_zalloc((size_t)C + 1, 1);
	if (!wch || !ord || !bb || !p0 || !rot)
		goto done;
	for (i = 0; i < W; i++)
		wch[i] = PM_NONE;
	for (i = 0; i < n; i++)
		wch[wedge[i]] = chart[i / 3];
	for (c = 0; c < C; c++)
	{
		bb[5 * (size_t)c] = bb[5 * (size_t)c + 1] = HUGE_VAL;
		bb[5 * (size_t)c + 2] = bb[5 * (size_t)c + 3] = -HUGE_VAL;
		bb[5 * (size_t)c + 4] = 0;
		ord[c] = c;
	}
	for (i = 0; i < W; i++)
		if (wch[i] != PM_NONE)
		{
			double* b = bb + 5 * (size_t)wch[i];
			b[0] = fmin(b[0], uv[2 * (size_t)i]);
			b[1] = fmin(b[1], uv[2 * (size_t)i + 1]);
			b[2] = fmax(b[2], uv[2 * (size_t)i]);
			b[3] = fmax(b[3], uv[2 * (size_t)i + 1]);
		}
	for (i = 0; i < n; i += 3)
	{
		double const *a = uv + 2 * (size_t)wedge[i], *b = uv + 2 * (size_t)wedge[i + 1], *d = uv + 2 * (size_t)wedge[i + 2];
		double s = 0.5 * fabs((b[0] - a[0]) * (d[1] - a[1]) - (d[0] - a[0]) * (b[1] - a[1]));
		bb[5 * (size_t)chart[i / 3] + 4] += s;
		area += s;
	}
	if (!(area > 0))
		goto done;
	ps = sqrt(area / 0.7) / 192;
	/* largest first */
	for (i = 1; i < C; i++)
		for (k = i; k > 0 && bb[5 * (size_t)ord[k] + 4] > bb[5 * (size_t)ord[k - 1] + 4]; k--)
		{
			unsigned int s = ord[k];
			ord[k] = ord[k - 1];
			ord[k - 1] = s;
		}
	/* each chart's pixels, unturned and turned a quarter (x, y packed in one int each pair), the gap and a pixel grown */
	{
		int g = (int)ceil(0.5 * gap / ps);
		for (c = 0; c < C; c++)
			for (r = 0; r < 2; r++)
			{
				double const* b = bb + 5 * (size_t)c;
				int bw = (int)ceil((r ? b[3] - b[1] : b[2] - b[0]) / ps) + 1 + 2 * g, bh = (int)ceil((r ? b[2] - b[0] : b[3] - b[1]) / ps) + 1 + 2 * g, px, py;
				unsigned char* m = (unsigned char*)qaws_internal_zalloc((size_t)bw * (size_t)bh + 1, 1);
				if (!m)
					goto done;
				p0[2 * (size_t)c + r] = np;
				for (i = 0; i < n; i += 3)
				{
					double P[3][2];
					int x0, x1, y0, y1;
					if (chart[i / 3] != c)
						continue;
					for (k = 0; k < 3; k++)
					{
						double u = uv[2 * (size_t)wedge[i + k]] - b[0], v = uv[2 * (size_t)wedge[i + k] + 1] - b[1];
						P[k][0] = (r ? (b[3] - b[1]) - v : u) / ps + g;
						P[k][1] = (r ? u : v) / ps + g;
					}
					x0 = (int)floor(fmin(P[0][0], fmin(P[1][0], P[2][0])));
					x1 = (int)ceil(fmax(P[0][0], fmax(P[1][0], P[2][0])));
					y0 = (int)floor(fmin(P[0][1], fmin(P[1][1], P[2][1])));
					y1 = (int)ceil(fmax(P[0][1], fmax(P[1][1], P[2][1])));
					for (py = y0; py <= y1; py++)
						for (px = x0; px <= x1; px++)
						{
							/* the pixel touched: its centre within the triangle grown by a pixel's half diagonal */
							double cx = px + 0.5, cy = py + 0.5, ar = (P[1][0] - P[0][0]) * (P[2][1] - P[0][1]) - (P[2][0] - P[0][0]) * (P[1][1] - P[0][1]), in = 1;
							for (k = 0; k < 3 && in; k++)
							{
								double const *e0 = P[k], *e1 = P[(k + 1) % 3];
								double ex = e1[0] - e0[0], ey = e1[1] - e0[1], l = sqrt(ex * ex + ey * ey);
								double sd = ((cx - e0[0]) * ey - (cy - e0[1]) * ex) * (ar > 0 ? -1 : 1);
								in = l == 0 || sd / l > -0.71;
							}
							if (in && px >= 0 && py >= 0 && px < bw && py < bh)
								m[(size_t)py * bw + px] = 1;
						}
				}
				/* grown by g pixels (the gap), then the set pixels listed */
				{
					unsigned char* m2 = (unsigned char*)qaws_internal_zalloc((size_t)bw * (size_t)bh + 1, 1);
					int dx, dy;
					if (!m2)
					{
						qaws_internal_free(m);
						goto done;
					}
					for (py = 0; py < bh; py++)
						for (px = 0; px < bw; px++)
							if (m[(size_t)py * bw + px])
								for (dy = -g; dy <= g; dy++)
									for (dx = -g; dx <= g; dx++)
										if (px + dx >= 0 && py + dy >= 0 && px + dx < bw && py + dy < bh)
											m2[(size_t)(py + dy) * bw + px + dx] = 1;
					for (py = 0; py < bh; py++)
						for (px = 0; px < bw; px++)
							if (m2[(size_t)py * bw + px])
							{
								if (np + 1 > capp)
								{
									unsigned int nc2 = qaws_internal_grow(capp, (size_t)np + 1, 4096);
									unsigned int* g2 = nc2 ? (unsigned int*)qaws_internal_realloc(pix, qaws_internal_mul2(sizeof(unsigned int), nc2)) : NULL;
									if (!g2)
									{
										qaws_internal_free(m);
										qaws_internal_free(m2);
										goto done;
									}
									pix = g2;
									capp = nc2;
								}
								pix[np++] = ((unsigned int)py << 16) | (unsigned int)px;
							}
					qaws_internal_free(m2);
				}
				qaws_internal_free(m);
			}
		p0[2 * (size_t)C] = np;
	}
	/* the strip's width: 0.8 .. 1.4 of the square, the lowest square kept; then laid out at it */
	for (trial = 0; trial < 8; trial++)
	{
		unsigned int H = 0, Wu = 0, side, wmin = 0;
		gw = trial < 7 ? (unsigned int)(192 * (0.8 + 0.1 * trial)) : best_gw;
		for (c = 0; c < C; c++)
			for (r = 0; r < 2; r++)
			{
				unsigned int e = p0[2 * (size_t)c + r + 1], j, mx = 0;
				for (j = p0[2 * (size_t)c + r]; j < e; j++)
					mx = mx > (pix[j] & 0xFFFF) ? mx : (pix[j] & 0xFFFF);
				if (!r)
					wmin = wmin > mx + 1 ? wmin : mx + 1;
			}
		if (gw < wmin)
			gw = wmin;
		qaws_internal_free(occ);
		occ = (unsigned char*)qaws_internal_zalloc((size_t)gw * 4096 + 1, 1);
		if (!occ)
			goto done;
		for (i = 0; i < C; i++)
		{
			unsigned int bx = 0, byy = 0, br = 0, bt = 0xFFFFFFFFu, j, bpw = 0;
			c = ord[i];
			for (r = 0; r < 2; r++)
			{
				unsigned int s0 = p0[2 * (size_t)c + r], e = p0[2 * (size_t)c + r + 1], pw = 0, ph = 0, x, y, found = 0;
				for (j = s0; j < e; j++)
				{
					pw = pw > (pix[j] & 0xFFFF) + 1 ? pw : (pix[j] & 0xFFFF) + 1;
					ph = ph > (pix[j] >> 16) + 1 ? ph : (pix[j] >> 16) + 1;
				}
				if (pw > gw)
					continue;
				for (y = 0; y + ph < 4096 && !found && y < bt; y++)
					for (x = 0; x + pw <= gw && !found; x++)
					{
						for (j = s0; j < e && !occ[(size_t)(y + (pix[j] >> 16)) * gw + x + (pix[j] & 0xFFFF)]; j++)
							;
						if (j == e)
						{
							found = 1;
							if (y + ph < bt)
							{
								bt = y + ph;
								bx = x;
								byy = y;
								br = r;
								bpw = pw;
							}
						}
					}
			}
			if (bt == 0xFFFFFFFFu)
				goto done;   /* (no room in 4096 rows: give up, the skyline instead) */
			{
				unsigned int s0 = p0[2 * (size_t)c + br], e = p0[2 * (size_t)c + br + 1];
				for (j = s0; j < e; j++)
					occ[(size_t)(byy + (pix[j] >> 16)) * gw + bx + (pix[j] & 0xFFFF)] = 1;
			}
			H = H > bt ? H : bt;
			Wu = Wu > bx + bpw ? Wu : bx + bpw;
			if (trial == 7)
			{
				rot[c] = (unsigned char)br;
				xy[2 * (size_t)c] = bx * ps;
				xy[2 * (size_t)c + 1] = byy * ps;
			}
		}
		side = H > Wu ? H : Wu;
		if (trial < 7 && side * ps < best_side)
		{
			best_side = side * ps;
			best_gw = gw;
		}
		if (trial == 7)
			best_side = side * ps;
	}
	/* the charts turned a quarter: u' = (box height) - (v - min v), v' = u - min u, then their box from 0; the others from 0 */
	for (i = 0; i < W; i++)
		if (wch[i] != PM_NONE)
		{
			double const* b = bb + 5 * (size_t)wch[i];
			double u = uv[2 * (size_t)i] - b[0], v = uv[2 * (size_t)i + 1] - b[1];
			int g = (int)ceil(0.5 * gap / ps);
			uv[2 * (size_t)i] = (rot[wch[i]] ? (b[3] - b[1]) - v : u) + g * ps;
			uv[2 * (size_t)i + 1] = (rot[wch[i]] ? u : v) + g * ps;
		}
	*out_side = best_side;
	ok = 1;
done:
	qaws_internal_free(wch); qaws_internal_free(ord); qaws_internal_free(bb); qaws_internal_free(p0); qaws_internal_free(rot);
	qaws_internal_free(pix); qaws_internal_free(occ);
	return ok;
}

/* Boxes (wh: width and height per box) packed bottom-left on a skyline in a strip W wide, in the order ord: each where its
   top is lowest (then leftmost), the skyline's segments (sk: 2 per segment, its start and height, 2 C + 2 of them)
   raised under it. xy: each box's corner, or NULL; *out_w, *out_h: the extent used */
static void pm_skyline(double const* wh, unsigned int const* ord, unsigned int C, double W, double* sk, double* xy, double* out_w, double* out_h)
{
	unsigned int i, s, nsk = 1;
	double used_w = 0, used_h = 0;
	sk[0] = 0;
	sk[1] = 0;
	for (i = 0; i < C; i++)
	{
		unsigned int c = ord[i], bs = 0, a, b, j, n2;
		double w = wh[2 * (size_t)c], hgt = wh[2 * (size_t)c + 1], best = HUGE_VAL, bx = 0;
		for (s = 0; s < nsk; s++)
		{
			double x = sk[2 * (size_t)s], top = 0;
			if (s > 0 && x + w > W * (1 + 1e-12))
				break;
			for (j = s; j < nsk && sk[2 * (size_t)j] < x + w - 1e-12; j++)
				top = fmax(top, sk[2 * (size_t)j + 1]);
			if (top < best - 1e-12)
			{
				best = top;
				bx = x;
				bs = s;
			}
		}
		if (xy)
		{
			xy[2 * (size_t)c] = bx;
			xy[2 * (size_t)c + 1] = best;
		}
		used_w = fmax(used_w, bx + w);
		used_h = fmax(used_h, best + hgt);
		/* the skyline: the segments before bx, the box's top over [bx, bx + w), the rest from bx + w on */
		for (a = bs; a + 1 < nsk && sk[2 * (size_t)(a + 1)] < bx + w - 1e-12; a++)
			;
		{
			double resty = sk[2 * (size_t)a + 1], tmp[8];
			unsigned int k = 0;
			tmp[k++] = bx;
			tmp[k++] = best + hgt;
			b = a + 1;
			if (b >= nsk || sk[2 * (size_t)b] > bx + w + 1e-12)
			{
				tmp[k++] = bx + w;
				tmp[k++] = resty;   /* (the segment the box ends in goes on past it) */
			}
			n2 = bs + k / 2 + (nsk - b);
			memmove(sk + 2 * ((size_t)bs + k / 2), sk + 2 * (size_t)b, sizeof(double) * 2 * (nsk - b));
			memcpy(sk + 2 * (size_t)bs, tmp, sizeof(double) * k);
			nsk = n2;
		}
	}
	*out_w = used_w;
	*out_h = used_h;
}

/* each chart turned to its smallest bounding box (16 angles in a half turn), then wider than tall; wch: W scratch */
static void pm_turn_charts(double* uv, unsigned int* wch, double* bb, unsigned int const* wedge, unsigned int const* chart, unsigned int n,
	unsigned int W, unsigned int C)
{
	unsigned int i, k;
	for (i = 0; i < W; i++)
		wch[i] = PM_NONE;
	for (i = 0; i < n; i++)
		wch[wedge[i]] = chart[i / 3];
	for (i = 0; i < 64 * C; i++)
		bb[i] = i % 4 < 2 ? HUGE_VAL : -HUGE_VAL;
	for (i = 0; i < W; i++)
		for (k = 0; wch[i] != PM_NONE && k < 16; k++)
		{
			double a = k * 3.14159265358979323846 / 16, c = cos(a), s = sin(a), *b = bb + 64 * (size_t)wch[i] + 4 * k;
			double X = c * uv[2 * (size_t)i] - s * uv[2 * (size_t)i + 1], Y = s * uv[2 * (size_t)i] + c * uv[2 * (size_t)i + 1];
			b[0] = fmin(b[0], X);
			b[1] = fmin(b[1], Y);
			b[2] = fmax(b[2], X);
			b[3] = fmax(b[3], Y);
		}
	for (i = 0; i < C; i++)
	{
		unsigned int best = 0;
		double* b = bb + 64 * (size_t)i;
		for (k = 1; k < 16; k++)
			if ((b[4 * k + 2] - b[4 * k]) * (b[4 * k + 3] - b[4 * k + 1]) < (b[4 * best + 2] - b[4 * best]) * (b[4 * best + 3] - b[4 * best + 1]))
				best = k;
		b[0] = best * 3.14159265358979323846 / 16 + (b[4 * best + 3] - b[4 * best + 1] > b[4 * best + 2] - b[4 * best] ? 1.5707963267948966 : 0);
	}
	for (i = 0; i < W; i++)
		if (wch[i] != PM_NONE)
		{
			double a = bb[64 * (size_t)wch[i]], c = cos(a), s = sin(a), X = uv[2 * (size_t)i], Y = uv[2 * (size_t)i + 1];
			uv[2 * (size_t)i] = c * X - s * Y;
			uv[2 * (size_t)i + 1] = s * X + c * Y;
		}
}

/* the parameterization with the given seams (and the axis groups' with auto_seams); tchart: NULL, or each triangle's chart */
static qaws_status pm_parameterize(qaws_scalar* out_uv, qaws_param_report* report, qaws_mesh_3d const* mesh, unsigned int const* seams,
	unsigned int ns, qaws_param_desc const* desc, unsigned int* tchart)
{
	qaws_param_method method = desc ? desc->method : QAWS_PARAM_LSCM;
	unsigned int iters = desc && desc->iterations ? desc->iterations : 10, nv, nt, n, i, j, k, W = 0, C = 0, nb = 0, ntr = 0;
	double *x = NULL, *xw = NULL, *uv = NULL, *mass = NULL, *rhs = NULL;
	unsigned int *t = NULL, *wedge = NULL, *par = NULL, *tpar = NULL, *tw = NULL, *chart = NULL, *pin_a = NULL, *pin_b = NULL;
	unsigned long long *he = NULL, *sk = NULL;
	unsigned char *bord = NULL, *move = NULL, *cls = NULL;
	qaws_scalar* pos = NULL;
	pm_trip* T = NULL;
	sm_csr K, Q;
	sm_system S;
	qaws_param_report local;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(&K, 0, sizeof K);
	memset(&Q, 0, sizeof Q);
	memset(&S, 0, sizeof S);
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out_uv || !mesh || (ns && !seams) || (unsigned int)method > (unsigned int)QAWS_PARAM_BFF)
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	n = 3 * nt;
	pos = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)nv + 1));
	x = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nv + 1));
	t = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	he = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul2(2 * sizeof(unsigned long long), (size_t)n + 1));
	sk = (unsigned long long*)qaws_internal_malloc(sizeof(unsigned long long) * ((size_t)ns + 1));
	par = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	tpar = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	bord = (unsigned char*)qaws_internal_zalloc((size_t)n + 1, 1);
	wedge = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
	if (!pos || !x || !t || !he || !sk || !par || !tpar || !bord || !wedge)
		goto done;
	if ((st = qaws_mesh_3d_get_positions(pos, mesh)) != QAWS_STATUS_OK || (st = qaws_mesh_3d_get_triangles(t, mesh)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 * nv; i++)
		x[i] = (double)pos[i];
	/* half-edges by their edge; the seams' keys */
	for (i = 0; i < n; i++)
	{
		he[2 * (size_t)i] = pm_key(t[i], t[3 * (i / 3) + (i % 3 + 1) % 3]);
		he[2 * (size_t)i + 1] = i;
		par[i] = i;
	}
	for (i = 0; i < nt; i++)
		tpar[i] = i;
	qsort(he, n, 2 * sizeof(unsigned long long), pm_cmp_u64);
	for (i = 0; i < ns; i++)
	{
		unsigned long long key = pm_key(seams[2 * (size_t)i], seams[2 * (size_t)i + 1]);
		unsigned int lo = 0, hi = n;
		while (lo < hi)
		{
			unsigned int mid = (lo + hi) / 2;
			if (he[2 * (size_t)mid] < key)
				lo = mid + 1;
			else
				hi = mid;
		}
		if (seams[2 * (size_t)i] >= nv || seams[2 * (size_t)i + 1] >= nv || lo == n || he[2 * (size_t)lo] != key)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		sk[i] = key;
	}
	qsort(sk, ns, sizeof *sk, pm_cmp_u64);
	if (desc && desc->auto_seams)
	{
		/* each triangle's group: the axis direction nearest its normal (6), an edge between two groups a seam */
		cls = (unsigned char*)qaws_internal_malloc((size_t)nt + 1);
		if (!cls)
			goto done;
		for (i = 0; i < nt; i++)
		{
			double const *a = x + 3 * (size_t)t[3 * i], *b = x + 3 * (size_t)t[3 * i + 1], *c = x + 3 * (size_t)t[3 * i + 2];
			double u[3], w[3], nr[3], best = -1;
			for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; w[k] = c[k] - a[k]; }
			nr[0] = u[1] * w[2] - u[2] * w[1];
			nr[1] = u[2] * w[0] - u[0] * w[2];
			nr[2] = u[0] * w[1] - u[1] * w[0];
			cls[i] = 0;
			for (k = 0; k < 3; k++)
				if (fabs(nr[k]) > best)
				{
					best = fabs(nr[k]);
					cls[i] = (unsigned char)(2 * k + (nr[k] < 0));
				}
		}
		/* the groups by a graph cut (sm_label_cut, as the polycube labels): a triangle's cost for a group its area times 1 - its
		   normal's component along the group's axis (three times its area where that is under 0.3), a shared edge's its
		   length times a quarter of the mean edge between two groups (twice between the two sides of one axis) */
		{
			unsigned int* tn = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
			double* cost = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, (size_t)nt + 1));
			double* ew = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)n + 1));
			double V[36], le = 0;
			int okc = tn && cost && ew;
			for (i = 0; okc && i < n; i++)
				tn[i] = PM_NONE;
			for (i = 0; okc && i + 1 < n; i++)
				if (he[2 * (size_t)i] == he[2 * (size_t)i + 2] && (i + 2 >= n || he[2 * (size_t)i + 4] != he[2 * (size_t)i]) &&
					(i == 0 || he[2 * (size_t)i - 2] != he[2 * (size_t)i]))
				{
					unsigned int c0 = (unsigned int)he[2 * (size_t)i + 1], c1 = (unsigned int)he[2 * (size_t)i + 3];
					tn[c0] = c1 / 3;
					tn[c1] = c0 / 3;
				}
			for (i = 0; okc && i < n; i++)
			{
				double const *a = x + 3 * (size_t)t[i], *b = x + 3 * (size_t)t[3 * (i / 3) + (i % 3 + 1) % 3];
				ew[i] = sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2]));
				le += ew[i];
			}
			le = n ? le / n : 1;
			for (i = 0; okc && i < n; i++)
				ew[i] *= 0.25 * le;
			for (i = 0; okc && i < nt; i++)
			{
				double const *a = x + 3 * (size_t)t[3 * i], *b = x + 3 * (size_t)t[3 * i + 1], *c = x + 3 * (size_t)t[3 * i + 2];
				double u[3], w[3], nr[3], l;
				unsigned int g;
				for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; w[k] = c[k] - a[k]; }
				nr[0] = u[1] * w[2] - u[2] * w[1];
				nr[1] = u[2] * w[0] - u[0] * w[2];
				nr[2] = u[0] * w[1] - u[1] * w[0];
				l = sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
				for (g = 0; g < 6; g++)
				{
					double s = l > 0 ? (g & 1 ? -nr[g / 2] : nr[g / 2]) / l : 0;
					cost[6 * (size_t)i + g] = 0.5 * l * (s > 0.3 ? 1 - s : 3);
				}
			}
			for (i = 0; i < 36; i++)
				V[i] = i / 6 == i % 6 ? 0 : i / 12 == (i % 6) / 2 ? 2 : 1;
			okc = okc && sm_label_cut(cls, cost, nt, 6, tn, ew, V, 4);
			qaws_internal_free(tn);
			qaws_internal_free(cost);
			qaws_internal_free(ew);
			if (!okc)
				goto done;
		}
	}
	/* edges of exactly two triangles turned opposite ways, not seams, join the corners at their ends and the triangles; the
	   others are the cut mesh's border */
	for (i = 0; i < n; i = j)
	{
		for (j = i + 1; j < n && he[2 * (size_t)j] == he[2 * (size_t)i]; j++)
			;
		{
			unsigned long long key = he[2 * (size_t)i];
			unsigned int lo = 0, hi = ns;
			int seam;
			while (lo < hi)
			{
				unsigned int mid = (lo + hi) / 2;
				if (sk[mid] < key)
					lo = mid + 1;
				else
					hi = mid;
			}
			seam = lo < ns && sk[lo] == key;
			if (j == i + 2 && cls && cls[(unsigned int)he[2 * (size_t)i + 1] / 3] != cls[(unsigned int)he[2 * (size_t)i + 3] / 3])
				seam = 1;
			if (j == i + 2 && !seam)
			{
				unsigned int c0 = (unsigned int)he[2 * (size_t)i + 1], c1 = (unsigned int)he[2 * (size_t)i + 3];
				unsigned int f = c0 / 3, g = c1 / 3, a0 = c0, b0 = 3 * f + (c0 % 3 + 1) % 3, a1 = c1, b1 = 3 * g + (c1 % 3 + 1) % 3;
				if (t[a0] == t[b1] && t[b0] == t[a1])
				{
					pm_join(par, a0, b1);
					pm_join(par, b0, a1);
					pm_join(tpar, f, g);
					continue;
				}
			}
			for (k = i; k < j; k++)
				bord[(unsigned int)he[2 * (size_t)k + 1]] = 1;
		}
	}
	qaws_internal_free(he);
	he = NULL;
	/* wedges (the corners' classes) and charts (the triangles') */
	for (i = 0; i < n; i++)
		wedge[i] = pm_find(par, i) == i ? W++ : PM_NONE;
	for (i = 0; i < n; i++)
		wedge[i] = wedge[pm_find(par, i)];
	chart = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	tw = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	xw = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)W + 1));
	uv = (double*)qaws_internal_zalloc(qaws_internal_mul2(2, (size_t)W + 1), sizeof(double));
	if (!chart || !tw || !xw || !uv)
		goto done;
	for (i = 0; i < nt; i++)
		tw[i] = pm_find(tpar, i) == i ? C++ : PM_NONE;
	for (i = 0; i < nt; i++)
		chart[i] = tw[pm_find(tpar, i)];
	for (i = 0; i < n; i++)
		memcpy(xw + 3 * (size_t)wedge[i], x + 3 * (size_t)t[i], 3 * sizeof(double));
	report->charts = C;
	/* pins: per chart a border wedge and the border wedge farthest from it */
	pin_a = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)C + 1));
	pin_b = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)C + 1));
	if (!pin_a || !pin_b)
		goto done;
	for (i = 0; i < C; i++)
		pin_a[i] = pin_b[i] = PM_NONE;
	for (i = 0; i < n; i++)
		if (bord[i] && pin_a[chart[i / 3]] == PM_NONE)
			pin_a[chart[i / 3]] = wedge[i];
	for (i = 0; i < C; i++)
		if (pin_a[i] == PM_NONE)
		{
			st = QAWS_STATUS_INVALID_TOPOLOGY;   /* (a closed chart: nothing to open it) */
			goto done;
		}
	{
		double* far = (double*)qaws_internal_zalloc((size_t)C + 1, sizeof(double));
		if (!far)
			goto done;
		for (i = 0; i < n; i++)
			if (bord[i])
			{
				unsigned int c = chart[i / 3], w = wedge[i];
				double const *a = xw + 3 * (size_t)pin_a[c], *b = xw + 3 * (size_t)w;
				double d = (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
				if (d > far[c])
				{
					far[c] = d;
					pin_b[c] = w;
				}
			}
		for (i = 0; i < C; i++)
		{
			uv[2 * (size_t)pin_a[i]] = uv[2 * (size_t)pin_a[i] + 1] = 0;
			if (pin_b[i] != PM_NONE)
			{
				uv[2 * (size_t)pin_b[i]] = sqrt(far[i]);
				uv[2 * (size_t)pin_b[i] + 1] = 0;
			}
		}
		qaws_internal_free(far);
	}
	/* LSCM: Q = (K on u and on v) - B, 1/2 x^T B x the signed area the border encloses; u and v of a wedge adjacent rows */
	mass = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)W + 1));
	if (!mass)
		goto done;
	{
		unsigned int* twc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)n + 1));
		if (!twc)
			goto done;
		for (i = 0; i < n; i++)
			twc[i] = wedge[i];
		k = (unsigned int)sm_laplacian(&K, mass, xw, W, twc, nt, 0, 0);
		qaws_internal_free(twc);
		if (!k)
			goto done;
	}
	if (method == QAWS_PARAM_BFF)
	{
		st = pm_bff(uv, xw, W, wedge, nt, bord, chart, C, &K, pin_a, desc && desc->bff_disk, t, desc ? desc->bff_turning : NULL);
		if (st != QAWS_STATUS_OK)
			goto done;
		goto laid;
	}
	for (i = 0; i < n; i++)
		nb += bord[i];
	T = (pm_trip*)qaws_internal_malloc(sizeof(pm_trip) * (2 * (size_t)K.rp[W] + 4 * (size_t)nb + 1));
	move = (unsigned char*)qaws_internal_malloc(2 * (size_t)W + 1);
	rhs = (double*)qaws_internal_zalloc(qaws_internal_mul2(2, (size_t)W + 1), sizeof(double));
	if (!T || !move || !rhs)
		goto done;
	for (i = 0; i < W; i++)
	{
		unsigned int p;
		for (p = K.rp[i]; p < K.rp[i + 1]; p++)
		{
			T[ntr].r = 2 * i;
			T[ntr].c = 2 * K.ci[p];
			T[ntr++].v = K.v[p];
			T[ntr].r = 2 * i + 1;
			T[ntr].c = 2 * K.ci[p] + 1;
			T[ntr++].v = K.v[p];
		}
	}
	for (i = 0; i < n; i++)
		if (bord[i])
		{
			unsigned int a = wedge[i], b = wedge[3 * (i / 3) + (i % 3 + 1) % 3];
			pm_trip const e[4] = { { 2 * a, 2 * b + 1, -0.5 }, { 2 * b + 1, 2 * a, -0.5 }, { 2 * b, 2 * a + 1, 0.5 }, { 2 * a + 1, 2 * b, 0.5 } };
			memcpy(T + ntr, e, sizeof e);
			ntr += 4;
		}
	if (!pm_csr(&Q, T, ntr, 2 * W))
		goto done;
	{
		double* pxyz = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, (size_t)W + 1));
		if (!pxyz)
			goto done;
		for (i = 0; i < W; i++)
			for (k = 0; k < 3; k++)
				pxyz[6 * (size_t)i + k] = pxyz[6 * (size_t)i + 3 + k] = xw[3 * (size_t)i + k];   /* (u and v of a wedge stay together) */
		for (i = 0; i < 2 * W; i++)
			move[i] = 1;
		for (i = 0; i < C; i++)
		{
			move[2 * (size_t)pin_a[i]] = move[2 * (size_t)pin_a[i] + 1] = 0;
			if (pin_b[i] != PM_NONE)
				move[2 * (size_t)pin_b[i]] = move[2 * (size_t)pin_b[i] + 1] = 0;
		}
		st = sm_system_factor(&S, &Q, pxyz, move);
		qaws_internal_free(pxyz);
		if (st == QAWS_STATUS_OK)
			st = sm_system_solve(&S, uv, rhs, 1);
		sm_system_free(&S);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	if (method == QAWS_PARAM_ARAP)
	{
		/* each triangle's map nearest a rotation of its own layout; the wedges' Laplacian with one pin per chart, factored once */
		double* X = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 6, (size_t)nt + 1));
		unsigned int it;
		if (!X)
			goto done;
		for (i = 0; i < nt; i++)
			pm_layout(x, t, i, (double(*)[2])(X + 6 * (size_t)i));
		for (i = 0; i < W; i++)
			move[i] = 1;
		for (i = 0; i < C; i++)
			move[pin_a[i]] = 0;
		st = sm_system_factor(&S, &K, xw, move);
		for (it = 0; it < iters && st == QAWS_STATUS_OK; it++)
		{
			memset(rhs, 0, sizeof(double) * 2 * W);
			for (i = 0; i < nt; i++)
			{
				double const* P = X + 6 * (size_t)i;
				double s00 = 0, s01 = 0, s10 = 0, s11 = 0, th, c, s;
				unsigned int e;
				/* the covariance over the edges (opposite corner e: from e + 1 to e + 2), weighted by half the cotangent */
				for (e = 0; e < 3; e++)
				{
					unsigned int p = (e + 1) % 3, r = (e + 2) % 3;
					double const *up = uv + 2 * (size_t)wedge[3 * i + p], *ur = uv + 2 * (size_t)wedge[3 * i + r];
					double w = 0.5 * pm_cot2(P, e), du[2] = { up[0] - ur[0], up[1] - ur[1] }, dx[2] = { P[2 * p] - P[2 * r], P[2 * p + 1] - P[2 * r + 1] };
					s00 += w * du[0] * dx[0];
					s01 += w * du[0] * dx[1];
					s10 += w * du[1] * dx[0];
					s11 += w * du[1] * dx[1];
				}
				th = atan2(s10 - s01, s00 + s11);
				c = cos(th);
				s = sin(th);
				for (e = 0; e < 3; e++)
				{
					unsigned int p = (e + 1) % 3, r = (e + 2) % 3, a = wedge[3 * i + p], b = wedge[3 * i + r];
					double w = 0.5 * pm_cot2(P, e), dx[2] = { P[2 * p] - P[2 * r], P[2 * p + 1] - P[2 * r + 1] };
					double rx = c * dx[0] - s * dx[1], ry = s * dx[0] + c * dx[1];
					rhs[2 * (size_t)a] += w * rx;
					rhs[2 * (size_t)a + 1] += w * ry;
					rhs[2 * (size_t)b] -= w * rx;
					rhs[2 * (size_t)b + 1] -= w * ry;
				}
			}
			st = sm_system_solve(&S, uv, rhs, 2);
		}
		sm_system_free(&S);
		qaws_internal_free(X);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
laid:
	{
		/* the charts side by side along u; the report */
		double *lo = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 4, (size_t)C + 1)), off = 0, a3 = 0, a2 = 0, dang = 0, gap;
		double bmin[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, bmax[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
		if (!lo)
			goto done;
		if (desc && desc->pack)
		{
			unsigned int* wch = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)W + 1));
			double* bb = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 64, (size_t)C + 1));
			if (wch && bb)
				pm_turn_charts(uv, wch, bb, wedge, chart, n, W, C);
			qaws_internal_free(wch);
			qaws_internal_free(bb);
			if (!wch || !bb)
			{
				qaws_internal_free(lo);
				goto done;
			}
		}
		for (i = 0; i < C; i++)
		{
			lo[4 * (size_t)i] = lo[4 * (size_t)i + 1] = HUGE_VAL;
			lo[4 * (size_t)i + 2] = lo[4 * (size_t)i + 3] = -HUGE_VAL;
		}
		for (i = 0; i < n; i++)
		{
			unsigned int c = chart[i / 3];
			double const* p = uv + 2 * (size_t)wedge[i];
			lo[4 * (size_t)c] = fmin(lo[4 * (size_t)c], p[0]);
			lo[4 * (size_t)c + 1] = fmin(lo[4 * (size_t)c + 1], p[1]);
			lo[4 * (size_t)c + 2] = fmax(lo[4 * (size_t)c + 2], p[0]);
			for (k = 0; k < 3; k++)
			{
				bmin[k] = fmin(bmin[k], x[3 * (size_t)t[i] + k]);
				bmax[k] = fmax(bmax[k], x[3 * (size_t)t[i] + k]);
			}
		}
		gap = 0.02 * sqrt((bmax[0] - bmin[0]) * (bmax[0] - bmin[0]) + (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]) + (bmax[2] - bmin[2]) * (bmax[2] - bmin[2]));
		if (desc && desc->pack)
		{
			/* shelves: the charts by height, tallest first, left to right in rows as wide as the square their boxes fill; then
			   all scaled into the unit square at once (the charts keep their relative sizes) */
			double* hi = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)C + 1));
			unsigned int* ord = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)C + 1));
			double sum = 0, row = 0, wmax = 0, side = 1, sc;
			if (!hi || !ord)
			{
				qaws_internal_free(hi);
				qaws_internal_free(ord);
				qaws_internal_free(lo);
				goto done;
			}
			for (i = 0; i < C; i++)
			{
				hi[2 * (size_t)i] = -HUGE_VAL;
				hi[2 * (size_t)i + 1] = -HUGE_VAL;
				ord[i] = i;
			}
			for (i = 0; i < n; i++)
			{
				unsigned int c = chart[i / 3];
				double const* p = uv + 2 * (size_t)wedge[i];
				hi[2 * (size_t)c] = fmax(hi[2 * (size_t)c], p[0]);
				hi[2 * (size_t)c + 1] = fmax(hi[2 * (size_t)c + 1], p[1]);
			}
			for (i = 0; i < C; i++)
			{
				double w = hi[2 * (size_t)i] - lo[4 * (size_t)i] + gap, hh = hi[2 * (size_t)i + 1] - lo[4 * (size_t)i + 1] + gap;
				hi[2 * (size_t)i] = w;
				hi[2 * (size_t)i + 1] = hh;
				sum += w * hh;
				wmax = fmax(wmax, w);
			}
			for (i = 1; i < C; i++)
				for (j = i; j > 0 && hi[2 * (size_t)ord[j] + 1] > hi[2 * (size_t)ord[j - 1] + 1]; j--)
				{
					unsigned int s = ord[j];
					ord[j] = ord[j - 1];
					ord[j - 1] = s;
				}
			/* by the charts' shapes (pm_pack_bitmap); on a skyline by their boxes when that fails */
			{
				double* xyb = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)C + 1));
				double sideb = 0;
				int packed = xyb && pm_pack_bitmap(uv, wedge, chart, n, W, C, gap, xyb, &sideb);
				for (i = 0; packed && i < C; i++)
				{
					lo[4 * (size_t)i + 2] = xyb[2 * (size_t)i];
					lo[4 * (size_t)i + 3] = xyb[2 * (size_t)i + 1];
				}
				if (packed)
					side = sideb;
				qaws_internal_free(xyb);
				if (!packed)
				/* on a skyline (pm_skyline), the strip's width the one of sqrt(sum) x 1, 1.1, .. 1.6 (at least the widest chart)
				   whose square is least */
				{
					double* sky = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 4, (size_t)C + 2));
					double* xy = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)C + 1));
					double uw, uh, bf = 1;
					if (!sky || !xy)
					{
						qaws_internal_free(sky);
						qaws_internal_free(xy);
						qaws_internal_free(hi);
						qaws_internal_free(ord);
						qaws_internal_free(lo);
						goto done;
					}
					for (k = 0, side = HUGE_VAL; k < 7; k++)
					{
						row = fmax(sqrt(sum) * (1 + 0.1 * k), wmax);
						pm_skyline(hi, ord, C, row, sky, NULL, &uw, &uh);
						if (fmax(uw, uh) < side)
						{
							side = fmax(uw, uh);
							bf = 1 + 0.1 * k;
						}
					}
					pm_skyline(hi, ord, C, fmax(sqrt(sum) * bf, wmax), sky, xy, &uw, &uh);
					for (i = 0; i < C; i++)
					{
						lo[4 * (size_t)i + 2] = xy[2 * (size_t)i] - lo[4 * (size_t)i];
						lo[4 * (size_t)i + 3] = xy[2 * (size_t)i + 1] - lo[4 * (size_t)i + 1];
					}
					side = fmax(uw, uh);
					qaws_internal_free(sky);
					qaws_internal_free(xy);
				}
			}
			sc = side > 0 ? 1 / side : 1;
			for (i = 0; i < n; i++)
			{
				unsigned int c = chart[i / 3];
				out_uv[2 * (size_t)i] = (qaws_scalar)((uv[2 * (size_t)wedge[i]] + lo[4 * (size_t)c + 2]) * sc);
				out_uv[2 * (size_t)i + 1] = (qaws_scalar)((uv[2 * (size_t)wedge[i] + 1] + lo[4 * (size_t)c + 3]) * sc);
			}
			/* the fill: the charts' area in the plane over the square's */
			for (i = 0; i < nt; i++)
			{
				qaws_scalar const *p0 = out_uv + 6 * (size_t)i, *p1 = p0 + 2, *p2 = p0 + 4;
				report->fill += (qaws_scalar)fabs(0.5 * (((double)p1[0] - p0[0]) * ((double)p2[1] - p0[1]) - ((double)p2[0] - p0[0]) * ((double)p1[1] - p0[1])));
			}
			qaws_internal_free(hi);
			qaws_internal_free(ord);
		}
		else
		{
			for (i = 0; i < C; i++)
			{
				double w = lo[4 * (size_t)i + 2] - lo[4 * (size_t)i];
				lo[4 * (size_t)i + 2] = off - lo[4 * (size_t)i];   /* (the chart's shift in u) */
				lo[4 * (size_t)i + 3] = -lo[4 * (size_t)i + 1];
				off += w + gap;
			}
			for (i = 0; i < n; i++)
			{
				unsigned int c = chart[i / 3];
				out_uv[2 * (size_t)i] = (qaws_scalar)(uv[2 * (size_t)wedge[i]] + lo[4 * (size_t)c + 2]);
				out_uv[2 * (size_t)i + 1] = (qaws_scalar)(uv[2 * (size_t)wedge[i] + 1] + lo[4 * (size_t)c + 3]);
			}
		}
		qaws_internal_free(lo);
		/* distortion: areas' shares, corners' angles */
		for (i = 0; i < nt; i++)
		{
			double const *p0 = uv + 2 * (size_t)wedge[3 * i], *p1 = uv + 2 * (size_t)wedge[3 * i + 1], *p2 = uv + 2 * (size_t)wedge[3 * i + 2];
			double X[3][2], au = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
			pm_layout(x, t, i, X);
			a3 += 0.5 * X[1][0] * X[2][1];
			a2 += au;
			report->flipped += au <= 0;
		}
		for (i = 0; i < nt; i++)
		{
			double const *p0 = uv + 2 * (size_t)wedge[3 * i], *p1 = uv + 2 * (size_t)wedge[3 * i + 1], *p2 = uv + 2 * (size_t)wedge[3 * i + 2];
			double X[3][2], au = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
			double const* q[3] = { p0, p1, p2 };
			pm_layout(x, t, i, X);
			if (a2 != 0 && a3 > 0)
				report->area_distortion += (qaws_scalar)fabs(au / a2 - 0.5 * X[1][0] * X[2][1] / a3);
			for (k = 0; k < 3; k++)
				dang += fabs(pm_angle(q[k], q[(k + 1) % 3], q[(k + 2) % 3], 2) - pm_angle(X[k], X[(k + 1) % 3], X[(k + 2) % 3], 2));
		}
		report->angle_distortion = (qaws_scalar)(n ? dang / n : 0);
	}
	if (tchart)
		memcpy(tchart, chart, sizeof(unsigned int) * nt);
	st = QAWS_STATUS_OK;
done:
	sm_system_free(&S);
	sm_csr_free(&K);
	sm_csr_free(&Q);
	qaws_internal_free(pos); qaws_internal_free(x); qaws_internal_free(t); qaws_internal_free(he); qaws_internal_free(sk);
	qaws_internal_free(par); qaws_internal_free(tpar); qaws_internal_free(bord); qaws_internal_free(wedge); qaws_internal_free(chart);
	qaws_internal_free(tw); qaws_internal_free(xw); qaws_internal_free(uv); qaws_internal_free(mass); qaws_internal_free(rhs);
	qaws_internal_free(pin_a); qaws_internal_free(pin_b); qaws_internal_free(T); qaws_internal_free(move); qaws_internal_free(cls);
	return st;
}

/* the triangles' area scales (plane over mesh) within one chart differ by more than max_stretch (or one turned): the chart
   split across its longest extent (the triangles by their centroids' median), the edges between the halves new seams;
   until none does, 32 splits at most */
qaws_status qaws_mesh_3d_parameterize(qaws_scalar* out_uv, qaws_param_report* report, qaws_mesh_3d const* mesh, unsigned int const* seams,
	unsigned int ns, qaws_param_desc const* desc)
{
	unsigned int nt, *tc = NULL, *all = NULL, na = ns, i, k, round;
	double *sc = NULL, *cen = NULL, *rng = NULL;
	unsigned long long* ek = NULL;
	qaws_status st;
	if (!desc || !desc->auto_seams || !(desc->max_stretch > 0) || !mesh || !out_uv || (ns && !seams))
		return pm_parameterize(out_uv, report, mesh, seams, ns, desc, NULL);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	tc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	all = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (2 * (size_t)ns + 6 * (size_t)nt + 2));
	sc = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 4, (size_t)nt + 1));
	cen = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)nt + 1));
	rng = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (size_t)nt + 1));
	ek = (unsigned long long*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned long long), 6, (size_t)nt + 1));
	if (!tc || !all || !sc || !cen || !rng || !ek)
		goto done;
	if (ns)
		memcpy(all, seams, sizeof(unsigned int) * 2 * ns);
	for (i = 0; i < nt; i++)
	{
		unsigned int v[3];
		qaws_vec3 a, b, c;
		qaws_mesh_3d_get_triangle(v, mesh, i);
		a = qaws_mesh_3d_get_vertex(mesh, v[0]);
		b = qaws_mesh_3d_get_vertex(mesh, v[1]);
		c = qaws_mesh_3d_get_vertex(mesh, v[2]);
		{
			double u[3] = { (double)b.x - a.x, (double)b.y - a.y, (double)b.z - a.z }, w[3] = { (double)c.x - a.x, (double)c.y - a.y, (double)c.z - a.z };
			double nx = u[1] * w[2] - u[2] * w[1], ny = u[2] * w[0] - u[0] * w[2], nz = u[0] * w[1] - u[1] * w[0];
			sc[4 * (size_t)i] = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
		}
		cen[3 * (size_t)i] = ((double)a.x + b.x + c.x) / 3;
		cen[3 * (size_t)i + 1] = ((double)a.y + b.y + c.y) / 3;
		cen[3 * (size_t)i + 2] = ((double)a.z + b.z + c.z) / 3;
	}
	for (round = 0;; round++)
	{
		unsigned int worst = PM_NONE, cnt = 0, ax = 0, ne = 0;
		double wr = desc->max_stretch, mid, ext = -1;
		st = pm_parameterize(out_uv, report, mesh, all, na, desc, tc);
		if (st != QAWS_STATUS_OK || round == 32)
			break;
		/* each chart's least and largest scale (sc[4 c + 1], sc[4 c + 2]; turned: -1) */
		for (i = 0; i < report->charts; i++)
		{
			sc[4 * (size_t)i + 1] = HUGE_VAL;
			sc[4 * (size_t)i + 2] = 0;
		}
		for (i = 0; i < nt; i++)
		{
			qaws_scalar const* q = out_uv + 6 * (size_t)i;
			double au = 0.5 * (((double)q[2] - q[0]) * ((double)q[5] - q[1]) - ((double)q[4] - q[0]) * ((double)q[3] - q[1]));
			double s = sc[4 * (size_t)i] > 0 ? au / sc[4 * (size_t)i] : 0;
			sc[4 * (size_t)tc[i] + 1] = fmin(sc[4 * (size_t)tc[i] + 1], s);
			sc[4 * (size_t)tc[i] + 2] = fmax(sc[4 * (size_t)tc[i] + 2], s);
		}
		for (i = 0; i < report->charts; i++)
		{
			double r = sc[4 * (size_t)i + 1] <= 0 ? HUGE_VAL : sc[4 * (size_t)i + 2] / sc[4 * (size_t)i + 1];
			if (r > wr)
			{
				wr = r;
				worst = i;
			}
		}
		if (worst == PM_NONE)
			break;
		/* the worst chart's longest extent, its centroids' median along it */
		for (k = 0; k < 3; k++)
		{
			double lo = HUGE_VAL, hi = -HUGE_VAL;
			for (i = 0; i < nt; i++)
				if (tc[i] == worst)
				{
					lo = fmin(lo, cen[3 * (size_t)i + k]);
					hi = fmax(hi, cen[3 * (size_t)i + k]);
				}
			if (hi - lo > ext)
			{
				ext = hi - lo;
				ax = k;
			}
		}
		for (i = 0; i < nt; i++)
			if (tc[i] == worst)
				rng[cnt++] = cen[3 * (size_t)i + ax];
		qsort(rng, cnt, sizeof *rng, pm_cmp_double);
		mid = rng[cnt / 2];
		/* the edges of the worst chart's triangles keyed with their side; an edge on both sides a seam */
		for (i = 0; i < nt; i++)
			if (tc[i] == worst)
			{
				unsigned int v[3];
				qaws_mesh_3d_get_triangle(v, mesh, i);
				for (k = 0; k < 3; k++)
				{
					ek[2 * (size_t)ne] = pm_key(v[k], v[(k + 1) % 3]);
					ek[2 * (size_t)ne + 1] = (cen[3 * (size_t)i + ax] < mid) | ((unsigned long long)v[k] << 1) | ((unsigned long long)v[(k + 1) % 3] << 33);
					ne++;
				}
			}
		qsort(ek, ne, 2 * sizeof *ek, pm_cmp_u64);
		cnt = 0;
		for (i = 0; i + 1 < ne; i++)
			if (ek[2 * (size_t)i] == ek[2 * (size_t)i + 2] && ((ek[2 * (size_t)i + 1] ^ ek[2 * (size_t)i + 3]) & 1))
			{
				all[2 * (size_t)na] = (unsigned int)(ek[2 * (size_t)i + 1] >> 1);
				all[2 * (size_t)na + 1] = (unsigned int)(ek[2 * (size_t)i + 1] >> 33);
				na++;
				cnt++;
			}
		if (!cnt || 2 * (size_t)na + 6 > 2 * (size_t)ns + 6 * (size_t)nt)
			break;   /* (nothing to cut along, or no room: as it is) */
	}
done:
	qaws_internal_free(tc);
	qaws_internal_free(all);
	qaws_internal_free(sc);
	qaws_internal_free(cen);
	qaws_internal_free(rng);
	qaws_internal_free(ek);
	return st;
}
