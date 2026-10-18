/*
 * Measures of surface meshes: shape (quality, angles, scaled Jacobians,
 * valence against an ideal read from the corners' angles) and the deviation
 * from a reference, the sup of the distance bounded from both sides by
 * branch and bound over triangles.
 *
 * Both directions run one driver: the surface measured (S) and the surface
 * measured against (R) are references, the mesh as its LINEAR one. A
 * triangle of S (a root, then its quarters) has corners projected onto R. Its
 * upper bound comes from a cubic map: its S-domain points to R-domain points
 * matched affinely (through a flat unfolding of R round a vertex shared by
 * the corners' closest triangles, cut where the unfolding changes triangle),
 * the difference of the two patches over each piece a cubic Bezier
 * triangle whose control points (blossoms) bound it by the convex hull
 * property; also the distance's Lipschitz bound from the corners (on a
 * curved S widened by how far its net strays from the corners' plane),
 * which still shrinks where the closest place on R jumps. The lower bound is
 * the largest corner distance projected. Any domain point of R is a point of
 * the surface, so a poor match only widens a bound, never breaks it.
 */

#include "qaws_mesh_internal.h"
#include "../qaws_mesh_measure.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MS_PI 3.14159265358979323846
#define MS_DEG (180.0 / MS_PI)
#define MS_FAN 64

static double ms_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void ms_cross(double const* a, double const* b, double* out)
{
	out[0] = a[1] * b[2] - a[2] * b[1];
	out[1] = a[2] * b[0] - a[0] * b[2];
	out[2] = a[0] * b[1] - a[1] * b[0];
}

/* the area of a patch's flat triangle (its corners) */
static double ms_area(double const (*p)[3])
{
	double u[3], w[3], x[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = p[1][k] - p[0][k];
		w[k] = p[2][k] - p[0][k];
	}
	ms_cross(u, w, x);
	return 0.5 * sqrt(ms_dot(x, x));
}

static double ms_dist(double const* a, double const* b)
{
	double d[3];
	d[0] = a[0] - b[0];
	d[1] = a[1] - b[1];
	d[2] = a[2] - b[2];
	return sqrt(ms_dot(d, d));
}

/* ======================================================================== */
/*  Shape                                                                   */
/* ======================================================================== */

/* the angle at corner c of its face, degrees */
static double ms_angle(qaws_hm const* m, unsigned int c)
{
	double const *p = m->x + 3 * (size_t)m->cv[c], *a = m->x + 3 * (size_t)m->cv[hm_prev(m, c)], *b = m->x + 3 * (size_t)hm_to(m, c);
	double u[3], w[3], x[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = a[k] - p[k];
		w[k] = b[k] - p[k];
	}
	ms_cross(u, w, x);
	return atan2(sqrt(ms_dot(x, x)), ms_dot(u, w)) * MS_DEG;
}

qaws_status qaws_mesh_3d_measure_shape(qaws_mesh_3d const* mesh, qaws_shape_quality* out, unsigned int* out_ideal)
{
	qaws_hm const* m;
	unsigned int f, c, v, ntri = 0, nquad = 0, ne = 0, nsmooth = 0, n567 = 0, nvert = 0, tri_major;
	double qmin = HUGE_VAL, qsum = 0, amin = HUGE_VAL, amax = 0, minsum = 0, below = 0, above = 0, jmin = HUGE_VAL, jsum = 0, ad = 0;
	double esum = 0, evar = 0, mean;
	if (!mesh || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	memset(out, 0, sizeof *out);
	for (f = 0; f < m->nf; f++)
	{
		unsigned int c0 = m->fc[f], n = m->fn[f];
		double lo = HUGE_VAL;
		for (c = c0; c < c0 + n; c++)
		{
			double a = ms_angle(m, c);
			amin = a < amin ? a : amin;
			amax = a > amax ? a : amax;
			lo = a < lo ? a : lo;
			if (n == 3)
			{
				below += a < 30;
				above += a > 90;
			}
			else if (n == 4)
				ad += fabs(a - 90);
		}
		if (n == 3)
		{
			/* Q = 2 sqrt(3) A / (s h): the inradius A / s over the longest edge h */
			double const *p0 = m->x + 3 * (size_t)m->cv[c0], *p1 = m->x + 3 * (size_t)m->cv[c0 + 1], *p2 = m->x + 3 * (size_t)m->cv[c0 + 2];
			double e1[3], e2[3], x[3], la = ms_dist(p0, p1), lb = ms_dist(p1, p2), lc = ms_dist(p2, p0), s, h, q;
			unsigned int k;
			for (k = 0; k < 3; k++)
			{
				e1[k] = p1[k] - p0[k];
				e2[k] = p2[k] - p0[k];
			}
			ms_cross(e1, e2, x);
			s = 0.5 * (la + lb + lc);
			h = la > lb ? (la > lc ? la : lc) : (lb > lc ? lb : lc);
			q = s > 0 && h > 0 ? sqrt(3.0) * sqrt(ms_dot(x, x)) / (s * h) : 0;
			qmin = q < qmin ? q : qmin;
			qsum += q;
			minsum += lo;
			ntri++;
		}
		else if (n == 4)
		{
			/* the smallest corner scaled Jacobian along the unit Newell normal */
			double nn[3], l, j = HUGE_VAL;
			hm_face_normal(m, f, nn);
			l = sqrt(ms_dot(nn, nn));
			for (c = c0; c < c0 + 4; c++)
			{
				double const *p = m->x + 3 * (size_t)m->cv[c], *a = m->x + 3 * (size_t)hm_to(m, c), *b = m->x + 3 * (size_t)m->cv[hm_prev(m, c)];
				double en[3], ep[3], x[3], d, sj;
				unsigned int k;
				for (k = 0; k < 3; k++)
				{
					en[k] = a[k] - p[k];
					ep[k] = b[k] - p[k];
				}
				ms_cross(en, ep, x);
				d = sqrt(ms_dot(en, en) * ms_dot(ep, ep)) * l;
				sj = d > 0 ? ms_dot(x, nn) / d : 0;
				j = sj < j ? sj : j;
			}
			jmin = j < jmin ? j : jmin;
			jsum += j;
			nquad++;
		}
	}
	/* edge lengths: their mean, then their spread */
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			esum += hm_length(m, c);
			ne++;
		}
	mean = ne ? esum / ne : 0;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			double d = hm_length(m, c) - mean;
			evar += d * d;
		}
	/* valences against their ideals */
	tri_major = !(nquad > ntri);
	for (v = 0; v < m->nv; v++)
	{
		unsigned int val, ideal, border, irr;
		if (out_ideal)
			out_ideal[v] = 0;
		if (m->vc[v] == HM_NONE)
			continue;
		nvert++;
		val = hm_valence(m, v);
		border = (unsigned int)hm_is_border_vertex(m, v);
		if (!border && hm_feature_degree(m, v) == 0)
		{
			ideal = tri_major ? 6 : 4;
			irr = val != ideal;
			nsmooth++;
			n567 += val >= 5 && val <= 7;
		}
		else
		{
			/* x: the faces the corners' angles hold at their ideal corner angles; regular within a half (so a tie at a half
			   rounds to the vertex's own count), one face at the least */
			double x = 0;
			unsigned int nf = 0;
			for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
			{
				unsigned int n = m->fn[m->cf[c]];
				x += ms_angle(m, c) / (n == 3 ? 60.0 : n == 4 ? 90.0 : 180.0 * (n - 2) / n);
				nf++;
			}
			ideal = fabs((double)nf - x) <= 0.5 + 1e-9 ? nf : (unsigned int)floor(x + 0.5);
			ideal = (ideal < 1 ? 1 : ideal) + border;
			irr = val != ideal;
		}
		out->irregular += irr;
		out->valence_excess += val > ideal ? val - ideal : ideal - val;
		if (out_ideal)
			out_ideal[v] = ideal;
	}
	out->min_quality = (qaws_scalar)(ntri ? qmin : 0);
	out->mean_quality = (qaws_scalar)(ntri ? qsum / ntri : 0);
	out->min_angle = (qaws_scalar)(m->nf ? amin : 0);
	out->max_angle = (qaws_scalar)amax;
	out->mean_min_angle = (qaws_scalar)(ntri ? minsum / ntri : 0);
	out->below_30 = (qaws_scalar)(ntri ? below / (3.0 * ntri) : 0);
	out->above_90 = (qaws_scalar)(ntri ? above / (3.0 * ntri) : 0);
	out->min_jacobian = (qaws_scalar)(nquad ? jmin : 0);
	out->mean_jacobian = (qaws_scalar)(nquad ? jsum / nquad : 0);
	out->angle_deviation = (qaws_scalar)(nquad ? ad / (4.0 * nquad) : 0);
	out->edge_deviation = (qaws_scalar)(mean > 0 ? sqrt(evar / ne) / mean : 0);
	out->regular_share = (qaws_scalar)(tri_major ? (nsmooth ? (double)n567 / nsmooth : 0) : (nvert ? (double)(nvert - out->irregular) / nvert : 0));
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Bezier nets                                                             */
/* ======================================================================== */

/* The net of a cubic Bezier triangle p (mi_e3 order) over the domain
   triangle of barycentric points W: its blossom at (W0^i, W1^j, W2^k), the
   same order; three de Casteljau levels with the arguments in turn. */
static void ms_blossom(double const (*p)[3], double const (*W)[3], double (*out)[3])
{
	/* the cubic net's points under each quadratic one (110 order: 200 020 002 110 011 101), the quadratic's under each linear one */
	static unsigned char const from3[6][3] = { { 0, 3, 8 }, { 4, 1, 5 }, { 7, 6, 2 }, { 3, 4, 9 }, { 9, 5, 6 }, { 8, 9, 7 } };
	static unsigned char const from2[3][3] = { { 0, 3, 5 }, { 3, 1, 4 }, { 5, 4, 2 } };
	/* each output's arguments, sorted */
	static unsigned char const arg[10][3] = { { 0, 0, 0 }, { 1, 1, 1 }, { 2, 2, 2 }, { 0, 0, 1 }, { 0, 1, 1 }, { 1, 1, 2 }, { 1, 2, 2 }, { 0, 2, 2 },
		{ 0, 0, 2 }, { 0, 1, 2 } };
	double l2[3][6][3], l1[3][3][3][3];
	unsigned int a, b, i, k;
	for (a = 0; a < 3; a++)
		for (i = 0; i < 6; i++)
			for (k = 0; k < 3; k++)
				l2[a][i][k] = W[a][0] * p[from3[i][0]][k] + W[a][1] * p[from3[i][1]][k] + W[a][2] * p[from3[i][2]][k];
	for (a = 0; a < 3; a++)
		for (b = a; b < 3; b++)
			for (i = 0; i < 3; i++)
				for (k = 0; k < 3; k++)
					l1[a][b][i][k] = W[b][0] * l2[a][from2[i][0]][k] + W[b][1] * l2[a][from2[i][1]][k] + W[b][2] * l2[a][from2[i][2]][k];
	for (i = 0; i < 10; i++)
	{
		double(*q)[3] = l1[arg[i][0]][arg[i][1]];
		double const* w = W[arg[i][2]];
		for (k = 0; k < 3; k++)
			out[i][k] = w[0] * q[0][k] + w[1] * q[1][k] + w[2] * q[2][k];
	}
}

/* The largest control point of a cubic net G (of a difference map), the net
   split into its four quarters while above goal, depth levels at most:
   quarters' points are convex combinations of G's, so splitting only
   tightens. */
static double ms_net_max(double const (*G)[3], double goal, unsigned int depth)
{
	static double const quarter[4][3][3] = { { { 1, 0, 0 }, { 0.5, 0.5, 0 }, { 0.5, 0, 0.5 } }, { { 0.5, 0.5, 0 }, { 0, 1, 0 }, { 0, 0.5, 0.5 } },
		{ { 0.5, 0, 0.5 }, { 0, 0.5, 0.5 }, { 0, 0, 1 } }, { { 0.5, 0.5, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 } } };
	double mx = 0, r = 0;
	unsigned int i;
	for (i = 0; i < 10; i++)
	{
		double d = ms_dot(G[i], G[i]);
		mx = d > mx ? d : mx;
	}
	mx = sqrt(mx);
	if (mx <= goal || depth == 0)
		return mx;
	for (i = 0; i < 4; i++)
	{
		double H[10][3], h;
		ms_blossom(G, quarter[i], H);
		h = ms_net_max((double const(*)[3])H, goal, depth - 1);
		r = h > r ? h : r;
	}
	return r < mx ? r : mx;
}

/* ======================================================================== */
/*  The driver                                                              */
/* ======================================================================== */

typedef struct ms_point
{
	double y[3];                    /* a point of S */
	double d;                       /* its distance to R (as found: at least the true one) */
	unsigned int t;                 /* its closest triangle of R; HM_NONE for a vertex S does not use */
	double b[3];                    /* its barycentrics there */
} ms_point;

typedef struct ms_tri
{
	double U;                       /* an upper bound of the distance to R over it */
	double w[3][3];                 /* its corners' barycentrics in its root */
	unsigned int p[3];              /* its corners' points */
	unsigned int root;
} ms_tri;

typedef struct ms_side
{
	qaws_mesh_reference const *S, *R;
	mr_scratch s;
	ms_point* pt;
	unsigned int np, cap;
	ms_tri* heap;                   /* a max-heap on U */
	unsigned int nh, caph;
	unsigned int* splits;           /* per root */
	double lower, worst[3];
	double tol, goal;               /* the gap wanted; lower + tol */
	int failed;
} ms_side;

/* (s, t) into the triangle s >= 0, t >= 0, s + t <= 1 */
static void ms_clamp(double* s, double* t)
{
	if (*s < 0) *s = 0;
	if (*t < 0) *t = 0;
	if (*s + *t > 1)
	{
		double l = *s + *t;
		*s /= l;
		*t /= l;
	}
}

/* the closest point of triangle (A, B, C) to p (Ericson): its barycentrics b, its squared distance */
static double ms_flat_closest(double const* A, double const* B, double const* C, double const* p, double* b)
{
	double ab[3], ac[3], ap[3], bp[3], cp[3], d1, d2, d3, d4, d5, d6, va, vb, vc, v, w, x[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		ab[k] = B[k] - A[k];
		ac[k] = C[k] - A[k];
		ap[k] = p[k] - A[k];
		bp[k] = p[k] - B[k];
		cp[k] = p[k] - C[k];
	}
	d1 = ms_dot(ab, ap); d2 = ms_dot(ac, ap); d3 = ms_dot(ab, bp); d4 = ms_dot(ac, bp); d5 = ms_dot(ab, cp); d6 = ms_dot(ac, cp);
	vc = d1 * d4 - d3 * d2; vb = d5 * d2 - d1 * d6; va = d3 * d6 - d5 * d4;
	if (d1 <= 0 && d2 <= 0) { v = 0; w = 0; }
	else if (d3 >= 0 && d4 <= d3) { v = 1; w = 0; }
	else if (vc <= 0 && d1 >= 0 && d3 <= 0) { v = d1 / (d1 - d3); w = 0; }
	else if (d6 >= 0 && d5 <= d6) { v = 0; w = 1; }
	else if (vb <= 0 && d2 >= 0 && d6 <= 0) { v = 0; w = d2 / (d2 - d6); }
	else if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); v = 1 - w; }
	else { double den = 1 / (va + vb + vc); v = vb * den; w = vc * den; }
	b[0] = 1 - v - w;
	b[1] = v;
	b[2] = w;
	for (k = 0; k < 3; k++)
		x[k] = A[k] + v * ab[k] + w * ac[k] - p[k];
	return ms_dot(x, x);
}

/* the exact closest point of a flat R: the grid's rings outward, Ericson per triangle, until a ring's lower bound passes the best */
typedef struct ms_query
{
	qaws_mesh_reference const* R;
	mr_scratch* s;
	double p[3], best, b[3];
	unsigned int t;
} ms_query;

static void ms_visit(void* user, unsigned int t)
{
	ms_query* q = (ms_query*)user;
	qaws_bp_box const* x = &q->R->box[t];
	double d = 0, b[3];
	unsigned int k;
	if (q->s->tried[t] == q->s->gen)
		return;
	q->s->tried[t] = q->s->gen;
	for (k = 0; k < 3; k++)
	{
		double e = q->p[k] < x->lo[k] ? x->lo[k] - q->p[k] : q->p[k] > x->hi[k] ? q->p[k] - x->hi[k] : 0;
		d += e * e;
	}
	if (d >= q->best)
		return;
	d = ms_flat_closest(q->R->patch[t].p[0], q->R->patch[t].p[1], q->R->patch[t].p[2], q->p, b);
	if (d < q->best)
	{
		q->best = d;
		q->t = t;
		memcpy(q->b, b, sizeof b);
	}
}

static double ms_flat_project(qaws_mesh_reference const* R, mr_scratch* s, double const* y, unsigned int* out_t, double* out_b)
{
	ms_query q;
	unsigned int ring;
	q.R = R;
	q.s = s;
	memcpy(q.p, y, sizeof q.p);
	q.best = HUGE_VAL;
	q.t = 0;
	q.b[0] = 1;
	q.b[1] = q.b[2] = 0;
	if (++s->gen == 0)
	{
		memset(s->tried, 0, sizeof(unsigned int) * R->nt);
		s->gen = 1;
	}
	for (ring = 0;; ring++)
	{
		double lb = qaws_internal_grid_ring(R->grid, y, ring, ms_visit, &q);
		if (lb == HUGE_VAL || lb * lb > q.best)
			break;
	}
	*out_t = q.t;
	memcpy(out_b, q.b, sizeof q.b);
	return q.best;
}

/* a cubic patch at barycentrics w by de Casteljau, its derivatives along a = w1 and b = w2 (xa, xb may be NULL; mi_eval's
   power form, with its second derivatives, costs several times as much in the projections' inner loop) */
static void ms_eval(double const (*p)[3], double const* w, double* x, double* xa, double* xb)
{
	static unsigned char const from3[6][3] = { { 0, 3, 8 }, { 4, 1, 5 }, { 7, 6, 2 }, { 3, 4, 9 }, { 9, 5, 6 }, { 8, 9, 7 } };
	static unsigned char const from2[3][3] = { { 0, 3, 5 }, { 3, 1, 4 }, { 5, 4, 2 } };
	double l2[6][3], l1[3][3];
	unsigned int i, k;
	for (i = 0; i < 6; i++)
		for (k = 0; k < 3; k++)
			l2[i][k] = w[0] * p[from3[i][0]][k] + w[1] * p[from3[i][1]][k] + w[2] * p[from3[i][2]][k];
	for (i = 0; i < 3; i++)
		for (k = 0; k < 3; k++)
			l1[i][k] = w[0] * l2[from2[i][0]][k] + w[1] * l2[from2[i][1]][k] + w[2] * l2[from2[i][2]][k];
	for (k = 0; k < 3; k++)
	{
		x[k] = w[0] * l1[0][k] + w[1] * l1[1][k] + w[2] * l1[2][k];
		if (xa)
			xa[k] = 3 * (l1[1][k] - l1[0][k]);
		if (xb)
			xb[k] = 3 * (l1[2][k] - l1[0][k]);
	}
}

/* Descent on R's curved patch P to y from barycentrics b (in, out), its squared distance: projected Gauss-Newton whose step, at
   a side the gradient leads out of, runs along that side alone (the active set; the reference's plain projected step stalls
   there short of the side's closest point), halved while it does not improve. */
static double ms_descend(mi_tri const* P, double const* y, double* b)
{
	double const(*p)[3] = (double const(*)[3])P->p;
	double x[3], xa[3], xb[3], best;
	unsigned int it, k;
	ms_eval(p, b, x, xa, xb);
	best = (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) + (x[2] - y[2]) * (x[2] - y[2]);
	for (it = 0; it < 32; it++)
	{
		double r[3], g0, g1, da, db, s, a = b[1], c = b[2];
		int improved = 0, out0, out1, out2;
		for (k = 0; k < 3; k++)
			r[k] = x[k] - y[k];
		g0 = ms_dot(xa, r);
		g1 = ms_dot(xb, r);
		out0 = a <= 0 && g0 > 0;
		out1 = c <= 0 && g1 > 0;
		out2 = b[0] <= 0 && g0 + g1 < 0;
		if (out0 + out1 + out2 > 1)
			break;   /* (a corner the descent leads out of both ways) */
		if (out0 || out1 || out2)
		{
			/* along the side: a fixed at 0, c at 0, or a + c at 1 */
			double e0 = out0 ? 0 : 1, e1 = out1 ? 0 : out2 ? -1 : 1, xe[3], h;
			for (k = 0; k < 3; k++)
				xe[k] = e0 * xa[k] + e1 * xb[k];
			h = ms_dot(xe, xe);
			if (h <= 1e-300)
				break;
			s = -ms_dot(xe, r) / h;
			da = s * e0;
			db = s * e1;
		}
		else
		{
			double h00 = ms_dot(xa, xa), h01 = ms_dot(xa, xb), h11 = ms_dot(xb, xb), det = h00 * h11 - h01 * h01;
			if (det <= 1e-300)
				break;
			da = -(h11 * g0 - h01 * g1) / det;
			db = -(h00 * g1 - h01 * g0) / det;
		}
		if (fabs(da) + fabs(db) <= 1e-14)
			break;   /* (converged) */
		for (s = 1; s > 1e-4; s *= 0.5)
		{
			double w[3], z[3], d;
			w[1] = a + s * da;
			w[2] = c + s * db;
			ms_clamp(&w[1], &w[2]);
			w[0] = 1 - w[1] - w[2];
			ms_eval(p, w, z, NULL, NULL);
			d = (z[0] - y[0]) * (z[0] - y[0]) + (z[1] - y[1]) * (z[1] - y[1]) + (z[2] - y[2]) * (z[2] - y[2]);
			if (d < best)
			{
				improved = fabs(w[1] - a) + fabs(w[2] - c) > 1e-14;
				memcpy(b, w, sizeof w);
				best = d;
				break;
			}
		}
		if (!improved)
			break;
		ms_eval(p, b, x, xa, xb);
	}
	return best;
}

/* the closest point of R's patch t to y: its barycentrics b, its squared distance; LINEAR exact (Ericson), curved the descent from
   the flat triangle's */
static double ms_closest(qaws_mesh_reference const* R, unsigned int t, double const* y, double* b)
{
	double d = ms_flat_closest(R->patch[t].p[0], R->patch[t].p[1], R->patch[t].p[2], y, b);
	return R->kind == QAWS_INTERPOLANT_LINEAR ? d : ms_descend(&R->patch[t], y, b);
}

/* The closest point of R's patch *t to y, its barycentrics b, its squared distance. LINEAR: exact (Ericson). Curved: the
   descent from b (in), then on across the side it ends on into the neighbour while that is closer (*t follows). */
static double ms_polish(qaws_mesh_reference const* R, unsigned int* t, double const* y, double* b)
{
	double best;
	unsigned int hop, i;
	if (R->kind == QAWS_INTERPOLANT_LINEAR)
		return ms_closest(R, *t, y, b);
	best = ms_descend(&R->patch[*t], y, b);
	for (hop = 0; hop < 4; hop++)
	{
		double w[3], d;
		unsigned int e, nb, u, k;
		for (i = 0; i < 3 && b[i] > 1e-12; i++)
			;
		if (i == 3)
			break;
		/* the side from corner e to e + 1 (b[i] = 0 there); across it corner k is our e + 1, k + 1 our e (the weights summing to 1) */
		e = (i + 1) % 3;
		nb = R->ta[3 * *t + e];
		if (nb == HM_NONE)
			break;
		u = nb / 3;
		k = nb % 3;
		w[k] = b[(e + 1) % 3];
		w[(k + 1) % 3] = 1 - w[k];
		w[(k + 2) % 3] = 0;
		d = ms_descend(&R->patch[u], y, w);
		if (!(d < best))
			break;
		best = d;
		*t = u;
		memcpy(b, w, sizeof w);
	}
	return best;
}

/* A point of S projected onto R (project 0: a vertex no triangle uses, kept out of the bounds). First onto the hint triangles
   (HM_NONE: none; on a flat R also their neighbours): a distance there no larger than the lower bound cannot raise it, and
   any point of R serves the upper bounds. Else the global search (the grid's rings over all of R): a hinted one could stop
   at a sheet nearby, too far for a lower bound. */
static unsigned int ms_point_add(ms_side* q, double const* y, int project, unsigned int h0, unsigned int h1)
{
	qaws_mesh_reference const* R = q->R;
	ms_point* p;
	qaws_mesh_location l;
	unsigned int h, e, hint[2];
	double best = HUGE_VAL;
	if (q->np == q->cap)
	{
		unsigned int cap = q->cap ? 2 * q->cap : 1024;
		ms_point* g = (ms_point*)realloc(q->pt, sizeof(ms_point) * cap);
		if (!g)
		{
			q->failed = 1;
			return HM_NONE;
		}
		q->pt = g;
		q->cap = cap;
	}
	p = &q->pt[q->np];
	memcpy(p->y, y, sizeof p->y);
	p->t = HM_NONE;
	p->b[0] = p->b[1] = p->b[2] = 0;
	p->d = 0;
	if (!project)
		return q->np++;
	hint[0] = h0;
	hint[1] = h1 != h0 ? h1 : HM_NONE;
	for (h = 0; h < 2; h++)
		for (e = 0; e < (R->kind == QAWS_INTERPOLANT_LINEAR ? 4u : 1u) && hint[h] < R->nt; e++)
		{
			unsigned int t = e == 0 ? hint[h] : R->ta[3 * hint[h] + e - 1];
			double b[3], d;
			if (t == HM_NONE)
				continue;
			t = e == 0 ? t : t / 3;
			d = ms_closest(R, t, y, b);
			if (d < best)
			{
				best = d;
				p->t = t;
				memcpy(p->b, b, sizeof b);
			}
		}
	p->d = sqrt(best);
	if (p->d <= q->lower && p->b[0] > 1e-12 && p->b[1] > 1e-12 && p->b[2] > 1e-12)
		return q->np++;   /* (inside its triangle: not clamped short of a closer neighbour) */
	/* the global search: on a curved R the reference's, polished from its location (read back from qaws_scalar into the
	   triangle) since its patch projection can stall on a side; the closest point found kept */
	{
		unsigned int t = HM_NONE;
		double b[3], d = HUGE_VAL;
		if (R->kind == QAWS_INTERPOLANT_LINEAR)
			d = ms_flat_project(R, &q->s, y, &t, b);
		else if (mr_project(R, y, HM_NONE, &q->s, &l) < HUGE_VAL && l.patch < R->nt)
		{
			t = l.patch;
			b[1] = (double)l.a;
			b[2] = (double)l.b;
			ms_clamp(&b[1], &b[2]);
			b[0] = 1 - b[1] - b[2];
			d = ms_polish(R, &t, y, b);
		}
		if (d < best)
		{
			best = d;
			p->t = t;
			memcpy(p->b, b, sizeof b);
		}
	}
	p->d = sqrt(best);
	if (p->d > q->lower)
	{
		q->lower = p->d;
		q->goal = q->lower + q->tol;
		memcpy(q->worst, y, sizeof q->worst);
	}
	return q->np++;
}

/* the bound over S's triangle W against R's triangle t at barycentrics u: the net of S at W minus R's at u, refined while above goal */
static double ms_pair(ms_side const* q, unsigned int root, double const (*W)[3], unsigned int t, double const (*u)[3], double goal)
{
	double A[10][3], B[10][3];
	unsigned int i, k;
	ms_blossom((double const(*)[3])q->S->patch[root].p, W, A);
	ms_blossom((double const(*)[3])q->R->patch[t].p, u, B);
	for (i = 0; i < 10; i++)
		for (k = 0; k < 3; k++)
			A[i][k] -= B[i][k];
	return ms_net_max((double const(*)[3])A, goal, 2);
}

/* ---- the unfolding of R round a vertex ---------------------------------- */

typedef struct ms_fan
{
	unsigned int n, closed;
	unsigned int t[MS_FAN];         /* the triangles round the centre, counter-clockwise */
	unsigned char k[MS_FAN];        /* the centre's corner in each */
	double th[MS_FAN], al[MS_FAN];  /* each one's start angle in the chart, its angle at the centre */
	double a[MS_FAN][2], b[MS_FAN][2];   /* chart points of its corners k + 1 and k + 2 */
} ms_fan;

/* the triangles of R round its vertex v, from triangle t (one of them); 0 when they are more than MS_FAN */
static int ms_fan_walk(qaws_mesh_reference const* R, unsigned int t, unsigned int v, ms_fan* F)
{
	unsigned int ring[2 * MS_FAN], k, n = 0, cw = 0, u = t, j;
	for (k = 0; k < 3 && R->tv[3 * t + k] != v; k++)
		;
	if (k == 3)
		return 0;
	F->closed = 0;
	/* counter-clockwise: across the edge into the centre (k + 2) */
	for (;;)
	{
		unsigned int nb;
		if (n == MS_FAN)
			return 0;
		ring[n++] = 3 * u + k;
		nb = R->ta[3 * u + (k + 2) % 3];
		if (nb == HM_NONE)
			break;
		u = nb / 3;
		k = nb % 3;
		if (u == t)
		{
			F->closed = 1;
			break;
		}
	}
	if (!F->closed)
	{
		/* clockwise from t: across the edge out of the centre (k), the centre at the neighbour's next corner */
		unsigned int c = ring[0];
		for (;;)
		{
			unsigned int nb = R->ta[c];
			if (nb == HM_NONE)
				break;
			if (n + cw == MS_FAN)
				return 0;
			c = 3 * (nb / 3) + (nb % 3 + 1) % 3;
			ring[MS_FAN + cw++] = c;
		}
	}
	F->n = n + cw;
	for (j = 0; j < cw; j++)
	{
		F->t[j] = ring[MS_FAN + cw - 1 - j] / 3;
		F->k[j] = (unsigned char)(ring[MS_FAN + cw - 1 - j] % 3);
	}
	for (j = 0; j < n; j++)
	{
		F->t[cw + j] = ring[j] / 3;
		F->k[cw + j] = (unsigned char)(ring[j] % 3);
	}
	return 1;
}

/* the fan laid flat from triangle s outward both ways, isometrically */
static void ms_fan_unfold(qaws_mesh_reference const* R, ms_fan* F, unsigned int s)
{
	unsigned int j, n = F->n, ccw = F->closed ? (n - 1) / 2 + 1 : n - s, cw = F->closed ? n - ccw : s;
	double la[MS_FAN], lb[MS_FAN];
	for (j = 0; j < n; j++)
	{
		double const *P = R->patch[F->t[j]].p[F->k[j]], *A = R->patch[F->t[j]].p[(F->k[j] + 1) % 3], *B = R->patch[F->t[j]].p[(F->k[j] + 2) % 3];
		double u[3], w[3], x[3];
		unsigned int k;
		for (k = 0; k < 3; k++)
		{
			u[k] = A[k] - P[k];
			w[k] = B[k] - P[k];
		}
		ms_cross(u, w, x);
		F->al[j] = atan2(sqrt(ms_dot(x, x)), ms_dot(u, w));
		la[j] = sqrt(ms_dot(u, u));
		lb[j] = sqrt(ms_dot(w, w));
	}
	F->th[s] = 0;
	for (j = 1; j < ccw; j++)
		F->th[(s + j) % n] = F->th[(s + j - 1) % n] + F->al[(s + j - 1) % n];
	for (j = 1; j <= cw; j++)
		F->th[(s + n - j) % n] = F->th[(s + n - j + 1) % n] - F->al[(s + n - j) % n];
	for (j = 0; j < n; j++)
	{
		double e = F->th[j] + F->al[j];
		F->a[j][0] = la[j] * cos(F->th[j]);
		F->a[j][1] = la[j] * sin(F->th[j]);
		F->b[j][0] = lb[j] * cos(e);
		F->b[j][1] = lb[j] * sin(e);
	}
}

/* a point of R (triangle t, barycentrics b) as a fan triangle and its barycentrics there (in its corner order): t itself, or one
   holding the point's support (its corners of weight above 1e-12: a vertex or an edge); 0 when none does */
static int ms_fan_find(qaws_mesh_reference const* R, ms_fan const* F, unsigned int t, double const* b, unsigned int* out_j, double* out_b)
{
	unsigned int j, i, k;
	for (j = 0; j < F->n; j++)
		if (F->t[j] == t)
		{
			*out_j = j;
			memcpy(out_b, b, sizeof(double) * 3);
			return 1;
		}
	for (j = 0; j < F->n; j++)
	{
		double s = 0;
		int ok = 1;
		out_b[0] = out_b[1] = out_b[2] = 0;
		for (i = 0; i < 3 && ok; i++)
		{
			if (b[i] <= 1e-12)
				continue;
			for (k = 0; k < 3 && R->tv[3 * F->t[j] + k] != R->tv[3 * t + i]; k++)
				;
			if (k == 3)
				ok = 0;
			else
			{
				out_b[k] = b[i];
				s += b[i];
			}
		}
		if (ok && s > 0)
		{
			out_b[0] /= s;
			out_b[1] /= s;
			out_b[2] /= s;
			*out_j = j;
			return 1;
		}
	}
	return 0;
}

typedef struct ms_vtx
{
	double c[2], l[3];              /* chart point; barycentrics in the triangle measured */
} ms_vtx;

/* the part of polygon in on the side of the line through the origin along d where side * cross(d, p) >= 0 */
static unsigned int ms_clip(ms_vtx const* in, unsigned int n, double const* d, double side, ms_vtx* out)
{
	unsigned int i, m = 0, k;
	for (i = 0; i < n; i++)
	{
		ms_vtx const *a = &in[i], *b = &in[(i + 1) % n];
		double sa = side * (d[0] * a->c[1] - d[1] * a->c[0]), sb = side * (d[0] * b->c[1] - d[1] * b->c[0]);
		if (sa >= 0)
			out[m++] = *a;
		if ((sa >= 0) != (sb >= 0))
		{
			double s = sa / (sa - sb);
			for (k = 0; k < 2; k++)
				out[m].c[k] = a->c[k] + s * (b->c[k] - a->c[k]);
			for (k = 0; k < 3; k++)
				out[m].l[k] = a->l[k] + s * (b->l[k] - a->l[k]);
			m++;
		}
	}
	return m;
}

/* The bound through the unfolding round centre v: the triangle's chart image
   cut by the wedges of the fan's triangles (the far side, where the fan
   does not close or overlaps itself, split between its two ends), each
   piece matched affinely to its triangle. HUGE_VAL when a corner is not
   round v or the image is flat. */
static double ms_chart(ms_side const* q, ms_tri const* T, ms_point const* const* P, unsigned int v)
{
	qaws_mesh_reference const* R = q->R;
	ms_fan F;
	ms_vtx tau[3];
	unsigned int j, k, jj[3], nw = 0, wj[2 * MS_FAN + 8];
	double bb[3][3], lo[2 * MS_FAN + 8], hi[2 * MS_FAN + 8], mid, area, size = 0, U = 0, cover = 0;
	for (k = 0; k < 3 && R->tv[3 * P[k]->t] != v && R->tv[3 * P[k]->t + 1] != v && R->tv[3 * P[k]->t + 2] != v; k++)
		;
	if (k == 3 || !ms_fan_walk(R, P[k]->t, v, &F))
		return HUGE_VAL;
	for (k = 0; k < 3; k++)
		if (!ms_fan_find(R, &F, P[k]->t, P[k]->b, &jj[k], bb[k]))
			return HUGE_VAL;
	ms_fan_unfold(R, &F, jj[0]);
	for (k = 0; k < 3; k++)
	{
		unsigned int c = F.k[jj[k]];
		double s = bb[k][(c + 1) % 3], t = bb[k][(c + 2) % 3];
		tau[k].c[0] = s * F.a[jj[k]][0] + t * F.b[jj[k]][0];
		tau[k].c[1] = s * F.a[jj[k]][1] + t * F.b[jj[k]][1];
		tau[k].l[0] = k == 0;
		tau[k].l[1] = k == 1;
		tau[k].l[2] = k == 2;
	}
	for (k = 0; k < 3; k++)
	{
		double dx = tau[(k + 1) % 3].c[0] - tau[k].c[0], dy = tau[(k + 1) % 3].c[1] - tau[k].c[1];
		size = dx * dx + dy * dy > size ? dx * dx + dy * dy : size;
	}
	area = (tau[1].c[0] - tau[0].c[0]) * (tau[2].c[1] - tau[0].c[1]) - (tau[1].c[1] - tau[0].c[1]) * (tau[2].c[0] - tau[0].c[0]);
	if (!(fabs(area) > 1e-12 * size))
		return HUGE_VAL;
	/* the wedges in angular order, clipped to [mid - 2 pi, mid] and stretched to it at the ends, wide ones halved; a zero-area
	   triangle's wedge maps nothing (every point to the centre), its range goes to its neighbours */
	{
		unsigned int first = F.closed ? (jj[0] + (F.n - 1) / 2 + 1) % F.n : 0, i, last = F.n;
		double a0 = F.th[first], a1 = F.th[(first + F.n - 1) % F.n] + F.al[(first + F.n - 1) % F.n];
		mid = 0.5 * (a0 + a1 + 2 * MS_PI);
		for (i = 0; i < F.n; i++)
			if (sin(F.al[(first + i) % F.n]) > 1e-12)
				last = i;
		for (i = 0; i < F.n; i++)
		{
			unsigned int jw = (first + i) % F.n;
			double l = F.th[jw] > mid - 2 * MS_PI ? F.th[jw] : mid - 2 * MS_PI, h = F.th[jw] + F.al[jw] < mid ? F.th[jw] + F.al[jw] : mid;
			if (!(sin(F.al[jw]) > 1e-12))
				continue;
			if (nw == 0)
				l = mid - 2 * MS_PI;
			if (i == last)
				h = mid;
			while (h > l)
			{
				double e = h - l > 3 ? l + 0.5 * (h - l < 6 ? h - l : 6) : h;
				if (nw == 2 * MS_FAN + 8)
					return HUGE_VAL;
				lo[nw] = l;
				hi[nw] = e;
				wj[nw++] = jw;
				l = e;
			}
		}
	}
	for (j = 0; j < nw; j++)
	{
		ms_vtx p1[8], p2[8];
		double d0[2], d1[2];
		unsigned int n, i, t = F.t[wj[j]], c = F.k[wj[j]];
		double const *A = F.a[wj[j]], *B = F.b[wj[j]], det = A[0] * B[1] - A[1] * B[0];
		d0[0] = cos(lo[j]);
		d0[1] = sin(lo[j]);
		d1[0] = cos(hi[j]);
		d1[1] = sin(hi[j]);
		n = ms_clip(tau, 3, d0, 1, p1);
		n = n ? ms_clip(p1, n, d1, -1, p2) : 0;
		for (i = 1; i + 1 < n; i++)
		{
			ms_vtx const* v3[3];
			double W[3][3], u[3][3], b, piece;
			unsigned int e, z;
			v3[0] = &p2[0];
			v3[1] = &p2[i];
			v3[2] = &p2[i + 1];
			piece = fabs(v3[0]->l[0] * (v3[1]->l[1] * v3[2]->l[2] - v3[1]->l[2] * v3[2]->l[1])
				- v3[0]->l[1] * (v3[1]->l[0] * v3[2]->l[2] - v3[1]->l[2] * v3[2]->l[0])
				+ v3[0]->l[2] * (v3[1]->l[0] * v3[2]->l[1] - v3[1]->l[1] * v3[2]->l[0]));
			if (piece == 0)
				continue;   /* (a sliver on a ray: the closed pieces round it hold it) */
			cover += piece;
			for (e = 0; e < 3; e++)
			{
				double s = (v3[e]->c[0] * B[1] - v3[e]->c[1] * B[0]) / det, r = (A[0] * v3[e]->c[1] - A[1] * v3[e]->c[0]) / det;
				for (z = 0; z < 3; z++)
					W[e][z] = v3[e]->l[0] * T->w[0][z] + v3[e]->l[1] * T->w[1][z] + v3[e]->l[2] * T->w[2][z];
				ms_clamp(&s, &r);
				u[e][c] = 1 - s - r;
				u[e][(c + 1) % 3] = s;
				u[e][(c + 2) % 3] = r;
			}
			b = ms_pair(q, T->root, (double const(*)[3])W, t, (double const(*)[3])u, U > q->goal ? U : q->goal);   /* (no use refining a piece below the max so far) */
			U = b > U ? b : U;
		}
	}
	/* the pieces tile the triangle (their barycentric areas sum to 1), or the bound is not taken */
	return fabs(cover - 1) <= 1e-9 ? U : HUGE_VAL;
}

/* the bound over T against one triangle of R at a time, its corners' closest points there */
static double ms_candidates(ms_side const* q, ms_tri const* T, ms_point const* const* P)
{
	double U = HUGE_VAL;
	unsigned int i, k;
	for (i = 0; i < 3; i++)
	{
		unsigned int t = P[i]->t;
		double u[3][3], b;
		if ((i > 0 && t == P[0]->t) || (i > 1 && t == P[1]->t))
			continue;
		for (k = 0; k < 3; k++)
			if (P[k]->t == t)
				memcpy(u[k], P[k]->b, sizeof u[k]);
			else
				ms_closest(q->R, t, P[k]->y, u[k]);
		b = ms_pair(q, T->root, (double const(*)[3])T->w, t, (double const(*)[3])u, q->goal);
		U = b < U ? b : U;
	}
	return U;
}

/* the distance's Lipschitz bound over a flat triangle: from the corners, out to the farthest point from them all (the
   circumcentre, or the longest edge's midpoint past a right angle), or each corner's farthest corner */
static double ms_lipschitz(ms_point const* const* P)
{
	double e0 = ms_dist(P[0]->y, P[1]->y), e1 = ms_dist(P[1]->y, P[2]->y), e2 = ms_dist(P[2]->y, P[0]->y), u[3], w[3], x[3], A2, g, h, dm, ua, ub;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = P[1]->y[k] - P[0]->y[k];
		w[k] = P[2]->y[k] - P[0]->y[k];
	}
	ms_cross(u, w, x);
	A2 = sqrt(ms_dot(x, x));
	h = e0 > e1 ? (e0 > e2 ? e0 : e2) : (e1 > e2 ? e1 : e2);
	if (2 * h * h >= e0 * e0 + e1 * e1 + e2 * e2 || A2 <= 0)
		g = 0.5 * h;
	else
		g = e0 * e1 * e2 / (2 * A2);
	dm = P[0]->d > P[1]->d ? (P[0]->d > P[2]->d ? P[0]->d : P[2]->d) : (P[1]->d > P[2]->d ? P[1]->d : P[2]->d);
	ua = dm + g;
	ub = P[0]->d + (e0 > e2 ? e0 : e2);
	ub = fmin(ub, P[1]->d + (e0 > e1 ? e0 : e1));
	ub = fmin(ub, P[2]->d + (e1 > e2 ? e1 : e2));
	return ua < ub ? ua : ub;
}

static double ms_bound(ms_side const* q, ms_tri const* T)
{
	ms_point const* P[3];
	double U;
	unsigned int k;
	for (k = 0; k < 3; k++)
		P[k] = &q->pt[T->p[k]];
	/* Lipschitz from the corners: over their flat triangle, on a curved S plus the farthest the sub-patch's net strays from that
	   triangle's (degree-elevated) net, since a point of the sub-patch and its flat twin are the same combination of the two
	   nets. It holds where the closest place on R jumps (the matched maps below do not shrink there), and on any R. */
	U = ms_lipschitz(P);
	if (q->S->kind != QAWS_INTERPOLANT_LINEAR)
	{
		static unsigned char const pw[10][3] = { { 3, 0, 0 }, { 0, 3, 0 }, { 0, 0, 3 }, { 2, 1, 0 }, { 1, 2, 0 }, { 0, 2, 1 }, { 0, 1, 2 }, { 1, 0, 2 },
			{ 2, 0, 1 }, { 1, 1, 1 } };
		double N[10][3], e = 0;
		unsigned int i;
		ms_blossom((double const(*)[3])q->S->patch[T->root].p, T->w, N);
		for (i = 0; i < 10; i++)
		{
			double d[3];
			for (k = 0; k < 3; k++)
				d[k] = N[i][k] - (pw[i][0] * P[0]->y[k] + pw[i][1] * P[1]->y[k] + pw[i][2] * P[2]->y[k]) / 3;
			e = fmax(e, ms_dot(d, d));
		}
		U += sqrt(e);
	}
	if (U <= q->goal)
		return U;
	if (P[0]->t == P[1]->t && P[1]->t == P[2]->t)
	{
		double u[3][3];
		for (k = 0; k < 3; k++)
			memcpy(u[k], P[k]->b, sizeof u[k]);
		return fmin(U, ms_pair(q, T->root, (double const(*)[3])T->w, P[0]->t, (double const(*)[3])u, q->goal));
	}
	{
		/* centres: the vertices of the corners' triangles, those most of them hold first; each chart that applies, until one
		   is within the goal (a poor match, as round a zero-area triangle, only widens its own) */
		unsigned int cand[9], score[9], nc = 0, i, j;
		double b = HUGE_VAL;
		for (k = 0; k < 3; k++)
			for (i = 0; i < 3; i++)
			{
				unsigned int v = q->R->tv[3 * P[k]->t + i];
				for (j = 0; j < nc && cand[j] != v; j++)
					;
				if (j == nc)
				{
					cand[nc] = v;
					score[nc++] = 0;
				}
				score[j]++;
			}
		for (i = 0; i < nc && b > q->goal; i++)
		{
			unsigned int best = 0;
			for (j = 1; j < nc; j++)
				if (score[j] > score[best])
					best = j;
			if (score[best] == 0)
				break;
			score[best] = 0;
			/* (a corner off the centre's triangles must sit on an edge or a vertex to be read in its fan) */
			for (k = 0; k < 3; k++)
			{
				unsigned int const* tv = q->R->tv + 3 * P[k]->t;
				if (tv[0] != cand[best] && tv[1] != cand[best] && tv[2] != cand[best] && P[k]->b[0] > 1e-12 && P[k]->b[1] > 1e-12 && P[k]->b[2] > 1e-12)
					break;
			}
			if (k == 3)
				b = fmin(b, ms_chart(q, T, P, cand[best]));
		}
		/* shortcut: corners that share no vertex star (a triangle over a ring of R's, or coarser) are matched to one of their
		   triangles at a time, which converges only linearly until the splits bring them into a star; a chart round a
		   triangle and its neighbours would cut a level where R is much finer than S */
		if (b > q->goal)
			b = fmin(b, ms_candidates(q, T, P));
		return fmin(U, b);
	}
}

static int ms_push(ms_side* q, ms_tri const* T)
{
	unsigned int i;
	if (q->nh == q->caph)
	{
		unsigned int cap = q->caph ? 2 * q->caph : 1024;
		ms_tri* g = (ms_tri*)realloc(q->heap, sizeof(ms_tri) * cap);
		if (!g)
		{
			q->failed = 1;
			return 0;
		}
		q->heap = g;
		q->caph = cap;
	}
	for (i = q->nh++; i > 0 && q->heap[(i - 1) / 2].U < T->U; i = (i - 1) / 2)
		q->heap[i] = q->heap[(i - 1) / 2];
	q->heap[i] = *T;
	return 1;
}

static void ms_pop(ms_side* q, ms_tri* out)
{
	ms_tri last;
	unsigned int i = 0;
	*out = q->heap[0];
	last = q->heap[--q->nh];
	for (;;)
	{
		unsigned int c = 2 * i + 1;
		if (c >= q->nh)
			break;
		if (c + 1 < q->nh && q->heap[c + 1].U > q->heap[c].U)
			c++;
		if (q->heap[c].U <= last.U)
			break;
		q->heap[i] = q->heap[c];
		i = c;
	}
	if (q->nh)
		q->heap[i] = last;
}

/* the sup over S of the distance to R: a lower and an upper bound; 1 when within tol of each other */
static qaws_status ms_run(qaws_mesh_reference const* S, qaws_mesh_reference const* R, double tol, unsigned int budget, double* out_rms,
	double* bounds, double* worst, int* converged)
{
	ms_side q;
	qaws_hm const* m = &S->mesh->hm;
	unsigned int t, v, k;
	double aside = 0, upper, scale = 0, area = 0;
	for (t = 0; t < R->nt; t++)
		area += ms_area(R->patch[t].p) / R->nt;
	area = area > 0 ? area : HUGE_VAL;
	memset(&q, 0, sizeof q);
	q.S = S;
	q.R = R;
	q.tol = tol;
	q.goal = tol;
	q.splits = (unsigned int*)calloc((size_t)S->nt + 1, sizeof(unsigned int));
	if (!q.splits || !mr_scratch_init(&q.s, R->nt))
	{
		free(q.splits);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* S's vertices are its first points; those of no triangle are left out of the bounds */
	for (v = 0; v < m->nv && !q.failed; v++)
		ms_point_add(&q, m->x + 3 * (size_t)v, S->vloc[v] != HM_NONE, m->foot && m->foot[v].patch < R->nt ? m->foot[v].patch : HM_NONE, HM_NONE);
	if (out_rms && !q.failed)
	{
		/* rms: per triangle the 7-point degree-5 rule (Dunavant) of the squared distance, projected from the nearest corner's
		   closest triangle (not certified) */
		static double const g[7][4] = { { 0.225, 1.0 / 3, 1.0 / 3, 1.0 / 3 },
			{ 0.132394152788506, 0.059715871789770, 0.470142064105115, 0.470142064105115 },
			{ 0.132394152788506, 0.470142064105115, 0.059715871789770, 0.470142064105115 },
			{ 0.132394152788506, 0.470142064105115, 0.470142064105115, 0.059715871789770 },
			{ 0.125939180544827, 0.797426985353087, 0.101286507323456, 0.101286507323456 },
			{ 0.125939180544827, 0.101286507323456, 0.797426985353087, 0.101286507323456 },
			{ 0.125939180544827, 0.101286507323456, 0.101286507323456, 0.797426985353087 } };
		double sum = 0, total = 0;
		for (t = 0; t < S->nt; t++)
		{
			double const *p0 = q.pt[S->tv[3 * t]].y, *p1 = q.pt[S->tv[3 * t + 1]].y, *p2 = q.pt[S->tv[3 * t + 2]].y;
			double u[3], w[3], x[3], A, s = 0;
			unsigned int i;
			for (k = 0; k < 3; k++)
			{
				u[k] = p1[k] - p0[k];
				w[k] = p2[k] - p0[k];
			}
			ms_cross(u, w, x);
			A = 0.5 * sqrt(ms_dot(x, x));
			for (i = 0; i < 7; i++)
			{
				ms_point const *c0 = &q.pt[S->tv[3 * t]], *c1 = &q.pt[S->tv[3 * t + 1]], *c2 = &q.pt[S->tv[3 * t + 2]];
				qaws_mesh_location l;
				double b[3], d;
				unsigned int near = S->tv[3 * t + (g[i][1] >= g[i][2] && g[i][1] >= g[i][3] ? 0 : g[i][2] >= g[i][3] ? 1 : 2)], tn;
				for (k = 0; k < 3; k++)
				{
					x[k] = g[i][1] * p0[k] + g[i][2] * p1[k] + g[i][3] * p2[k];
					b[k] = g[i][1] * c0->b[k] + g[i][2] * c1->b[k] + g[i][3] * c2->b[k];
				}
				/* the corners on one patch: from their barycentrics interpolated there; else hinted by the nearest corner's */
				tn = c0->t;
				if (c0->t == c1->t && c1->t == c2->t)
					d = ms_polish(R, &tn, x, b);
				else
					d = R->kind == QAWS_INTERPOLANT_LINEAR ? ms_flat_project(R, &q.s, x, &tn, b) : mr_project(R, x, q.pt[near].t, &q.s, &l);
				s += g[i][0] * d;
			}
			sum += A * s;
			total += A;
		}
		*out_rms = total > 0 ? sqrt(sum / total) : 0;
	}
	/* the roots */
	for (t = 0; t < S->nt && !q.failed; t++)
	{
		ms_tri T;
		memset(T.w, 0, sizeof T.w);
		for (k = 0; k < 3; k++)
		{
			T.w[k][k] = 1;
			T.p[k] = S->tv[3 * t + k];
		}
		T.root = t;
		T.U = ms_bound(&q, &T);
		if (T.U > q.lower)
			ms_push(&q, &T);
	}
	/* branch and bound: the largest upper bound split into four while it is more than tol above the lower bound */
	while (q.nh && !q.failed && q.heap[0].U - q.lower > tol)
	{
		ms_tri T, C;
		unsigned int mp[3], e;
		double mw[3][3];
		static unsigned char const child[4][3] = { { 0, 3, 5 }, { 3, 1, 4 }, { 5, 4, 2 }, { 3, 4, 5 } };
		ms_pop(&q, &T);
		/* the budget, past the splits that bring the root down to R's triangles' size */
		if (q.splits[T.root] >= budget + (unsigned int)fmin(1e6, 2 * ms_area(S->patch[T.root].p) / area))
		{
			aside = T.U > aside ? T.U : aside;
			continue;
		}
		q.splits[T.root]++;
		for (e = 0; e < 3; e++)
		{
			double y[3];
			for (k = 0; k < 3; k++)
				mw[e][k] = 0.5 * (T.w[e][k] + T.w[(e + 1) % 3][k]);
			ms_eval((double const(*)[3])S->patch[T.root].p, mw[e], y, NULL, NULL);
			mp[e] = ms_point_add(&q, y, 1, q.pt[T.p[e]].t, q.pt[T.p[(e + 1) % 3]].t);
		}
		if (q.failed)
			break;
		for (e = 0; e < 4; e++)
		{
			for (k = 0; k < 3; k++)
			{
				unsigned int i = child[e][k];
				C.p[k] = i < 3 ? T.p[i] : mp[i - 3];
				memcpy(C.w[k], i < 3 ? T.w[i] : mw[i - 3], sizeof C.w[k]);
			}
			C.root = T.root;
			C.U = ms_bound(&q, &C);
			if (C.U > q.lower && !ms_push(&q, &C))
				break;
		}
	}
	upper = q.lower;
	if (q.nh && q.heap[0].U > upper)
		upper = q.heap[0].U;
	upper = aside > upper ? aside : upper;
	/* rounding: the bounds widened by a few units in the last place of the coordinates */
	for (t = 0; t < S->nt + R->nt; t++)
		for (v = 0; v < 3; v++)
			for (k = 0; k < 3; k++)
			{
				double x = fabs(t < S->nt ? S->patch[t].p[v][k] : R->patch[t - S->nt].p[v][k]);
				scale = x > scale ? x : scale;
			}
	scale = 1e-14 * (scale + S->extent + R->extent);
	bounds[0] = q.lower - scale > 0 ? q.lower - scale : 0;
	bounds[1] = upper + scale;
	memcpy(worst, q.worst, sizeof q.worst);
	*converged = bounds[1] - bounds[0] <= tol + 2 * scale;
	mr_scratch_free(&q.s);
	free(q.pt);
	free(q.heap);
	free(q.splits);
	return q.failed ? QAWS_STATUS_ALLOCATION_FAILURE : QAWS_STATUS_OK;
}

/* stored outward: a lower bound rounded down, an upper bound up */
static qaws_scalar ms_store(double x, int up)
{
	qaws_scalar f = (qaws_scalar)x;
#if QAWS_SCALAR_IS_FLOAT
	if (up ? (double)f < x : (double)f > x)
		f = nextafterf(f, up ? HUGE_VALF : -HUGE_VALF);
#else
	(void)up;
#endif
	return f;
}

qaws_status qaws_mesh_3d_measure_deviation(qaws_mesh_3d const* mesh, qaws_mesh_reference const* reference, qaws_deviation_desc const* desc,
	qaws_mesh_deviation* out)
{
	qaws_mesh_reference_desc ld;
	qaws_mesh_reference* lin = NULL;
	qaws_status st;
	double prec, lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, diag, to[2], from[2], rms = 0, wt[3], wf[3];
	unsigned int budget, t, k, i;
	int c1 = 0, c2 = 0;
	if (!mesh || !reference || !out || !reference->nt || !mesh->hm.nf || (desc && !(desc->precision >= 0)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	prec = desc && desc->precision > 0 ? (double)desc->precision : 1e-4;
	budget = desc && desc->budget ? desc->budget : 256;
	for (t = 0; t < reference->nt; t++)
		for (i = 0; i < 3; i++)
			for (k = 0; k < 3; k++)
			{
				double x = reference->patch[t].p[i][k];
				lo[k] = x < lo[k] ? x : lo[k];
				hi[k] = x > hi[k] ? x : hi[k];
			}
	diag = sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
	/* the mesh as its flat triangles: what is measured one way, and measured against the other */
	memset(&ld, 0, sizeof ld);
	ld.interpolant = QAWS_INTERPOLANT_LINEAR;
	st = qaws_mesh_reference_create(mesh, &ld, &lin);
	if (st == QAWS_STATUS_OK)
		st = ms_run(lin, reference, prec * diag, budget, &rms, to, wt, &c1);
	if (st == QAWS_STATUS_OK)
		st = ms_run(reference, lin, prec * diag, budget, NULL, from, wf, &c2);
	qaws_mesh_reference_destroy(lin);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(out, 0, sizeof *out);
	out->to_reference[0] = ms_store(to[0], 0);
	out->to_reference[1] = ms_store(to[1], 1);
	out->from_reference[0] = ms_store(from[0], 0);
	out->from_reference[1] = ms_store(from[1], 1);
	out->hausdorff[0] = out->to_reference[0] > out->from_reference[0] ? out->to_reference[0] : out->from_reference[0];
	out->hausdorff[1] = out->to_reference[1] > out->from_reference[1] ? out->to_reference[1] : out->from_reference[1];
	out->rms = (qaws_scalar)rms;
	out->diagonal = (qaws_scalar)diag;
	for (k = 0; k < 3; k++)
	{
		out->to_reference_at[k] = (qaws_scalar)wt[k];
		out->from_reference_at[k] = (qaws_scalar)wf[k];
	}
	out->converged = c1 && c2;
	return QAWS_STATUS_OK;
}
