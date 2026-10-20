/*
 * Remeshing by a centroidal Voronoi tessellation of the reference (Du,
 * Faber and Gunzburger 1999), on its restricted Voronoi diagram (Yan, Levy
 * et al. 2009), Lloyd then L-BFGS (Liu et al. 2009, geogram's 5 and 30).
 *
 * The seeds are the working mesh's vertices after LOCAL: their number set
 * from the size (each takes a hexagon of the target), each feature line
 * given as many vertices as its length over the size and spread along it
 * by the size (a one-dimensional CVT; lines alike, as the two sides of a
 * part thinner than the size, come out alike, their vertices facing each
 * other across it) and locked, the others moved on the energy of the
 * Voronoi diagram restricted to the reference, density size^-4. The
 * reference is sampled into flat sub-triangles once; a cell is clipped
 * from them (Sutherland-Hodgman, the bisectors of its nearest seeds up to
 * Levy and Bonneel's radius of security), found breadth first from its
 * seed, region by region, so that no cell crosses a crease. The connectivity stays the working mesh's,
 * kept extrinsic-Delaunay by flips: on a well-sampled surface this is the
 * restricted Delaunay triangulation, and it stays a manifold whatever the
 * seeds do.
 *
 * Every move keeps the faces round the seed from folding: none turned by
 * more than 60 degrees in one step, none turned to face against the
 * reference's normals at its corners (steps that each turn a little fold a
 * part curved tighter than the size over a few of them), and a seed walks
 * along the reference rather than jumping to its other side.
 *
 * Lloyd and L-BFGS anneal a lattice's dislocations (5-7 pairs) slowly: 30
 * steps leave a fifth of the vertices irregular. So the seeds are solved
 * twice: a coarse level first (a quarter of them at twice the size, the
 * feature lines at twice the spacing), each triangle then split in four,
 * which keeps every vertex's valence and doubles the lines (kept within 5%
 * of the size), and the full level relaxed from there; the irregular
 * share drops about fourfold for about twice the time. Where the split
 * folds faces (a part curved tighter than twice the size) the coarse level
 * is dropped. Below 200 seeds, a single level.
 */

#include "qaws_remesh_cvt.h"
#include "../internal/qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RV_KMAX 16u                 /* sub-triangles along a reference edge at most */
#define RV_POLY 96u                 /* a clipped polygon's vertices at most */
#define RV_HISTORY 7u               /* L-BFGS pairs */
#define RV_PI 3.14159265358979323846

rv_probe* rv_trace = NULL;

typedef struct rv_pt
{
	double x[3], rho;
} rv_pt;

typedef struct rv_cand
{
	double d2;
	unsigned int j;
} rv_cand;

typedef struct rv_cvt
{
	qaws_hm* m;
	qaws_mesh_reference const* r;
	qaws_size_field const* size;
	double L;
	unsigned int sa, ra;
	mr_scratch* s;
	qaws_mesh_report* rep;
	/* the integration mesh: per reference triangle its grid, its first sub-vertex and sub-triangle */
	unsigned int *k, *v0, *t0;
	double* sv;                     /* per sub-vertex: its point and density (4 each) */
	unsigned int* tri;              /* per sub-triangle: its sub-vertices */
	unsigned int *a0, *adj;         /* its neighbours (CSR), across the reference's edges but its creases */
	unsigned int nst;
	/* the seeds (the vertices; their count fixed once set) */
	unsigned int n;
	unsigned char* locked;
	unsigned int *e0, *ereg, *estart;   /* per seed its entries (CSR): a region, the sub-triangle a locked seed's cell starts from */
	/* neighbour search */
	qaws_bp_box* box;
	qaws_bp_grid* grid;
	rv_cand* cand;
	unsigned int ncand;
	double reach2;                  /* every seed not among the candidates is at least this far (squared) */
	double hcell;                   /* the grid's smallest cell side */
	/* breadth first search */
	unsigned int *stamp, *queue, gen;
	/* flips: edges as vertex pairs */
	unsigned int *eq, ceq;
	unsigned int evaluations;
	double scale;                   /* the size attribute over the size: 2 on the coarse level */
	double mass, covered;            /* the integral of the density over the reference, over the cells at the last evaluation */
	double slack;                   /* a line keeps its count while within this share of its length over the size (or half a segment) */
} rv_cvt;

static double rv_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void rv_cross(double const* a, double const* b, double* out)
{
	out[0] = a[1] * b[2] - a[2] * b[1];
	out[1] = a[2] * b[0] - a[0] * b[2];
	out[2] = a[0] * b[1] - a[1] * b[0];
}

static double rv_d2(double const* a, double const* b)
{
	return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
}

static double* rv_x(qaws_hm const* m, unsigned int v)
{
	return m->x + 3 * (size_t)v;
}

/* the size at a location */
static double rv_size(rv_cvt const* w, qaws_mesh_location const* l, double const* x)
{
	double h = w->size ? sf_eval(w->size, l, x) : w->L;
	return h > 1e-12 * w->L ? h : 1e-12 * w->L;
}

static double rv_h(rv_cvt const* w, unsigned int v)
{
	return hm_attr_at(w->m, w->sa, v)[0];
}

static unsigned int rv_region(rv_cvt const* w, unsigned int f)
{
	return (unsigned int)hm_attr_at(w->m, w->ra, f)[0];
}

/* the reference's normal at v's footpoint */
static void rv_normal(rv_cvt const* w, unsigned int v, double* n)
{
	double x[3];
	mr_eval(w->r, &w->m->foot[v], x, n);
}

/* the reference's normal at v's footpoint added to s, unless v is on a feature (where it depends on the side) */
static void rv_side(rv_cvt const* w, unsigned int v, double* s)
{
	double n[3];
	unsigned int k;
	if (hm_feature_degree(w->m, v))
		return;
	rv_normal(w, v, n);
	for (k = 0; k < 3; k++)
		s[k] += n[k];
}

/* the faces facing against the reference's normals at their corners off the features */
static unsigned int rv_folded(rv_cvt const* w)
{
	qaws_hm const* m = w->m;
	unsigned int f, q, bad = 0;
	for (f = 0; f < m->nf; f++)
	{
		double n[3], s[3] = { 0, 0, 0 };
		if (m->fflag[f] & HM_DELETED)
			continue;
		hm_face_normal(m, f, n);
		for (q = 0; q < 3; q++)
			rv_side(w, m->cv[m->fc[f] + q], s);
		bad += rv_dot(n, s) < 0;
	}
	return bad;
}

/* would moving v to y (the reference's normal there ny) turn a face round v (faces holding `skip` aside) by more than acos(c),
   flatten it, or turn it to face against the reference's normals at its corners off the features? */
static int rv_folds(rv_cvt const* w, unsigned int v, double const* y, double const* ny, unsigned int skip, double c)
{
	qaws_hm const* m = w->m;
	unsigned int e, k;
	double nv[3] = { 0, 0, 0 }, nw[3] = { 0, 0, 0 };
	if (!hm_feature_degree(m, v))
	{
		rv_normal(w, v, nv);
		memcpy(nw, ny, sizeof nw);
	}
	for (e = m->vc[v]; e != HM_NONE; e = hm_round(m, v, e))
	{
		unsigned int a = hm_to(m, e), b = m->cv[hm_prev(m, e)];
		double const *p = rv_x(m, v), *A = rv_x(m, a), *B = rv_x(m, b);
		double u[3], t[3], n0[3], n1[3], s0[3], s1[3] = { 0, 0, 0 }, d0, d1;
		if (a == skip || b == skip)
			continue;
		for (k = 0; k < 3; k++) { u[k] = A[k] - p[k]; t[k] = B[k] - p[k]; }
		rv_cross(u, t, n0);
		for (k = 0; k < 3; k++) { u[k] = A[k] - y[k]; t[k] = B[k] - y[k]; }
		rv_cross(u, t, n1);
		d0 = rv_dot(n0, n0);
		d1 = rv_dot(n1, n1);
		if (rv_dot(n0, n1) <= c * sqrt(d0 * d1) || d1 < 1e-12 * d0)
			return 1;
		/* turned over against the surface: a fold that small steps build up on a part curved tighter than the size; judged at
		   the corners only where the face leans more than 60 degrees from the surface at y */
		if (rv_dot(n1, nw) > 0.5 * sqrt(d1))
			continue;
		rv_side(w, a, s1);
		rv_side(w, b, s1);
		for (k = 0; k < 3; k++)
		{
			s0[k] = s1[k] + nv[k];
			s1[k] += nw[k];
		}
		if (rv_dot(n1, s1) <= 0 && rv_dot(n0, s0) > 0)
			return 1;
	}
	return 0;
}

/* v's neighbours along its feature lines, at most cap; their count */
static unsigned int rv_line(qaws_hm const* m, unsigned int v, unsigned int* f, unsigned int cap)
{
	unsigned int nf = 0, c, last = HM_NONE;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
	{
		if ((m->cflag[c] & HM_FEATURE) && nf < cap)
			f[nf++] = hm_to(m, c);
		last = c;
	}
	if (nf < cap && last != HM_NONE && hm_is_border_vertex(m, v) && (m->cflag[hm_prev(m, last)] & HM_FEATURE))
		f[nf++] = m->cv[hm_prev(m, last)];
	return nf;
}

/* is v on a feature but not inside a line (a corner, or where lines meet or stop)? */
static int rv_line_end(qaws_hm const* m, unsigned int v)
{
	unsigned int fd = hm_feature_degree(m, v);
	return fd && (fd != 2 || (m->vflag[v] & HM_CORNER));
}

/* ---- the integration mesh ----------------------------------------------- */

/* within a reference triangle of grid k: the upward sub-triangle (i, j) has corners (i, j), (i + 1, j), (i, j + 1) in units
   of 1 / k of (a, b), the downward one (i + 1, j), (i + 1, j + 1), (i, j + 1) */
#define RV_UP(k, i, j) (2 * (k) * (j) - (j) * (j) + 2 * (i))
#define RV_VERT(k, i, j) ((j) * (2 * (k) + 3 - (j)) / 2 + (i))

/* the sub-triangle on segment q of edge e (from corner e to e + 1, q of k) */
static unsigned int rv_edge_tri(unsigned int k, unsigned int e, unsigned int q)
{
	return e == 0 ? RV_UP(k, q, 0) : e == 1 ? RV_UP(k, k - 1 - q, q) : RV_UP(k, 0, k - 1 - q);
}

/* the sub-triangles across segment q of t's edge e (none across a crease or the border); their count */
static unsigned int rv_across(rv_cvt const* w, unsigned int t, unsigned int e, unsigned int q, unsigned int* out)
{
	unsigned int nb = w->r->ta[3 * t + e], u, ue, ku, k = w->k[t], lo, hi, p, n = 0;
	if (nb == HM_NONE || w->r->tcrease[3 * t + e])
		return 0;
	u = nb / 3;
	ue = nb % 3;
	ku = w->k[u];
	/* u's edge runs the other way: the segment [q, q + 1] / k is [k - q - 1, k - q] / k there */
	lo = (k - q - 1) * ku / k;
	hi = ((k - q) * ku + k - 1) / k;
	for (p = lo; p < hi && p < ku; p++)
		out[n++] = w->t0[u] + rv_edge_tri(ku, ue, p);
	return n;
}

/* sub-triangle (i, j, down) of t: its sub-vertices and neighbours; the neighbours' count */
static unsigned int rv_sub(rv_cvt const* w, unsigned int t, unsigned int i, unsigned int j, int down, unsigned int* sv, unsigned int* nb)
{
	unsigned int k = w->k[t], v = w->v0[t], b = w->t0[t], n = 0;
	if (down)
	{
		sv[0] = v + RV_VERT(k, i + 1, j);
		sv[1] = v + RV_VERT(k, i + 1, j + 1);
		sv[2] = v + RV_VERT(k, i, j + 1);
		nb[n++] = b + RV_UP(k, i + 1, j);
		nb[n++] = b + RV_UP(k, i, j + 1);
		nb[n++] = b + RV_UP(k, i, j);
		return n;
	}
	sv[0] = v + RV_VERT(k, i, j);
	sv[1] = v + RV_VERT(k, i + 1, j);
	sv[2] = v + RV_VERT(k, i, j + 1);
	if (j > 0)
		nb[n++] = b + RV_UP(k, i, j - 1) + 1;
	else
		n += rv_across(w, t, 0, i, nb + n);
	if (i + j + 2 <= k)
		nb[n++] = b + RV_UP(k, i, j) + 1;
	else
		n += rv_across(w, t, 1, j, nb + n);
	if (i > 0)
		nb[n++] = b + RV_UP(k, i - 1, j) + 1;
	else
		n += rv_across(w, t, 2, k - 1 - j, nb + n);
	return n;
}

/* the sub-triangle holding location l */
static unsigned int rv_locate(rv_cvt const* w, qaws_mesh_location const* l)
{
	unsigned int t = l->patch, k = w->k[t], i, j;
	double a = (double)l->a * k, b = (double)l->b * k;
	i = a <= 0 ? 0 : a >= k - 1 ? k - 1 : (unsigned int)a;
	j = b <= 0 ? 0 : b >= k - 1 ? k - 1 : (unsigned int)b;
	if (i + j > k - 1)
		return w->t0[t] + RV_UP(k, k - 1 - j, j);
	if (a - i + b - j > 1 && i + j + 2 <= k)
		return w->t0[t] + RV_UP(k, i, j) + 1;
	return w->t0[t] + RV_UP(k, i, j);
}

/* each reference triangle split k^2 ways, k from its control net's longest edge over its smallest size; the seeds wanted
   (the integral of 2 / (sqrt(3) h^2): a hexagon of side h per seed) */
static int rv_build(rv_cvt* w, double* nstar)
{
	static unsigned char const net[9][2] = { { 0, 3 }, { 3, 4 }, { 4, 1 }, { 1, 5 }, { 5, 6 }, { 6, 2 }, { 2, 7 }, { 7, 8 }, { 8, 0 } };
	qaws_mesh_reference const* r = w->r;
	unsigned int t, nsv = 0, nst = 0, nadj = 0, pass, i, j, q;
	double total = 0;
	w->k = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (r->nt + 1)));
	w->v0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (r->nt + 1)));
	w->t0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (r->nt + 1)));
	if (!w->k || !w->v0 || !w->t0)
		return 0;
	for (t = 0; t < r->nt; t++)
	{
		static double const ca[4] = { 0, 1, 0, 1.0 / 3 }, cb[4] = { 0, 0, 1, 1.0 / 3 };
		mi_tri const* P = &r->patch[t];
		double hmin = HUGE_VAL, len = 0, kk;
		for (q = 0; q < 4; q++)
		{
			qaws_mesh_location l;
			double x[3], h;
			l.patch = t;
			l.a = (qaws_scalar)ca[q];
			l.b = (qaws_scalar)cb[q];
			mi_eval(P, ca[q], cb[q], x, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
			h = rv_size(w, &l, x);
			hmin = h < hmin ? h : hmin;
		}
		for (q = 0; q < 9; q++)
			len = fmax(len, sqrt(rv_d2(P->p[net[q][0]], P->p[net[q][1]])));
		kk = ceil(3 * len / (w->scale * hmin));
		{
			/* a sliver: no more sub-triangles than 8 per size squared of its area, which its long side alone would ask */
			double e1[3], e2[3], n[3];
			for (q = 0; q < 3; q++)
			{
				e1[q] = P->p[1][q] - P->p[0][q];
				e2[q] = P->p[2][q] - P->p[0][q];
			}
			rv_cross(e1, e2, n);
			kk = fmin(kk, ceil(sqrt(4 * sqrt(rv_dot(n, n))) / (w->scale * hmin)));
		}
		w->k[t] = kk < 1 ? 1 : kk > RV_KMAX ? RV_KMAX : (unsigned int)kk;
		w->v0[t] = nsv;
		w->t0[t] = nst;
		nsv += (w->k[t] + 1) * (w->k[t] + 2) / 2;
		nst += w->k[t] * w->k[t];
	}
	w->nst = nst;
	w->sv = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 4, ((size_t)nsv + 1)));
	w->tri = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, ((size_t)nst + 1)));
	w->a0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nst + 1)));
	w->stamp = (unsigned int*)qaws_internal_zalloc((size_t)nst + 1, sizeof(unsigned int));
	w->queue = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nst + 1)));
	if (!w->sv || !w->tri || !w->a0 || !w->stamp || !w->queue)
		return 0;
	for (t = 0; t < r->nt; t++)
	{
		unsigned int k = w->k[t];
		for (j = 0; j <= k; j++)
			for (i = 0; i + j <= k; i++)
			{
				double* o = w->sv + 4 * ((size_t)w->v0[t] + RV_VERT(k, i, j));
				qaws_mesh_location l;
				double h;
				l.patch = t;
				l.a = (qaws_scalar)((double)i / k);
				l.b = (qaws_scalar)((double)j / k);
				mi_eval(&r->patch[t], (double)l.a, (double)l.b, o, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
				h = rv_size(w, &l, o);
				o[3] = pow(w->L / h, 4);
			}
	}
	/* the sub-triangles and their neighbours: counted, then stored */
	for (pass = 0; pass < 2; pass++)
	{
		nadj = 0;
		for (t = 0; t < r->nt; t++)
		{
			unsigned int k = w->k[t];
			for (j = 0; j < k; j++)
				for (i = 0; i + j < k; i++)
				{
					int down;
					for (down = 0; down < 2 && (!down || i + j + 2 <= k); down++)
					{
						unsigned int u = w->t0[t] + RV_UP(k, i, j) + (unsigned int)down, sv[3], nb[3 * (RV_KMAX + 2)], nn, c;
						nn = rv_sub(w, t, i, j, down, sv, nb);
						if (pass == 0)
						{
							memcpy(w->tri + 3 * (size_t)u, sv, sizeof sv);
							w->a0[u] = nadj;
						}
						else
							for (c = 0; c < nn; c++)
								w->adj[nadj + c] = nb[c];
						nadj += nn;
					}
				}
		}
		if (pass == 0)
		{
			w->a0[nst] = nadj;
			w->adj = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nadj + 1)));
			if (!w->adj)
				return 0;
		}
	}
	w->mass = 0;
	for (t = 0; t < nst; t++)
	{
		double const *A = w->sv + 4 * (size_t)w->tri[3 * t], *B = w->sv + 4 * (size_t)w->tri[3 * t + 1], *C = w->sv + 4 * (size_t)w->tri[3 * t + 2];
		double u[3], v[3], n[3], a;
		unsigned int c;
		for (c = 0; c < 3; c++) { u[c] = B[c] - A[c]; v[c] = C[c] - A[c]; }
		rv_cross(u, v, n);
		a = 0.5 * sqrt(rv_dot(n, n));
		total += a * (sqrt(A[3]) + sqrt(B[3]) + sqrt(C[3])) / 3;
		w->mass += a * (A[3] + B[3] + C[3]) / 3;
	}
	*nstar = total * 2 / (sqrt(3.0) * w->L * w->L);
	return 1;
}

/* ---- the seed count ----------------------------------------------------- */

static int rv_cmp_key(void const* a, void const* b)
{
	rv_cand const *p = (rv_cand const*)a, *q = (rv_cand const*)b;
	return p->d2 < q->d2 ? -1 : p->d2 > q->d2 ? 1 : p->j < q->j ? -1 : p->j > q->j;
}

/* would a new edge from p to q run along a crease (both ends on features, its middle on one): a chord over a feature line? */
static int rv_chord(rv_cvt* w, unsigned int p, unsigned int q)
{
	qaws_hm const* m = w->m;
	double mid[3], tol = 1e-3 * 0.5 * (rv_h(w, p) + rv_h(w, q));
	qaws_mesh_location l;
	unsigned int k;
	if (!hm_feature_degree(m, p) || !hm_feature_degree(m, q))
		return 0;
	for (k = 0; k < 3; k++)
		mid[k] = 0.5 * (rv_x(m, p)[k] + rv_x(m, q)[k]);
	return mr_project_crease(w->r, mid, m->foot[p].patch, w->s, &l) < tol * tol;
}

/* may x, an end of c's edge, go into the other end y (which stays)? Off the features, or from a feature line along it; never a
   chord along a crease, never a face folded or turned by more than 60 degrees */
static int rv_can_remove(rv_cvt* w, unsigned int c, unsigned int x, unsigned int y)
{
	qaws_hm const* m = w->m;
	unsigned int ring[256], n, i, fd = hm_feature_degree(m, x), f[2] = { HM_NONE, HM_NONE };
	double ny[3];
	if (fd && (fd != 2 || (m->vflag[x] & HM_CORNER) || !(m->cflag[c] & HM_FEATURE) || rv_line(m, x, f, 2) != 2))
		return 0;
	rv_normal(w, y, ny);
	if (!hm_can_collapse(m, c, x) || rv_folds(w, x, rv_x(m, y), ny, y, 0.5))
		return 0;
	n = hm_ring(m, x, ring, 256);
	for (i = 0; i < n && i < 256; i++)
		if (ring[i] != y && ring[i] != f[0] && ring[i] != f[1] && hm_find_edge(m, y, ring[i]) == HM_NONE && rv_chord(w, y, ring[i]))
			return 0;   /* (x's other neighbour on its line is y's next one there) */
	return n <= 256;
}

/* as many vertices as the size asks, within 1%: the shortest edges collapsed or the longest split, a few rounds; the edges
   off the features with both ends free, or with `lines` (a coarse level) those with an end on a feature too, the free end going */
static int rv_count(rv_cvt* w, double nstar, int lines)
{
	qaws_hm* m = w->m;
	unsigned int round;
	for (round = 0; round < (lines ? 16u : 5u); round++)
	{
		unsigned int c, ne = 0, i, done = 0, want, *e;
		int more = m->nv < nstar;
		rv_cand* key;
		unsigned char* touched;
		if (fabs(m->nv - nstar) <= 0.01 * nstar)
			break;
		want = (unsigned int)fabs(m->nv - nstar);
		key = (rv_cand*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(rv_cand), (m->nc + 1)));
		e = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (m->nc + 1)));
		touched = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, 1);
		if (!key || !e || !touched)
		{
			qaws_internal_free(key); qaws_internal_free(e); qaws_internal_free(touched);
			return 0;
		}
		for (c = 0; c < m->nc; c++)
			if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
			{
				unsigned int a = m->cv[c], b = hm_to(m, c);
				double l = sqrt(rv_d2(rv_x(m, a), rv_x(m, b))) / (0.5 * (rv_h(w, a) + rv_h(w, b)));
				if ((m->cflag[c] & HM_FEATURE) || (!lines && (hm_feature_degree(m, a) || hm_feature_degree(m, b))))
					continue;
				e[2 * ne] = a;
				e[2 * ne + 1] = b;
				key[ne].d2 = more ? -l : l;
				key[ne].j = ne;
				ne++;
			}
		qsort(key, ne, sizeof *key, rv_cmp_key);
		for (i = 0; i < ne && done < want; i++)
		{
			unsigned int a = e[2 * key[i].j], b = e[2 * key[i].j + 1], cc, x, y, side;
			if (touched[a] || touched[b] || (cc = hm_find_edge(m, a, b)) == HM_NONE)
				continue;
			if (more)
			{
				double mid[3], p[3];
				qaws_mesh_location l;
				unsigned int reg = rv_region(w, m->cf[cc]), hint = m->foot[a].patch, v, k;
				for (k = 0; k < 3; k++)
					mid[k] = 0.5 * (rv_x(m, a)[k] + rv_x(m, b)[k]);
				if (mr_project_region(w->r, mid, hint < w->r->nt && w->r->treg[hint] == reg ? hint : HM_NONE, reg, w->s, &l) == HUGE_VAL
					|| l.patch == HM_NONE)
					continue;
				mr_eval(w->r, &l, p, NULL);
				v = hm_split_edge(m, cc, 0.5, p, 1);
				if (v == HM_NONE)
				{
					qaws_internal_free(key); qaws_internal_free(e); qaws_internal_free(touched);
					return 0;
				}
				m->foot[v] = l;
				hm_attr_at(m, w->sa, v)[0] = w->scale * rv_size(w, &l, p);
				touched[a] = touched[b] = 1;
				w->rep->splits++;
				done++;
				continue;
			}
			/* a collapse: one end into the other, which stays; the rings round both left alone this round */
			x = y = HM_NONE;
			for (side = 0; side < 2; side++)
			{
				x = side ? a : b;
				y = side ? b : a;
				if (!hm_feature_degree(m, x) && !(m->vflag[x] & HM_CORNER) && rv_can_remove(w, cc, x, y))
					break;
			}
			if (side == 2)
				continue;
			{
				unsigned int ring[256], nr = hm_ring(m, y, ring, 128), q;
				nr = nr > 128 ? 128 : nr;
				nr += hm_ring(m, x, ring + nr, 128);
				if (!hm_collapse(m, cc, x, 1, rv_x(m, y)))
				{
					qaws_internal_free(key); qaws_internal_free(e); qaws_internal_free(touched);
					return 0;
				}
				touched[x] = touched[y] = 1;
				for (q = 0; q < nr && q < 256; q++)
					touched[ring[q]] = 1;
			}
			w->rep->collapses++;
			done++;
		}
		qaws_internal_free(key); qaws_internal_free(e); qaws_internal_free(touched);
		hm_compact(m);
		if (!done)
			break;
	}
	return 1;
}

/* every edge split at its middle on the reference (on the creases for a feature edge), then the new vertices' edges to the
   old ones across the faces flipped: each triangle in four, every old vertex keeping its valence; 0 on allocation failure.
   Left whole: an edge between two features that is not one (across a thin part), a feature edge no longer than 4/3 of the
   size (a coarse line kept short: one too short to split, or where a collapse was refused) */
static int rv_subdivide(rv_cvt* w)
{
	qaws_hm* m = w->m;
	unsigned int c, ne = 0, i, nv0 = m->nv, *e = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (m->nc + 1))), *src;
	if (!e)
		return 0;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			e[2 * ne] = m->cv[c];
			e[2 * ne + 1] = hm_to(m, c);
			ne++;
		}
	src = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, ((size_t)ne + 1)));
	if (!src)
	{
		qaws_internal_free(e);
		return 0;
	}
	for (i = 0; i < ne; i++)
	{
		unsigned int a = e[2 * i], b = e[2 * i + 1], cc = hm_find_edge(m, a, b), v, k;
		double mid[3], p[3], d;
		qaws_mesh_location l;
		if (cc == HM_NONE || (!(m->cflag[cc] & HM_FEATURE) && hm_feature_degree(m, a) && hm_feature_degree(m, b)))
			continue;   /* (across a part narrower than twice the size: no vertex inside it) */
		if ((m->cflag[cc] & HM_FEATURE) && rv_d2(rv_x(m, a), rv_x(m, b)) <= (4.0 / 9) * (rv_h(w, a) + rv_h(w, b)) * (rv_h(w, a) + rv_h(w, b)))
			continue;   /* (a feature edge no longer than 4/3 of its size: one the coarse line left short) */
		for (k = 0; k < 3; k++)
			mid[k] = 0.5 * (rv_x(m, a)[k] + rv_x(m, b)[k]);
		if (m->cflag[cc] & HM_FEATURE)
			d = mr_project_crease(w->r, mid, m->foot[a].patch, w->s, &l);
		else
		{
			unsigned int reg = rv_region(w, m->cf[cc]), hint = m->foot[a].patch;
			d = mr_project_region(w->r, mid, hint < w->r->nt && w->r->treg[hint] == reg ? hint : HM_NONE, reg, w->s, &l);
		}
		if (d == HUGE_VAL || l.patch == HM_NONE)
			continue;
		mr_eval(w->r, &l, p, NULL);
		v = hm_split_edge(m, cc, 0.5, p, 1);
		if (v == HM_NONE)
		{
			qaws_internal_free(e); qaws_internal_free(src);
			return 0;
		}
		m->foot[v] = l;
		hm_attr_at(m, w->sa, v)[0] = w->scale * rv_size(w, &l, p);
		src[2 * (v - nv0)] = a;
		src[2 * (v - nv0) + 1] = b;
	}
	/* the edges from a new vertex to an old one that is not an end of its own edge, flipped where that joins two new ones (not
	   where an edge was left whole) */
	ne = 0;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && m->co[c] != HM_NONE && c < m->co[c])
		{
			unsigned int p = m->cv[c], q = hm_to(m, c), nw, old;
			if ((p >= nv0) == (q >= nv0))
				continue;
			nw = p >= nv0 ? p : q;
			old = p >= nv0 ? q : p;
			if (old == src[2 * (nw - nv0)] || old == src[2 * (nw - nv0) + 1])
				continue;
			e[2 * ne] = p;
			e[2 * ne + 1] = q;
			ne++;
		}
	for (i = 0; i < ne; i++)
	{
		unsigned int cc = hm_find_edge(m, e[2 * i], e[2 * i + 1]);
		if (cc != HM_NONE && !(m->cflag[cc] & HM_FEATURE) && m->co[cc] != HM_NONE && m->cv[hm_prev(m, cc)] >= nv0
			&& m->cv[hm_prev(m, m->co[cc])] >= nv0 && hm_can_flip(m, cc) && !hm_flip(m, cc))
		{
			qaws_internal_free(e); qaws_internal_free(src);
			return 0;
		}
	}
	qaws_internal_free(e);
	qaws_internal_free(src);
	hm_compact(m);
	return 1;
}

/* ---- creases ------------------------------------------------------------ */

/* the vertex of chain ch (n segments) at metric parameter t of the polyline (cumulative arc over the size in s) */
static void rv_at(rv_cvt const* w, unsigned int const* ch, double const* s, unsigned int n, double t, double* y)
{
	unsigned int i = 0, k;
	double u;
	while (i + 1 < n && s[i + 1] < t)
		i++;
	u = s[i + 1] > s[i] ? (t - s[i]) / (s[i + 1] - s[i]) : 0;
	u = u < 0 ? 0 : u > 1 ? 1 : u;
	for (k = 0; k < 3; k++)
		y[k] = rv_x(w->m, ch[i])[k] + u * (rv_x(w->m, ch[i + 1])[k] - rv_x(w->m, ch[i])[k]);
}

/* chain ch's cumulative arc over the size, s[0..n]; its total */
static double rv_arc(rv_cvt const* w, unsigned int const* ch, unsigned int n, double* s)
{
	unsigned int i;
	s[0] = 0;
	for (i = 0; i < n; i++)
		s[i + 1] = s[i] + sqrt(rv_d2(rv_x(w->m, ch[i]), rv_x(w->m, ch[i + 1]))) / (0.5 * (rv_h(w, ch[i]) + rv_h(w, ch[i + 1])));
	return s[n];
}

/* one chain of a feature line, ch[0..n] (its ends fixed; a loop's end its start): as many segments as its length over the
   size (rounded; a loop at least 3), the longest split or the vertex between the shortest pair merged along it, then each
   vertex moved, a few sweeps, to its share of the length over the size (a 1D CVT of density size^-3, which spaces the
   vertices by the size), back onto the creases, unless a face round it would fold. Lines alike (the two sides of a part
   thinner than the size, between the same corners) come out alike, so their vertices face each other across it. ch holds
   room for cap vertices; the chain's segments at the end */
static unsigned int rv_chain(rv_cvt* w, unsigned int* ch, unsigned int n, unsigned int cap, double* s, int loop)
{
	qaws_hm* m = w->m;
	unsigned int want, i, k, sweep;
	double total = rv_arc(w, ch, n, s);
	want = fabs(n - total) <= fmax(0.5, w->slack * total) ? n : (unsigned int)floor(total + 0.5);
	want = want < (loop ? 3u : 1u) ? (loop ? 3u : 1u) : want;
	while (n < want && n + 1 < cap)
	{
		unsigned int best = 0, c, v;
		double mid[3], p[3];
		qaws_mesh_location l;
		for (i = 1; i < n; i++)
			if (s[i + 1] - s[i] > s[best + 1] - s[best])
				best = i;
		c = hm_find_edge(m, ch[best], ch[best + 1]);
		for (k = 0; k < 3; k++)
			mid[k] = 0.5 * (rv_x(m, ch[best])[k] + rv_x(m, ch[best + 1])[k]);
		if (c == HM_NONE || mr_project_crease(w->r, mid, m->foot[ch[best]].patch, w->s, &l) == HUGE_VAL || l.patch == HM_NONE)
			break;
		mr_eval(w->r, &l, p, NULL);
		if ((v = hm_split_edge(m, c, 0.5, p, 1)) == HM_NONE)
			return 0;
		m->foot[v] = l;
		hm_attr_at(m, w->sa, v)[0] = w->scale * rv_size(w, &l, p);
		w->rep->splits++;
		memmove(ch + best + 2, ch + best + 1, sizeof(unsigned int) * (n - best));
		ch[best + 1] = v;
		n++;
		rv_arc(w, ch, n, s);
	}
	while (n > want)
	{
		/* the inner vertex whose two segments are the shortest that may go into a neighbour along the line */
		unsigned int best = 0, side = 0, tries, last = 0;
		double used = -1;
		for (tries = 0; tries < n && !best; tries++)
		{
			double low = HUGE_VAL, key;
			unsigned int q;
			best = 0;
			for (i = 1; i < n; i++)
				if ((key = s[i + 1] - s[i - 1]) < low && (key > used || (key == used && i > last)))
				{
					low = key;
					best = i;
				}
			if (!best)
				break;
			used = low;
			last = best;
			for (q = 0; q < 2; q++)
			{
				/* (into the nearer one first) */
				unsigned int y = (q == 0) == (s[best + 1] - s[best] < s[best] - s[best - 1]) ? ch[best + 1] : ch[best - 1];
				unsigned int c = hm_find_edge(m, ch[best], y);
				if (c != HM_NONE && rv_can_remove(w, c, ch[best], y))
				{
					if (!hm_collapse(m, c, ch[best], 1, rv_x(m, y)))
						return 0;
					w->rep->collapses++;
					side = 1;
					break;
				}
			}
			if (!side)
				best = 0;
		}
		if (!side)
			break;
		memmove(ch + best, ch + best + 1, sizeof(unsigned int) * (n - best));
		n--;
		rv_arc(w, ch, n, s);
	}
	for (sweep = 0; sweep < 4 && n > 1; sweep++)
	{
		double* t = s + cap;
		total = rv_arc(w, ch, n, s);
		for (i = 1; i < n; i++)
			rv_at(w, ch, s, n, total * i / n, t + 3 * (size_t)i);   /* (the targets from the line as it was) */
		for (i = 1; i < n; i++)
		{
			unsigned int v = ch[i], half;
			double y[3], x[3], nr[3], p[3];
			qaws_mesh_location l;
			memcpy(y, t + 3 * (size_t)i, sizeof y);
			memcpy(p, rv_x(m, v), sizeof p);
			for (half = 0; half < 3; half++, y[0] = 0.5 * (y[0] + p[0]), y[1] = 0.5 * (y[1] + p[1]), y[2] = 0.5 * (y[2] + p[2]))
			{
				if (mr_project_crease(w->r, y, m->foot[v].patch, w->s, &l) == HUGE_VAL || l.patch == HM_NONE)
					continue;
				mr_eval(w->r, &l, x, nr);
				if (rv_folds(w, v, x, nr, HM_NONE, 0.2))
					continue;
				memcpy(rv_x(m, v), x, sizeof x);
				m->foot[v] = l;
				hm_attr_at(m, w->sa, v)[0] = w->scale * rv_size(w, &l, x);
				w->rep->moves++;
				break;
			}
		}
	}
	return n;
}

/* every feature line's vertices from the size (rv_chain), chain by chain: from the lines' ends and corners, then the closed
   loops; corners stay; 0 on allocation failure */
static int rv_lines(rv_cvt* w)
{
	qaws_hm* m = w->m;
	unsigned int v, pass, nv = m->nv, cap = 0, *ch = NULL;
	unsigned char* seen = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nv + 1, 1);
	double* s = NULL;
	int ok = seen != NULL;
	for (pass = 0; ok && pass < 2; pass++)
		for (v = 0; ok && v < nv; v++)
		{
			unsigned int f[16], nf, q;
			if ((m->vflag[v] & HM_DELETED) || m->vc[v] == HM_NONE || seen[v] || !hm_feature_degree(m, v) || rv_line_end(m, v) != (pass == 0))
				continue;
			if (pass == 1)
				seen[v] = 1;
			nf = rv_line(m, v, f, 16);
			for (q = 0; ok && q < nf && (pass == 0 || q < 1); q++)
			{
				unsigned int prev = v, cur = f[q], n = 0, g[2];
				if (cur >= nv || seen[cur])
					continue;   /* (a vertex added on a chain done) */
				if (rv_line_end(m, cur) && cur < v)
					continue;   /* (a single segment between two ends: once) */
				for (;;)
				{
					if ((size_t)n + 2 >= cap)
					{
						unsigned int* c2;
						double* s2;
						unsigned int nc = qaws_internal_grow(cap, (size_t)n + 3, 64);
						if (!nc)
						{
							ok = 0;
							break;
						}
						cap = nc;
						c2 = (unsigned int*)qaws_internal_realloc(ch, qaws_internal_mul2(sizeof(unsigned int), cap));
						s2 = (double*)qaws_internal_realloc(s, qaws_internal_mul3(sizeof(double), 4, (size_t)cap));
						if (c2)
							ch = c2;
						if (s2)
							s = s2;
						if (!c2 || !s2)
						{
							ok = 0;
							break;
						}
					}
					if (!n)
						ch[n++] = v;
					ch[n++] = cur;
					if (cur == v || rv_line_end(m, cur) || rv_line(m, cur, g, 2) != 2)
						break;
					seen[cur] = 1;
					{
						unsigned int next = g[0] == prev ? g[1] : g[0];
						prev = cur;
						cur = next;
					}
				}
				if (!ok)
					break;
				/* room for the splits: the length over the size, twice */
				{
					double need = 2 * rv_arc(w, ch, n - 1, s) + 8;
					if (need > cap)
					{
						unsigned int* c2;
						double* s2;
						if (!(need < 4294967295.0))   /* (also NaN) */
						{
							ok = 0;
							break;
						}
						cap = (unsigned int)need;
						c2 = (unsigned int*)qaws_internal_realloc(ch, qaws_internal_mul2(sizeof(unsigned int), cap));
						s2 = (double*)qaws_internal_realloc(s, qaws_internal_mul3(sizeof(double), 4, (size_t)cap));
						if (c2)
							ch = c2;
						if (s2)
							s = s2;
						if (!c2 || !s2)
						{
							ok = 0;
							break;
						}
					}
				}
				ok = rv_chain(w, ch, n - 1, cap, s, ch[n - 1] == v) != 0;
			}
		}
	qaws_internal_free(seen);
	qaws_internal_free(ch);
	qaws_internal_free(s);
	hm_compact(m);
	return ok;
}

/* ---- flips -------------------------------------------------------------- */

/* the angle at P of triangle (P, A, B) */
static double rv_angle(double const* P, double const* A, double const* B)
{
	double u[3], t[3], n[3];
	unsigned int k;
	for (k = 0; k < 3; k++) { u[k] = A[k] - P[k]; t[k] = B[k] - P[k]; }
	rv_cross(u, t, n);
	return atan2(sqrt(rv_dot(n, n)), rv_dot(u, t));
}

/* Lawson's flips to extrinsic Delaunay over the edges off the features (at most 10 per vertex): never a chord along a crease,
   never a face turned by more than 60 degrees from the pair it replaces; 0 on allocation failure */
static int rv_delaunay(rv_cvt* w)
{
	qaws_hm* m = w->m;
	unsigned int c, head = 0, tail = 0, flips = 0, cap = 10 * m->nv, need = m->nc + 4 * cap + 8;
	if (m->nv > 0xFFFFFFFFu / 64 || m->nc > 0xFFFFFFFFu / 4)   /* need = nc + 40 nv + 8 counts in unsigned int */
		return 0;
	if (need > w->ceq)
	{
		unsigned int* g = (unsigned int*)qaws_internal_realloc(w->eq, qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)need));
		if (!g)
			return 0;
		w->eq = g;
		w->ceq = need;
	}
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && m->co[c] != HM_NONE && c < m->co[c] && !(m->cflag[c] & HM_FEATURE))
		{
			w->eq[2 * tail] = m->cv[c];
			w->eq[2 * tail + 1] = hm_to(m, c);
			tail++;
		}
	while (head < tail && flips < cap)
	{
		unsigned int a = w->eq[2 * head], b = w->eq[2 * head + 1], d, p, q, k;
		double const *A, *B, *P, *Q;
		double n0[3], n1[3], n2[3], g[3], u[3], t[3];
		head++;
		c = hm_find(m, a, b);
		if (c == HM_NONE)
			c = hm_find(m, b, a);
		if (c == HM_NONE || (d = m->co[c]) == HM_NONE || (m->cflag[c] & HM_FEATURE) || !hm_can_flip(m, c))
			continue;
		a = m->cv[c]; b = hm_to(m, c); p = m->cv[hm_prev(m, c)]; q = m->cv[hm_prev(m, d)];
		A = rv_x(m, a); B = rv_x(m, b); P = rv_x(m, p); Q = rv_x(m, q);
		if (rv_angle(P, A, B) + rv_angle(Q, B, A) <= RV_PI + 1e-12 || rv_chord(w, p, q))
			continue;
		hm_face_normal(m, m->cf[c], n0);
		hm_face_normal(m, m->cf[d], g);
		for (k = 0; k < 3; k++) { n0[k] += g[k]; u[k] = Q[k] - A[k]; t[k] = P[k] - A[k]; }
		rv_cross(u, t, n1);
		for (k = 0; k < 3; k++) { u[k] = B[k] - Q[k]; t[k] = P[k] - Q[k]; }
		rv_cross(u, t, n2);
		if (rv_dot(n1, n0) < 0.5 * sqrt(rv_dot(n1, n1) * rv_dot(n0, n0)) || rv_dot(n2, n0) < 0.5 * sqrt(rv_dot(n2, n2) * rv_dot(n0, n0)))
			continue;
		if (!hm_flip(m, c))
			return 0;
		flips++;
		/* the quad's sides may have turned non-Delaunay */
		w->eq[2 * tail] = a; w->eq[2 * tail + 1] = p;
		w->eq[2 * tail + 2] = p; w->eq[2 * tail + 3] = b;
		w->eq[2 * tail + 4] = b; w->eq[2 * tail + 5] = q;
		w->eq[2 * tail + 6] = q; w->eq[2 * tail + 7] = a;
		tail += 4;
	}
	hm_compact(m);
	w->rep->flips += flips;
	return 1;
}

/* ---- the restricted Voronoi diagram ------------------------------------- */

/* is seed j in region g? */
static int rv_in(rv_cvt const* w, unsigned int j, unsigned int g)
{
	unsigned int e;
	for (e = w->e0[j]; e < w->e0[j + 1]; e++)
		if (w->ereg[e] == g)
			return 1;
	return 0;
}

typedef struct rv_near
{
	rv_cvt* w;
	double const* z;
	double const* p;
	unsigned int self, region;
} rv_near;

static void rv_near_visit(void* user, unsigned int j)
{
	rv_near* q = (rv_near*)user;
	rv_cvt* w = q->w;
	if (j == q->self || !rv_in(w, j, q->region))
		return;
	w->cand[w->ncand].d2 = rv_d2(q->z + 3 * (size_t)j, q->p);
	w->cand[w->ncand].j = j;
	w->ncand++;
}

static void rv_nothing(void* user, unsigned int j)
{
	(void)user;
	(void)j;
}

/* at least the k seeds of region g nearest seed i, by distance; reach2: how far every other one is at least */
static void rv_neighbours(rv_cvt* w, double const* z, unsigned int i, unsigned int g, unsigned int k)
{
	rv_near q;
	unsigned int ring, a, b;
	q.w = w;
	q.z = z;
	q.p = z + 3 * (size_t)i;
	q.self = i;
	q.region = g;
	w->ncand = 0;
	w->reach2 = 0;
	for (ring = 0;; ring++)
	{
		unsigned int in = 0;
		if (qaws_internal_grid_ring(w->grid, q.p, ring, rv_near_visit, &q) == HUGE_VAL)
		{
			w->reach2 = HUGE_VAL;
			break;
		}
		/* the seeds not visited are in the rings past this one, at least ring cells away */
		w->reach2 = (ring * w->hcell) * (ring * w->hcell);
		for (a = 0; a < w->ncand; a++)
			in += w->cand[a].d2 <= w->reach2;
		if (in >= k)
			break;   /* (the k nearest are all found) */
	}
	/* by distance, ties by index (insertion: the lists are short) */
	for (a = 1; a < w->ncand; a++)
	{
		rv_cand c = w->cand[a];
		for (b = a; b > 0 && rv_cmp_key(&c, &w->cand[b - 1]) < 0; b--)
			w->cand[b] = w->cand[b - 1];
		w->cand[b] = c;
	}
}

/* polygon `in` cut by the half-space nearer zi than zj */
static unsigned int rv_clip(rv_pt const* in, unsigned int n, rv_pt* out, double const* zi, double const* zj)
{
	double mid[3], dir[3], dp, dq;
	unsigned int a, k, o = 0;
	for (k = 0; k < 3; k++)
	{
		mid[k] = 0.5 * (zi[k] + zj[k]);
		dir[k] = zj[k] - zi[k];
	}
	dq = (in[0].x[0] - mid[0]) * dir[0] + (in[0].x[1] - mid[1]) * dir[1] + (in[0].x[2] - mid[2]) * dir[2];
	for (a = 0; a < n; a++)
	{
		rv_pt const *P = &in[a], *Q = &in[a + 1 < n ? a + 1 : 0];
		dp = dq;
		dq = (Q->x[0] - mid[0]) * dir[0] + (Q->x[1] - mid[1]) * dir[1] + (Q->x[2] - mid[2]) * dir[2];
		if (dp <= 0 && o < RV_POLY)
			out[o++] = *P;
		if ((dp <= 0) != (dq <= 0) && o < RV_POLY)
		{
			double t = dp / (dp - dq);
			for (k = 0; k < 3; k++)
				out[o].x[k] = P->x[k] + t * (Q->x[k] - P->x[k]);
			out[o].rho = P->rho + t * (Q->rho - P->rho);
			o++;
		}
	}
	return o;
}

/* a cut polygon's mass, moment about z and energy (exact for a density linear over it) */
static void rv_accumulate(rv_pt const* p, unsigned int n, double const* z, double* acc)
{
	unsigned int f, k;
	for (f = 1; f + 1 < n; f++)
	{
		rv_pt const* v[3];
		double u[3][3], e1[3], e2[3], nn[3], A, G[3][3], r[3];
		unsigned int a, b;
		v[0] = &p[0]; v[1] = &p[f]; v[2] = &p[f + 1];
		for (k = 0; k < 3; k++)
		{
			e1[k] = v[1]->x[k] - v[0]->x[k];
			e2[k] = v[2]->x[k] - v[0]->x[k];
		}
		rv_cross(e1, e2, nn);
		A = 0.5 * sqrt(rv_dot(nn, nn));
		if (A <= 0)
			continue;
		for (a = 0; a < 3; a++)
		{
			r[a] = v[a]->rho;
			for (k = 0; k < 3; k++)
				u[a][k] = v[a]->x[k] - z[k];
		}
		for (a = 0; a < 3; a++)
			for (b = 0; b < 3; b++)
				G[a][b] = rv_dot(u[a], u[b]);
		acc[0] += A * (r[0] + r[1] + r[2]) / 3;
		for (k = 0; k < 3; k++)
			acc[1 + k] += A / 12 * ((2 * r[0] + r[1] + r[2]) * u[0][k] + (r[0] + 2 * r[1] + r[2]) * u[1][k] + (r[0] + r[1] + 2 * r[2]) * u[2][k]);
		for (a = 0; a < 3; a++)
		{
			unsigned int b1 = (a + 1) % 3, b2 = (a + 2) % 3;
			acc[4] += A / 30 * r[a] * (3 * G[a][a] + G[b1][b1] + G[b2][b2] + 2 * G[a][b1] + 2 * G[a][b2] + G[b1][b2]);
		}
	}
}

/* seed i's cell in region g = ereg[e], breadth first over the sub-triangles it reaches from where it starts (a free seed's
   footpoint; a locked one's entries in g, either side of its creases): acc gets its mass, moment about z_i and energy */
static void rv_cell(rv_cvt* w, double const* z, qaws_mesh_location const* foot, unsigned int i, unsigned int e, double* acc)
{
	rv_pt buf[2][RV_POLY];
	double const* zi = z + 3 * (size_t)i;
	unsigned int head = 0, tail = 0, k = 8, g = w->ereg[e], nstart, q;
	rv_neighbours(w, z, i, g, k);
	if (++w->gen == 0)
	{
		memset(w->stamp, 0, sizeof(unsigned int) * w->nst);
		w->gen = 1;
	}
	for (q = e; q < w->e0[i + 1]; q++)
	{
		unsigned int u = w->locked[i] ? w->estart[q] : rv_locate(w, &foot[i]);
		if (w->ereg[q] == g && u != HM_NONE && w->stamp[u] != w->gen)
		{
			w->stamp[u] = w->gen;
			w->queue[tail++] = u;
		}
	}
	nstart = tail;
	while (head < tail)
	{
		unsigned int u = w->queue[head], n = 3, c, cur = 0, a;
		double R2 = 0;
		int sure = 0;
		for (a = 0; a < 3; a++)
		{
			double const* s = w->sv + 4 * (size_t)w->tri[3 * u + a];
			memcpy(buf[0][a].x, s, sizeof(double) * 3);
			buf[0][a].rho = s[3];
			R2 = fmax(R2, rv_d2(s, zi));
		}
		for (c = 0; c < w->ncand && n >= 3; c++)
		{
			if (w->cand[c].d2 > 4 * R2)
			{
				sure = w->reach2 > 4 * R2;
				break;
			}
			n = rv_clip(buf[cur], n, buf[1 - cur], zi, z + 3 * (size_t)w->cand[c].j);
			cur = 1 - cur;
			R2 = 0;
			for (a = 0; a < n; a++)
				R2 = fmax(R2, rv_d2(buf[cur][a].x, zi));
		}
		if (n >= 3 && c == w->ncand)
			sure = w->reach2 > 4 * R2;
		if (n >= 3 && !sure)
		{
			/* the candidates may miss a bisector: more of them, the sub-triangle again */
			k *= 2;
			rv_neighbours(w, z, i, g, k);
			continue;
		}
		head++;
		if (n < 3 && head > nstart)
			continue;
		if (n >= 3)
			rv_accumulate(buf[cur], n, zi, acc);
		for (a = w->a0[u]; a < w->a0[u + 1]; a++)
			if (w->stamp[w->adj[a]] != w->gen)
			{
				w->stamp[w->adj[a]] = w->gen;
				w->queue[tail++] = w->adj[a];
			}
	}
}

/* every seed's cell at positions z (footpoints foot): mass, centroid and the total energy; 0 on allocation failure */
static int rv_eval(rv_cvt* w, double const* z, qaws_mesh_location const* foot, double* mass, double* cen, double* energy)
{
	unsigned int i, e, q, k;
	double E = 0;
	w->covered = 0;
	for (i = 0; i < w->n; i++)
		for (k = 0; k < 3; k++)
			w->box[i].lo[k] = w->box[i].hi[k] = z[3 * (size_t)i + k];
	if (qaws_internal_grid_create(w->box, w->n, 3, &w->grid) != QAWS_STATUS_OK)
		return 0;
	{
		double lb = qaws_internal_grid_ring(w->grid, z, 2, rv_nothing, NULL);   /* (ring 2 is at least a cell away) */
		w->hcell = lb == HUGE_VAL ? 0 : lb;
	}
	for (i = 0; i < w->n; i++)
	{
		double acc[5] = { 0, 0, 0, 0, 0 };
		/* once per region: its first entry */
		for (e = w->e0[i]; e < w->e0[i + 1]; e++)
		{
			for (q = w->e0[i]; q < e && w->ereg[q] != w->ereg[e]; q++)
			{
			}
			if (q == e)
				rv_cell(w, z, foot, i, e, acc);
		}
		mass[i] = acc[0];
		for (k = 0; k < 3; k++)
			cen[3 * (size_t)i + k] = z[3 * (size_t)i + k] + (acc[0] > 0 ? acc[1 + k] / acc[0] : 0);
		E += acc[4];
		w->covered += acc[0];
	}
	qaws_internal_grid_destroy(w->grid);
	w->grid = NULL;
	*energy = E;
	w->evaluations++;
	return 1;
}

/* the energy's gradient 2 m (z - centroid), in each free seed's tangent plane; its norm */
static double rv_gradient(rv_cvt const* w, double const* z, qaws_mesh_location const* foot, double const* mass, double const* cen, double* g)
{
	unsigned int i, k;
	double s = 0;
	for (i = 0; i < w->n; i++)
	{
		double* gi = g + 3 * (size_t)i;
		double x[3], nr[3], gn;
		if (w->locked[i])
		{
			gi[0] = gi[1] = gi[2] = 0;
			continue;
		}
		for (k = 0; k < 3; k++)
			gi[k] = 2 * mass[i] * (z[3 * (size_t)i + k] - cen[3 * (size_t)i + k]);
		mr_eval(w->r, &foot[i], x, nr);
		gn = rv_dot(gi, nr);
		for (k = 0; k < 3; k++)
			gi[k] -= gn * nr[k];
		s += rv_dot(gi, gi);
	}
	return sqrt(s);
}

/* ---- optimisation ------------------------------------------------------- */

/* Gauss-Newton for the point of patch t nearest y, from (a, b): 1 and its location when it converges inside the patch */
static int rv_newton(qaws_mesh_reference const* r, unsigned int t, double a, double b, double const* y, qaws_mesh_location* l)
{
	unsigned int it, k;
	for (it = 0; it < 12; it++)
	{
		double p[3], pa[3], pb[3], d[3], g0, g1, h00, h01, h11, det, da, db;
		mi_eval(&r->patch[t], a, b, p, pa, pb, NULL, NULL, NULL, NULL, NULL, NULL);
		for (k = 0; k < 3; k++)
			d[k] = p[k] - y[k];
		g0 = rv_dot(pa, d); g1 = rv_dot(pb, d);
		h00 = rv_dot(pa, pa); h01 = rv_dot(pa, pb); h11 = rv_dot(pb, pb);
		det = h00 * h11 - h01 * h01;
		if (det <= 1e-300)
			return 0;
		da = -(h11 * g0 - h01 * g1) / det;
		db = -(h00 * g1 - h01 * g0) / det;
		a += da;
		b += db;
		if (a < 0 || b < 0 || a + b > 1)
			return 0;
		if (fabs(da) + fabs(db) < 1e-13)
		{
			l->patch = t;
			l->a = (qaws_scalar)a;
			l->b = (qaws_scalar)b;
			return 1;
		}
	}
	return 0;
}

/* a free seed onto the reference near y: its position, footpoint and normal there. Nearest y on the footpoint's patch (a
   step's foot is near the last one); off it, walked along the reference by the step's tangent part, within the region, then
   nearest y on the patch reached. 0 when the normal would turn by more than 90 degrees: onto the other side of a thin part */
static int rv_place(rv_cvt* w, unsigned int i, double const* y, double* x, double* n, qaws_mesh_location* l)
{
	qaws_mesh_location const* f = &w->m->foot[i];
	double n0[3], d[3], dn;
	unsigned int k;
	rv_normal(w, i, n0);
	if (!rv_newton(w->r, f->patch, (double)f->a, (double)f->b, y, l))
	{
		qaws_mesh_location q;
		for (k = 0; k < 3; k++)
			d[k] = y[k] - rv_x(w->m, i)[k];
		dn = rv_dot(d, n0);
		for (k = 0; k < 3; k++)
			d[k] -= dn * n0[k];
		mr_walk(w->r, *f, d, 1, l);
		if (rv_newton(w->r, l->patch, (double)l->a, (double)l->b, y, &q))
			*l = q;
	}
	mr_eval(w->r, l, x, n);
	return rv_dot(n, n0) > 0;
}

/* Lloyd: each free seed to its cell's centroid on the reference, unless a face round it would fold */
static void rv_lloyd(rv_cvt* w, double const* mass, double const* cen)
{
	qaws_hm* m = w->m;
	unsigned int i, k;
	for (i = 0; i < w->n; i++)
	{
		double x[3], n[3];
		qaws_mesh_location l;
		if (w->locked[i] || mass[i] <= 0 || !rv_place(w, i, cen + 3 * (size_t)i, x, n, &l) || rv_folds(w, i, x, n, HM_NONE, 0.5))
			continue;
		for (k = 0; k < 3; k++)
			rv_x(m, i)[k] = x[k];
		m->foot[i] = l;
		hm_attr_at(m, w->sa, i)[0] = w->scale * rv_size(w, &l, x);
		w->rep->moves++;
	}
}

/* each seed's entries: a free one its footpoint's region; one on the features, per face round it, the face's region and the
   sub-triangle a little way into the face where its cell starts there (both sides of a crease that has one region round it),
   repeats dropped; 0 on allocation failure */
static int rv_seeds(rv_cvt* w)
{
	qaws_hm* m = w->m;
	unsigned int i, ne = 0, pass;
	w->n = m->nv;
	w->locked = (unsigned char*)qaws_internal_zalloc((size_t)w->n + 1, 1);
	w->e0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)w->n + 1)));
	w->box = (qaws_bp_box*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_bp_box), ((size_t)w->n + 1)));
	w->cand = (rv_cand*)qaws_internal_malloc(qaws_internal_mul2(sizeof(rv_cand), ((size_t)w->n + 1)));
	if (!w->locked || !w->e0 || !w->box || !w->cand)
		return 0;
	for (pass = 0; pass < 2; pass++)
	{
		ne = 0;
		for (i = 0; i < w->n; i++)
		{
			unsigned int c, q;
			w->e0[i] = ne;
			w->locked[i] = 1;   /* (an isolated vertex, or one off the reference, holds no cell) */
			if (m->vc[i] == HM_NONE || m->foot[i].patch >= w->r->nt)
				continue;
			w->locked[i] = hm_feature_degree(m, i) != 0 || (m->vflag[i] & HM_CORNER);   /* (a pinned corner off the lines too) */
			if (!w->locked[i])
			{
				if (pass == 1)
				{
					w->ereg[ne] = w->r->treg[m->foot[i].patch];
					w->estart[ne] = HM_NONE;
				}
				ne++;
				continue;
			}
			for (c = m->vc[i]; c != HM_NONE; c = hm_round(m, i, c))
			{
				unsigned int g = rv_region(w, m->cf[c]), h = m->foot[i].patch, f = m->cf[c], k, s = HM_NONE, d;
				double p[3] = { 0, 0, 0 };
				qaws_mesh_location l;
				if (pass == 0)
				{
					ne++;
					continue;
				}
				for (d = m->fc[f]; d < m->fc[f] + m->fn[f]; d++)
					for (k = 0; k < 3; k++)
						p[k] += rv_x(m, m->cv[d])[k] / m->fn[f];
				for (k = 0; k < 3; k++)
					p[k] = rv_x(m, i)[k] + 0.1 * (p[k] - rv_x(m, i)[k]);
				if (mr_project_region(w->r, p, h < w->r->nt && w->r->treg[h] == g ? h : HM_NONE, g, w->s, &l) != HUGE_VAL && l.patch != HM_NONE)
					s = rv_locate(w, &l);
				for (q = w->e0[i]; q < ne && (w->ereg[q] != g || w->estart[q] != s); q++)
				{
				}
				if (q < ne)
					continue;
				w->ereg[ne] = g;
				w->estart[ne] = s;
				ne++;
			}
		}
		w->e0[w->n] = ne;
		if (pass == 0)
		{
			w->ereg = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)ne + 1)));
			w->estart = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)ne + 1)));
			if (!w->ereg || !w->estart)
				return 0;
		}
	}
	return 1;
}

/* Lloyd steps, then L-BFGS (the two-loop recursion of Nocedal and Wright, H0 the cells' 1 / (2 m), a Lloyd step, scaled by
   s.y / y.H0.y); each seed's step at most half its size, no seed moved where a face round it would fold, Armijo backtracking
   on the energy after projection, the history dropped when a projection bends a step, Lloyd when the search fails; flips to
   Delaunay along the way; 0 on allocation failure */
static int rv_solve(rv_cvt* w, unsigned int steps)
{
	qaws_hm* m = w->m;
	unsigned int n3 = 3 * w->n, i, k, it, hn = 0, hhead = 0;
	size_t sz = sizeof(double) * ((size_t)n3 + 3);
	double *mass = (double*)qaws_internal_scratch_malloc(sz), *cen = (double*)qaws_internal_malloc(sz), *mt = (double*)qaws_internal_malloc(sz), *ct = (double*)qaws_internal_malloc(sz);
	double *g = (double*)qaws_internal_malloc(sz), *gt = (double*)qaws_internal_malloc(sz), *d = (double*)qaws_internal_malloc(sz), *zt = (double*)qaws_internal_malloc(sz);
	double *S = (double*)qaws_internal_malloc(qaws_internal_mul2(sz, RV_HISTORY)), *Y = (double*)qaws_internal_malloc(qaws_internal_mul2(sz, RV_HISTORY)), rh[RV_HISTORY], al[RV_HISTORY];
	qaws_mesh_location* ft = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), ((size_t)w->n + 1)));
	double E = 0, Et = 0, gnorm, gam = 1, stop = 1e-9 * w->r->extent * w->r->extent * w->r->extent;
	int ok = mass && cen && mt && ct && g && gt && d && zt && S && Y && ft && w->n <= 0xFFFFFFFFu / 3;   /* n3 counts in unsigned int */
	for (it = 0; ok && it < 5; it++)
	{
		ok = rv_eval(w, m->x, m->foot, mass, cen, &E);
		if (ok)
		{
			rv_lloyd(w, mass, cen);
			ok = rv_delaunay(w);
		}
	}
	ok = ok && rv_eval(w, m->x, m->foot, mass, cen, &E);
	if (ok && rv_trace)
		rv_trace->lloyd_energy = E;
	gnorm = ok ? rv_gradient(w, m->x, m->foot, mass, cen, g) : 0;
	for (it = 0; ok && it < steps && gnorm >= stop; it++)
	{
		double gd = 0, alpha = 1, far = 0;
		unsigned int tries, moved = 0;
		int accepted = 0;
		/* the direction */
		memcpy(d, g, sizeof(double) * n3);
		for (k = 0; k < hn; k++)
		{
			unsigned int h = (hhead + RV_HISTORY - 1 - k) % RV_HISTORY;
			double s = 0;
			for (i = 0; i < n3; i++)
				s += S[(size_t)h * n3 + i] * d[i];
			al[h] = rh[h] * s;
			for (i = 0; i < n3; i++)
				d[i] -= al[h] * Y[(size_t)h * n3 + i];
		}
		for (i = 0; i < n3; i++)
			d[i] = mass[i / 3] > 0 ? (hn ? gam : 1) * d[i] / (2 * mass[i / 3]) : 0;
		for (k = hn; k > 0; k--)
		{
			unsigned int h = (hhead + RV_HISTORY - k) % RV_HISTORY;
			double s = 0;
			for (i = 0; i < n3; i++)
				s += Y[(size_t)h * n3 + i] * d[i];
			s *= rh[h];
			for (i = 0; i < n3; i++)
				d[i] += S[(size_t)h * n3 + i] * (al[h] - s);
		}
		for (i = 0; i < n3; i++)
		{
			d[i] = w->locked[i / 3] ? 0 : -d[i];
			gd += g[i] * d[i];
		}
		if (gd >= 0)
		{
			/* not a descent: the history dropped, a Lloyd direction */
			hn = 0;
			for (i = 0; i < n3; i++)
				d[i] = w->locked[i / 3] || mass[i / 3] <= 0 ? 0 : -g[i] / (2 * mass[i / 3]);
		}
		/* each seed's step at most half its size */
		gd = 0;
		for (i = 0; i < w->n; i++)
		{
			double l = sqrt(rv_dot(d + 3 * (size_t)i, d + 3 * (size_t)i)), hm = 0.5 * rv_h(w, i);
			for (k = 0; k < 3; k++)
			{
				if (l > hm)
					d[3 * (size_t)i + k] *= hm / l;
				gd += g[3 * (size_t)i + k] * d[3 * (size_t)i + k];
			}
		}
		/* the trial in place, the seeds moved one by one, each kept where it is when a face round it would fold; zt and ft
		   hold where they were */
		memcpy(zt, m->x, sizeof(double) * n3);
		memcpy(ft, m->foot, sizeof(qaws_mesh_location) * w->n);
		for (tries = 0; tries < 10 && gd < 0; tries++, alpha *= 0.5)
		{
			far = 0;
			moved = 0;
			memcpy(m->x, zt, sizeof(double) * n3);
			memcpy(m->foot, ft, sizeof(qaws_mesh_location) * w->n);
			for (i = 0; i < w->n; i++)
			{
				double y[3], x[3], nr[3];
				qaws_mesh_location l;
				if (w->locked[i] || (d[3 * (size_t)i] == 0 && d[3 * (size_t)i + 1] == 0 && d[3 * (size_t)i + 2] == 0))
					continue;
				for (k = 0; k < 3; k++)
					y[k] = zt[3 * (size_t)i + k] + alpha * d[3 * (size_t)i + k];
				if (!rv_place(w, i, y, x, nr, &l) || rv_folds(w, i, x, nr, HM_NONE, 0.5))
					continue;
				far = fmax(far, sqrt(rv_d2(y, x)) / rv_h(w, i));
				memcpy(rv_x(m, i), x, sizeof x);
				m->foot[i] = l;
				moved++;
			}
			if (!rv_eval(w, m->x, m->foot, mt, ct, &Et))
			{
				ok = 0;
				break;
			}
			if (Et <= E + 1e-4 * alpha * gd)
			{
				accepted = 1;
				break;
			}
		}
		if (!ok)
			break;
		if (accepted)
		{
			double sy = 0, ss = 0, yy = 0;
			rv_gradient(w, m->x, m->foot, mt, ct, gt);
			for (i = 0; i < n3; i++)
			{
				double si = m->x[i] - zt[i], yi = gt[i] - g[i];
				sy += si * yi;
				ss += si * si;
				yy += mt[i / 3] > 0 ? yi * yi / (2 * mt[i / 3]) : 0;
			}
			if (far > 0.1)
				hn = 0;   /* (the projection bent the step: the history no longer holds) */
			else if (sy > 1e-12 * ss)
			{
				/* the pair into the oldest slot, only now: with the history full it is still in use */
				for (i = 0; i < n3; i++)
				{
					S[(size_t)hhead * n3 + i] = m->x[i] - zt[i];
					Y[(size_t)hhead * n3 + i] = gt[i] - g[i];
				}
				rh[hhead] = 1 / sy;
				gam = yy > 0 ? fmin(50, fmax(0.2, sy / yy)) : 1;
				hhead = (hhead + 1) % RV_HISTORY;
				hn = hn < RV_HISTORY ? hn + 1 : RV_HISTORY;
			}
			for (i = 0; i < w->n; i++)
				if (!w->locked[i])
					hm_attr_at(m, w->sa, i)[0] = w->scale * rv_size(w, &m->foot[i], rv_x(m, i));
			w->rep->moves += moved;
			E = Et;
			memcpy(mass, mt, sizeof(double) * w->n);
			memcpy(cen, ct, sizeof(double) * n3);
			memcpy(g, gt, sizeof(double) * n3);
		}
		else
		{
			hn = 0;
			memcpy(m->x, zt, sizeof(double) * n3);
			memcpy(m->foot, ft, sizeof(qaws_mesh_location) * w->n);
			rv_lloyd(w, mass, cen);
			ok = rv_eval(w, m->x, m->foot, mass, cen, &E);
			if (!ok)
				break;
			rv_gradient(w, m->x, m->foot, mass, cen, g);
		}
		gnorm = 0;
		for (i = 0; i < n3; i++)
			gnorm += g[i] * g[i];
		gnorm = sqrt(gnorm);
		if ((it + 1) % 5 == 0)
			ok = rv_delaunay(w);
	}
	ok = ok && rv_delaunay(w);
	if (ok && rv_trace)
	{
		rv_trace->energy = E;
		rv_trace->coverage = w->mass > 0 ? w->covered / w->mass : 0;
		rv_trace->evaluations = w->evaluations;
	}
	qaws_internal_free(mass); qaws_internal_free(cen); qaws_internal_free(mt); qaws_internal_free(ct); qaws_internal_free(g); qaws_internal_free(gt); qaws_internal_free(d); qaws_internal_free(zt); qaws_internal_free(S); qaws_internal_free(Y); qaws_internal_free(ft);
	return ok;
}

static void rv_seeds_free(rv_cvt* w)
{
	qaws_internal_free(w->locked); qaws_internal_free(w->e0); qaws_internal_free(w->ereg); qaws_internal_free(w->estart); qaws_internal_free(w->box); qaws_internal_free(w->cand);
	w->locked = NULL; w->e0 = w->ereg = w->estart = NULL; w->box = NULL; w->cand = NULL;
}

/* the size attribute from the footpoints, times s */
static void rv_resize(rv_cvt* w, double s)
{
	qaws_hm* m = w->m;
	unsigned int v;
	w->scale = s;
	for (v = 0; v < m->nv; v++)
		if (m->foot[v].patch < w->r->nt)
			hm_attr_at(m, w->sa, v)[0] = s * rv_size(w, &m->foot[v], rv_x(m, v));
}

/* one level: the feature lines from the size, the count, Delaunay, the seeds solved */
static int rv_level(rv_cvt* w, double nstar, int coarse, unsigned int steps)
{
	int ok = rv_lines(w) && rv_count(w, nstar, coarse);
	ok = ok && rv_delaunay(w) && rv_seeds(w) && rv_solve(w, steps);
	rv_seeds_free(w);
	return ok;
}

/* the integration mesh freed */
static void rv_unbuild(rv_cvt* w)
{
	qaws_internal_free(w->k); qaws_internal_free(w->v0); qaws_internal_free(w->t0); qaws_internal_free(w->sv); qaws_internal_free(w->tri); qaws_internal_free(w->a0); qaws_internal_free(w->adj); qaws_internal_free(w->stamp); qaws_internal_free(w->queue);
	w->k = w->v0 = w->t0 = w->tri = w->a0 = w->adj = w->stamp = w->queue = NULL;
	w->sv = NULL;
	w->gen = 0;
}

int rv_optimize(qaws_hm* m, qaws_mesh_reference const* r, qaws_size_field const* size, double L, unsigned int size_attr, unsigned int region_attr,
	mr_scratch* s, unsigned int steps, qaws_mesh_report* rep)
{
	rv_cvt w;
	double nstar = 0;
	int ok;
	memset(&w, 0, sizeof w);
	w.m = m;
	w.r = r;
	w.size = size;
	w.L = L;
	w.sa = size_attr;
	w.ra = region_attr;
	w.s = s;
	w.rep = rep;
	w.scale = 2;
	ok = rv_build(&w, &nstar);
	/* a coarse level first (a quarter of the seeds at twice the size, on sub-triangles twice as large: dislocations anneal four
	   times sooner), each triangle then split in four, which keeps every vertex's valence; dropped, the working mesh put back,
	   where the split folds faces: a part curved tighter than twice the size */
	if (ok && nstar >= 200)
	{
		qaws_hm keep;
		unsigned int before = rv_folded(&w);
		ok = hm_copy(&keep, m) == QAWS_STATUS_OK;
		if (ok)
		{
			rv_resize(&w, 2);
			ok = rv_level(&w, nstar / 4, 1, steps);
			rv_resize(&w, 1);
			ok = ok && rv_subdivide(&w);
			if (ok && rv_folded(&w) > before)
			{
				hm_free(m);
				*m = keep;
			}
			else
			{
				hm_free(&keep);
				w.slack = 0.05;   /* (the split doubled the coarse lines: kept within 5% of the size, the split's structure with them) */
			}
		}
	}
	rv_unbuild(&w);
	w.scale = 1;
	ok = ok && rv_build(&w, &nstar) && rv_level(&w, nstar, 0, steps);
	rv_unbuild(&w);
	qaws_internal_free(w.eq);
	return ok;
}
