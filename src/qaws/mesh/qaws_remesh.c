/*
 * Local remeshing on a reference (Botsch and Kobbelt 2004): long edges
 * split at their midpoint on the reference, short ones collapsed, flips
 * toward valence 6, vertices relaxed along the surface by walking; feature
 * corners fixed, feature lines and the border slid along; footpoints kept
 * throughout; the reference mesh's attributes transferred at the end.
 */

#include "qaws_mesh_internal.h"
#include "../qaws_remesh.h"
#include "qaws_remesh_cvt.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct rmsh
{
	qaws_hm* m;                     /* the working mesh (the output's store) */
	qaws_mesh_reference const* r;
	qaws_size_field const* size;
	double L;                       /* the uniform target */
	unsigned int sa;                /* the size attribute (per vertex) */
	unsigned int ra;                /* the region attribute (per face): the reference's triangles between creases the face is on */
	mr_scratch s;
	qaws_mesh_report rep;
	unsigned int* ep;               /* the edges as vertex pairs at a pass's start */
	unsigned int nep, cep;
	int failed;
} rmsh;

static double rm_d2(double const* a, double const* b)
{
	return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
}

static double* rm_x(qaws_hm const* m, unsigned int v)
{
	return m->x + 3 * (size_t)v;
}

static double rm_size_of(rmsh const* w, unsigned int v)
{
	return hm_attr_at(w->m, w->sa, v)[0];
}

static double rm_target(rmsh const* w, unsigned int a, unsigned int b)
{
	return 0.5 * (rm_size_of(w, a) + rm_size_of(w, b));
}

/* the size at v from its footpoint */
static void rm_set_size(rmsh* w, unsigned int v)
{
	hm_attr_at(w->m, w->sa, v)[0] = w->size ? sf_eval(w->size, &w->m->foot[v], rm_x(w->m, v)) : w->L;
}

/* a triangle's shape: 4 sqrt(3) area over the sum of its squared edges (1 equilateral, 0 flat) */
static double rm_quality(double const* A, double const* B, double const* C)
{
	double u[3], t[3], n[3], s = 0;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = B[k] - A[k];
		t[k] = C[k] - A[k];
		s += u[k] * u[k] + t[k] * t[k] + (C[k] - B[k]) * (C[k] - B[k]);
	}
	n[0] = u[1] * t[2] - u[2] * t[1]; n[1] = u[2] * t[0] - u[0] * t[2]; n[2] = u[0] * t[1] - u[1] * t[0];
	return s > 0 ? 2 * sqrt(3.0) * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) / s : 0;
}

#define RM_FAIR 0.3   /* below this a triangle is poor: no operation may leave the poorest one round it poorer */

/* would moving v to y turn over a face round v (faces holding `skip` aside), or
   leave the poorest face poorer than before and than fair? */
static int rm_folds(qaws_hm const* m, unsigned int v, double const* y, unsigned int skip)
{
	unsigned int c;
	double qb = 1, qa = 1;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
	{
		unsigned int a = hm_to(m, c), b = m->cv[hm_prev(m, c)], k;
		double const *p = rm_x(m, v), *A = rm_x(m, a), *B = rm_x(m, b);
		double u[3], t[3], n0[3], n1[3], d0, d1, dd;
		if (a == skip || b == skip)
			continue;
		for (k = 0; k < 3; k++) { u[k] = A[k] - p[k]; t[k] = B[k] - p[k]; }
		n0[0] = u[1] * t[2] - u[2] * t[1]; n0[1] = u[2] * t[0] - u[0] * t[2]; n0[2] = u[0] * t[1] - u[1] * t[0];
		for (k = 0; k < 3; k++) { u[k] = A[k] - y[k]; t[k] = B[k] - y[k]; }
		n1[0] = u[1] * t[2] - u[2] * t[1]; n1[1] = u[2] * t[0] - u[0] * t[2]; n1[2] = u[0] * t[1] - u[1] * t[0];
		d0 = n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2];
		d1 = n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2];
		dd = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2];
		if (dd <= 0.2 * sqrt(d0 * d1) || d1 < 1e-12 * d0)
			return 1;   /* turned by more than about 78 degrees, or flattened */
		qb = fmin(qb, rm_quality(p, A, B));
		qa = fmin(qa, rm_quality(y, A, B));
	}
	return qa < qb && qa < RM_FAIR;
}

/* the region of face f */
static unsigned int rm_region(rmsh const* w, unsigned int f)
{
	return (unsigned int)hm_attr_at(w->m, w->ra, f)[0];
}

/* an end's footpoint triangle in the region as a hint (HM_NONE: none is) */
static unsigned int rm_hint(rmsh const* w, unsigned int a, unsigned int b, unsigned int region)
{
	unsigned int ta = w->m->foot[a].patch, tb = w->m->foot[b].patch;
	if (ta < w->r->nt && w->r->treg[ta] == region)
		return ta;
	return tb < w->r->nt && w->r->treg[tb] == region ? tb : HM_NONE;
}

/* the edges as vertex pairs, taken at a pass's start: an edit moves the corners of the faces it
   touches past the end of the store, so walking the corners would skip their edges */
static int rm_snapshot(rmsh* w)
{
	qaws_hm const* m = w->m;
	unsigned int c;
	if (m->nc > w->cep)
	{
		unsigned int* g = (unsigned int*)realloc(w->ep, sizeof(unsigned int) * 2 * m->nc);
		if (!g)
		{
			w->failed = 1;
			return 0;
		}
		w->ep = g;
		w->cep = m->nc;
	}
	w->nep = 0;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			w->ep[2 * w->nep] = m->cv[c];
			w->ep[2 * w->nep + 1] = hm_to(m, c);
			w->nep++;
		}
	return 1;
}

/* the corner of snapshot edge i now, or HM_NONE when it is gone */
static unsigned int rm_edge(rmsh const* w, unsigned int i)
{
	unsigned int a = w->ep[2 * i], b = w->ep[2 * i + 1];
	if ((w->m->vflag[a] | w->m->vflag[b]) & HM_DELETED)
		return HM_NONE;
	return hm_find_edge(w->m, a, b);
}

/* ---- splits ------------------------------------------------------------- */

static void rm_split(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int a, b, v, k, c = rm_edge(w, i);
		double mid[3], x[3], d;
		qaws_mesh_location l;
		if (c == HM_NONE)
			continue;
		a = m->cv[c];
		b = hm_to(m, c);
		if (rm_d2(rm_x(m, a), rm_x(m, b)) <= (16.0 / 9) * rm_target(w, a, b) * rm_target(w, a, b))
			continue;
		for (k = 0; k < 3; k++)
			mid[k] = 0.5 * (rm_x(m, a)[k] + rm_x(m, b)[k]);
		/* the midpoint on the reference: on the creases for a feature edge, on the hint's side otherwise */
		if (m->cflag[c] & HM_FEATURE)
			d = mr_project_crease(w->r, mid, m->foot[a].patch, &w->s, &l);
		else
		{
			unsigned int reg = rm_region(w, m->cf[c]);
			d = mr_project_region(w->r, mid, rm_hint(w, a, b, reg), reg, &w->s, &l);
		}
		if (d == HUGE_VAL || l.patch == HM_NONE)
			continue;
		mr_eval(w->r, &l, x, NULL);
		v = hm_split_edge(m, c, 0.5, x, 1);
		if (v == HM_NONE)
		{
			w->failed = 1;
			break;
		}
		m->foot[v] = l;
		rm_set_size(w, v);
		w->rep.splits++;
	}
}

/* ---- collapses ---------------------------------------------------------- */

/* may x, an end of c's edge, go into the other end y? */
static int rm_along_crease(rmsh* w, unsigned int p, unsigned int q);

static int rm_can_remove(rmsh* w, unsigned int c, unsigned int x, unsigned int y)
{
	qaws_hm const* m = w->m;
	unsigned int ring[256], n, i;
	double hi;
	if (m->vflag[x] & HM_CORNER)
		return 0;
	if (hm_feature_degree(m, x) && !(m->cflag[c] & HM_FEATURE))
		return 0;
	if (!hm_can_collapse(m, c, x))
		return 0;
	n = hm_ring(m, x, ring, 256);
	if (n > 256)
		return 0;
	for (i = 0; i < n; i++)
	{
		if (ring[i] == y)
			continue;
		hi = (4.0 / 3) * rm_target(w, y, ring[i]);
		if (rm_d2(rm_x(m, y), rm_x(m, ring[i])) > hi * hi)
			return 0;
		/* no new edge along a crease over its own vertices */
		if (hm_find_edge(m, y, ring[i]) == HM_NONE && rm_along_crease(w, y, ring[i]))
			return 0;
	}
	return !rm_folds(m, x, rm_x(m, y), y);
}

static void rm_collapse(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int a, b, x, y, c = rm_edge(w, i);
		double lo;
		if (c == HM_NONE)
			continue;
		a = m->cv[c];
		b = hm_to(m, c);
		lo = 0.8 * rm_target(w, a, b);
		if (rm_d2(rm_x(m, a), rm_x(m, b)) >= lo * lo)
			continue;
		/* keep the end that is on a feature, else the one of higher valence */
		if (hm_feature_degree(m, a) > hm_feature_degree(m, b) || (!hm_feature_degree(m, a) && !hm_feature_degree(m, b) && hm_valence(m, a) > hm_valence(m, b)))
		{
			x = b;
			y = a;
		}
		else
		{
			x = a;
			y = b;
		}
		if (!rm_can_remove(w, c, x, y))
		{
			unsigned int t = x;
			x = y;
			y = t;
			if (!rm_can_remove(w, c, x, y))
				continue;
		}
		if (!hm_collapse(m, c, x, 1, rm_x(m, y)))
		{
			w->failed = 1;
			break;
		}
		w->rep.collapses++;
	}
}

/* ---- flips -------------------------------------------------------------- */

static int rm_ideal(qaws_hm const* m, unsigned int v)
{
	return hm_is_border_vertex(m, v) ? 4 : 6;
}

/* would a new edge from p to q run along a crease (both ends on features, its
   middle on one): a chord over a feature line's own vertices? */
static int rm_along_crease(rmsh* w, unsigned int p, unsigned int q)
{
	qaws_hm const* m = w->m;
	double mid[3], d2 = 0;
	qaws_mesh_location l;
	unsigned int k;
	if (!hm_feature_degree(m, p) || !hm_feature_degree(m, q))
		return 0;
	for (k = 0; k < 3; k++)
	{
		mid[k] = 0.5 * (rm_x(m, p)[k] + rm_x(m, q)[k]);
		d2 += (rm_x(m, p)[k] - rm_x(m, q)[k]) * (rm_x(m, p)[k] - rm_x(m, q)[k]);
	}
	return mr_project_crease(w->r, mid, m->foot[p].patch, &w->s, &l) < 1e-6 * d2;
}

static void rm_flip(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int d, a, b, p, q, k, c = rm_edge(w, i);
		int va, vb, vp, vq, before, after;
		double const *A, *B, *P, *Q;
		double n1[3], n2[3], n0[3], u[3], t[3];
		if (c == HM_NONE || (d = m->co[c]) == HM_NONE || (m->cflag[c] & HM_FEATURE) || !hm_can_flip(m, c))
			continue;
		a = m->cv[c]; b = hm_to(m, c); p = m->cv[hm_prev(m, c)]; q = m->cv[hm_prev(m, d)];
		if (rm_along_crease(w, p, q))
			continue;
		va = (int)hm_valence(m, a) - rm_ideal(m, a);
		vb = (int)hm_valence(m, b) - rm_ideal(m, b);
		vp = (int)hm_valence(m, p) - rm_ideal(m, p);
		vq = (int)hm_valence(m, q) - rm_ideal(m, q);
		before = abs(va) + abs(vb) + abs(vp) + abs(vq);
		after = abs(va - 1) + abs(vb - 1) + abs(vp + 1) + abs(vq + 1);
		if (after >= before)
			continue;
		/* the new pair (a, q, p) and (q, b, p) facing the old pair's way */
		A = rm_x(m, a); B = rm_x(m, b); P = rm_x(m, p); Q = rm_x(m, q);
		for (k = 0; k < 3; k++) { u[k] = Q[k] - A[k]; t[k] = P[k] - A[k]; }
		n1[0] = u[1] * t[2] - u[2] * t[1]; n1[1] = u[2] * t[0] - u[0] * t[2]; n1[2] = u[0] * t[1] - u[1] * t[0];
		for (k = 0; k < 3; k++) { u[k] = B[k] - Q[k]; t[k] = P[k] - Q[k]; }
		n2[0] = u[1] * t[2] - u[2] * t[1]; n2[1] = u[2] * t[0] - u[0] * t[2]; n2[2] = u[0] * t[1] - u[1] * t[0];
		hm_face_normal(m, m->cf[c], n0);
		{
			double g[3];
			hm_face_normal(m, m->cf[d], g);
			for (k = 0; k < 3; k++) n0[k] += g[k];
		}
		if (n1[0] * n0[0] + n1[1] * n0[1] + n1[2] * n0[2] <= 0 || n2[0] * n0[0] + n2[1] * n0[1] + n2[2] * n0[2] <= 0
			|| n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] <= 0)
			continue;
		{
			double qo = fmin(rm_quality(A, B, P), rm_quality(B, A, Q)), qn = fmin(rm_quality(A, Q, P), rm_quality(Q, B, P));
			if (qn < qo && qn < RM_FAIR)
				continue;
		}
		if (!hm_flip(m, c))
		{
			w->failed = 1;
			break;
		}
		w->rep.flips++;
	}
}

/* the cosine of the angle at P of triangle (P, A, B) */
static double rm_cos(double const* P, double const* A, double const* B)
{
	double u[3], t[3], uu = 0, tt = 0, ut = 0;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = A[k] - P[k];
		t[k] = B[k] - P[k];
		uu += u[k] * u[k];
		tt += t[k] * t[k];
		ut += u[k] * t[k];
	}
	return uu > 0 && tt > 0 ? ut / sqrt(uu * tt) : 1;
}

/* caps (an angle over 150 degrees): the edge across it flipped when that makes the pair better (as pmp's remove_caps) */
static void rm_caps(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int i;
	double lim = cos(150 * 3.14159265358979323846 / 180);
	for (i = 0; (i > 0 || rm_snapshot(w)) && i < w->nep && !w->failed; i++)
	{
		unsigned int d, c = rm_edge(w, i);
		double const *A, *B, *P, *Q;
		if (c == HM_NONE || (d = m->co[c]) == HM_NONE || (m->cflag[c] & HM_FEATURE) || !hm_can_flip(m, c)
			|| rm_along_crease(w, m->cv[hm_prev(m, c)], m->cv[hm_prev(m, d)]))
			continue;
		A = rm_x(m, m->cv[c]); B = rm_x(m, hm_to(m, c)); P = rm_x(m, m->cv[hm_prev(m, c)]); Q = rm_x(m, m->cv[hm_prev(m, d)]);
		if (rm_cos(P, A, B) > lim && rm_cos(Q, B, A) > lim)
			continue;
		if (fmin(rm_quality(A, Q, P), rm_quality(Q, B, P)) <= fmin(rm_quality(A, B, P), rm_quality(B, A, Q)))
			continue;
		if (!hm_flip(m, c))
		{
			w->failed = 1;
			break;
		}
		w->rep.flips++;
	}
}

/* ---- tangential relaxation ---------------------------------------------- */

static void rm_smooth(rmsh* w)
{
	qaws_hm* m = w->m;
	unsigned int v, k;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int ring[256], n, i, fd;
		double cen[3] = { 0, 0, 0 }, d[3], x[3], nrm[3], dn;
		qaws_mesh_location l;
		if ((m->vflag[v] & (HM_DELETED | HM_CORNER)) || m->vc[v] == HM_NONE)
			continue;
		fd = hm_feature_degree(m, v);
		if (fd != 0 && fd != 2)
			continue;
		n = hm_ring(m, v, ring, 256);
		if (n == 0 || n > 256)
			continue;
		if (fd == 0)
		{
			for (i = 0; i < n; i++)
				for (k = 0; k < 3; k++)
					cen[k] += rm_x(m, ring[i])[k] / n;
			mr_eval(w->r, &m->foot[v], x, nrm);
			for (k = 0; k < 3; k++)
				d[k] = 0.5 * (cen[k] - rm_x(m, v)[k]);
			dn = d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2];
			for (k = 0; k < 3; k++)
				d[k] -= dn * nrm[k];
			if (mr_walk(w->r, m->foot[v], d, 1, &l))
				continue;   /* (it would reach a crease: stay off it) */
		}
		else
		{
			/* along the feature line: toward the middle of its two feature neighbours */
			unsigned int f[2], nf = 0, c, last = HM_NONE;
			double t[3], dt, y[3];
			for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
			{
				if ((m->cflag[c] & HM_FEATURE) && nf < 2)
					f[nf++] = hm_to(m, c);
				last = c;
			}
			if (nf < 2 && last != HM_NONE && hm_is_border_vertex(m, v) && (m->cflag[hm_prev(m, last)] & HM_FEATURE))
				f[nf++] = m->cv[hm_prev(m, last)];
			if (nf != 2)
				continue;
			for (k = 0; k < 3; k++)
			{
				t[k] = rm_x(m, f[1])[k] - rm_x(m, f[0])[k];
				d[k] = 0.5 * (0.5 * (rm_x(m, f[0])[k] + rm_x(m, f[1])[k]) - rm_x(m, v)[k]);
			}
			dt = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
			if (dt <= 0)
				continue;
			dt = (d[0] * t[0] + d[1] * t[1] + d[2] * t[2]) / dt;
			for (k = 0; k < 3; k++)
				y[k] = rm_x(m, v)[k] + dt * t[k];
			if (mr_project_crease(w->r, y, m->foot[v].patch, &w->s, &l) == HUGE_VAL || l.patch == HM_NONE)
				continue;
		}
		mr_eval(w->r, &l, x, NULL);
		if (rm_folds(m, v, x, HM_NONE))
			continue;
		for (k = 0; k < 3; k++)
			rm_x(m, v)[k] = x[k];
		m->foot[v] = l;
		rm_set_size(w, v);
		w->rep.moves++;
	}
}

/* ---- set-up and the end ------------------------------------------------- */

/* the working mesh: the input's faces fanned, its features and the border as features, footpoints */
static qaws_status rm_begin(rmsh* w, qaws_mesh_3d const* in, qaws_mesh_3d** out)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(in), nt = qaws_mesh_3d_get_triangle_count(in), nfe = qaws_mesh_3d_get_feature_edge_count(in), i, k;
	unsigned int *tri = (unsigned int*)malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1)), *fe = (unsigned int*)malloc(sizeof(unsigned int) * (2 * (size_t)nfe + 1));
	qaws_hm* m;
	qaws_status st;
	if (!tri || !fe)
	{
		free(tri);
		free(fe);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	qaws_mesh_3d_get_triangles(in, tri);
	qaws_mesh_3d_get_feature_edges(in, fe);
	st = hm_mesh_3d_create(in->hm.x, nv, tri, NULL, nt, out);
	free(tri);
	if (st != QAWS_STATUS_OK)
	{
		free(fe);
		return st;
	}
	m = &(*out)->hm;
	w->m = m;
	for (i = 0; i < nfe; i++)
	{
		unsigned int c = hm_find_edge(m, fe[2 * i], fe[2 * i + 1]);
		if (c == HM_NONE)
			continue;
		m->cflag[c] |= HM_FEATURE;
		if (m->co[c] != HM_NONE)
			m->cflag[m->co[c]] |= HM_FEATURE;
	}
	free(fe);
	for (i = 0; i < m->nc; i++)
		if (m->cf[i] != HM_NONE && m->co[i] == HM_NONE)
			m->cflag[i] |= HM_FEATURE;   /* (the border) */
	/* corners: the input's, where feature lines meet or end, and where one turns by more than 60 degrees */
	for (i = 0; i < m->nv; i++)
	{
		unsigned int fd = hm_feature_degree(m, i), f[2], nf = 0, c, last = HM_NONE;
		if ((i < in->hm.nv && (in->hm.vflag[i] & HM_CORNER)) || (fd != 0 && fd != 2))
		{
			m->vflag[i] |= HM_CORNER;
			continue;
		}
		if (fd != 2)
			continue;
		for (c = m->vc[i]; c != HM_NONE; c = hm_round(m, i, c))
		{
			if ((m->cflag[c] & HM_FEATURE) && nf < 2)
				f[nf++] = hm_to(m, c);
			last = c;
		}
		if (nf < 2 && last != HM_NONE && hm_is_border_vertex(m, i) && (m->cflag[hm_prev(m, last)] & HM_FEATURE))
			f[nf++] = m->cv[hm_prev(m, last)];
		if (nf == 2)
		{
			double u[3], v[3], uu = 0, vv = 0, uv = 0;
			for (k = 0; k < 3; k++)
			{
				u[k] = rm_x(m, i)[k] - rm_x(m, f[0])[k];
				v[k] = rm_x(m, f[1])[k] - rm_x(m, i)[k];
				uu += u[k] * u[k];
				vv += v[k] * v[k];
				uv += u[k] * v[k];
			}
			if (uv < 0.5 * sqrt(uu * vv))
				m->vflag[i] |= HM_CORNER;
		}
	}
	if (!hm_foot_enable(m) || !mr_scratch_init(&w->s, w->r->nt))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	w->sa = hm_attr_add(m, "qaws.size", QAWS_MESH_VERTEX, QAWS_ATTRIBUTE_LINEAR, 1);
	w->ra = hm_attr_add(m, "qaws.region", QAWS_MESH_FACE, QAWS_ATTRIBUTE_NEAREST, 1);
	if (w->sa == HM_NONE || w->ra == HM_NONE)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* each face's region: its own triangle's when the reference is the input's fan, else under its centre */
	for (i = 0; i < m->nf; i++)
	{
		unsigned int reg = HM_NONE;
		if (w->r->mesh->hm.nv == in->hm.nv && w->r->nt == m->nf && w->r->tv[3 * i] == m->cv[m->fc[i]] && w->r->tv[3 * i + 1] == m->cv[m->fc[i] + 1]
			&& w->r->tv[3 * i + 2] == m->cv[m->fc[i] + 2])
			reg = w->r->treg[i];
		else
		{
			double cen[3] = { 0, 0, 0 };
			qaws_mesh_location l;
			unsigned int c;
			for (c = m->fc[i]; c < m->fc[i] + 3; c++)
				for (k = 0; k < 3; k++)
					cen[k] += rm_x(m, m->cv[c])[k] / 3;
			mr_project(w->r, cen, HM_NONE, &w->s, &l);
			reg = l.patch < w->r->nt ? w->r->treg[l.patch] : 0;
		}
		hm_attr_at(m, w->ra, i)[0] = reg;
	}
	/* footpoints: the input's when it has them, else projected (feature vertices onto the creases) */
	for (i = 0; i < m->nv; i++)
	{
		qaws_mesh_location l;
		double x[3];
		unsigned int src = m->vsrc[i];
		if (in->hm.foot && src < in->hm.nv && in->hm.foot[src].patch < w->r->nt)
			l = in->hm.foot[src];
		else if (hm_feature_degree(m, i) && mr_project_crease(w->r, rm_x(m, i), HM_NONE, &w->s, &l) != HUGE_VAL && l.patch != HM_NONE)
		{
		}
		else if (m->vc[i] != HM_NONE)
			mr_project_region(w->r, rm_x(m, i), HM_NONE, rm_region(w, m->cf[m->vc[i]]), &w->s, &l);
		else
			mr_project(w->r, rm_x(m, i), HM_NONE, &w->s, &l);
		if (l.patch == HM_NONE)
			return QAWS_STATUS_INVALID_ARGUMENT;
		mr_eval(w->r, &l, x, NULL);
		for (k = 0; k < 3; k++)
			rm_x(m, i)[k] = x[k];
		m->foot[i] = l;
		rm_set_size(w, i);
	}
	return QAWS_STATUS_OK;
}

/* barycentric coordinates of p in triangle t's flat corners, clamped into it */
static void rm_bary(qaws_mesh_reference const* r, unsigned int t, double const* p, double* a, double* b)
{
	double const *P0 = r->patch[t].p[0], *P1 = r->patch[t].p[1], *P2 = r->patch[t].p[2];
	double d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0, den, u, v;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		double e1 = P1[k] - P0[k], e2 = P2[k] - P0[k], q = p[k] - P0[k];
		d00 += e1 * e1; d01 += e1 * e2; d11 += e2 * e2; d20 += q * e1; d21 += q * e2;
	}
	den = d00 * d11 - d01 * d01;
	u = den > 0 ? (d11 * d20 - d01 * d21) / den : 1.0 / 3;
	v = den > 0 ? (d00 * d21 - d01 * d20) / den : 1.0 / 3;
	if (u < 0) u = 0;
	if (v < 0) v = 0;
	if (u + v > 1)
	{
		double s = u + v;
		u /= s;
		v /= s;
	}
	*a = u;
	*b = v;
}

/* the reference mesh's attributes onto m by its footpoints and its faces' centre locations fl */
qaws_status mx_transfer(qaws_hm* m, qaws_mesh_reference const* r, qaws_mesh_location const* fl)
{
	qaws_hm const* src = &r->mesh->hm;
	unsigned int f, k, v, c;
	for (k = 0; k < src->na; k++)
	{
		hm_attr const* at = &src->attr[k];
		unsigned int a = hm_attr_add(m, at->name, at->element, at->rule, at->n);
		qaws_scalar vals[16];
		if (a == HM_NONE)
		{
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		if (at->element == QAWS_MESH_VERTEX)
			for (v = 0; v < m->nv; v++)
			{
				unsigned int q;
				qaws_mesh_reference_evaluate_attribute(r, k, &m->foot[v], 1, vals);
				for (q = 0; q < at->n; q++)
					hm_attr_at(m, a, v)[q] = (double)vals[q];
			}
		else
			for (f = 0; f < m->nf; f++)
			{
				unsigned int q;
				if (fl[f].patch == HM_NONE)
					continue;
				if (at->element == QAWS_MESH_FACE)
				{
					qaws_mesh_reference_evaluate_attribute(r, k, &fl[f], 1, vals);
					for (q = 0; q < at->n; q++)
						hm_attr_at(m, a, f)[q] = (double)vals[q];
					continue;
				}
				/* corners: the vertex in the reference triangle under the face's centre */
				for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
				{
					qaws_mesh_location l;
					double ba, bb;
					rm_bary(r, fl[f].patch, rm_x(m, m->cv[c]), &ba, &bb);
					l.patch = fl[f].patch;
					l.a = (qaws_scalar)ba;
					l.b = (qaws_scalar)bb;
					qaws_mesh_reference_evaluate_attribute(r, k, &l, 1, vals);
					for (q = 0; q < at->n; q++)
						hm_attr_at(m, a, c)[q] = (double)vals[q];
				}
			}
	}
	return QAWS_STATUS_OK;
}

/* the reference mesh's attributes onto the result; the faces' deviations */
static qaws_status rm_end(rmsh* w)
{
	qaws_hm* m = w->m;
	qaws_mesh_location* fl;
	unsigned int f, k, c;
	double dmax = 0, dsum = 0;
	hm_compact(m);
	fl = (qaws_mesh_location*)malloc(sizeof(qaws_mesh_location) * (m->nf + 1));
	if (!fl)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* each face's centre on the reference, from its first vertex's side */
	for (f = 0; f < m->nf; f++)
	{
		double cen[3] = { 0, 0, 0 }, d;
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			for (k = 0; k < 3; k++)
				cen[k] += rm_x(m, m->cv[c])[k] / m->fn[f];
		{
			unsigned int reg = rm_region(w, f), h = m->foot[m->cv[m->fc[f]]].patch;
			d = sqrt(mr_project_region(w->r, cen, h < w->r->nt && w->r->treg[h] == reg ? h : HM_NONE, reg, &w->s, &fl[f]));
		}
		dmax = d > dmax ? d : dmax;
		dsum += d;
	}
	hm_attr_remove(m, w->ra > w->sa ? w->ra : w->sa);
	hm_attr_remove(m, w->ra > w->sa ? w->sa : w->ra);
	w->rep.max_deviation = (qaws_scalar)dmax;
	w->rep.mean_deviation = (qaws_scalar)(m->nf ? dsum / m->nf : 0);
	{
		qaws_status st = mx_transfer(m, w->r, fl);
		free(fl);
		return st;
	}
}

qaws_status qaws_mesh_remesh(qaws_mesh_3d const* mesh, qaws_remesh_desc const* desc, qaws_mesh_3d** out_mesh, qaws_mesh_report* out_report)
{
	rmsh w;
	qaws_mesh_reference* own = NULL;
	qaws_mesh_3d* out = NULL;
	qaws_status st;
	unsigned int it, iterations;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh || !desc || (unsigned int)desc->method > QAWS_REMESH_FIELD_ALIGNED || desc->target_length < 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->method == QAWS_REMESH_FIELD_ALIGNED)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(&w, 0, sizeof w);
	w.r = desc->reference;
	if (!w.r)
	{
		st = qaws_mesh_reference_create(mesh, NULL, &own);
		if (st != QAWS_STATUS_OK)
			return st;
		w.r = own;
	}
	if (desc->size && qaws_size_field_get_reference(desc->size) && qaws_size_field_get_reference(desc->size) != w.r)
	{
		qaws_mesh_reference_destroy(own);
		return QAWS_STATUS_INVALID_ARGUMENT;   /* (a field sampled on another reference) */
	}
	if (!w.r->nt)
	{
		qaws_mesh_reference_destroy(own);
		return qaws_mesh_3d_clone(mesh, out_mesh);
	}
	w.size = desc->size;
	w.L = desc->target_length > 0 ? (double)desc->target_length : 0;
	if (w.L == 0)
	{
		qaws_mesh_quality q;
		qaws_mesh_3d_measure_quality(mesh, &q);
		w.L = q.mean_edge > 0 ? (double)q.mean_edge : 1;
	}
	iterations = desc->method == QAWS_REMESH_CVT ? 5 : desc->iterations ? desc->iterations : 10;
	st = rm_begin(&w, mesh, &out);
	for (it = 0; st == QAWS_STATUS_OK && it < iterations && !w.failed; it++)
	{
		rm_split(&w);
		hm_compact(w.m);
		rm_collapse(&w);
		hm_compact(w.m);
		rm_flip(&w);
		hm_compact(w.m);
		rm_smooth(&w);
		rm_caps(&w);
		hm_compact(w.m);
		w.rep.iterations++;
	}
	if (st == QAWS_STATUS_OK && !w.failed && desc->method == QAWS_REMESH_CVT && !rv_optimize(w.m, w.r, w.size, w.L, w.sa, w.ra, &w.s, desc->iterations ? desc->iterations : 30, &w.rep))
		w.failed = 1;
	if (w.failed && st == QAWS_STATUS_OK)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (st == QAWS_STATUS_OK)
		st = rm_end(&w);
	mr_scratch_free(&w.s);
	free(w.ep);
	qaws_mesh_reference_destroy(own);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(out);
		return st;
	}
	qaws_mesh_3d_measure_quality(out, &w.rep.quality);
	if (out_report)
		*out_report = w.rep;
	*out_mesh = out;
	return QAWS_STATUS_OK;
}
