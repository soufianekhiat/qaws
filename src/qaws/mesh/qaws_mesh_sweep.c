/*
 * Swept solids (qaws_mesh_sweep.h): time cut into slabs by the chords'
 * deviation, each slab's front faces extruded, all united by the Boolean.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_sweep.h"
#include "../qaws_mesh_decimate.h"
#include "../qaws_clip64.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct swp_ctx
{
	qaws_motion const* m;
	unsigned int nv;
	qaws_scalar const* body;   /* a rigid motion's points (the mesh's positions) */
	qaws_scalar* buf;          /* scratch: nv points */
} swp_ctx;

/* the vertices at t, in double */
static qaws_status swp_at(swp_ctx* c, double t, double* X)
{
	unsigned int i;
	qaws_status st = qaws_motion_transform_points(c->buf, NULL, c->m, (qaws_scalar)t, qaws_motion_is_per_vertex(c->m) ? NULL : c->body, c->nv);
	for (i = 0; st == QAWS_STATUS_OK && i < 3 * c->nv; i++) X[i] = (double)c->buf[i];
	return st;
}

/* the soup as a mesh carrying the swept mesh's attributes: per soup face its triangle (sof; bit 31: turned, 30: a wall), per
   soup corner the triangle's corner it repeats (soc), on the layer of its vertex; vertex attributes become corner ones,
   directions (NORMAL, and NORMALIZED of 3 components) turn with a rigid motion, NORMAL ones are negated on turned faces;
   a wall's NORMAL face value is its own normal (its face's is not: a wall stands across it); a soup edge repeating an
   edge of the mesh at one layer (a cap's, a wall's top or bottom) takes its edge values and feature flag, the others
   zeros */
/* o = the rotation R's quaternion times q (w x y z) */
static void swp_compose(qaws_scalar* o, qaws_scalar const* q, double const* R)
{
	double tr = R[0] + R[4] + R[8], r[4], s;
	if (tr > 0)
	{
		s = 2 * sqrt(tr + 1);
		r[0] = s / 4; r[1] = (R[7] - R[5]) / s; r[2] = (R[2] - R[6]) / s; r[3] = (R[3] - R[1]) / s;
	}
	else if (R[0] > R[4] && R[0] > R[8])
	{
		s = 2 * sqrt(1 + R[0] - R[4] - R[8]);
		r[0] = (R[7] - R[5]) / s; r[1] = s / 4; r[2] = (R[1] + R[3]) / s; r[3] = (R[2] + R[6]) / s;
	}
	else if (R[4] > R[8])
	{
		s = 2 * sqrt(1 + R[4] - R[0] - R[8]);
		r[0] = (R[2] - R[6]) / s; r[1] = (R[1] + R[3]) / s; r[2] = s / 4; r[3] = (R[5] + R[7]) / s;
	}
	else
	{
		s = 2 * sqrt(1 + R[8] - R[0] - R[4]);
		r[0] = (R[3] - R[1]) / s; r[1] = (R[2] + R[6]) / s; r[2] = (R[5] + R[7]) / s; r[3] = s / 4;
	}
	{
		double w = (double)q[0], x = (double)q[1], y = (double)q[2], z = (double)q[3];
		o[0] = (qaws_scalar)(r[0] * w - r[1] * x - r[2] * y - r[3] * z);
		o[1] = (qaws_scalar)(r[0] * x + r[1] * w + r[2] * z - r[3] * y);
		o[2] = (qaws_scalar)(r[0] * y - r[1] * z + r[2] * w + r[3] * x);
		o[3] = (qaws_scalar)(r[0] * z + r[1] * y - r[2] * x + r[3] * w);
	}
}

/* a rail's copy in the soup: its given vertices' key, its soup edge, its wall's normal */
typedef struct swp_rail
{
	unsigned long long key;
	unsigned int a, b;
	double n[3];
} swp_rail;

static int swp_rail_cmp(void const* p, void const* q)
{
	unsigned long long x = ((swp_rail const*)p)->key, y = ((swp_rail const*)q)->key;
	return x < y ? -1 : x > y;
}

/* a triangle edge's record: its vertices (smaller first), its triangle and corner */
typedef struct swp_erec
{
	unsigned int a, b, t, c;
} swp_erec;

static int swp_erec_cmp(void const* p, void const* q)
{
	swp_erec const *x = (swp_erec const*)p, *y = (swp_erec const*)q;
	if (x->a != y->a)
		return x->a < y->a ? -1 : 1;
	if (x->b != y->b)
		return x->b < y->b ? -1 : 1;
	return x->t < y->t ? -1 : x->t > y->t;
}

/* the triangles t (nt) as flat faces: each paired with a neighbour across an edge in its plane (positions x, the plane
   within eps) into a quad, else alone; per face its vertices gv (4 each, in turn), their corners gc, its size gs; per
   triangle its face tg. The face count (0: out of memory) */
static unsigned int swp_flat_faces(unsigned int* gv, unsigned int* gc, unsigned int* gs, unsigned int* tg, unsigned int const* t, unsigned int nt,
	qaws_scalar const* x, double eps)
{
	swp_erec* e = (swp_erec*)qaws_internal_scratch_malloc(sizeof(swp_erec) * (3 * (size_t)nt + 1));
	unsigned int i, k, ng = 0;
	if (!e)
		return 0;
	for (i = 0; i < 3 * nt; i++)
	{
		unsigned int a = t[i], b = t[3 * (i / 3) + (i % 3 + 1) % 3];
		e[i].a = a < b ? a : b;
		e[i].b = a < b ? b : a;
		e[i].t = i / 3;
		e[i].c = i;
	}
	qsort(e, 3 * (size_t)nt, sizeof(swp_erec), swp_erec_cmp);
	for (i = 0; i < nt; i++)
		tg[i] = 0xFFFFFFFFu;
	for (i = 0; i < 3 * nt; i++)
	{
		unsigned int ti, tj, ci, cj, s;
		double u[3], w[3], n1[3], n2[3], l1, l2, d;
		if (i + 1 >= 3 * nt || e[i].a != e[i + 1].a || e[i].b != e[i + 1].b || (i + 2 < 3 * nt && e[i + 2].a == e[i].a && e[i + 2].b == e[i].b))
			continue;
		ti = e[i].t;
		tj = e[i + 1].t;
		ci = e[i].c;
		cj = e[i + 1].c;
		if (ti == tj || tg[ti] != 0xFFFFFFFFu || tg[tj] != 0xFFFFFFFFu || t[ci] == t[cj])
			continue;   /* (the same edge run the same way: not a pair) */
		for (k = 0; k < 3; k++)
		{
			u[k] = (double)x[3 * (size_t)t[3 * ti + 1] + k] - (double)x[3 * (size_t)t[3 * ti] + k];
			w[k] = (double)x[3 * (size_t)t[3 * ti + 2] + k] - (double)x[3 * (size_t)t[3 * ti] + k];
		}
		n1[0] = u[1] * w[2] - u[2] * w[1]; n1[1] = u[2] * w[0] - u[0] * w[2]; n1[2] = u[0] * w[1] - u[1] * w[0];
		for (k = 0; k < 3; k++)
		{
			u[k] = (double)x[3 * (size_t)t[3 * tj + 1] + k] - (double)x[3 * (size_t)t[3 * tj] + k];
			w[k] = (double)x[3 * (size_t)t[3 * tj + 2] + k] - (double)x[3 * (size_t)t[3 * tj] + k];
		}
		n2[0] = u[1] * w[2] - u[2] * w[1]; n2[1] = u[2] * w[0] - u[0] * w[2]; n2[2] = u[0] * w[1] - u[1] * w[0];
		l1 = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
		l2 = sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
		if (!(l1 > 0) || !(l2 > 0) || (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]) / (l1 * l2) < 1 - 1e-12)
			continue;
		s = t[3 * tj + (cj % 3 + 2) % 3];   /* (tj's third vertex, off the shared edge) */
		d = 0;
		for (k = 0; k < 3; k++)
			d += ((double)x[3 * (size_t)s + k] - (double)x[3 * (size_t)t[3 * ti] + k]) * n1[k] / l1;
		if (fabs(d) > eps)
			continue;
		/* ti runs a -> b (corner ci), tj b -> a then s: the quad a, s, b, ti's third */
		gv[4 * ng] = t[ci]; gc[4 * ng] = ci;
		gv[4 * ng + 1] = s; gc[4 * ng + 1] = 3 * tj + (cj % 3 + 2) % 3;
		gv[4 * ng + 2] = t[3 * ti + (ci % 3 + 1) % 3]; gc[4 * ng + 2] = 3 * ti + (ci % 3 + 1) % 3;
		gv[4 * ng + 3] = t[3 * ti + (ci % 3 + 2) % 3]; gc[4 * ng + 3] = 3 * ti + (ci % 3 + 2) % 3;
		gs[ng] = 4;
		tg[ti] = tg[tj] = ng++;
	}
	for (i = 0; i < nt; i++)
		if (tg[i] == 0xFFFFFFFFu)
		{
			for (k = 0; k < 3; k++)
			{
				gv[4 * ng + k] = t[3 * i + k];
				gc[4 * ng + k] = 3 * i + k;
			}
			gs[ng] = 3;
			tg[i] = ng++;
		}
	qaws_internal_free(e);
	return ng;
}

/* the convex hull of a flat face's two positions (its gn vertices tv, 3 or 4, at layers oa and ob of positions x; tc their
   corners), when a corner moves out of
   the solid by more than out_min (less is within the tolerance: a face sliding along itself): its triangles (3 soup vertices each) into out, their face flags (f; bit 31 a cap turned against the face, bit 30
   a side) into of, their corners (the face's) into oc; their count, at most 12 (0: none, or flat). caps: the caps it has (1 the
   start's, 2 the end's; 4 and 8 when turned), those in skip left out. Coplanar hull points (a translated face's sides) make one polygon, fanned
   once; eps: the plane distance below which points are on it */
static unsigned int swp_hull(unsigned int* out, unsigned int* of, unsigned int* oc, unsigned int* caps, qaws_scalar const* x, unsigned int const* tv,
	unsigned int const* tc, unsigned int gn, unsigned int f, unsigned int oa, unsigned int ob, double eps, double out_min, unsigned int skip)
{
	double P[8][3], fn[3], u[3], w[3], fl, smax = 0;
	unsigned int id[8], cor[8], np = 0, i, j, k, m, n = 0, masks[56], nm = 0;
	int side[56];
	*caps = 0;
	for (i = 0; i < 2 * gn; i++)
	{
		unsigned int v = (i < gn ? oa : ob) + tv[i % gn], dup = 0;
		double q[3];
		for (k = 0; k < 3; k++)
			q[k] = (double)x[3 * (size_t)v + k];
		for (j = 0; j < np && !dup; j++)
			dup = fabs(P[j][0] - q[0]) + fabs(P[j][1] - q[1]) + fabs(P[j][2] - q[2]) <= eps;   /* (a corner that does not move) */
		if (dup)
			continue;
		memcpy(P[np], q, sizeof q);
		id[np] = v;
		cor[np++] = tc[i % gn];
	}
	/* a corner moving out: its displacement against the face's normal at the start */
	for (k = 0; k < 3; k++)
	{
		u[k] = (double)x[3 * ((size_t)oa + tv[1]) + k] - (double)x[3 * ((size_t)oa + tv[0]) + k];
		w[k] = (double)x[3 * ((size_t)oa + tv[2]) + k] - (double)x[3 * ((size_t)oa + tv[0]) + k];
	}
	fn[0] = u[1] * w[2] - u[2] * w[1]; fn[1] = u[2] * w[0] - u[0] * w[2]; fn[2] = u[0] * w[1] - u[1] * w[0];
	fl = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);
	if (!(fl > 0) || np < 4)
		return 0;
	for (i = 0; i < gn; i++)
	{
		double s = 0;
		for (k = 0; k < 3; k++)
			s += ((double)x[3 * ((size_t)ob + tv[i]) + k] - (double)x[3 * ((size_t)oa + tv[i]) + k]) * fn[k] / fl;
		smax = fmax(smax, s);
	}
	if (!(smax > out_min))
		return 0;
	/* the supporting planes: every triple with all the points on one side, by the points on it (a mask) */
	for (i = 0; i < np; i++)
		for (j = i + 1; j < np; j++)
			for (k = j + 1; k < np; k++)
			{
				double a[3], b[3], nn[3], l, lo = 0, hi = 0;
				unsigned int mask = 0, q;
				int sd;
				for (m = 0; m < 3; m++)
				{
					a[m] = P[j][m] - P[i][m];
					b[m] = P[k][m] - P[i][m];
				}
				nn[0] = a[1] * b[2] - a[2] * b[1]; nn[1] = a[2] * b[0] - a[0] * b[2]; nn[2] = a[0] * b[1] - a[1] * b[0];
				l = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
				if (!(l > 0))
					continue;
				for (m = 0; m < np; m++)
				{
					double d = ((P[m][0] - P[i][0]) * nn[0] + (P[m][1] - P[i][1]) * nn[1] + (P[m][2] - P[i][2]) * nn[2]) / l;
					if (fabs(d) <= eps)
						mask |= 1u << m;
					else
					{
						lo = fmin(lo, d);
						hi = fmax(hi, d);
					}
				}
				if (lo < 0 && hi > 0)
					continue;
				if (lo == 0 && hi == 0)
					return 0;   /* (all on one plane: flat) */
				sd = hi > 0 ? -1 : 1;   /* (the others behind: outward along nn; in front: against it) */
				for (q = 0; q < nm && masks[q] != mask; q++) {}
				if (q < nm)
					continue;
				masks[nm] = mask;
				side[nm++] = sd;
				/* the plane's points in turn round their centroid, seen from outside, fanned */
				{
					unsigned int pts[8], c = 0, r, s2;
					double cen[3] = { 0, 0, 0 }, e1[3], e2[3], on[3], ang[8], el;
					for (m = 0; m < np; m++)
						if (mask >> m & 1)
							pts[c++] = m;
					for (m = 0; m < c; m++)
						for (r = 0; r < 3; r++)
							cen[r] += P[pts[m]][r] / c;
					for (r = 0; r < 3; r++)
						on[r] = sd * nn[r] / l;
					for (r = 0; r < 3; r++)
						e1[r] = P[pts[0]][r] - cen[r];
					el = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
					for (r = 0; r < 3; r++)
						e1[r] /= el;
					e2[0] = on[1] * e1[2] - on[2] * e1[1]; e2[1] = on[2] * e1[0] - on[0] * e1[2]; e2[2] = on[0] * e1[1] - on[1] * e1[0];
					for (m = 0; m < c; m++)
					{
						double dv[3];
						for (r = 0; r < 3; r++)
							dv[r] = P[pts[m]][r] - cen[r];
						ang[m] = atan2(dv[0] * e2[0] + dv[1] * e2[1] + dv[2] * e2[2], dv[0] * e1[0] + dv[1] * e1[1] + dv[2] * e1[2]);
					}
					for (m = 1; m < c; m++)
						for (r = m; r > 0 && ang[r - 1] > ang[r]; r--)
						{
							double ta = ang[r];
							unsigned int tp = pts[r];
							ang[r] = ang[r - 1]; ang[r - 1] = ta;
							pts[r] = pts[r - 1]; pts[r - 1] = tp;
						}
					/* (a point inside the plane's polygon, or on one of its edges, dropped: the fan covers the polygon once) */
					for (r = 0; c > 3 && r < c; )
					{
						unsigned int pm = pts[(r + c - 1) % c], pc = pts[r], pn = pts[(r + 1) % c], q2;
						double a1[3], a2[3], cr;
						for (q2 = 0; q2 < 3; q2++)
						{
							a1[q2] = P[pc][q2] - P[pm][q2];
							a2[q2] = P[pn][q2] - P[pc][q2];
						}
						cr = on[0] * (a1[1] * a2[2] - a1[2] * a2[1]) + on[1] * (a1[2] * a2[0] - a1[0] * a2[2]) + on[2] * (a1[0] * a2[1] - a1[1] * a2[0]);
						if (cr > eps * (sqrt(a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2]) + sqrt(a2[0] * a2[0] + a2[1] * a2[1] + a2[2] * a2[2])))
						{
							r++;
							continue;
						}
						for (q2 = r; q2 + 1 < c; q2++)
							pts[q2] = pts[q2 + 1];
						c--;
						r = 0;
					}
					/* (fanned from its lowest soup vertex: two hulls sharing the polygon fan it alike, their triangles cancel exactly) */
					{
						unsigned int low = 0, rt[8];
						for (r = 1; r < c; r++)
							if (id[pts[r]] < id[pts[low]])
								low = r;
						for (r = 0; r < c; r++)
							rt[r] = pts[(low + r) % c];
						memcpy(pts, rt, sizeof(unsigned int) * c);
					}
					for (s2 = 1; s2 + 1 < c && n < 12; s2++)
					{
						unsigned int tri[3], lay = 0;
						double tn;
						tri[0] = pts[0]; tri[1] = pts[s2]; tri[2] = pts[s2 + 1];
						for (r = 0; r < 3; r++)
						{
							out[3 * n + r] = id[tri[r]];
							oc[3 * n + r] = cor[tri[r]];
							lay |= id[tri[r]] >= ob ? 2u : 1u;
						}
						/* a cap (all three at one layer, bit 1 the start's, 2 the end's; 4 and 8 turned): turned when it faces against the face;
						   else a side */
						tn = on[0] * fn[0] + on[1] * fn[1] + on[2] * fn[2];
						if (lay != 3)
						{
							*caps |= lay | (tn < 0 ? lay << 2 : 0u);
							if (skip & lay)
								continue;
						}
						of[n] = lay == 3 ? f | 0x40000000u : tn < 0 ? f | 0x80000000u : f;
						n++;
					}
				}
			}
	return n;
}

/* a soup triangle's key: its vertices from the lowest, in turn */
typedef struct swp_tkey
{
	unsigned int v[3], i;
} swp_tkey;

static int swp_tkey_cmp(void const* p, void const* q)
{
	swp_tkey const *a = (swp_tkey const*)p, *b = (swp_tkey const*)q;
	unsigned int k;
	for (k = 0; k < 3; k++)
		if (a->v[k] != b->v[k])
			return a->v[k] < b->v[k] ? -1 : 1;
	return a->i < b->i ? -1 : a->i > b->i;
}

static void swp_tkey_set(swp_tkey* k, unsigned int a, unsigned int b, unsigned int c)
{
	if (a < b && a < c) { k->v[0] = a; k->v[1] = b; k->v[2] = c; }
	else if (b < c) { k->v[0] = b; k->v[1] = c; k->v[2] = a; }
	else { k->v[0] = c; k->v[1] = a; k->v[2] = b; }
}

/* the soup's triangles that meet one turned the other way (the same three vertices) left out in pairs: they cancel, and
   the union resolves coincident faces less well than none */
static qaws_status swp_cancel(unsigned int* sf, unsigned int* sof, unsigned int* soc, unsigned int* snf)
{
	unsigned int n = *snf, i, j, m = 0;
	swp_tkey* key = (swp_tkey*)qaws_internal_scratch_malloc(sizeof(swp_tkey) * ((size_t)n + 1));
	unsigned char* gone = (unsigned char*)qaws_internal_scratch_zalloc((size_t)n + 1, 1);
	if (!key || !gone)
	{
		qaws_internal_free(key);
		qaws_internal_free(gone);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < n; i++)
	{
		swp_tkey_set(&key[i], sf[3 * (size_t)i], sf[3 * (size_t)i + 1], sf[3 * (size_t)i + 2]);
		key[i].i = i;
	}
	qsort(key, n, sizeof(swp_tkey), swp_tkey_cmp);
	for (i = 0; i < n; i++)
	{
		swp_tkey want;
		unsigned int lo = 0, hi = n, t0 = key[i].i;
		if (gone[t0])
			continue;
		swp_tkey_set(&want, sf[3 * (size_t)t0], sf[3 * (size_t)t0 + 2], sf[3 * (size_t)t0 + 1]);
		want.i = 0;
		while (lo < hi)
		{
			unsigned int mid = lo + (hi - lo) / 2;
			if (swp_tkey_cmp(&key[mid], &want) < 0) lo = mid + 1; else hi = mid;
		}
		for (j = lo; j < n && !memcmp(key[j].v, want.v, sizeof want.v); j++)
			if (!gone[key[j].i] && key[j].i != t0)
			{
				gone[key[j].i] = gone[t0] = 1;
				break;
			}
	}
	for (i = 0; i < n; i++)
		if (!gone[i])
		{
			memmove(sf + 3 * (size_t)m, sf + 3 * (size_t)i, sizeof(unsigned int) * 3);
			memmove(soc + 3 * (size_t)m, soc + 3 * (size_t)i, sizeof(unsigned int) * 3);
			sof[m++] = sof[i];
		}
	*snf = m;
	qaws_internal_free(key);
	qaws_internal_free(gone);
	return QAWS_STATUS_OK;
}

/* soup vertex v's given index (creation splits the soup's non-manifold vertices: copies past its count) */
static unsigned int swp_src(qaws_mesh_3d const* m, unsigned int v, unsigned int n)
{
	unsigned int s = qaws_mesh_3d_get_vertex_source(m, v);
	return s < n ? s : v;
}

/* each edge's two first faces' unit normals (6 per edge) and how many were seen (rc), by m's triangles; NULL out of memory */
static double* swp_edge_normals(qaws_mesh_3d const* m, unsigned int const* e, unsigned int ne, unsigned long long const* key, unsigned char* rc)
{
	unsigned int nt = qaws_mesh_3d_get_triangle_count(m), i, j, k;
	unsigned int* t = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	double* rn = (double*)qaws_internal_scratch_zalloc(6 * (size_t)ne + 1, sizeof(double));
	(void)e;
	if (!t || !rn)
	{
		qaws_internal_free(t);
		qaws_internal_free(rn);
		return NULL;
	}
	qaws_mesh_3d_get_triangles(t, m);
	for (i = 0; i < nt; i++)
	{
		qaws_vec3 a = qaws_mesh_3d_get_vertex(m, t[3 * (size_t)i]), b = qaws_mesh_3d_get_vertex(m, t[3 * (size_t)i + 1]), c = qaws_mesh_3d_get_vertex(m, t[3 * (size_t)i + 2]);
		double u[3], w[3], n[3], l;
		u[0] = (double)b.x - (double)a.x; u[1] = (double)b.y - (double)a.y; u[2] = (double)b.z - (double)a.z;
		w[0] = (double)c.x - (double)a.x; w[1] = (double)c.y - (double)a.y; w[2] = (double)c.z - (double)a.z;
		n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		for (j = 0; j < 3 && l > 0; j++)
		{
			unsigned int x = cy_find(key, ne, t[3 * (size_t)i + j], t[3 * (size_t)i + (j + 1) % 3]);
			if (x == QAWS_MESH_NONE || rc[x] >= 2)
				continue;
			for (k = 0; k < 3; k++)
				rn[6 * (size_t)x + 3 * rc[x] + k] = n[k] / l;
			rc[x]++;
		}
	}
	qaws_internal_free(t);
	return rn;
}

/* the smallest angle between the faces at a feature edge of m (radians), pi / 6 without features; 0 out of memory */
static double swp_feature_angle(qaws_mesh_3d const* m)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), nfe = qaws_mesh_3d_get_feature_edge_count(m), i;
	unsigned int *e = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)(ne > nfe ? ne : nfe) + 1));
	unsigned char* rc = (unsigned char*)qaws_internal_scratch_zalloc((size_t)ne + 1, 1);
	unsigned long long* key = NULL;
	double* rn = NULL, best = HUGE_VAL;
	if (e && rc)
	{
		qaws_mesh_3d_get_edges(e, m);
		key = cy_keys(e, ne);
	}
	if (key && (rn = swp_edge_normals(m, e, ne, key, rc)) != NULL)
	{
		qaws_mesh_3d_get_feature_edges(e, m);
		for (i = 0; i < nfe; i++)
		{
			unsigned int x = cy_find(key, ne, e[2 * (size_t)i], e[2 * (size_t)i + 1]);
			double const* r = rn + 6 * (size_t)(x == QAWS_MESH_NONE ? 0 : x);
			double d;
			if (x == QAWS_MESH_NONE || rc[x] != 2)
				continue;
			d = r[0] * r[3] + r[1] * r[4] + r[2] * r[5];
			best = fmin(best, acos(d > 1 ? 1 : d < -1 ? -1 : d));
		}
	}
	qaws_internal_free(e); qaws_internal_free(rc); qaws_internal_free(key); qaws_internal_free(rn);
	return !rn ? 0 : best == HUGE_VAL ? 0.52359877559829887 : best;
}

/* q within 1e-9 of segment ab (relative to its length) */
static int swp_on_segment(qaws_scalar const* a, qaws_scalar const* b, qaws_vec3 q)
{
	double u[3], w[3], uu = 0, uw = 0, t, d2 = 0;
	unsigned int k;
	u[0] = (double)b[0] - (double)a[0]; u[1] = (double)b[1] - (double)a[1]; u[2] = (double)b[2] - (double)a[2];
	w[0] = (double)q.x - (double)a[0]; w[1] = (double)q.y - (double)a[1]; w[2] = (double)q.z - (double)a[2];
	for (k = 0; k < 3; k++)
	{
		uu += u[k] * u[k];
		uw += u[k] * w[k];
	}
	if (!(uu > 0))
		return 0;
	t = uw / uu;
	if (t < -1e-9 || t > 1 + 1e-9)
		return 0;
	for (k = 0; k < 3; k++)
		d2 += (w[k] - t * u[k]) * (w[k] - t * u[k]);
	return d2 <= 1e-18 * uu;
}

/* the swept solid's features: those turning by more than angle (radians) there (a copy of an input crease left inside a
   face of the solid, or on a wall bent only by the slabs' chords, dropped); those on the mesh's own features at the first
   and last poses (positions pos, layers 0 and last of nv each) need only turn at all */
static qaws_status swp_crease_features(qaws_mesh_3d* m, double angle, qaws_mesh_3d const* mesh, qaws_scalar const* pos, unsigned int nv, unsigned int last)
{
	unsigned int nin = qaws_mesh_3d_get_feature_edge_count(mesh), *ife = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nin + 1));
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), nfe = qaws_mesh_3d_get_feature_edge_count(m), i, kept = 0;
	double cmax = cos(angle);
	unsigned int *e = NULL, *fe = NULL;
	unsigned char* rc = NULL;
	double* rn = NULL;
	unsigned long long* key = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!nfe)
		return QAWS_STATUS_OK;
	e = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)ne + 1));
	fe = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)nfe + 1));
	rc = (unsigned char*)qaws_internal_scratch_zalloc((size_t)ne + 1, 1);
	if (!e || !fe || !rc || !ife)
		goto done;
	qaws_mesh_3d_get_feature_edges(ife, mesh);
	qaws_mesh_3d_get_edges(e, m);
	qaws_mesh_3d_get_feature_edges(fe, m);
	if (!(key = cy_keys(e, ne)) || !(rn = swp_edge_normals(m, e, ne, key, rc)))
		goto done;
	for (i = 0; i < nfe; i++)
	{
		unsigned int x = cy_find(key, ne, fe[2 * (size_t)i], fe[2 * (size_t)i + 1]);
		double const* r = rn + 6 * (size_t)(x == QAWS_MESH_NONE ? 0 : x);
		double lim = cmax;
		qaws_vec3 pa = qaws_mesh_3d_get_vertex(m, fe[2 * (size_t)i]), pb = qaws_mesh_3d_get_vertex(m, fe[2 * (size_t)i + 1]);
		unsigned int j, l;
		for (l = 0; l < 2 && lim == cmax; l++)
			for (j = 0; j < nin; j++)
			{
				qaws_scalar const *sa = pos + 3 * ((size_t)(l ? last : 0) * nv + ife[2 * (size_t)j]), *sb = pos + 3 * ((size_t)(l ? last : 0) * nv + ife[2 * (size_t)j + 1]);
				if (swp_on_segment(sa, sb, pa) && swp_on_segment(sa, sb, pb))
				{
					lim = 1 - 1e-9;
					break;
				}
			}
		if (x == QAWS_MESH_NONE || (rc[x] == 2 && r[0] * r[3] + r[1] * r[4] + r[2] * r[5] > lim))
			continue;
		fe[2 * (size_t)kept] = fe[2 * (size_t)i];
		fe[2 * (size_t)kept + 1] = fe[2 * (size_t)i + 1];
		kept++;
	}
	st = kept == nfe ? QAWS_STATUS_OK : qaws_mesh_3d_set_feature_edges(m, fe, kept);
done:
	qaws_internal_free(e); qaws_internal_free(fe); qaws_internal_free(rc); qaws_internal_free(rn); qaws_internal_free(key); qaws_internal_free(ife);
	return st;
}

static qaws_status swp_soup(qaws_mesh_3d** out, qaws_mesh_3d const* mesh, qaws_motion const* motion, qaws_mesh_3d_desc const* md,
	unsigned int const* sof, unsigned int const* soc, unsigned int nv, double const* layer_t, unsigned int layers)
{
	unsigned int na = qaws_mesh_3d_get_attribute_count(mesh), nf = qaws_mesh_3d_get_face_count(mesh), nc = qaws_mesh_3d_get_corner_count(mesh);
	unsigned int nt = qaws_mesh_3d_get_triangle_count(mesh), a, i, j, k, *fv = NULL, *fs = NULL, *tf = NULL, *tc = NULL;
	unsigned int ine = qaws_mesh_3d_get_edge_count(mesh), sne, *ie = NULL, *se = NULL, *smap = NULL, nfe = 0;
	unsigned long long* ik = NULL;
	unsigned char* ifeat = NULL;
	qaws_scalar *val = NULL, *sv = NULL;
	double* rot = NULL;
	qaws_status st;
	int rigid = !qaws_motion_is_per_vertex(motion);
	if ((st = qaws_mesh_3d_create(out, md)) != QAWS_STATUS_OK)
		return st;
	fv = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nc + 1));
	fs = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	tf = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
	tc = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)nt + 1));
	val = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 16, (size_t)(nv > nc ? (nv > nf ? nv : nf) : (nc > nf ? nc : nf)) + 1));
	sv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 16, 3 * (size_t)md->face_count + 1));
	rot = (double*)qaws_internal_scratch_malloc(sizeof(double) * (9 * (size_t)layers + 1));
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!fv || !fs || !tf || !tc || !val || !sv || !rot)
		goto done;
	/* the soup's edges on the mesh's (one layer, an edge of the mesh), their features */
	sne = qaws_mesh_3d_get_edge_count(*out);
	{
		unsigned int nfeat = qaws_mesh_3d_get_feature_edge_count(mesh), m = ine > nfeat ? ine : nfeat;
		ie = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)m + 1));
		se = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)sne + 1));
		smap = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)sne + 1));
		ifeat = (unsigned char*)qaws_internal_zalloc((size_t)ine + 1, 1);
		if (!ie || !se || !smap || !ifeat)
			goto done;
		qaws_mesh_3d_get_edges(ie, mesh);
		if (!(ik = cy_keys(ie, ine)))
			goto done;
		qaws_mesh_3d_get_feature_edges(ie, mesh);
		for (i = 0; i < nfeat; i++)
		{
			unsigned int x = cy_find(ik, ine, ie[2 * (size_t)i], ie[2 * (size_t)i + 1]);
			if (x != QAWS_MESH_NONE)
				ifeat[x] = 1;
		}
		qaws_mesh_3d_get_edges(se, *out);
		for (i = 0; i < sne; i++)
		{
			unsigned int x = se[2 * (size_t)i], y = se[2 * (size_t)i + 1], sx = swp_src(*out, x, md->vertex_count), sy = swp_src(*out, y, md->vertex_count);
			smap[i] = sx / nv == sy / nv ? cy_find(ik, ine, sx % nv, sy % nv) : QAWS_MESH_NONE;
			if (smap[i] != QAWS_MESH_NONE && ifeat[smap[i]])
			{
				se[2 * (size_t)nfe] = x;
				se[2 * (size_t)nfe + 1] = y;
				nfe++;
			}
		}
		/* the rails a feature corner runs along (its positions at two layers: a swept box's long edges) and a feature edge's
		   chords across two layers (where it slides along itself: its envelope), unless two hull faces meet flat there; each
		   hull holds its own copy of one (the soup's vertices split there), so the copies are grouped by their given vertices;
		   the result keeps those that are creases of it */
		{
			swp_rail* rl = (swp_rail*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(swp_rail), 3, (size_t)md->face_count + 1));
			unsigned char* corner = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nv + 1, 1);
			unsigned int nr = 0, r0;
			if (!rl || !corner)
			{
				qaws_internal_free(rl); qaws_internal_free(corner);
				goto done;
			}
			for (i = 0; i < nv; i++)
				corner[i] = (unsigned char)qaws_mesh_3d_is_feature_corner(mesh, i);
			for (i = 0; i < md->face_count; i++)
			{
				unsigned int t[3], g[3];
				double u[3], w[3], n[3], l;
				qaws_mesh_3d_get_triangle(t, *out, i);
				for (j = 0; j < 3; j++)
					g[j] = swp_src(*out, t[j], md->vertex_count);
				for (k = 0; k < 3; k++)
				{
					u[k] = (double)md->positions[3 * (size_t)g[1] + k] - (double)md->positions[3 * (size_t)g[0] + k];
					w[k] = (double)md->positions[3 * (size_t)g[2] + k] - (double)md->positions[3 * (size_t)g[0] + k];
				}
				n[0] = u[1] * w[2] - u[2] * w[1]; n[1] = u[2] * w[0] - u[0] * w[2]; n[2] = u[0] * w[1] - u[1] * w[0];
				l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
				for (j = 0; j < 3 && l > 0; j++)
				{
					unsigned int ga = g[j], gb = g[(j + 1) % 3];
					if (ga / nv == gb / nv)
						continue;
					if (ga % nv == gb % nv ? !corner[ga % nv] : cy_find(ik, ine, ga % nv, gb % nv) == QAWS_MESH_NONE || !ifeat[cy_find(ik, ine, ga % nv, gb % nv)])
						continue;   /* (a corner's rail, or a feature edge's chord across two layers: where it slides along itself, its envelope) */
					rl[nr].key = ((unsigned long long)(ga < gb ? ga : gb) << 32) | (ga < gb ? gb : ga);
					rl[nr].a = t[j];
					rl[nr].b = t[(j + 1) % 3];
					for (k = 0; k < 3; k++)
						rl[nr].n[k] = n[k] / l;
					nr++;
				}
			}
			qsort(rl, nr, sizeof(swp_rail), swp_rail_cmp);
			for (r0 = 0; r0 < nr; )
			{
				unsigned int r1 = r0 + 1, q;
				int sharp;
				while (r1 < nr && rl[r1].key == rl[r0].key)
					r1++;
				sharp = r1 - r0 == 1;   /* (one wall: the rail lies between it and a face of the solid, the result tells) */
				for (q = r0 + 1; q < r1; q++)
					sharp |= fabs(rl[r0].n[0] * rl[q].n[0] + rl[r0].n[1] * rl[q].n[1] + rl[r0].n[2] * rl[q].n[2]) < 1 - 1e-9;
				for (q = r0; sharp && q < r1; q++)
				{
					unsigned int q2, seen = 0;
					for (q2 = r0; q2 < q; q2++)
						seen |= (rl[q2].a == rl[q].a && rl[q2].b == rl[q].b) || (rl[q2].a == rl[q].b && rl[q2].b == rl[q].a);
					if (seen || nfe >= sne)
						continue;
					/* (every copy: the union reads each result edge against its own source face's) */
					se[2 * (size_t)nfe] = rl[q].a;
					se[2 * (size_t)nfe + 1] = rl[q].b;
					nfe++;
				}
				r0 = r1;
			}
			qaws_internal_free(rl); qaws_internal_free(corner);
		}
		if (nfe && qaws_mesh_3d_set_feature_edges(*out, se, nfe) != QAWS_STATUS_OK)
			goto done;
	}
	/* the triangles (fanned from each face's first corner) back to their faces' corners */
	qaws_mesh_3d_get_faces(fv, fs, mesh);
	for (i = 0, j = 0, k = 0; i < nf; k += fs[i++])
		for (a = 0; a + 2 < fs[i]; a++, j++)
		{
			tf[j] = i;
			tc[3 * (size_t)j] = k;
			tc[3 * (size_t)j + 1] = k + a + 1;
			tc[3 * (size_t)j + 2] = k + a + 2;
		}
	for (i = 0; i < layers; i++)
	{
		qaws_motion_state s;
		if (rigid && qaws_motion_evaluate(&s, motion, (qaws_scalar)layer_t[i]) == QAWS_STATUS_OK)
			for (k = 0; k < 9; k++)
				rot[9 * (size_t)i + k] = (double)s.rotation[k];
		else
			for (k = 0; k < 9; k++)
				rot[9 * (size_t)i + k] = k % 4 == 0;
	}
	st = QAWS_STATUS_OK;
	for (a = 0; a < na && st == QAWS_STATUS_OK; a++)
	{
		qaws_attribute_desc d;
		unsigned int n, id, cnt;
		int dir;
		qaws_mesh_3d_get_attribute_desc(&d, mesh, a);
		if (d.name && !strcmp(d.name, "source"))
			continue;
		n = d.components;
		dir = d.rule == QAWS_ATTRIBUTE_NORMAL || d.rule == QAWS_ATTRIBUTE_TANGENT || d.rule == QAWS_ATTRIBUTE_VECTOR || (d.rule == QAWS_ATTRIBUTE_NORMALIZED && n == 3);
		qaws_mesh_3d_get_attribute(val, mesh, a);
		if (d.element == QAWS_MESH_EDGE)
		{
			for (i = 0; i < sne; i++)
				for (k = 0; k < n; k++)
					sv[(size_t)n * i + k] = smap[i] != QAWS_MESH_NONE ? val[(size_t)n * smap[i] + k] : 0;
			st = qaws_mesh_3d_add_attribute(&id, *out, sv, &d);
			continue;
		}
		cnt = d.element == QAWS_MESH_FACE ? md->face_count : 3 * md->face_count;
		for (i = 0; i < cnt; i++)
		{
			unsigned int g = d.element == QAWS_MESH_FACE ? sof[i] : sof[i / 3], f = g & 0x3FFFFFFFu, turned = g >> 31;
			unsigned int c = d.element == QAWS_MESH_FACE ? 0 : tc[soc[i]], e = d.element == QAWS_MESH_FACE ? tf[f] : d.element == QAWS_MESH_VERTEX ? fv[c] : c;
			unsigned int layer = d.element == QAWS_MESH_FACE ? 0 : md->face_vertices[i] / nv;
			qaws_scalar const* x = val + (size_t)n * e;
			qaws_scalar* o = sv + (size_t)n * i;
			for (k = 0; k < n; k++)
				o[k] = x[k];
			if (dir && d.element != QAWS_MESH_FACE)
			{
				double const* R = rot + 9 * (size_t)layer;
				for (k = 0; k < 3; k++)
					o[k] = (qaws_scalar)(R[3 * k] * (double)x[0] + R[3 * k + 1] * (double)x[1] + R[3 * k + 2] * (double)x[2]);
			}
			if (d.rule == QAWS_ATTRIBUTE_QUATERNION && d.element != QAWS_MESH_FACE)
				swp_compose(o, x, rot + 9 * (size_t)layer);
			if (d.rule == QAWS_ATTRIBUTE_NORMAL && turned)
				for (k = 0; k < 3; k++)
					o[k] = -o[k];
			if (d.rule == QAWS_ATTRIBUTE_TANGENT && turned)
				o[3] = -o[3];
		}
		for (i = 0; d.rule == QAWS_ATTRIBUTE_NORMAL && d.element == QAWS_MESH_FACE && i < cnt; i++)
			if (sof[i] >> 30 & 1)
			{
				qaws_scalar const *p0 = md->positions + 3 * (size_t)md->face_vertices[3 * (size_t)i], *p1 = md->positions + 3 * (size_t)md->face_vertices[3 * (size_t)i + 1],
					*p2 = md->positions + 3 * (size_t)md->face_vertices[3 * (size_t)i + 2];
				double u[3], w[3], nn[3], l;
				for (k = 0; k < 3; k++)
				{
					u[k] = (double)p1[k] - (double)p0[k];
					w[k] = (double)p2[k] - (double)p0[k];
				}
				nn[0] = u[1] * w[2] - u[2] * w[1];
				nn[1] = u[2] * w[0] - u[0] * w[2];
				nn[2] = u[0] * w[1] - u[1] * w[0];
				l = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
				for (k = 0; k < 3; k++)
					sv[3 * (size_t)i + k] = (qaws_scalar)(l > 0 ? nn[k] / l : 0);
			}
		/* (a face attribute's directions: the motion at 0, as the solid at 0's) */
		if (dir && d.element == QAWS_MESH_FACE)
			for (i = 0; i < cnt; i++)
			{
				qaws_scalar* o = sv + (size_t)n * i;
				double y[3];
				for (k = 0; k < 3; k++)
					y[k] = rot[3 * k] * (double)o[0] + rot[3 * k + 1] * (double)o[1] + rot[3 * k + 2] * (double)o[2];
				for (k = 0; k < 3; k++)
					o[k] = (qaws_scalar)y[k];
			}
		if (d.element == QAWS_MESH_VERTEX)
			d.element = QAWS_MESH_CORNER;
		st = qaws_mesh_3d_add_attribute(&id, *out, sv, &d);
	}
done:
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
	}
	qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(tf); qaws_internal_free(tc); qaws_internal_free(val);
	qaws_internal_free(sv); qaws_internal_free(rot); qaws_internal_free(ie); qaws_internal_free(se); qaws_internal_free(smap);
	qaws_internal_free(ik); qaws_internal_free(ifeat);
	return st;
}

static int swp_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : x > y;
}

/* the angle at p of the triangle p a b */
static double swp_angle(double const* p, double const* a, double const* b)
{
	double u[3], v[3], w[3];
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = a[k] - p[k];
		v[k] = b[k] - p[k];
	}
	w[0] = u[1] * v[2] - u[2] * v[1];
	w[1] = u[2] * v[0] - u[0] * v[2];
	w[2] = u[0] * v[1] - u[1] * v[0];
	return atan2(sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]), u[0] * v[0] + u[1] * v[1] + u[2] * v[2]);
}

/* edges inside flat regions (two triangles in one plane, not a feature or seam) turned toward Delaunay: the collapses
   leave fans of long thin triangles; the shape does not change. 0 on allocation failure. */
static int swp_flat_delaunay(qaws_hm* m, double eps)
{
	unsigned int pass, c, flips = 1;
	for (pass = 0; pass < 64 && flips; pass++)
	{
		flips = 0;
		for (c = 0; c < m->nc; c++)
		{
			unsigned int d = m->co[c], k;
			double const *a, *b, *pp, *q;
			double n1[3], n2[3], l1, l2, s = 0, h = 0, e1[3], e2[3];
			if (m->cf[c] == HM_NONE || d == HM_NONE || d < c || (m->cflag[c] & HM_FEATURE) || !hm_can_flip(m, c))
				continue;
			a = m->x + 3 * (size_t)m->cv[c];
			b = m->x + 3 * (size_t)hm_to(m, c);
			pp = m->x + 3 * (size_t)m->cv[hm_prev(m, c)];
			q = m->x + 3 * (size_t)m->cv[hm_prev(m, d)];
			hm_face_normal(m, m->cf[c], n1);
			hm_face_normal(m, m->cf[d], n2);
			l1 = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
			l2 = sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
			if (!(l1 > 0) || !(l2 > 0))
				continue;
			for (k = 0; k < 3; k++)
			{
				s += n1[k] * n2[k];
				h += n1[k] * (q[k] - a[k]);
			}
			if (s < (1 - 1e-12) * l1 * l2 || fabs(h) > eps * l1 || swp_angle(pp, a, b) + swp_angle(q, b, a) <= 3.14159265358979323846 + 1e-9)
				continue;
			/* the new faces (a q p) and (q b p) turned as the old */
			for (k = 0; k < 3; k++)
			{
				e1[k] = q[k] - a[k];
				e2[k] = pp[k] - a[k];
			}
			if ((e1[1] * e2[2] - e1[2] * e2[1]) * n1[0] + (e1[2] * e2[0] - e1[0] * e2[2]) * n1[1] + (e1[0] * e2[1] - e1[1] * e2[0]) * n1[2] <= 0)
				continue;
			for (k = 0; k < 3; k++)
			{
				e1[k] = b[k] - q[k];
				e2[k] = pp[k] - q[k];
			}
			if ((e1[1] * e2[2] - e1[2] * e2[1]) * n1[0] + (e1[2] * e2[0] - e1[0] * e2[2]) * n1[1] + (e1[0] * e2[1] - e1[1] * e2[0]) * n1[2] <= 0)
				continue;
			if (!hm_flip(m, c))
				return 0;
			flips++;
		}
	}
	hm_compact(m);
	return 1;
}

/* the overlapping hulls and cuts leave many vertices in flat regions and on straight creases: those whose removal moves
   nothing (an error within a thousandth of the tolerance) are collapsed, features and attributes kept (seams are lines
   the collapses keep), then the flat regions turned toward Delaunay. The Boolean's
   "source" (the soup's faces, internal) is dropped first: every cut would be a seam */
static qaws_status swp_tidy(qaws_mesh_3d* m, double tol)
{
	qaws_decimate_desc d;
	qaws_attribute_desc ad;
	unsigned int a, i, ne, nf = 0, *ev = NULL;
	qaws_scalar* fv = NULL;
	qaws_status st;
	if (!m)
		return QAWS_STATUS_OK;
	if (qaws_mesh_3d_find_attribute(&a, m, "source") == QAWS_STATUS_OK && (st = qaws_mesh_3d_remove_attribute(m, a)) != QAWS_STATUS_OK)
		return st;
	/* the features carried as an edge value (pieces keep it): the decimation makes every seam a feature, they are put back after */
	ne = qaws_mesh_3d_get_edge_count(m);
	ev = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (2 * (size_t)ne + 1));
	fv = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * ((size_t)ne + 1));
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!ev || !fv)
		goto done;
	qaws_mesh_3d_get_edges(ev, m);
	for (i = 0; i < ne; i++)
	{
		unsigned int c = hm_find_edge(&m->hm, ev[2 * (size_t)i], ev[2 * (size_t)i + 1]);
		fv[i] = (qaws_scalar)(c != HM_NONE && (m->hm.cflag[c] & HM_FEATURE) ? 1 : 0);
	}
	memset(&ad, 0, sizeof ad);
	ad.name = "qaws_sweep_feature";
	ad.element = QAWS_MESH_EDGE;
	ad.components = 1;
	ad.rule = QAWS_ATTRIBUTE_NEAREST;
	if ((st = qaws_mesh_3d_add_attribute(&a, m, fv, &ad)) != QAWS_STATUS_OK)
		goto done;
	memset(&d, 0, sizeof d);
	d.max_error = (qaws_scalar)(1e-3 * tol);
	if ((st = qaws_mesh_3d_decimate(NULL, m, &d)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!swp_flat_delaunay(&m->hm, 1e-3 * tol))
		goto done;
	ne = qaws_mesh_3d_get_edge_count(m);
	qaws_internal_free(ev);
	qaws_internal_free(fv);
	ev = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (2 * (size_t)ne + 1));
	fv = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * ((size_t)ne + 1));
	if (!ev || !fv || qaws_mesh_3d_find_attribute(&a, m, "qaws_sweep_feature") != QAWS_STATUS_OK)
		goto done;
	qaws_mesh_3d_get_edges(ev, m);
	qaws_mesh_3d_get_attribute(fv, m, a);
	for (i = 0; i < ne; i++)
		if (fv[i] != 0)
		{
			ev[2 * (size_t)nf] = ev[2 * (size_t)i];
			ev[2 * (size_t)nf + 1] = ev[2 * (size_t)i + 1];
			nf++;
		}
	if ((st = qaws_mesh_3d_remove_attribute(m, a)) == QAWS_STATUS_OK)
		st = qaws_mesh_3d_set_feature_edges(m, ev, nf);
done:
	qaws_internal_free(ev);
	qaws_internal_free(fv);
	return st;
}

/* ---- spins: a body turning a whole turn or more about one axis while drifting at most the tolerance a turn passes over
   the same places again and again. Its sweep is its solid of revolution about the axis (each turn reaches every angle)
   swept along the drift; the solid of revolution is made from the body's profile in the half plane of the axis (height,
   distance), the union over heights of each slice's pieces' distance ranges ---- */

/* the motion a spin: rigid, its angular velocity keeping one direction u and sense, a whole turn or more; its turns (0: not
   a spin), *slow whether it moves at most tol along any turn */
static double swp_spin(qaws_motion const* m, double tol, double* u, int* slow)
{
	double total = 0, prev = 0, l;
	unsigned int i, k, ns = 512;
	qaws_motion_state s;
	*slow = 1;
	if (qaws_motion_is_per_vertex(m) || qaws_motion_evaluate(&s, m, (qaws_scalar)0.5) != QAWS_STATUS_OK)
		return 0;
	u[0] = (double)s.angular_velocity.x;
	u[1] = (double)s.angular_velocity.y;
	u[2] = (double)s.angular_velocity.z;
	l = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
	if (!(l > 0))
		return 0;
	for (k = 0; k < 3; k++)
		u[k] /= l;
	for (i = 0; i <= ns; i++)
	{
		double w[3], c[3], a, lw, v;
		if (qaws_motion_evaluate(&s, m, (qaws_scalar)((double)i / ns)) != QAWS_STATUS_OK)
			return 0;
		w[0] = (double)s.angular_velocity.x;
		w[1] = (double)s.angular_velocity.y;
		w[2] = (double)s.angular_velocity.z;
		lw = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
		a = w[0] * u[0] + w[1] * u[1] + w[2] * u[2];
		c[0] = w[1] * u[2] - w[2] * u[1];
		c[1] = w[2] * u[0] - w[0] * u[2];
		c[2] = w[0] * u[1] - w[1] * u[0];
		v = sqrt((double)s.velocity.x * s.velocity.x + (double)s.velocity.y * s.velocity.y + (double)s.velocity.z * s.velocity.z);
		if (!(a > 0) || sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) > 1e-6 * lw)
			return 0;
		if (v * 2 * 3.14159265358979323846 / a > tol)
			*slow = 0;
		if (i)
			total += 0.5 * (prev + a) / ns;
		prev = a;
	}
	return total >= 2 * 3.14159265358979323846 ? total / (2 * 3.14159265358979323846) : 0;
}

/* a slice's piece: the edge key it is known by (the least of its outer loop's), its distance range from the axis */
typedef struct swp_rv_piece
{
	unsigned long long key;
	double rmin, rmax;
} swp_rv_piece;

static int swp_rv_piece_cmp(void const* p, void const* q)
{
	unsigned long long a = ((swp_rv_piece const*)p)->key, b = ((swp_rv_piece const*)q)->key;
	return a < b ? -1 : a > b;
}

static int swp_rv_u64_cmp(void const* p, void const* q)
{
	unsigned long long a = *(unsigned long long const*)p, b = *(unsigned long long const*)q;
	return a < b ? -1 : a > b;
}

/* point (x, y) in the polygon P (n points), even-odd */
static int swp_rv_inside(double const* P, unsigned int n, double x, double y)
{
	unsigned int i, j;
	int in = 0;
	for (i = 0, j = n - 1; i < n; j = i++)
		if ((P[2 * i + 1] > y) != (P[2 * j + 1] > y) && x < (P[2 * j] - P[2 * i]) * (y - P[2 * i + 1]) / (P[2 * j + 1] - P[2 * i + 1]) + P[2 * i])
			in = !in;
	return in;
}

/* the distance from the origin to segment p q */
static double swp_rv_segdist(double const* p, double const* q)
{
	double d[2] = { q[0] - p[0], q[1] - p[1] }, l = d[0] * d[0] + d[1] * d[1], s = l > 0 ? -(p[0] * d[0] + p[1] * d[1]) / l : 0;
	s = s < 0 ? 0 : s > 1 ? 1 : s;
	return hypot(p[0] + s * d[0], p[1] + s * d[1]);
}

/* the slice at height h (no vertex there): its pieces into *out (grown), their count; Z, A, B the vertices' height and plane
   coordinates. Scratch arrays grown as needed. -1 out of memory */
static int swp_rv_slice(swp_rv_piece** out, unsigned int* capo, double h, double const* Z, double const* A, double const* B, unsigned int const* t,
	unsigned int nt, unsigned long long** ends, unsigned int* cape, double** pts, unsigned int* capp, unsigned int** loop, unsigned int* capl)
{
	unsigned int f, ne = 0, i, nl = 0, np = 0, npc = 0, j;
	unsigned char* used;
	/* each crossing triangle's segment: its two crossed edges (keys), 2 ends a segment */
	for (f = 0; f < nt; f++)
	{
		unsigned int k, c = 0;
		for (k = 0; k < 3; k++)
			c += Z[t[3 * (size_t)f + k]] > h;
		if (c == 0 || c == 3)
			continue;
		if (ne + 2 > *cape)
		{
			unsigned int cap = 2 * *cape + 64;
			unsigned long long* g = (unsigned long long*)qaws_internal_realloc(*ends, sizeof(unsigned long long) * 2 * (size_t)cap);
			if (!g)
				return -1;
			*ends = g;
			*cape = cap;
		}
		for (k = 0; k < 3; k++)
		{
			unsigned int a = t[3 * (size_t)f + k], b = t[3 * (size_t)f + (k + 1) % 3];
			if ((Z[a] > h) != (Z[b] > h))
			{
				unsigned long long key = a < b ? ((unsigned long long)a << 32 | b) : ((unsigned long long)b << 32 | a);
				(*ends)[2 * (size_t)ne] = key;
				(*ends)[2 * (size_t)ne + 1] = ne;   /* (the segment is ne / 2, its end ne % 2) */
				ne++;
			}
		}
	}
	if (!ne)
		return 0;
	/* the loops: segments chained by their shared edge keys */
	{
		unsigned long long* sorted = (unsigned long long*)qaws_internal_malloc(sizeof(unsigned long long) * 2 * (size_t)ne);
		unsigned int* mate = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (size_t)ne);
		used = (unsigned char*)qaws_internal_zalloc((size_t)ne / 2 + 1, 1);
		if (!sorted || !mate || !used)
		{
			qaws_internal_free(sorted); qaws_internal_free(mate); qaws_internal_free(used);
			return -1;
		}
		memcpy(sorted, *ends, sizeof(unsigned long long) * 2 * (size_t)ne);
		qsort(sorted, ne, 2 * sizeof(unsigned long long), swp_rv_u64_cmp);
		for (i = 0; i < ne; i++)
			mate[i] = 0xFFFFFFFFu;
		for (i = 0; i + 1 < ne; i++)
			if (sorted[2 * (size_t)i] == sorted[2 * (size_t)i + 2])
			{
				mate[(unsigned int)sorted[2 * (size_t)i + 1]] = (unsigned int)sorted[2 * (size_t)i + 3];
				mate[(unsigned int)sorted[2 * (size_t)i + 3]] = (unsigned int)sorted[2 * (size_t)i + 1];
				i++;
			}
		for (i = 0; i < ne / 2; i++)
		{
			unsigned int e, start = np;
			unsigned long long kmin = ~0ull;
			if (used[i])
				continue;
			/* a loop from segment i: leave by its end 1, enter the mate's segment, leave by its other end */
			e = 2 * i + 1;
			for (;;)
			{
				unsigned long long key = (*ends)[2 * (size_t)e];
				unsigned int a = (unsigned int)(key >> 32), b = (unsigned int)key, m2;
				double s = (h - Z[a]) / (Z[b] - Z[a]);
				used[e / 2] = 1;
				if (np + 1 > *capp)
				{
					unsigned int cap = 2 * *capp + 64;
					double* g = (double*)qaws_internal_realloc(*pts, sizeof(double) * 2 * (size_t)cap);
					if (!g)
					{
						qaws_internal_free(sorted); qaws_internal_free(mate); qaws_internal_free(used);
						return -1;
					}
					*pts = g;
					*capp = cap;
				}
				(*pts)[2 * (size_t)np] = A[a] + s * (A[b] - A[a]);
				(*pts)[2 * (size_t)np + 1] = B[a] + s * (B[b] - B[a]);
				np++;
				kmin = key < kmin ? key : kmin;
				m2 = mate[e];
				if (m2 == 0xFFFFFFFFu || used[m2 / 2])
					break;   /* (closed, or an open chain on a mesh not closed) */
				e = m2 ^ 1u;
			}
			if (nl + 1 > *capl)
			{
				unsigned int cap = 2 * *capl + 16;
				unsigned int* g = (unsigned int*)qaws_internal_realloc(*loop, sizeof(unsigned int) * 4 * (size_t)cap);
				if (!g)
				{
					qaws_internal_free(sorted); qaws_internal_free(mate); qaws_internal_free(used);
					return -1;
				}
				*loop = g;
				*capl = cap;
			}
			(*loop)[4 * (size_t)nl] = start;
			(*loop)[4 * (size_t)nl + 1] = np - start;
			(*loop)[4 * (size_t)nl + 2] = (unsigned int)(kmin >> 32);
			(*loop)[4 * (size_t)nl + 3] = (unsigned int)kmin;
			nl++;
		}
		qaws_internal_free(sorted); qaws_internal_free(mate); qaws_internal_free(used);
	}
	/* the pieces: each loop's depth (loops round it); an even one an outer loop, its holes the odd ones directly in it */
	for (i = 0; i < nl; i++)
	{
		double const* P = *pts + 2 * (size_t)(*loop)[4 * (size_t)i];
		unsigned int n = (*loop)[4 * (size_t)i + 1], depth = 0, k;
		double rmin, rmax = 0;
		int axis;
		if (n < 3)
			continue;
		for (j = 0; j < nl; j++)
			if (j != i && (*loop)[4 * (size_t)j + 1] >= 3 && swp_rv_inside(*pts + 2 * (size_t)(*loop)[4 * (size_t)j], (*loop)[4 * (size_t)j + 1], P[0], P[1]))
				depth++;
		if (depth & 1)
			continue;
		axis = swp_rv_inside(P, n, 0, 0);
		rmin = HUGE_VAL;
		for (k = 0; k < n; k++)
		{
			rmax = fmax(rmax, hypot(P[2 * k], P[2 * k + 1]));
			rmin = fmin(rmin, swp_rv_segdist(P + 2 * k, P + 2 * ((k + 1) % n)));
		}
		/* its holes: the loops inside it one level deeper */
		for (j = 0; j < nl; j++)
		{
			double const* H = *pts + 2 * (size_t)(*loop)[4 * (size_t)j];
			unsigned int nh = (*loop)[4 * (size_t)j + 1], dh = 0, q;
			if (j == i || nh < 3 || !swp_rv_inside(P, n, H[0], H[1]))
				continue;
			for (q = 0; q < nl; q++)
				if (q != j && (*loop)[4 * (size_t)q + 1] >= 3 && swp_rv_inside(*pts + 2 * (size_t)(*loop)[4 * (size_t)q], (*loop)[4 * (size_t)q + 1], H[0], H[1]))
					dh++;
			if (dh != depth + 1)
				continue;
			if (swp_rv_inside(H, nh, 0, 0))
				axis = 0;
			for (q = 0; q < nh; q++)
				rmin = fmin(rmin, swp_rv_segdist(H + 2 * q, H + 2 * ((q + 1) % nh)));
		}
		if (npc + 1 > *capo)
		{
			unsigned int cap = 2 * *capo + 16;
			swp_rv_piece* g = (swp_rv_piece*)qaws_internal_realloc(*out, sizeof(swp_rv_piece) * (size_t)cap);
			if (!g)
				return -1;
			*out = g;
			*capo = cap;
		}
		(*out)[npc].key = (unsigned long long)(*loop)[4 * (size_t)i + 2] << 32 | (*loop)[4 * (size_t)i + 3];
		(*out)[npc].rmin = axis ? 0 : rmin;
		(*out)[npc].rmax = rmax;
		npc++;
	}
	qsort(*out, npc, sizeof(swp_rv_piece), swp_rv_piece_cmp);
	return (int)npc;
}

/* where a point of the solid of revolution comes from: per profile point the body's surface points on its circle (angle,
   triangle, barycentric weights), flat triangles across it kept whole (any angle inside them is its own) */
typedef struct swp_rv_src
{
	double th, w[3], slack;
	unsigned int f;
	int flat;
} swp_rv_src;

/* the body's points on the circle of height z0, radius r0 about the axis: into *src (grown), their count. -1 out of memory */
static int swp_rv_sources(swp_rv_src** src, unsigned int* cap, double z0, double r0, double const* Z, double const* A, double const* B,
	unsigned int const* t, unsigned int nt, double tz, double tol)
{
	unsigned int f, n = 0, k;
	for (f = 0; f < nt; f++)
	{
		unsigned int v[3] = { t[3 * (size_t)f], t[3 * (size_t)f + 1], t[3 * (size_t)f + 2] };
		double zl = fmin(Z[v[0]], fmin(Z[v[1]], Z[v[2]])), zh = fmax(Z[v[0]], fmax(Z[v[1]], Z[v[2]])), P[3][2], np = 0;
		double cand[8][3];   /* (angle, slack, the point's s along the slice: unused) */
		unsigned int nc = 0, q;
		if (z0 < zl - tz || z0 > zh + tz)
			continue;
		if (n + 9 > *cap)
		{
			unsigned int c2 = 2 * *cap + 32;
			swp_rv_src* g = (swp_rv_src*)qaws_internal_realloc(*src, sizeof(swp_rv_src) * (size_t)c2);
			if (!g)
				return -1;
			*src = g;
			*cap = c2;
		}
		if (zh - zl <= tz)
		{
			(*src)[n].f = f;
			(*src)[n].flat = 1;
			(*src)[n].slack = 0;
			n++;
		}
		/* the slice's segment: where the edges meet the height (a vertex at it counts) */
		for (k = 0; k < 3; k++)
		{
			unsigned int a = v[k], b = v[(k + 1) % 3];
			double za = Z[a], zb = Z[b];
			if (fabs(za - z0) <= tz)
			{
				P[(unsigned int)np][0] = A[a];
				P[(unsigned int)np][1] = B[a];
				np += np < 3;
			}
			else if ((za - z0) * (zb - z0) < 0 && fabs(zb - z0) > tz)
			{
				double s = (z0 - za) / (zb - za);
				P[(unsigned int)np][0] = A[a] + s * (A[b] - A[a]);
				P[(unsigned int)np][1] = B[a] + s * (B[b] - B[a]);
				np += np < 3;
			}
		}
		/* the points of each edge of what was found at distance r0 (the nearest within tol when none is) */
		for (k = 0; k < (unsigned int)np; k++)
		{
			double const *p = P[k], *qp = P[(k + 1) % (unsigned int)np];
			double d[2] = { qp[0] - p[0], qp[1] - p[1] }, aa = d[0] * d[0] + d[1] * d[1], bb = 2 * (p[0] * d[0] + p[1] * d[1]),
				cc = p[0] * p[0] + p[1] * p[1] - r0 * r0, disc = bb * bb - 4 * aa * cc;
			int found = 0;
			if (aa > 0 && disc >= 0)
			{
				double sq = sqrt(disc), s1 = (-bb - sq) / (2 * aa), s2 = (-bb + sq) / (2 * aa);
				if (s1 >= -1e-12 && s1 <= 1 + 1e-12 && nc < 8)
				{
					cand[nc][0] = atan2(p[1] + s1 * d[1], p[0] + s1 * d[0]);
					cand[nc++][1] = 0;
					found = 1;
				}
				if (s2 >= -1e-12 && s2 <= 1 + 1e-12 && s2 - s1 > 1e-12 && nc < 8)
				{
					cand[nc][0] = atan2(p[1] + s2 * d[1], p[0] + s2 * d[0]);
					cand[nc++][1] = 0;
					found = 1;
				}
			}
			if (!found && nc < 8)
			{
				double s = aa > 0 ? -bb / (2 * aa) : 0, x, y, sl;
				s = s < 0 ? 0 : s > 1 ? 1 : s;
				x = p[0] + s * d[0];
				y = p[1] + s * d[1];
				/* (the nearest of the segment to the circle: its foot, or an end) */
				{
					double r1 = hypot(p[0], p[1]), r2 = hypot(qp[0], qp[1]), rs = hypot(x, y);
					double e1 = fabs(r1 - r0), e2 = fabs(r2 - r0), es = fabs(rs - r0);
					if (e1 <= e2 && e1 <= es) { x = p[0]; y = p[1]; }
					else if (e2 <= es) { x = qp[0]; y = qp[1]; }
					sl = fmin(es, fmin(e1, e2));
				}
				if (sl <= tol)
				{
					cand[nc][0] = atan2(y, x);
					cand[nc++][1] = sl;
				}
			}
			if (np == 1)
				break;
		}
		for (q = 0; q < nc; q++)
		{
			(*src)[n].f = f;
			(*src)[n].flat = 0;
			(*src)[n].th = cand[q][0];
			(*src)[n].slack = cand[q][1];
			n++;
		}
	}
	return (int)n;
}

/* barycentric weights of x in the triangle a b c (3D, projected on its plane) */
static void swp_rv_bary(double const* x, double const* a, double const* b, double const* c, double* w)
{
	double u[3], v[3], p[3], d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0, den;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = b[k] - a[k];
		v[k] = c[k] - a[k];
		p[k] = x[k] - a[k];
		d00 += u[k] * u[k]; d01 += u[k] * v[k]; d11 += v[k] * v[k]; d20 += p[k] * u[k]; d21 += p[k] * v[k];
	}
	den = d00 * d11 - d01 * d01;
	w[1] = den != 0 ? (d11 * d20 - d01 * d21) / den : 0;
	w[2] = den != 0 ? (d00 * d21 - d01 * d20) / den : 0;
	w[0] = 1 - w[1] - w[2];
	for (k = 0; k < 3; k++)
		w[k] = w[k] < 0 ? 0 : w[k];
	den = w[0] + w[1] + w[2];
	for (k = 0; k < 3; k++)
		w[k] = den > 0 ? w[k] / den : 1.0 / 3;
}

/* the solid of revolution of the body (X: its vertices at 0, t its triangles) about the axis through o along u, its
   attributes from the body's surface where the turn first brings it (the least angle behind each point: passes coincide,
   the first stands), directions turned by that angle after the motion's rotation at 0 (R0) */
static qaws_status swp_revolve(qaws_mesh_3d** out, qaws_mesh_3d const* mesh, double const* X, unsigned int nv, unsigned int const* t,
	unsigned int nt, double const* o, double const* u, double const* R0, double tol, int symmetric_only)
{
	double e1[3], e2[3], l, ext = 0, tz, q, rmaxall = 0;
	double *Z = NULL, *A = NULL, *Bc = NULL, *hs = NULL, *pts = NULL, *prof = NULL, *vx = NULL, *alpha = NULL, *wsrc = NULL;
	unsigned int i, k, nh = 0, *loop = NULL, capl = 0, capp = 0, cape = 0, capo[2] = { 0, 0 }, nquad = 0, capq = 0, nprof = 0, npath = 0;
	unsigned int *pathoff = NULL, nvo = 0, nto = 0, *tri = NULL, *ring = NULL, N, *fsrc = NULL, capsrc = 0;
	unsigned long long* ends = NULL;
	swp_rv_piece* pc[2] = { NULL, NULL };
	int64_t* quads = NULL;
	qaws_path64* paths = NULL;
	qaws_clip64_result* cr = NULL;
	swp_rv_src* src = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	*out = NULL;
	/* the frame: u, e1, e2 */
	e1[0] = fabs(u[0]) < 0.9 ? 1 : 0;
	e1[1] = fabs(u[0]) < 0.9 ? 0 : 1;
	e1[2] = 0;
	l = e1[0] * u[0] + e1[1] * u[1] + e1[2] * u[2];
	for (k = 0; k < 3; k++)
		e1[k] -= l * u[k];
	l = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
	for (k = 0; k < 3; k++)
		e1[k] /= l;
	e2[0] = u[1] * e1[2] - u[2] * e1[1];
	e2[1] = u[2] * e1[0] - u[0] * e1[2];
	e2[2] = u[0] * e1[1] - u[1] * e1[0];
	Z = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	A = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	Bc = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	hs = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nv + 1));
	if (!Z || !A || !Bc || !hs)
		goto done;
	for (i = 0; i < nv; i++)
	{
		double d[3];
		for (k = 0; k < 3; k++)
			d[k] = X[3 * (size_t)i + k] - o[k];
		Z[i] = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
		A[i] = d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2];
		Bc[i] = d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2];
		ext = fmax(ext, fmax(fabs(Z[i]), hypot(A[i], Bc[i])));
		rmaxall = fmax(rmaxall, hypot(A[i], Bc[i]));
		hs[i] = Z[i];
	}
	if (!(ext > 0))
		goto bad;
	tz = 1e-12 * ext;
	q = 1e-7 * ext;   /* (the profile's integer grid) */
	qsort(hs, nv, sizeof(double), swp_cmp_double);
	for (i = 0; i < nv; i++)
		if (!nh || hs[i] - hs[nh - 1] > tz)
			hs[nh++] = hs[i];
	/* the profile's pieces: per gap between vertex heights its slices (the topology holds within one), each piece's distance
	   range joined to the next slice's as a quad (height, distance), the gap's ends at the vertex heights */
	for (i = 0; i + 1 < nh; i++)
	{
		double g = hs[i + 1] - hs[i];
		unsigned int m = (unsigned int)ceil(g / (2 * tol)), s, cur = 0, np0 = 0;
		double zprev = hs[i];
		m = m < 1 ? 1 : m > 64 ? 64 : m;
		for (s = 0; s <= m; s++)
		{
			double h = s == 0 ? hs[i] + 1e-7 * g : s == m ? hs[i + 1] - 1e-7 * g : hs[i] + g * s / m, zz = s == 0 ? hs[i] : s == m ? hs[i + 1] : h;
			int npc = swp_rv_slice(&pc[cur], &capo[cur], h, Z, A, Bc, t, nt, &ends, &cape, &pts, &capp, &loop, &capl);
			unsigned int a2 = 0, b2 = 0;
			if (npc < 0)
				goto done;
			if (s > 0)
				while (a2 < np0 && b2 < (unsigned int)npc)
				{
					swp_rv_piece const *p0 = pc[cur ^ 1] + a2, *p1 = pc[cur] + b2;
					if (p0->key < p1->key)
						a2++;
					else if (p1->key < p0->key)
						b2++;
					else
					{
						int64_t* Q;
						if (nquad + 1 > capq)
						{
							unsigned int c2 = 2 * capq + 64;
							int64_t* gq = (int64_t*)qaws_internal_realloc(quads, sizeof(int64_t) * 8 * (size_t)c2);
							if (!gq)
								goto done;
							quads = gq;
							capq = c2;
						}
						Q = quads + 8 * (size_t)nquad++;
						Q[0] = (int64_t)llround(zprev / q); Q[1] = (int64_t)llround(p0->rmin / q);
						Q[2] = (int64_t)llround(zz / q);    Q[3] = (int64_t)llround(p1->rmin / q);
						Q[4] = (int64_t)llround(zz / q);    Q[5] = (int64_t)llround(p1->rmax / q);
						Q[6] = (int64_t)llround(zprev / q); Q[7] = (int64_t)llround(p0->rmax / q);
						a2++;
						b2++;
					}
				}
			zprev = zz;
			np0 = (unsigned int)npc;
			cur ^= 1;
		}
	}
	if (!nquad)
		goto bad;
	/* their union: the profile (outer paths and holes) */
	paths = (qaws_path64*)qaws_internal_malloc(sizeof(qaws_path64) * (size_t)nquad);
	if (!paths)
		goto done;
	for (i = 0; i < nquad; i++)
	{
		paths[i].points = quads + 8 * (size_t)i;
		paths[i].point_count = 4;
	}
	{
		qaws_clip64_desc cd;
		unsigned int nr, tot = 0;
		memset(&cd, 0, sizeof cd);
		cd.subjects = paths;
		cd.subject_count = nquad;
		cd.clip_type = QAWS_CLIP_UNION;
		cd.fill_rule = QAWS_FILL_NON_ZERO;
		if ((st = qaws_clip64_execute(&cr, &cd)) != QAWS_STATUS_OK)
			goto done;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		nr = qaws_clip64_result_get_path_count(cr);
		for (i = 0; i < nr; i++)
		{
			double const* p;
			unsigned int n;
			if (qaws_clip64_result_get_path_exact(&p, &n, cr, i) == QAWS_STATUS_OK)
				tot += n;
		}
		prof = (double*)qaws_internal_malloc(sizeof(double) * 2 * ((size_t)tot + 1));
		pathoff = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nr + 2));
		if (!prof || !pathoff)
			goto done;
		for (i = 0; i < nr; i++)
		{
			double const* p;
			unsigned int n, j;
			if (qaws_clip64_result_get_path_exact(&p, &n, cr, i) != QAWS_STATUS_OK || n < 3)
				continue;
			pathoff[npath++] = nprof;
			for (j = 0; j < n; j++)
			{
				prof[2 * (size_t)nprof] = p[2 * j] * q;
				prof[2 * (size_t)nprof + 1] = p[2 * j + 1] * q <= 1e-6 * ext ? 0 : p[2 * j + 1] * q;   /* (on the axis) */
				nprof++;
			}
			{
				/* the points on a straight line between their neighbours (each slice leaves one on a straight side) left out */
				unsigned int s0 = pathoff[npath - 1], cnt = nprof - s0, changed = 1;
				while (changed && cnt > 3)
				{
					changed = 0;
					for (j = 0; j < cnt && cnt > 3;)
					{
						double const *a = prof + 2 * ((size_t)s0 + (j + cnt - 1) % cnt), *b = prof + 2 * ((size_t)s0 + j), *c = prof + 2 * ((size_t)s0 + (j + 1) % cnt);
						double dx = c[0] - a[0], dy = c[1] - a[1], ll = hypot(dx, dy);
						if (ll > 0 && fabs(dx * (b[1] - a[1]) - dy * (b[0] - a[0])) / ll <= 2 * q)
						{
							memmove(prof + 2 * ((size_t)s0 + j), prof + 2 * ((size_t)s0 + j + 1), sizeof(double) * 2 * (cnt - j - 1));
							cnt--;
							changed = 1;
						}
						else
							j++;
					}
				}
				nprof = s0 + cnt;
			}
		}
		pathoff[npath] = nprof;
	}
	if (symmetric_only)
	{
		/* a fast drift: the revolution stands for the body only if the body is one (every vertex, edge middle and face centre
		   within tol of the profile's border: a drill, an end mill, a wheel), else the slabs (UNSUPPORTED_OPERATION) */
		unsigned int f2, s2;
		for (f2 = 0; f2 < nt; f2++)
			for (s2 = 0; s2 < 7; s2++)
			{
				/* (the corners, the edges' middles, the centre) */
				static double const wt[7][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0.5, 0.5, 0 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 }, { 1.0 / 3, 1.0 / 3, 1.0 / 3 } };
				double zz = 0, aa = 0, bb = 0, rr, best = HUGE_VAL;
				unsigned int pth, j;
				for (k = 0; k < 3; k++)
				{
					unsigned int v = t[3 * (size_t)f2 + k];
					zz += wt[s2][k] * Z[v];
					aa += wt[s2][k] * A[v];
					bb += wt[s2][k] * Bc[v];
				}
				rr = hypot(aa, bb);
				for (pth = 0; pth < npath && best > tol; pth++)
					for (j = pathoff[pth]; j < pathoff[pth + 1] && best > tol; j++)
					{
						unsigned int j2 = j + 1 < pathoff[pth + 1] ? j + 1 : pathoff[pth];
						double const *pa = prof + 2 * (size_t)j, *pb = prof + 2 * (size_t)j2;
						double dz = pb[0] - pa[0], dr = pb[1] - pa[1], ll = dz * dz + dr * dr, s3 = ll > 0 ? ((zz - pa[0]) * dz + (rr - pa[1]) * dr) / ll : 0;
						s3 = s3 < 0 ? 0 : s3 > 1 ? 1 : s3;
						best = fmin(best, hypot(zz - pa[0] - s3 * dz, rr - pa[1] - s3 * dr));
					}
				if (best > tol)
				{
					st = QAWS_STATUS_UNSUPPORTED_OPERATION;
					goto done;
				}
			}
	}
	/* the revolution: N segments, the chords within tol of the widest circle; a profile point on the axis one vertex */
	N = tol < rmaxall ? (unsigned int)ceil(3.14159265358979323846 / acos(1 - tol / rmaxall)) : 8;
	N = N < 8 ? 8 : N > 2048 ? 2048 : N;
	ring = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nprof + 1));
	vx = (double*)qaws_internal_malloc(sizeof(double) * 3 * ((size_t)nprof * N + 1));
	tri = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * 6 * ((size_t)nprof * N + 1));
	alpha = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)nprof * N + 1));
	wsrc = (double*)qaws_internal_malloc(sizeof(double) * 3 * ((size_t)nprof * N + 1));
	fsrc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nprof * N + 1));
	if (!ring || !vx || !tri || !alpha || !wsrc || !fsrc)
		goto done;
	for (i = 0; i < nprof; i++)
	{
		double z = prof[2 * (size_t)i], r = prof[2 * (size_t)i + 1];
		unsigned int nk = r > 0 ? N : 1, kk;
		int ns = swp_rv_sources(&src, &capsrc, z, r, Z, A, Bc, t, nt, 1e-9 * ext, tol);
		if (ns < 0)
			goto done;
		ring[i] = nvo;
		for (kk = 0; kk < nk; kk++)
		{
			double th = 2 * 3.14159265358979323846 * kk / N, ct = cos(th), sn = sin(th), best = HUGE_VAL, bth = th, x[3];
			unsigned int bf = 0xFFFFFFFFu, j;
			for (k = 0; k < 3; k++)
				vx[3 * (size_t)nvo + k] = o[k] + z * u[k] + r * (ct * e1[k] + sn * e2[k]);
			/* the source: the body's point on this circle the turn brings here first (the least angle behind; a flat
			   triangle holding this point itself) */
			for (j = 0; j < (unsigned int)ns; j++)
			{
				double dlt;
				if (src[j].flat)
				{
					unsigned int const* v = t + 3 * (size_t)src[j].f;
					double P[3][2], pa = r * ct, pb = r * sn, s0, s1, s2;
					unsigned int z2;
					for (z2 = 0; z2 < 3; z2++)
					{
						P[z2][0] = A[v[z2]];
						P[z2][1] = Bc[v[z2]];
					}
					s0 = (P[1][0] - P[0][0]) * (pb - P[0][1]) - (P[1][1] - P[0][1]) * (pa - P[0][0]);
					s1 = (P[2][0] - P[1][0]) * (pb - P[1][1]) - (P[2][1] - P[1][1]) * (pa - P[1][0]);
					s2 = (P[0][0] - P[2][0]) * (pb - P[2][1]) - (P[0][1] - P[2][1]) * (pa - P[2][0]);
					if (!((s0 >= 0 && s1 >= 0 && s2 >= 0) || (s0 <= 0 && s1 <= 0 && s2 <= 0)))
						continue;
					dlt = 0;
					if (dlt < best)
					{
						best = dlt;
						bf = j;
						bth = th;
					}
					continue;
				}
				dlt = fmod(th - src[j].th, 2 * 3.14159265358979323846);
				if (dlt < 0)
					dlt += 2 * 3.14159265358979323846;
				if (dlt > 2 * 3.14159265358979323846 - 1e-9)
					dlt = 0;
				dlt += 1e-3 * src[j].slack / (r > 0 ? r : 1);   /* (the nearer of two at one angle) */
				if (dlt < best)
				{
					best = dlt;
					bf = j;
					bth = src[j].th;
				}
			}
			if (bf == 0xFFFFFFFFu)
			{
				/* (none: the nearest vertex of the body by height and distance) */
				double bd = HUGE_VAL;
				unsigned int bv = 0, ff;
				for (j = 0; j < nv; j++)
				{
					double dd = (Z[j] - z) * (Z[j] - z) + (hypot(A[j], Bc[j]) - r) * (hypot(A[j], Bc[j]) - r);
					if (dd < bd)
					{
						bd = dd;
						bv = j;
					}
				}
				for (ff = 0; ff < nt && t[3 * (size_t)ff] != bv && t[3 * (size_t)ff + 1] != bv && t[3 * (size_t)ff + 2] != bv; ff++)
					;
				fsrc[nvo] = ff < nt ? ff : 0;
				for (k = 0; k < 3; k++)
					wsrc[3 * (size_t)nvo + k] = ff < nt ? (t[3 * (size_t)ff + k] == bv) : 1.0 / 3;
				alpha[nvo] = th - atan2(Bc[bv], A[bv]);
			}
			else
			{
				unsigned int f = src[bf].f;
				for (k = 0; k < 3; k++)
					x[k] = o[k] + z * u[k] + r * (cos(bth) * e1[k] + sin(bth) * e2[k]);
				swp_rv_bary(x, X + 3 * (size_t)t[3 * (size_t)f], X + 3 * (size_t)t[3 * (size_t)f + 1], X + 3 * (size_t)t[3 * (size_t)f + 2], wsrc + 3 * (size_t)nvo);
				fsrc[nvo] = f;
				alpha[nvo] = th - bth;
			}
			nvo++;
		}
	}
	/* the faces round each profile edge */
	for (k = 0; k < npath; k++)
	{
		unsigned int a0 = pathoff[k], a1 = pathoff[k + 1];
		for (i = a0; i < a1; i++)
		{
			unsigned int j = i + 1 < a1 ? i + 1 : a0, kk;
			int ai = prof[2 * (size_t)i + 1] <= 0, aj = prof[2 * (size_t)j + 1] <= 0;
			if (ai && aj)
				continue;
			for (kk = 0; kk < N; kk++)
			{
				unsigned int k1 = (kk + 1) % N, p0 = ring[i] + (ai ? 0 : kk), p1 = ring[i] + (ai ? 0 : k1), q0 = ring[j] + (aj ? 0 : kk), q1 = ring[j] + (aj ? 0 : k1);
				if (!ai)
				{
					tri[3 * (size_t)nto] = p0; tri[3 * (size_t)nto + 1] = q0; tri[3 * (size_t)nto + 2] = p1;
					nto++;
				}
				if (!aj)
				{
					tri[3 * (size_t)nto] = p1; tri[3 * (size_t)nto + 1] = q0; tri[3 * (size_t)nto + 2] = q1;
					nto++;
				}
			}
		}
	}
	{
		/* outward: the signed volume */
		double vol = 0;
		for (i = 0; i < nto; i++)
		{
			double const *a = vx + 3 * (size_t)tri[3 * (size_t)i], *b = vx + 3 * (size_t)tri[3 * (size_t)i + 1], *c = vx + 3 * (size_t)tri[3 * (size_t)i + 2];
			vol += a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
		}
		if (vol < 0)
			for (i = 0; i < nto; i++)
			{
				unsigned int z2 = tri[3 * (size_t)i + 1];
				tri[3 * (size_t)i + 1] = tri[3 * (size_t)i + 2];
				tri[3 * (size_t)i + 2] = z2;
			}
	}
	{
		qaws_mesh_3d_desc md;
		qaws_scalar* xs = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * 3 * ((size_t)nvo + 1));
		if (!xs)
			goto done;
		for (i = 0; i < 3 * nvo; i++)
			xs[i] = (qaws_scalar)vx[i];
		memset(&md, 0, sizeof md);
		md.positions = xs;
		md.vertex_count = nvo;
		md.face_vertices = tri;
		md.face_count = nto;
		st = qaws_mesh_3d_create(out, &md);
		qaws_internal_free(xs);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	/* the attributes from their sources (vertex ones per corner, face ones from the first corner's source) */
	{
		unsigned int na = qaws_mesh_3d_get_attribute_count(mesh), nf = qaws_mesh_3d_get_face_count(mesh), nc = qaws_mesh_3d_get_corner_count(mesh), a;
		unsigned int *fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nc + 1)), *fs = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
		unsigned int *tf = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)nt + 1)), *tc = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * 3 * ((size_t)nt + 1));
		size_t cmax = (size_t)(nv > nc ? (nv > nf ? nv : nf) : (nc > nf ? nc : nf)) + 1;
		qaws_scalar* val = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * 16 * cmax);
		double* dv = (double*)qaws_internal_malloc(sizeof(double) * 16 * cmax);
		qaws_scalar* sv = (qaws_scalar*)qaws_internal_malloc(sizeof(qaws_scalar) * 16 * (3 * (size_t)nto + (size_t)qaws_mesh_3d_get_edge_count(*out) + 1));
		unsigned int j, c2;
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		if (!fv || !fs || !tf || !tc || !val || !dv || !sv)
		{
			qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(tf); qaws_internal_free(tc); qaws_internal_free(val); qaws_internal_free(dv); qaws_internal_free(sv);
			goto done;
		}
		qaws_mesh_3d_get_faces(fv, fs, mesh);
		for (i = 0, j = 0, k = 0; i < nf; k += fs[i++])
			for (c2 = 0; c2 + 2 < fs[i]; c2++, j++)
			{
				tf[j] = i;
				tc[3 * (size_t)j] = k;
				tc[3 * (size_t)j + 1] = k + c2 + 1;
				tc[3 * (size_t)j + 2] = k + c2 + 2;
			}
		st = QAWS_STATUS_OK;
		for (a = 0; a < na && st == QAWS_STATUS_OK; a++)
		{
			qaws_attribute_desc d;
			unsigned int n, id, cnt, e;
			int dir;
			qaws_mesh_3d_get_attribute_desc(&d, mesh, a);
			if (d.name && !strcmp(d.name, "source"))
				continue;
			n = d.components;
			dir = d.rule == QAWS_ATTRIBUTE_NORMAL || d.rule == QAWS_ATTRIBUTE_TANGENT || d.rule == QAWS_ATTRIBUTE_VECTOR || (d.rule == QAWS_ATTRIBUTE_NORMALIZED && n == 3);
			qaws_mesh_3d_get_attribute(val, mesh, a);
			if (d.element == QAWS_MESH_EDGE)
			{
				/* (edge values do not survive the revolution: zeros) */
				memset(sv, 0, sizeof(qaws_scalar) * n * (size_t)qaws_mesh_3d_get_edge_count(*out));
				st = qaws_mesh_3d_add_attribute(&id, *out, sv, &d);
				continue;
			}
			cnt = d.element == QAWS_MESH_VERTEX ? nv : d.element == QAWS_MESH_CORNER ? nc : nf;
			for (e = 0; e < n * cnt; e++)
				dv[e] = (double)val[e];
			for (i = 0; i < (d.element == QAWS_MESH_FACE ? nto : 3 * nto); i++)
			{
				unsigned int v = d.element == QAWS_MESH_FACE ? tri[3 * (size_t)i] : tri[i], f = fsrc[v], pos[3];
				double R[9], ca = cos(alpha[v]), sa = sin(alpha[v]), Ru[9], y[16];
				qaws_scalar* o2 = sv + (size_t)n * i;
				/* the rotation by alpha about u (Rodrigues), after R0 */
				for (k = 0; k < 9; k++)
				{
					unsigned int r = k / 3, c = k % 3;
					double kx[3][3] = { { 0, -u[2], u[1] }, { u[2], 0, -u[0] }, { -u[1], u[0], 0 } };
					Ru[k] = (r == c) * ca + sa * kx[r][c] + (1 - ca) * u[r] * u[c];
				}
				for (k = 0; k < 9; k++)
					R[k] = Ru[3 * (k / 3)] * R0[k % 3] + Ru[3 * (k / 3) + 1] * R0[3 + k % 3] + Ru[3 * (k / 3) + 2] * R0[6 + k % 3];
				if (d.element == QAWS_MESH_FACE)
					for (k = 0; k < n; k++)
						y[k] = dv[(size_t)n * tf[f] + k];
				else
				{
					for (k = 0; k < 3; k++)
						pos[k] = d.element == QAWS_MESH_VERTEX ? fv[tc[3 * (size_t)f + k]] : tc[3 * (size_t)f + k];
					hm_blend(y, dv, n, pos, wsrc + 3 * (size_t)v, 3, n, d.rule);
				}
				for (k = 0; k < n; k++)
					o2[k] = (qaws_scalar)y[k];
				if (dir)
					for (k = 0; k < 3; k++)
						o2[k] = (qaws_scalar)(R[3 * k] * y[0] + R[3 * k + 1] * y[1] + R[3 * k + 2] * y[2]);
				if (d.rule == QAWS_ATTRIBUTE_NORMAL)
				{
					/* a source on a crease (a corner or edge tracing the surface) has no normal of it: the face's own where the turned
					   one strays past 30 degrees */
					unsigned int const* tv = tri + 3 * (size_t)(d.element == QAWS_MESH_FACE ? i : i / 3);
					double const *pa = vx + 3 * (size_t)tv[0], *pb = vx + 3 * (size_t)tv[1], *pq = vx + 3 * (size_t)tv[2];
					double ea[3], eb[3], fn[3], lf, lo = 0, dt = 0;
					for (k = 0; k < 3; k++)
					{
						ea[k] = pb[k] - pa[k];
						eb[k] = pq[k] - pa[k];
					}
					fn[0] = ea[1] * eb[2] - ea[2] * eb[1];
					fn[1] = ea[2] * eb[0] - ea[0] * eb[2];
					fn[2] = ea[0] * eb[1] - ea[1] * eb[0];
					lf = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);
					for (k = 0; k < 3; k++)
					{
						lo += (double)o2[k] * o2[k];
						dt += (double)o2[k] * fn[k];
					}
					if (lf > 0 && !(dt > 0.866 * lf * sqrt(lo)))
						for (k = 0; k < 3; k++)
							o2[k] = (qaws_scalar)(fn[k] / lf);
				}
				if (d.rule == QAWS_ATTRIBUTE_QUATERNION)
				{
					qaws_scalar x4[4];
					for (k = 0; k < 4; k++)
						x4[k] = (qaws_scalar)y[k];
					swp_compose(o2, x4, R);
				}
			}
			if (d.element == QAWS_MESH_VERTEX)
				d.element = QAWS_MESH_CORNER;
			st = qaws_mesh_3d_add_attribute(&id, *out, sv, &d);
		}
		qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(tf); qaws_internal_free(tc); qaws_internal_free(val); qaws_internal_free(dv); qaws_internal_free(sv);
	}
	goto done;
bad:
	st = QAWS_STATUS_INVALID_ARGUMENT;
done:
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(*out);
		*out = NULL;
	}
	qaws_clip64_result_destroy(cr);
	qaws_internal_free(Z); qaws_internal_free(A); qaws_internal_free(Bc); qaws_internal_free(hs); qaws_internal_free(pts); qaws_internal_free(prof);
	qaws_internal_free(vx); qaws_internal_free(alpha); qaws_internal_free(wsrc); qaws_internal_free(loop); qaws_internal_free(pathoff); qaws_internal_free(tri);
	qaws_internal_free(ring); qaws_internal_free(fsrc); qaws_internal_free(ends); qaws_internal_free(pc[0]); qaws_internal_free(pc[1]); qaws_internal_free(quads);
	qaws_internal_free(paths); qaws_internal_free(src);
	return st;
}

/* a spin's sweep: the solid of revolution at 0 (features its creases sharper than half the body's softest), swept along
   the drift (the translation from 0, no turn) */
static qaws_status swp_spun(qaws_mesh_3d** out, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion, double const* X,
	unsigned int nv, unsigned int const* t, unsigned int nt, double const* u, double tol, qaws_sweep_desc const* desc, double turns, int slow)
{
	qaws_motion_state s0, s1;
	qaws_motion_pose p0, p1;
	qaws_motion* drift = NULL;
	qaws_mesh_3d* rev = NULL;
	qaws_sweep_desc d2;
	double o[3], R0[9];
	unsigned int k;
	qaws_status st;
	if ((st = qaws_motion_evaluate(&s0, motion, 0)) != QAWS_STATUS_OK || (st = qaws_motion_evaluate(&s1, motion, 1)) != QAWS_STATUS_OK)
		return st;
	o[0] = (double)s0.position.x;
	o[1] = (double)s0.position.y;
	o[2] = (double)s0.position.z;
	for (k = 0; k < 9; k++)
		R0[k] = (double)s0.rotation[k];
	if ((st = swp_revolve(&rev, mesh, X, nv, t, nt, o, u, R0, tol, !slow)) != QAWS_STATUS_OK)
		return st;
	if ((st = qaws_mesh_3d_detect_features(rev, (qaws_scalar)(0.5 * swp_feature_angle(mesh)), 0, 0)) != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(rev);
		return st;
	}
	memset(&p0, 0, sizeof p0);
	memset(&p1, 0, sizeof p1);
	p0.rotation[0] = p1.rotation[0] = 1;
	p0.velocity = s0.velocity;
	p1.velocity = s1.velocity;
	p1.position.x = (qaws_scalar)((double)s1.position.x - o[0]);
	p1.position.y = (qaws_scalar)((double)s1.position.y - o[1]);
	p1.position.z = (qaws_scalar)((double)s1.position.z - o[2]);
	if ((st = qaws_motion_create_poses(&drift, &p0, &p1, 0)) == QAWS_STATUS_OK)
	{
		if (desc)
			d2 = *desc;
		else
			memset(&d2, 0, sizeof d2);
		d2.tolerance = (qaws_scalar)tol;
		st = qaws_mesh_3d_sweep(out, report, rev, drift, &d2);
		if (report)
			report->turns = (qaws_scalar)turns;
	}
	qaws_motion_destroy(drift);
	qaws_mesh_3d_destroy(rev);
	return st;
}

qaws_status qaws_mesh_3d_sweep(qaws_mesh_3d** out, qaws_sweep_report* report, qaws_mesh_3d const* mesh, qaws_motion const* motion,
	qaws_sweep_desc const* desc)
{
	qaws_sweep_report local;
	swp_ctx c;
	unsigned int nv, nt, i, k, *t = NULL, ntimes = 0, captimes = 0, maxs, ns;
	double *Xa = NULL, *Xb = NULL, *Xm = NULL, *times = NULL, tol, lo[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL }, hi[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL }, diag = 0;
	qaws_scalar *px = NULL, *sx = NULL;
	unsigned int *sf = NULL, *sof = NULL, *soc = NULL, snv = 0, snf = 0, capv = 0, capf = 0;
	unsigned char* fr = NULL;
	unsigned int *gv = NULL, *gc = NULL, *gs = NULL, *tg = NULL, ng = 0;
	qaws_intersect_desc id;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!mesh || !motion)
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nt = qaws_mesh_3d_get_triangle_count(mesh);
	if (!nt || (qaws_motion_is_per_vertex(motion) && qaws_motion_get_vertex_count(motion) != nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&c, 0, sizeof c);
	c.m = motion;
	c.nv = nv;
	px = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (3 * (size_t)nv + 1));
	c.buf = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (3 * (size_t)nv + 1));
	t = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (3 * (size_t)nt + 1));
	Xa = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nv + 1));
	Xb = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nv + 1));
	Xm = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nv + 1));
	if (!px || !c.buf || !t || !Xa || !Xb || !Xm)
		goto done;
	qaws_mesh_3d_get_positions(px, mesh);
	qaws_mesh_3d_get_triangles(t, mesh);
	c.body = px;
	maxs = desc && desc->max_slabs ? desc->max_slabs : 256;
	/* the tolerance from the solid's box at 0 */
	if ((st = swp_at(&c, 0, Xa)) != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < nv; i++)
		for (k = 0; k < 3; k++)
		{
			lo[k] = fmin(lo[k], Xa[3 * (size_t)i + k]);
			hi[k] = fmax(hi[k], Xa[3 * (size_t)i + k]);
		}
	for (k = 0; k < 3; k++) diag += (hi[k] - lo[k]) * (hi[k] - lo[k]);
	diag = sqrt(diag);
	tol = desc && desc->tolerance > 0 ? (double)desc->tolerance : 1e-3 * diag;
	{
		/* a spin (a whole turn or more about one axis, drifting under tol a turn): its solid of revolution swept along the drift */
		/* (drifting faster: only a body symmetric about the axis; UNSUPPORTED_OPERATION otherwise, and the slabs below) */
		int slow;
		double ax[3], turns = swp_spin(motion, tol, ax, &slow);
		if (turns > 0)
		{
			st = swp_spun(out, report, mesh, motion, Xa, nv, t, nt, ax, tol, desc, turns, slow);
			if (st != QAWS_STATUS_UNSUPPORTED_OPERATION)
				goto done;
			memset(report, 0, sizeof *report);
		}
	}
	/* the slabs: [0, 1] halved until the chords' deviation (at the middle and quarters) is within tol */
	{
		double* stack = (double*)qaws_internal_scratch_malloc(sizeof(double) * 2 * (2 * (size_t)maxs + 2));
		unsigned int nsk = 0;
		if (!stack)
			goto oom;
		stack[nsk++] = 0;
		stack[nsk++] = 1;
		ns = 1;
		while (nsk)
		{
			double b = stack[--nsk], a = stack[--nsk], dev = 0;
			unsigned int q;
			if ((st = swp_at(&c, a, Xa)) != QAWS_STATUS_OK || (st = swp_at(&c, b, Xb)) != QAWS_STATUS_OK)
				break;
			for (q = 1; q <= 3 && st == QAWS_STATUS_OK; q++)
			{
				double f = 0.25 * q;
				if ((st = swp_at(&c, a + f * (b - a), Xm)) != QAWS_STATUS_OK)
					break;
				for (i = 0; i < nv; i++)
				{
					double d2 = 0;
					for (k = 0; k < 3; k++)
					{
						double e = Xm[3 * (size_t)i + k] - ((1 - f) * Xa[3 * (size_t)i + k] + f * Xb[3 * (size_t)i + k]);
						d2 += e * e;
					}
					dev = fmax(dev, sqrt(d2));
				}
			}
			if (st != QAWS_STATUS_OK)
				break;
			if (dev > tol && ns < maxs)
			{
				stack[nsk++] = 0.5 * (a + b);
				stack[nsk++] = b;
				stack[nsk++] = a;
				stack[nsk++] = 0.5 * (a + b);
				ns++;
				continue;
			}
			report->deviation = (qaws_scalar)fmax((double)report->deviation, dev);
			if (ntimes + 2 > captimes)
			{
				unsigned int cap = captimes ? 2 * captimes : 64;
				double* g = (double*)qaws_internal_realloc(times, sizeof(double) * cap);
				if (!g)
				{
					qaws_internal_free(stack);
					goto oom;
				}
				times = g;
				captimes = cap;
			}
			times[ntimes++] = a;
			times[ntimes++] = b;
		}
		qaws_internal_free(stack);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	/* the slabs in time order (pairs sorted by their start) */
	qsort(times, ntimes / 2, 2 * sizeof(double), swp_cmp_double);
	report->slabs = ntimes / 2;
	/* the soup on layers of vertices, one per slab end: the solid at 0, and per slab the convex hull of each face's two
	   positions where any of its corners moves out of the solid (its normal at the start against its displacement): the hull
	   holds every position the face passes through between them (its vertices move along their chords), never turns inside
	   out as a twisted extrusion would, and is the face's exact prism under a translation; every point the solid reaches
	   was crossed by a face moving out. Caps that coincide turned (a hull's end and the next one's start, the solid at 0's
	   face and its first hull's start) are left out in pairs, as they cancel */
	{
		size_t S = report->slabs, need_v = (size_t)nv * (S + 1), need_f = (size_t)nt * (1 + 8 * S);   /* (8 hull triangles a face and slab at most) */
		if (need_v > 0xFFFFFFF0u || need_f > 0xFFFFFFF0u)
			goto oom;
		capv = (unsigned int)need_v;
		capf = (unsigned int)need_f;
		sx = (qaws_scalar*)qaws_internal_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)capv + 1));
		sf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)capf + 1));
		fr = (unsigned char*)qaws_internal_scratch_zalloc(qaws_internal_mul2((size_t)nt, S) + 1, 1);
		gv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, (size_t)nt + 1));
		gc = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, (size_t)nt + 1));
		gs = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
		tg = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nt + 1));
		sof = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)capf + 1));
		soc = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)capf + 1));
		if (!sx || !sf || !fr || !sof || !soc || !gv || !gc || !gs || !tg)
			goto oom;
	}
	for (ns = 0; ns <= report->slabs; ns++)
	{
		/* layer ns: the vertices at the slab ends */
		if ((st = swp_at(&c, ns < report->slabs ? times[2 * (size_t)ns] : times[2 * (size_t)ns - 1], Xa)) != QAWS_STATUS_OK)
			goto done;
		for (i = 0; i < 3 * nv; i++) sx[3 * (size_t)ns * nv + i] = (qaws_scalar)Xa[i];
	}
	snv = nv * (report->slabs + 1);
	/* the flat faces: two triangles sharing an edge in one plane make a quad, hulled as one (the hulls of two halves of one
	   face would overlap in its plane, and the union resolves coincident pieces less well than none) */
	if (!(ng = swp_flat_faces(gv, gc, gs, tg, t, nt, sx, 1e-12 * diag)))
		goto oom;
	/* each face's caps per slab (fr: 1 start, 2 end), then the soup */
	for (ns = 0; ns < report->slabs; ns++)
	{
		unsigned int f, cp;
		for (f = 0; f < ng; f++)
		{
			swp_hull(sf, sof, soc, &cp, sx, gv + 4 * (size_t)f, gc + 4 * (size_t)f, gs[f], gc[4 * (size_t)f] / 3, ns * nv, (ns + 1) * nv, 1e-12 * diag, tol, 3);
			fr[(size_t)ns * nt + f] = (unsigned char)cp;
		}
	}
	snf = 0;
	for (i = 0; i < nt; i++)
	{
		if (report->slabs && (fr[tg[i]] & 1) && (fr[tg[i]] & 4))
			continue;   /* (cancelled by its face's first hull's start, turned) */
		for (k = 0; k < 3; k++)
		{
			sf[3 * (size_t)snf + k] = t[3 * (size_t)i + k];
			soc[3 * (size_t)snf + k] = 3 * i + k;
		}
		sof[snf++] = i;
	}
	for (ns = 0; ns < report->slabs; ns++)
	{
		unsigned int f, cp;
		for (f = 0; f < ng; f++)
		{
			/* a cap left out with the one it meets turned the other way (the solid at 0's face: an end cap not turned) */
			unsigned char here = fr[(size_t)ns * nt + f], next = ns + 1 < report->slabs ? fr[((size_t)ns + 1) * nt + f] : 0, prev = ns ? fr[((size_t)ns - 1) * nt + f] : 2;
			unsigned int skip = ((here & 1) && (prev & 2) && !(here & 4) != !(prev & 8) ? 1u : 0u) | ((here & 2) && (next & 1) && !(here & 8) != !(next & 4) ? 2u : 0u), n0;
			n0 = swp_hull(sf + 3 * (size_t)snf, sof + snf, soc + 3 * (size_t)snf, &cp, sx, gv + 4 * (size_t)f, gc + 4 * (size_t)f, gs[f], gc[4 * (size_t)f] / 3,
				ns * nv, (ns + 1) * nv, 1e-12 * diag, tol, skip);
			snf += n0;
			report->front_faces += n0 > 0;
		}
	}
	if ((st = swp_cancel(sf, sof, soc, &snf)) != QAWS_STATUS_OK)
		goto done;
	{
		qaws_mesh_3d_desc md;
		memset(&md, 0, sizeof md);
		md.positions = sx;
		md.vertex_count = snv;
		md.face_vertices = sf;
		md.face_count = snf;
		memset(&id, 0, sizeof id);
		id.executor = desc ? desc->executor : NULL;
		{
			/* the soup as a mesh carrying the attributes and features, its self-union carrying them on (the cuts features too) */
			qaws_mesh_3d* soup = NULL;
			for (ns = 0; ns <= report->slabs; ns++)
				times[ns] = ns < report->slabs ? times[2 * (size_t)ns] : times[2 * (size_t)ns - 1];   /* (each layer's time, in place) */
			st = swp_soup(&soup, mesh, motion, &md, sof, soc, nv, times, report->slabs + 1);
			if (st == QAWS_STATUS_OK)
				st = cy_boolean(out, &report->boolean, (qaws_mesh_3d const* const*)&soup, 1, NULL, 0, QAWS_BOOLEAN_UNION, &id);
			if (st == QAWS_STATUS_OK && *out)
				st = swp_crease_features(*out, 0.5 * swp_feature_angle(mesh), mesh, sx, nv, report->slabs);   /* (half the input's softest crease: the slabs' chords bend far less) */
			if (st == QAWS_STATUS_OK && *out && desc && desc->simplify)
				st = swp_tidy(*out, tol);
			qaws_mesh_3d_destroy(soup);
		}
		report->closed = report->boolean.closed;
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(px); qaws_internal_free(c.buf); qaws_internal_free(t);
	qaws_internal_free(Xa); qaws_internal_free(Xb); qaws_internal_free(Xm); qaws_internal_free(times);
	qaws_internal_free(gv); qaws_internal_free(gc); qaws_internal_free(gs); qaws_internal_free(tg);
	qaws_internal_free(sx); qaws_internal_free(sf); qaws_internal_free(fr); qaws_internal_free(sof); qaws_internal_free(soc);
	return st;
}

/* a closed mesh's signed volume (its triangles fanned from the origin) */
static double swp_volume(qaws_mesh_3d const* m)
{
	unsigned int nt = m ? qaws_mesh_3d_get_triangle_count(m) : 0, i;
	double v = 0;
	for (i = 0; i < nt; i++)
	{
		unsigned int c[3];
		qaws_vec3 a, b, e;
		qaws_mesh_3d_get_triangle(c, m, i);
		a = qaws_mesh_3d_get_vertex(m, c[0]);
		b = qaws_mesh_3d_get_vertex(m, c[1]);
		e = qaws_mesh_3d_get_vertex(m, c[2]);
		v += ((double)a.x * ((double)b.y * e.z - (double)b.z * e.y) + (double)a.y * ((double)b.z * e.x - (double)b.x * e.z) +
			(double)a.z * ((double)b.x * e.y - (double)b.y * e.x)) / 6;
	}
	return v;
}

qaws_status qaws_mesh_sweep_boolean(qaws_mesh_3d** out, qaws_sweep_boolean_report* report, qaws_mesh_3d const* tool, qaws_motion const* const* moves,
	qaws_mesh_3d const* stock, unsigned int move_count, qaws_boolean_op op, qaws_sweep_desc const* desc)
{
	qaws_sweep_boolean_report local;
	qaws_mesh_3d** sw = NULL;
	unsigned int i, ns = 0;
	qaws_intersect_desc id;
	qaws_status st = QAWS_STATUS_OK;
	if (!report)
		report = &local;
	memset(report, 0, sizeof *report);
	if (!out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out = NULL;
	if (!tool || !moves || !move_count || !stock ||
		(op != QAWS_BOOLEAN_DIFFERENCE && op != QAWS_BOOLEAN_UNION && op != QAWS_BOOLEAN_INTERSECTION))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < move_count; i++)
		if (!moves[i])
			return QAWS_STATUS_INVALID_ARGUMENT;
	sw = (qaws_mesh_3d**)qaws_internal_zalloc(move_count, sizeof *sw);
	if (!sw)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	/* each move's sweep (its attributes the tool's), all of them one operand */
	for (i = 0; i < move_count && st == QAWS_STATUS_OK; i++)
	{
		qaws_sweep_report r;
		st = qaws_mesh_3d_sweep(&sw[ns], &r, tool, moves[i], desc);
		report->slabs += r.slabs;
		if (st == QAWS_STATUS_OK && sw[ns])
		{
			ns++;
			report->moves++;
		}
	}
	if (st != QAWS_STATUS_OK)
		goto done;
	report->stock_volume = (qaws_scalar)swp_volume(stock);
	memset(&id, 0, sizeof id);
	id.executor = desc ? desc->executor : NULL;
	if (!ns)
	{
		/* (the tool swept nothing: the stock, or nothing, as the op says) */
		if (op == QAWS_BOOLEAN_INTERSECTION)
			goto done;
		st = cy_boolean(out, &report->boolean, &stock, 1, NULL, 0, QAWS_BOOLEAN_UNION, &id);
	}
	else
		st = cy_boolean(out, &report->boolean, &stock, 1, (qaws_mesh_3d const* const*)sw, ns, op, &id);   /* (sweep op stock: the stock first, taken from for the difference) */
	report->closed = report->boolean.closed;
	report->result_volume = (qaws_scalar)swp_volume(*out);
	report->changed_volume = (qaws_scalar)(op == QAWS_BOOLEAN_DIFFERENCE ? (double)report->stock_volume - (double)report->result_volume :
		op == QAWS_BOOLEAN_UNION ? (double)report->result_volume - (double)report->stock_volume : (double)report->result_volume);
done:
	for (i = 0; i < ns; i++)
		qaws_mesh_3d_destroy(sw[i]);
	qaws_internal_free(sw);
	return st;
}
