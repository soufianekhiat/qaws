/*
 * Swept solids (qaws_mesh_sweep.h): time cut into slabs by the chords'
 * deviation, each slab's front faces extruded, all united by the Boolean.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_sweep.h"
#include "../qaws_mesh_decimate.h"
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
