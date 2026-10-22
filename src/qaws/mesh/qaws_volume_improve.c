/*
 * Improvement of tetrahedral meshes (qaws_volume_improve.h): edge removal by
 * the ring's best triangulation, face removal, and smoothing up the worst
 * quality's gradient, on a working copy of the tetrahedra with their
 * adjacency; the border and the other cells untouched.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_volume_improve.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define VIM_RING 10u          /* the longest edge ring removed */
#define VIM_ACTIVE 8u         /* the qualities followed by a smoothing step */
#define VIM_EPS 1e-9          /* the least gain kept */

static unsigned char const vim_face[4][3] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };

typedef struct vim_mesh
{
	double* x;
	unsigned char* fixed;     /* per vertex: on the border or another cell */
	unsigned int* v2t;        /* per vertex: a live tetrahedron holding it */
	unsigned int nv;
	unsigned int *tv, *tn;    /* 4 per tetrahedron: vertices; the tetrahedron across facet i (HM_NONE: none) */
	unsigned char* dead;
	unsigned int *freel, nfree;
	unsigned int nt, cap;
} vim_mesh;

typedef struct vim_fkey
{
	unsigned int s[3], t, f;
} vim_fkey;

static int vim_fkey_cmp(void const* p, void const* q)
{
	vim_fkey const *a = (vim_fkey const*)p, *b = (vim_fkey const*)q;
	unsigned int i;
	for (i = 0; i < 3; i++)
		if (a->s[i] != b->s[i])
			return a->s[i] < b->s[i] ? -1 : 1;
	return 0;
}

static void vim_key(unsigned int const* v, unsigned int f, unsigned int* s)
{
	unsigned int t;
	s[0] = v[vim_face[f][0]];
	s[1] = v[vim_face[f][1]];
	s[2] = v[vim_face[f][2]];
	if (s[0] > s[1]) { t = s[0]; s[0] = s[1]; s[1] = t; }
	if (s[1] > s[2]) { t = s[1]; s[1] = s[2]; s[2] = t; }
	if (s[0] > s[1]) { t = s[0]; s[0] = s[1]; s[1] = t; }
}

/* the smallest sine of a tetrahedron's dihedral angles (sin = 3 V |e| / (2 A A'), negative when inverted) */
static double vim_quality(double const* a, double const* b, double const* c, double const* d)
{
	double const* p[4];
	double area[4], u[3], w[3], e[3], cr[3], vol, q = HUGE_VAL;
	unsigned int i, j, k;
	p[0] = a; p[1] = b; p[2] = c; p[3] = d;
	for (k = 0; k < 3; k++) { u[k] = b[k] - a[k]; w[k] = c[k] - a[k]; e[k] = d[k] - a[k]; }
	cr[0] = u[1] * w[2] - u[2] * w[1];
	cr[1] = u[2] * w[0] - u[0] * w[2];
	cr[2] = u[0] * w[1] - u[1] * w[0];
	vol = (cr[0] * e[0] + cr[1] * e[1] + cr[2] * e[2]) / 6;
	for (i = 0; i < 4; i++)
	{
		double const *x0 = p[vim_face[i][0]], *x1 = p[vim_face[i][1]], *x2 = p[vim_face[i][2]];
		for (k = 0; k < 3; k++) { u[k] = x1[k] - x0[k]; w[k] = x2[k] - x0[k]; }
		cr[0] = u[1] * w[2] - u[2] * w[1];
		cr[1] = u[2] * w[0] - u[0] * w[2];
		cr[2] = u[0] * w[1] - u[1] * w[0];
		area[i] = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
	}
	/* the edge between vertices i and j joins the facets opposite the other two */
	for (i = 0; i < 4; i++)
		for (j = i + 1; j < 4; j++)
		{
			unsigned int f0 = 4, f1 = 4;
			double l = 0;
			for (k = 0; k < 4; k++)
				if (k != i && k != j)
				{
					if (f0 == 4) f0 = k;
					else f1 = k;
				}
			for (k = 0; k < 3; k++) l += (p[j][k] - p[i][k]) * (p[j][k] - p[i][k]);
			q = fmin(q, area[f0] > 0 && area[f1] > 0 ? 1.5 * vol * sqrt(l) / (area[f0] * area[f1]) : -1);
		}
	return q;
}

static double vim_q(vim_mesh const* m, unsigned int const* v)
{
	return vim_quality(m->x + 3 * (size_t)v[0], m->x + 3 * (size_t)v[1], m->x + 3 * (size_t)v[2], m->x + 3 * (size_t)v[3]);
}

static void vim_free(vim_mesh* m)
{
	qaws_internal_free(m->x);
	qaws_internal_free(m->fixed);
	qaws_internal_free(m->v2t);
	qaws_internal_free(m->tv);
	qaws_internal_free(m->tn);
	qaws_internal_free(m->dead);
	qaws_internal_free(m->freel);
	memset(m, 0, sizeof *m);
}

/* room for k more tetrahedra; 0 on allocation failure (nothing changed) */
static int vim_reserve(vim_mesh* m, unsigned int k)
{
	unsigned int c;
	void *a, *b, *d, *f;
	if (m->nfree + (m->cap - m->nt) >= k)
		return 1;
	c = m->cap * 2 + k + 16;
	a = qaws_internal_realloc(m->tv, sizeof(unsigned int) * 4 * (size_t)c);
	if (a) m->tv = (unsigned int*)a;
	b = qaws_internal_realloc(m->tn, sizeof(unsigned int) * 4 * (size_t)c);
	if (b) m->tn = (unsigned int*)b;
	d = qaws_internal_realloc(m->dead, c);
	if (d) m->dead = (unsigned char*)d;
	f = qaws_internal_realloc(m->freel, sizeof(unsigned int) * (size_t)c);
	if (f) m->freel = (unsigned int*)f;
	if (!a || !b || !d || !f)
		return 0;
	m->cap = c;
	return 1;
}

/* The tetrahedra old (a group) replaced by nw (positive, filling the same region): the group's outer facets matched to
   the new tetrahedra's. 1 done; 0 when the facets do not match (nothing changed); -1 on allocation failure. */
static int vim_replace(vim_mesh* m, unsigned int const* old, unsigned int nold, unsigned int const (*nw)[4], unsigned int nnew)
{
	unsigned int ext[4 * VIM_RING][5], link[2 * VIM_RING][4], made[2 * VIM_RING], next = 0, i, f, k, j, g;
	for (i = 0; i < nold; i++)
		for (f = 0; f < 4; f++)
		{
			unsigned int u = m->tn[4 * (size_t)old[i] + f];
			int inner = 0;
			for (k = 0; k < nold; k++) inner |= u == old[k];
			if (inner)
				continue;
			vim_key(m->tv + 4 * (size_t)old[i], f, ext[next]);
			ext[next][3] = old[i];
			ext[next][4] = u;
			next++;
		}
	/* link[i][g]: a new tetrahedron j (j), or an outer facet k (0x80000000 | k) */
	for (i = 0; i < nnew; i++)
		for (g = 0; g < 4; g++)
		{
			unsigned int s[3], s2[3], h;
			link[i][g] = HM_NONE;
			vim_key(nw[i], g, s);
			for (j = 0; j < nnew && link[i][g] == HM_NONE; j++)
				for (h = 0; h < 4 && j != i; h++)
				{
					vim_key(nw[j], h, s2);
					if (s[0] == s2[0] && s[1] == s2[1] && s[2] == s2[2])
					{
						link[i][g] = j;
						break;
					}
				}
			for (k = 0; k < next && link[i][g] == HM_NONE; k++)
				if (s[0] == ext[k][0] && s[1] == ext[k][1] && s[2] == ext[k][2])
					link[i][g] = 0x80000000u | k;
			if (link[i][g] == HM_NONE)
				return 0;
		}
	if (!vim_reserve(m, nnew))
		return -1;
	for (i = 0; i < nnew; i++)
	{
		made[i] = m->nfree ? m->freel[--m->nfree] : m->nt++;
		memcpy(m->tv + 4 * (size_t)made[i], nw[i], sizeof(unsigned int) * 4);
		m->dead[made[i]] = 0;
		for (k = 0; k < 4; k++) m->v2t[nw[i][k]] = made[i];
	}
	for (i = 0; i < nnew; i++)
		for (g = 0; g < 4; g++)
		{
			unsigned int l = link[i][g], h;
			if (!(l & 0x80000000u))
			{
				m->tn[4 * (size_t)made[i] + g] = made[l];
				continue;
			}
			l &= 0x7FFFFFFFu;
			m->tn[4 * (size_t)made[i] + g] = ext[l][4];
			if (ext[l][4] != HM_NONE)
				for (h = 0; h < 4; h++)
					if (m->tn[4 * (size_t)ext[l][4] + h] == ext[l][3])
					{
						m->tn[4 * (size_t)ext[l][4] + h] = made[i];
						break;
					}
		}
	for (i = 0; i < nold; i++)
	{
		m->dead[old[i]] = 1;
		m->freel[m->nfree++] = old[i];
	}
	return 1;
}

/* The tetrahedra round edge {a, b} from t in turn and the ring's vertices (r[i] in cells i and i + 1), turned so that
   (a, b, r[0], r[1]) is positive; 0 when the ring is open (the edge on the border) or longer than VIM_RING. */
static unsigned int vim_ring(vim_mesh const* m, unsigned int t, unsigned int a, unsigned int b, unsigned int* cells, unsigned int* r)
{
	unsigned int n = 0, cur = t, prev = HM_NONE, i;
	for (;;)
	{
		unsigned int p[2], np = 0, k, nx = HM_NONE;
		if (n == VIM_RING)
			return 0;
		for (k = 0; k < 4; k++)
			if (m->tv[4 * (size_t)cur + k] != a && m->tv[4 * (size_t)cur + k] != b && np < 2)
				p[np++] = k;
		if (np != 2)
			return 0;
		for (k = 0; k < 2; k++)
		{
			unsigned int u = m->tn[4 * (size_t)cur + p[k]];
			if (u == HM_NONE)
				return 0;
			if (nx == HM_NONE && u != prev)
			{
				nx = u;
				r[n] = m->tv[4 * (size_t)cur + p[1 - k]];
			}
		}
		if (nx == HM_NONE)
			return 0;
		cells[n++] = cur;
		prev = cur;
		cur = nx;
		if (cur == t)
			break;
	}
	if (n < 3)
		return 0;
	{
		double const *xa = m->x + 3 * (size_t)a, *xb = m->x + 3 * (size_t)b, *x0 = m->x + 3 * (size_t)r[0], *x1 = m->x + 3 * (size_t)r[1];
		double u[3], v[3], w[3];
		for (i = 0; i < 3; i++) { u[i] = xb[i] - xa[i]; v[i] = x0[i] - xa[i]; w[i] = x1[i] - xa[i]; }
		if (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]) < 0)
			for (i = 0; i < n / 2; i++)
			{
				unsigned int s = r[i];
				r[i] = r[n - 1 - i];
				r[n - 1 - i] = s;
			}
	}
	return n;
}

/* The ring polygon's triangulation whose worst tetrahedron (each triangle joined to a and to b) is best (Shewchuk's
   dynamic program): that worst quality; the 2 (n - 2) tetrahedra in out when not NULL. */
static double vim_ring_best(vim_mesh const* m, unsigned int a, unsigned int b, unsigned int const* r, unsigned int n, unsigned int (*out)[4])
{
	double best[VIM_RING][VIM_RING];
	unsigned char arg[VIM_RING][VIM_RING];
	unsigned int len, i, j, k;
	for (i = 0; i + 1 < n; i++) best[i][i + 1] = HUGE_VAL;
	for (len = 2; len < n; len++)
		for (i = 0; i + len < n; i++)
		{
			j = i + len;
			best[i][j] = -HUGE_VAL;
			arg[i][j] = (unsigned char)(i + 1);
			for (k = i + 1; k < j; k++)
			{
				double q = fmin(best[i][k], best[k][j]);
				if (q > best[i][j])
				{
					unsigned int ta[4], tb[4];
					ta[0] = r[i]; ta[1] = r[j]; ta[2] = r[k]; ta[3] = a;
					tb[0] = r[i]; tb[1] = r[k]; tb[2] = r[j]; tb[3] = b;
					q = fmin(q, fmin(vim_q(m, ta), vim_q(m, tb)));
					if (q > best[i][j])
					{
						best[i][j] = q;
						arg[i][j] = (unsigned char)k;
					}
				}
			}
		}
	if (out)
	{
		unsigned int st[2 * VIM_RING][2], ns = 0, no = 0;
		st[ns][0] = 0;
		st[ns++][1] = n - 1;
		while (ns)
		{
			--ns;
			i = st[ns][0];
			j = st[ns][1];
			if (j - i < 2)
				continue;
			k = arg[i][j];
			out[no][0] = r[i]; out[no][1] = r[j]; out[no][2] = r[k]; out[no][3] = a; no++;
			out[no][0] = r[i]; out[no][1] = r[k]; out[no][2] = r[j]; out[no][3] = b; no++;
			st[ns][0] = i; st[ns++][1] = k;
			st[ns][0] = k; st[ns++][1] = j;
		}
	}
	return best[0][n - 1];
}

/* The best edge or face removal at tetrahedron t, applied when it raises the worst quality of what it replaces:
   1 (edge), 2 (face), 0 none, -1 on allocation failure. */
static int vim_flip(vim_mesh* m, unsigned int t)
{
	unsigned int v[4], i, j, kind = 0, ba = 0, bb = 0;
	double gain = VIM_EPS;
	memcpy(v, m->tv + 4 * (size_t)t, sizeof v);
	for (i = 0; i < 4; i++)
		for (j = i + 1; j < 4; j++)
		{
			unsigned int cells[VIM_RING], r[VIM_RING], n = vim_ring(m, t, v[i], v[j], cells, r), k;
			double qo = HUGE_VAL, qn;
			if (!n)
				continue;
			for (k = 0; k < n; k++) qo = fmin(qo, vim_q(m, m->tv + 4 * (size_t)cells[k]));
			qn = vim_ring_best(m, v[i], v[j], r, n, NULL);
			if (qn > 0 && qn - qo > gain)
			{
				gain = qn - qo;
				kind = 1;
				ba = v[i];
				bb = v[j];
			}
		}
	for (i = 0; i < 4; i++)
	{
		unsigned int u = m->tn[4 * (size_t)t + i], e = HM_NONE, k, nw[3][4];
		double qn;
		if (u == HM_NONE)
			continue;
		for (k = 0; k < 4; k++)
		{
			unsigned int w = m->tv[4 * (size_t)u + k];
			if (w != v[vim_face[i][0]] && w != v[vim_face[i][1]] && w != v[vim_face[i][2]])
				e = w;
		}
		for (k = 0; k < 3; k++)
		{
			nw[k][0] = v[vim_face[i][k]];
			nw[k][1] = v[vim_face[i][(k + 1) % 3]];
			nw[k][2] = v[i];
			nw[k][3] = e;
		}
		qn = fmin(vim_q(m, nw[0]), fmin(vim_q(m, nw[1]), vim_q(m, nw[2])));
		if (qn > 0 && qn - fmin(vim_q(m, v), vim_q(m, m->tv + 4 * (size_t)u)) > gain)
		{
			gain = qn - fmin(vim_q(m, v), vim_q(m, m->tv + 4 * (size_t)u));
			kind = 2;
			ba = i;
		}
	}
	if (kind == 1)
	{
		unsigned int cells[VIM_RING], r[VIM_RING], nw[2 * VIM_RING][4], n = vim_ring(m, t, ba, bb, cells, r);
		vim_ring_best(m, ba, bb, r, n, nw);
		return vim_replace(m, cells, n, (unsigned int const(*)[4])nw, 2 * (n - 2));
	}
	if (kind == 2)
	{
		unsigned int u = m->tn[4 * (size_t)t + ba], e = HM_NONE, k, nw[3][4], old[2];
		int res;
		for (k = 0; k < 4; k++)
		{
			unsigned int w = m->tv[4 * (size_t)u + k];
			if (w != v[vim_face[ba][0]] && w != v[vim_face[ba][1]] && w != v[vim_face[ba][2]])
				e = w;
		}
		for (k = 0; k < 3; k++)
		{
			nw[k][0] = v[vim_face[ba][k]];
			nw[k][1] = v[vim_face[ba][(k + 1) % 3]];
			nw[k][2] = v[ba];
			nw[k][3] = e;
		}
		old[0] = t;
		old[1] = u;
		res = vim_replace(m, old, 2, (unsigned int const(*)[4])nw, 3);
		return res == 1 ? 2 : res;
	}
	return 0;
}

/* the worst quality of the tetrahedra cl[0 .. n) */
static double vim_worst(vim_mesh const* m, unsigned int const* cl, unsigned int n)
{
	double q = HUGE_VAL;
	unsigned int i;
	for (i = 0; i < n; i++) q = fmin(q, vim_q(m, m->tv + 4 * (size_t)cl[i]));
	return q;
}

/* Vertex v moved to raise the worst quality of its tetrahedra cl: its neighbours' centre (and halfway there) tried,
   then steps up the nearest point to zero of the active gradients' hull, each kept only when the worst rises. 1 moved. */
static int vim_smooth_vertex(vim_mesh* m, unsigned int v, unsigned int const* cl, unsigned int n)
{
	double* x = m->x + 3 * (size_t)v;
	double keep[3], cen[3] = { 0, 0, 0 }, q0 = vim_worst(m, cl, n), L = 0, q;
	unsigned int i, k, it, cnt = 0;
	int moved = 0;
	memcpy(keep, x, sizeof keep);
	for (i = 0; i < n; i++)
		for (k = 0; k < 4; k++)
		{
			unsigned int w = m->tv[4 * (size_t)cl[i] + k];
			if (w == v)
				continue;
			cen[0] += m->x[3 * (size_t)w];
			cen[1] += m->x[3 * (size_t)w + 1];
			cen[2] += m->x[3 * (size_t)w + 2];
			L += sqrt((m->x[3 * (size_t)w] - x[0]) * (m->x[3 * (size_t)w] - x[0]) + (m->x[3 * (size_t)w + 1] - x[1]) * (m->x[3 * (size_t)w + 1] - x[1]) +
				(m->x[3 * (size_t)w + 2] - x[2]) * (m->x[3 * (size_t)w + 2] - x[2]));
			cnt++;
		}
	if (!cnt)
		return 0;
	L /= cnt;
	for (k = 0; k < 2; k++)
	{
		double t = k == 0 ? 1 : 0.5;
		unsigned int d;
		for (d = 0; d < 3; d++) x[d] = keep[d] + t * (cen[d] / cnt - keep[d]);
		q = vim_worst(m, cl, n);
		if (q > q0 + VIM_EPS)
		{
			q0 = q;
			memcpy(keep, x, sizeof keep);
			moved = 1;
			break;
		}
		memcpy(x, keep, sizeof keep);
	}
	for (it = 0; it < 8; it++)
	{
		double g[VIM_ACTIVE][3], dir[3], h = 1e-7 * L, nd, s;
		unsigned int na = 0, d, j;
		for (i = 0; i < n && na < VIM_ACTIVE; i++)
		{
			unsigned int const* t = m->tv + 4 * (size_t)cl[i];
			if (vim_q(m, t) > q0 + 1e-3)
				continue;
			for (d = 0; d < 3; d++)
			{
				double qp, qm;
				x[d] = keep[d] + h;
				qp = vim_q(m, t);
				x[d] = keep[d] - h;
				qm = vim_q(m, t);
				x[d] = keep[d];
				g[na][d] = (qp - qm) / (2 * h);
			}
			na++;
		}
		if (!na)
			break;
		/* the hull's point nearest zero (Frank-Wolfe) */
		memcpy(dir, g[0], sizeof dir);
		for (j = 0; j < 40; j++)
		{
			double best = HUGE_VAL, dd = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2], e[3], ee, gam;
			unsigned int bi = 0;
			for (i = 0; i < na; i++)
			{
				double p = g[i][0] * dir[0] + g[i][1] * dir[1] + g[i][2] * dir[2];
				if (p < best) { best = p; bi = i; }
			}
			if (dd - best <= 1e-14 * dd)
				break;
			for (d = 0; d < 3; d++) e[d] = dir[d] - g[bi][d];
			ee = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
			gam = ee > 0 ? (e[0] * dir[0] + e[1] * dir[1] + e[2] * dir[2]) / ee : 0;
			gam = gam < 0 ? 0 : (gam > 1 ? 1 : gam);
			for (d = 0; d < 3; d++) dir[d] -= gam * e[d];
		}
		nd = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
		if (nd < 1e-12)
			break;
		for (s = 0.25 * L; s > 1e-4 * L; s *= 0.5)
		{
			for (d = 0; d < 3; d++) x[d] = keep[d] + s * dir[d] / nd;
			q = vim_worst(m, cl, n);
			if (q > q0 + VIM_EPS)
				break;
		}
		if (s <= 1e-4 * L)
		{
			memcpy(x, keep, sizeof keep);
			break;
		}
		q0 = q;
		memcpy(keep, x, sizeof keep);
		moved = 1;
	}
	return moved;
}

/* one smoothing sweep over the free vertices of tetrahedra below goal: the moves made, or -1 on allocation failure */
static int vim_smooth(vim_mesh* m, double goal)
{
	unsigned int *v0 = (unsigned int*)qaws_internal_scratch_zalloc((size_t)m->nv + 2, sizeof(unsigned int));
	unsigned int *vt = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * (4 * (size_t)m->nt + 1));
	unsigned int* fill = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)m->nv + 1));
	unsigned int t, k, v;
	int moves = 0;
	if (!v0 || !vt || !fill)
	{
		qaws_internal_free(v0);
		qaws_internal_free(vt);
		qaws_internal_free(fill);
		return -1;
	}
	for (t = 0; t < m->nt; t++)
		for (k = 0; !m->dead[t] && k < 4; k++) v0[m->tv[4 * (size_t)t + k] + 1]++;
	for (v = 0; v < m->nv; v++) v0[v + 1] += v0[v];
	memcpy(fill, v0, sizeof(unsigned int) * m->nv);
	for (t = 0; t < m->nt; t++)
		for (k = 0; !m->dead[t] && k < 4; k++) vt[fill[m->tv[4 * (size_t)t + k]]++] = t;
	for (v = 0; v < m->nv; v++)
	{
		if (m->fixed[v] || v0[v] == v0[v + 1] || vim_worst(m, vt + v0[v], v0[v + 1] - v0[v]) >= goal)
			continue;
		moves += vim_smooth_vertex(m, v, vt + v0[v], v0[v + 1] - v0[v]);
	}
	qaws_internal_free(v0);
	qaws_internal_free(vt);
	qaws_internal_free(fill);
	return moves;
}

static double vim_min_quality(vim_mesh const* m)
{
	double q = HUGE_VAL;
	unsigned int t;
	for (t = 0; t < m->nt; t++)
		if (!m->dead[t])
			q = fmin(q, vim_q(m, m->tv + 4 * (size_t)t));
	return q;
}

/* the working copy of a mesh's tetrahedra: adjacency by facets' vertices, the border's and the other cells' vertices fixed,
   each vertex's cell; INVALID_ARGUMENT without tetrahedra */
static qaws_status vim_load(vim_mesh* m, qaws_volume_mesh const* mesh)
{
	vim_fkey* keys;
	unsigned int nt = 0, c, i, k, nk = 0;
	memset(m, 0, sizeof *m);
	for (c = 0; c < mesh->nc; c++) nt += mesh->type[c] == QAWS_CELL_TETRAHEDRON;
	if (!nt)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m->nv = mesh->nv;
	m->cap = 2 * nt + 16;
	m->x = (double*)qaws_internal_malloc(sizeof(double) * 3 * ((size_t)m->nv + 1));
	m->fixed = (unsigned char*)qaws_internal_zalloc((size_t)m->nv + 1, 1);
	m->v2t = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)m->nv + 1));
	m->tv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * 4 * (size_t)m->cap);
	m->tn = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * 4 * (size_t)m->cap);
	m->dead = (unsigned char*)qaws_internal_zalloc(m->cap, 1);
	m->freel = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (size_t)m->cap);
	keys = (vim_fkey*)qaws_internal_scratch_malloc(sizeof(vim_fkey) * 4 * (size_t)nt);
	if (!m->x || !m->fixed || !m->v2t || !m->tv || !m->tn || !m->dead || !m->freel || !keys)
	{
		qaws_internal_free(keys);
		vim_free(m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	memcpy(m->x, mesh->x, sizeof(double) * 3 * (size_t)m->nv);
	for (i = 0; i < m->nv; i++) m->v2t[i] = HM_NONE;
	for (c = 0; c < mesh->nc; c++)
	{
		unsigned int const* v = mesh->cv + mesh->c0[c];
		if (mesh->type[c] != QAWS_CELL_TETRAHEDRON)
		{
			for (k = mesh->c0[c]; k < mesh->c0[c + 1]; k++) m->fixed[mesh->cv[k]] = 1;
			continue;
		}
		memcpy(m->tv + 4 * (size_t)m->nt, v, sizeof(unsigned int) * 4);
		for (k = 0; k < 4; k++)
		{
			vim_key(v, k, keys[nk].s);
			keys[nk].t = m->nt;
			keys[nk++].f = k;
			m->tn[4 * (size_t)m->nt + k] = HM_NONE;
			m->v2t[v[k]] = m->nt;
		}
		m->nt++;
	}
	/* facets with the same vertices paired; the unpaired ones the border, their vertices fixed */
	qsort(keys, nk, sizeof(vim_fkey), vim_fkey_cmp);
	for (i = 0; i < nk; i++)
	{
		if (i + 1 < nk && !vim_fkey_cmp(&keys[i], &keys[i + 1]))
		{
			m->tn[4 * (size_t)keys[i].t + keys[i].f] = keys[i + 1].t;
			m->tn[4 * (size_t)keys[i + 1].t + keys[i + 1].f] = keys[i].t;
			i++;
			continue;
		}
		for (k = 0; k < 3; k++) m->fixed[keys[i].s[k]] = 1;
	}
	qaws_internal_free(keys);
	return QAWS_STATUS_OK;
}

/* the result: the input's other cells first, in their order, then the live tetrahedra */
static qaws_status vim_store(qaws_volume_mesh** out_mesh, vim_mesh const* m, qaws_volume_mesh const* mesh)
{
	unsigned int* cv = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)mesh->c0[mesh->nc] + 4 * (size_t)m->nt + 1));
	unsigned char* type = (unsigned char*)qaws_internal_scratch_malloc((size_t)mesh->nc + m->nt + 1);
	unsigned int c, k, ncv = 0, nc = 0;
	qaws_status st;
	if (!cv || !type)
	{
		qaws_internal_free(cv);
		qaws_internal_free(type);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (c = 0; c < mesh->nc; c++)
		if (mesh->type[c] != QAWS_CELL_TETRAHEDRON)
		{
			for (k = mesh->c0[c]; k < mesh->c0[c + 1]; k++) cv[ncv++] = mesh->cv[k];
			type[nc++] = mesh->type[c];
		}
	for (c = 0; c < m->nt; c++)
		if (!m->dead[c])
		{
			memcpy(cv + ncv, m->tv + 4 * (size_t)c, sizeof(unsigned int) * 4);
			ncv += 4;
			type[nc++] = QAWS_CELL_TETRAHEDRON;
		}
	st = vm_create(m->x, m->nv, cv, type, nc, out_mesh);
	qaws_internal_free(cv);
	qaws_internal_free(type);
	return st;
}

qaws_status qaws_volume_mesh_improve(qaws_volume_mesh** out_mesh, qaws_volume_improve_report* report, qaws_volume_mesh const* mesh,
	qaws_volume_improve_desc const* desc)
{
	vim_mesh m;
	unsigned int c, it, iters = desc && desc->iterations ? desc->iterations : 10;
	double goal = sin((desc && desc->target_dihedral > 0 ? (double)desc->target_dihedral : 40.0) * 3.14159265358979323846 / 180);
	qaws_volume_improve_report rep;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&rep, 0, sizeof rep);
	st = vim_load(&m, mesh);
	if (st != QAWS_STATUS_OK)
		return st;
	rep.min_quality_before = (qaws_scalar)vim_min_quality(&m);
	for (it = 0; it < iters; it++)
	{
		unsigned int changed = 0, snap = m.nt;
		double q0 = vim_min_quality(&m);
		int stalled = 1;
		if (!(desc && desc->no_smoothing))
		{
			int r = vim_smooth(&m, goal);
			if (r < 0)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto done;
			}
			rep.moves += (unsigned int)r;
			changed += (unsigned int)r;
			stalled = vim_min_quality(&m) < q0 + 1e-3;
		}
		/* flips where smoothing stalls (smoothing first: a flip can leave cells of border vertices only, beyond its reach) */
		if (stalled && !(desc && desc->no_flips))
			for (c = 0; c < snap; c++)
			{
				int r;
				if (m.dead[c] || vim_q(&m, m.tv + 4 * (size_t)c) >= goal)
					continue;
				r = vim_flip(&m, c);
				if (r < 0)
				{
					st = QAWS_STATUS_ALLOCATION_FAILURE;
					goto done;
				}
				rep.edge_removals += r == 1;
				rep.face_removals += r == 2;
				changed += r > 0;
			}
		rep.iterations = it + 1;
		if (!changed)
			break;
	}
	rep.min_quality_after = (qaws_scalar)vim_min_quality(&m);
	st = vim_store(out_mesh, &m, mesh);
	if (st == QAWS_STATUS_OK)
		qaws_volume_mesh_measure_quality(&rep.quality, *out_mesh);
done:
	if (report)
		*report = rep;
	vim_free(&m);
	return st;
}

/* ======================================================================== */
/*  Border points moved inside                                              */
/* ======================================================================== */

/* the cells round vertex v from cell t0 (holding it), across their facets through v: their count, 0 past cap */
static unsigned int vim_star(vim_mesh const* m, unsigned int v, unsigned int t0, unsigned int* out, unsigned int cap)
{
	unsigned int n = 1, i, k, j;
	out[0] = t0;
	for (i = 0; i < n; i++)
		for (k = 0; k < 4; k++)
		{
			unsigned int u = m->tn[4 * (size_t)out[i] + k];
			if (m->tv[4 * (size_t)out[i] + k] == v || u == HM_NONE)
				continue;   /* (facet k, opposite vertex k, holds v unless vertex k is v) */
			for (j = 0; j < n && out[j] != u; j++) {}
			if (j < n)
				continue;
			if (n == cap)
				return 0;
			out[n++] = u;
		}
	return n;
}

/* The tetrahedra old replaced by nw (positive, filling the same region): new facets matched to each other first, then to
   the old group's outer facets (their cells across, or none on the border); the nb left unmatched are the new border. 1 done; 0 when the facets do not match (nothing changed); -1 on
   allocation failure. */
static int vim_replace_open(vim_mesh* m, unsigned int const* old, unsigned int nold, unsigned int const (*nw)[4], unsigned int nnew, unsigned int nb)
{
	unsigned int (*ext)[5] = (unsigned int(*)[5])qaws_internal_scratch_malloc(sizeof(unsigned int) * 5 * 4 * ((size_t)nold + 1));
	unsigned int (*link)[4] = (unsigned int(*)[4])qaws_internal_scratch_malloc(sizeof(unsigned int) * 4 * ((size_t)nnew + 1));
	unsigned int* made = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nnew + 1));
	unsigned int next = 0, open = 0, i, f, k, j, g;
	int res = 1;
	if (!ext || !link || !made)
	{
		res = -1;
		goto out;
	}
	for (i = 0; i < nold; i++)
		for (f = 0; f < 4; f++)
		{
			unsigned int u = m->tn[4 * (size_t)old[i] + f];
			int inner = 0;
			for (k = 0; k < nold; k++) inner |= u == old[k];
			if (inner)
				continue;
			vim_key(m->tv + 4 * (size_t)old[i], f, ext[next]);
			ext[next][3] = old[i];
			ext[next][4] = u;
			next++;
		}
	for (i = 0; i < nnew; i++)
		for (g = 0; g < 4; g++)
		{
			unsigned int s[3], s2[3], h;
			link[i][g] = HM_NONE;
			vim_key(nw[i], g, s);
			for (j = 0; j < nnew && link[i][g] == HM_NONE; j++)
				for (h = 0; h < 4 && j != i; h++)
				{
					vim_key(nw[j], h, s2);
					if (s[0] == s2[0] && s[1] == s2[1] && s[2] == s2[2])
					{
						link[i][g] = j;
						break;
					}
				}
			for (k = 0; k < next && link[i][g] == HM_NONE; k++)
				if (s[0] == ext[k][0] && s[1] == ext[k][1] && s[2] == ext[k][2])
					link[i][g] = 0x80000000u | k;
			open += link[i][g] == HM_NONE;
		}
	if (open != nb)
	{
		res = 0;
		goto out;
	}
	if (!vim_reserve(m, nnew))
	{
		res = -1;
		goto out;
	}
	for (i = 0; i < nnew; i++)
	{
		made[i] = m->nfree ? m->freel[--m->nfree] : m->nt++;
		memcpy(m->tv + 4 * (size_t)made[i], nw[i], sizeof(unsigned int) * 4);
		m->dead[made[i]] = 0;
		for (k = 0; k < 4; k++) m->v2t[nw[i][k]] = made[i];
	}
	for (i = 0; i < nnew; i++)
		for (g = 0; g < 4; g++)
		{
			unsigned int l = link[i][g], h;
			if (l == HM_NONE)
			{
				m->tn[4 * (size_t)made[i] + g] = HM_NONE;
				continue;
			}
			if (!(l & 0x80000000u))
			{
				m->tn[4 * (size_t)made[i] + g] = made[l];
				continue;
			}
			l &= 0x7FFFFFFFu;
			m->tn[4 * (size_t)made[i] + g] = ext[l][4];
			for (h = 0; h < 4 && ext[l][4] != HM_NONE; h++)
				if (m->tn[4 * (size_t)ext[l][4] + h] == ext[l][3])
				{
					m->tn[4 * (size_t)ext[l][4] + h] = made[i];
					break;
				}
		}
	for (i = 0; i < nold; i++)
	{
		m->dead[old[i]] = 1;
		m->freel[m->nfree++] = old[i];
	}
out:
	qaws_internal_free(ext);
	qaws_internal_free(link);
	qaws_internal_free(made);
	return res;
}

/* A planar polygon (ids p[0 .. n), counter-clockwise seen along -nrm) cut into triangles by ear clipping, appended to tri
   (3 each, the same turn); their count, 0 when it gets stuck */
static unsigned int vim_ears(vim_mesh const* m, unsigned int* p, unsigned int n, double const* nrm, unsigned int* tri)
{
	double e1[3], e2[3], l = 0, (*q)[2], scale = 0;
	unsigned int i, k, nt = 0;
	q = (double(*)[2])qaws_internal_scratch_malloc(sizeof(double) * 2 * ((size_t)n + 1));
	if (!q)
		return 0;
	for (k = 0; k < 3; k++) e1[k] = m->x[3 * (size_t)p[1] + k] - m->x[3 * (size_t)p[0] + k];
	l = e1[0] * nrm[0] + e1[1] * nrm[1] + e1[2] * nrm[2];
	for (k = 0; k < 3; k++) e1[k] -= l * nrm[k];
	l = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
	if (!(l > 0))
	{
		qaws_internal_free(q);
		return 0;
	}
	for (k = 0; k < 3; k++) e1[k] /= l;
	e2[0] = nrm[1] * e1[2] - nrm[2] * e1[1];
	e2[1] = nrm[2] * e1[0] - nrm[0] * e1[2];
	e2[2] = nrm[0] * e1[1] - nrm[1] * e1[0];
	for (i = 0; i < n; i++)
	{
		double const* x = m->x + 3 * (size_t)p[i];
		q[i][0] = x[0] * e1[0] + x[1] * e1[1] + x[2] * e1[2];
		q[i][1] = x[0] * e2[0] + x[1] * e2[1] + x[2] * e2[2];
		scale = fmax(scale, fabs(q[i][0] - q[0][0]) + fabs(q[i][1] - q[0][1]));
	}
	while (n > 3)
	{
		unsigned int cut = HM_NONE;
		for (i = 0; i < n && cut == HM_NONE; i++)
		{
			unsigned int a = (i + n - 1) % n, b = (i + 1) % n, j;
			double cr = (q[i][0] - q[a][0]) * (q[b][1] - q[a][1]) - (q[i][1] - q[a][1]) * (q[b][0] - q[a][0]);
			int clear = cr > 1e-12 * scale * scale;
			for (j = 0; clear && j < n; j++)
			{
				double d0, d1, d2;
				if (j == a || j == i || j == b)
					continue;
				d0 = (q[i][0] - q[a][0]) * (q[j][1] - q[a][1]) - (q[i][1] - q[a][1]) * (q[j][0] - q[a][0]);
				d1 = (q[b][0] - q[i][0]) * (q[j][1] - q[i][1]) - (q[b][1] - q[i][1]) * (q[j][0] - q[i][0]);
				d2 = (q[a][0] - q[b][0]) * (q[j][1] - q[b][1]) - (q[a][1] - q[b][1]) * (q[j][0] - q[b][0]);
				clear = !(d0 >= 0 && d1 >= 0 && d2 >= 0);
			}
			if (clear)
				cut = i;
		}
		if (cut == HM_NONE)
		{
			qaws_internal_free(q);
			return 0;
		}
		tri[3 * nt] = p[(cut + n - 1) % n];
		tri[3 * nt + 1] = p[cut];
		tri[3 * nt + 2] = p[(cut + 1) % n];
		nt++;
		for (i = cut; i + 1 < n; i++)
		{
			p[i] = p[i + 1];
			q[i][0] = q[i + 1][0];
			q[i][1] = q[i + 1][1];
		}
		n--;
	}
	if ((q[1][0] - q[0][0]) * (q[2][1] - q[0][1]) - (q[1][1] - q[0][1]) * (q[2][0] - q[0][0]) <= 1e-12 * scale * scale)
	{
		qaws_internal_free(q);
		return 0;
	}
	tri[3 * nt] = p[0];
	tri[3 * nt + 1] = p[1];
	tri[3 * nt + 2] = p[2];
	qaws_internal_free(q);
	return nt + 1;
}

#define VIM_STAR 256u   /* the largest star of a border point moved inside */

/* Border point v moved inside: its border triangles merged (the link re-triangulated in its plane, or in each of its two
   planes along a crease), its star coned from the point inside where the worst cell is best. 1 done, 0 not, -1 memory. */
static int vim_unsplit_point(vim_mesh* m, unsigned int v)
{
	unsigned int star[VIM_STAR], ns, i, k, ne = 0, g, ng = 0, nb = 0, nn = 0, cnt = 0;
	unsigned int ea[VIM_STAR][2], eg[VIM_STAR], poly[VIM_STAR + 2], tri[3 * VIM_STAR];
	double nrm[2][3], cen[3] = { 0, 0, 0 }, keep[3], best = -HUGE_VAL, bt = 0;
	unsigned int (*nw)[4];
	int res;
	static double const ts[5] = { 1, 0.7, 0.5, 0.3, 0.15 };
	if (m->v2t[v] == HM_NONE || m->dead[m->v2t[v]])
		return 0;
	ns = vim_star(m, v, m->v2t[v], star, VIM_STAR);
	if (!ns)
		return 0;
	/* the border triangles at v: their sides opposite v (turned as the triangle), grouped by plane */
	for (i = 0; i < ns; i++)
		for (k = 0; k < 4; k++)
		{
			unsigned int const* t = m->tv + 4 * (size_t)star[i];
			unsigned int a, b, c, r;
			double u[3], w[3], n[3], l;
			if (t[k] == v || m->tn[4 * (size_t)star[i] + k] != HM_NONE)
				continue;
			a = t[vim_face[k][0]];
			b = t[vim_face[k][1]];
			c = t[vim_face[k][2]];
			while (a != v) { r = a; a = b; b = c; c = r; }
			for (r = 0; r < 3; r++)
			{
				u[r] = m->x[3 * (size_t)b + r] - m->x[3 * (size_t)a + r];
				w[r] = m->x[3 * (size_t)c + r] - m->x[3 * (size_t)a + r];
			}
			n[0] = u[1] * w[2] - u[2] * w[1];
			n[1] = u[2] * w[0] - u[0] * w[2];
			n[2] = u[0] * w[1] - u[1] * w[0];
			l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
			if (!(l > 0) || ne == VIM_STAR)
				return 0;
			for (r = 0; r < 3; r++) n[r] /= l;
			for (g = 0; g < ng && nrm[g][0] * n[0] + nrm[g][1] * n[1] + nrm[g][2] * n[2] < 1 - 1e-9; g++) {}
			if (g == ng)
			{
				if (ng == 2)
					return 0;
				memcpy(nrm[ng++], n, sizeof n);
			}
			ea[ne][0] = b;
			ea[ne][1] = c;
			eg[ne++] = g;
		}
	if (!ne)
		return 0;
	/* each plane's polygon: its sides chained (a cycle, or a path closed through v's place along the crease) */
	for (g = 0; g < ng; g++)
	{
		unsigned int start = HM_NONE, cur, np = 0, ming = 0, e;
		for (e = 0; e < ne; e++)
			if (eg[e] == g)
			{
				unsigned int f, in = 0;
				ming++;
				for (f = 0; f < ne; f++) in |= eg[f] == g && ea[f][1] == ea[e][0];
				if (!in)
					start = ea[e][0];
			}
		if ((ng == 1) != (start == HM_NONE))
			return 0;
		if (start == HM_NONE)
			for (e = 0; e < ne && start == HM_NONE; e++)
				if (eg[e] == g)
					start = ea[e][0];
		cur = start;
		for (;;)
		{
			unsigned int f;
			if (np > ming + 1)
				return 0;
			poly[np++] = cur;
			for (f = 0; f < ne && !(eg[f] == g && ea[f][0] == cur); f++) {}
			if (f == ne)
				break;
			cur = ea[f][1];
			if (cur == start)
				break;
		}
		if (np != ming + (ng == 1 ? 0u : 1u))
			return 0;
		if (ng == 2)
		{
			/* v on the closing side, between its ends */
			double const *s = m->x + 3 * (size_t)poly[0], *t = m->x + 3 * (size_t)poly[np - 1], *x = m->x + 3 * (size_t)v;
			double d[3], p[3], cr[3], dd, pd;
			for (k = 0; k < 3; k++) { d[k] = t[k] - s[k]; p[k] = x[k] - s[k]; }
			cr[0] = d[1] * p[2] - d[2] * p[1];
			cr[1] = d[2] * p[0] - d[0] * p[2];
			cr[2] = d[0] * p[1] - d[1] * p[0];
			dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
			pd = p[0] * d[0] + p[1] * d[1] + p[2] * d[2];
			if (!(dd > 0) || cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2] > 1e-18 * dd * dd || pd <= 0 || pd >= dd)
				return 0;
		}
		k = vim_ears(m, poly, np, nrm[g], tri + 3 * nb);
		if (!k)
			return 0;
		nb += k;
	}
	/* the new cells: the star's with v at its new place, and each new border triangle coned from it */
	nw = (unsigned int(*)[4])qaws_internal_scratch_malloc(sizeof(unsigned int) * 4 * ((size_t)ns + nb + 1));
	if (!nw)
		return -1;
	for (i = 0; i < ns; i++) memcpy(nw[nn++], m->tv + 4 * (size_t)star[i], sizeof(unsigned int) * 4);
	for (i = 0; i < nb; i++)
	{
		nw[nn][0] = tri[3 * i];
		nw[nn][1] = tri[3 * i + 2];
		nw[nn][2] = tri[3 * i + 1];
		nw[nn++][3] = v;
	}
	for (i = 0; i < ns; i++)
		for (k = 0; k < 4; k++)
		{
			unsigned int w = m->tv[4 * (size_t)star[i] + k], r;
			if (w == v)
				continue;
			for (r = 0; r < 3; r++) cen[r] += m->x[3 * (size_t)w + r];
			cnt++;
		}
	for (k = 0; k < 3; k++) cen[k] /= cnt;
	memcpy(keep, m->x + 3 * (size_t)v, sizeof keep);
	for (i = 0; i < 5; i++)
	{
		double q = HUGE_VAL;
		unsigned int j;
		for (k = 0; k < 3; k++) m->x[3 * (size_t)v + k] = keep[k] + ts[i] * (cen[k] - keep[k]);
		for (j = 0; j < nn && q > best; j++) q = fmin(q, vim_q(m, nw[j]));
		if (q > best)
		{
			best = q;
			bt = ts[i];
		}
	}
	if (!(best > 1e-3))
	{
		memcpy(m->x + 3 * (size_t)v, keep, sizeof keep);
		qaws_internal_free(nw);
		return 0;
	}
	for (k = 0; k < 3; k++) m->x[3 * (size_t)v + k] = keep[k] + bt * (cen[k] - keep[k]);
	res = vim_replace_open(m, star, ns, (unsigned int const(*)[4])nw, nn, nb);
	qaws_internal_free(nw);
	if (res != 1)
		memcpy(m->x + 3 * (size_t)v, keep, sizeof keep);
	if (res == 1)
		m->fixed[v] = 0;
	return res;
}

qaws_status vim_unsplit_border(qaws_volume_mesh** out_mesh, unsigned int* out_moved, unsigned int* out_left, qaws_volume_mesh const* mesh,
	unsigned char const* steiner)
{
	vim_mesh m;
	unsigned int v, pass, moved = 0, left = 0;
	unsigned char* done;
	qaws_status st;
	*out_mesh = NULL;
	st = vim_load(&m, mesh);
	if (st != QAWS_STATUS_OK)
		return st;
	done = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m.nv + 1, 1);
	if (!done)
	{
		vim_free(&m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	/* (one point moved can free its neighbours along a crease: a few passes) */
	for (pass = 0; pass < 4; pass++)
	{
		unsigned int now = 0;
		for (v = 0; v < m.nv; v++)
		{
			int r;
			if (!steiner[v] || done[v])
				continue;
			r = vim_unsplit_point(&m, v);
			if (r < 0)
			{
				st = QAWS_STATUS_ALLOCATION_FAILURE;
				goto out;
			}
			if (r == 1)
			{
				done[v] = 1;
				moved++;
				now++;
			}
		}
		if (!now)
			break;
	}
	for (v = 0; v < m.nv; v++) left += steiner[v] && !done[v];
	st = vim_store(out_mesh, &m, mesh);
out:
	if (out_moved)
		*out_moved = moved;
	if (out_left)
		*out_left = left;
	qaws_internal_free(done);
	vim_free(&m);
	return st;
}
