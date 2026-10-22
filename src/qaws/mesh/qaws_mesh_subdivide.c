/*
 * Subdivision of meshes (qaws_mesh_subdivide.h): Loop and Catmull-Clark
 * levels on flat arrays (a row per vertex: the position, then the vertex
 * attributes; a row per corner and per face for theirs), the creases
 * semi-sharp, the result built once after the last level.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_subdivide.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SDV_INF 1e30          /* an infinitely sharp crease (borders, no sharpness given) */
#define SDV_MAX 0x7FFFFFFFu   /* the largest count a level may reach */

typedef struct sdv_attr
{
	unsigned int element, rule, off, n;   /* off: the first column in its element's row */
} sdv_attr;

typedef struct sdv_cols
{
	unsigned int kv, kc, kf, na;          /* row widths: vertex (3 + attributes), corner, face */
	sdv_attr const* at;
} sdv_cols;

typedef struct sdv_mesh
{
	unsigned int nv, nf, nc, nfe;
	unsigned int *fo, *fv;                /* per face its first corner (nf + 1); per corner its vertex */
	double *X, *C, *F;                    /* the rows */
	unsigned char* vk;                    /* per vertex: a corner, it stays */
	unsigned int* fe;                     /* the feature edges, 2 vertices each */
	double* fs;                           /* their sharpness */
} sdv_mesh;

typedef struct sdv_key
{
	unsigned int a, b, c;                 /* an edge's smaller and larger vertex; a corner on it */
} sdv_key;

static int sdv_key_cmp(void const* p, void const* q)
{
	sdv_key const *x = (sdv_key const*)p, *y = (sdv_key const*)q;
	if (x->a != y->a)
		return x->a < y->a ? -1 : 1;
	if (x->b != y->b)
		return x->b < y->b ? -1 : 1;
	return x->c < y->c ? -1 : (x->c > y->c ? 1 : 0);
}

static void sdv_free(sdv_mesh* m)
{
	qaws_internal_free(m->fo);
	qaws_internal_free(m->fv);
	qaws_internal_free(m->X);
	qaws_internal_free(m->C);
	qaws_internal_free(m->F);
	qaws_internal_free(m->vk);
	qaws_internal_free(m->fe);
	qaws_internal_free(m->fs);
	memset(m, 0, sizeof *m);
}

static int sdv_alloc(sdv_mesh* m, sdv_cols const* k)
{
	m->fo = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)m->nf + 1));
	m->fv = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * ((size_t)m->nc + 1));
	m->X = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)m->nv * k->kv + 1));
	m->C = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)m->nc * k->kc + 1));
	m->F = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)m->nf * k->kf + 1));
	m->vk = (unsigned char*)qaws_internal_zalloc((size_t)m->nv + 1, 1);
	m->fe = (unsigned int*)qaws_internal_malloc(sizeof(unsigned int) * (2 * (size_t)m->nfe + 1));
	m->fs = (double*)qaws_internal_malloc(sizeof(double) * ((size_t)m->nfe + 1));
	return m->fo && m->fv && m->X && m->C && m->F && m->vk && m->fe && m->fs;
}

/* y = sum of w[i] * x[i] over n rows of k columns */
static void sdv_mix(double* y, double const* const* x, double const* w, unsigned int n, unsigned int k)
{
	unsigned int i, j;
	for (j = 0; j < k; j++)
	{
		double s = 0;
		for (i = 0; i < n; i++) s += w[i] * x[i][j];
		y[j] = s;
	}
}

/* the attributes of one element in a new row by their rules: NEAREST from pick, NORMALIZED rescaled */
static void sdv_rules(double* row, double const* pick, sdv_cols const* k, unsigned int element)
{
	unsigned int a, j;
	for (a = 0; a < k->na; a++)
	{
		sdv_attr const* t = &k->at[a];
		double l = 0;
		if (t->element != element)
			continue;
		if (t->rule == QAWS_ATTRIBUTE_NEAREST)
			memcpy(row + t->off, pick + t->off, sizeof(double) * t->n);
		else if (t->rule == QAWS_ATTRIBUTE_NORMALIZED)
		{
			for (j = 0; j < t->n; j++) l += row[t->off + j] * row[t->off + j];
			l = sqrt(l);
			for (j = 0; l > 1e-300 && j < t->n; j++) row[t->off + j] /= l;
		}
	}
}

/* the edge (a, b) among the ne edges sorted by their two vertices (smaller first), HM_NONE if none */
static unsigned int sdv_find_edge(unsigned int const* ev, unsigned int ne, unsigned int a, unsigned int b)
{
	unsigned int lo = 0, hi = ne;
	if (a > b)
	{
		unsigned int s = a;
		a = b;
		b = s;
	}
	while (lo < hi)
	{
		unsigned int mid = lo + (hi - lo) / 2;
		if (ev[2 * mid] < a || (ev[2 * mid] == a && ev[2 * mid + 1] < b)) lo = mid + 1;
		else hi = mid;
	}
	return lo < ne && ev[2 * lo] == a && ev[2 * lo + 1] == b ? lo : HM_NONE;
}

/* the triangles' corner c's next vertex (d 1) or previous one (d 2) in its face */
#define SDV_TRI(m, cf, c, d) ((m)->fv[(m)->fo[(cf)[c]] + ((c) - (m)->fo[(cf)[c]] + (d)) % 3])

/* v's ring of neighbours in order from start into out (corners of v from co[v] in cv, their faces in cf); its size,
   0 when the ring is not closed (a border, not a manifold) or has fewer than 3 */
static unsigned int sdv_ring(unsigned int* out, sdv_mesh const* m, unsigned int const* co, unsigned int const* cv, unsigned int const* cf,
	unsigned int v, unsigned int start)
{
	unsigned int n = co[v + 1] - co[v], c = HM_NONE, c0, q, i;
	for (q = co[v]; q < co[v + 1]; q++)
		if (SDV_TRI(m, cf, cv[q], 1) == start)
			c = cv[q];
	if (c == HM_NONE || n < 3)
		return 0;
	c0 = c;
	for (i = 0; i < n; i++)
	{
		unsigned int pv = SDV_TRI(m, cf, c, 2), b = HM_NONE;
		out[i] = SDV_TRI(m, cf, c, 1);
		for (q = co[v]; q < co[v + 1]; q++)
			if (SDV_TRI(m, cf, cv[q], 1) == pv)
				b = cv[q];
		if (b == HM_NONE)
			return 0;
		c = b;
	}
	return c == c0 ? n : 0;
}

/* v's sector holding the smooth edge to w: its neighbours from one sharp edge (es > 0, a border's too) round to the other,
   into out (p_0 .. p_n, n triangles), *iw w's place; 0 when none (a ring closed without a sharp edge, not a manifold) */
static unsigned int sdv_sector(unsigned int* out, sdv_mesh const* m, unsigned int const* co, unsigned int const* cv, unsigned int const* cf,
	unsigned int const* ev, unsigned int ne, double const* es, unsigned int v, unsigned int w, unsigned int* iw)
{
	unsigned int n = co[v + 1] - co[v], nb = 0, nf = 0, d, q, x, guard;
	unsigned int* back = out + n + 2;   /* (the backward walk, reversed in after) */
	for (d = 0; d < 2; d++)
	{
		/* d 0: forward (w's corner's previous vertex), 1: backward */
		x = w;
		for (guard = 0; guard <= n; guard++)
		{
			unsigned int c = HM_NONE, y, e;
			for (q = co[v]; q < co[v + 1]; q++)
				if (SDV_TRI(m, cf, cv[q], d ? 2 : 1) == x)
					c = cv[q];
			if (c == HM_NONE)
				return 0;   /* (a border without its edge marked) */
			y = SDV_TRI(m, cf, c, d ? 1 : 2);
			if (y == w)
				return 0;   /* (round the whole ring: no sharp edge) */
			if (d) back[nb++] = y; else out[1 + nf++] = y;
			e = sdv_find_edge(ev, ne, v, y);
			if (e == HM_NONE)
				return 0;
			if (es[e] > 0)
				break;
			x = y;
		}
		if (guard > n)
			return 0;
	}
	/* p_0 .. p_n: the backward walk reversed, w, the forward walk */
	memmove(out + nb + 1, out + 1, sizeof(unsigned int) * nf);
	for (q = 0; q < nb; q++) out[q] = back[nb - 1 - q];
	out[nb] = w;
	*iw = nb;
	return nb + nf + 1;
}

/* one level of m into o (loop: 1 Loop, 2 the modified butterfly, 0 Catmull-Clark); limit: o is m with its vertices moved
   onto the limit surface instead */
static qaws_status sdv_step(sdv_mesh* o, sdv_mesh const* m, int loop, sdv_cols const* k, int limit)
{
	unsigned int nv = m->nv, nf = m->nf, nc = m->nc, kv = k->kv, kc = k->kc, kf = k->kf, ne = 0, i, j, f, e, t;
	size_t nv2, nf2, nc2;
	sdv_key* key = (sdv_key*)qaws_internal_scratch_malloc(sizeof(sdv_key) * ((size_t)nc + 1));
	unsigned int* ui = (unsigned int*)qaws_internal_scratch_zalloc(9 * (size_t)nc + 6 * (size_t)nv + m->nfe + 1, sizeof(unsigned int));
	double* dd = (double*)qaws_internal_scratch_zalloc((size_t)nc + (size_t)nv + 2 * (size_t)nv * kv + (loop ? 0 : (size_t)nf * kv) + 4 * (size_t)kv + 1,
		sizeof(double));
	unsigned int *ev, *ef, *eo, *en, *ce, *cf, *val, *nfc, *nsh, *sh2, *fidx, *co = NULL, *cv = NULL, *ring = NULL;
	double *es, *ssum, *NS, *QS, *FP, *sm, *rr;
	qaws_status st = QAWS_STATUS_OK;
	memset(o, 0, sizeof *o);
	if (!key || !ui || !dd)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	ev = ui;                    /* per edge its two vertices (smaller first) */
	ef = ev + 2 * (size_t)nc;   /* its faces (HM_NONE: none) */
	eo = ef + 2 * (size_t)nc;   /* Loop: the vertex opposite in each */
	en = eo + 2 * (size_t)nc;   /* its faces' count */
	ce = en + nc;               /* per corner its edge */
	cf = ce + nc;               /* per corner its face */
	val = cf + nc;              /* per vertex: edges */
	nfc = val + nv;             /* faces */
	nsh = nfc + nv;             /* sharp edges */
	sh2 = nsh + nv;             /* the first two sharp edges' far ends */
	fidx = sh2 + 2 * (size_t)nv;
	es = dd;                    /* per edge its sharpness */
	ssum = es + nc;             /* per vertex its sharp edges' sharpness, summed */
	NS = ssum + nv;             /* its neighbours' rows, summed */
	QS = NS + (size_t)nv * kv;  /* Catmull-Clark: its faces' points, summed */
	sm = QS + (size_t)nv * kv;  /* two scratch rows */
	rr = sm + kv;
	FP = rr + kv;               /* Catmull-Clark: the face points */

	/* the edges: corners sorted by their edge */
	for (f = 0; f < nf; f++)
		for (j = m->fo[f]; j < m->fo[f + 1]; j++)
		{
			unsigned int a = m->fv[j], b = m->fv[j + 1 == m->fo[f + 1] ? m->fo[f] : j + 1];
			key[j].a = a < b ? a : b;
			key[j].b = a < b ? b : a;
			key[j].c = j;
			cf[j] = f;
		}
	qsort(key, nc, sizeof(sdv_key), sdv_key_cmp);
	for (i = 0; i < nc; i++)
	{
		unsigned int c = key[i].c, p;
		if (i == 0 || key[i].a != key[i - 1].a || key[i].b != key[i - 1].b)
		{
			e = ne++;
			ev[2 * e] = key[i].a;
			ev[2 * e + 1] = key[i].b;
			ef[2 * e] = ef[2 * e + 1] = HM_NONE;
		}
		e = ne - 1;
		ce[c] = e;
		p = c == m->fo[cf[c]] ? m->fo[cf[c] + 1] - 1 : c - 1;
		if (en[e] < 2)
		{
			ef[2 * e + en[e]] = cf[c];
			eo[2 * e + en[e]] = m->fv[p];
		}
		en[e]++;
	}
	for (e = 0; e < ne; e++) es[e] = en[e] < 2 ? SDV_INF : 0;
	for (t = 0; t < m->nfe; t++)
	{
		unsigned int lo = sdv_find_edge(ev, ne, m->fe[2 * t], m->fe[2 * t + 1]);
		fidx[t] = lo;
		if (lo != HM_NONE && es[lo] < SDV_INF)
			es[lo] = m->fs[t];
	}
	if (loop == 2 && !limit)
	{
		/* the butterfly's rings: per vertex its corners */
		unsigned int* fill;
		co = (unsigned int*)qaws_internal_scratch_zalloc(2 * (size_t)nv + 2 * (size_t)nc + 2, sizeof(unsigned int));
		if (!co)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		cv = co + nv + 1;
		ring = cv + nc;
		fill = ring + nc + 1;
		for (j = 0; j < nc; j++) co[m->fv[j] + 1]++;
		for (i = 0; i < nv; i++) co[i + 1] += co[i];
		for (j = 0; j < nc; j++) cv[co[m->fv[j]] + fill[m->fv[j]]++] = j;
	}

	nv2 = limit ? nv : (size_t)nv + ne + (loop ? 0 : nf);
	nf2 = limit ? nf : (loop ? 4 * (size_t)nf : nc);
	nc2 = limit ? nc : (loop ? 12 * (size_t)nf : 4 * (size_t)nc);
	if (nv2 > SDV_MAX || nf2 > SDV_MAX || nc2 > SDV_MAX)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	o->nv = (unsigned int)nv2;
	o->nf = (unsigned int)nf2;
	o->nc = (unsigned int)nc2;
	o->nfe = limit ? m->nfe : 2 * m->nfe;
	if (!sdv_alloc(o, k))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if (limit)
	{
		memcpy(o->fo, m->fo, sizeof(unsigned int) * ((size_t)nf + 1));
		memcpy(o->fv, m->fv, sizeof(unsigned int) * nc);
		memcpy(o->C, m->C, sizeof(double) * (size_t)nc * kc);
		memcpy(o->F, m->F, sizeof(double) * (size_t)nf * kf);
		memcpy(o->vk, m->vk, nv);
		memcpy(o->fe, m->fe, sizeof(unsigned int) * 2 * (size_t)m->nfe);
		memcpy(o->fs, m->fs, sizeof(double) * m->nfe);
	}

	/* the vertices' sums */
	for (e = 0; e < ne; e++)
		for (i = 0; i < 2; i++)
		{
			unsigned int v = ev[2 * e + i], w = ev[2 * e + 1 - i];
			val[v]++;
			for (j = 0; j < kv; j++) NS[(size_t)v * kv + j] += m->X[(size_t)w * kv + j];
			if (es[e] > 0)
			{
				if (nsh[v] < 2)
					sh2[2 * v + nsh[v]] = w;
				nsh[v]++;
				ssum[v] += es[e];
			}
		}
	for (f = 0; limit && !loop && f < nf; f++)
		for (j = m->fo[f]; m->fo[f + 1] - m->fo[f] == 4 && j < m->fo[f + 1]; j++)
		{
			/* the limit's diagonal neighbours (quads, as every face is after a level), summed in QS */
			unsigned int opp = m->fv[m->fo[f] + (j - m->fo[f] + 2) % 4];
			nfc[m->fv[j]]++;
			for (i = 0; i < kv; i++) QS[(size_t)m->fv[j] * kv + i] += m->X[(size_t)opp * kv + i];
		}
	for (f = 0; !limit && f < nf; f++)
	{
		double* y = FP + (size_t)f * kv;
		unsigned int n = m->fo[f + 1] - m->fo[f];
		for (j = m->fo[f]; !loop && j < m->fo[f + 1]; j++)
			for (i = 0; i < kv; i++) y[i] += m->X[(size_t)m->fv[j] * kv + i] / n;
		for (j = m->fo[f]; j < m->fo[f + 1]; j++)
		{
			nfc[m->fv[j]]++;
			for (i = 0; !loop && i < kv; i++) QS[(size_t)m->fv[j] * kv + i] += y[i];
		}
	}

	if (limit)
	{
		/* the limit positions: Loop's (1 - n chi, chi each neighbour, chi = 1 / (3 / (8 beta) + n)), Catmull-Clark's
		   ((n^2 V + 4 E + F) / (n (n + 5)), E and F the edge and diagonal neighbours summed), a crease's (1, 4, 1) / 6 */
		if (loop)
			for (f = 0; f < nf; f++)
				for (j = m->fo[f]; j < m->fo[f + 1]; j++) nfc[m->fv[j]]++;
		for (i = 0; i < nv; i++)
		{
			double const* x = m->X + (size_t)i * kv;
			double* y = o->X + (size_t)i * kv;
			double n = (double)val[i];
			memcpy(y, x, sizeof(double) * kv);
			if (val[i] > 0 && loop)
			{
				double c = 0.375 + 0.25 * cos(6.283185307179586 / n), beta = (0.625 - c * c) / n, chi = 1 / (0.375 / beta + n);
				for (j = 0; j < kv; j++) y[j] = (1 - n * chi) * x[j] + chi * NS[(size_t)i * kv + j];
			}
			else if (val[i] > 0 && nfc[i] == val[i])
				for (j = 0; j < kv; j++) y[j] = (n * n * x[j] + 4 * NS[(size_t)i * kv + j] + QS[(size_t)i * kv + j]) / (n * (n + 5));
			if (nfc[i] == 1)
				memcpy(y, x, sizeof(double) * kv);
			else if (nsh[i] >= 2 || (m->vk[i] && nsh[i] >= 1))
			{
				double s = ssum[i] / nsh[i];
				if (nsh[i] == 2 && !m->vk[i])
					for (j = 0; j < kv; j++)
						rr[j] = (m->X[(size_t)sh2[2 * i] * kv + j] + 4 * x[j] + m->X[(size_t)sh2[2 * i + 1] * kv + j]) / 6;
				else
					memcpy(rr, x, sizeof(double) * kv);
				for (j = 0; j < kv; j++) y[j] = s >= 1 ? rr[j] : (1 - s) * y[j] + s * rr[j];
			}
			sdv_rules(y, x, k, QAWS_MESH_VERTEX);
		}
		goto done;
	}

	/* vertex points */
	for (i = 0; i < nv; i++)
	{
		double const* x = m->X + (size_t)i * kv;
		double* y = o->X + (size_t)i * kv;
		double n = (double)val[i];
		memcpy(y, x, sizeof(double) * kv);
		o->vk[i] = m->vk[i];
		if (loop == 2)
			continue;   /* (the butterfly interpolates: the vertices stay) */
		if (val[i] > 0 && loop)
		{
			double c = 0.375 + 0.25 * cos(6.283185307179586 / n), beta = (0.625 - c * c) / n;
			for (j = 0; j < kv; j++) y[j] = (1 - n * beta) * x[j] + beta * NS[(size_t)i * kv + j];
		}
		else if (val[i] > 0 && nfc[i] > 0)
			for (j = 0; j < kv; j++) y[j] = (QS[(size_t)i * kv + j] / nfc[i] + NS[(size_t)i * kv + j] / n + (n - 2) * x[j]) / n;
		if (nfc[i] == 1)
			memcpy(y, x, sizeof(double) * kv);
		else if (nsh[i] >= 2 || (m->vk[i] && nsh[i] >= 1))
		{
			/* a crease or a corner (a feature corner while one of its creases is sharp), blended below sharpness 1 */
			double s = ssum[i] / nsh[i];
			if (nsh[i] == 2 && !m->vk[i])
			{
				double const* r3[3];
				static double const w3[3] = { 0.125, 0.75, 0.125 };
				r3[0] = m->X + (size_t)sh2[2 * i] * kv;
				r3[1] = x;
				r3[2] = m->X + (size_t)sh2[2 * i + 1] * kv;
				sdv_mix(rr, r3, w3, 3, kv);
			}
			else
				memcpy(rr, x, sizeof(double) * kv);
			for (j = 0; j < kv; j++) y[j] = s >= 1 ? rr[j] : (1 - s) * y[j] + s * rr[j];
		}
		sdv_rules(y, x, k, QAWS_MESH_VERTEX);
		o->vk[i] = m->vk[i];
	}
	/* edge points */
	for (e = 0; e < ne; e++)
	{
		double const* a = m->X + (size_t)ev[2 * e] * kv;
		double const* b = m->X + (size_t)ev[2 * e + 1] * kv;
		double* y = o->X + ((size_t)nv + e) * kv;
		for (j = 0; j < kv; j++) rr[j] = 0.5 * (a[j] + b[j]);
		if (loop == 2)
		{
			/* the butterfly: a crease's four points, -1/16, 9/16, 9/16, -1/16 (an end not inside a crease reflected, 2 a - b) */
			unsigned int u = ev[2 * e], w = ev[2 * e + 1];
			unsigned int pu = nsh[u] == 2 && !m->vk[u] ? (sh2[2 * u] == w ? sh2[2 * u + 1] : sh2[2 * u]) : HM_NONE;
			unsigned int pw = nsh[w] == 2 && !m->vk[w] ? (sh2[2 * w] == u ? sh2[2 * w + 1] : sh2[2 * w]) : HM_NONE;
			for (j = 0; j < kv; j++)
				rr[j] = (9 * (a[j] + b[j]) - (pu != HM_NONE ? m->X[(size_t)pu * kv + j] : 2 * a[j] - b[j]) -
					(pw != HM_NONE ? m->X[(size_t)pw * kv + j] : 2 * b[j] - a[j])) / 16;
		}
		if (loop == 2 && es[e] < 1)
		{
			/* smooth: an extraordinary end's ring (Zorin, Schroeder and Sweldens 1996: 3/4 the end, s_j its j-th
			   neighbour from the other end), both ends' averaged; between regular ends the ten points (1/2 the ends, 1/8
			   the opposite vertices, -1/16 the wings across their edges; a wing beyond a border virtual, reflected across it, as
			   Zorin, Schroeder and Sweldens do) */
			unsigned int ends[2], ext = 0, s, r;
			ends[0] = ev[2 * e];
			ends[1] = ev[2 * e + 1];
			for (j = 0; j < kv; j++) sm[j] = 0;
			for (s = 0; s < 2; s++)
			{
				unsigned int nr = sdv_ring(ring, m, co, cv, cf, ends[s], ends[1 - s]), np, iw = 0;
				double const* xv = m->X + (size_t)ends[s] * kv;
				if (!nr && es[e] == 0 && nsh[ends[s]] == 2 && !m->vk[ends[s]] &&
					(np = sdv_sector(ring, m, co, cv, cf, ev, ne, es, ends[s], ends[1 - s], &iw)) >= 3 && np != 4)
				{
					/* an extraordinary crease or border end (n triangles in its sector, not 3): Zorin's rule (Zorin and
					   Schroeder, SIGGRAPH 2000 course notes, figure 4.7), the valence 2 n ring reflected oddly across the
					   crease: theta = pi / n, the other end p_i, c_j = (sin i theta sin j theta + sin 2 i theta sin 2 j
					   theta / 2) / n inside, the end 1 - sin theta sin i theta / (n (1 - cos theta)), p_0 and p_n
					   +-(cos i theta / 4 - sin 2 theta sin 2 i theta / (4 n (cos theta - cos 2 theta))) */
					unsigned int nt = np - 1;
					double th = 3.14159265358979323846 / nt, si = sin(iw * th), s2 = sin(2 * iw * th);
					double c0 = 1 - sin(th) * si / (nt * (1 - cos(th)));
					double cb = 0.25 * cos(iw * th) - sin(2 * th) * s2 / (4 * nt * (cos(th) - cos(2 * th)));
					for (j = 0; j < kv; j++)
						sm[j] += c0 * xv[j] + cb * (m->X[(size_t)ring[0] * kv + j] - m->X[(size_t)ring[nt] * kv + j]);
					for (r = 1; r < nt; r++)
					{
						double cj = (si * sin(r * th) + 0.5 * s2 * sin(2 * r * th)) / nt;
						for (j = 0; j < kv; j++) sm[j] += cj * m->X[(size_t)ring[r] * kv + j];
					}
					ext++;
					continue;
				}
				if (!nr || nr == 6)
					continue;
				for (j = 0; j < kv; j++) sm[j] += 0.75 * xv[j];
				for (r = 0; r < nr; r++)
				{
					double sj = nr == 3 ? (r ? -1.0 / 12 : 5.0 / 12) : nr == 4 ? (r == 0 ? 0.375 : (r == 2 ? -0.125 : 0))
						: (0.25 + cos(6.283185307179586 * r / nr) + 0.5 * cos(12.566370614359172 * r / nr)) / nr;
					for (j = 0; j < kv; j++) sm[j] += sj * m->X[(size_t)ring[r] * kv + j];
				}
				ext++;
			}
			if (ext)
				for (j = 0; j < kv; j++) sm[j] /= ext;
			else
			{
				unsigned int cd[2];
				cd[0] = eo[2 * e];
				cd[1] = eo[2 * e + 1];
				for (j = 0; j < kv; j++)
					sm[j] = 0.5 * (a[j] + b[j]) + 0.125 * (m->X[(size_t)cd[0] * kv + j] + m->X[(size_t)cd[1] * kv + j]);
				for (r = 0; r < 4; r++)
				{
					unsigned int p = ends[r & 1], q = cd[r >> 1], z = ends[1 - (r & 1)], e2 = sdv_find_edge(ev, ne, p, q);
					unsigned int g = e2 != HM_NONE && en[e2] == 2 ? (eo[2 * e2] == z ? eo[2 * e2 + 1] : eo[2 * e2]) : HM_NONE;
					for (j = 0; j < kv; j++)
						sm[j] -= (g != HM_NONE ? m->X[(size_t)g * kv + j] : m->X[(size_t)p * kv + j] + m->X[(size_t)q * kv + j] - m->X[(size_t)z * kv + j]) / 16;
				}
			}
			for (j = 0; j < kv; j++) y[j] = (1 - es[e]) * sm[j] + es[e] * rr[j];
		}
		else if (loop == 2)
			memcpy(y, rr, sizeof(double) * kv);
		else if (es[e] < 1)
		{
			double const* r4[4];
			static double const wl[4] = { 0.375, 0.375, 0.125, 0.125 }, wc[4] = { 0.25, 0.25, 0.25, 0.25 };
			r4[0] = a;
			r4[1] = b;
			r4[2] = loop ? m->X + (size_t)eo[2 * e] * kv : FP + (size_t)ef[2 * e] * kv;
			r4[3] = loop ? m->X + (size_t)eo[2 * e + 1] * kv : FP + (size_t)ef[2 * e + 1] * kv;
			sdv_mix(sm, r4, loop ? wl : wc, 4, kv);
			for (j = 0; j < kv; j++) y[j] = (1 - es[e]) * sm[j] + es[e] * rr[j];
		}
		else
			memcpy(y, rr, sizeof(double) * kv);
		sdv_rules(y, a, k, QAWS_MESH_VERTEX);
	}
	/* face points */
	for (f = 0; !loop && f < nf; f++)
	{
		double* y = o->X + ((size_t)nv + ne + f) * kv;
		memcpy(y, FP + (size_t)f * kv, sizeof(double) * kv);
		sdv_rules(y, m->X + (size_t)m->fv[m->fo[f]] * kv, k, QAWS_MESH_VERTEX);
	}

	/* faces, their corners' rows (a child corner the mean of parent corners pa and pb of its face; pb 4: the centroid) */
	{
		static unsigned char const lc[12][2] = { { 0, 0 }, { 0, 1 }, { 2, 0 }, { 1, 1 }, { 1, 2 }, { 0, 1 }, { 2, 2 }, { 2, 0 }, { 1, 2 },
			{ 0, 1 }, { 1, 2 }, { 2, 0 } };
		unsigned int q = 0, g = 0;
		for (f = 0; f < nf; f++)
		{
			unsigned int c0 = m->fo[f], n = m->fo[f + 1] - c0, nk = loop ? 4 : n, r;
			for (r = 0; r < nk; r++, g++)
			{
				unsigned int cs = loop ? 3 : 4, s;
				o->fo[g] = q;
				memcpy(o->F + (size_t)g * kf, m->F + (size_t)f * kf, sizeof(double) * kf);
				for (s = 0; s < cs; s++, q++)
				{
					unsigned int pa, pb, v;
					double* y = o->C + (size_t)q * kc;
					if (loop)
					{
						static unsigned int const lv[4][3] = { { 0, 3, 5 }, { 1, 4, 3 }, { 2, 5, 4 }, { 3, 4, 5 } };
						unsigned int id = lv[r][s];
						pa = lc[3 * r + s][0];
						pb = lc[3 * r + s][1];
						v = id < 3 ? m->fv[c0 + id] : nv + ce[c0 + id - 3];
					}
					else
					{
						pa = s == 3 ? (r + n - 1) % n : r;
						pb = s == 0 ? r : (s == 2 ? 4 : (s == 1 ? (r + 1) % n : r));
						v = s == 0 ? m->fv[c0 + r] : (s == 1 ? nv + ce[c0 + r] : (s == 2 ? nv + ne + f : nv + ce[c0 + (r + n - 1) % n]));
						if (s == 2)
							pa = 0;
					}
					o->fv[q] = v;
					if (pb == 4 && !loop)
					{
						for (j = 0; j < kc; j++) y[j] = 0;
						for (i = 0; i < n; i++)
							for (j = 0; j < kc; j++) y[j] += m->C[((size_t)c0 + i) * kc + j] / n;
					}
					else
						for (j = 0; j < kc; j++) y[j] = 0.5 * (m->C[((size_t)c0 + pa) * kc + j] + m->C[((size_t)c0 + pb) * kc + j]);
					sdv_rules(y, m->C + ((size_t)c0 + pa) * kc, k, QAWS_MESH_CORNER);
				}
			}
		}
		o->fo[g] = q;
	}

	/* the creases' halves, while sharp */
	o->nfe = 0;
	for (t = 0; t < m->nfe; t++)
	{
		double s = m->fs[t] >= SDV_INF ? SDV_INF : m->fs[t] - 1;
		unsigned int mid;
		if (fidx[t] == HM_NONE || s <= 0)
			continue;
		mid = nv + fidx[t];
		o->fe[2 * o->nfe] = m->fe[2 * t];
		o->fe[2 * o->nfe + 1] = mid;
		o->fs[o->nfe++] = s;
		o->fe[2 * o->nfe] = mid;
		o->fe[2 * o->nfe + 1] = m->fe[2 * t + 1];
		o->fs[o->nfe++] = s;
	}
done:
	if (st != QAWS_STATUS_OK)
		sdv_free(o);
	qaws_internal_free(key);
	qaws_internal_free(ui);
	qaws_internal_free(dd);
	qaws_internal_free(co);
	return st;
}

/* the limit surface's unit normals at m's vertices, 3 per vertex in N, by the tangent masks on m's control points (Loop:
   sum cos(2 pi i / n) e_i and sum sin(2 pi i / n) e_i over the ring in order; Catmull-Clark's of Halstead, Kass and DeRose,
   the edge neighbours by A_n cos, the diagonal ones by cos + the next cos); a vertex on a border, a crease or a corner by
   its faces' normals at the limit positions L, rows of kv (shortcut: not the creases' own masks, Biermann et al. 2000,
   for exact normals there). 0 on allocation failure. */
static int sdv_limit_normals(double* N, sdv_mesh const* m, double const* L, int loop, unsigned int kv)
{
	unsigned int nv = m->nv, nc = m->nc, *co = (unsigned int*)qaws_internal_scratch_zalloc(2 * (size_t)nv + 2 * (size_t)nc + 2, sizeof(unsigned int));
	unsigned char* odd = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nv + 1, 1);
	unsigned int *cv, *fc, *cf, v, f, j, i;
	if (!co || !odd)
	{
		qaws_internal_free(co);
		qaws_internal_free(odd);
		return 0;
	}
	cv = co + nv + 1;   /* per vertex its corners, from co[v] */
	cf = cv + nc;       /* per corner its face */
	fc = cf + nc;       /* per vertex its corners placed */
	for (j = 0; j < nc; j++) co[m->fv[j] + 1]++;
	for (v = 0; v < nv; v++) co[v + 1] += co[v];
	for (f = 0; f < m->nf; f++)
		for (j = m->fo[f]; j < m->fo[f + 1]; j++)
		{
			cf[j] = f;
			cv[co[m->fv[j]] + fc[m->fv[j]]++] = j;
		}
	for (j = 0; j < m->nfe; j++)
		odd[m->fe[2 * j]] = odd[m->fe[2 * j + 1]] = 1;
	for (v = 0; v < nv; v++)
	{
		unsigned int n = co[v + 1] - co[v], c = n ? cv[co[v]] : 0, k;
		double t1[3] = { 0, 0, 0 }, t2[3] = { 0, 0, 0 }, *y = N + 3 * (size_t)v, l, w = 6.283185307179586 / (n ? n : 1),
			An = 1 + cos(w) + cos(w / 2) * sqrt(2 * (9 + cos(w)));
		odd[v] |= m->vk[v] | !n;
		/* the ring in order: from a corner, the next one around v is the face whose edge out of v is this one's edge in */
		for (k = 0; !odd[v] && k < n; k++)
		{
			unsigned int s = m->fo[cf[c]], fn = m->fo[cf[c] + 1] - s, nx = m->fv[s + (c - s + 1) % fn], pv = m->fv[s + (c - s + fn - 1) % fn],
				dg = m->fv[s + (c - s + 2) % fn], q, b = HM_NONE;
			double ca = cos(w * k), sa = sin(w * k), cb = cos(w * (k + 1)), sb = sin(w * (k + 1));
			for (i = 0; i < 3; i++)
			{
				double e = m->X[(size_t)nx * kv + i], d = m->X[(size_t)dg * kv + i];
				t1[i] += loop ? ca * e : An * ca * e + (ca + cb) * d;
				t2[i] += loop ? sa * e : An * sa * e + (sa + sb) * d;
			}
			for (q = co[v]; q < co[v + 1]; q++)
			{
				unsigned int h0 = m->fo[cf[cv[q]]], hn = m->fo[cf[cv[q]] + 1] - h0;
				if (m->fv[h0 + (cv[q] - h0 + 1) % hn] == pv)
					b = cv[q];
			}
			odd[v] = b == HM_NONE || (!loop && fn != 4);   /* a border; a face that is not a quad */
			c = b;
		}
		if (!odd[v] && c != cv[co[v]])
			odd[v] = 1;   /* (the ring not closed after n faces: not a manifold vertex) */
		y[0] = t1[1] * t2[2] - t1[2] * t2[1];
		y[1] = t1[2] * t2[0] - t1[0] * t2[2];
		y[2] = t1[0] * t2[1] - t1[1] * t2[0];
		for (k = co[v]; odd[v] && k < co[v + 1]; k++)
		{
			unsigned int s = m->fo[cf[cv[k]]], fn = m->fo[cf[cv[k]] + 1] - s;
			double const *p = L + (size_t)v * kv, *a = L + (size_t)m->fv[s + (cv[k] - s + 1) % fn] * kv, *b = L + (size_t)m->fv[s + (cv[k] - s + fn - 1) % fn] * kv;
			if (k == co[v])
				y[0] = y[1] = y[2] = 0;
			y[0] += (a[1] - p[1]) * (b[2] - p[2]) - (a[2] - p[2]) * (b[1] - p[1]);
			y[1] += (a[2] - p[2]) * (b[0] - p[0]) - (a[0] - p[0]) * (b[2] - p[2]);
			y[2] += (a[0] - p[0]) * (b[1] - p[1]) - (a[1] - p[1]) * (b[0] - p[0]);
		}
		l = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
		for (i = 0; l > 1e-300 && i < 3; i++) y[i] /= l;
	}
	qaws_internal_free(co);
	qaws_internal_free(odd);
	return 1;
}

/* one level of Kobbelt's sqrt(3) on triangles m into o: a point at each triangle's centroid, each old vertex relaxed
   ((1 - a) v + a times its neighbours' mean, a = (4 - 2 cos(2 pi / n)) / 9), each old edge between two triangles flipped
   to join their centroids: a triangle per half-edge, (a, centroid across, own centroid) where flipped, (a, b, own
   centroid) where not. Border and feature edges stay unflipped and their vertices fixed (shortcut: Kobbelt's ternary
   rule for the border every second level not used; the border keeps its polyline) */
static qaws_status sdv_sqrt3(sdv_mesh* o, sdv_mesh const* m, sdv_cols const* k)
{
	unsigned int nv = m->nv, nf = m->nf, nc = m->nc, kv = k->kv, kc = k->kc, kf = k->kf, ne = 0, i, j, e, f, t;
	sdv_key* key = (sdv_key*)qaws_internal_scratch_malloc(sizeof(sdv_key) * ((size_t)nc + 1));
	unsigned int* ui = (unsigned int*)qaws_internal_scratch_zalloc(6 * (size_t)nc + (size_t)nv + 1, sizeof(unsigned int));
	unsigned char* fixed = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nv + 1, 1);
	double* NS = (double*)qaws_internal_scratch_zalloc((size_t)nv * kv + 1, sizeof(double));
	unsigned int *ev, *ec, *en, *ce, *val;
	qaws_status st = QAWS_STATUS_OK;
	memset(o, 0, sizeof *o);
	if (!key || !ui || !fixed || !NS)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	ev = ui;                    /* per edge its two vertices (smaller first) */
	ec = ev + 2 * (size_t)nc;   /* its first two corners */
	en = ec + 2 * (size_t)nc;   /* its corners' count (3: not flipped) */
	ce = en + nc;               /* per corner its edge */
	val = ce + nc;              /* per vertex its edges */
	for (f = 0; f < nf; f++)
		for (j = m->fo[f]; j < m->fo[f + 1]; j++)
		{
			unsigned int a = m->fv[j], b = m->fv[j + 1 == m->fo[f + 1] ? m->fo[f] : j + 1];
			key[j].a = a < b ? a : b;
			key[j].b = a < b ? b : a;
			key[j].c = j;
		}
	qsort(key, nc, sizeof(sdv_key), sdv_key_cmp);
	for (i = 0; i < nc; i++)
	{
		if (i == 0 || key[i].a != key[i - 1].a || key[i].b != key[i - 1].b)
		{
			e = ne++;
			ev[2 * e] = key[i].a;
			ev[2 * e + 1] = key[i].b;
		}
		e = ne - 1;
		ce[key[i].c] = e;
		if (en[e] < 2)
			ec[2 * e + en[e]] = key[i].c;
		en[e]++;
	}
	/* fixed: on the border, on a feature, a corner */
	for (e = 0; e < ne; e++)
		if (en[e] != 2)
			fixed[ev[2 * e]] = fixed[ev[2 * e + 1]] = 1;
	for (t = 0; t < m->nfe; t++)
	{
		unsigned int x = sdv_find_edge(ev, ne, m->fe[2 * t], m->fe[2 * t + 1]);
		if (x != HM_NONE)
		{
			en[x] = 3;   /* (unflipped like a border) */
			fixed[ev[2 * x]] = fixed[ev[2 * x + 1]] = 1;
		}
	}
	for (i = 0; i < nv; i++) fixed[i] |= m->vk[i];
	for (e = 0; e < ne; e++)
	{
		/* (an edge between two faces flips only when they run it opposite ways) */
		unsigned int a = ev[2 * e], b = ev[2 * e + 1];
		if (en[e] == 2 && m->fv[ec[2 * e]] == m->fv[ec[2 * e + 1]])
			en[e] = 3;
		val[a]++;
		val[b]++;
		for (j = 0; j < kv; j++)
		{
			NS[(size_t)a * kv + j] += m->X[(size_t)b * kv + j];
			NS[(size_t)b * kv + j] += m->X[(size_t)a * kv + j];
		}
	}
	if ((size_t)nv + nf > SDV_MAX || 3 * (size_t)nf > SDV_MAX || 9 * (size_t)nf > SDV_MAX)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	o->nv = nv + nf;
	o->nf = 3 * nf;
	o->nc = 9 * nf;
	o->nfe = m->nfe;
	if (!sdv_alloc(o, k))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	memcpy(o->fe, m->fe, sizeof(unsigned int) * 2 * (size_t)m->nfe);
	memcpy(o->fs, m->fs, sizeof(double) * m->nfe);
	/* the old vertices relaxed, the centroids */
	for (i = 0; i < nv; i++)
	{
		double const* x = m->X + (size_t)i * kv;
		double* y = o->X + (size_t)i * kv;
		double n = (double)val[i], a = (4 - 2 * cos(6.283185307179586 / (val[i] ? n : 1))) / 9;
		memcpy(y, x, sizeof(double) * kv);
		o->vk[i] = m->vk[i];
		if (fixed[i] || !val[i])
			continue;
		for (j = 0; j < kv; j++) y[j] = (1 - a) * x[j] + a * NS[(size_t)i * kv + j] / n;
		sdv_rules(y, x, k, QAWS_MESH_VERTEX);
	}
	for (f = 0; f < nf; f++)
	{
		double* y = o->X + ((size_t)nv + f) * kv;
		for (j = 0; j < kv; j++) y[j] = 0;
		for (i = m->fo[f]; i < m->fo[f + 1]; i++)
			for (j = 0; j < kv; j++) y[j] += m->X[(size_t)m->fv[i] * kv + j] / 3;
		sdv_rules(y, m->X + (size_t)m->fv[m->fo[f]] * kv, k, QAWS_MESH_VERTEX);
	}
	/* a triangle per half-edge, its parent's rows; corners: an old vertex its corner in the face, a centroid its face's corners' mean */
	for (f = 0; f < nf; f++)
	{
		unsigned int c0 = m->fo[f];
		for (i = 0; i < 3; i++)
		{
			unsigned int c = c0 + i, cn = c0 + (i + 1) % 3, g = 3 * f + i, q = 3 * g, x = ce[c], tw, gf, s;
			unsigned int src[3];   /* per new corner: its corner in m, or the face (as nc + face) whose mean it takes */
			o->fo[g] = q;
			memcpy(o->F + (size_t)g * kf, m->F + (size_t)f * kf, sizeof(double) * kf);
			tw = en[x] == 2 ? (ec[2 * x] == c ? ec[2 * x + 1] : ec[2 * x]) : HM_NONE;
			if (tw != HM_NONE)
			{
				/* flipped: (a, the twin's centroid, own centroid) */
				gf = tw / 3;   /* (triangles only) */
				o->fv[q] = m->fv[c];
				o->fv[q + 1] = nv + gf;
				o->fv[q + 2] = nv + f;
				src[0] = c;
				src[1] = nc + gf;
				src[2] = nc + f;
			}
			else
			{
				o->fv[q] = m->fv[c];
				o->fv[q + 1] = m->fv[cn];
				o->fv[q + 2] = nv + f;
				src[0] = c;
				src[1] = cn;
				src[2] = nc + f;
			}
			for (s = 0; s < 3; s++)
			{
				double* y = o->C + ((size_t)q + s) * kc;
				if (src[s] < nc)
					memcpy(y, m->C + (size_t)src[s] * kc, sizeof(double) * kc);
				else
				{
					unsigned int ff = src[s] - nc, r;
					for (j = 0; j < kc; j++) y[j] = 0;
					for (r = m->fo[ff]; r < m->fo[ff + 1]; r++)
						for (j = 0; j < kc; j++) y[j] += m->C[(size_t)r * kc + j] / 3;
					sdv_rules(y, m->C + (size_t)m->fo[ff] * kc, k, QAWS_MESH_CORNER);
				}
			}
		}
	}
	o->fo[3 * nf] = 9 * nf;
done:
	if (st != QAWS_STATUS_OK)
		sdv_free(o);
	qaws_internal_free(key);
	qaws_internal_free(ui);
	qaws_internal_free(fixed);
	qaws_internal_free(NS);
	return st;
}

/* one adaptive level of Loop or the butterfly on triangles m into o: the faces at a turn (a neighbour's normal more than
   angle away, or a feature edge on them) split, and every edge of theirs; a face with one or two split edges cut into two
   or three (no T-junction), with three into four. The positions are the uniform level's (full, the rules all in one
   place): a split edge's point, and an old vertex's moved by its rule where it touches a split edge, kept elsewhere. */
static qaws_status sdv_adaptive(sdv_mesh* o, sdv_mesh const* m, int loop, sdv_cols const* k, double angle)
{
	sdv_mesh full;
	unsigned int nv = m->nv, nf = m->nf, nc = m->nc, kv = k->kv, kc = k->kc, kf = k->kf, ne = 0, i, j, e, f, t, ns = 0, nfo = 0, nco = 0, nfe2 = 0;
	sdv_key* key = (sdv_key*)qaws_internal_scratch_malloc(sizeof(sdv_key) * ((size_t)nc + 1));
	unsigned int* ui = (unsigned int*)qaws_internal_scratch_zalloc(8 * (size_t)nc + 2 * (size_t)nv + 2, sizeof(unsigned int));
	unsigned char* mark = (unsigned char*)qaws_internal_scratch_zalloc((size_t)nf + 1, 1);
	double* nrm = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)nf + 1));
	unsigned int *ev, *ef, *en, *ce, *eid, *touch;
	double cosa = cos(angle);
	qaws_status st;
	memset(o, 0, sizeof *o);
	memset(&full, 0, sizeof full);
	if (!key || !ui || !mark || !nrm)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	st = sdv_step(&full, m, loop, k, 0);
	if (st != QAWS_STATUS_OK)
		goto done;
	ev = ui;                    /* per edge its two vertices (smaller first) */
	ef = ev + 2 * (size_t)nc;   /* its first two faces */
	en = ef + 2 * (size_t)nc;   /* its faces' count */
	ce = en + nc;               /* per corner its edge */
	eid = ce + nc;              /* per edge its new vertex (HM_NONE: not split) */
	touch = eid + nc;           /* per vertex: a split edge on it */
	/* the edges, as the uniform level numbers them (sorted by their vertices) */
	for (f = 0; f < nf; f++)
		for (j = m->fo[f]; j < m->fo[f + 1]; j++)
		{
			unsigned int a = m->fv[j], b = m->fv[j + 1 == m->fo[f + 1] ? m->fo[f] : j + 1];
			key[j].a = a < b ? a : b;
			key[j].b = a < b ? b : a;
			key[j].c = j;
		}
	qsort(key, nc, sizeof(sdv_key), sdv_key_cmp);
	for (i = 0; i < nc; i++)
	{
		if (i == 0 || key[i].a != key[i - 1].a || key[i].b != key[i - 1].b)
		{
			e = ne++;
			ev[2 * e] = key[i].a;
			ev[2 * e + 1] = key[i].b;
			ef[2 * e] = ef[2 * e + 1] = HM_NONE;
		}
		e = ne - 1;
		ce[key[i].c] = e;
		if (en[e] < 2)
			ef[2 * e + en[e]] = key[i].c / 3;
		en[e]++;
	}
	/* the faces at a turn */
	for (f = 0; f < nf; f++)
	{
		double const *a = m->X + (size_t)m->fv[3 * f] * kv, *b = m->X + (size_t)m->fv[3 * f + 1] * kv, *c = m->X + (size_t)m->fv[3 * f + 2] * kv;
		double u[3], w[3], *nn = nrm + 3 * (size_t)f, l;
		for (j = 0; j < 3; j++) { u[j] = b[j] - a[j]; w[j] = c[j] - a[j]; }
		nn[0] = u[1] * w[2] - u[2] * w[1];
		nn[1] = u[2] * w[0] - u[0] * w[2];
		nn[2] = u[0] * w[1] - u[1] * w[0];
		l = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
		for (j = 0; j < 3; j++) nn[j] = l > 0 ? nn[j] / l : 0;
	}
	for (e = 0; e < ne; e++)
		if (en[e] == 2)
		{
			double const *n0 = nrm + 3 * (size_t)ef[2 * e], *n1 = nrm + 3 * (size_t)ef[2 * e + 1];
			if (n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] < cosa)
				mark[ef[2 * e]] = mark[ef[2 * e + 1]] = 1;
		}
	for (t = 0; t < m->nfe; t++)
	{
		unsigned int x = sdv_find_edge(ev, ne, m->fe[2 * t], m->fe[2 * t + 1]);
		if (x != HM_NONE && m->fs[t] > 0)
			for (j = 0; j < en[x] && j < 2; j++) mark[ef[2 * x + j]] = 1;
	}
	for (e = 0; e < ne; e++) eid[e] = HM_NONE;
	for (f = 0; f < nf; f++)
		if (mark[f])
			for (j = 0; j < 3; j++) eid[ce[3 * f + j]] = 0;
	for (e = 0; e < ne; e++)
		if (eid[e] == 0)
		{
			eid[e] = nv + ns++;
			touch[ev[2 * e]] = touch[ev[2 * e + 1]] = 1;
		}
	for (f = 0; f < nf; f++)
	{
		unsigned int s = (eid[ce[3 * f]] != HM_NONE) + (eid[ce[3 * f + 1]] != HM_NONE) + (eid[ce[3 * f + 2]] != HM_NONE);
		nfo += s == 3 ? 4 : s + 1;
	}
	nco = 3 * nfo;
	for (t = 0; t < m->nfe; t++)
	{
		unsigned int x = sdv_find_edge(ev, ne, m->fe[2 * t], m->fe[2 * t + 1]);
		nfe2 += x != HM_NONE && eid[x] != HM_NONE ? 2 : 1;
	}
	o->nv = nv + ns;
	o->nf = nfo;
	o->nc = nco;
	o->nfe = nfe2;
	if (!sdv_alloc(o, k))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	/* the positions: the uniform level's where the mesh is split, the old ones elsewhere */
	for (i = 0; i < nv; i++)
	{
		memcpy(o->X + (size_t)i * kv, (touch[i] ? full.X : m->X) + (size_t)i * kv, sizeof(double) * kv);
		o->vk[i] = m->vk[i];
	}
	for (e = 0; e < ne; e++)
		if (eid[e] != HM_NONE)
			memcpy(o->X + (size_t)eid[e] * kv, full.X + ((size_t)nv + e) * kv, sizeof(double) * kv);
	/* the faces: by how many of their edges split (corner rows: an old corner's own, an edge point's the mean of its two) */
	{
		unsigned int g = 0, q = 0;
		for (f = 0; f < nf; f++)
		{
			unsigned int c0 = 3 * f, v[3], mid[3], s = 0, tri[4][3], ntr = 0, r, z;
			for (j = 0; j < 3; j++)
			{
				v[j] = m->fv[c0 + j];
				mid[j] = eid[ce[c0 + j]];
				s += mid[j] != HM_NONE;
			}
			if (s == 0)
			{
				tri[0][0] = 0; tri[0][1] = 1; tri[0][2] = 2;
				ntr = 1;
			}
			else if (s == 3)
			{
				/* (Loop's four: the corners' then the middle one) */
				static unsigned int const lv[4][3] = { { 0, 3, 5 }, { 1, 4, 3 }, { 2, 5, 4 }, { 3, 4, 5 } };
				for (r = 0; r < 4; r++)
					for (z = 0; z < 3; z++) tri[r][z] = lv[r][z];
				ntr = 4;
			}
			else if (s == 1)
			{
				for (j = 0; mid[j] == HM_NONE; j++)
					;
				/* edge j (corner j to j + 1) split at 3 + j */
				tri[0][0] = j; tri[0][1] = 3 + j; tri[0][2] = (j + 2) % 3;
				tri[1][0] = 3 + j; tri[1][1] = (j + 1) % 3; tri[1][2] = (j + 2) % 3;
				ntr = 2;
			}
			else
			{
				unsigned int a1, a2;
				for (j = 0; mid[j] != HM_NONE; j++)
					;
				/* edge j kept, j + 1 and j + 2 split */
				a1 = (j + 1) % 3;
				a2 = (j + 2) % 3;
				tri[0][0] = a2; tri[0][1] = 3 + a2; tri[0][2] = 3 + a1;
				tri[1][0] = j; tri[1][1] = a1; tri[1][2] = 3 + a1;
				tri[2][0] = j; tri[2][1] = 3 + a1; tri[2][2] = 3 + a2;
				ntr = 3;
			}
			for (r = 0; r < ntr; r++, g++)
			{
				o->fo[g] = q;
				memcpy(o->F + (size_t)g * kf, m->F + (size_t)f * kf, sizeof(double) * kf);
				for (z = 0; z < 3; z++, q++)
				{
					unsigned int id = tri[r][z];
					double* y = o->C + (size_t)q * kc;
					if (id < 3)
					{
						o->fv[q] = v[id];
						memcpy(y, m->C + ((size_t)c0 + id) * kc, sizeof(double) * kc);
					}
					else
					{
						unsigned int a = id - 3, b = (id - 2) % 3;
						o->fv[q] = mid[a];
						for (j = 0; j < kc; j++) y[j] = 0.5 * (m->C[((size_t)c0 + a) * kc + j] + m->C[((size_t)c0 + b) * kc + j]);
						sdv_rules(y, m->C + ((size_t)c0 + a) * kc, k, QAWS_MESH_CORNER);
					}
				}
			}
		}
		o->fo[g] = q;
	}
	/* the creases: a split one's halves (sharp one level less), the others as they were */
	o->nfe = 0;
	for (t = 0; t < m->nfe; t++)
	{
		unsigned int x = sdv_find_edge(ev, ne, m->fe[2 * t], m->fe[2 * t + 1]);
		if (x != HM_NONE && eid[x] != HM_NONE)
		{
			double s = m->fs[t] >= SDV_INF ? SDV_INF : m->fs[t] - 1;
			if (s <= 0)
				continue;
			o->fe[2 * o->nfe] = m->fe[2 * t];
			o->fe[2 * o->nfe + 1] = eid[x];
			o->fs[o->nfe++] = s;
			o->fe[2 * o->nfe] = eid[x];
			o->fe[2 * o->nfe + 1] = m->fe[2 * t + 1];
			o->fs[o->nfe++] = s;
		}
		else
		{
			o->fe[2 * o->nfe] = m->fe[2 * t];
			o->fe[2 * o->nfe + 1] = m->fe[2 * t + 1];
			o->fs[o->nfe++] = m->fs[t];
		}
	}
done:
	if (st != QAWS_STATUS_OK)
		sdv_free(o);
	sdv_free(&full);
	qaws_internal_free(key);
	qaws_internal_free(ui);
	qaws_internal_free(mark);
	qaws_internal_free(nrm);
	return st;
}

qaws_status qaws_mesh_3d_subdivide(qaws_mesh_3d** out_mesh, qaws_mesh_3d const* mesh, qaws_mesh_subdivide_desc const* desc)
{
	sdv_mesh m, o;
	sdv_cols k;
	sdv_attr* at = NULL;
	qaws_scalar* sv = NULL;
	double* ln = NULL;
	unsigned int* fsz = NULL;
	unsigned int nv, nf, nc, na, i, j, a, lv, levels;
	int loop;
	size_t big;
	qaws_mesh_3d* out = NULL;
	qaws_status st = QAWS_STATUS_OK;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh || !desc || (desc->scheme != QAWS_MESH_SUBDIVIDE_CATMULL_CLARK && desc->scheme != QAWS_MESH_SUBDIVIDE_LOOP &&
		desc->scheme != QAWS_MESH_SUBDIVIDE_BUTTERFLY && desc->scheme != QAWS_MESH_SUBDIVIDE_SQRT3) ||
		(desc->scheme == QAWS_MESH_SUBDIVIDE_BUTTERFLY && desc->limit_normals) || (desc->scheme == QAWS_MESH_SUBDIVIDE_SQRT3 && (desc->limit || desc->limit_normals)) ||
		!(desc->adaptive_angle >= 0) || (desc->adaptive_angle > 0 && (desc->scheme == QAWS_MESH_SUBDIVIDE_CATMULL_CLARK || desc->scheme == QAWS_MESH_SUBDIVIDE_SQRT3)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nf = qaws_mesh_3d_get_face_count(mesh);
	nc = qaws_mesh_3d_get_corner_count(mesh);
	na = qaws_mesh_3d_get_attribute_count(mesh);
	if (nf == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	loop = desc->scheme == QAWS_MESH_SUBDIVIDE_LOOP ? 1 : desc->scheme == QAWS_MESH_SUBDIVIDE_BUTTERFLY ? 2 : desc->scheme == QAWS_MESH_SUBDIVIDE_SQRT3 ? 3 : 0;
	levels = desc->levels ? desc->levels : 1;
	memset(&m, 0, sizeof m);
	memset(&o, 0, sizeof o);
	memset(&k, 0, sizeof k);

	/* the columns */
	at = (sdv_attr*)qaws_internal_scratch_zalloc((size_t)na + 1, sizeof(sdv_attr));
	if (!at)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	k.kv = 3;
	for (a = 0; a < na; a++)
	{
		qaws_attribute_desc d;
		qaws_mesh_3d_get_attribute_desc(&d, mesh, a);
		at[a].element = (unsigned int)d.element;
		at[a].rule = (unsigned int)d.rule;
		at[a].n = d.components;
		at[a].off = d.element == QAWS_MESH_VERTEX ? k.kv : (d.element == QAWS_MESH_CORNER ? k.kc : k.kf);
		*(d.element == QAWS_MESH_VERTEX ? &k.kv : (d.element == QAWS_MESH_CORNER ? &k.kc : &k.kf)) += d.components;
	}
	k.na = na;
	k.at = at;

	/* the input in rows */
	m.nv = nv;
	m.nf = nf;
	m.nc = nc;
	m.nfe = qaws_mesh_3d_get_feature_edge_count(mesh);
	big = (size_t)nv * k.kv + (size_t)nc * k.kc + (size_t)nf * k.kf + 3 * (size_t)nv;
	sv = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (big + 1));
	fsz = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)nf + 1));
	if (!sv || !fsz || !sdv_alloc(&m, &k))
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	qaws_mesh_3d_get_faces(m.fv, fsz, mesh);
	m.fo[0] = 0;
	for (i = 0; i < nf; i++)
	{
		if (loop && fsz[i] != 3)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		m.fo[i + 1] = m.fo[i] + fsz[i];
	}
	qaws_mesh_3d_get_feature_edges(m.fe, mesh);
	for (i = 0; i < m.nfe; i++)
	{
		double s = desc->sharpness ? (double)desc->sharpness[i] : SDV_INF;
		if (!(s >= 0))
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto done;
		}
		m.fs[i] = s;
	}
	qaws_mesh_3d_get_positions(sv, mesh);
	for (i = 0; i < nv; i++)
	{
		for (j = 0; j < 3; j++) m.X[(size_t)i * k.kv + j] = (double)sv[3 * (size_t)i + j];
		m.vk[i] = (unsigned char)(qaws_mesh_3d_is_feature_corner(mesh, i) && hm_feature_degree(&mesh->hm, i) != 1);
	}
	for (a = 0; a < na; a++)
	{
		unsigned int count = at[a].element == QAWS_MESH_VERTEX ? nv : (at[a].element == QAWS_MESH_CORNER ? nc : nf);
		unsigned int w = at[a].element == QAWS_MESH_VERTEX ? k.kv : (at[a].element == QAWS_MESH_CORNER ? k.kc : k.kf);
		double* rows = at[a].element == QAWS_MESH_VERTEX ? m.X : (at[a].element == QAWS_MESH_CORNER ? m.C : m.F);
		qaws_mesh_3d_get_attribute(sv, mesh, a);
		for (i = 0; i < count; i++)
			for (j = 0; j < at[a].n; j++) rows[(size_t)i * w + at[a].off + j] = (double)sv[(size_t)i * at[a].n + j];
	}
	qaws_internal_free(sv);
	qaws_internal_free(fsz);
	sv = NULL;
	fsz = NULL;

	for (lv = 0; lv < levels; lv++)
	{
		st = loop == 3 ? sdv_sqrt3(&o, &m, &k) : desc->adaptive_angle > 0 ? sdv_adaptive(&o, &m, loop, &k, (double)desc->adaptive_angle) : sdv_step(&o, &m, loop, &k, 0);
		sdv_free(&m);
		if (st != QAWS_STATUS_OK)
			goto done;
		m = o;
		memset(&o, 0, sizeof o);
	}
	if ((desc->limit || desc->limit_normals) && loop != 2)
	{
		st = sdv_step(&o, &m, loop, &k, 1);
		if (st == QAWS_STATUS_OK && desc->limit_normals)
		{
			ln = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)m.nv + 1));
			if (!ln || !sdv_limit_normals(ln, &m, o.X, loop, k.kv))
				st = QAWS_STATUS_ALLOCATION_FAILURE;
		}
		sdv_free(&m);
		if (st != QAWS_STATUS_OK)
			goto done;
		m = o;
		memset(&o, 0, sizeof o);
	}

	/* the result */
	{
		double* x = (double*)qaws_internal_scratch_malloc(sizeof(double) * (3 * (size_t)m.nv + 1));
		big = (size_t)m.nv * k.kv + (size_t)m.nc * k.kc + (size_t)m.nf * k.kf;
		sv = (qaws_scalar*)qaws_internal_scratch_malloc(sizeof(qaws_scalar) * (big + 1));
		fsz = (unsigned int*)qaws_internal_scratch_malloc(sizeof(unsigned int) * ((size_t)m.nf + 1));
		if (!x || !sv || !fsz)
		{
			qaws_internal_free(x);
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (i = 0; i < m.nv; i++)
			for (j = 0; j < 3; j++) x[3 * (size_t)i + j] = m.X[(size_t)i * k.kv + j];
		for (i = 0; i < m.nf; i++) fsz[i] = m.fo[i + 1] - m.fo[i];
		st = hm_mesh_3d_create(x, m.nv, m.fv, fsz, m.nf, &out);
		qaws_internal_free(x);
		if (st != QAWS_STATUS_OK)
			goto done;
		if (qaws_mesh_3d_get_vertex_count(out) != m.nv)
		{
			st = QAWS_STATUS_INTERNAL_ERROR;
			goto done;
		}
	}
	if (m.nfe > 0)
	{
		st = qaws_mesh_3d_set_feature_edges(out, m.fe, m.nfe);
		if (st != QAWS_STATUS_OK)
			goto done;
	}
	for (i = 0; i < m.nv; i++)
		if (m.vk[i])
			out->hm.vflag[i] |= HM_CORNER;
	for (a = 0; a < na && st == QAWS_STATUS_OK; a++)
	{
		qaws_attribute_desc d;
		unsigned int count = at[a].element == QAWS_MESH_VERTEX ? m.nv : (at[a].element == QAWS_MESH_CORNER ? m.nc : m.nf);
		unsigned int w = at[a].element == QAWS_MESH_VERTEX ? k.kv : (at[a].element == QAWS_MESH_CORNER ? k.kc : k.kf), id;
		double const* rows = at[a].element == QAWS_MESH_VERTEX ? m.X : (at[a].element == QAWS_MESH_CORNER ? m.C : m.F);
		qaws_mesh_3d_get_attribute_desc(&d, mesh, a);
		for (i = 0; i < count; i++)
			for (j = 0; j < at[a].n; j++) sv[(size_t)i * at[a].n + j] = (qaws_scalar)rows[(size_t)i * w + at[a].off + j];
		st = qaws_mesh_3d_add_attribute(&id, out, sv, &d);
	}
	if (ln && st == QAWS_STATUS_OK)
	{
		qaws_attribute_desc d;
		unsigned int id;
		d.name = "limit_normal";
		d.element = QAWS_MESH_VERTEX;
		d.components = 3;
		d.rule = QAWS_ATTRIBUTE_NORMALIZED;
		for (i = 0; i < 3 * m.nv; i++) sv[i] = (qaws_scalar)ln[i];
		st = qaws_mesh_3d_add_attribute(&id, out, sv, &d);
	}
done:
	if (st != QAWS_STATUS_OK)
		qaws_mesh_3d_destroy(out);
	else
		*out_mesh = out;
	sdv_free(&m);
	sdv_free(&o);
	qaws_internal_free(at);
	qaws_internal_free(sv);
	qaws_internal_free(ln);
	qaws_internal_free(fsz);
	return st;
}
