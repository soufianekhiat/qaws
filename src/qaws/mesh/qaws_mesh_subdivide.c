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

/* one level of m into o; limit: o is m with its vertices moved onto the limit surface instead */
static qaws_status sdv_step(sdv_mesh* o, sdv_mesh const* m, int loop, sdv_cols const* k, int limit)
{
	unsigned int nv = m->nv, nf = m->nf, nc = m->nc, kv = k->kv, kc = k->kc, kf = k->kf, ne = 0, i, j, f, e, t;
	size_t nv2, nf2, nc2;
	sdv_key* key = (sdv_key*)qaws_internal_scratch_malloc(sizeof(sdv_key) * ((size_t)nc + 1));
	unsigned int* ui = (unsigned int*)qaws_internal_scratch_zalloc(9 * (size_t)nc + 6 * (size_t)nv + m->nfe + 1, sizeof(unsigned int));
	double* dd = (double*)qaws_internal_scratch_zalloc((size_t)nc + (size_t)nv + 2 * (size_t)nv * kv + (loop ? 0 : (size_t)nf * kv) + 4 * (size_t)kv + 1,
		sizeof(double));
	unsigned int *ev, *ef, *eo, *en, *ce, *cf, *val, *nfc, *nsh, *sh2, *fidx;
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
		unsigned int a = m->fe[2 * t], b = m->fe[2 * t + 1], lo = 0, hi = ne;
		if (a > b) { unsigned int s = a; a = b; b = s; }
		fidx[t] = HM_NONE;
		while (lo < hi)
		{
			unsigned int mid = lo + (hi - lo) / 2;
			if (ev[2 * mid] < a || (ev[2 * mid] == a && ev[2 * mid + 1] < b)) lo = mid + 1;
			else hi = mid;
		}
		if (lo < ne && ev[2 * lo] == a && ev[2 * lo + 1] == b)
		{
			fidx[t] = lo;
			if (es[lo] < SDV_INF)
				es[lo] = m->fs[t];
		}
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
		if (es[e] < 1)
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
	return st;
}

qaws_status qaws_mesh_3d_subdivide(qaws_mesh_3d** out_mesh, qaws_mesh_3d const* mesh, qaws_mesh_subdivide_desc const* desc)
{
	sdv_mesh m, o;
	sdv_cols k;
	sdv_attr* at = NULL;
	qaws_scalar* sv = NULL;
	unsigned int* fsz = NULL;
	unsigned int nv, nf, nc, na, i, j, a, lv, levels;
	int loop;
	size_t big;
	qaws_mesh_3d* out = NULL;
	qaws_status st = QAWS_STATUS_OK;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh || !desc || (desc->scheme != QAWS_MESH_SUBDIVIDE_CATMULL_CLARK && desc->scheme != QAWS_MESH_SUBDIVIDE_LOOP))
		return QAWS_STATUS_INVALID_ARGUMENT;
	nv = qaws_mesh_3d_get_vertex_count(mesh);
	nf = qaws_mesh_3d_get_face_count(mesh);
	nc = qaws_mesh_3d_get_corner_count(mesh);
	na = qaws_mesh_3d_get_attribute_count(mesh);
	if (nf == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;
	loop = desc->scheme == QAWS_MESH_SUBDIVIDE_LOOP;
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
		st = sdv_step(&o, &m, loop, &k, 0);
		sdv_free(&m);
		if (st != QAWS_STATUS_OK)
			goto done;
		m = o;
		memset(&o, 0, sizeof o);
	}
	if (desc->limit)
	{
		st = sdv_step(&o, &m, loop, &k, 1);
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
done:
	if (st != QAWS_STATUS_OK)
		qaws_mesh_3d_destroy(out);
	else
		*out_mesh = out;
	sdv_free(&m);
	sdv_free(&o);
	qaws_internal_free(at);
	qaws_internal_free(sv);
	qaws_internal_free(fsz);
	return st;
}
