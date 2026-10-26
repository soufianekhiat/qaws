/*
 * The corner store: building from faces (input that is not a manifold split
 * into one), navigation, local operations (a patch of faces replaced by new
 * ones and re-stitched to its surroundings) and compaction.
 */

#include "qaws_mesh_core.h"
#include "../internal/qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64))
#include <xmmintrin.h>
#define HM_PREFETCH(p) _mm_prefetch((char const*)(p), _MM_HINT_T0)
#elif defined(__GNUC__) || defined(__clang__)
#define HM_PREFETCH(p) __builtin_prefetch(p)
#else
#define HM_PREFETCH(p) ((void)(p))
#endif

static int hm_re(void** p, size_t bytes)
{
	void* q = qaws_internal_realloc(*p, bytes ? bytes : 1);
	if (!q)
		return 0;
	*p = q;
	return 1;
}

int hm_reserve(qaws_hm* m, unsigned int dv, unsigned int dc, unsigned int df)
{
	unsigned int i, c;
	if (qaws_internal_add2(m->nv, dv) > m->capv)
	{
		c = qaws_internal_grow(m->capv, qaws_internal_add2(m->nv, dv), 16);
		if (!c)
			return 0;
		if (!hm_re((void**)&m->x, qaws_internal_mul3(sizeof(double), m->dim, c)) || !hm_re((void**)&m->vc, qaws_internal_mul2(sizeof(unsigned int), c))
			|| !hm_re((void**)&m->vsrc, qaws_internal_mul2(sizeof(unsigned int), c)) || !hm_re((void**)&m->vflag, c)
			|| !hm_re((void**)&m->hole[0], qaws_internal_mul2(sizeof(unsigned int), c))
			|| (m->foot && !hm_re((void**)&m->foot, qaws_internal_mul2(sizeof(qaws_mesh_location), c))))
			return 0;
		for (i = 0; i < m->na; i++)
			if (m->attr[i].element == QAWS_MESH_VERTEX && !hm_re((void**)&m->attr[i].v, qaws_internal_mul3(sizeof(double), m->attr[i].n, c)))
				return 0;
		m->capv = c;
	}
	if (qaws_internal_add2(m->nc, dc) > m->capc)
	{
		c = qaws_internal_grow(m->capc, qaws_internal_add2(m->nc, dc), 64);
		if (!c)
			return 0;
		if (!hm_re((void**)&m->cv, qaws_internal_mul2(sizeof(unsigned int), c)) || !hm_re((void**)&m->co, qaws_internal_mul2(sizeof(unsigned int), c))
			|| !hm_re((void**)&m->cf, qaws_internal_mul2(sizeof(unsigned int), c)) || !hm_re((void**)&m->cflag, c))
			return 0;
		for (i = 0; i < m->na; i++)
			if (HM_ON_CORNERS(m->attr[i].element) && !hm_re((void**)&m->attr[i].v, qaws_internal_mul3(sizeof(double), m->attr[i].n, c)))
				return 0;
		m->capc = c;
	}
	if (qaws_internal_add2(m->nf, df) > m->capf)
	{
		c = qaws_internal_grow(m->capf, qaws_internal_add2(m->nf, df), 16);
		if (!c)
			return 0;
		if (!hm_re((void**)&m->fc, qaws_internal_mul2(sizeof(unsigned int), c)) || !hm_re((void**)&m->fn, qaws_internal_mul2(sizeof(unsigned int), c))
			|| !hm_re((void**)&m->fflag, c) || !hm_re((void**)&m->hole[1], qaws_internal_mul2(sizeof(unsigned int), c)))
			return 0;
		for (i = 0; i < m->na; i++)
			if (m->attr[i].element == QAWS_MESH_FACE && !hm_re((void**)&m->attr[i].v, qaws_internal_mul3(sizeof(double), m->attr[i].n, c)))
				return 0;
		m->capf = c;
	}
	return 1;
}

void hm_free(qaws_hm* m)
{
	unsigned int i;
	hm_patch* p;
	if (!m)
		return;
	p = &m->pb;
	qaws_internal_free(m->x); qaws_internal_free(m->vc); qaws_internal_free(m->vsrc); qaws_internal_free(m->vflag);
	qaws_internal_free(m->cv); qaws_internal_free(m->co); qaws_internal_free(m->cf); qaws_internal_free(m->cflag);
	qaws_internal_free(m->fc); qaws_internal_free(m->fn); qaws_internal_free(m->fflag);
	for (i = 0; i < m->na; i++)
		qaws_internal_free(m->attr[i].v);
	qaws_internal_free(m->attr);
	qaws_internal_free(m->foot);
	qaws_internal_free(m->hole[0]); qaws_internal_free(m->hole[1]);
	qaws_internal_free(p->old); qaws_internal_free(p->off); qaws_internal_free(p->fs); qaws_internal_free(p->fa); qaws_internal_free(p->fv); qaws_internal_free(p->s0); qaws_internal_free(p->sn); qaws_internal_free(p->hint); qaws_internal_free(p->esrc);
	qaws_internal_free(p->sp); qaws_internal_free(p->sw);
	memset(m, 0, sizeof *m);
}

double* hm_attr_at(qaws_hm const* m, unsigned int a, unsigned int i)
{
	return m->attr[a].v + (size_t)m->attr[a].n * i;
}

unsigned int hm_new_vertex(qaws_hm* m, double const* x)
{
	unsigned int v = m->nv++, k;
	for (k = 0; k < m->dim; k++)
		m->x[(size_t)m->dim * v + k] = x[k];
	m->vc[v] = HM_NONE;
	m->vsrc[v] = HM_NONE;
	m->vflag[v] = 0;
	if (m->foot)
	{
		m->foot[v].patch = HM_NONE;
		m->foot[v].a = m->foot[v].b = 0;
	}
	for (k = 0; k < m->na; k++)
		if (m->attr[k].element == QAWS_MESH_VERTEX)
			memset(hm_attr_at(m, k, v), 0, sizeof(double) * m->attr[k].n);
	return v;
}

void hm_kill_vertex(qaws_hm* m, unsigned int v)
{
	m->vflag[v] |= HM_DELETED;
	m->vc[v] = HM_NONE;
	m->hole[0][m->nhole[0]++] = v;
}

int hm_rule_unit(unsigned int rule)
{
	return rule == QAWS_ATTRIBUTE_NORMALIZED || rule == QAWS_ATTRIBUTE_NORMAL || rule == QAWS_ATTRIBUTE_TANGENT;
}

/* the sRGB transfer function and its inverse (odd beyond [0, 1]) */
static double hm_srgb_decode(double c)
{
	double a = fabs(c), l = a <= 0.04045 ? a / 12.92 : pow((a + 0.055) / 1.055, 2.4);
	return c < 0 ? -l : l;
}

static double hm_srgb_encode(double l)
{
	double a = fabs(l), c = a <= 0.0031308 ? 12.92 * a : 1.055 * pow(a, 1 / 2.4) - 0.055;
	return l < 0 ? -c : c;
}

/* Ottosson's linear sRGB to LMS, and the cube roots of LMS to Oklab; back by their exact inverses (the published
   inverses are rounded: a colour would come back off by 1e-8) */
static double const hm_lms[9] = { 0.4122214708, 0.5363325363, 0.0514459929, 0.2119034982, 0.6806995451, 0.1073969566,
	0.0883024619, 0.2817188376, 0.6299787005 };
static double const hm_lab[9] = { 0.2104542553, 0.7936177850, -0.0040720468, 1.9779984951, -2.4285922050, 0.4505937099,
	0.0259040371, 0.7827717662, -0.8086757660 };

static void hm_mul3(double const* M, double* v)
{
	double x = v[0], y = v[1], z = v[2];
	v[0] = M[0] * x + M[1] * y + M[2] * z;
	v[1] = M[3] * x + M[4] * y + M[5] * z;
	v[2] = M[6] * x + M[7] * y + M[8] * z;
}

/* v = M^-1 v (Cramer) */
static void hm_solve3(double const* M, double* v)
{
	double c[9], d;
	unsigned int i;
	c[0] = M[4] * M[8] - M[5] * M[7]; c[1] = M[2] * M[7] - M[1] * M[8]; c[2] = M[1] * M[5] - M[2] * M[4];
	c[3] = M[5] * M[6] - M[3] * M[8]; c[4] = M[0] * M[8] - M[2] * M[6]; c[5] = M[2] * M[3] - M[0] * M[5];
	c[6] = M[3] * M[7] - M[4] * M[6]; c[7] = M[1] * M[6] - M[0] * M[7]; c[8] = M[0] * M[4] - M[1] * M[3];
	d = M[0] * c[0] + M[1] * c[3] + M[2] * c[6];
	for (i = 0; i < 9; i++)
		c[i] /= d;
	hm_mul3(c, v);
}

void hm_space_in(double* v, unsigned int n, unsigned int rule)
{
	unsigned int k;
	if (rule != QAWS_ATTRIBUTE_SRGB && rule != QAWS_ATTRIBUTE_OKLAB)
		return;
	for (k = 0; k < 3; k++)
		v[k] = hm_srgb_decode(v[k]);
	if (rule == QAWS_ATTRIBUTE_OKLAB)
	{
		hm_mul3(hm_lms, v);
		for (k = 0; k < 3; k++)
			v[k] = cbrt(v[k]);
		hm_mul3(hm_lab, v);
	}
	if (n > 3)
		for (k = 0; k < 3; k++)
			v[k] *= v[3];
}

void hm_space_out(double* v, unsigned int n, unsigned int rule)
{
	unsigned int k;
	if (hm_rule_unit(rule))
	{
		/* (a tangent: its direction made unit, its handedness a sign) */
		unsigned int d = rule == QAWS_ATTRIBUTE_TANGENT ? 3 : n;
		double l = 0;
		for (k = 0; k < d; k++)
			l += v[k] * v[k];
		if (l > 0)
			for (l = 1 / sqrt(l), k = 0; k < d; k++)
				v[k] *= l;
		if (rule == QAWS_ATTRIBUTE_TANGENT)
			v[3] = v[3] < 0 ? -1 : 1;
		return;
	}
	if (rule != QAWS_ATTRIBUTE_SRGB && rule != QAWS_ATTRIBUTE_OKLAB)
		return;
	if (n > 3 && v[3] != 0)
		for (k = 0; k < 3; k++)
			v[k] /= v[3];
	if (rule == QAWS_ATTRIBUTE_OKLAB)
	{
		hm_solve3(hm_lab, v);
		for (k = 0; k < 3; k++)
			v[k] = v[k] * v[k] * v[k];
		hm_solve3(hm_lms, v);
	}
	for (k = 0; k < 3; k++)
		v[k] = hm_srgb_encode(v[k]);
}

void hm_lerp(double* dst, double const* a, double const* b, double t, unsigned int n, unsigned int rule)
{
	double x[16], y[16];
	unsigned int k;
	if (rule == QAWS_ATTRIBUTE_NEAREST)
	{
		memmove(dst, t < 0.5 ? a : b, sizeof(double) * n);
		return;
	}
	memcpy(x, a, sizeof(double) * n);
	memcpy(y, b, sizeof(double) * n);
	hm_space_in(x, n, rule);
	hm_space_in(y, n, rule);
	for (k = 0; k < n; k++)
		x[k] = (1 - t) * x[k] + t * y[k];
	hm_space_out(x, n, rule);
	memcpy(dst, x, sizeof(double) * n);
}

/* the weighted mean of s values (vals + stride * pos[i]) by rule; NEAREST: the heaviest */
static void hm_blend(double* out, double const* vals, unsigned int stride, unsigned int const* pos, double const* w,
	unsigned int s, unsigned int n, unsigned int rule)
{
	unsigned int i, j;
	double best = -1, acc[16], v[16];
	memset(acc, 0, sizeof(double) * n);
	for (i = 0; i < s; i++)
	{
		memcpy(v, vals + (size_t)stride * pos[i], sizeof(double) * n);
		if (rule == QAWS_ATTRIBUTE_NEAREST)
		{
			if (w[i] > best)
			{
				best = w[i];
				memcpy(acc, v, sizeof(double) * n);
			}
			continue;
		}
		hm_space_in(v, n, rule);
		for (j = 0; j < n; j++)
			acc[j] += w[i] * v[j];
	}
	if (rule != QAWS_ATTRIBUTE_NEAREST)
		hm_space_out(acc, n, rule);
	memcpy(out, acc, sizeof(double) * n);
}

/* ======================================================================== */
/*  Navigation                                                              */
/* ======================================================================== */

unsigned int hm_next(qaws_hm const* m, unsigned int c)
{
	unsigned int f = m->cf[c];
	return c + 1 < m->fc[f] + m->fn[f] ? c + 1 : m->fc[f];
}

unsigned int hm_prev(qaws_hm const* m, unsigned int c)
{
	unsigned int f = m->cf[c];
	return c > m->fc[f] ? c - 1 : m->fc[f] + m->fn[f] - 1;
}

unsigned int hm_to(qaws_hm const* m, unsigned int c)
{
	return m->cv[hm_next(m, c)];
}

unsigned int hm_ccw(qaws_hm const* m, unsigned int c)
{
	return m->co[hm_prev(m, c)];
}

unsigned int hm_cw(qaws_hm const* m, unsigned int c)
{
	unsigned int o = m->co[c];
	return o == HM_NONE ? HM_NONE : hm_next(m, o);
}

/* the corner after c round vertex v counter-clockwise, HM_NONE past the last */
unsigned int hm_round(qaws_hm const* m, unsigned int v, unsigned int c)
{
	unsigned int n = hm_ccw(m, c);
	return n == m->vc[v] ? HM_NONE : n;
}

unsigned int hm_find(qaws_hm const* m, unsigned int a, unsigned int b)
{
	unsigned int c;
	for (c = m->vc[a]; c != HM_NONE; c = hm_round(m, a, c))
		if (hm_to(m, c) == b)
			return c;
	return HM_NONE;
}

unsigned int hm_find_edge(qaws_hm const* m, unsigned int a, unsigned int b)
{
	unsigned int c = hm_find(m, a, b);
	return c != HM_NONE ? c : hm_find(m, b, a);
}

int hm_is_border_vertex(qaws_hm const* m, unsigned int v)
{
	return m->vc[v] != HM_NONE && m->co[m->vc[v]] == HM_NONE;
}

unsigned int hm_valence(qaws_hm const* m, unsigned int v)
{
	unsigned int c, n = 0;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		n++;
	return n + (unsigned int)hm_is_border_vertex(m, v);
}

unsigned int hm_ring(qaws_hm const* m, unsigned int v, unsigned int* out, unsigned int cap)
{
	unsigned int c, last = HM_NONE, n = 0;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
	{
		if (n < cap)
			out[n] = hm_to(m, c);
		n++;
		last = c;
	}
	if (last != HM_NONE && hm_is_border_vertex(m, v))
	{
		/* the neighbour across the border edge into v */
		if (n < cap)
			out[n] = m->cv[hm_prev(m, last)];
		n++;
	}
	return n;
}

void hm_fix_vertex(qaws_hm* m, unsigned int v, unsigned int c)
{
	unsigned int k = c, n = c;
	if (c == HM_NONE)
	{
		m->vc[v] = HM_NONE;
		return;
	}
	for (;;)
	{
		n = hm_cw(m, k);
		if (n == HM_NONE || n == c)
			break;
		k = n;
	}
	m->vc[v] = n == HM_NONE ? k : c;
}

double hm_length(qaws_hm const* m, unsigned int c)
{
	double const *a = m->x + (size_t)m->dim * m->cv[c], *b = m->x + (size_t)m->dim * hm_to(m, c);
	double s = 0;
	unsigned int k;
	for (k = 0; k < m->dim; k++)
		s += (b[k] - a[k]) * (b[k] - a[k]);
	return sqrt(s);
}

void hm_face_normal(qaws_hm const* m, unsigned int f, double* out)
{
	unsigned int c, e = m->fc[f] + m->fn[f];
	out[0] = out[1] = out[2] = 0;
	for (c = m->fc[f]; c < e; c++)
	{
		double const *p = m->x + 3 * (size_t)m->cv[c], *q = m->x + 3 * (size_t)m->cv[c + 1 < e ? c + 1 : m->fc[f]];
		out[0] += (p[1] - q[1]) * (p[2] + q[2]);
		out[1] += (p[2] - q[2]) * (p[0] + q[0]);
		out[2] += (p[0] - q[0]) * (p[1] + q[1]);
	}
	out[0] *= 0.5; out[1] *= 0.5; out[2] *= 0.5;
}

unsigned int hm_edge_count(qaws_hm const* m)
{
	unsigned int c, n = 0;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
			n++;
	return n;
}

/* ======================================================================== */
/*  Building                                                                */
/* ======================================================================== */

/* a half-edge in the table: its corner and the vertices it leaves and reaches (read without going through the corner) */
typedef struct hm_slot
{
	unsigned int v, w, c;
} hm_slot;

/* both directions of an edge hash alike, so one probe finds both */
static unsigned int hm_hkey(unsigned int a, unsigned int b)
{
	unsigned long long k = (unsigned long long)(a < b ? a : b) << 32 | (a < b ? b : a);
	k ^= k >> 33;
	k *= 0xff51afd7ed558ccdull;
	k ^= k >> 33;
	return (unsigned int)k;
}

/* the corners of half-edges a -> b and b -> a in the table, with the free slot where a -> b would go */
static void hm_lookup(hm_slot const* slot, unsigned int mask, unsigned int a, unsigned int b, unsigned int* out_ab, unsigned int* out_ba,
	unsigned int* out_pos)
{
	unsigned int i = hm_hkey(a, b) & mask;
	*out_ab = *out_ba = HM_NONE;
	for (; slot[i].c != HM_NONE; i = (i + 1) & mask)
	{
		if (slot[i].v == a && slot[i].w == b)
			*out_ab = slot[i].c;
		else if (slot[i].v == b && slot[i].w == a)
			*out_ba = slot[i].c;
	}
	*out_pos = i;
}

/* every face's vertices in range and distinct */
static qaws_status hm_faces_ok(unsigned int const* fv, unsigned int const* fs, unsigned int nf, unsigned int nv)
{
	unsigned int f, i, off, *stamp = (unsigned int*)qaws_internal_scratch_zalloc(nv ? nv : 1, sizeof(unsigned int));
	if (!stamp)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (f = 0, off = 0; f < nf; f++)
	{
		unsigned int s = fs ? fs[f] : 3;
		for (i = 0; i < s; i++)
		{
			unsigned int v = fv[off + i];
			if (v >= nv || stamp[v] == f + 1)
			{
				qaws_internal_free(stamp);
				return QAWS_STATUS_INVALID_ARGUMENT;
			}
			stamp[v] = f + 1;
		}
		off += s;
	}
	qaws_internal_free(stamp);
	return QAWS_STATUS_OK;
}

/* each vertex's fans: a vertex on several gets a copy per fan past the first */
static qaws_status hm_fans(qaws_hm* m)
{
	unsigned int v, c, j, n = m->nv;
	unsigned int* vn = (unsigned int*)qaws_internal_scratch_zalloc((size_t)n + 1, sizeof(unsigned int));
	unsigned int* vl = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nc + 1)));
	unsigned char* seen = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nc + 1, 1);
	qaws_status st = QAWS_STATUS_OK;
	if (!vn || !vl || !seen)
		goto oom;
	for (c = 0; c < m->nc; c++)
		vn[m->cv[c] + 1]++;
	for (v = 0; v < n; v++)
		vn[v + 1] += vn[v];
	for (c = 0; c < m->nc; c++)
		vl[vn[m->cv[c]]++] = c;
	for (v = n; v > 0; v--)
		vn[v] = vn[v - 1];
	vn[0] = 0;
	for (v = 0; v < n; v++)
	{
		unsigned int fans = 0;
		for (j = vn[v]; j < vn[v + 1]; j++)
		{
			unsigned int k0 = vl[j], k, nk = HM_NONE, start, w = v;
			if (seen[k0])
				continue;
			if (fans++)
			{
				if (!hm_reserve(m, 1, 0, 0))
					goto oom;
				w = hm_new_vertex(m, m->x + (size_t)m->dim * v);
				m->vsrc[w] = m->vsrc[v];
			}
			/* the fan from its clockwise end, counter-clockwise (navigation
			   does not read the corners' vertices: renaming on the way is safe) */
			for (k = k0;;)
			{
				nk = hm_cw(m, k);
				if (nk == HM_NONE || nk == k0)
					break;
				k = nk;
			}
			start = nk == HM_NONE ? k : k0;
			for (k = start;;)
			{
				seen[k] = 1;
				m->cv[k] = w;
				nk = hm_ccw(m, k);
				if (nk == HM_NONE || nk == start)
					break;
				k = nk;
			}
			m->vc[w] = start;
		}
	}
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(vn); qaws_internal_free(vl); qaws_internal_free(seen);
	return st;
}

qaws_status hm_build(qaws_hm* m, unsigned int dim, double const* x, unsigned int nv, unsigned int const* fv,
	unsigned int const* fs, unsigned int nf)
{
	size_t total = 0;
	unsigned int f, i, off, mask, smax = 3, fa = 0, offa = 0;
	unsigned int *rep = NULL, *hab, *hba, *hpos;
	int marks;
	hm_slot* slot = NULL;
	qaws_status st = QAWS_STATUS_OK;
	memset(m, 0, sizeof *m);
	m->dim = dim;
	for (f = 0; f < nf; f++)
	{
		unsigned int s = fs ? fs[f] : 3;
		if (s < 3)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (s > smax)
			smax = s;
		total += s;
	}
	if (total > 0x3FFFFFFFu)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = hm_faces_ok(fv, fs, nf, nv);
	if (st != QAWS_STATUS_OK)
		return st;
	for (mask = 1; mask < 2 * total; mask <<= 1) {}
	slot = (hm_slot*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(hm_slot), mask));
	rep = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 4, smax));
	if (!slot || !rep || !hm_reserve(m, nv, (unsigned int)total, nf))
		goto oom;
	hab = rep + smax;   /* (a face's lookups: the corners of a -> b and b -> a, the free slot) */
	hba = hab + smax;
	hpos = hba + smax;
	memset(slot, 0xFF, sizeof(hm_slot) * mask);
	mask -= 1;
	for (i = 0; i < nv; i++)
	{
		hm_new_vertex(m, x + (size_t)dim * i);
		m->vsrc[i] = i;
	}
	for (f = 0, off = 0; f < nf; f++)
	{
		unsigned int s = fs ? fs[f] : 3, first = m->nc;
		/* the table is read at random: its slots for a face ahead are fetched while this one is linked */
		for (; fa < nf && fa < f + 8; fa++)
		{
			unsigned int sa = fs ? fs[fa] : 3;
			for (i = 0; i < sa; i++)
				HM_PREFETCH(&slot[hm_hkey(fv[offa + i], fv[offa + (i + 1) % sa]) & mask]);
			offa += sa;
		}
		m->fc[f] = first;
		m->fn[f] = s;
		m->fflag[f] = 0;
		for (i = 0; i < s; i++)
		{
			m->cv[first + i] = fv[off + i];
			m->co[first + i] = HM_NONE;
			m->cf[first + i] = f;
			m->cflag[first + i] = 0;
		}
		off += s;
		m->nc += s;
		m->nf++;
		/* an edge used the same way already, or by two faces already: this
		   face gets its own copies of the edge's ends */
		marks = 0;
		for (i = 0; i < s; i++)
		{
			unsigned int a = m->cv[first + i], b = m->cv[first + (i + 1) % s];
			hm_lookup(slot, mask, a, b, &hab[i], &hba[i], &hpos[i]);
			if (hab[i] != HM_NONE || (hba[i] != HM_NONE && m->co[hba[i]] != HM_NONE))
			{
				m->cflag[first + i] |= HM_MARK;
				marks = 1;
			}
		}
		for (i = 0; marks && i < s; i++)
		{
			unsigned int src = m->cv[first + i];
			rep[i] = src;
			if ((m->cflag[first + i] | m->cflag[first + (i + s - 1) % s]) & HM_MARK)
			{
				if (!hm_reserve(m, 1, 0, 0))
					goto oom;
				rep[i] = hm_new_vertex(m, m->x + (size_t)dim * src);
				m->vsrc[rep[i]] = m->vsrc[src];
			}
		}
		for (i = 0; marks && i < s; i++)
		{
			m->cv[first + i] = rep[i];
			m->cflag[first + i] = 0;
		}
		for (i = 0; i < s; i++)
		{
			unsigned int a = m->cv[first + i], b = m->cv[first + (i + 1) % s], pos, ab, o;
			if (marks)
				hm_lookup(slot, mask, a, b, &ab, &o, &pos);
			else
			{
				/* the face's vertices as looked up: the same answers (its edges are distinct), the free slot past
				   those its earlier edges took */
				ab = hab[i];
				o = hba[i];
				for (pos = hpos[i]; slot[pos].c != HM_NONE; pos = (pos + 1) & mask) {}
			}
			if (ab != HM_NONE)
				continue;   /* (not after the copies) */
			slot[pos].v = a;
			slot[pos].w = b;
			slot[pos].c = first + i;
			if (o != HM_NONE && m->co[o] == HM_NONE)
			{
				m->co[first + i] = o;
				m->co[o] = first + i;
			}
		}
	}
	st = hm_fans(m);
	goto done;
oom:
	st = QAWS_STATUS_ALLOCATION_FAILURE;
done:
	qaws_internal_free(slot); qaws_internal_free(rep);
	if (st != QAWS_STATUS_OK)
		hm_free(m);
	return st;
}

qaws_status hm_build_paired(qaws_hm* m, unsigned int dim, double const* x, unsigned int nv, unsigned int const* tri,
	unsigned int const* opposite, unsigned int nt)
{
	unsigned int i;
	qaws_status st;
	memset(m, 0, sizeof *m);
	m->dim = dim;
	if (3 * (size_t)nt > 0x3FFFFFFFu)
		return QAWS_STATUS_INVALID_ARGUMENT;
	st = hm_faces_ok(tri, NULL, nt, nv);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!hm_reserve(m, nv, 3 * nt, nt))
	{
		hm_free(m);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < nv; i++)
	{
		hm_new_vertex(m, x + (size_t)dim * i);
		m->vsrc[i] = i;
	}
	for (i = 0; i < nt; i++)
	{
		m->fc[i] = 3 * i;
		m->fn[i] = 3;
		m->fflag[i] = 0;
	}
	for (i = 0; i < 3 * nt; i++)
	{
		m->cv[i] = tri[i];
		m->co[i] = opposite[i];
		m->cf[i] = i / 3;
		m->cflag[i] = 0;
	}
	m->nc = 3 * nt;
	m->nf = nt;
	st = hm_fans(m);
	if (st != QAWS_STATUS_OK)
		hm_free(m);
	return st;
}

qaws_status hm_copy(qaws_hm* d, qaws_hm const* s)
{
	unsigned int i;
	memset(d, 0, sizeof *d);
	d->dim = s->dim;
	if (s->foot && !hm_foot_enable(d))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < s->na; i++)
		if (hm_attr_add(d, s->attr[i].name, s->attr[i].element, s->attr[i].rule, s->attr[i].n) == HM_NONE)
		{
			hm_free(d);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
	if (!hm_reserve(d, s->nv, s->nc, s->nf))
	{
		hm_free(d);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	d->nv = s->nv; d->nc = s->nc; d->nf = s->nf; d->cfree = s->cfree;
	memcpy(d->x, s->x, sizeof(double) * s->dim * s->nv);
	memcpy(d->vc, s->vc, sizeof(unsigned int) * s->nv);
	memcpy(d->vsrc, s->vsrc, sizeof(unsigned int) * s->nv);
	memcpy(d->vflag, s->vflag, s->nv);
	if (s->foot)
		memcpy(d->foot, s->foot, sizeof(qaws_mesh_location) * s->nv);
	memcpy(d->cv, s->cv, sizeof(unsigned int) * s->nc);
	memcpy(d->co, s->co, sizeof(unsigned int) * s->nc);
	memcpy(d->cf, s->cf, sizeof(unsigned int) * s->nc);
	memcpy(d->cflag, s->cflag, s->nc);
	memcpy(d->fc, s->fc, sizeof(unsigned int) * s->nf);
	memcpy(d->fn, s->fn, sizeof(unsigned int) * s->nf);
	memcpy(d->fflag, s->fflag, s->nf);
	for (i = 0; i < s->na; i++)
	{
		unsigned int count = s->attr[i].element == QAWS_MESH_VERTEX ? s->nv : HM_ON_CORNERS(s->attr[i].element) ? s->nc : s->nf;
		memcpy(d->attr[i].v, s->attr[i].v, sizeof(double) * s->attr[i].n * count);
	}
	for (i = 0; i < 2; i++)
	{
		d->nhole[i] = s->nhole[i];
		if (s->nhole[i])
			memcpy(d->hole[i], s->hole[i], sizeof(unsigned int) * s->nhole[i]);
	}
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Checking                                                                */
/* ======================================================================== */

qaws_status hm_check(qaws_hm const* m)
{
	unsigned int c, v, f, nfree = 0, *deg;
	qaws_status st = QAWS_STATUS_INTERNAL_ERROR;
	deg = (unsigned int*)qaws_internal_zalloc((size_t)m->nv + 1, sizeof(unsigned int));
	if (!deg)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (f = 0; f < m->nf; f++)
	{
		if (m->fflag[f] & HM_DELETED)
			continue;
		if (m->fn[f] < 3 || m->fc[f] + m->fn[f] > m->nc)
			goto done;
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			if (m->cf[c] != f)
				goto done;
	}
	for (c = 0; c < m->nc; c++)
	{
		unsigned int o = m->co[c], fa = m->cf[c];
		if (fa == HM_NONE)
		{
			nfree++;
			continue;
		}
		if (fa >= m->nf || (m->fflag[fa] & HM_DELETED) || c < m->fc[fa] || c >= m->fc[fa] + m->fn[fa])
			goto done;
		if (m->cv[c] >= m->nv || (m->vflag[m->cv[c]] & HM_DELETED) || m->cv[c] == hm_to(m, c))
			goto done;
		if (o != HM_NONE)
		{
			if (o >= m->nc || m->cf[o] == HM_NONE || m->co[o] != c || m->cv[o] != hm_to(m, c) || hm_to(m, o) != m->cv[c])
				goto done;
			if ((m->cflag[o] ^ m->cflag[c]) & HM_FEATURE)
				goto done;
		}
		deg[m->cv[c]]++;
	}
	if (nfree != m->cfree)
		goto done;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int n = 0, k, nk;
		if (m->vflag[v] & HM_DELETED)
			continue;
		c = m->vc[v];
		if (c == HM_NONE)
		{
			if (deg[v])
				goto done;
			continue;
		}
		if (c >= m->nc || m->cf[c] == HM_NONE || m->cv[c] != v)
			goto done;
		/* round the vertex: a border start ends at the border, an interior one closes */
		for (k = c;;)
		{
			if (m->cv[k] != v || ++n > deg[v])
				goto done;
			nk = hm_ccw(m, k);
			if (nk == HM_NONE)
			{
				if (m->co[c] != HM_NONE)
					goto done;
				break;
			}
			if (nk == c)
				break;
			k = nk;
		}
		if (n != deg[v])
			goto done;
	}
	st = QAWS_STATUS_OK;
done:
	qaws_internal_free(deg);
	return st;
}

/* ======================================================================== */
/*  Patches: the operations                                                 */
/* ======================================================================== */

static int hm_grow(void** p, unsigned int* cap, size_t need, size_t size)
{
	unsigned int c;
	if (need <= *cap)
		return 1;
	c = qaws_internal_grow(*cap, need, 16);
	if (!c || !hm_re(p, qaws_internal_mul2(c, size)))
		return 0;
	*cap = c;
	return 1;
}

static void pb_reset(qaws_hm* m)
{
	hm_patch* p = &m->pb;
	p->nold = p->nnew = p->ncv = p->nsrc = p->total = 0;
	p->failed = 0;
}

/* face f into the patch; its index there */
static unsigned int pb_old(qaws_hm* m, unsigned int f)
{
	hm_patch* p = &m->pb;
	unsigned int c = p->cold;
	if (p->failed || !hm_grow((void**)&p->old, &c, (size_t)p->nold + 1, sizeof(unsigned int))
		|| !hm_grow((void**)&p->off, &p->cold, (size_t)p->nold + 1, sizeof(unsigned int)))
	{
		p->failed = 1;
		return 0;
	}
	p->old[p->nold] = f;
	p->off[p->nold] = p->total;
	p->total += m->fn[f];
	return p->nold++;
}

/* a new face taking old face k's face attributes (HM_NONE: zeros) */
static void pb_face(qaws_hm* m, unsigned int k)
{
	hm_patch* p = &m->pb;
	unsigned int c = p->cnew;
	if (p->failed || !hm_grow((void**)&p->fs, &c, (size_t)p->nnew + 1, sizeof(unsigned int))
		|| !hm_grow((void**)&p->fa, &p->cnew, (size_t)p->nnew + 1, sizeof(unsigned int)))
	{
		p->failed = 1;
		return;
	}
	p->fs[p->nnew] = 0;
	p->fa[p->nnew] = k;
	p->nnew++;
}

/* a corner of the last new face at v, its edge's flags or'ed with hint */
static void pb_corner(qaws_hm* m, unsigned int v, unsigned char hint)
{
	hm_patch* p = &m->pb;
	unsigned int c1 = p->ccv, c2 = p->ccv, c3 = p->ccv, c4 = p->ccv;
	if (p->failed || !p->nnew || !hm_grow((void**)&p->fv, &c1, (size_t)p->ncv + 1, sizeof(unsigned int))
		|| !hm_grow((void**)&p->s0, &c2, (size_t)p->ncv + 1, sizeof(unsigned int)) || !hm_grow((void**)&p->sn, &c3, (size_t)p->ncv + 1, sizeof(unsigned int))
		|| !hm_grow((void**)&p->esrc, &c4, (size_t)p->ncv + 1, sizeof(unsigned int))
		|| !hm_grow((void**)&p->hint, &p->ccv, (size_t)p->ncv + 1, 1))
	{
		p->failed = 1;
		return;
	}
	p->fv[p->ncv] = v;
	p->s0[p->ncv] = p->nsrc;
	p->sn[p->ncv] = 0;
	p->hint[p->ncv] = hint;
	p->esrc[p->ncv] = HM_NONE;
	p->ncv++;
	p->fs[p->nnew - 1]++;
}

/* the last corner's edge a piece of the edge of the old corner at patch position pos */
static void pb_edge(qaws_hm* m, unsigned int pos)
{
	if (!m->pb.failed && m->pb.ncv)
		m->pb.esrc[m->pb.ncv - 1] = pos;
}

/* an attribute source of the last corner: an old corner's position, weighted */
static void pb_src(qaws_hm* m, unsigned int pos, double w)
{
	hm_patch* p = &m->pb;
	unsigned int c = p->csrc;
	if (p->failed || pos == HM_NONE || !hm_grow((void**)&p->sp, &c, (size_t)p->nsrc + 1, sizeof(unsigned int))
		|| !hm_grow((void**)&p->sw, &p->csrc, (size_t)p->nsrc + 1, sizeof(double)))
	{
		p->failed = 1;
		return;
	}
	p->sp[p->nsrc] = pos;
	p->sw[p->nsrc] = w;
	p->nsrc++;
	p->sn[p->ncv - 1]++;
}

/* the patch position of old face k's corner c */
static unsigned int pb_cpos(qaws_hm const* m, unsigned int k, unsigned int c)
{
	return m->pb.off[k] + (c - m->fc[m->pb.old[k]]);
}

/* the patch position of old face k's corner at vertex v */
static unsigned int pb_pos(qaws_hm const* m, unsigned int k, unsigned int v)
{
	unsigned int f = m->pb.old[k], j;
	for (j = 0; j < m->fn[f]; j++)
		if (m->cv[m->fc[f] + j] == v)
			return m->pb.off[k] + j;
	return HM_NONE;
}

/* a new corner at v copying old face k's corner there */
static void pb_copy(qaws_hm* m, unsigned int v, unsigned int k, unsigned char hint)
{
	pb_corner(m, v, hint);
	pb_src(m, pb_pos(m, k, v), 1);
}

/* a live corner at v among the patch's new corners and outside neighbours */
static unsigned int hm_any_corner(qaws_hm const* m, unsigned int v, unsigned int const* nc, unsigned int nnc,
	unsigned int const* otw, unsigned int nt)
{
	unsigned int i, o;
	if (m->vc[v] != HM_NONE && m->cf[m->vc[v]] != HM_NONE)
		return m->vc[v];
	for (i = 0; i < nnc; i++)
		if (m->cv[nc[i]] == v)
			return nc[i];
	for (i = 0; i < nt; i++)
		if ((o = otw[i]) != HM_NONE)
		{
			if (m->cv[o] == v)
				return o;
			if (m->cv[hm_next(m, o)] == v)
				return hm_next(m, o);
		}
	return HM_NONE;
}

/* The old faces replaced by the new ones, old vertex `from` read as `to` (a
   collapse): new corners stitched to each other and to the old faces'
   outside neighbours by their edges (an outside neighbour left alone is on
   the border now), feature flags kept by edge, corner attributes blended
   from their sources, the touched vertices' corners fixed. */
static int hm_replace(qaws_hm* m, unsigned int from, unsigned int to)
{
	hm_patch* p = &m->pb;
	unsigned int nt = p->total, ext = p->nnew > p->nold ? p->nnew - p->nold : 0, i, j, k, c, t;
	unsigned int ca = 0, fa = 0, ea = 0, *oc, *ou, *ow, *otw, *oflag, *nc, *sel, *buf;
	double *cvals = NULL, *fvals = NULL;
	if (p->failed)
		return 0;
	for (k = 0; k < m->na; k++)
	{
		if (m->attr[k].element == QAWS_MESH_CORNER)
			ca += m->attr[k].n;
		if (m->attr[k].element == QAWS_MESH_FACE)
			fa += m->attr[k].n;
		if (m->attr[k].element == QAWS_MESH_EDGE)
			ea += m->attr[k].n;
	}
	if (!hm_reserve(m, 0, p->ncv, ext))
		return 0;
	buf = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), (5 * (size_t)nt + 2 * (size_t)p->ncv + 1)));
	if (ca)
		cvals = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), ca, (nt ? nt : 1)));
	if (fa)
		fvals = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), fa, (p->nold ? p->nold : 1)));
	if (!buf || (ca && !cvals) || (fa && !fvals))
	{
		qaws_internal_free(buf); qaws_internal_free(cvals); qaws_internal_free(fvals);
		return 0;
	}
	oc = buf; ou = oc + nt; ow = ou + nt; otw = ow + nt; oflag = otw + nt; nc = oflag + nt; sel = nc + p->ncv;
	for (i = 0; i < p->nold; i++)
		m->fflag[p->old[i]] |= HM_MARK;
	/* the old corners: edges (read with the renaming), outside neighbours, flags, attributes */
	for (i = 0; i < p->nold; i++)
	{
		unsigned int f = p->old[i];
		for (j = 0; j < m->fn[f]; j++)
		{
			unsigned int pos = p->off[i] + j, o, a = 0;
			c = m->fc[f] + j;
			o = m->co[c];
			oc[pos] = c;
			ou[pos] = m->cv[c] == from ? to : m->cv[c];
			ow[pos] = hm_to(m, c) == from ? to : hm_to(m, c);
			otw[pos] = o != HM_NONE && !(m->fflag[m->cf[o]] & HM_MARK) ? o : HM_NONE;
			oflag[pos] = m->cflag[c] & HM_FEATURE;
			for (k = 0; k < m->na; k++)
				if (m->attr[k].element == QAWS_MESH_CORNER)
				{
					memcpy(cvals + (size_t)ca * pos + a, hm_attr_at(m, k, c), sizeof(double) * m->attr[k].n);
					a += m->attr[k].n;
				}
		}
		for (k = 0, j = 0; k < m->na; k++)
			if (m->attr[k].element == QAWS_MESH_FACE)
			{
				memcpy(fvals + (size_t)fa * i + j, hm_attr_at(m, k, f), sizeof(double) * m->attr[k].n);
				j += m->attr[k].n;
			}
	}
	for (i = 0; i < p->nold; i++)
	{
		unsigned int f = p->old[i];
		m->fflag[f] &= (unsigned char)(0xFFu & ~HM_MARK);
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			m->cf[c] = HM_NONE;   /* (freed; its vertex stays readable) */
		m->cfree += m->fn[f];
	}
	/* the new faces: the old faces' indices first, then appended */
	for (i = 0, t = 0; i < p->nnew; i++)
	{
		unsigned int f = i < p->nold ? p->old[i] : m->nf++, a;
		m->fc[f] = m->nc;
		m->fn[f] = p->fs[i];
		m->fflag[f] = 0;
		for (j = 0; j < p->fs[i]; j++, t++)
		{
			c = m->nc++;
			m->cv[c] = p->fv[t];
			m->co[c] = HM_NONE;
			m->cf[c] = f;
			m->cflag[c] = p->hint[t];
			nc[t] = c;
		}
		for (k = 0, a = 0; k < m->na; k++)
			if (m->attr[k].element == QAWS_MESH_FACE)
			{
				if (p->fa[i] == HM_NONE)
					memset(hm_attr_at(m, k, f), 0, sizeof(double) * m->attr[k].n);
				else
					memcpy(hm_attr_at(m, k, f), fvals + (size_t)fa * p->fa[i] + a, sizeof(double) * m->attr[k].n);
				a += m->attr[k].n;
			}
	}
	for (i = p->nnew; i < p->nold; i++)
	{
		m->fflag[p->old[i]] |= HM_DELETED;
		m->hole[1][m->nhole[1]++] = p->old[i];
	}
	/* stitching by edges, flags by edge */
	for (t = 0; t < p->ncv; t++)
	{
		unsigned int u, w, s, fl = 0;
		c = nc[t];
		u = m->cv[c];
		w = hm_to(m, c);
		sel[t] = p->esrc[t];
		for (j = 0; j < nt; j++)
			if ((ou[j] == u && ow[j] == w) || (ou[j] == w && ow[j] == u))
			{
				fl |= oflag[j];
				if (sel[t] == HM_NONE)
					sel[t] = j;
			}
		m->cflag[c] |= (unsigned char)fl;
		if (m->co[c] != HM_NONE)
			continue;
		for (s = t + 1; s < p->ncv; s++)
		{
			unsigned int d = nc[s];
			if (m->co[d] == HM_NONE && m->cv[d] == w && hm_to(m, d) == u)
			{
				m->co[c] = d;
				m->co[d] = c;
				break;
			}
		}
		if (m->co[c] == HM_NONE)
			for (j = 0; j < nt; j++)
				if (otw[j] != HM_NONE && ou[j] == u && ow[j] == w)
				{
					m->co[c] = otw[j];
					m->co[otw[j]] = c;
					break;
				}
	}
	for (j = 0; j < nt; j++)
		if (otw[j] != HM_NONE && m->co[otw[j]] != HM_NONE && m->cf[m->co[otw[j]]] == HM_NONE)
			m->co[otw[j]] = HM_NONE;   /* (its edge has nothing on the patch's side now) */
	for (t = 0; t < p->ncv; t++)
	{
		unsigned int o = m->co[nc[t]];
		unsigned char fl;
		if (o == HM_NONE)
			continue;
		fl = (unsigned char)((m->cflag[nc[t]] | m->cflag[o]) & HM_FEATURE);
		m->cflag[nc[t]] = (unsigned char)((m->cflag[nc[t]] & ~HM_FEATURE) | fl);
		m->cflag[o] = (unsigned char)((m->cflag[o] & ~HM_FEATURE) | fl);
	}
	/* edge attributes: an edge kept, or a piece of one (a split's halves), keeps its value; a new edge stitched to an
	   outside corner takes that one's, its twin among the new ones the twin's; the others zeros */
	for (t = 0; ea && t < p->ncv; t++)
	{
		unsigned int o = m->co[nc[t]], s = sel[t] != HM_NONE ? oc[sel[t]] : o != HM_NONE && o < nc[0] ? o : HM_NONE;
		for (k = 0; k < m->na; k++)
			if (m->attr[k].element == QAWS_MESH_EDGE)
			{
				if (s == HM_NONE)
					memset(hm_attr_at(m, k, nc[t]), 0, sizeof(double) * m->attr[k].n);
				else
					memcpy(hm_attr_at(m, k, nc[t]), hm_attr_at(m, k, s), sizeof(double) * m->attr[k].n);
			}
	}
	for (t = 0; ea && t < p->ncv; t++)
	{
		unsigned int o = m->co[nc[t]];
		if (sel[t] == HM_NONE && o != HM_NONE && o >= nc[0])
			for (k = 0; k < m->na; k++)
				if (m->attr[k].element == QAWS_MESH_EDGE)
					memcpy(hm_attr_at(m, k, nc[t]), hm_attr_at(m, k, o), sizeof(double) * m->attr[k].n);
	}
	/* corner attributes from their sources */
	for (t = 0; ca && t < p->ncv; t++)
	{
		unsigned int a = 0;
		for (k = 0; k < m->na; k++)
			if (m->attr[k].element == QAWS_MESH_CORNER)
			{
				hm_blend(hm_attr_at(m, k, nc[t]), cvals + a, ca, p->sp + p->s0[t], p->sw + p->s0[t], p->sn[t], m->attr[k].n, m->attr[k].rule);
				a += m->attr[k].n;
			}
	}
	/* every touched vertex's corner */
	for (t = 0; t < p->ncv + 3 * nt; t++)
	{
		unsigned int v = HM_NONE;
		if (t < p->ncv)
			v = m->cv[nc[t]];
		else
		{
			j = (t - p->ncv) / 3;
			k = (t - p->ncv) % 3;
			if (k == 0)
				v = m->cv[oc[j]];
			else if (otw[j] != HM_NONE)
				v = k == 1 ? m->cv[otw[j]] : hm_to(m, otw[j]);
		}
		if (v == HM_NONE || (m->vflag[v] & HM_DELETED))
			continue;
		hm_fix_vertex(m, v, hm_any_corner(m, v, nc, p->ncv, otw, nt));
	}
	qaws_internal_free(buf); qaws_internal_free(cvals); qaws_internal_free(fvals);
	return 1;
}

unsigned int hm_split_edge(qaws_hm* m, unsigned int c, double t, double const* x, int triangulate)
{
	unsigned int d = m->co[c], f = m->cf[c], g = d == HM_NONE ? HM_NONE : m->cf[d], a = m->cv[c], b = hm_to(m, c);
	unsigned int v = m->nv, kf, kg = 0, k, p = 0, q = 0;
	unsigned char feat = (unsigned char)(m->cflag[c] & HM_FEATURE);
	int ftri = triangulate && m->fn[f] == 3, gtri = triangulate && g != HM_NONE && m->fn[g] == 3;
	if (!hm_reserve(m, 1, 0, 0))
		return HM_NONE;
	pb_reset(m);
	kf = pb_old(m, f);
	if (g != HM_NONE)
		kg = pb_old(m, g);
	/* the faces' first parts keep their indices, the second halves follow */
	pb_face(m, kf);
	if (ftri)
	{
		p = m->cv[hm_prev(m, c)];
		pb_copy(m, a, kf, feat);
		pb_edge(m, pb_cpos(m, kf, c));
		pb_corner(m, v, 0); pb_src(m, pb_pos(m, kf, a), 1 - t); pb_src(m, pb_pos(m, kf, b), t);
		pb_copy(m, p, kf, 0);
	}
	else
		for (k = c;; k = hm_next(m, k))
		{
			pb_copy(m, m->cv[k], kf, k == c ? feat : 0);
			if (k == c)
			{
				pb_edge(m, pb_cpos(m, kf, c));
				pb_corner(m, v, feat); pb_src(m, pb_pos(m, kf, a), 1 - t); pb_src(m, pb_pos(m, kf, b), t);
				pb_edge(m, pb_cpos(m, kf, c));
			}
			if (hm_next(m, k) == c)
				break;
		}
	if (g != HM_NONE)
	{
		pb_face(m, kg);
		if (gtri)
		{
			q = m->cv[hm_prev(m, d)];
			pb_copy(m, b, kg, feat);
			pb_edge(m, pb_cpos(m, kg, d));
			pb_corner(m, v, 0); pb_src(m, pb_pos(m, kg, a), 1 - t); pb_src(m, pb_pos(m, kg, b), t);
			pb_copy(m, q, kg, 0);
		}
		else
			for (k = d;; k = hm_next(m, k))
			{
				pb_copy(m, m->cv[k], kg, k == d ? feat : 0);
				if (k == d)
				{
					pb_edge(m, pb_cpos(m, kg, d));
					pb_corner(m, v, feat); pb_src(m, pb_pos(m, kg, a), 1 - t); pb_src(m, pb_pos(m, kg, b), t);
					pb_edge(m, pb_cpos(m, kg, d));
				}
				if (hm_next(m, k) == d)
					break;
			}
	}
	if (ftri)
	{
		pb_face(m, kf);
		pb_corner(m, v, feat); pb_src(m, pb_pos(m, kf, a), 1 - t); pb_src(m, pb_pos(m, kf, b), t);
		pb_edge(m, pb_cpos(m, kf, c));
		pb_copy(m, b, kf, 0);
		pb_copy(m, p, kf, 0);
	}
	if (gtri)
	{
		pb_face(m, kg);
		pb_corner(m, v, feat); pb_src(m, pb_pos(m, kg, a), 1 - t); pb_src(m, pb_pos(m, kg, b), t);
		pb_edge(m, pb_cpos(m, kg, d));
		pb_copy(m, a, kg, 0);
		pb_copy(m, q, kg, 0);
	}
	if (m->pb.failed)
		return HM_NONE;
	hm_new_vertex(m, x);
	for (k = 0; k < m->na; k++)
		if (m->attr[k].element == QAWS_MESH_VERTEX)
			hm_lerp(hm_attr_at(m, k, v), hm_attr_at(m, k, a), hm_attr_at(m, k, b), t, m->attr[k].n, m->attr[k].rule);
	if (!hm_replace(m, HM_NONE, HM_NONE))
	{
		m->nv--;
		return HM_NONE;
	}
	return v;
}

unsigned int hm_split_face(qaws_hm* m, unsigned int f, double const* x)
{
	unsigned int s = m->fn[f], i, j, k, v = m->nv, kf, first = m->fc[f], vx3[3], *vx;
	double w3[3], *w, cx[16];
	int centroid = x == NULL;
	if (!hm_reserve(m, 1, 0, 0))
		return HM_NONE;
	w = s == 3 ? w3 : (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), s));
	vx = s == 3 ? vx3 : (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), s));
	if (!w || !vx)
	{
		if (w != w3)
			qaws_internal_free(w);
		if (vx != vx3)
			qaws_internal_free(vx);
		return HM_NONE;
	}
	for (i = 0; i < s; i++)
		w[i] = 1.0 / s;
	if (centroid)
	{
		/* the centroid: equal weights, exactly */
		for (k = 0; k < m->dim; k++)
			for (cx[k] = 0, i = 0; i < s; i++)
				cx[k] += m->x[(size_t)m->dim * m->cv[first + i] + k] / s;
		x = cx;
	}
	else if (s == 3)
	{
		/* barycentric in the triangle (equal weights when it is too thin to tell) */
		double const *p0 = m->x + (size_t)m->dim * m->cv[first], *p1 = m->x + (size_t)m->dim * m->cv[first + 1],
			*p2 = m->x + (size_t)m->dim * m->cv[first + 2];
		double d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0, den;
		for (k = 0; k < m->dim; k++)
		{
			double u0 = p1[k] - p0[k], u1 = p2[k] - p0[k], q = x[k] - p0[k];
			d00 += u0 * u0; d01 += u0 * u1; d11 += u1 * u1; d20 += q * u0; d21 += q * u1;
		}
		den = d00 * d11 - d01 * d01;
		if (den > 1e-12 * d00 * d11)
		{
			w[1] = (d11 * d20 - d01 * d21) / den;
			w[2] = (d00 * d21 - d01 * d20) / den;
			w[0] = 1 - w[1] - w[2];
		}
	}
	pb_reset(m);
	kf = pb_old(m, f);
	for (i = 0; i < s; i++)
	{
		pb_face(m, kf);
		pb_corner(m, m->cv[first + i], 0); pb_src(m, m->pb.off[kf] + i, 1);
		pb_corner(m, m->cv[first + (i + 1) % s], 0); pb_src(m, m->pb.off[kf] + (i + 1) % s, 1);
		pb_corner(m, v, 0);
		for (j = 0; j < s; j++)
			pb_src(m, m->pb.off[kf] + j, w[j]);
	}
	if (!m->pb.failed)
	{
		hm_new_vertex(m, x);
		for (j = 0; j < s; j++)
			vx[j] = m->cv[first + j];
		for (k = 0; k < m->na; k++)
			if (m->attr[k].element == QAWS_MESH_VERTEX)
				hm_blend(hm_attr_at(m, k, v), m->attr[k].v, m->attr[k].n, vx, w, s, m->attr[k].n, m->attr[k].rule);
	}
	if (w != w3)
		qaws_internal_free(w);
	if (vx != vx3)
		qaws_internal_free(vx);
	if (m->pb.failed)
		return HM_NONE;
	if (!hm_replace(m, HM_NONE, HM_NONE))
	{
		m->nv--;
		return HM_NONE;
	}
	return v;
}

int hm_cut_face(qaws_hm* m, unsigned int c0, unsigned int c1)
{
	unsigned int f = m->cf[c0], k, kf;
	pb_reset(m);
	kf = pb_old(m, f);
	pb_face(m, kf);
	for (k = c0;; k = hm_next(m, k))
	{
		pb_corner(m, m->cv[k], 0);
		pb_src(m, pb_cpos(m, kf, k), 1);
		if (k == c1)
			break;
	}
	pb_face(m, kf);
	for (k = c1;; k = hm_next(m, k))
	{
		pb_corner(m, m->cv[k], 0);
		pb_src(m, pb_cpos(m, kf, k), 1);
		if (k == c0)
			break;
	}
	return hm_replace(m, HM_NONE, HM_NONE);
}

int hm_can_collapse(qaws_hm const* m, unsigned int c, unsigned int a)
{
	unsigned int d = m->co[c], f = m->cf[c], g = d == HM_NONE ? HM_NONE : m->cf[d], b = a == m->cv[c] ? hm_to(m, c) : m->cv[c];
	unsigned int vl = HM_NONE, vr = HM_NONE, k, last = HM_NONE;
	if (m->fn[f] == 3)
	{
		vl = m->cv[hm_prev(m, c)];
		if (m->co[hm_next(m, c)] == HM_NONE && m->co[hm_prev(m, c)] == HM_NONE)
			return 0;
	}
	if (g != HM_NONE && m->fn[g] == 3)
	{
		vr = m->cv[hm_prev(m, d)];
		if (m->co[hm_next(m, d)] == HM_NONE && m->co[hm_prev(m, d)] == HM_NONE)
			return 0;
	}
	if (vl != HM_NONE && vl == vr)
		return 0;
	if (d != HM_NONE && hm_is_border_vertex(m, a) && hm_is_border_vertex(m, b))
		return 0;
	/* a triangle's third corner left with two edges would fold (a tetrahedron's edge) */
	if ((vl != HM_NONE && !hm_is_border_vertex(m, vl) && hm_valence(m, vl) <= 3)
		|| (vr != HM_NONE && !hm_is_border_vertex(m, vr) && hm_valence(m, vr) <= 3))
		return 0;
	/* the one-rings meet only at the edge's triangles' third corners, and no
	   other face round a holds b */
	for (k = m->vc[a]; k != HM_NONE; k = hm_round(m, a, k))
	{
		unsigned int w = hm_to(m, k), e, h = m->cf[k];
		if (w != b && w != vl && w != vr && hm_find_edge(m, w, b) != HM_NONE)
			return 0;
		if (h != f && h != g)
			for (e = m->fc[h]; e < m->fc[h] + m->fn[h]; e++)
				if (m->cv[e] == b)
					return 0;
		last = k;
	}
	if (last != HM_NONE && hm_is_border_vertex(m, a))
	{
		unsigned int w = m->cv[hm_prev(m, last)];
		if (w != b && w != vl && w != vr && hm_find_edge(m, w, b) != HM_NONE)
			return 0;
	}
	return 1;
}

int hm_collapse(qaws_hm* m, unsigned int c, unsigned int a, double t, double const* x)
{
	unsigned int d = m->co[c], f = m->cf[c], g = d == HM_NONE ? HM_NONE : m->cf[d], b = a == m->cv[c] ? hm_to(m, c) : m->cv[c], k, e, i;
	pb_reset(m);
	/* a's faces: those off the edge keep their indices, polygons on the edge shrink, its triangles go */
	for (k = m->vc[a]; k != HM_NONE; k = hm_round(m, a, k))
		if (m->cf[k] != f && m->cf[k] != g)
			pb_old(m, m->cf[k]);
	if (m->fn[f] > 3)
		pb_old(m, f);
	if (g != HM_NONE && m->fn[g] > 3)
		pb_old(m, g);
	if (m->fn[f] == 3)
		pb_old(m, f);
	if (g != HM_NONE && m->fn[g] == 3)
		pb_old(m, g);
	for (i = 0; i < m->pb.nold && !m->pb.failed; i++)
	{
		unsigned int h = m->pb.old[i], edge = h == f || h == g;
		if (edge && m->fn[h] == 3)
			break;
		pb_face(m, i);
		for (e = m->fc[h]; e < m->fc[h] + m->fn[h]; e++)
		{
			if (edge && m->cv[e] == a)
				continue;   /* (the polygon loses a, b stays) */
			pb_corner(m, m->cv[e] == a ? b : m->cv[e], 0);
			pb_src(m, pb_cpos(m, i, e), 1);
		}
	}
	if (m->pb.failed)
		return 0;
	if (!hm_replace(m, a, b))
		return 0;
	for (k = 0; k < m->na; k++)
		if (m->attr[k].element == QAWS_MESH_VERTEX)
			hm_lerp(hm_attr_at(m, k, b), hm_attr_at(m, k, a), hm_attr_at(m, k, b), t, m->attr[k].n, m->attr[k].rule);
	for (k = 0; k < m->dim; k++)
		m->x[(size_t)m->dim * b + k] = x[k];
	m->vflag[b] |= m->vflag[a] & HM_CORNER;
	hm_kill_vertex(m, a);
	return 1;
}

int hm_can_flip(qaws_hm const* m, unsigned int c)
{
	unsigned int d = m->co[c], p, q, a, b;
	if (d == HM_NONE || m->fn[m->cf[c]] != 3 || m->fn[m->cf[d]] != 3)
		return 0;
	a = m->cv[c];
	b = hm_to(m, c);
	p = m->cv[hm_prev(m, c)];
	q = m->cv[hm_prev(m, d)];
	if (p == q || hm_find_edge(m, p, q) != HM_NONE)
		return 0;
	/* an interior end left with two edges would fold its faces together */
	if ((!hm_is_border_vertex(m, a) && hm_valence(m, a) <= 3) || (!hm_is_border_vertex(m, b) && hm_valence(m, b) <= 3))
		return 0;
	return 1;
}

int hm_flip(qaws_hm* m, unsigned int c)
{
	unsigned int d = m->co[c], f = m->cf[c], g = m->cf[d], a = m->cv[c], b = hm_to(m, c);
	unsigned int p = m->cv[hm_prev(m, c)], q = m->cv[hm_prev(m, d)], kf, kg;
	pb_reset(m);
	kf = pb_old(m, f);
	kg = pb_old(m, g);
	/* (a, b, p) and (b, a, q) become (a, q, p) and (q, b, p) */
	pb_face(m, kf); pb_copy(m, a, kf, 0); pb_copy(m, q, kg, 0); pb_copy(m, p, kf, 0);
	pb_face(m, kg); pb_copy(m, q, kg, 0); pb_copy(m, b, kg, 0); pb_copy(m, p, kf, 0);
	return hm_replace(m, HM_NONE, HM_NONE);
}

/* ======================================================================== */
/*  Compaction                                                              */
/* ======================================================================== */

static void hm_move_vertex(qaws_hm* m, unsigned int j, unsigned int i)
{
	unsigned int c, k;
	for (c = m->vc[j]; c != HM_NONE; c = hm_round(m, j, c))
		m->cv[c] = i;
	memcpy(m->x + (size_t)m->dim * i, m->x + (size_t)m->dim * j, sizeof(double) * m->dim);
	m->vc[i] = m->vc[j];
	m->vsrc[i] = m->vsrc[j];
	m->vflag[i] = m->vflag[j];
	if (m->foot)
		m->foot[i] = m->foot[j];
	for (k = 0; k < m->na; k++)
		if (m->attr[k].element == QAWS_MESH_VERTEX)
			memcpy(hm_attr_at(m, k, i), hm_attr_at(m, k, j), sizeof(double) * m->attr[k].n);
}

static void hm_move_face(qaws_hm* m, unsigned int j, unsigned int i)
{
	unsigned int c, k;
	for (c = m->fc[j]; c < m->fc[j] + m->fn[j]; c++)
		m->cf[c] = i;
	m->fc[i] = m->fc[j];
	m->fn[i] = m->fn[j];
	m->fflag[i] = m->fflag[j];
	for (k = 0; k < m->na; k++)
		if (m->attr[k].element == QAWS_MESH_FACE)
			memcpy(hm_attr_at(m, k, i), hm_attr_at(m, k, j), sizeof(double) * m->attr[k].n);
}

/* the corner pool without its free slots (runs keep their order) */
static void hm_pool(qaws_hm* m)
{
	unsigned int* map = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nc + 1)));
	unsigned int c, n = 0, k, v;
	if (!map)
		return;   /* (the garbage stays) */
	for (c = 0; c < m->nc; c++)
		map[c] = m->cf[c] == HM_NONE ? HM_NONE : n++;
	for (c = 0; c < m->nc; c++)
	{
		unsigned int d = map[c];
		if (d == HM_NONE)
			continue;
		if (m->fc[m->cf[c]] == c)
			m->fc[m->cf[c]] = d;
		m->cv[d] = m->cv[c];
		m->co[d] = m->co[c];
		m->cf[d] = m->cf[c];
		m->cflag[d] = m->cflag[c];
		for (k = 0; k < m->na; k++)
			if (HM_ON_CORNERS(m->attr[k].element) && d != c)
				memcpy(hm_attr_at(m, k, d), hm_attr_at(m, k, c), sizeof(double) * m->attr[k].n);
	}
	for (c = 0; c < n; c++)
		if (m->co[c] != HM_NONE)
			m->co[c] = map[m->co[c]];
	for (v = 0; v < m->nv; v++)
		if (m->vc[v] != HM_NONE)
			m->vc[v] = map[m->vc[v]];
	m->nc = n;
	m->cfree = 0;
	qaws_internal_free(map);
}

void hm_compact(qaws_hm* m)
{
	unsigned int kind, k;
	for (kind = 0; kind < 2; kind++)
	{
		unsigned int* n = kind == 0 ? &m->nv : &m->nf;
		unsigned char* fl = kind == 0 ? m->vflag : m->fflag;
		for (k = 0; k < m->nhole[kind]; k++)
		{
			unsigned int i = m->hole[kind][k];
			while (*n > 0 && (fl[*n - 1] & HM_DELETED))
				(*n)--;
			if (i >= *n)
				continue;
			if (kind == 0)
				hm_move_vertex(m, *n - 1, i);
			else
				hm_move_face(m, *n - 1, i);
			(*n)--;
		}
		while (*n > 0 && (fl[*n - 1] & HM_DELETED))
			(*n)--;
		m->nhole[kind] = 0;
	}
	if (m->cfree && 2 * m->cfree >= m->nc)
		hm_pool(m);
}

/* ======================================================================== */
/*  Attributes and footpoints                                               */
/* ======================================================================== */

unsigned int hm_attr_add(qaws_hm* m, char const* name, unsigned int element, unsigned int rule, unsigned int n)
{
	hm_attr* g = (hm_attr*)qaws_internal_realloc(m->attr, qaws_internal_mul2(sizeof(hm_attr), (m->na + 1)));
	hm_attr* a;
	size_t slots;
	if (!g)
		return HM_NONE;
	m->attr = g;
	a = &g[m->na];
	memset(a, 0, sizeof *a);
	slots = element == QAWS_MESH_VERTEX ? m->capv : HM_ON_CORNERS(element) ? m->capc : m->capf;
	a->v = (double*)qaws_internal_zalloc(slots ? qaws_internal_mul2(slots, n) : 1, sizeof(double));
	if (!a->v)
		return HM_NONE;
	if (name)
	{
		strncpy(a->name, name, sizeof a->name - 1);
		a->name[sizeof a->name - 1] = 0;
	}
	a->element = element;
	a->rule = rule;
	a->n = n;
	return m->na++;
}

void hm_attr_remove(qaws_hm* m, unsigned int a)
{
	qaws_internal_free(m->attr[a].v);
	memmove(m->attr + a, m->attr + a + 1, sizeof(hm_attr) * (m->na - a - 1));
	m->na--;
}

int hm_foot_enable(qaws_hm* m)
{
	unsigned int v;
	if (m->foot)
		return 1;
	m->foot = (qaws_mesh_location*)qaws_internal_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (m->capv ? m->capv : 1)));
	if (!m->foot)
		return 0;
	for (v = 0; v < m->nv; v++)
	{
		m->foot[v].patch = HM_NONE;
		m->foot[v].a = m->foot[v].b = 0;
	}
	return 1;
}

unsigned int hm_feature_degree(qaws_hm const* m, unsigned int v)
{
	unsigned int c, last = HM_NONE, n = 0;
	for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
	{
		n += (m->cflag[c] & HM_FEATURE) != 0;
		last = c;
	}
	if (last != HM_NONE && hm_is_border_vertex(m, v))
		n += (m->cflag[hm_prev(m, last)] & HM_FEATURE) != 0;
	return n;
}

/* ======================================================================== */
/*  Quad operations (Tarini et al. 2010)                                    */
/* ======================================================================== */

/* the hexagon round the two quads at c's edge (a, b): b, f2, f3, a, g2, g3 counter-clockwise (f = a b f2 f3, g = b a g2 g3);
   0 when the faces are not two quads sharing only that edge */
static int hm_hexagon(qaws_hm const* m, unsigned int c, unsigned int* h)
{
	unsigned int d = m->co[c], i, j;
	if (d == HM_NONE || m->cf[c] == m->cf[d] || m->fn[m->cf[c]] != 4 || m->fn[m->cf[d]] != 4)
		return 0;
	h[0] = hm_to(m, c);
	h[1] = m->cv[hm_next(m, hm_next(m, c))];
	h[2] = m->cv[hm_prev(m, c)];
	h[3] = m->cv[c];
	h[4] = m->cv[hm_next(m, hm_next(m, d))];
	h[5] = m->cv[hm_prev(m, d)];
	for (i = 0; i < 6; i++)
		for (j = i + 1; j < 6; j++)
			if (h[i] == h[j])
				return 0;
	return 1;
}

/* an end keeping three edges (two on the border) */
static int hm_keeps(qaws_hm const* m, unsigned int v)
{
	return hm_valence(m, v) >= 4u - (unsigned int)hm_is_border_vertex(m, v);
}

int hm_can_rotate(qaws_hm const* m, unsigned int c, int dir)
{
	unsigned int h[6], p, q;
	if (c >= m->nc || m->cf[c] == HM_NONE || (m->cflag[c] & HM_FEATURE) || !hm_hexagon(m, c, h))
		return 0;
	p = dir > 0 ? h[1] : h[5];
	q = dir > 0 ? h[4] : h[2];
	if (((m->vflag[h[0]] | m->vflag[h[3]] | m->vflag[p] | m->vflag[q]) & HM_CORNER) || hm_find_edge(m, p, q) != HM_NONE)
		return 0;
	return hm_keeps(m, h[0]) && hm_keeps(m, h[3]);
}

int hm_rotate(qaws_hm* m, unsigned int c, int dir)
{
	unsigned int h[6], kf, kg, s = dir > 0 ? 1 : 5, i, j;
	if (!hm_can_rotate(m, c, dir))
		return 0;
	hm_hexagon(m, c, h);
	pb_reset(m);
	kf = pb_old(m, m->cf[c]);
	kg = pb_old(m, m->cf[m->co[c]]);
	/* the hexagon cut from h[s] to h[s + 3]; each corner from the old face it was in (f: h0..h3) */
	for (i = 0; i < 2; i++)
	{
		pb_face(m, i ? kg : kf);
		for (j = 0; j < 4; j++)
		{
			unsigned int at = (s + 3 * i + j) % 6;
			pb_copy(m, h[at], at <= 3 ? kf : kg, 0);
		}
	}
	return hm_replace(m, HM_NONE, HM_NONE);
}

int hm_can_collapse_diagonal(qaws_hm const* m, unsigned int f, unsigned int k)
{
	unsigned int c0, A, B, C, D, c, e;
	if (f >= m->nf || (m->fflag[f] & HM_DELETED) || m->fn[f] != 4 || k > 3)
		return 0;
	c0 = m->fc[f];
	A = m->cv[c0 + k];
	B = m->cv[c0 + (k + 1) % 4];
	C = m->cv[c0 + (k + 2) % 4];
	D = m->cv[c0 + (k + 3) % 4];
	if (A == C || B == D || ((m->vflag[A] | m->vflag[C]) & HM_CORNER))
		return 0;
	/* C goes: inside, off the features (A keeps the lines it is on) */
	if (hm_is_border_vertex(m, C) || hm_feature_degree(m, C) || !hm_keeps(m, B) || !hm_keeps(m, D) || hm_valence(m, A) + hm_valence(m, C) < 5)
		return 0;
	/* C's neighbours but B and D not A's, and no face but f holds both */
	for (c = m->vc[C]; c != HM_NONE; c = hm_round(m, C, c))
	{
		unsigned int w = hm_to(m, c), h = m->cf[c];
		if (w != B && w != D && (w == A || hm_find_edge(m, w, A) != HM_NONE))
			return 0;
		if (h != f)
			for (e = m->fc[h]; e < m->fc[h] + m->fn[h]; e++)
				if (m->cv[e] == A)
					return 0;
	}
	return 1;
}

int hm_collapse_diagonal(qaws_hm* m, unsigned int f, unsigned int k, double const* x)
{
	unsigned int A, C, c, e, i;
	double t = 0, l = 0;
	if (!hm_can_collapse_diagonal(m, f, k))
		return 0;
	A = m->cv[m->fc[f] + k];
	C = m->cv[m->fc[f] + (k + 2) % 4];
	pb_reset(m);
	/* C's faces but f keep their indices, A for C; f goes */
	for (c = m->vc[C]; c != HM_NONE; c = hm_round(m, C, c))
		if (m->cf[c] != f)
			pb_old(m, m->cf[c]);
	pb_old(m, f);
	for (i = 0; i + 1 < m->pb.nold && !m->pb.failed; i++)
	{
		unsigned int h = m->pb.old[i];
		pb_face(m, i);
		for (e = m->fc[h]; e < m->fc[h] + m->fn[h]; e++)
		{
			pb_corner(m, m->cv[e] == C ? A : m->cv[e], 0);
			pb_src(m, pb_cpos(m, i, e), 1);
		}
	}
	if (m->pb.failed || !hm_replace(m, C, A))
		return 0;
	/* A's vertex attributes from A's and C's, at x's place along A to C */
	for (i = 0; i < m->dim; i++)
	{
		double u = m->x[(size_t)m->dim * C + i] - m->x[(size_t)m->dim * A + i];
		t += (x[i] - m->x[(size_t)m->dim * A + i]) * u;
		l += u * u;
	}
	t = l > 0 ? t / l : 0;
	t = t < 0 ? 0 : t > 1 ? 1 : t;
	for (i = 0; i < m->na; i++)
		if (m->attr[i].element == QAWS_MESH_VERTEX)
			hm_lerp(hm_attr_at(m, i, A), hm_attr_at(m, i, A), hm_attr_at(m, i, C), t, m->attr[i].n, m->attr[i].rule);
	for (i = 0; i < m->dim; i++)
		m->x[(size_t)m->dim * A + i] = x[i];
	hm_kill_vertex(m, C);
	return 1;
}

int hm_remove_doublet(qaws_hm* m, unsigned int v)
{
	unsigned int c, d, a, b, x, y, kf, kg;
	if (v >= m->nv || m->vc[v] == HM_NONE || (m->vflag[v] & (HM_DELETED | HM_CORNER)) || hm_is_border_vertex(m, v) || hm_valence(m, v) != 2
		|| hm_feature_degree(m, v))
		return 0;
	/* f = v a x b and g = v b y a become a x b y */
	c = m->vc[v];
	d = hm_ccw(m, c);
	if (d == HM_NONE || m->cf[c] == m->cf[d] || m->fn[m->cf[c]] != 4 || m->fn[m->cf[d]] != 4)
		return 0;
	a = hm_to(m, c);
	x = m->cv[hm_next(m, hm_next(m, c))];
	b = m->cv[hm_prev(m, c)];
	y = m->cv[hm_next(m, hm_next(m, d))];
	if (hm_to(m, d) != b || m->cv[hm_prev(m, d)] != a || a == b || x == y || hm_valence(m, a) < 3 || hm_valence(m, b) < 3)
		return 0;
	pb_reset(m);
	kf = pb_old(m, m->cf[c]);
	kg = pb_old(m, m->cf[d]);
	pb_face(m, kf);
	pb_copy(m, a, kf, 0);
	pb_copy(m, x, kf, 0);
	pb_copy(m, b, kf, 0);
	pb_copy(m, y, kg, 0);
	if (!hm_replace(m, HM_NONE, HM_NONE))
		return 0;
	hm_kill_vertex(m, v);
	return 1;
}
