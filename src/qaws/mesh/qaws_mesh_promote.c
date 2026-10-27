/*
 * Attribute promotion (qaws_mesh_promote.h): the sources meeting each target
 * element from the mesh's faces and edges, combined by the method.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_promote.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int pm_cmp_u64(void const* a, void const* b)
{
	unsigned long long x = *(unsigned long long const*)a, y = *(unsigned long long const*)b;
	return x < y ? -1 : x > y;
}

#define PM_REC 18u   /* a mode record: its component count, up to 16 values, its source */

static int pm_cmp_double(void const* a, void const* b)
{
	double x = *(double const*)a, y = *(double const*)b;
	return x < y ? -1 : x > y;
}

/* records by value, then by source (so a run starts at its first source) */
static int pm_cmp_rec(void const* a, void const* b)
{
	double const *x = (double const*)a, *y = (double const*)b;
	unsigned int k, n = (unsigned int)x[0];
	for (k = 1; k <= n + 1; k++)
		if (x[k] != y[k])
			return x[k] < y[k] ? -1 : 1;
	return 0;
}

/* the method over the sources src (ns of them, values vals of n components each) into out, the giving source into *idx;
   w and tmp: ns scratch each */
static void pm_combine(double* out, int* idx, double const* vals, unsigned int const* src, unsigned int ns, unsigned int n, unsigned int rule,
	qaws_promote_method method, double* w, double* tmp)
{
	unsigned int i, k;
	*idx = -1;
	for (k = 0; k < n; k++)
		out[k] = 0;
	if (!ns)
		return;
	switch (method)
	{
	case QAWS_PROMOTE_MAXIMUM:
	case QAWS_PROMOTE_MINIMUM:
		for (k = 0; k < n; k++)
		{
			unsigned int b = 0;
			for (i = 1; i < ns; i++)
			{
				double v = vals[(size_t)n * src[i] + k], bv = vals[(size_t)n * src[b] + k];
				if (method == QAWS_PROMOTE_MAXIMUM ? v > bv : v < bv)
					b = i;
			}
			out[k] = vals[(size_t)n * src[b] + k];
			if (k == 0)
				*idx = (int)src[b];
		}
		break;
	case QAWS_PROMOTE_AVERAGE:
		for (i = 0; i < ns; i++)
			w[i] = 1.0 / ns;
		hm_blend(out, vals, n, src, w, ns, n, rule);
		break;
	case QAWS_PROMOTE_MODE:
	{
		/* (records [n, values, source] sorted by value, then the longest run: ties the least value, sorted first; the first
		   source holding it) */
		double* r = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), PM_REC, (size_t)ns));
		unsigned int best = 0, bc = 0, run = 0;
		if (!r)
			break;
		for (i = 0; i < ns; i++)
		{
			r[(size_t)PM_REC * i] = n;
			for (k = 0; k < n; k++)
				r[(size_t)PM_REC * i + 1 + k] = vals[(size_t)n * src[i] + k];
			r[(size_t)PM_REC * i + 1 + n] = src[i];
		}
		qsort(r, ns, sizeof(double) * PM_REC, pm_cmp_rec);
		for (i = 0; i < ns; i++)
		{
			int same = i > 0;
			for (k = 0; k < n && same; k++)
				same = r[(size_t)PM_REC * i + 1 + k] == r[(size_t)PM_REC * (i - 1) + 1 + k];
			run = same ? run + 1 : 1;
			if (run > bc)
			{
				bc = run;
				best = i + 1 - run;
			}
		}
		for (k = 0; k < n; k++)
			out[k] = r[(size_t)PM_REC * best + 1 + k];
		*idx = (int)r[(size_t)PM_REC * best + 1 + n];
		qaws_internal_free(r);
		break;
	}
	case QAWS_PROMOTE_MEDIAN:
		for (k = 0; k < n; k++)
		{
			for (i = 0; i < ns; i++)
				tmp[i] = vals[(size_t)n * src[i] + k];
			qsort(tmp, ns, sizeof(double), pm_cmp_double);
			out[k] = tmp[ns / 2];
		}
		break;
	case QAWS_PROMOTE_SUM:
	case QAWS_PROMOTE_SUM_OF_SQUARES:
	case QAWS_PROMOTE_RMS:
		for (k = 0; k < n; k++)
		{
			double s = 0;
			for (i = 0; i < ns; i++)
			{
				double v = vals[(size_t)n * src[i] + k];
				s += method == QAWS_PROMOTE_SUM ? v : v * v;
			}
			out[k] = method == QAWS_PROMOTE_RMS ? sqrt(s / ns) : s;
		}
		break;
	case QAWS_PROMOTE_FIRST:
	case QAWS_PROMOTE_LAST:
	{
		unsigned int b = 0;
		for (i = 1; i < ns; i++)
			if (method == QAWS_PROMOTE_FIRST ? src[i] < src[b] : src[i] > src[b])
				b = i;
		for (k = 0; k < n; k++)
			out[k] = vals[(size_t)n * src[b] + k];
		*idx = (int)src[b];
		break;
	}
	}
}

static unsigned int pm_count(qaws_mesh_3d const* m, unsigned int el)
{
	return el == QAWS_MESH_VERTEX ? qaws_mesh_3d_get_vertex_count(m) : el == QAWS_MESH_CORNER ? qaws_mesh_3d_get_corner_count(m)
		: el == QAWS_MESH_FACE ? qaws_mesh_3d_get_face_count(m) : qaws_mesh_3d_get_edge_count(m);
}

/* the (target, source) pairs of the relation to <- from, sorted and unique, into *out (target << 32 | source), their count;
   ~0u on allocation failure */
static unsigned int pm_pairs(unsigned long long** out, qaws_mesh_3d const* m, unsigned int to, unsigned int from)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), nc = qaws_mesh_3d_get_corner_count(m), ne = qaws_mesh_3d_get_edge_count(m), i, c, np = 0, w;
	unsigned int *fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
	unsigned int *fs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nf + 1));
	unsigned int *cf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
	unsigned int *cn = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
	unsigned int *cp = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
	unsigned int *ce = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)nc + 1));
	unsigned int *e = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)ne + 1));
	unsigned long long* key = NULL;
	unsigned long long* p = (unsigned long long*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned long long), 4, (size_t)nc + 2 * (size_t)ne + 1));
	*out = NULL;
	if (!fv || !fs || !cf || !cn || !cp || !ce || !e || !p)
		goto fail;
	qaws_mesh_3d_get_faces(fv, fs, m);
	qaws_mesh_3d_get_edges(e, m);
	if (ne && !(key = cy_keys(e, ne)))
		goto fail;
	for (i = 0, c = 0; i < nf; c += fs[i++])
	{
		unsigned int j;
		for (j = 0; j < fs[i]; j++)
		{
			cf[c + j] = i;
			cn[c + j] = c + (j + 1) % fs[i];
			cp[c + j] = c + (j + fs[i] - 1) % fs[i];
		}
	}
	for (c = 0; c < nc; c++)
		ce[c] = ne ? cy_find(key, ne, fv[c], fv[cn[c]]) : 0;
#define PM_PAIR(t, s) (p[np++] = (unsigned long long)(t) << 32 | (unsigned long long)(s))
	if (to == from)
	{
		unsigned int cnt = pm_count(m, to);
		for (i = 0; i < cnt; i++)
			PM_PAIR(i, i);
	}
	else if (to == QAWS_MESH_VERTEX && from == QAWS_MESH_EDGE)
		for (i = 0; i < ne; i++)
		{
			PM_PAIR(e[2 * (size_t)i], i);
			PM_PAIR(e[2 * (size_t)i + 1], i);
		}
	else if (to == QAWS_MESH_EDGE && from == QAWS_MESH_VERTEX)
		for (i = 0; i < ne; i++)
		{
			PM_PAIR(i, e[2 * (size_t)i]);
			PM_PAIR(i, e[2 * (size_t)i + 1]);
		}
	else
		for (c = 0; c < nc; c++)
		{
			unsigned int v = fv[c], f = cf[c], eo = ce[c], ei = ce[cp[c]];
			if (to == QAWS_MESH_VERTEX)
				PM_PAIR(v, from == QAWS_MESH_CORNER ? c : f);
			else if (to == QAWS_MESH_CORNER)
			{
				if (from == QAWS_MESH_EDGE)
				{
					PM_PAIR(c, eo);
					PM_PAIR(c, ei);
				}
				else
					PM_PAIR(c, from == QAWS_MESH_VERTEX ? v : f);
			}
			else if (to == QAWS_MESH_FACE)
				PM_PAIR(f, from == QAWS_MESH_VERTEX ? v : from == QAWS_MESH_CORNER ? c : eo);
			else
			{
				if (from == QAWS_MESH_CORNER)
				{
					PM_PAIR(eo, c);
					PM_PAIR(eo, cn[c]);
				}
				else
					PM_PAIR(eo, f);
			}
		}
#undef PM_PAIR
	qsort(p, np, sizeof(unsigned long long), pm_cmp_u64);
	for (i = 0, w = 0; i < np; i++)
		if (!w || p[i] != p[w - 1])
			p[w++] = p[i];
	qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(cf); qaws_internal_free(cn); qaws_internal_free(cp); qaws_internal_free(ce);
	qaws_internal_free(e); qaws_internal_free(key);
	*out = p;
	return w;
fail:
	qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(cf); qaws_internal_free(cn); qaws_internal_free(cp); qaws_internal_free(ce);
	qaws_internal_free(e); qaws_internal_free(key); qaws_internal_free(p);
	return ~0u;
}

static int pm_cmp_piece(void const* a, void const* b)
{
	double const *x = (double const*)a, *y = (double const*)b;
	return x[0] < y[0] ? -1 : x[0] > y[0] ? 1 : x[1] < y[1] ? -1 : x[1] > y[1];
}

qaws_status qaws_mesh_3d_promote_attribute(qaws_mesh_3d* mesh, qaws_attribute_promote_desc const* desc)
{
	qaws_attribute_desc sd, od;
	unsigned int id, pid = 0, from, to, n, cs, ct, np, i, k, j;
	char const* out_name;
	qaws_scalar *sv = NULL, *ov = NULL, *pv = NULL, *iv = NULL;
	double *dv = NULL, *w = NULL, *tmp = NULL, *pc = NULL, y[16];
	unsigned long long* pairs = NULL;
	unsigned int *off = NULL, *gather = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!mesh || !desc || !desc->name || qaws_mesh_3d_find_attribute(&id, mesh, desc->name) != QAWS_STATUS_OK ||
		(unsigned int)desc->element > QAWS_MESH_EDGE || (unsigned int)desc->method > QAWS_PROMOTE_LAST)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_mesh_3d_get_attribute_desc(&sd, mesh, id);
	from = (unsigned int)sd.element;
	to = (unsigned int)desc->element;
	n = sd.components;
	out_name = desc->new_name ? desc->new_name : desc->name;
	if ((from == to && !desc->piece) || (strcmp(out_name, desc->name) && qaws_mesh_3d_find_attribute(&k, mesh, out_name) == QAWS_STATUS_OK))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (desc->piece)
	{
		qaws_attribute_desc pd;
		if (qaws_mesh_3d_find_attribute(&pid, mesh, desc->piece) != QAWS_STATUS_OK)
			return QAWS_STATUS_INVALID_ARGUMENT;
		qaws_mesh_3d_get_attribute_desc(&pd, mesh, pid);
		if ((unsigned int)pd.element != to || pd.components != 1)
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	cs = pm_count(mesh, from);
	ct = pm_count(mesh, to);
	sv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), n, (size_t)cs + 1));
	dv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, (size_t)cs + 1));
	ov = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), n, (size_t)ct + 1));
	iv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)ct + 1));
	off = (unsigned int*)qaws_internal_scratch_zalloc((size_t)ct + 2, sizeof(unsigned int));
	if (!sv || !dv || !ov || !iv || !off)
		goto done;
	qaws_mesh_3d_get_attribute(sv, mesh, id);
	for (i = 0; i < n * cs; i++)
		dv[i] = (double)sv[i];
	if ((np = pm_pairs(&pairs, mesh, to, from)) == ~0u)
		goto done;
	/* (the pairs' sources per target: a CSR over the sorted pairs) */
	gather = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)np + 1));
	w = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)np + 1));
	tmp = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)np + 1));
	if (!gather || !w || !tmp)
		goto done;
	for (i = 0; i < np; i++)
	{
		off[(unsigned int)(pairs[i] >> 32) + 1]++;
		gather[i] = (unsigned int)pairs[i];
	}
	for (i = 0; i < ct; i++)
		off[i + 1] += off[i];
	if (!desc->piece)
		for (i = 0; i < ct; i++)
		{
			int ix;
			pm_combine(y, &ix, dv, gather + off[i], off[i + 1] - off[i], n, (unsigned int)sd.rule, desc->method, w, tmp);
			for (k = 0; k < n; k++)
				ov[(size_t)n * i + k] = (qaws_scalar)y[k];
			iv[i] = (qaws_scalar)ix;
		}
	else
	{
		/* pieces: the targets by their piece value; each piece's sources together */
		unsigned int* src = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)np + 1));
		pv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), (size_t)ct + 1));
		pc = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)ct + 1));
		if (!src || !pv || !pc)
		{
			qaws_internal_free(src);
			goto done;
		}
		qaws_mesh_3d_get_attribute(pv, mesh, pid);
		for (i = 0; i < ct; i++)
		{
			pc[2 * (size_t)i] = (double)pv[i];
			pc[2 * (size_t)i + 1] = i;
		}
		qsort(pc, ct, 2 * sizeof(double), pm_cmp_piece);
		for (i = 0; i < ct; i = j)
		{
			unsigned int ns = 0, u, m2;
			int ix;
			for (j = i; j < ct && pc[2 * (size_t)j] == pc[2 * (size_t)i]; j++)
			{
				unsigned int t = (unsigned int)pc[2 * (size_t)j + 1];
				for (u = off[t]; u < off[t + 1]; u++)
					src[ns++] = gather[u];
			}
			/* (each source once, in order) */
			for (u = 1; u < ns; u++)
				for (m2 = u; m2 > 0 && src[m2 - 1] > src[m2]; m2--)
				{
					unsigned int s2 = src[m2];
					src[m2] = src[m2 - 1];
					src[m2 - 1] = s2;
				}
			for (u = 0, m2 = 0; u < ns; u++)
				if (!m2 || src[u] != src[m2 - 1])
					src[m2++] = src[u];
			pm_combine(y, &ix, dv, src, m2, n, (unsigned int)sd.rule, desc->method, w, tmp);
			for (u = i; u < j; u++)
			{
				unsigned int t = (unsigned int)pc[2 * (size_t)u + 1];
				for (k = 0; k < n; k++)
					ov[(size_t)n * t + k] = (qaws_scalar)y[k];
				iv[t] = (qaws_scalar)ix;
			}
		}
		qaws_internal_free(src);
	}
	/* the result: the original removed first when the name is one, or when asked */
	if (!strcmp(out_name, desc->name) || desc->delete_original)
	{
		if ((st = qaws_mesh_3d_find_attribute(&id, mesh, desc->name)) != QAWS_STATUS_OK || (st = qaws_mesh_3d_remove_attribute(mesh, id)) != QAWS_STATUS_OK)
			goto done;
	}
	memset(&od, 0, sizeof od);
	od.name = out_name;
	od.element = (qaws_mesh_element)to;
	od.components = n;
	od.rule = sd.rule;
	if ((st = qaws_mesh_3d_add_attribute(&k, mesh, ov, &od)) != QAWS_STATUS_OK)
		goto done;
	if (desc->index_name)
	{
		unsigned int xid;
		if (qaws_mesh_3d_find_attribute(&xid, mesh, desc->index_name) == QAWS_STATUS_OK)
			qaws_mesh_3d_remove_attribute(mesh, xid);
		od.name = desc->index_name;
		od.components = 1;
		od.rule = QAWS_ATTRIBUTE_NEAREST;
		st = qaws_mesh_3d_add_attribute(&k, mesh, iv, &od);
	}
done:
	qaws_internal_free(sv); qaws_internal_free(dv); qaws_internal_free(ov); qaws_internal_free(iv); qaws_internal_free(off); qaws_internal_free(pairs);
	qaws_internal_free(gather); qaws_internal_free(w); qaws_internal_free(tmp); qaws_internal_free(pv); qaws_internal_free(pc);
	return st;
}

qaws_status qaws_mesh_3d_reduce_attribute(qaws_scalar* out_values, int* out_index, qaws_mesh_3d const* mesh, char const* name, qaws_promote_method method)
{
	qaws_attribute_desc sd;
	unsigned int id, cs, n, i;
	qaws_scalar* sv;
	double *dv, *w, *tmp, y[16];
	unsigned int* all;
	int ix;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out_values || !mesh || !name || qaws_mesh_3d_find_attribute(&id, mesh, name) != QAWS_STATUS_OK || (unsigned int)method > QAWS_PROMOTE_LAST)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_mesh_3d_get_attribute_desc(&sd, mesh, id);
	n = sd.components;
	cs = pm_count(mesh, (unsigned int)sd.element);
	sv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), n, (size_t)cs + 1));
	dv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, (size_t)cs + 1));
	w = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)cs + 1));
	tmp = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (size_t)cs + 1));
	all = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)cs + 1));
	if (sv && dv && w && tmp && all)
	{
		qaws_mesh_3d_get_attribute(sv, mesh, id);
		for (i = 0; i < n * cs; i++)
			dv[i] = (double)sv[i];
		for (i = 0; i < cs; i++)
			all[i] = i;
		pm_combine(y, &ix, dv, all, cs, n, (unsigned int)sd.rule, method, w, tmp);
		for (i = 0; i < n; i++)
			out_values[i] = (qaws_scalar)y[i];
		if (out_index)
			*out_index = ix;
		st = QAWS_STATUS_OK;
	}
	qaws_internal_free(sv); qaws_internal_free(dv); qaws_internal_free(w); qaws_internal_free(tmp); qaws_internal_free(all);
	return st;
}
