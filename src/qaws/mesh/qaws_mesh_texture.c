/*
 * Texture attributes (qaws_mesh_texture.h): triangles rasterized in (u, v)
 * at the texels' centres (each texel its triangle and weights), valued by
 * the attribute's rule from the triangle's corners or from a source mesh's
 * closest point; reading bilinear by the rule.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_texture.h"
#include "../qaws_mesh_bvh.h"
#include "../qaws_inspect.h"
#include "qaws_mesh_internal.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* the mesh's faces, each corner's (u, v) and which faces count (desc's id) */
typedef struct tx_faces
{
	unsigned int nf, nc, *fv, *fs, *fo;   /* faces' vertices and sizes, each face's first corner */
	double* uv;                           /* 2 per corner */
	unsigned char* use;                   /* per face */
} tx_faces;

/* a texel in a triangle: the texel, the face, its three corners and the centre's weights on them */
typedef struct tx_hit
{
	size_t tx;
	unsigned int f, c[3];
	double w[3];
} tx_hit;

static void tx_faces_free(tx_faces* t)
{
	qaws_internal_free(t->fv); qaws_internal_free(t->fs); qaws_internal_free(t->fo); qaws_internal_free(t->uv); qaws_internal_free(t->use);
}

static qaws_status tx_faces_read(tx_faces* t, qaws_mesh_3d const* m, qaws_texture_desc const* desc)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), id = 0, i, c;
	qaws_attribute_desc ad;
	qaws_scalar* v = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	qaws_mesh_location* loc = NULL;
	int foot = desc && desc->footpoints;
	memset(t, 0, sizeof *t);
	t->nf = qaws_mesh_3d_get_face_count(m);
	t->nc = qaws_mesh_3d_get_corner_count(m);
	if (foot)
	{
		/* (the patch's own (u, v), from the footpoints, over its domain) */
		if (!((double)desc->domain[1] != (double)desc->domain[0]) || !((double)desc->domain[3] != (double)desc->domain[2]))
			return QAWS_STATUS_INVALID_ARGUMENT;
		loc = (qaws_mesh_location*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_mesh_location), (size_t)nv + 1));
		if (!loc)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		if (qaws_mesh_3d_get_footpoints(loc, m) != QAWS_STATUS_OK)
		{
			qaws_internal_free(loc);
			return QAWS_STATUS_INVALID_ARGUMENT;
		}
		ad.element = QAWS_MESH_VERTEX;
	}
	else
	{
		if (qaws_mesh_3d_find_attribute(&id, m, desc && desc->uv ? desc->uv : "uv") != QAWS_STATUS_OK)
			return QAWS_STATUS_INVALID_ARGUMENT;
		qaws_mesh_3d_get_attribute_desc(&ad, m, id);
		if (ad.components != 2 || (ad.element != QAWS_MESH_CORNER && ad.element != QAWS_MESH_VERTEX))
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	t->fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)t->nc + 1));
	t->fs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)t->nf + 1));
	t->fo = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)t->nf + 1));
	t->uv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 2, (size_t)t->nc + 1));
	t->use = (unsigned char*)qaws_internal_scratch_malloc((size_t)t->nf + 1);
	v = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 2, (size_t)(nv > t->nc ? nv : t->nc) + 1));
	if (!t->fv || !t->fs || !t->fo || !t->uv || !t->use || !v)
		goto fail;
	qaws_mesh_3d_get_faces(t->fv, t->fs, m);
	for (i = 0, c = 0; i < t->nf; c += t->fs[i++])
		t->fo[i] = c;
	memset(t->use, 1, t->nf);
	if (foot)
	{
		double u0 = (double)desc->domain[0], du = (double)desc->domain[1] - u0, v0 = (double)desc->domain[2], dv = (double)desc->domain[3] - v0;
		for (c = 0; c < t->nc; c++)
		{
			qaws_mesh_location const* l = loc + t->fv[c];
			t->uv[2 * (size_t)c] = ((double)l->a - u0) / du;
			t->uv[2 * (size_t)c + 1] = ((double)l->b - v0) / dv;
		}
		/* (a face counts when its corners all lie on the patch; given its surface, when one does, the others projected) */
		for (i = 0; i < t->nf; i++)
		{
			unsigned int on = 0;
			for (c = t->fo[i]; c < t->fo[i] + t->fs[i]; c++)
				on += loc[t->fv[c]].patch == desc->patch;
			if (on == t->fs[i])
				continue;
			if (!on || !desc->surface)
			{
				t->use[i] = 0;
				continue;
			}
			for (c = t->fo[i]; c < t->fo[i] + t->fs[i]; c++)
				if (loc[t->fv[c]].patch != desc->patch)
				{
					qaws_scalar pu, pv;
					qaws_vec3 cp;
					if (qaws_surface_find_closest_point(&pu, &pv, &cp, desc->surface, qaws_mesh_3d_get_vertex(m, t->fv[c])) != QAWS_STATUS_OK)
					{
						t->use[i] = 0;
						break;
					}
					t->uv[2 * (size_t)c] = ((double)pu - u0) / du;
					t->uv[2 * (size_t)c + 1] = ((double)pv - v0) / dv;
				}
		}
	}
	else
	{
		qaws_mesh_3d_get_attribute(v, m, id);
		for (c = 0; c < t->nc; c++)
		{
			unsigned int e = ad.element == QAWS_MESH_CORNER ? c : t->fv[c];
			t->uv[2 * (size_t)c] = (double)v[2 * (size_t)e];
			t->uv[2 * (size_t)c + 1] = (double)v[2 * (size_t)e + 1];
		}
	}
	if (desc && desc->texture_id)
	{
		if (qaws_mesh_3d_find_attribute(&id, m, desc->texture_id) != QAWS_STATUS_OK)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto fail;
		}
		qaws_mesh_3d_get_attribute_desc(&ad, m, id);
		if (ad.components != 1 || ad.element != QAWS_MESH_FACE)
		{
			st = QAWS_STATUS_INVALID_ARGUMENT;
			goto fail;
		}
		qaws_mesh_3d_get_attribute(v, m, id);
		for (i = 0; i < t->nf; i++)
			t->use[i] &= (unsigned char)((double)v[i] == (double)desc->id);
	}
	qaws_internal_free(v);
	qaws_internal_free(loc);
	return QAWS_STATUS_OK;
fail:
	qaws_internal_free(v);
	qaws_internal_free(loc);
	tx_faces_free(t);
	return st;
}

static int tx_check(qaws_texture const* x)
{
	return x && x->texels && x->width && x->height && x->components >= 1 && x->components <= 16;
}

/* a texel index along an axis of n: wrapped when repeated, else ~0u outside */
static unsigned int tx_index(long long i, unsigned int n, int repeat)
{
	if (repeat)
	{
		long long r = i % (long long)n;
		return (unsigned int)(r < 0 ? r + n : r);
	}
	return i < 0 || i >= (long long)n ? ~0u : (unsigned int)i;
}

/* texture's value at (u, v) into out (components): the four nearest centres' by the rule (clamped or repeated) */
static void tx_read(double* out, qaws_texture const* x, double u, double v)
{
	double fx = u * x->width - 0.5, fy = v * x->height - 0.5, vals[64], w[4];
	long long i0 = (long long)floor(fx), j0 = (long long)floor(fy);
	unsigned int pos[4] = { 0, 1, 2, 3 }, k, q, n = x->components;
	fx -= (double)i0;
	fy -= (double)j0;
	w[0] = (1 - fx) * (1 - fy); w[1] = fx * (1 - fy); w[2] = (1 - fx) * fy; w[3] = fx * fy;
	for (q = 0; q < 4; q++)
	{
		long long ii = i0 + (q & 1), jj = j0 + (q >> 1);
		unsigned int a, b;
		if (x->repeat_u)
			a = tx_index(ii, x->width, 1);
		else
			a = ii < 0 ? 0 : ii >= (long long)x->width ? x->width - 1 : (unsigned int)ii;
		if (x->repeat_v)
			b = tx_index(jj, x->height, 1);
		else
			b = jj < 0 ? 0 : jj >= (long long)x->height ? x->height - 1 : (unsigned int)jj;
		for (k = 0; k < n; k++)
			vals[(size_t)n * q + k] = (double)x->texels[((size_t)b * x->width + a) * n + k];
	}
	hm_blend(out, vals, n, pos, w, 4, n, (unsigned int)x->rule);
}

/* the texels whose centres lie in the faces' triangles (fanned from their first corner; a triangle across a repeated
   axis' edge unwrapped to its first corner's copy), each once (a later triangle over it replaces the earlier) */
static qaws_status tx_raster(tx_hit** out, size_t* out_count, tx_faces const* t, qaws_texture const* x)
{
	unsigned int W = x->width, H = x->height, f;
	size_t nt = (size_t)W * H, n = 0, cap = 0;
	size_t* at = (size_t*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(size_t), nt + 1));
	tx_hit* h = NULL;
	*out = NULL;
	*out_count = 0;
	if (!at)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	memset(at, 0xFF, sizeof(size_t) * nt);
	for (f = 0; f < t->nf; f++)
	{
		unsigned int j;
		if (!t->use[f])
			continue;
		for (j = 1; j + 1 < t->fs[f]; j++)
		{
			unsigned int c[3] = { t->fo[f], t->fo[f] + j, t->fo[f] + j + 1 }, q;
			double U[3][2], lo[2], hi[2], den;
			long long i0, i1, j0, j1, ii, jj;
			for (q = 0; q < 3; q++)
			{
				U[q][0] = t->uv[2 * (size_t)c[q]];
				U[q][1] = t->uv[2 * (size_t)c[q] + 1];
				if (q && x->repeat_u)
					U[q][0] -= floor(U[q][0] - U[0][0] + 0.5);
				if (q && x->repeat_v)
					U[q][1] -= floor(U[q][1] - U[0][1] + 0.5);
			}
			lo[0] = fmin(U[0][0], fmin(U[1][0], U[2][0])); hi[0] = fmax(U[0][0], fmax(U[1][0], U[2][0]));
			lo[1] = fmin(U[0][1], fmin(U[1][1], U[2][1])); hi[1] = fmax(U[0][1], fmax(U[1][1], U[2][1]));
			den = (U[1][0] - U[0][0]) * (U[2][1] - U[0][1]) - (U[2][0] - U[0][0]) * (U[1][1] - U[0][1]);
			if (!(fabs(den) > 0) || !isfinite(den))
				continue;
			i0 = (long long)ceil(lo[0] * W - 0.5); i1 = (long long)floor(hi[0] * W - 0.5);
			j0 = (long long)ceil(lo[1] * H - 0.5); j1 = (long long)floor(hi[1] * H - 0.5);
			if (i1 - i0 > 4 * (long long)W || j1 - j0 > 4 * (long long)H)
				continue;   /* (a triangle over many copies of the square: a broken map) */
			for (jj = j0; jj <= j1; jj++)
				for (ii = i0; ii <= i1; ii++)
				{
					unsigned int a = tx_index(ii, W, x->repeat_u), b = tx_index(jj, H, x->repeat_v);
					double px = (ii + 0.5) / W, py = (jj + 0.5) / H, w1, w2;
					size_t tx, k;
					if (a == ~0u || b == ~0u)
						continue;
					w1 = ((px - U[0][0]) * (U[2][1] - U[0][1]) - (U[2][0] - U[0][0]) * (py - U[0][1])) / den;
					w2 = ((U[1][0] - U[0][0]) * (py - U[0][1]) - (px - U[0][0]) * (U[1][1] - U[0][1])) / den;
					if (w1 < -1e-12 || w2 < -1e-12 || 1 - w1 - w2 < -1e-12)
						continue;
					tx = (size_t)b * W + a;
					k = at[tx];
					if (k == (size_t)-1)
					{
						if (n + 1 > cap)
						{
							size_t c2 = cap ? 2 * cap : 1024;
							tx_hit* g = (tx_hit*)(h ? qaws_internal_realloc(h, qaws_internal_mul2(sizeof(tx_hit), c2))
							                        : qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(tx_hit), c2)));
							if (!g)
							{
								qaws_internal_free(h);
								qaws_internal_free(at);
								return QAWS_STATUS_ALLOCATION_FAILURE;
							}
							h = g;
							cap = c2;
						}
						k = at[tx] = n++;
					}
					h[k].tx = tx;
					h[k].f = f;
					memcpy(h[k].c, c, sizeof c);
					h[k].w[0] = 1 - w1 - w2;
					h[k].w[1] = w1;
					h[k].w[2] = w2;
				}
		}
	}
	qaws_internal_free(at);
	*out = h;
	*out_count = n;
	return QAWS_STATUS_OK;
}

/* gutter rings grown past the written texels (ring: 1 written, else 0), each new texel the rule's mean of its written
   neighbours (8 about it); ring left 2 on them */
static qaws_status tx_gutter(qaws_texture* x, unsigned char* ring, unsigned int gutter)
{
	unsigned int W = x->width, H = x->height, n = x->components, r, k;
	size_t nt = (size_t)W * H;
	unsigned int* grow;
	double* gval;
	if (!gutter)
		return QAWS_STATUS_OK;
	grow = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), nt + 1));
	gval = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, nt + 1));
	if (!grow || !gval)
	{
		qaws_internal_free(grow);
		qaws_internal_free(gval);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (r = 0; r < gutter; r++)
	{
		size_t tx, ng = 0;
		for (tx = 0; tx < nt; tx++)
		{
			unsigned int a = (unsigned int)(tx % W), b = (unsigned int)(tx / W), pos[8], m = 0, q;
			double w[8], vals[8 * 16];
			int di, dj;
			if (ring[tx])
				continue;
			for (dj = -1; dj <= 1; dj++)
				for (di = -1; di <= 1; di++)
				{
					unsigned int a2 = tx_index((long long)a + di, W, x->repeat_u), b2 = tx_index((long long)b + dj, H, x->repeat_v);
					size_t t2;
					if ((!di && !dj) || a2 == ~0u || b2 == ~0u)
						continue;
					t2 = (size_t)b2 * W + a2;
					if (!ring[t2])
						continue;   /* (this ring's new texels are written after the scan) */
					for (q = 0; q < n; q++)
						vals[(size_t)n * m + q] = (double)x->texels[t2 * n + q];
					pos[m] = m;
					m++;
				}
			if (!m)
				continue;
			for (q = 0; q < m; q++)
				w[q] = 1.0 / m;
			hm_blend(gval + (size_t)n * ng, vals, n, pos, w, m, n, (unsigned int)x->rule);
			grow[ng++] = (unsigned int)tx;
		}
		for (tx = 0; tx < ng; tx++)
		{
			for (k = 0; k < n; k++)
				x->texels[(size_t)grow[tx] * n + k] = (qaws_scalar)gval[tx * n + k];
			ring[grow[tx]] = 2;
		}
		if (!ng)
			break;
	}
	qaws_internal_free(grow);
	qaws_internal_free(gval);
	return QAWS_STATUS_OK;
}

/* an attribute's values as doubles (scratch; NULL on failure), *count its elements */
static double* tx_values(qaws_mesh_3d const* m, unsigned int id, qaws_attribute_desc const* ad, unsigned int* count)
{
	unsigned int cnt = ad->element == QAWS_MESH_VERTEX ? qaws_mesh_3d_get_vertex_count(m) : ad->element == QAWS_MESH_CORNER
		? qaws_mesh_3d_get_corner_count(m) : qaws_mesh_3d_get_face_count(m), k;
	double* dv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), ad->components, (size_t)cnt + 1));
	qaws_scalar* sv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), ad->components, (size_t)cnt + 1));
	if (!dv || !sv)
	{
		qaws_internal_free(dv);
		qaws_internal_free(sv);
		return NULL;
	}
	qaws_mesh_3d_get_attribute(sv, m, id);
	for (k = 0; k < ad->components * cnt; k++)
		dv[k] = (double)sv[k];
	qaws_internal_free(sv);
	*count = cnt;
	return dv;
}

qaws_status tx_write(unsigned char* out_coverage, qaws_texture* x, size_t const* tx, double const* vals, unsigned char const* keep, size_t n,
	unsigned int gutter)
{
	size_t nt = (size_t)x->width * x->height, i;
	unsigned int k, nc = x->components;
	unsigned char* ring = (unsigned char*)qaws_internal_scratch_zalloc(nt + 1, 1);
	qaws_status st;
	if (!ring)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < n; i++)
	{
		if (keep && !keep[i])
			continue;
		for (k = 0; k < nc; k++)
			x->texels[tx[i] * nc + k] = (qaws_scalar)vals[i * nc + k];
		ring[tx[i]] = 1;
	}
	st = tx_gutter(x, ring, gutter);
	if (st == QAWS_STATUS_OK && out_coverage)
		memcpy(out_coverage, ring, nt);
	qaws_internal_free(ring);
	return st;
}

/* the hits' values written (tx_write over their texels) */
static qaws_status tx_hits_write(unsigned char* out_coverage, qaws_texture* x, tx_hit const* h, size_t nh, unsigned char const* keep, double const* vals,
	unsigned int gutter)
{
	size_t* tx = (size_t*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(size_t), nh + 1)), i;
	qaws_status st;
	if (!tx)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < nh; i++)
		tx[i] = h[i].tx;
	st = tx_write(out_coverage, x, tx, vals, keep, nh, gutter);
	qaws_internal_free(tx);
	return st;
}

qaws_status tx_texel_points(double** out_points, size_t** out_tx, size_t* out_count, qaws_texture const* texture, qaws_mesh_3d const* mesh,
	qaws_texture_desc const* desc)
{
	tx_faces t;
	tx_hit* h = NULL;
	size_t nh = 0, i;
	unsigned int k, q;
	double* P = NULL;
	size_t* tx = NULL;
	qaws_scalar* x = NULL;
	qaws_status st;
	*out_points = NULL;
	*out_tx = NULL;
	*out_count = 0;
	if (!tx_check(texture) || !mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = tx_faces_read(&t, mesh, desc)) != QAWS_STATUS_OK)
		return st;
	x = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)qaws_mesh_3d_get_vertex_count(mesh) + 1));
	if (!x)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	if ((st = tx_raster(&h, &nh, &t, texture)) != QAWS_STATUS_OK)
		goto done;
	qaws_mesh_3d_get_positions(x, mesh);
	P = (double*)qaws_internal_malloc(qaws_internal_mul3(sizeof(double), 9, nh + 1));
	tx = (size_t*)qaws_internal_malloc(qaws_internal_mul2(sizeof(size_t), nh + 1));
	if (!P || !tx)
	{
		qaws_internal_free(P);
		qaws_internal_free(tx);
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < nh; i++)
	{
		double V[3][3];
		for (q = 0; q < 3; q++)
			for (k = 0; k < 3; k++)
				V[q][k] = (double)x[3 * (size_t)t.fv[h[i].c[q]] + k];
		for (k = 0; k < 3; k++)
		{
			P[9 * i + k] = h[i].w[0] * V[0][k] + h[i].w[1] * V[1][k] + h[i].w[2] * V[2][k];
			P[9 * i + 3 + k] = V[1][k] - V[0][k];
			P[9 * i + 6 + k] = V[2][k] - V[0][k];
		}
		tx[i] = h[i].tx;
	}
	*out_points = P;
	*out_tx = tx;
	*out_count = nh;
	st = QAWS_STATUS_OK;
done:
	tx_faces_free(&t);
	qaws_internal_free(h);
	qaws_internal_free(x);
	return st;
}

qaws_status qaws_mesh_3d_bake_texture(unsigned char* out_coverage, qaws_texture* texture, qaws_mesh_3d const* mesh, char const* name,
	qaws_texture_desc const* desc)
{
	tx_faces t;
	qaws_attribute_desc ad;
	unsigned int id, n, cnt, q;
	size_t nh = 0, i;
	double *dv = NULL, *vals = NULL;
	tx_hit* h = NULL;
	qaws_status st;
	if (!tx_check(texture) || !mesh || !name || qaws_mesh_3d_find_attribute(&id, mesh, name) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_mesh_3d_get_attribute_desc(&ad, mesh, id);
	if (ad.element == QAWS_MESH_EDGE || ad.components != texture->components)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = tx_faces_read(&t, mesh, desc)) != QAWS_STATUS_OK)
		return st;
	n = ad.components;
	texture->rule = ad.rule;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!(dv = tx_values(mesh, id, &ad, &cnt)) || (st = tx_raster(&h, &nh, &t, texture)) != QAWS_STATUS_OK)
		goto done;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!(vals = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, nh + 1))))
		goto done;
	for (i = 0; i < nh; i++)
	{
		unsigned int pos[3];
		for (q = 0; q < 3; q++)
			pos[q] = ad.element == QAWS_MESH_VERTEX ? t.fv[h[i].c[q]] : ad.element == QAWS_MESH_CORNER ? h[i].c[q] : h[i].f;
		hm_blend(vals + i * n, dv, n, pos, h[i].w, 3, n, (unsigned int)ad.rule);
	}
	st = tx_hits_write(out_coverage, texture, h, nh, NULL, vals, desc ? desc->gutter : 0);
done:
	tx_faces_free(&t);
	qaws_internal_free(dv); qaws_internal_free(vals); qaws_internal_free(h);
	return st;
}

qaws_status qaws_mesh_3d_bake_texture_from(unsigned char* out_coverage, qaws_texture* texture, qaws_mesh_3d const* mesh, qaws_mesh_3d const* source,
	char const* name, qaws_texture_desc const* desc)
{
	tx_faces t;
	qaws_attribute_desc ad;
	unsigned int id, n, cnt, q, snv, snt, snf, snc, f, j, k;
	size_t nh = 0, i;
	double *dv = NULL, *vals = NULL;
	tx_hit* h = NULL;
	qaws_scalar *x = NULL, *sx = NULL, *pts = NULL;
	unsigned int *stri = NULL, *sfv = NULL, *sfs = NULL, *tf = NULL, *tc = NULL;
	unsigned char* keep = NULL;
	qaws_mesh_bvh* bvh = NULL;
	qaws_mesh_hit* hits = NULL;
	qaws_status st;
	if (!tx_check(texture) || !mesh || !source || !name || qaws_mesh_3d_find_attribute(&id, source, name) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_mesh_3d_get_attribute_desc(&ad, source, id);
	if (ad.element == QAWS_MESH_EDGE || ad.components != texture->components)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = tx_faces_read(&t, mesh, desc)) != QAWS_STATUS_OK)
		return st;
	n = ad.components;
	texture->rule = ad.rule;
	snv = qaws_mesh_3d_get_vertex_count(source);
	snt = qaws_mesh_3d_get_triangle_count(source);
	snf = qaws_mesh_3d_get_face_count(source);
	snc = qaws_mesh_3d_get_corner_count(source);
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	x = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)qaws_mesh_3d_get_vertex_count(mesh) + 1));
	sx = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, (size_t)snv + 1));
	stri = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)snt + 1));
	sfv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)snc + 1));
	sfs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)snf + 1));
	tf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (size_t)snt + 1));
	tc = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)snt + 1));
	if (!x || !sx || !stri || !sfv || !sfs || !tf || !tc || !(dv = tx_values(source, id, &ad, &cnt)))
		goto done;
	qaws_mesh_3d_get_positions(x, mesh);
	qaws_mesh_3d_get_positions(sx, source);
	qaws_mesh_3d_get_triangles(stri, source);
	qaws_mesh_3d_get_faces(sfv, sfs, source);
	/* (the source's triangles back to their faces' corners: faces fanned from their first corner) */
	for (f = 0, j = 0, k = 0; f < snf; k += sfs[f++])
		for (q = 1; q + 1 < sfs[f] && j < snt; q++, j++)
		{
			tf[j] = f;
			tc[3 * (size_t)j] = k;
			tc[3 * (size_t)j + 1] = k + q;
			tc[3 * (size_t)j + 2] = k + q + 1;
		}
	{
		qaws_mesh_bvh_desc bd;
		memset(&bd, 0, sizeof bd);
		bd.positions = sx;
		bd.vertex_count = snv;
		bd.dimension = 3;
		bd.simplices = stri;
		bd.simplex_size = 3;
		bd.simplex_count = snt;
		if ((st = qaws_mesh_bvh_create(&bvh, &bd)) != QAWS_STATUS_OK || (st = tx_raster(&h, &nh, &t, texture)) != QAWS_STATUS_OK)
			goto done;
	}
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	pts = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), 3, nh + 1));
	hits = (qaws_mesh_hit*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_mesh_hit), nh + 1));
	vals = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, nh + 1));
	keep = (unsigned char*)qaws_internal_scratch_zalloc(nh + 1, 1);
	if (!pts || !hits || !vals || !keep)
		goto done;
	/* each texel's point on the mesh, the source's closest point to it, the source's value there */
	for (i = 0; i < nh; i++)
		for (k = 0; k < 3; k++)
		{
			double p = 0;
			for (q = 0; q < 3; q++)
				p += h[i].w[q] * (double)x[3 * (size_t)t.fv[h[i].c[q]] + k];
			pts[3 * i + k] = (qaws_scalar)p;
		}
	if (nh > 0xFFFFFFFFu)
		goto done;
	if ((st = qaws_mesh_bvh_closest(hits, bvh, pts, (unsigned int)nh, desc ? desc->max_distance : 0, NULL)) != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < nh; i++)
	{
		unsigned int s = hits[i].simplex, pos[3];
		double w[3];
		if (hits[i].element == QAWS_MESH_NONE || s >= snt)
			continue;
		for (q = 0; q < 3; q++)
		{
			unsigned int cn = tc[3 * (size_t)s + q];
			pos[q] = ad.element == QAWS_MESH_VERTEX ? sfv[cn] : ad.element == QAWS_MESH_CORNER ? cn : tf[s];
			/* (the hit's weights are on the triangle's vertices as stri lists them: tc follows the same fan) */
			w[q] = (double)hits[i].weights[q];
		}
		hm_blend(vals + i * n, dv, n, pos, w, 3, n, (unsigned int)ad.rule);
		keep[i] = 1;
	}
	st = tx_hits_write(out_coverage, texture, h, nh, keep, vals, desc ? desc->gutter : 0);
done:
	tx_faces_free(&t);
	qaws_mesh_bvh_destroy(bvh);
	qaws_internal_free(dv); qaws_internal_free(vals); qaws_internal_free(h); qaws_internal_free(x); qaws_internal_free(sx); qaws_internal_free(stri);
	qaws_internal_free(sfv); qaws_internal_free(sfs); qaws_internal_free(tf); qaws_internal_free(tc); qaws_internal_free(pts); qaws_internal_free(hits);
	qaws_internal_free(keep);
	return st;
}

qaws_status qaws_texture_sample(qaws_scalar* out_values, qaws_texture const* texture, qaws_scalar const* uv, unsigned int count)
{
	unsigned int i, k;
	if (!out_values || !tx_check(texture) || (count && !uv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < count; i++)
	{
		double y[16];
		tx_read(y, texture, (double)uv[2 * (size_t)i], (double)uv[2 * (size_t)i + 1]);
		for (k = 0; k < texture->components; k++)
			out_values[(size_t)texture->components * i + k] = (qaws_scalar)y[k];
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_sample_texture(qaws_mesh_3d* mesh, qaws_texture const* texture, char const* name, qaws_texture_desc const* desc)
{
	tx_faces t;
	unsigned int id, f, k, n;
	int have;
	qaws_scalar* sv = NULL;
	qaws_attribute_desc ad;
	qaws_status st;
	if (!mesh || !tx_check(texture) || !name)
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = texture->components;
	have = qaws_mesh_3d_find_attribute(&id, mesh, name) == QAWS_STATUS_OK;
	if (have)
	{
		qaws_mesh_3d_get_attribute_desc(&ad, mesh, id);
		if (ad.element != QAWS_MESH_CORNER || ad.components != n)
			return QAWS_STATUS_INVALID_ARGUMENT;
	}
	if ((st = tx_faces_read(&t, mesh, desc)) != QAWS_STATUS_OK)
		return st;
	sv = (qaws_scalar*)qaws_internal_scratch_zalloc(qaws_internal_mul2((size_t)n, (size_t)t.nc + 1), sizeof(qaws_scalar));
	if (!sv)
	{
		tx_faces_free(&t);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	if (have)
		qaws_mesh_3d_get_attribute(sv, mesh, id);
	for (f = 0; f < t.nf; f++)
	{
		unsigned int c;
		if (!t.use[f])
			continue;
		for (c = t.fo[f]; c < t.fo[f] + t.fs[f]; c++)
		{
			double y[16];
			tx_read(y, texture, t.uv[2 * (size_t)c], t.uv[2 * (size_t)c + 1]);
			for (k = 0; k < n; k++)
				sv[(size_t)n * c + k] = (qaws_scalar)y[k];
		}
	}
	if (have)
		st = qaws_mesh_3d_set_attribute(mesh, id, sv);
	else
	{
		memset(&ad, 0, sizeof ad);
		ad.name = name;
		ad.element = QAWS_MESH_CORNER;
		ad.components = n;
		ad.rule = texture->rule;
		st = qaws_mesh_3d_add_attribute(&id, mesh, sv, &ad);
	}
	qaws_internal_free(sv);
	tx_faces_free(&t);
	return st;
}
