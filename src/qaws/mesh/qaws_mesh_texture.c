/*
 * Texture attributes (qaws_mesh_texture.h): triangles rasterized in (u, v)
 * at the texels' centres, values blended by the attribute's rule; reading
 * bilinear by the rule.
 */

#include "../internal/qaws_internal_memory.h"
#include "../qaws_mesh_texture.h"
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

static void tx_faces_free(tx_faces* t)
{
	qaws_internal_free(t->fv); qaws_internal_free(t->fs); qaws_internal_free(t->fo); qaws_internal_free(t->uv); qaws_internal_free(t->use);
}

static qaws_status tx_faces_read(tx_faces* t, qaws_mesh_3d const* m, qaws_texture_desc const* desc)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), id, i, c;
	qaws_attribute_desc ad;
	qaws_scalar* v = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	memset(t, 0, sizeof *t);
	t->nf = qaws_mesh_3d_get_face_count(m);
	t->nc = qaws_mesh_3d_get_corner_count(m);
	if (qaws_mesh_3d_find_attribute(&id, m, desc && desc->uv ? desc->uv : "uv") != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_mesh_3d_get_attribute_desc(&ad, m, id);
	if (ad.components != 2 || (ad.element != QAWS_MESH_CORNER && ad.element != QAWS_MESH_VERTEX))
		return QAWS_STATUS_INVALID_ARGUMENT;
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
	qaws_mesh_3d_get_attribute(v, m, id);
	for (c = 0; c < t->nc; c++)
	{
		unsigned int e = ad.element == QAWS_MESH_CORNER ? c : t->fv[c];
		t->uv[2 * (size_t)c] = (double)v[2 * (size_t)e];
		t->uv[2 * (size_t)c + 1] = (double)v[2 * (size_t)e + 1];
	}
	memset(t->use, 1, t->nf);
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
			t->use[i] = (unsigned char)((double)v[i] == (double)desc->id);
	}
	qaws_internal_free(v);
	return QAWS_STATUS_OK;
fail:
	qaws_internal_free(v);
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

qaws_status qaws_mesh_3d_bake_texture(unsigned char* out_coverage, qaws_texture* texture, qaws_mesh_3d const* mesh, char const* name,
	qaws_texture_desc const* desc)
{
	tx_faces t;
	qaws_attribute_desc ad;
	unsigned int id, f, k, n, W, H, cnt, r, gutter = desc ? desc->gutter : 0;
	size_t nt;
	double* dv = NULL;
	qaws_scalar* sv = NULL;
	unsigned char* ring = NULL;
	unsigned int* grow = NULL;
	double* gval = NULL;
	qaws_status st;
	if (!tx_check(texture) || !mesh || !name || qaws_mesh_3d_find_attribute(&id, mesh, name) != QAWS_STATUS_OK)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_mesh_3d_get_attribute_desc(&ad, mesh, id);
	if (ad.element == QAWS_MESH_EDGE || ad.components != texture->components)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((st = tx_faces_read(&t, mesh, desc)) != QAWS_STATUS_OK)
		return st;
	n = ad.components;
	W = texture->width;
	H = texture->height;
	nt = (size_t)W * H;
	texture->rule = ad.rule;
	cnt = ad.element == QAWS_MESH_VERTEX ? qaws_mesh_3d_get_vertex_count(mesh) : ad.element == QAWS_MESH_CORNER ? t.nc : t.nf;
	st = QAWS_STATUS_ALLOCATION_FAILURE;
	dv = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, (size_t)cnt + 1));
	sv = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), n, (size_t)cnt + 1));
	ring = (unsigned char*)qaws_internal_scratch_zalloc(nt + 1, 1);
	if (!dv || !sv || !ring)
		goto done;
	qaws_mesh_3d_get_attribute(sv, mesh, id);
	for (k = 0; k < n * cnt; k++)
		dv[k] = (double)sv[k];
	/* the triangles (faces fanned from their first corner) at the texels' centres */
	for (f = 0; f < t.nf; f++)
	{
		unsigned int j;
		if (!t.use[f])
			continue;
		for (j = 1; j + 1 < t.fs[f]; j++)
		{
			unsigned int c[3] = { t.fo[f], t.fo[f] + j, t.fo[f] + j + 1 }, pos[3], q;
			double U[3][2], lo[2], hi[2], den;
			long long i0, i1, j0, j1, ii, jj;
			for (q = 0; q < 3; q++)
			{
				U[q][0] = t.uv[2 * (size_t)c[q]];
				U[q][1] = t.uv[2 * (size_t)c[q] + 1];
				/* (across a repeated axis' edge: the copy nearest the first corner) */
				if (q && texture->repeat_u)
					U[q][0] -= floor(U[q][0] - U[0][0] + 0.5);
				if (q && texture->repeat_v)
					U[q][1] -= floor(U[q][1] - U[0][1] + 0.5);
				pos[q] = ad.element == QAWS_MESH_VERTEX ? t.fv[c[q]] : ad.element == QAWS_MESH_CORNER ? c[q] : f;
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
					unsigned int a = tx_index(ii, W, texture->repeat_u), b = tx_index(jj, H, texture->repeat_v);
					double px = (ii + 0.5) / W, py = (jj + 0.5) / H, w[3], y[16];
					size_t tx;
					if (a == ~0u || b == ~0u)
						continue;
					w[1] = ((px - U[0][0]) * (U[2][1] - U[0][1]) - (U[2][0] - U[0][0]) * (py - U[0][1])) / den;
					w[2] = ((U[1][0] - U[0][0]) * (py - U[0][1]) - (px - U[0][0]) * (U[1][1] - U[0][1])) / den;
					w[0] = 1 - w[1] - w[2];
					if (w[0] < -1e-12 || w[1] < -1e-12 || w[2] < -1e-12)
						continue;
					hm_blend(y, dv, n, pos, w, 3, n, (unsigned int)ad.rule);
					tx = (size_t)b * W + a;
					for (q = 0; q < n; q++)
						texture->texels[tx * n + q] = (qaws_scalar)y[q];
					ring[tx] = 1;
				}
		}
	}
	/* the gutter: ring after ring, each new texel the mean of its written neighbours (8 about it) */
	if (gutter)
	{
		grow = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), nt + 1));
		gval = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), n, nt + 1));
		if (!grow || !gval)
			goto done;
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
					unsigned int a2 = tx_index((long long)a + di, W, texture->repeat_u), b2 = tx_index((long long)b + dj, H, texture->repeat_v);
					size_t t2;
					if ((!di && !dj) || a2 == ~0u || b2 == ~0u)
						continue;
					t2 = (size_t)b2 * W + a2;
					if (!ring[t2])
						continue;   /* (this ring's new texels are written after the scan) */
					for (q = 0; q < n; q++)
						vals[(size_t)n * m + q] = (double)texture->texels[t2 * n + q];
					pos[m] = m;
					m++;
				}
			if (!m)
				continue;
			for (q = 0; q < m; q++)
				w[q] = 1.0 / m;
			hm_blend(gval + (size_t)n * ng, vals, n, pos, w, m, n, (unsigned int)ad.rule);
			grow[ng++] = (unsigned int)tx;
		}
		for (tx = 0; tx < ng; tx++)
		{
			for (k = 0; k < n; k++)
				texture->texels[(size_t)grow[tx] * n + k] = (qaws_scalar)gval[tx * n + k];
			ring[grow[tx]] = (unsigned char)(r + 2 > 255 ? 255 : r + 2);
		}
		if (!ng)
			break;
	}
	if (out_coverage)
		for (k = 0; k < nt; k++)
			out_coverage[k] = (unsigned char)(ring[k] > 1 ? 2 : ring[k]);
	st = QAWS_STATUS_OK;
done:
	tx_faces_free(&t);
	qaws_internal_free(dv); qaws_internal_free(sv); qaws_internal_free(ring); qaws_internal_free(grow); qaws_internal_free(gval);
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
