#include "qaws_diff.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include <string.h>

/* Jet slot i holds d^(a+b) S / du^a dv^b with (a, b) below. */
static unsigned char const g_surface_jet_a[QAWS_SURFACE_JET_COUNT] = { 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 };
static unsigned char const g_surface_jet_b[QAWS_SURFACE_JET_COUNT] = { 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 };

/* ================================================================== */
/*  Schema                                                            */
/* ================================================================== */

static qaws_surface_diff_vtable const* surface_diff(qaws_surface const* surface)
{
	return (surface && surface->vtable) ? surface->vtable->diff : NULL;
}

unsigned int qaws_surface_get_diff_capabilities(qaws_surface const* surface)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	return d ? d->capabilities : 0u;
}

qaws_diff_class qaws_surface_get_diff_class(qaws_surface const* surface)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	return d ? d->diff_class : QAWS_DIFF_UNSUPPORTED;
}

qaws_status qaws_surface_describe_fields(
	qaws_surface const* surface,
	qaws_field_desc* out_fields,
	unsigned int capacity,
	unsigned int* out_count)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	unsigned int n;
	if (!surface || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (!d || !d->describe_fields)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	n = d->describe_fields(surface, out_fields, out_fields ? capacity : 0u);
	*out_count = n;
	return (out_fields && n > capacity) ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}

qaws_status qaws_surface_read_field(
	qaws_surface const* surface,
	qaws_diff_field field,
	qaws_scalar* out_values,
	unsigned int capacity,
	unsigned int* out_scalar_count)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	qaws_scalar const* data = NULL;
	unsigned int count = 0, components = 0, total;
	qaws_status s;

	if (!surface || !out_scalar_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_scalar_count = 0;
	if (!d || !d->primal_field)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	s = d->primal_field(surface, field, &data, &count, &components);
	if (s != QAWS_STATUS_OK)
		return s;
	total = count * components;
	*out_scalar_count = total;
	if (!out_values)
		return QAWS_STATUS_OK;
	if (capacity < total)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	memcpy(out_values, data, sizeof(qaws_scalar) * (size_t)total);
	return QAWS_STATUS_OK;
}

/* ================================================================== */
/*  Linear tensor-product engine                                      */
/* ================================================================== */

typedef struct surface_sample
{
	qaws_surface_support support;
	qaws_scalar const* primal;
	qaws_scalar const* weights;   /* rational weights, NULL when polynomial */
} surface_sample;

static void clamp_uv(qaws_surface const* surface, qaws_scalar* u, qaws_scalar* v)
{
	if (*u < surface->u_range.min_value) *u = surface->u_range.min_value;
	if (*u > surface->u_range.max_value) *u = surface->u_range.max_value;
	if (*v < surface->v_range.min_value) *v = surface->v_range.min_value;
	if (*v > surface->v_range.max_value) *v = surface->v_range.max_value;
}

/* Highest single-direction derivative order needed by the channels. */
static unsigned int channel_order(unsigned int channels)
{
	unsigned int i, order = 0;
	for (i = 0; i < QAWS_SURFACE_JET_COUNT; i++)
	{
		if (!(channels & (1u << i)))
			continue;
		if (g_surface_jet_a[i] > order) order = g_surface_jet_a[i];
		if (g_surface_jet_b[i] > order) order = g_surface_jet_b[i];
	}
	return order;
}

/* base_order: highest single-direction order of the requested channels;
   extra: additional orders for coordinate tangents (1) or second
   tangents (2). Rational surfaces always need every partial up to order
   three because the quotient recurrence couples them. */
static qaws_status surface_prepare(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int base_order,
	unsigned int extra,
	surface_sample* s)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	unsigned int count = 0, components = 0, last, order = base_order + extra;
	qaws_surface_support const* sp = &s->support;
	qaws_status st;

	if (!d || !d->linear_support || !d->primal_field)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_OUT_OF_RANGE;

	clamp_uv(surface, &u, &v);
	st = d->linear_support(surface, u, v, order, &s->support);
	if (st != QAWS_STATUS_OK)
		return st;
	if (sp->weight_field != QAWS_FIELD_NONE && order < 3 + extra)
	{
		order = 3 + extra;
		st = d->linear_support(surface, u, v, order, &s->support);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	if (!sp->has_weights || sp->order < order || !sp->u_count || !sp->v_count)
		return QAWS_STATUS_INTERNAL_ERROR;

	st = d->primal_field(surface, sp->field, &s->primal, &count, &components);
	if (st != QAWS_STATUS_OK)
		return st;
	last = (sp->u_first + sp->u_count - 1) * sp->u_stride + (sp->v_first + sp->v_count - 1) * sp->v_stride;
	if (components != 3 || last >= count)
		return QAWS_STATUS_INTERNAL_ERROR;

	s->weights = NULL;
	if (sp->weight_field != QAWS_FIELD_NONE)
	{
		unsigned int wcount = 0, wcomp = 0;
		st = d->primal_field(surface, sp->weight_field, &s->weights, &wcount, &wcomp);
		if (st != QAWS_STATUS_OK)
			return st;
		if (wcomp != 1 || last >= wcount)
			return QAWS_STATUS_INTERNAL_ERROR;
	}
	return QAWS_STATUS_OK;
}

static unsigned int support_element(qaws_surface_support const* sp, unsigned int i, unsigned int j)
{
	return (sp->u_first + i) * sp->u_stride + (sp->v_first + j) * sp->v_stride;
}

/* d^(a+b) S / du^a dv^b from the primal control net. */
static void surface_partial(surface_sample const* s, unsigned int a, unsigned int b, qaws_scalar* out)
{
	qaws_surface_support const* sp = &s->support;
	unsigned int i, j;
	out[0] = out[1] = out[2] = QAWS_ZERO;
	for (i = 0; i < sp->u_count; i++)
	{
		qaws_scalar wu = sp->u_weights[a][i];
		if (wu == QAWS_ZERO)
			continue;
		for (j = 0; j < sp->v_count; j++)
		{
			qaws_scalar w = wu * sp->v_weights[b][j];
			qaws_scalar const* p = s->primal + (size_t)support_element(sp, i, j) * 3;
			out[0] += w * p[0];
			out[1] += w * p[1];
			out[2] += w * p[2];
		}
	}
}

/* Same partial applied to the parameter tangent (masked). */
static void surface_partial_tangent(
	surface_sample const* s,
	unsigned int a,
	unsigned int b,
	qaws_field_view const* view,
	qaws_scalar* out)
{
	qaws_surface_support const* sp = &s->support;
	unsigned int i, j;
	out[0] = out[1] = out[2] = QAWS_ZERO;
	if (!view)
		return;
	for (i = 0; i < sp->u_count; i++)
	{
		qaws_scalar wu = sp->u_weights[a][i];
		if (wu == QAWS_ZERO)
			continue;
		for (j = 0; j < sp->v_count; j++)
		{
			qaws_scalar w = wu * sp->v_weights[b][j];
			qaws_scalar p[3];
			qaws_internal_view_read(view, support_element(sp, i, j), 3, p);
			out[0] += w * p[0];
			out[1] += w * p[1];
			out[2] += w * p[2];
		}
	}
}

/* ------------------------------------------------------------------ */
/*  Rational surfaces                                                 */
/*                                                                    */
/*  A_ab = sum U_a V_b w P,  W_ab = sum U_a V_b w,  S = A / W         */
/*  S_kl = (A_kl - sum_{(i,j) != (0,0)} C(k,i) C(l,j) W_ij S_{k-i,l-j})*/
/*         / W_00                                                     */
/* ------------------------------------------------------------------ */

static qaws_scalar const g_surface_binom[4][4] = {
	{ 1, 0, 0, 0 },
	{ 1, 1, 0, 0 },
	{ 1, 2, 1, 0 },
	{ 1, 3, 3, 1 }
};

static unsigned int jet_index(unsigned int a, unsigned int b)
{
	unsigned int t = a + b;
	return t * (t + 1) / 2 + b;
}

typedef struct surface_rational_sums
{
	qaws_vec3 A, Ad, Add;
	qaws_scalar W, Wd;
} surface_rational_sums;

/* Homogeneous sums of partial (a, b): primal, parameter tangent and the
   parameter-only second derivative 2 w' P'. */
static void surface_hom(
	surface_sample const* s,
	unsigned int a,
	unsigned int b,
	qaws_field_view const* pv,
	qaws_field_view const* wv,
	surface_rational_sums* out)
{
	qaws_surface_support const* sp = &s->support;
	unsigned int i, j;
	memset(out, 0, sizeof(*out));
	for (i = 0; i < sp->u_count; i++)
	{
		qaws_scalar wu = sp->u_weights[a][i];
		if (wu == QAWS_ZERO)
			continue;
		for (j = 0; j < sp->v_count; j++)
		{
			unsigned int e = support_element(sp, i, j);
			qaws_scalar n = wu * sp->v_weights[b][j];
			qaws_scalar w = s->weights[e];
			qaws_scalar const* p = s->primal + (size_t)e * 3;
			qaws_scalar pd[3] = { 0, 0, 0 }, wd = 0;
			if (pv) qaws_internal_view_read(pv, e, 3, pd);
			if (wv) qaws_internal_view_read(wv, e, 1, &wd);
			out->W += n * w;
			out->Wd += n * wd;
			out->A = qaws_v3_add(out->A, qaws_v3(n * w * p[0], n * w * p[1], n * w * p[2]));
			out->Ad = qaws_v3_add(out->Ad, qaws_v3(n * (wd * p[0] + w * pd[0]),
			                                      n * (wd * p[1] + w * pd[1]),
			                                      n * (wd * p[2] + w * pd[2])));
			out->Add = qaws_v3_add(out->Add, qaws_v3(n * QAWS_LITERAL(2.0) * wd * pd[0],
			                                        n * QAWS_LITERAL(2.0) * wd * pd[1],
			                                        n * QAWS_LITERAL(2.0) * wd * pd[2]));
		}
	}
}

/* Rational tangent along (u', v', parameter tangent). All ten partials
   are formed because the recurrence couples them. */
static qaws_status surface_rational_tangent(
	surface_sample const* s,
	qaws_scalar u_dot,
	qaws_scalar v_dot,
	qaws_diff_views const* views,
	int want_second,
	qaws_dual3* out_S)
{
	qaws_field_view const* pv = views ? qaws_diff_views_find(views, s->support.field) : NULL;
	qaws_field_view const* wv = views ? qaws_diff_views_find(views, s->support.weight_field) : NULL;
	qaws_dual3 A[QAWS_SURFACE_JET_COUNT];
	qaws_dual1 W[QAWS_SURFACE_JET_COUNT];
	unsigned int ch, k, l, i, j;
	int has_coord = (u_dot != QAWS_ZERO) || (v_dot != QAWS_ZERO);

	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		unsigned int a = g_surface_jet_a[ch], b = g_surface_jet_b[ch];
		surface_rational_sums h, hu, hv, huu, huv, hvv;
		surface_hom(s, a, b, pv, wv, &h);
		A[ch].v = h.A;
		A[ch].t = h.Ad;
		A[ch].tt = h.Add;
		W[ch].v = h.W;
		W[ch].t = h.Wd;
		W[ch].tt = QAWS_ZERO;
		if (!has_coord)
			continue;

		surface_hom(s, a + 1, b, pv, wv, &hu);
		surface_hom(s, a, b + 1, pv, wv, &hv);
		A[ch].t = qaws_v3_axpy(qaws_v3_axpy(A[ch].t, hu.A, u_dot), hv.A, v_dot);
		W[ch].t += u_dot * hu.W + v_dot * hv.W;
		if (!want_second)
			continue;

		surface_hom(s, a + 2, b, pv, wv, &huu);
		surface_hom(s, a + 1, b + 1, pv, wv, &huv);
		surface_hom(s, a, b + 2, pv, wv, &hvv);
		A[ch].tt = qaws_v3_axpy(A[ch].tt, hu.Ad, QAWS_LITERAL(2.0) * u_dot);
		A[ch].tt = qaws_v3_axpy(A[ch].tt, hv.Ad, QAWS_LITERAL(2.0) * v_dot);
		A[ch].tt = qaws_v3_axpy(A[ch].tt, huu.A, u_dot * u_dot);
		A[ch].tt = qaws_v3_axpy(A[ch].tt, huv.A, QAWS_LITERAL(2.0) * u_dot * v_dot);
		A[ch].tt = qaws_v3_axpy(A[ch].tt, hvv.A, v_dot * v_dot);
		W[ch].tt = QAWS_LITERAL(2.0) * u_dot * hu.Wd + QAWS_LITERAL(2.0) * v_dot * hv.Wd
		         + u_dot * u_dot * huu.W + QAWS_LITERAL(2.0) * u_dot * v_dot * huv.W + v_dot * v_dot * hvv.W;
	}

	if (W[0].v < QAWS_LITERAL(1e-15) && W[0].v > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;

	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		qaws_dual3 num = A[ch];
		k = g_surface_jet_a[ch];
		l = g_surface_jet_b[ch];
		for (i = 0; i <= k; i++)
			for (j = 0; j <= l; j++)
			{
				qaws_scalar bc;
				if (i == 0 && j == 0)
					continue;
				bc = g_surface_binom[k][i] * g_surface_binom[l][j];
				num = qaws_dual3_sub(num, qaws_dual3_mul_const(
					qaws_dual3_scale(out_S[jet_index(k - i, l - j)], W[jet_index(i, j)]), bc));
			}
		out_S[ch] = qaws_dual3_div(num, W[0]);
	}
	return QAWS_STATUS_OK;
}

/* Reverse of the rational recurrence: control point and weight entries,
   plus coordinate adjoints when requested. */
static qaws_status surface_rational_adjoint(
	surface_sample const* s,
	unsigned int channels,
	qaws_surface_jet const* ybar,
	int want_coordinate,
	qaws_diff_entry* entries,
	unsigned int capacity,
	unsigned int* out_count,
	qaws_scalar* u_adj,
	qaws_scalar* v_adj)
{
	qaws_surface_support const* sp = &s->support;
	surface_rational_sums h[QAWS_SURFACE_JET_COUNT];
	qaws_vec3 S[QAWS_SURFACE_JET_COUNT], Sbar[QAWS_SURFACE_JET_COUNT], Abar[QAWS_SURFACE_JET_COUNT];
	qaws_scalar Wbar[QAWS_SURFACE_JET_COUNT];
	unsigned int ch, k, l, i, j, n = 0;
	qaws_scalar W0;

	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		surface_hom(s, g_surface_jet_a[ch], g_surface_jet_b[ch], NULL, NULL, &h[ch]);
	W0 = h[0].W;
	if (W0 < QAWS_LITERAL(1e-15) && W0 > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;

	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		qaws_vec3 num = h[ch].A;
		k = g_surface_jet_a[ch];
		l = g_surface_jet_b[ch];
		for (i = 0; i <= k; i++)
			for (j = 0; j <= l; j++)
				if (i || j)
					num = qaws_v3_sub(num, qaws_v3_scale(S[jet_index(k - i, l - j)],
						g_surface_binom[k][i] * g_surface_binom[l][j] * h[jet_index(i, j)].W));
		S[ch] = qaws_v3_scale(num, QAWS_ONE / W0);
		Sbar[ch] = (channels & (1u << ch)) ? ybar->d[ch] : qaws_v3_zero();
		Abar[ch] = qaws_v3_zero();
		Wbar[ch] = QAWS_ZERO;
	}

	for (ch = QAWS_SURFACE_JET_COUNT; ch-- > 0;)
	{
		qaws_vec3 g = qaws_v3_scale(Sbar[ch], QAWS_ONE / W0);
		k = g_surface_jet_a[ch];
		l = g_surface_jet_b[ch];
		Abar[ch] = qaws_v3_add(Abar[ch], g);
		Wbar[0] -= qaws_v3_dot(g, S[ch]);
		for (i = 0; i <= k; i++)
			for (j = 0; j <= l; j++)
			{
				qaws_scalar bc;
				unsigned int lower = jet_index(k - i, l - j), wij = jet_index(i, j);
				if (i == 0 && j == 0)
					continue;
				bc = g_surface_binom[k][i] * g_surface_binom[l][j];
				Wbar[wij] -= bc * qaws_v3_dot(g, S[lower]);
				Sbar[lower] = qaws_v3_axpy(Sbar[lower], g, -bc * h[wij].W);
			}
	}

	for (i = 0; i < sp->u_count; i++)
		for (j = 0; j < sp->v_count; j++)
		{
			unsigned int e = support_element(sp, i, j);
			qaws_scalar w = s->weights[e];
			qaws_vec3 p = qaws_v3(s->primal[e * 3 + 0], s->primal[e * 3 + 1], s->primal[e * 3 + 2]);
			qaws_vec3 gp = qaws_v3_zero();
			qaws_scalar gw = QAWS_ZERO;
			if (n + 2 > capacity)
				return QAWS_STATUS_INTERNAL_ERROR;
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				qaws_scalar nab = sp->u_weights[g_surface_jet_a[ch]][i] * sp->v_weights[g_surface_jet_b[ch]][j];
				gp = qaws_v3_axpy(gp, Abar[ch], nab * w);
				gw += nab * (qaws_v3_dot(p, Abar[ch]) + Wbar[ch]);
			}
			entries[n].field = sp->field;
			entries[n].element = e;
			entries[n].g[0] = gp.x; entries[n].g[1] = gp.y; entries[n].g[2] = gp.z;
			n++;
			entries[n].field = sp->weight_field;
			entries[n].element = e;
			entries[n].g[0] = gw; entries[n].g[1] = QAWS_ZERO; entries[n].g[2] = QAWS_ZERO;
			n++;
		}
	*out_count = n;

	if (want_coordinate)
	{
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			surface_rational_sums hu, hv;
			unsigned int a = g_surface_jet_a[ch], b = g_surface_jet_b[ch];
			if (u_adj)
			{
				surface_hom(s, a + 1, b, NULL, NULL, &hu);
				*u_adj += qaws_v3_dot(Abar[ch], hu.A) + Wbar[ch] * hu.W;
			}
			if (v_adj)
			{
				surface_hom(s, a, b + 1, NULL, NULL, &hv);
				*v_adj += qaws_v3_dot(Abar[ch], hv.A) + Wbar[ch] * hv.W;
			}
		}
	}
	return QAWS_STATUS_OK;
}

static void store_vec3(qaws_vec3* dst, qaws_scalar const* src)
{
	dst->x = src[0];
	dst->y = src[1];
	dst->z = src[2];
}

static void report_surface(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	surface_sample const* s,
	int has_coordinate,
	unsigned int index)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	qaws_internal_diff_report_note(ctx,
		has_coordinate ? d->diff_class : QAWS_DIFF_SMOOTH,
		(has_coordinate && s->support.on_boundary) ? QAWS_DIFF_AT_BOUNDARY : QAWS_DIFF_VALID,
		(has_coordinate && d->diff_class == QAWS_DIFF_PIECEWISE_SMOOTH) ? (unsigned int)QAWS_FREEZE_SPAN : 0u,
		index);
}

static qaws_status surface_tangent_sample(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	qaws_scalar u_dot,
	qaws_scalar v_dot,
	unsigned int channels,
	qaws_diff_views const* views,
	unsigned int index,
	qaws_surface_jet* primal,
	qaws_surface_jet* tangent,
	qaws_surface_jet* tangent2)
{
	surface_sample s;
	int has_coord = (u_dot != QAWS_ZERO) || (v_dot != QAWS_ZERO);
	unsigned int extra = has_coord ? (tangent2 ? 2u : 1u) : 0u;
	unsigned int ch, c;
	qaws_field_view const* view;
	qaws_status st = surface_prepare(surface, u, v, channel_order(channels), extra, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	view = views ? qaws_diff_views_find(views, s.support.field) : NULL;

	if (primal) memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2) memset(tangent2, 0, sizeof(*tangent2));

	if (s.weights)
	{
		qaws_dual3 S[QAWS_SURFACE_JET_COUNT];
		st = surface_rational_tangent(&s, u_dot, v_dot, views, tangent2 != NULL, S);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
		{
			if (!(channels & (1u << ch)))
				continue;
			if (primal) primal->d[ch] = S[ch].v;
			tangent->d[ch] = S[ch].t;
			if (tangent2) tangent2->d[ch] = S[ch].tt;
		}
	}
	else
	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		unsigned int a = g_surface_jet_a[ch], b = g_surface_jet_b[ch];
		qaws_scalar p[3], tg[3], q[3];
		if (!(channels & (1u << ch)))
			continue;

		if (primal)
		{
			surface_partial(&s, a, b, p);
			store_vec3(&primal->d[ch], p);
		}

		surface_partial_tangent(&s, a, b, view, tg);
		if (has_coord)
		{
			surface_partial(&s, a + 1, b, q);
			for (c = 0; c < 3; c++) tg[c] += u_dot * q[c];
			surface_partial(&s, a, b + 1, q);
			for (c = 0; c < 3; c++) tg[c] += v_dot * q[c];
		}
		store_vec3(&tangent->d[ch], tg);

		if (tangent2 && has_coord)
		{
			/* 2 u' S'_{a+1,b}[P'] + 2 v' S'_{a,b+1}[P'] + u'^2 S_{a+2,b}
			   + 2 u' v' S_{a+1,b+1} + v'^2 S_{a,b+2} */
			qaws_scalar acc[3] = { 0, 0, 0 };
			surface_partial_tangent(&s, a + 1, b, view, q);
			for (c = 0; c < 3; c++) acc[c] += QAWS_LITERAL(2.0) * u_dot * q[c];
			surface_partial_tangent(&s, a, b + 1, view, q);
			for (c = 0; c < 3; c++) acc[c] += QAWS_LITERAL(2.0) * v_dot * q[c];
			surface_partial(&s, a + 2, b, q);
			for (c = 0; c < 3; c++) acc[c] += u_dot * u_dot * q[c];
			surface_partial(&s, a + 1, b + 1, q);
			for (c = 0; c < 3; c++) acc[c] += QAWS_LITERAL(2.0) * u_dot * v_dot * q[c];
			surface_partial(&s, a, b + 2, q);
			for (c = 0; c < 3; c++) acc[c] += v_dot * v_dot * q[c];
			store_vec3(&tangent2->d[ch], acc);
		}
	}

	if (primal) primal->channels = channels;
	tangent->channels = channels;
	if (tangent2) tangent2->channels = channels;

	report_surface(ctx, surface, &s, has_coord, index);
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Primal jet                                                        */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_eval_jet(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int channels,
	qaws_surface_jet* out_jet)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	surface_sample s;
	unsigned int ch;
	qaws_status st;

	if (!surface || !out_jet)
		return QAWS_STATUS_INVALID_ARGUMENT;
	channels &= (unsigned int)QAWS_SJET_ORDER3;
	if (d && d->eval_jet)
	{
		clamp_uv(surface, &u, &v);
		return d->eval_jet(surface, u, v, channels, out_jet);
	}

	st = surface_prepare(surface, u, v, channel_order(channels), 0, &s);
	if (st != QAWS_STATUS_OK)
		return st;
	memset(out_jet, 0, sizeof(*out_jet));
	if (s.weights)
	{
		qaws_dual3 S[QAWS_SURFACE_JET_COUNT];
		st = surface_rational_tangent(&s, QAWS_ZERO, QAWS_ZERO, NULL, 0, S);
		if (st != QAWS_STATUS_OK)
			return st;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			if (channels & (1u << ch))
				out_jet->d[ch] = S[ch].v;
		out_jet->channels = channels;
		return QAWS_STATUS_OK;
	}
	for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
	{
		qaws_scalar p[3];
		if (!(channels & (1u << ch)))
			continue;
		surface_partial(&s, g_surface_jet_a[ch], g_surface_jet_b[ch], p);
		store_vec3(&out_jet->d[ch], p);
	}
	out_jet->channels = channels;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Tangent                                                           */
/* ------------------------------------------------------------------ */

static qaws_status surface_batch_tangent(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	qaws_scalar const* u_tangent,
	qaws_scalar const* v_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal,
	qaws_surface_jet* out_tangent,
	qaws_surface_jet* out_tangent2)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	unsigned int i;
	qaws_status st;

	if (!surface || !u || !v || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!d)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (out_tangent2 && !(d->capabilities & QAWS_CAP_TANGENT2))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	st = qaws_internal_check_views(param_tangent, 3);
	if (st != QAWS_STATUS_OK)
		return st;
	channels &= (unsigned int)QAWS_SJET_ORDER3;

	for (i = 0; i < count; i++)
	{
		st = surface_tangent_sample(ctx, surface, u[i], v[i],
			u_tangent ? u_tangent[i] : QAWS_ZERO,
			v_tangent ? v_tangent[i] : QAWS_ZERO,
			channels, param_tangent, i,
			out_primal ? &out_primal[i] : NULL,
			&out_tangent[i],
			out_tangent2 ? &out_tangent2[i] : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_eval_batch_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar const* u, qaws_scalar const* v,
	qaws_scalar const* u_tangent, qaws_scalar const* v_tangent,
	unsigned int count, unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal, qaws_surface_jet* out_tangent)
{
	return surface_batch_tangent(ctx, surface, u, v, u_tangent, v_tangent, count,
		channels, param_tangent, out_primal, out_tangent, NULL);
}

qaws_status qaws_surface_eval_batch_tangent2(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar const* u, qaws_scalar const* v,
	qaws_scalar const* u_tangent, qaws_scalar const* v_tangent,
	unsigned int count, unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal, qaws_surface_jet* out_tangent,
	qaws_surface_jet* out_tangent2)
{
	if (!out_tangent2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return surface_batch_tangent(ctx, surface, u, v, u_tangent, v_tangent, count,
		channels, param_tangent, out_primal, out_tangent, out_tangent2);
}

qaws_status qaws_surface_eval_tangent(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, qaws_scalar u_tangent, qaws_scalar v_tangent,
	unsigned int channels, qaws_diff_views const* param_tangent,
	qaws_surface_jet* out_primal, qaws_surface_jet* out_tangent)
{
	return surface_batch_tangent(ctx, surface, &u, &v, &u_tangent, &v_tangent, 1,
		channels, param_tangent, out_primal, out_tangent, NULL);
}

/* ------------------------------------------------------------------ */
/*  Adjoint                                                           */
/* ------------------------------------------------------------------ */

typedef struct surface_adjoint_job
{
	qaws_diff_context const* ctx;
	qaws_surface const* surface;
	qaws_scalar const* u;
	qaws_scalar const* v;
	unsigned int channels;
	unsigned int order;
	qaws_surface_jet const* jets;
	qaws_scalar* u_adjoint;
	qaws_scalar* v_adjoint;
} surface_adjoint_job;

static qaws_status surface_collect(
	void const* user,
	unsigned int sample,
	int first_pass,
	qaws_diff_entry* entries,
	unsigned int capacity,
	unsigned int* out_count)
{
	surface_adjoint_job const* job = (surface_adjoint_job const*)user;
	qaws_surface_jet const* ybar = &job->jets[sample];
	surface_sample s;
	qaws_surface_support const* sp;
	unsigned int i, j, ch, n = 0;
	qaws_status st = surface_prepare(job->surface, job->u[sample], job->v[sample], job->order, (job->u_adjoint || job->v_adjoint) ? 1u : 0u, &s);

	*out_count = 0;
	if (st != QAWS_STATUS_OK)
		return st;
	sp = &s.support;

	if (s.weights)
	{
		int want = first_pass && (job->u_adjoint || job->v_adjoint);
		st = surface_rational_adjoint(&s, job->channels, ybar, want, entries, capacity, out_count,
			(want && job->u_adjoint) ? &job->u_adjoint[sample] : NULL,
			(want && job->v_adjoint) ? &job->v_adjoint[sample] : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		if (first_pass)
			report_surface(job->ctx, job->surface, &s, job->u_adjoint || job->v_adjoint, sample);
		return QAWS_STATUS_OK;
	}

	if (sp->u_count * sp->v_count > capacity)
		return QAWS_STATUS_INTERNAL_ERROR;

	for (i = 0; i < sp->u_count; i++)
	{
		for (j = 0; j < sp->v_count; j++)
		{
			qaws_diff_entry* e = &entries[n++];
			e->field = sp->field;
			e->element = support_element(sp, i, j);
			e->g[0] = e->g[1] = e->g[2] = QAWS_ZERO;
			for (ch = 0; ch < QAWS_SURFACE_JET_COUNT; ch++)
			{
				qaws_scalar w;
				if (!(job->channels & (1u << ch)))
					continue;
				w = sp->u_weights[g_surface_jet_a[ch]][i] * sp->v_weights[g_surface_jet_b[ch]][j];
				e->g[0] += w * ybar->d[ch].x;
				e->g[1] += w * ybar->d[ch].y;
				e->g[2] += w * ybar->d[ch].z;
			}
		}
	}
	*out_count = n;

	if (first_pass)
	{
		int has_coord = job->u_adjoint || job->v_adjoint;
		for (ch = 0; ch < QAWS_SURFACE_JET_COUNT && has_coord; ch++)
		{
			unsigned int a = g_surface_jet_a[ch], b = g_surface_jet_b[ch];
			qaws_scalar q[3];
			if (!(job->channels & (1u << ch)))
				continue;
			if (job->u_adjoint)
			{
				surface_partial(&s, a + 1, b, q);
				job->u_adjoint[sample] += ybar->d[ch].x * q[0] + ybar->d[ch].y * q[1] + ybar->d[ch].z * q[2];
			}
			if (job->v_adjoint)
			{
				surface_partial(&s, a, b + 1, q);
				job->v_adjoint[sample] += ybar->d[ch].x * q[0] + ybar->d[ch].y * q[1] + ybar->d[ch].z * q[2];
			}
		}
		report_surface(job->ctx, job->surface, &s, has_coord, sample);
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_surface_eval_batch_adjoint(
	qaws_diff_context const* ctx,
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	unsigned int count,
	unsigned int channels,
	qaws_surface_jet const* out_adjoint,
	qaws_diff_views* param_adjoint,
	qaws_scalar* u_adjoint,
	qaws_scalar* v_adjoint)
{
	surface_adjoint_job job;
	qaws_status st;

	if (!surface || !u || !v || !out_adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!surface_diff(surface))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	st = qaws_internal_check_views(param_adjoint, 3);
	if (st != QAWS_STATUS_OK)
		return st;

	job.ctx = ctx;
	job.surface = surface;
	job.u = u;
	job.v = v;
	job.channels = channels & (unsigned int)QAWS_SJET_ORDER3;
	job.order = channel_order(job.channels);
	job.jets = out_adjoint;
	job.u_adjoint = u_adjoint;
	job.v_adjoint = v_adjoint;

	return qaws_internal_diff_accumulate(ctx, 3, count,
		2 * QAWS_DIFF_MAX_SUPPORT * QAWS_DIFF_MAX_SUPPORT, surface_collect, &job, param_adjoint);
}

qaws_status qaws_surface_eval_adjoint(
	qaws_diff_context const* ctx, qaws_surface const* surface,
	qaws_scalar u, qaws_scalar v, unsigned int channels,
	qaws_surface_jet const* out_adjoint, qaws_diff_views* param_adjoint,
	qaws_scalar* u_adjoint, qaws_scalar* v_adjoint)
{
	return qaws_surface_eval_batch_adjoint(ctx, surface, &u, &v, 1, channels,
		out_adjoint, param_adjoint, u_adjoint, v_adjoint);
}

/* ------------------------------------------------------------------ */
/*  Local support and support index                                   */
/* ------------------------------------------------------------------ */

qaws_status qaws_surface_local_support(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int order,
	qaws_surface_support* out_support)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	if (!surface || !out_support)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_OUT_OF_RANGE;
	if (!d || !d->linear_support)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	clamp_uv(surface, &u, &v);
	return d->linear_support(surface, u, v, order, out_support);
}

qaws_status qaws_surface_build_support_index(
	qaws_surface const* surface,
	qaws_scalar const* u,
	qaws_scalar const* v,
	unsigned int count,
	qaws_diff_field field,
	unsigned int* out_offsets,
	unsigned int offset_capacity,
	unsigned int* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_entry_count)
{
	qaws_surface_diff_vtable const* d = surface_diff(surface);
	qaws_scalar const* data = NULL;
	unsigned int element_count = 0, components = 0, s, i, j, e, total, pass;
	unsigned int* cursor = NULL;
	qaws_status st;

	if (!surface || !u || !v || !out_offsets || !out_entry_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_entry_count = 0;
	if (!d || !d->linear_support || !d->primal_field)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	st = d->primal_field(surface, field, &data, &element_count, &components);
	if (st != QAWS_STATUS_OK)
		return st;
	(void)data;
	if (offset_capacity < element_count + 1)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	memset(out_offsets, 0, sizeof(unsigned int) * ((size_t)element_count + 1));

	/* pass 0 counts, pass 1 fills */
	for (pass = 0; pass < 2; pass++)
	{
		for (s = 0; s < count; s++)
		{
			qaws_surface_support sp;
			st = qaws_surface_local_support(surface, u[s], v[s], 0, &sp);
			if (st != QAWS_STATUS_OK)
				break;
			if (sp.field != field)
				continue;
			for (i = 0; i < sp.u_count; i++)
				for (j = 0; j < sp.v_count; j++)
				{
					e = support_element(&sp, i, j);
					if (e >= element_count)
						continue;
					if (pass == 0)
						out_offsets[e + 1]++;
					else
						out_samples[cursor[e]++] = s;
				}
		}
		if (st != QAWS_STATUS_OK)
			break;

		if (pass == 0)
		{
			for (e = 0; e < element_count; e++)
				out_offsets[e + 1] += out_offsets[e];
			total = out_offsets[element_count];
			*out_entry_count = total;
			if (!out_samples || sample_capacity < total)
				return QAWS_STATUS_BUFFER_TOO_SMALL;
			cursor = (unsigned int*)qaws_internal_alloc(surface->allocator,
				(unsigned long)(sizeof(unsigned int) * ((size_t)element_count + 1)));
			if (!cursor)
				return QAWS_STATUS_ALLOCATION_FAILURE;
			memcpy(cursor, out_offsets, sizeof(unsigned int) * ((size_t)element_count + 1));
		}
	}

	if (cursor)
		qaws_internal_dealloc(surface->allocator, cursor);
	return st;
}
