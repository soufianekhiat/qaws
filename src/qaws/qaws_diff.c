#include "qaws_diff.h"
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_span.h"
#include "internal/qaws_internal_diff.h"
#include "core/qaws_dual_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIFF_DEFAULT_TILE 256u

/* ================================================================== */
/*  Context and report                                                */
/* ================================================================== */

void qaws_diff_context_init(qaws_diff_context* ctx)
{
	if (!ctx)
		return;
	memset(ctx, 0, sizeof(*ctx));
	ctx->order = 1;
	ctx->frozen = QAWS_FREEZE_NONE;
	ctx->accumulation = QAWS_ACCUMULATE_SCATTER;
}

void qaws_diff_report_reset(qaws_diff_report* report)
{
	if (!report)
		return;
	memset(report, 0, sizeof(*report));
	report->diff_class = QAWS_DIFF_SMOOTH;
	report->validity = QAWS_DIFF_VALID;
	report->condition_number = QAWS_ONE;
	report->branch_gap = QAWS_DIFF_UNBOUNDED;
}

void qaws_internal_diff_report_note(
	qaws_diff_context const* ctx,
	qaws_diff_class diff_class,
	qaws_diff_validity validity,
	unsigned int frozen,
	unsigned int index)
{
	qaws_diff_report* r;
	if (!ctx || !ctx->report)
		return;
	r = ctx->report;
	if (diff_class > r->diff_class)
		r->diff_class = diff_class;
	if (validity > r->validity)
	{
		r->validity = validity;
		r->worst_index = index;
	}
	r->frozen_used |= frozen;
	r->evaluation_count++;
}

/* ================================================================== */
/*  Views                                                             */
/* ================================================================== */

qaws_field_view qaws_field_view_make(
	qaws_diff_field field,
	qaws_scalar* data,
	unsigned int count,
	unsigned int components)
{
	qaws_field_view v;
	memset(&v, 0, sizeof(v));
	v.field = field;
	v.data = data;
	v.count = count;
	v.components = components;
	v.stride = components;
	return v;
}

qaws_field_view* qaws_diff_views_find(
	qaws_diff_views const* views,
	qaws_diff_field field)
{
	unsigned int i;
	if (!views)
		return NULL;
	for (i = 0; i < views->field_count; i++)
		if (views->fields[i].field == field && views->fields[i].data)
			return &views->fields[i];
	return NULL;
}

unsigned int qaws_internal_view_stride(qaws_field_view const* v)
{
	return v->stride ? v->stride : v->components;
}

int qaws_internal_view_element_active(qaws_field_view const* v, unsigned int e)
{
	return e < v->count && (!v->active || v->active[e]);
}

int qaws_internal_view_component_active(qaws_field_view const* v, unsigned int c)
{
	return !v->component_mask || ((v->component_mask >> c) & 1u);
}

void qaws_diff_views_clear(qaws_diff_views* views)
{
	unsigned int i, e, c;
	if (!views)
		return;
	for (i = 0; i < views->field_count; i++)
	{
		qaws_field_view* v = &views->fields[i];
		unsigned int stride = qaws_internal_view_stride(v);
		if (!v->data)
			continue;
		for (e = 0; e < v->count; e++)
		{
			if (!qaws_internal_view_element_active(v, e))
				continue;
			for (c = 0; c < v->components; c++)
				if (qaws_internal_view_component_active(v, c))
					v->data[e * stride + c] = QAWS_ZERO;
		}
	}
	for (i = 0; i < views->child_count; i++)
		qaws_diff_views_clear(&views->children[i]);
}

/* ================================================================== */
/*  Parameter keys                                                    */
/* ================================================================== */

static char const* const g_field_names[QAWS_FIELD_COUNT] = {
	"none",
	"control_points",
	"weights",
	"knots",
	"u_knots",
	"v_knots",
	"points",
	"derivatives",
	"coefficients",
	"key_times",
	"center",
	"radius",
	"radius_b",
	"angle_start",
	"angle_end",
	"axis_u",
	"axis_v",
	"offset_distance",
	"scale",
	"curvature",
	"curvature_rate",
	"direction",
	"parameter"
};

char const* qaws_diff_field_name(qaws_diff_field field)
{
	if ((unsigned int)field >= (unsigned int)QAWS_FIELD_COUNT)
		return "unknown";
	return g_field_names[field];
}

qaws_diff_field qaws_diff_field_from_name(char const* name)
{
	unsigned int i;
	if (!name)
		return QAWS_FIELD_NONE;
	for (i = 1; i < (unsigned int)QAWS_FIELD_COUNT; i++)
		if (strcmp(name, g_field_names[i]) == 0)
			return (qaws_diff_field)i;
	return QAWS_FIELD_NONE;
}

qaws_param_key qaws_param_key_make(
	qaws_diff_field field,
	unsigned int element,
	unsigned int component,
	unsigned int components)
{
	qaws_param_key k;
	memset(&k, 0, sizeof(k));
	k.field = (unsigned char)field;
	k.element = element;
	k.component = (unsigned char)component;
	k.components = (unsigned char)components;
	return k;
}

qaws_status qaws_param_key_prepend_child(qaws_param_key* key, unsigned int child_index)
{
	int i;
	if (!key || child_index > 0xFFFFu)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (key->depth >= QAWS_PARAM_KEY_MAX_DEPTH)
		return QAWS_STATUS_OUT_OF_RANGE;
	for (i = (int)key->depth; i > 0; i--)
		key->child[i] = key->child[i - 1];
	key->child[0] = (unsigned short)child_index;
	key->depth++;
	return QAWS_STATUS_OK;
}

int qaws_param_key_compare(qaws_param_key const* a, qaws_param_key const* b)
{
	unsigned int i;
	if (a->depth != b->depth)
		return a->depth < b->depth ? -1 : 1;
	for (i = 0; i < a->depth; i++)
		if (a->child[i] != b->child[i])
			return a->child[i] < b->child[i] ? -1 : 1;
	if (a->field != b->field)
		return a->field < b->field ? -1 : 1;
	if (a->element != b->element)
		return a->element < b->element ? -1 : 1;
	if (a->component != b->component)
		return a->component < b->component ? -1 : 1;
	return 0;
}

unsigned int qaws_param_key_to_string(
	qaws_param_key const* key,
	char* buffer,
	unsigned int capacity)
{
	char tmp[160];
	int n = 0;
	unsigned int i;
	static char const comp_names[3] = { 'x', 'y', 'z' };

	if (!key)
		return 0;
	tmp[0] = 0;
	for (i = 0; i < key->depth && i < QAWS_PARAM_KEY_MAX_DEPTH; i++)
		n += snprintf(tmp + n, sizeof(tmp) - (size_t)n, "child[%u]/", (unsigned int)key->child[i]);
	n += snprintf(tmp + n, sizeof(tmp) - (size_t)n, "%s/%u",
		qaws_diff_field_name((qaws_diff_field)key->field), key->element);
	if (key->components != 1)
	{
		if (key->component < 3)
			n += snprintf(tmp + n, sizeof(tmp) - (size_t)n, "/%c", comp_names[key->component]);
		else
			n += snprintf(tmp + n, sizeof(tmp) - (size_t)n, "/%u", (unsigned int)key->component);
	}

	if (buffer && capacity > 0)
	{
		unsigned int copy = (unsigned int)n < capacity - 1 ? (unsigned int)n : capacity - 1;
		memcpy(buffer, tmp, copy);
		buffer[copy] = 0;
	}
	return (unsigned int)n;
}

qaws_status qaws_param_key_parse(char const* text, qaws_param_key* out_key)
{
	char token[64];
	char const* p = text;
	qaws_param_key k;
	int stage = 0; /* 0 = children/field, 1 = element, 2 = component, 3 = done */

	if (!text || !out_key)
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(&k, 0, sizeof(k));
	k.components = 1;

	while (*p)
	{
		size_t len = 0;
		while (p[len] && p[len] != '/')
			len++;
		if (len == 0 || len >= sizeof(token))
			return QAWS_STATUS_INVALID_ARGUMENT;
		memcpy(token, p, len);
		token[len] = 0;
		p += len;
		if (*p == '/')
			p++;

		if (stage == 0)
		{
			unsigned int idx;
			if (sscanf(token, "child[%u]", &idx) == 1)
			{
				if (k.depth >= QAWS_PARAM_KEY_MAX_DEPTH || idx > 0xFFFFu)
					return QAWS_STATUS_OUT_OF_RANGE;
				k.child[k.depth++] = (unsigned short)idx;
				continue;
			}
			k.field = (unsigned char)qaws_diff_field_from_name(token);
			if (k.field == QAWS_FIELD_NONE)
				return QAWS_STATUS_INVALID_ARGUMENT;
			stage = 1;
		}
		else if (stage == 1)
		{
			char* end = NULL;
			unsigned long v = strtoul(token, &end, 10);
			if (!end || *end)
				return QAWS_STATUS_INVALID_ARGUMENT;
			k.element = (unsigned int)v;
			stage = 2;
		}
		else if (stage == 2)
		{
			if (token[1] != 0)
				return QAWS_STATUS_INVALID_ARGUMENT;
			if (token[0] >= 'x' && token[0] <= 'z')
				k.component = (unsigned char)(token[0] - 'x');
			else if (token[0] >= '0' && token[0] <= '2')
				k.component = (unsigned char)(token[0] - '0');
			else
				return QAWS_STATUS_INVALID_ARGUMENT;
			k.components = 3;
			stage = 3;
		}
		else
			return QAWS_STATUS_INVALID_ARGUMENT;
	}

	if (stage < 2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_key = k;
	return QAWS_STATUS_OK;
}

/* ================================================================== */
/*  Curve schema                                                      */
/* ================================================================== */

static qaws_curve_diff_vtable const* curve_diff(qaws_curve const* curve)
{
	return (curve && curve->vtable) ? curve->vtable->diff : NULL;
}

unsigned int qaws_curve_get_diff_capabilities(qaws_curve const* curve)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	return d ? d->capabilities : 0u;
}

qaws_diff_class qaws_curve_get_diff_class(qaws_curve const* curve)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	return d ? d->diff_class : QAWS_DIFF_UNSUPPORTED;
}

qaws_coordinate_kind qaws_curve_get_coordinate_kind(qaws_curve const* curve)
{
	if (curve && curve->kind == QAWS_CURVE_KIND_REPARAMETERIZED)
		return QAWS_COORDINATE_ARC_LENGTH;
	return QAWS_COORDINATE_PARAMETRIC;
}

qaws_status qaws_curve_describe_fields(
	qaws_curve const* curve,
	qaws_field_desc* out_fields,
	unsigned int capacity,
	unsigned int* out_count)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	unsigned int n;
	if (!curve || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (!d || !d->describe_fields)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	n = d->describe_fields(curve, out_fields, out_fields ? capacity : 0u);
	*out_count = n;
	return (out_fields && n > capacity) ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}

qaws_status qaws_curve_read_field(
	qaws_curve const* curve,
	qaws_diff_field field,
	qaws_scalar* out_values,
	unsigned int capacity,
	unsigned int* out_scalar_count)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	qaws_scalar const* data = NULL;
	unsigned int count = 0, components = 0, total;
	qaws_status s;

	if (!curve || !out_scalar_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_scalar_count = 0;
	if (!d || !d->primal_field)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	s = d->primal_field(curve, field, &data, &count, &components);
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
/*  Linear curve engine                                               */
/* ================================================================== */

typedef struct curve_sample
{
	qaws_local_support support;
	qaws_scalar const* primal[QAWS_DIFF_MAX_RANGES];
	qaws_scalar const* weights;   /* rational weights (range 0), NULL when polynomial */
	int at_boundary;
} curve_sample;

static unsigned int highest_channel(unsigned int channels)
{
	unsigned int k, kmax = 0;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
		if (channels & (1u << k))
			kmax = k;
	return kmax;
}

static qaws_status curve_prepare(
	qaws_curve const* curve,
	unsigned int dim,
	qaws_scalar t,
	unsigned int order,
	curve_sample* s)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	qaws_scalar local_t;
	unsigned int span, r;
	qaws_status st;
	qaws_scalar const eps = QAWS_LITERAL(16.0) * QAWS_EPSILON;

	if (!d || !d->linear_support || !d->primal_field)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;

	span = qaws_internal_find_span(curve, t, &local_t);
	st = d->linear_support(curve, span, local_t, order, &s->support);
	if (st != QAWS_STATUS_OK)
		return st;
	if (!s->support.has_weights || s->support.order < order)
		return QAWS_STATUS_INTERNAL_ERROR;

	for (r = 0; r < s->support.range_count; r++)
	{
		unsigned int count = 0, components = 0;
		st = d->primal_field(curve, s->support.ranges[r].field, &s->primal[r], &count, &components);
		if (st != QAWS_STATUS_OK)
			return st;
		if (components != dim ||
		    s->support.ranges[r].first + s->support.ranges[r].count > count)
			return QAWS_STATUS_INTERNAL_ERROR;
	}

	s->weights = NULL;
	if (s->support.weight_field != QAWS_FIELD_NONE)
	{
		unsigned int count = 0, components = 0;
		st = d->primal_field(curve, s->support.weight_field, &s->weights, &count, &components);
		if (st != QAWS_STATUS_OK)
			return st;
		if (components != 1 || s->support.range_count != 1 ||
		    s->support.ranges[0].first + s->support.ranges[0].count > count)
			return QAWS_STATUS_INTERNAL_ERROR;
	}

	s->at_boundary = (span > 0 && local_t <= eps) ||
	                 (span + 1 < curve->span_count && local_t >= QAWS_ONE - eps);
	return QAWS_STATUS_OK;
}

/* sum_j w[r][k][j] * primal_r[first + j] over all ranges */
static void sample_primal(curve_sample const* s, unsigned int dim, unsigned int k, qaws_scalar* out)
{
	unsigned int r, j, c;
	for (c = 0; c < dim; c++)
		out[c] = QAWS_ZERO;
	for (r = 0; r < s->support.range_count; r++)
	{
		qaws_support_range const* rg = &s->support.ranges[r];
		for (j = 0; j < rg->count; j++)
		{
			qaws_scalar w = s->support.weights[r][k][j];
			qaws_scalar const* p = s->primal[r] + (size_t)(rg->first + j) * dim;
			for (c = 0; c < dim; c++)
				out[c] += w * p[c];
		}
	}
}

/* sum_j w[r][k][j] * tangent_r[first + j] over all ranges (masked) */
static qaws_status sample_param_tangent(
	curve_sample const* s,
	unsigned int dim,
	unsigned int k,
	qaws_diff_views const* views,
	qaws_scalar* out)
{
	unsigned int r, j, c;
	for (c = 0; c < dim; c++)
		out[c] = QAWS_ZERO;
	if (!views)
		return QAWS_STATUS_OK;
	for (r = 0; r < s->support.range_count; r++)
	{
		qaws_support_range const* rg = &s->support.ranges[r];
		qaws_field_view const* v = qaws_diff_views_find(views, rg->field);
		unsigned int stride;
		if (!v)
			continue;
		if (v->components != dim)
			return QAWS_STATUS_INVALID_ARGUMENT;
		stride = qaws_internal_view_stride(v);
		for (j = 0; j < rg->count; j++)
		{
			unsigned int e = rg->first + j;
			qaws_scalar w = s->support.weights[r][k][j];
			if (!qaws_internal_view_element_active(v, e))
				continue;
			for (c = 0; c < dim; c++)
				if (qaws_internal_view_component_active(v, c))
					out[c] += w * v->data[(size_t)e * stride + c];
		}
	}
	return QAWS_STATUS_OK;
}

typedef struct curve_jet_buf
{
	qaws_scalar d[QAWS_CURVE_JET_ORDER + 1][3];
} curve_jet_buf;

static qaws_status rational_tangent_sample(
	curve_sample const* s, unsigned int dim, qaws_scalar t_dot, unsigned int channels,
	qaws_diff_views const* views, int want_second,
	curve_jet_buf* primal, curve_jet_buf* tangent, curve_jet_buf* tangent2);

static qaws_status curve_tangent_sample(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	unsigned int dim,
	qaws_scalar t,
	qaws_scalar t_dot,
	unsigned int channels,
	qaws_diff_views const* views,
	int want_second,
	unsigned int index,
	curve_jet_buf* primal,
	curve_jet_buf* tangent,
	curve_jet_buf* tangent2)
{
	curve_sample s;
	unsigned int kmax = highest_channel(channels);
	int has_t = (t_dot != QAWS_ZERO);
	unsigned int order = kmax + (has_t ? (want_second ? 2u : 1u) : 0u);
	unsigned int k, c;
	qaws_status st;
	qaws_curve_diff_vtable const* d = curve_diff(curve);

	st = curve_prepare(curve, dim, t, order, &s);
	if (st != QAWS_STATUS_OK)
		return st;

	memset(primal, 0, sizeof(*primal));
	memset(tangent, 0, sizeof(*tangent));
	if (tangent2)
		memset(tangent2, 0, sizeof(*tangent2));

	if (s.weights)
	{
		st = rational_tangent_sample(&s, dim, t_dot, channels, views, want_second, primal, tangent, tangent2);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	else
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		qaws_scalar dp[3], p1[3], p2[3];
		if (!(channels & (1u << k)))
			continue;

		sample_primal(&s, dim, k, primal->d[k]);

		st = sample_param_tangent(&s, dim, k, views, tangent->d[k]);
		if (st != QAWS_STATUS_OK)
			return st;
		if (has_t)
		{
			sample_primal(&s, dim, k + 1, p1);
			for (c = 0; c < dim; c++)
				tangent->d[k][c] += t_dot * p1[c];
		}

		if (tangent2 && has_t)
		{
			/* d2/de2 C_k(t + e t', P + e P') = 2 t' C'_{k+1}[P'] + t'^2 C_{k+2} */
			st = sample_param_tangent(&s, dim, k + 1, views, dp);
			if (st != QAWS_STATUS_OK)
				return st;
			sample_primal(&s, dim, k + 2, p2);
			for (c = 0; c < dim; c++)
				tangent2->d[k][c] = QAWS_LITERAL(2.0) * t_dot * dp[c] + t_dot * t_dot * p2[c];
		}
	}

	qaws_internal_diff_report_note(ctx,
		has_t ? d->diff_class : QAWS_DIFF_SMOOTH,
		(has_t && s.at_boundary) ? QAWS_DIFF_AT_BOUNDARY : QAWS_DIFF_VALID,
		(has_t && d->diff_class == QAWS_DIFF_PIECEWISE_SMOOTH) ? (unsigned int)QAWS_FREEZE_SPAN : 0u,
		index);
	return QAWS_STATUS_OK;
}

/* Contribution g_j = sum_k w[r][k][j] * ybar_k for one range element. */
static void adjoint_contribution(
	curve_sample const* s,
	unsigned int r,
	unsigned int j,
	unsigned int dim,
	unsigned int channels,
	curve_jet_buf const* ybar,
	qaws_scalar* g)
{
	unsigned int k, c;
	for (c = 0; c < dim; c++)
		g[c] = QAWS_ZERO;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		qaws_scalar w;
		if (!(channels & (1u << k)))
			continue;
		w = s->support.weights[r][k][j];
		for (c = 0; c < dim; c++)
			g[c] += w * ybar->d[k][c];
	}
}

static qaws_scalar coordinate_adjoint(
	curve_sample const* s,
	unsigned int dim,
	unsigned int channels,
	curve_jet_buf const* ybar)
{
	unsigned int k, c;
	qaws_scalar acc = QAWS_ZERO;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		qaws_scalar p1[3];
		if (!(channels & (1u << k)))
			continue;
		sample_primal(s, dim, k + 1, p1);
		for (c = 0; c < dim; c++)
			acc += ybar->d[k][c] * p1[c];
	}
	return acc;
}

/* ------------------------------------------------------------------ */
/*  Rational curves: homogeneous sums and the quotient recurrence     */
/*                                                                    */
/*  A^(k) = sum N^(k) w P,  W^(k) = sum N^(k) w,  C = A / W           */
/*  C^(k) = (A^(k) - sum_{i=1..k} binom(k,i) W^(i) C^(k-i)) / W^(0)   */
/* ------------------------------------------------------------------ */

static qaws_scalar const g_binom[4][4] = {
	{ 1, 0, 0, 0 },
	{ 1, 1, 0, 0 },
	{ 1, 2, 1, 0 },
	{ 1, 3, 3, 1 }
};

/* Homogeneous sums of one basis row: primal, parameter tangent, and the
   parameter-only second derivative (2 w' P'). */
typedef struct rational_row
{
	qaws_scalar A[3], W;
	qaws_scalar Ad[3], Wd;
	qaws_scalar Add[3];
} rational_row;

static void rational_sums(
	curve_sample const* s,
	unsigned int dim,
	unsigned int r,
	qaws_field_view const* pv,
	qaws_field_view const* wv,
	rational_row* out)
{
	qaws_support_range const* rg = &s->support.ranges[0];
	unsigned int j, c;
	memset(out, 0, sizeof(*out));
	for (j = 0; j < rg->count; j++)
	{
		unsigned int e = rg->first + j;
		qaws_scalar n = s->support.weights[0][r][j];
		qaws_scalar w = s->weights[e];
		qaws_scalar const* p = s->primal[0] + (size_t)e * dim;
		qaws_scalar pd[3] = { 0, 0, 0 }, wd = 0;
		if (n == QAWS_ZERO)
			continue;
		if (pv) qaws_internal_view_read(pv, e, dim, pd);
		if (wv) qaws_internal_view_read(wv, e, 1, &wd);
		out->W += n * w;
		out->Wd += n * wd;
		for (c = 0; c < dim; c++)
		{
			out->A[c] += n * w * p[c];
			out->Ad[c] += n * (wd * p[c] + w * pd[c]);
			out->Add[c] += n * QAWS_LITERAL(2.0) * wd * pd[c];
		}
	}
}

static qaws_vec3 vec3_from(qaws_scalar const* a, unsigned int dim)
{
	return qaws_v3(a[0], a[1], dim > 2 ? a[2] : QAWS_ZERO);
}

static qaws_status rational_tangent_sample(
	curve_sample const* s,
	unsigned int dim,
	qaws_scalar t_dot,
	unsigned int channels,
	qaws_diff_views const* views,
	int want_second,
	curve_jet_buf* primal,
	curve_jet_buf* tangent,
	curve_jet_buf* tangent2)
{
	qaws_field_view const* pv = views ? qaws_diff_views_find(views, s->support.ranges[0].field) : NULL;
	qaws_field_view const* wv = views ? qaws_diff_views_find(views, s->support.weight_field) : NULL;
	rational_row rows[QAWS_DIFF_MAX_ORDER + 1];
	qaws_dual3 A[QAWS_CURVE_JET_ORDER + 1], C[QAWS_CURVE_JET_ORDER + 1];
	qaws_dual1 W[QAWS_CURVE_JET_ORDER + 1];
	unsigned int kmax = highest_channel(channels), k, i, r;
	unsigned int rows_needed = kmax + (t_dot != QAWS_ZERO ? (want_second ? 2u : 1u) : 0u);
	(void)tangent2;

	for (r = 0; r <= rows_needed; r++)
		rational_sums(s, dim, r, pv, wv, &rows[r]);

	for (k = 0; k <= kmax; k++)
	{
		qaws_vec3 a1 = qaws_v3_zero(), ad1 = qaws_v3_zero(), a2 = qaws_v3_zero();
		qaws_scalar w1 = 0, wd1 = 0, w2 = 0;
		if (k + 1 <= rows_needed) { a1 = vec3_from(rows[k + 1].A, dim); ad1 = vec3_from(rows[k + 1].Ad, dim); w1 = rows[k + 1].W; wd1 = rows[k + 1].Wd; }
		if (k + 2 <= rows_needed) { a2 = vec3_from(rows[k + 2].A, dim); w2 = rows[k + 2].W; }

		A[k].v = vec3_from(rows[k].A, dim);
		A[k].t = qaws_v3_axpy(vec3_from(rows[k].Ad, dim), a1, t_dot);
		A[k].tt = qaws_v3_axpy(qaws_v3_axpy(vec3_from(rows[k].Add, dim), ad1, QAWS_LITERAL(2.0) * t_dot), a2, t_dot * t_dot);
		W[k].v = rows[k].W;
		W[k].t = rows[k].Wd + t_dot * w1;
		W[k].tt = QAWS_LITERAL(2.0) * t_dot * wd1 + t_dot * t_dot * w2;
	}

	if (W[0].v < QAWS_LITERAL(1e-15) && W[0].v > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;

	for (k = 0; k <= kmax; k++)
	{
		qaws_dual3 num = A[k];
		for (i = 1; i <= k; i++)
			num = qaws_dual3_sub(num, qaws_dual3_mul_const(qaws_dual3_scale(C[k - i], W[i]), g_binom[k][i]));
		C[k] = qaws_dual3_div(num, W[0]);
	}

	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		if (!(channels & (1u << k)))
			continue;
		primal->d[k][0] = C[k].v.x; primal->d[k][1] = C[k].v.y; primal->d[k][2] = C[k].v.z;
		tangent->d[k][0] = C[k].t.x; tangent->d[k][1] = C[k].t.y; tangent->d[k][2] = C[k].t.z;
		if (tangent2)
		{
			tangent2->d[k][0] = C[k].tt.x; tangent2->d[k][1] = C[k].tt.y; tangent2->d[k][2] = C[k].tt.z;
		}
	}
	if (dim == 2)
		for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
		{
			primal->d[k][2] = tangent->d[k][2] = QAWS_ZERO;
			if (tangent2) tangent2->d[k][2] = QAWS_ZERO;
		}
	return QAWS_STATUS_OK;
}

/* Reverse of the quotient recurrence: entries for control points and
   weights plus the coordinate adjoint. */
static qaws_status rational_adjoint_sample(
	curve_sample const* s,
	unsigned int dim,
	unsigned int channels,
	curve_jet_buf const* ybar,
	int want_coordinate,
	qaws_diff_entry* entries,
	unsigned int capacity,
	unsigned int* out_count,
	qaws_scalar* out_t_adjoint)
{
	rational_row rows[QAWS_CURVE_JET_ORDER + 2];
	qaws_scalar C[QAWS_CURVE_JET_ORDER + 1][3];
	qaws_scalar Cbar[QAWS_CURVE_JET_ORDER + 1][3];
	qaws_scalar Abar[QAWS_CURVE_JET_ORDER + 1][3];
	qaws_scalar Wbar[QAWS_CURVE_JET_ORDER + 1];
	qaws_support_range const* rg = &s->support.ranges[0];
	unsigned int kmax = highest_channel(channels), k, i, c, j, n = 0;
	qaws_scalar W0;

	for (k = 0; k <= kmax + (want_coordinate ? 1u : 0u); k++)
		rational_sums(s, dim, k, NULL, NULL, &rows[k]);
	W0 = rows[0].W;
	if (W0 < QAWS_LITERAL(1e-15) && W0 > -QAWS_LITERAL(1e-15))
		return QAWS_STATUS_DEGENERATE_CURVE;

	for (k = 0; k <= kmax; k++)
	{
		for (c = 0; c < dim; c++)
		{
			qaws_scalar num = rows[k].A[c];
			for (i = 1; i <= k; i++)
				num -= g_binom[k][i] * rows[i].W * C[k - i][c];
			C[k][c] = num / W0;
		}
	}

	memset(Abar, 0, sizeof(Abar));
	memset(Wbar, 0, sizeof(Wbar));
	for (k = 0; k <= kmax; k++)
		for (c = 0; c < dim; c++)
			Cbar[k][c] = (channels & (1u << k)) ? ybar->d[k][c] : QAWS_ZERO;

	for (k = kmax + 1; k-- > 0;)
	{
		qaws_scalar g[3], gc = 0;
		for (c = 0; c < dim; c++)
		{
			g[c] = Cbar[k][c] / W0;
			Abar[k][c] += g[c];
			gc += g[c] * C[k][c];
		}
		Wbar[0] -= gc;
		for (i = 1; i <= k; i++)
		{
			qaws_scalar dot = 0;
			for (c = 0; c < dim; c++)
				dot += g[c] * C[k - i][c];
			Wbar[i] -= g_binom[k][i] * dot;
			for (c = 0; c < dim; c++)
				Cbar[k - i][c] -= g_binom[k][i] * rows[i].W * g[c];
		}
	}

	for (j = 0; j < rg->count; j++)
	{
		unsigned int e = rg->first + j;
		qaws_scalar w = s->weights[e];
		qaws_scalar const* p = s->primal[0] + (size_t)e * dim;
		qaws_diff_entry* ep;
		qaws_diff_entry* ew;
		if (n + 2 > capacity)
			return QAWS_STATUS_INTERNAL_ERROR;
		ep = &entries[n++];
		ew = &entries[n++];
		ep->field = rg->field;
		ep->element = e;
		ep->g[0] = ep->g[1] = ep->g[2] = QAWS_ZERO;
		ew->field = s->support.weight_field;
		ew->element = e;
		ew->g[0] = ew->g[1] = ew->g[2] = QAWS_ZERO;
		for (k = 0; k <= kmax; k++)
		{
			qaws_scalar nk = s->support.weights[0][k][j];
			qaws_scalar pa = 0;
			for (c = 0; c < dim; c++)
			{
				ep->g[c] += nk * w * Abar[k][c];
				pa += p[c] * Abar[k][c];
			}
			ew->g[0] += nk * (pa + Wbar[k]);
		}
	}
	*out_count = n;

	if (want_coordinate)
	{
		qaws_scalar acc = 0;
		for (k = 0; k <= kmax; k++)
		{
			for (c = 0; c < dim; c++)
				acc += Abar[k][c] * rows[k + 1].A[c];
			acc += Wbar[k] * rows[k + 1].W;
		}
		*out_t_adjoint += acc;
	}
	return QAWS_STATUS_OK;
}

/* Adds g into the view element (masked). */
void qaws_internal_view_add(qaws_field_view* v, unsigned int e, unsigned int dim, qaws_scalar const* g)
{
	unsigned int c, stride = qaws_internal_view_stride(v);
	if (!qaws_internal_view_element_active(v, e))
		return;
	for (c = 0; c < dim; c++)
		if (qaws_internal_view_component_active(v, c))
			v->data[(size_t)e * stride + c] += g[c];
}

/* ------------------------------------------------------------------ */
/*  Batch adjoint with three accumulation strategies                  */
/* ------------------------------------------------------------------ */

typedef qaws_status (*jet_reader_fn)(void const* jets, unsigned int i, curve_jet_buf* out);

static qaws_status read_jet_2d(void const* jets, unsigned int i, curve_jet_buf* out)
{
	qaws_curve_jet_2d const* j = (qaws_curve_jet_2d const*)jets + i;
	unsigned int k;
	memset(out, 0, sizeof(*out));
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		out->d[k][0] = j->d[k].x;
		out->d[k][1] = j->d[k].y;
	}
	return QAWS_STATUS_OK;
}

static qaws_status read_jet_3d(void const* jets, unsigned int i, curve_jet_buf* out)
{
	qaws_curve_jet_3d const* j = (qaws_curve_jet_3d const*)jets + i;
	unsigned int k;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		out->d[k][0] = j->d[k].x;
		out->d[k][1] = j->d[k].y;
		out->d[k][2] = j->d[k].z;
	}
	return QAWS_STATUS_OK;
}

static void* scratch_alloc(qaws_diff_context const* ctx, size_t size)
{
	return qaws_internal_alloc(ctx ? ctx->allocator : NULL, (unsigned long)size);
}

static void scratch_free(qaws_diff_context const* ctx, void* p)
{
	qaws_internal_dealloc(ctx ? ctx->allocator : NULL, p);
}

unsigned int qaws_internal_field_components(qaws_diff_field field, unsigned int dim)
{
	switch (field)
	{
	case QAWS_FIELD_WEIGHTS:
	case QAWS_FIELD_KNOTS:
	case QAWS_FIELD_U_KNOTS:
	case QAWS_FIELD_V_KNOTS:
	case QAWS_FIELD_KEY_TIMES:
	case QAWS_FIELD_RADIUS:
	case QAWS_FIELD_RADIUS_B:
	case QAWS_FIELD_ANGLE_START:
	case QAWS_FIELD_ANGLE_END:
	case QAWS_FIELD_OFFSET_DISTANCE:
	case QAWS_FIELD_SCALE:
	case QAWS_FIELD_CURVATURE:
	case QAWS_FIELD_CURVATURE_RATE:
	case QAWS_FIELD_PARAMETER:
		return 1;
	default:
		return dim;
	}
}

qaws_status qaws_internal_check_views(qaws_diff_views const* views, unsigned int dim)
{
	unsigned int i;
	if (!views)
		return QAWS_STATUS_OK;
	for (i = 0; i < views->field_count; i++)
		if (views->fields[i].data &&
		    views->fields[i].components != qaws_internal_field_components(views->fields[i].field, dim))
			return QAWS_STATUS_INVALID_ARGUMENT;
	return QAWS_STATUS_OK;
}

/* View read with masks: inactive elements/components read as zero. */
void qaws_internal_view_read(qaws_field_view const* v, unsigned int e,
	unsigned int components, qaws_scalar* out)
{
	unsigned int c, stride = qaws_internal_view_stride(v);
	for (c = 0; c < components; c++)
		out[c] = QAWS_ZERO;
	if (!qaws_internal_view_element_active(v, e))
		return;
	for (c = 0; c < components && c < v->components; c++)
		if (qaws_internal_view_component_active(v, c))
			out[c] = v->data[(size_t)e * stride + c];
}

/* ------------------------------------------------------------------ */
/*  Generic accumulation: scatter, tiled, gather                      */
/* ------------------------------------------------------------------ */

#define DIFF_MAX_VIEW_FIELDS 16u

static int views_have_data(qaws_diff_views const* views)
{
	unsigned int i;
	if (!views)
		return 0;
	for (i = 0; i < views->field_count; i++)
		if (views->fields[i].data && views->fields[i].count)
			return 1;
	return 0;
}

static int view_index(qaws_diff_views const* views, qaws_diff_field field)
{
	unsigned int i;
	for (i = 0; i < views->field_count; i++)
		if (views->fields[i].field == field && views->fields[i].data)
			return (int)i;
	return -1;
}

/* Strategy: scatter. Each sample adds into the parameter adjoints. */
static qaws_status accumulate_scatter(
	unsigned int count,
	qaws_diff_entry* entries, unsigned int capacity,
	qaws_diff_collect_fn collect, void const* user, qaws_diff_views* views)
{
	unsigned int i, n, q;
	for (i = 0; i < count; i++)
	{
		qaws_status st = collect(user, i, 1, entries, capacity, &n);
		if (st != QAWS_STATUS_OK)
			return st;
		for (q = 0; q < n; q++)
		{
			int f = view_index(views, entries[q].field);
			if (f >= 0)
				qaws_internal_view_add(&views->fields[f], entries[q].element, views->fields[f].components, entries[q].g);
		}
	}
	return QAWS_STATUS_OK;
}

/*
 * Strategy: tiled. Every tile accumulates into a dense local buffer (the
 * CPU analogue of workgroup shared memory) and flushes only the touched
 * window of each field, once per tile.
 */
static qaws_status accumulate_tiled(
	qaws_diff_context const* ctx,
	unsigned int count,
	qaws_diff_entry* entries, unsigned int capacity,
	qaws_diff_collect_fn collect, void const* user, qaws_diff_views* views)
{
	unsigned int tile = (ctx && ctx->tile_size) ? ctx->tile_size : DIFF_DEFAULT_TILE;
	qaws_scalar* local[DIFF_MAX_VIEW_FIELDS];
	unsigned int lo[DIFF_MAX_VIEW_FIELDS], hi[DIFF_MAX_VIEW_FIELDS];
	unsigned int nf = views->field_count, f, i, q, c, base, n;
	qaws_status st = QAWS_STATUS_OK;

	if (nf > DIFF_MAX_VIEW_FIELDS)
		return accumulate_scatter(count, entries, capacity, collect, user, views);

	for (f = 0; f < DIFF_MAX_VIEW_FIELDS; f++)
		local[f] = NULL;
	for (f = 0; f < nf; f++)
	{
		qaws_field_view const* v = &views->fields[f];
		size_t size = sizeof(qaws_scalar) * (size_t)v->count * v->components;
		if (!v->data || !v->count)
			continue;
		local[f] = (qaws_scalar*)scratch_alloc(ctx, size);
		if (!local[f])
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto cleanup;
		}
		memset(local[f], 0, size);
	}

	for (base = 0; base < count; base += tile)
	{
		unsigned int end = (base + tile < count) ? base + tile : count;
		for (f = 0; f < nf; f++)
		{
			lo[f] = ~0u;
			hi[f] = 0u;
		}

		for (i = base; i < end; i++)
		{
			st = collect(user, i, 1, entries, capacity, &n);
			if (st != QAWS_STATUS_OK)
				goto cleanup;
			for (q = 0; q < n; q++)
			{
				int fi = view_index(views, entries[q].field);
				unsigned int e = entries[q].element;
				if (fi < 0 || !local[fi] || e >= views->fields[fi].count)
					continue;
				for (c = 0; c < views->fields[fi].components; c++)
					local[fi][(size_t)e * views->fields[fi].components + c] += entries[q].g[c];
				if (e < lo[fi]) lo[fi] = e;
				if (e + 1 > hi[fi]) hi[fi] = e + 1;
			}
		}

		for (f = 0; f < nf; f++)
		{
			unsigned int e;
			if (!local[f] || lo[f] == ~0u)
				continue;
			for (e = lo[f]; e < hi[f]; e++)
			{
				qaws_internal_view_add(&views->fields[f], e, views->fields[f].components, &local[f][(size_t)e * views->fields[f].components]);
				for (c = 0; c < views->fields[f].components; c++)
					local[f][(size_t)e * views->fields[f].components + c] = QAWS_ZERO;
			}
		}
	}

cleanup:
	for (f = 0; f < nf; f++)
		if (local[f])
			scratch_free(ctx, local[f]);
	return st;
}

/*
 * Strategy: gather. Contributions are bucketed per parameter element
 * (counting sort), then every element sums its own bucket in sample
 * order. Deterministic and free of write conflicts: the layout a
 * control-point-centric GPU kernel consumes.
 */
static qaws_status accumulate_gather(
	qaws_diff_context const* ctx,
	unsigned int count,
	qaws_diff_entry* entries, unsigned int capacity,
	qaws_diff_collect_fn collect, void const* user, qaws_diff_views* views)
{
	unsigned int* offsets[DIFF_MAX_VIEW_FIELDS];
	unsigned int* cursor[DIFF_MAX_VIEW_FIELDS];
	qaws_scalar* contrib[DIFF_MAX_VIEW_FIELDS];
	unsigned int nf = views->field_count, f, i, q, c, n, e;
	qaws_status st = QAWS_STATUS_OK;

	if (nf > DIFF_MAX_VIEW_FIELDS)
		return accumulate_scatter(count, entries, capacity, collect, user, views);

	for (f = 0; f < DIFF_MAX_VIEW_FIELDS; f++)
	{
		offsets[f] = NULL;
		cursor[f] = NULL;
		contrib[f] = NULL;
	}
	for (f = 0; f < nf; f++)
	{
		qaws_field_view const* v = &views->fields[f];
		size_t size = sizeof(unsigned int) * ((size_t)v->count + 1);
		if (!v->data || !v->count)
			continue;
		offsets[f] = (unsigned int*)scratch_alloc(ctx, size);
		cursor[f] = (unsigned int*)scratch_alloc(ctx, size);
		if (!offsets[f] || !cursor[f])
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto cleanup;
		}
		memset(offsets[f], 0, size);
	}

	/* Pass 1: count entries per element (coordinate adjoints written here). */
	for (i = 0; i < count; i++)
	{
		st = collect(user, i, 1, entries, capacity, &n);
		if (st != QAWS_STATUS_OK)
			goto cleanup;
		for (q = 0; q < n; q++)
		{
			int fi = view_index(views, entries[q].field);
			if (fi >= 0 && offsets[fi] && entries[q].element < views->fields[fi].count)
				offsets[fi][entries[q].element + 1]++;
		}
	}

	for (f = 0; f < nf; f++)
	{
		unsigned int total;
		if (!offsets[f])
			continue;
		for (e = 0; e < views->fields[f].count; e++)
			offsets[f][e + 1] += offsets[f][e];
		total = offsets[f][views->fields[f].count];
		memcpy(cursor[f], offsets[f], sizeof(unsigned int) * ((size_t)views->fields[f].count + 1));
		contrib[f] = (qaws_scalar*)scratch_alloc(ctx, sizeof(qaws_scalar) * ((size_t)total * views->fields[f].components + 1));
		if (!contrib[f])
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto cleanup;
		}
	}

	/* Pass 2: store contributions bucketed by element, in sample order. */
	for (i = 0; i < count; i++)
	{
		st = collect(user, i, 0, entries, capacity, &n);
		if (st != QAWS_STATUS_OK)
			goto cleanup;
		for (q = 0; q < n; q++)
		{
			int fi = view_index(views, entries[q].field);
			if (fi < 0 || !contrib[fi] || entries[q].element >= views->fields[fi].count)
				continue;
			e = entries[q].element;
			for (c = 0; c < views->fields[fi].components; c++)
				contrib[fi][(size_t)cursor[fi][e] * views->fields[fi].components + c] = entries[q].g[c];
			cursor[fi][e]++;
		}
	}

	/* Pass 3: every element gathers its bucket. */
	for (f = 0; f < nf; f++)
	{
		if (!contrib[f])
			continue;
		for (e = 0; e < views->fields[f].count; e++)
		{
			qaws_scalar g[3] = { 0, 0, 0 };
			unsigned int k;
			if (offsets[f][e + 1] == offsets[f][e])
				continue;
			for (k = offsets[f][e]; k < offsets[f][e + 1]; k++)
				for (c = 0; c < views->fields[f].components; c++)
					g[c] += contrib[f][(size_t)k * views->fields[f].components + c];
			qaws_internal_view_add(&views->fields[f], e, views->fields[f].components, g);
		}
	}

cleanup:
	for (f = 0; f < nf; f++)
	{
		if (offsets[f]) scratch_free(ctx, offsets[f]);
		if (cursor[f]) scratch_free(ctx, cursor[f]);
		if (contrib[f]) scratch_free(ctx, contrib[f]);
	}
	return st;
}

qaws_status qaws_internal_diff_accumulate(
	qaws_diff_context const* ctx,
	unsigned int components,
	unsigned int sample_count,
	unsigned int entry_capacity,
	qaws_diff_collect_fn collect,
	void const* user,
	qaws_diff_views* views)
{
	qaws_diff_entry* entries;
	qaws_status st;
	unsigned int i, n;

	(void)components;

	entries = (qaws_diff_entry*)scratch_alloc(ctx, sizeof(qaws_diff_entry) * (size_t)(entry_capacity ? entry_capacity : 1u));
	if (!entries)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	if (!views_have_data(views))
	{
		/* Only coordinate adjoints are requested. */
		st = QAWS_STATUS_OK;
		for (i = 0; i < sample_count && st == QAWS_STATUS_OK; i++)
			st = collect(user, i, 1, entries, entry_capacity, &n);
	}
	else
	{
		switch (ctx ? ctx->accumulation : QAWS_ACCUMULATE_SCATTER)
		{
		case QAWS_ACCUMULATE_TILED:
			st = accumulate_tiled(ctx, sample_count, entries, entry_capacity, collect, user, views);
			break;
		case QAWS_ACCUMULATE_GATHER:
			st = accumulate_gather(ctx, sample_count, entries, entry_capacity, collect, user, views);
			break;
		case QAWS_ACCUMULATE_SCATTER:
		default:
			st = accumulate_scatter(sample_count, entries, entry_capacity, collect, user, views);
			break;
		}
	}

	scratch_free(ctx, entries);
	return st;
}

/* ------------------------------------------------------------------ */
/*  Curve adjoint collector                                           */
/* ------------------------------------------------------------------ */

typedef struct curve_adjoint_job
{
	qaws_diff_context const* ctx;
	qaws_curve const* curve;
	unsigned int dim;
	qaws_scalar const* t;
	unsigned int channels;
	unsigned int order;
	void const* jets;
	jet_reader_fn read_jet;
	qaws_scalar* t_adjoint;
} curve_adjoint_job;

static qaws_status curve_collect(
	void const* user,
	unsigned int i,
	int first_pass,
	qaws_diff_entry* entries,
	unsigned int capacity,
	unsigned int* out_count)
{
	curve_adjoint_job const* job = (curve_adjoint_job const*)user;
	curve_sample s;
	curve_jet_buf ybar;
	unsigned int r, j, n = 0;
	qaws_curve_diff_vtable const* d = curve_diff(job->curve);
	qaws_status st = curve_prepare(job->curve, job->dim, job->t[i], job->order, &s);

	*out_count = 0;
	if (st != QAWS_STATUS_OK)
		return st;
	job->read_jet(job->jets, i, &ybar);

	if (s.weights)
	{
		qaws_scalar dummy = 0;
		int want_t = first_pass && job->t_adjoint != NULL;
		st = rational_adjoint_sample(&s, job->dim, job->channels, &ybar, want_t,
			entries, capacity, &n, want_t ? &job->t_adjoint[i] : &dummy);
		if (st != QAWS_STATUS_OK)
			return st;
	}
	else
	{
		for (r = 0; r < s.support.range_count; r++)
		{
			for (j = 0; j < s.support.ranges[r].count; j++)
			{
				if (n >= capacity)
					return QAWS_STATUS_INTERNAL_ERROR;
				entries[n].field = s.support.ranges[r].field;
				entries[n].element = s.support.ranges[r].first + j;
				entries[n].g[2] = QAWS_ZERO;
				adjoint_contribution(&s, r, j, job->dim, job->channels, &ybar, entries[n].g);
				n++;
			}
		}
		if (first_pass && job->t_adjoint)
			job->t_adjoint[i] += coordinate_adjoint(&s, job->dim, job->channels, &ybar);
	}
	*out_count = n;

	if (first_pass)
	{
		int has_t = job->t_adjoint != NULL;
		qaws_internal_diff_report_note(job->ctx,
			has_t ? d->diff_class : QAWS_DIFF_SMOOTH,
			(has_t && s.at_boundary) ? QAWS_DIFF_AT_BOUNDARY : QAWS_DIFF_VALID,
			(has_t && d->diff_class == QAWS_DIFF_PIECEWISE_SMOOTH) ? (unsigned int)QAWS_FREEZE_SPAN : 0u,
			i);
	}
	return QAWS_STATUS_OK;
}

static qaws_status curve_batch_adjoint(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	unsigned int dim,
	qaws_scalar const* t,
	unsigned int count,
	unsigned int channels,
	void const* jets,
	jet_reader_fn read_jet,
	qaws_diff_views* views,
	qaws_scalar* t_adjoint)
{
	curve_adjoint_job job;
	qaws_status st;

	if (!curve || !t || !jets)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((unsigned int)curve->dimension != dim)
		return QAWS_STATUS_INVALID_DIMENSION;
	if (!curve_diff(curve))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	st = qaws_internal_check_views(views, dim);
	if (st != QAWS_STATUS_OK)
		return st;

	if (curve_diff(curve)->adjoint_span)
	{
		unsigned int i;
		for (i = 0; i < count; i++)
		{
			curve_jet_buf buf;
			qaws_curve_jet_3d j3;
			qaws_scalar local_t;
			unsigned int k, span = qaws_internal_find_span(curve, t[i], &local_t);
			read_jet(jets, i, &buf);
			for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
			{
				j3.d[k].x = buf.d[k][0];
				j3.d[k].y = buf.d[k][1];
				j3.d[k].z = dim == 3 ? buf.d[k][2] : QAWS_ZERO;
			}
			j3.channels = channels & 0xFu;
			st = curve_diff(curve)->adjoint_span(ctx, curve, span, local_t, channels & 0xFu, &j3, views,
				t_adjoint ? &t_adjoint[i] : NULL);
			if (st != QAWS_STATUS_OK)
				return st;
		}
		return QAWS_STATUS_OK;
	}

	job.ctx = ctx;
	job.curve = curve;
	job.dim = dim;
	job.t = t;
	job.channels = channels & 0xFu;
	job.order = highest_channel(job.channels) + (t_adjoint ? 1u : 0u);
	job.jets = jets;
	job.read_jet = read_jet;
	job.t_adjoint = t_adjoint;

	return qaws_internal_diff_accumulate(ctx, dim, count,
		QAWS_DIFF_MAX_RANGES * QAWS_DIFF_MAX_SUPPORT, curve_collect, &job, views);
}

/* ------------------------------------------------------------------ */
/*  Batch tangent                                                     */
/* ------------------------------------------------------------------ */

static void store_jet_2d(qaws_curve_jet_2d* out, curve_jet_buf const* b, unsigned int channels)
{
	unsigned int k;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		out->d[k].x = b->d[k][0];
		out->d[k].y = b->d[k][1];
	}
	out->channels = channels;
}

static void store_jet_3d(qaws_curve_jet_3d* out, curve_jet_buf const* b, unsigned int channels)
{
	unsigned int k;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		out->d[k].x = b->d[k][0];
		out->d[k].y = b->d[k][1];
		out->d[k].z = b->d[k][2];
	}
	out->channels = channels;
}

static void store_lifted_2d(qaws_curve_jet_2d* out, qaws_curve_jet_3d const* in)
{
	unsigned int k;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		out->d[k].x = in->d[k].x;
		out->d[k].y = in->d[k].y;
	}
	out->channels = in->channels;
}

static qaws_status curve_batch_tangent(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	unsigned int dim,
	qaws_scalar const* t,
	qaws_scalar const* t_tangent,
	unsigned int count,
	unsigned int channels,
	qaws_diff_views const* views,
	void* out_primal,
	void* out_tangent,
	void* out_tangent2)
{
	unsigned int i;
	qaws_status st;
	qaws_curve_diff_vtable const* d;

	if (!curve || !t || !out_tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((unsigned int)curve->dimension != dim)
		return QAWS_STATUS_INVALID_DIMENSION;
	d = curve_diff(curve);
	if (!d)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (out_tangent2 && !(d->capabilities & QAWS_CAP_TANGENT2))
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	st = qaws_internal_check_views(views, dim);
	if (st != QAWS_STATUS_OK)
		return st;
	channels &= 0xFu;

	if (d->tangent_span)
	{
		/* Derived curve: direct rule chaining into its children. */
		for (i = 0; i < count; i++)
		{
			qaws_curve_jet_3d p3, t3, tt3;
			qaws_scalar local_t;
			unsigned int span = qaws_internal_find_span(curve, t[i], &local_t);
			qaws_scalar t_dot = t_tangent ? t_tangent[i] : QAWS_ZERO;
			qaws_scalar const eps = QAWS_LITERAL(16.0) * QAWS_EPSILON;
			int boundary = (span > 0 && local_t <= eps) || (span + 1 < curve->span_count && local_t >= QAWS_ONE - eps);
			st = d->tangent_span(ctx, curve, span, local_t, t_dot, channels, views, &p3, &t3, out_tangent2 ? &tt3 : NULL);
			if (st != QAWS_STATUS_OK)
				return st;
			if (dim == 2)
			{
				if (out_primal) store_lifted_2d((qaws_curve_jet_2d*)out_primal + i, &p3);
				store_lifted_2d((qaws_curve_jet_2d*)out_tangent + i, &t3);
				if (out_tangent2) store_lifted_2d((qaws_curve_jet_2d*)out_tangent2 + i, &tt3);
			}
			else
			{
				if (out_primal) ((qaws_curve_jet_3d*)out_primal)[i] = p3;
				((qaws_curve_jet_3d*)out_tangent)[i] = t3;
				if (out_tangent2) ((qaws_curve_jet_3d*)out_tangent2)[i] = tt3;
			}
			qaws_internal_diff_report_note(ctx,
				t_dot != QAWS_ZERO ? d->diff_class : QAWS_DIFF_SMOOTH,
				(t_dot != QAWS_ZERO && boundary) ? QAWS_DIFF_AT_BOUNDARY : QAWS_DIFF_VALID,
				t_dot != QAWS_ZERO ? (unsigned int)QAWS_FREEZE_SPAN : 0u, i);
		}
		return QAWS_STATUS_OK;
	}

	for (i = 0; i < count; i++)
	{
		curve_jet_buf p, tg, tg2;
		qaws_scalar t_dot = t_tangent ? t_tangent[i] : QAWS_ZERO;
		st = curve_tangent_sample(ctx, curve, dim, t[i], t_dot, channels, views,
			out_tangent2 != NULL, i, &p, &tg, out_tangent2 ? &tg2 : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		if (dim == 2)
		{
			if (out_primal) store_jet_2d((qaws_curve_jet_2d*)out_primal + i, &p, channels);
			store_jet_2d((qaws_curve_jet_2d*)out_tangent + i, &tg, channels);
			if (out_tangent2) store_jet_2d((qaws_curve_jet_2d*)out_tangent2 + i, &tg2, channels);
		}
		else
		{
			if (out_primal) store_jet_3d((qaws_curve_jet_3d*)out_primal + i, &p, channels);
			store_jet_3d((qaws_curve_jet_3d*)out_tangent + i, &tg, channels);
			if (out_tangent2) store_jet_3d((qaws_curve_jet_3d*)out_tangent2 + i, &tg2, channels);
		}
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_eval_batch_tangent_2d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar const* t, qaws_scalar const* t_tangent, unsigned int count,
	unsigned int channels, qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal, qaws_curve_jet_2d* out_tangent)
{
	return curve_batch_tangent(ctx, curve, 2, t, t_tangent, count, channels,
		param_tangent, out_primal, out_tangent, NULL);
}

qaws_status qaws_curve_eval_batch_tangent_3d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar const* t, qaws_scalar const* t_tangent, unsigned int count,
	unsigned int channels, qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal, qaws_curve_jet_3d* out_tangent)
{
	return curve_batch_tangent(ctx, curve, 3, t, t_tangent, count, channels,
		param_tangent, out_primal, out_tangent, NULL);
}

qaws_status qaws_curve_eval_batch_tangent2_2d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar const* t, qaws_scalar const* t_tangent, unsigned int count,
	unsigned int channels, qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal, qaws_curve_jet_2d* out_tangent,
	qaws_curve_jet_2d* out_tangent2)
{
	if (!out_tangent2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return curve_batch_tangent(ctx, curve, 2, t, t_tangent, count, channels,
		param_tangent, out_primal, out_tangent, out_tangent2);
}

qaws_status qaws_curve_eval_batch_tangent2_3d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar const* t, qaws_scalar const* t_tangent, unsigned int count,
	unsigned int channels, qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal, qaws_curve_jet_3d* out_tangent,
	qaws_curve_jet_3d* out_tangent2)
{
	if (!out_tangent2)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return curve_batch_tangent(ctx, curve, 3, t, t_tangent, count, channels,
		param_tangent, out_primal, out_tangent, out_tangent2);
}

qaws_status qaws_curve_eval_batch_adjoint_2d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar const* t, unsigned int count, unsigned int channels,
	qaws_curve_jet_2d const* out_adjoint,
	qaws_diff_views* param_adjoint, qaws_scalar* t_adjoint)
{
	return curve_batch_adjoint(ctx, curve, 2, t, count, channels,
		out_adjoint, read_jet_2d, param_adjoint, t_adjoint);
}

qaws_status qaws_curve_eval_batch_adjoint_3d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar const* t, unsigned int count, unsigned int channels,
	qaws_curve_jet_3d const* out_adjoint,
	qaws_diff_views* param_adjoint, qaws_scalar* t_adjoint)
{
	return curve_batch_adjoint(ctx, curve, 3, t, count, channels,
		out_adjoint, read_jet_3d, param_adjoint, t_adjoint);
}

/* ------------------------------------------------------------------ */
/*  Single sample wrappers                                            */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_eval_tangent_2d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar t, qaws_scalar t_tangent, unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_2d* out_primal, qaws_curve_jet_2d* out_tangent)
{
	return curve_batch_tangent(ctx, curve, 2, &t, &t_tangent, 1, channels,
		param_tangent, out_primal, out_tangent, NULL);
}

qaws_status qaws_curve_eval_tangent_3d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar t, qaws_scalar t_tangent, unsigned int channels,
	qaws_diff_views const* param_tangent,
	qaws_curve_jet_3d* out_primal, qaws_curve_jet_3d* out_tangent)
{
	return curve_batch_tangent(ctx, curve, 3, &t, &t_tangent, 1, channels,
		param_tangent, out_primal, out_tangent, NULL);
}

qaws_status qaws_curve_eval_adjoint_2d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar t, unsigned int channels,
	qaws_curve_jet_2d const* out_adjoint,
	qaws_diff_views* param_adjoint, qaws_scalar* t_adjoint)
{
	return curve_batch_adjoint(ctx, curve, 2, &t, 1, channels,
		out_adjoint, read_jet_2d, param_adjoint, t_adjoint);
}

qaws_status qaws_curve_eval_adjoint_3d(
	qaws_diff_context const* ctx, qaws_curve const* curve,
	qaws_scalar t, unsigned int channels,
	qaws_curve_jet_3d const* out_adjoint,
	qaws_diff_views* param_adjoint, qaws_scalar* t_adjoint)
{
	return curve_batch_adjoint(ctx, curve, 3, &t, 1, channels,
		out_adjoint, read_jet_3d, param_adjoint, t_adjoint);
}

/* ------------------------------------------------------------------ */
/*  Local support and support index                                   */
/* ------------------------------------------------------------------ */

qaws_status qaws_curve_local_support(
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int order,
	qaws_local_support* out_support)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	qaws_scalar local_t;
	unsigned int span;

	if (!curve || !out_support)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (order > QAWS_DIFF_MAX_ORDER)
		return QAWS_STATUS_OUT_OF_RANGE;
	if (!d || !d->linear_support)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	span = qaws_internal_find_span(curve, t, &local_t);
	return d->linear_support(curve, span, local_t, order, out_support);
}

qaws_status qaws_curve_build_support_index(
	qaws_curve const* curve,
	qaws_scalar const* t,
	unsigned int count,
	qaws_diff_field field,
	unsigned int* out_offsets,
	unsigned int offset_capacity,
	unsigned int* out_samples,
	unsigned int sample_capacity,
	unsigned int* out_entry_count)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	qaws_scalar const* data = NULL;
	unsigned int element_count = 0, components = 0, i, r, j, e, total;
	qaws_status st;

	if (!curve || !t || !out_offsets || !out_entry_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_entry_count = 0;
	if (!d || !d->linear_support || !d->primal_field)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	st = d->primal_field(curve, field, &data, &element_count, &components);
	if (st != QAWS_STATUS_OK)
		return st;
	if (offset_capacity < element_count + 1)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	memset(out_offsets, 0, sizeof(unsigned int) * ((size_t)element_count + 1));
	for (i = 0; i < count; i++)
	{
		qaws_local_support s;
		st = qaws_curve_local_support(curve, t[i], 0, &s);
		if (st != QAWS_STATUS_OK)
			return st;
		for (r = 0; r < s.range_count; r++)
		{
			if (s.ranges[r].field != field)
				continue;
			for (j = 0; j < s.ranges[r].count; j++)
			{
				e = s.ranges[r].first + j;
				if (e < element_count)
					out_offsets[e + 1]++;
			}
		}
	}
	for (e = 0; e < element_count; e++)
		out_offsets[e + 1] += out_offsets[e];
	total = out_offsets[element_count];
	*out_entry_count = total;
	if (!out_samples || sample_capacity < total)
		return QAWS_STATUS_BUFFER_TOO_SMALL;

	{
		/* Fill in sample order; use a moving cursor stored past the end of
		   each bucket start by recomputing offsets afterwards. */
		unsigned int* cursor = (unsigned int*)qaws_internal_alloc(NULL,
			(unsigned long)(sizeof(unsigned int) * ((size_t)element_count + 1)));
		if (!cursor)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		memcpy(cursor, out_offsets, sizeof(unsigned int) * ((size_t)element_count + 1));
		for (i = 0; i < count; i++)
		{
			qaws_local_support s;
			st = qaws_curve_local_support(curve, t[i], 0, &s);
			if (st != QAWS_STATUS_OK)
			{
				qaws_internal_dealloc(NULL, cursor);
				return st;
			}
			for (r = 0; r < s.range_count; r++)
			{
				if (s.ranges[r].field != field)
					continue;
				for (j = 0; j < s.ranges[r].count; j++)
				{
					e = s.ranges[r].first + j;
					if (e < element_count)
						out_samples[cursor[e]++] = i;
				}
			}
		}
		qaws_internal_dealloc(NULL, cursor);
	}
	(void)data;
	return QAWS_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/*  Child evaluation helpers for derived objects                      */
/* ------------------------------------------------------------------ */

qaws_diff_views const* qaws_internal_child_views(qaws_diff_views const* views, unsigned int i)
{
	if (!views || i >= views->child_count || !views->children)
		return NULL;
	return &views->children[i];
}

static void lift_jet(qaws_curve_jet_2d const* in, qaws_curve_jet_3d* out)
{
	unsigned int k;
	for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
	{
		out->d[k].x = in->d[k].x;
		out->d[k].y = in->d[k].y;
		out->d[k].z = QAWS_ZERO;
	}
	out->channels = in->channels;
}

qaws_status qaws_internal_curve_tangent_any(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	qaws_scalar t_dot,
	unsigned int channels,
	qaws_diff_views const* views,
	qaws_curve_jet_3d* primal,
	qaws_curve_jet_3d* tangent,
	qaws_curve_jet_3d* tangent2)
{
	if (!curve || !primal || !tangent)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension == QAWS_DIMENSION_3D)
		return curve_batch_tangent(ctx, curve, 3, &t, &t_dot, 1, channels, views,
			primal, tangent, tangent2);
	{
		qaws_curve_jet_2d p, tg, tt;
		qaws_status st = curve_batch_tangent(ctx, curve, 2, &t, &t_dot, 1, channels, views,
			&p, &tg, tangent2 ? &tt : NULL);
		if (st != QAWS_STATUS_OK)
			return st;
		lift_jet(&p, primal);
		lift_jet(&tg, tangent);
		if (tangent2)
			lift_jet(&tt, tangent2);
		return QAWS_STATUS_OK;
	}
}

qaws_status qaws_internal_curve_adjoint_any(
	qaws_diff_context const* ctx,
	qaws_curve const* curve,
	qaws_scalar t,
	unsigned int channels,
	qaws_curve_jet_3d const* jet_adjoint,
	qaws_diff_views* views,
	qaws_scalar* t_adjoint)
{
	if (!curve || !jet_adjoint)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (curve->dimension == QAWS_DIMENSION_3D)
		return curve_batch_adjoint(ctx, curve, 3, &t, 1, channels, jet_adjoint, read_jet_3d, views, t_adjoint);
	{
		qaws_curve_jet_2d j;
		unsigned int k;
		for (k = 0; k <= QAWS_CURVE_JET_ORDER; k++)
		{
			j.d[k].x = jet_adjoint->d[k].x;
			j.d[k].y = jet_adjoint->d[k].y;
		}
		j.channels = jet_adjoint->channels;
		return curve_batch_adjoint(ctx, curve, 2, &t, 1, channels, &j, read_jet_2d, views, t_adjoint);
	}
}

qaws_status qaws_curve_diff_children(
	qaws_curve const* curve,
	qaws_diff_child* out_children,
	unsigned int capacity,
	unsigned int* out_count)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	unsigned int n;
	if (!curve || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = 0;
	if (!d)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (!d->children)
		return QAWS_STATUS_OK;
	n = d->children(curve, out_children, out_children ? capacity : 0u);
	*out_count = n;
	return (out_children && n > capacity) ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}

qaws_status qaws_internal_field_override(
	qaws_diff_views const* values,
	qaws_diff_field field,
	qaws_scalar const* fallback,
	unsigned int count,
	unsigned int components,
	qaws_scalar const** out_values,
	qaws_scalar** out_owned)
{
	qaws_field_view const* v = values ? qaws_diff_views_find(values, field) : NULL;
	unsigned int e, c, stride;
	qaws_scalar* buf;

	*out_owned = NULL;
	*out_values = fallback;
	if (!v)
		return QAWS_STATUS_OK;
	if (v->count != count || v->components != components)
		return QAWS_STATUS_INVALID_ARGUMENT;
	buf = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (size_t)count * components + 1));
	if (!buf)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	stride = qaws_internal_view_stride(v);
	for (e = 0; e < count; e++)
		for (c = 0; c < components; c++)
			buf[(size_t)e * components + c] = v->data[(size_t)e * stride + c];
	*out_owned = buf;
	*out_values = buf;
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_clone_with_fields(
	qaws_curve const* curve,
	qaws_diff_views const* values,
	qaws_curve** out_curve)
{
	qaws_curve_diff_vtable const* d = curve_diff(curve);
	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_curve = NULL;
	if (!d || !d->rebuild)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	return d->rebuild(curve, values, out_curve);
}
