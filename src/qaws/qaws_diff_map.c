#include "qaws_diff_map.h"
#include "qaws_diff.h"
#include "qaws_curve.h"
#include "qaws_operations.h"
#include "qaws_convert.h"
#include "qaws_export.h"
#include "internal/qaws_internal_basis.h"
#include "internal/qaws_internal_fit.h"
#include "qaws_diff_geometry.h"
#include <math.h>
#include "internal/qaws_internal_types.h"
#include "internal/qaws_internal_curve.h"
#include "internal/qaws_internal_diff.h"
#include <string.h>

#define MAP_MAX_OBJECTS 8
#define MAP_MAX_FIELDS 8

struct qaws_diff_map
{
	qaws_diff_map_kind kind;
	unsigned int in_count, out_count;
	qaws_diff_map_entry* entries;
	unsigned int entry_count, capacity;
	unsigned int frozen_inputs;   /* bit i: input i has no derivative */
};

/* ================================================================== */
/*  Map object                                                        */
/* ================================================================== */

static qaws_diff_map* map_create(qaws_diff_map_kind kind, unsigned int in_count, unsigned int out_count)
{
	qaws_diff_map* m = (qaws_diff_map*)qaws_internal_alloc(NULL, (unsigned long)sizeof(qaws_diff_map));
	if (!m)
		return NULL;
	memset(m, 0, sizeof(*m));
	m->kind = kind;
	m->in_count = in_count;
	m->out_count = out_count;
	return m;
}

void qaws_diff_map_destroy(qaws_diff_map* map)
{
	if (!map)
		return;
	qaws_internal_dealloc(NULL, map->entries);
	qaws_internal_dealloc(NULL, map);
}

static qaws_status map_add(qaws_diff_map* m, unsigned int out_obj, qaws_diff_field out_field, unsigned int out_elem,
	unsigned int out_comp, unsigned int in_obj, qaws_diff_field in_field, unsigned int in_elem, unsigned int in_comp,
	qaws_scalar weight)
{
	qaws_diff_map_entry* e;
	if (m->entry_count == m->capacity)
	{
		unsigned int cap = m->capacity ? m->capacity * 2 : 64;
		qaws_diff_map_entry* grown = (qaws_diff_map_entry*)qaws_internal_alloc(NULL,
			(unsigned long)(sizeof(qaws_diff_map_entry) * cap));
		if (!grown)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		if (m->entries)
		{
			memcpy(grown, m->entries, sizeof(qaws_diff_map_entry) * m->entry_count);
			qaws_internal_dealloc(NULL, m->entries);
		}
		m->entries = grown;
		m->capacity = cap;
	}
	e = &m->entries[m->entry_count++];
	e->out_object = out_obj;
	e->out_field = (unsigned char)out_field;
	e->out_element = out_elem;
	e->out_component = (unsigned char)out_comp;
	e->in_object = in_obj;
	e->in_field = (unsigned char)in_field;
	e->in_element = in_elem;
	e->in_component = (unsigned char)in_comp;
	e->weight = weight;
	return QAWS_STATUS_OK;
}

qaws_diff_map_kind qaws_diff_map_get_kind(qaws_diff_map const* map)
{
	return map ? map->kind : QAWS_DIFF_MAP_LINEAR_SPARSE;
}

unsigned int qaws_diff_map_input_count(qaws_diff_map const* map)
{
	return map ? map->in_count : 0u;
}

unsigned int qaws_diff_map_output_count(qaws_diff_map const* map)
{
	return map ? map->out_count : 0u;
}

qaws_status qaws_diff_map_get_entries(qaws_diff_map const* map, qaws_diff_map_entry* out_entries,
	unsigned int capacity, unsigned int* out_count)
{
	if (!map || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_count = map->entry_count;
	if (!out_entries)
		return QAWS_STATUS_OK;
	if (capacity < map->entry_count)
		return QAWS_STATUS_BUFFER_TOO_SMALL;
	memcpy(out_entries, map->entries, sizeof(qaws_diff_map_entry) * map->entry_count);
	return QAWS_STATUS_OK;
}

/* Apply one entry: dst += w * src for the selected components. */
static void entry_apply(qaws_diff_map_entry const* e, qaws_field_view const* src, unsigned int src_elem,
	unsigned int src_comp, qaws_field_view* dst, unsigned int dst_elem, unsigned int dst_comp)
{
	qaws_scalar in[3], add[3] = { 0, 0, 0 };
	unsigned int c;
	qaws_internal_view_read(src, src_elem, src->components, in);
	if (src_comp == QAWS_DIFF_MAP_ALL_COMPONENTS && dst_comp == QAWS_DIFF_MAP_ALL_COMPONENTS)
	{
		for (c = 0; c < dst->components && c < src->components; c++)
			add[c] = e->weight * in[c];
	}
	else
	{
		unsigned int sc = src_comp == QAWS_DIFF_MAP_ALL_COMPONENTS ? 0u : src_comp;
		unsigned int dc = dst_comp == QAWS_DIFF_MAP_ALL_COMPONENTS ? 0u : dst_comp;
		if (sc < src->components && dc < dst->components)
			add[dc] = e->weight * in[sc];
	}
	qaws_internal_view_add(dst, dst_elem, dst->components, add);
}

qaws_status qaws_diff_map_tangent(
	qaws_diff_map const* map,
	qaws_diff_context const* ctx,
	qaws_diff_views const* const* in_tangents,
	unsigned int in_count,
	qaws_diff_views* const* out_tangents,
	unsigned int out_count)
{
	unsigned int i;
	(void)ctx;
	if (!map || (in_count && !in_tangents) || (out_count && !out_tangents))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < in_count && i < map->in_count; i++)
		if (in_tangents[i] && (map->frozen_inputs & (1u << i)))
			return QAWS_STATUS_UNSUPPORTED_OPERATION;
	for (i = 0; i < out_count; i++)
		if (out_tangents[i])
			qaws_diff_views_clear(out_tangents[i]);
	for (i = 0; i < map->entry_count; i++)
	{
		qaws_diff_map_entry const* e = &map->entries[i];
		qaws_field_view const* src;
		qaws_field_view* dst;
		if (e->in_object >= in_count || e->out_object >= out_count)
			continue;
		if (!in_tangents[e->in_object] || !out_tangents[e->out_object])
			continue;
		src = qaws_diff_views_find(in_tangents[e->in_object], (qaws_diff_field)e->in_field);
		dst = qaws_diff_views_find(out_tangents[e->out_object], (qaws_diff_field)e->out_field);
		if (src && dst)
			entry_apply(e, src, e->in_element, e->in_component, dst, e->out_element, e->out_component);
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_diff_map_adjoint(
	qaws_diff_map const* map,
	qaws_diff_context const* ctx,
	qaws_diff_views const* const* out_adjoints,
	unsigned int out_count,
	qaws_diff_views* const* in_adjoints,
	unsigned int in_count)
{
	unsigned int i;
	(void)ctx;
	if (!map || (in_count && !in_adjoints) || (out_count && !out_adjoints))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < map->entry_count; i++)
	{
		qaws_diff_map_entry const* e = &map->entries[i];
		qaws_field_view const* src;
		qaws_field_view* dst;
		if (e->in_object >= in_count || e->out_object >= out_count)
			continue;
		if (!in_adjoints[e->in_object] || !out_adjoints[e->out_object])
			continue;
		src = qaws_diff_views_find(out_adjoints[e->out_object], (qaws_diff_field)e->out_field);
		dst = qaws_diff_views_find(in_adjoints[e->in_object], (qaws_diff_field)e->in_field);
		if (src && dst)
			entry_apply(e, src, e->out_element, e->out_component, dst, e->in_element, e->in_component);
	}
	return QAWS_STATUS_OK;
}

/* ================================================================== */
/*  Superposition builder for linear constructions                    */
/*                                                                    */
/*  For a construction that is affine in the vector fields of its     */
/*  inputs (polynomial families), column k of the Jacobian is         */
/*  F(e_k) - F(0): exact, no finite step. Constructions act on every  */
/*  coordinate alike, so one unit per element gives the weight for    */
/*  all components. NURBS inputs are handled in homogeneous form.      */
/* ================================================================== */

typedef qaws_status (*construct_fn)(void const* user, qaws_curve const* const* in, qaws_curve** out);

typedef struct field_info
{
	qaws_diff_field field;
	unsigned int count, components;
} field_info;

/* Vector fields (linear) and whether the object carries weights. */
static unsigned int vector_fields(qaws_curve const* c, field_info* out, int* rational)
{
	qaws_field_desc fields[MAP_MAX_FIELDS];
	unsigned int n = 0, i, k = 0;
	*rational = 0;
	if (qaws_curve_describe_fields(c, fields, MAP_MAX_FIELDS, &n) != QAWS_STATUS_OK)
		return 0;
	for (i = 0; i < n && i < MAP_MAX_FIELDS; i++)
	{
		if (fields[i].field == QAWS_FIELD_WEIGHTS)
			*rational = 1;
		if ((fields[i].value_type == QAWS_VALUE_VEC2 || fields[i].value_type == QAWS_VALUE_VEC3) && fields[i].capabilities)
		{
			out[k].field = fields[i].field;
			out[k].count = fields[i].count;
			out[k].components = (unsigned int)fields[i].value_type;
			k++;
		}
	}
	return k;
}

static qaws_status read_values(qaws_curve const* c, qaws_diff_field f, qaws_scalar** out, unsigned int* n)
{
	qaws_status st = qaws_curve_read_field(c, f, NULL, 0, n);
	if (st != QAWS_STATUS_OK)
		return st;
	*out = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (*n + 1)));
	if (!*out)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	return qaws_curve_read_field(c, f, *out, *n, n);
}

/* Clone of `c` with every vector field zero (weights one), and optionally
   component 0 of one element of one field set to one. */
static qaws_status make_probe(qaws_curve const* c, field_info const* fields, unsigned int nf, int rational,
	int unit_field, unsigned int unit_elem, qaws_curve** out)
{
	qaws_field_view views[MAP_MAX_FIELDS + 1];
	qaws_scalar* buffers[MAP_MAX_FIELDS + 1];
	qaws_diff_views vs;
	unsigned int i, n = 0;
	qaws_status st = QAWS_STATUS_OK;

	for (i = 0; i < nf; i++)
	{
		size_t size = sizeof(qaws_scalar) * (size_t)fields[i].count * fields[i].components;
		buffers[n] = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(size + 1));
		if (!buffers[n])
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		memset(buffers[n], 0, size);
		if ((int)i == unit_field)
			buffers[n][(size_t)unit_elem * fields[i].components] = QAWS_ONE;
		views[n] = qaws_field_view_make(fields[i].field, buffers[n], fields[i].count, fields[i].components);
		n++;
	}
	if (rational)
	{
		unsigned int wc = fields[0].count;
		buffers[n] = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (wc + 1)));
		if (!buffers[n])
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			goto done;
		}
		for (i = 0; i < wc; i++)
			buffers[n][i] = QAWS_ONE;
		views[n] = qaws_field_view_make(QAWS_FIELD_WEIGHTS, buffers[n], wc, 1);
		n++;
	}
	vs.fields = views;
	vs.field_count = n;
	vs.children = NULL;
	vs.child_count = 0;
	st = qaws_curve_clone_with_fields(c, &vs, out);
done:
	for (i = 0; i < n; i++)
		qaws_internal_dealloc(NULL, buffers[i]);
	return st;
}

static void destroy_curves(qaws_curve** c, unsigned int n)
{
	unsigned int i;
	for (i = 0; i < n; i++)
		if (c[i])
		{
			qaws_curve_destroy(c[i]);
			c[i] = NULL;
		}
}

/* Accumulates entries for the columns of one probe: (F(probe) - F0). */
static qaws_status add_column(qaws_diff_map* m, qaws_curve* const* f0, qaws_curve* const* fu, unsigned int n_out,
	unsigned int in_obj, qaws_diff_field in_field, unsigned int in_elem)
{
	unsigned int o, g, k;
	qaws_status st;
	for (o = 0; o < n_out; o++)
	{
		field_info fields[MAP_MAX_FIELDS];
		int rational;
		unsigned int nf = vector_fields(f0[o], fields, &rational);
		for (g = 0; g < nf; g++)
		{
			qaws_scalar *a = NULL, *b = NULL;
			unsigned int na = 0, nb = 0;
			st = read_values(f0[o], fields[g].field, &a, &na);
			if (st == QAWS_STATUS_OK)
				st = read_values(fu[o], fields[g].field, &b, &nb);
			if (st == QAWS_STATUS_OK && na == nb)
			{
				qaws_scalar const tol = QAWS_LITERAL(64.0) * QAWS_EPSILON;
				for (k = 0; k < fields[g].count; k++)
				{
					qaws_scalar w = b[k * fields[g].components] - a[k * fields[g].components];
					if (QAWS_FABS(w) > tol)
					{
						st = map_add(m, o, fields[g].field, k, QAWS_DIFF_MAP_ALL_COMPONENTS,
							in_obj, in_field, in_elem, QAWS_DIFF_MAP_ALL_COMPONENTS, w);
						if (st != QAWS_STATUS_OK)
							break;
					}
				}
			}
			qaws_internal_dealloc(NULL, a);
			qaws_internal_dealloc(NULL, b);
			if (st != QAWS_STATUS_OK)
				return st;
		}
	}
	return QAWS_STATUS_OK;
}

/* Rewrites position entries of rational inputs/outputs into the projected
   form: P'_j = sum_k A_jk w_k P_k / w'_j,  w'_j = sum_k A_jk w_k. */
static qaws_status rationalize(qaws_diff_map* m, qaws_curve const* const* in, unsigned int n_in, qaws_curve* const* out)
{
	qaws_diff_map_entry* old = m->entries;
	unsigned int old_count = m->entry_count, i;
	qaws_status st = QAWS_STATUS_OK;
	qaws_scalar *w_in[MAP_MAX_OBJECTS], *p_in[MAP_MAX_OBJECTS], *w_out[MAP_MAX_OBJECTS], *p_out[MAP_MAX_OBJECTS];
	unsigned int n;

	memset(w_in, 0, sizeof(w_in));
	memset(p_in, 0, sizeof(p_in));
	memset(w_out, 0, sizeof(w_out));
	memset(p_out, 0, sizeof(p_out));
	for (i = 0; i < n_in && i < MAP_MAX_OBJECTS; i++)
		if (in[i])
		{
			read_values(in[i], QAWS_FIELD_WEIGHTS, &w_in[i], &n);
			read_values(in[i], QAWS_FIELD_CONTROL_POINTS, &p_in[i], &n);
		}
	for (i = 0; i < m->out_count && i < MAP_MAX_OBJECTS; i++)
	{
		read_values(out[i], QAWS_FIELD_WEIGHTS, &w_out[i], &n);
		read_values(out[i], QAWS_FIELD_CONTROL_POINTS, &p_out[i], &n);
	}

	m->entries = NULL;
	m->entry_count = m->capacity = 0;
	for (i = 0; i < old_count && st == QAWS_STATUS_OK; i++)
	{
		qaws_diff_map_entry const* e = &old[i];
		unsigned int dim = (unsigned int)out[e->out_object]->dimension, c;
		qaws_scalar A = e->weight;
		qaws_scalar wk = w_in[e->in_object] ? w_in[e->in_object][e->in_element] : QAWS_ONE;
		qaws_scalar wj = w_out[e->out_object] ? w_out[e->out_object][e->out_element] : QAWS_ONE;
		if (e->in_field != QAWS_FIELD_CONTROL_POINTS || e->out_field != QAWS_FIELD_CONTROL_POINTS)
		{
			st = map_add(m, e->out_object, (qaws_diff_field)e->out_field, e->out_element, e->out_component,
				e->in_object, (qaws_diff_field)e->in_field, e->in_element, e->in_component, A);
			continue;
		}
		st = map_add(m, e->out_object, QAWS_FIELD_CONTROL_POINTS, e->out_element, QAWS_DIFF_MAP_ALL_COMPONENTS,
			e->in_object, QAWS_FIELD_CONTROL_POINTS, e->in_element, QAWS_DIFF_MAP_ALL_COMPONENTS, A * wk / wj);
		if (st != QAWS_STATUS_OK || !w_in[e->in_object])
			continue;
		for (c = 0; c < dim && st == QAWS_STATUS_OK; c++)
		{
			qaws_scalar pk = p_in[e->in_object][(size_t)e->in_element * dim + c];
			qaws_scalar pj = p_out[e->out_object][(size_t)e->out_element * dim + c];
			st = map_add(m, e->out_object, QAWS_FIELD_CONTROL_POINTS, e->out_element, c,
				e->in_object, QAWS_FIELD_WEIGHTS, e->in_element, 0, A * (pk - pj) / wj);
		}
		if (st == QAWS_STATUS_OK && w_out[e->out_object])
			st = map_add(m, e->out_object, QAWS_FIELD_WEIGHTS, e->out_element, 0,
				e->in_object, QAWS_FIELD_WEIGHTS, e->in_element, 0, A);
	}
	qaws_internal_dealloc(NULL, old);
	for (i = 0; i < MAP_MAX_OBJECTS; i++)
	{
		qaws_internal_dealloc(NULL, w_in[i]);
		qaws_internal_dealloc(NULL, p_in[i]);
		qaws_internal_dealloc(NULL, w_out[i]);
		qaws_internal_dealloc(NULL, p_out[i]);
	}
	return st;
}

static qaws_status build_linear_map(
	qaws_curve const* const* inputs, unsigned int n_in,
	qaws_curve* const* outputs, unsigned int n_out,
	unsigned int extra_inputs,
	construct_fn fn, void const* user,
	qaws_diff_map** out_map)
{
	qaws_curve* zero[MAP_MAX_OBJECTS];
	qaws_curve* f0[MAP_MAX_OBJECTS];
	qaws_curve* fu[MAP_MAX_OBJECTS];
	qaws_curve const* probe_in[MAP_MAX_OBJECTS];
	field_info fields[MAP_MAX_OBJECTS][MAP_MAX_FIELDS];
	unsigned int nf[MAP_MAX_OBJECTS];
	int rational[MAP_MAX_OBJECTS], any_rational = 0;
	qaws_diff_map* m;
	unsigned int i, f, e;
	qaws_status st = QAWS_STATUS_OK;

	*out_map = NULL;
	if (n_in > MAP_MAX_OBJECTS || n_out > MAP_MAX_OBJECTS)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	memset(zero, 0, sizeof(zero));
	memset(f0, 0, sizeof(f0));
	memset(fu, 0, sizeof(fu));

	m = map_create(QAWS_DIFF_MAP_LINEAR_SPARSE, n_in + extra_inputs, n_out);
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	for (i = 0; i < n_in; i++)
	{
		nf[i] = vector_fields(inputs[i], fields[i], &rational[i]);
		any_rational |= rational[i];
		st = make_probe(inputs[i], fields[i], nf[i], rational[i], -1, 0, &zero[i]);
		if (st != QAWS_STATUS_OK)
			goto fail;
		probe_in[i] = zero[i];
	}
	st = fn(user, probe_in, f0);
	if (st != QAWS_STATUS_OK)
		goto fail;

	for (i = 0; i < n_in && st == QAWS_STATUS_OK; i++)
		for (f = 0; f < nf[i] && st == QAWS_STATUS_OK; f++)
			for (e = 0; e < fields[i][f].count && st == QAWS_STATUS_OK; e++)
			{
				qaws_curve* unit = NULL;
				st = make_probe(inputs[i], fields[i], nf[i], rational[i], (int)f, e, &unit);
				if (st != QAWS_STATUS_OK)
					break;
				probe_in[i] = unit;
				st = fn(user, probe_in, fu);
				if (st == QAWS_STATUS_OK)
					st = add_column(m, f0, fu, n_out, i, fields[i][f].field, e);
				destroy_curves(fu, n_out);
				probe_in[i] = zero[i];
				qaws_curve_destroy(unit);
			}
	if (st == QAWS_STATUS_OK)
	{
		int out_rational = 0;
		for (i = 0; i < n_out; i++)
		{
			field_info tmp[MAP_MAX_FIELDS];
			int r;
			vector_fields(outputs[i], tmp, &r);
			out_rational |= r;
		}
		if (any_rational || out_rational)
			st = rationalize(m, inputs, n_in, outputs);
	}
	if (st != QAWS_STATUS_OK)
		goto fail;

	destroy_curves(f0, n_out);
	destroy_curves(zero, n_in);
	*out_map = m;
	return QAWS_STATUS_OK;

fail:
	destroy_curves(f0, n_out);
	destroy_curves(fu, n_out);
	destroy_curves(zero, n_in);
	qaws_diff_map_destroy(m);
	return st;
}

/* ================================================================== */
/*  Operations                                                        */
/* ================================================================== */

static qaws_status construct_split(void const* user, qaws_curve const* const* in, qaws_curve** out)
{
	return qaws_curve_split(in[0], *(qaws_scalar const*)user, &out[0], &out[1]);
}

static qaws_status construct_join(void const* user, qaws_curve const* const* in, qaws_curve** out)
{
	(void)user;
	return qaws_curve_join(in[0], in[1], &out[0]);
}

static qaws_status construct_hermite_to_bezier(void const* user, qaws_curve const* const* in, qaws_curve** out)
{
	return qaws_curve_convert_hermite_to_bezier(in[0], *(unsigned int const*)user, &out[0]);
}

static qaws_status construct_bezier_to_bspline(void const* user, qaws_curve const* const* in, qaws_curve** out)
{
	(void)user;
	return qaws_curve_convert_bezier_to_bspline(in[0], &out[0]);
}

static qaws_status construct_bspline_to_nurbs(void const* user, qaws_curve const* const* in, qaws_curve** out)
{
	(void)user;
	return qaws_curve_convert_bspline_to_nurbs(in[0], &out[0]);
}

static qaws_status construct_elevate(void const* user, qaws_curve const* const* in, qaws_curve** out)
{
	(void)user;
	return qaws_curve_elevate_degree(in[0], &out[0]);
}

static qaws_status construct_reduce(void const* user, qaws_curve const* const* in, qaws_curve** out)
{
	(void)user;
	return qaws_curve_reduce_degree(in[0], &out[0]);
}

/* Runs the construction on the real inputs, then builds its map. */
static qaws_status run_with_map(
	qaws_curve const* const* inputs, unsigned int n_in, unsigned int extra_inputs,
	qaws_curve** outputs, unsigned int n_out,
	construct_fn fn, void const* user, qaws_diff_map** out_map)
{
	qaws_status st;
	unsigned int i;
	for (i = 0; i < n_out; i++)
		outputs[i] = NULL;
	st = fn(user, inputs, outputs);
	if (st != QAWS_STATUS_OK || !out_map)
		return st;
	st = build_linear_map(inputs, n_in, outputs, n_out, extra_inputs, fn, user, out_map);
	if (st != QAWS_STATUS_OK)
		destroy_curves(outputs, n_out);
	return st;
}

/* Bezier split: left_i = sum_{j<=i} B^i_j(t) P_j, right_i = sum_{j>=i} B^{n-i}_{j-i}(t) P_j.
   Adds the split parameter column (input 1). */
static qaws_status bezier_split_parameter_column(qaws_diff_map* m, qaws_curve const* curve, qaws_scalar t)
{
	unsigned int n = curve->degree, dim = (unsigned int)curve->dimension, i, j, c;
	qaws_scalar b[2 * QAWS_DIFF_MAX_SUPPORT];
	qaws_scalar* p = NULL;
	unsigned int count = 0;
	qaws_status st = read_values(curve, QAWS_FIELD_CONTROL_POINTS, &p, &count);
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i <= n && st == QAWS_STATUS_OK; i++)
	{
		qaws_scalar dl[3] = { 0, 0, 0 }, dr[3] = { 0, 0, 0 };
		/* left_i uses degree-i Bernstein in t */
		qaws_internal_bernstein_derivs(i, t, 1, b);
		for (j = 0; j <= i; j++)
			for (c = 0; c < dim; c++)
				dl[c] += b[(i + 1) + j] * p[j * dim + c];
		/* right_i uses degree-(n-i) Bernstein in t over P_i..P_n */
		qaws_internal_bernstein_derivs(n - i, t, 1, b);
		for (j = i; j <= n; j++)
			for (c = 0; c < dim; c++)
				dr[c] += b[(n - i + 1) + (j - i)] * p[j * dim + c];
		for (c = 0; c < dim && st == QAWS_STATUS_OK; c++)
		{
			st = map_add(m, 0, QAWS_FIELD_CONTROL_POINTS, i, c, 1, QAWS_FIELD_PARAMETER, 0, 0, dl[c]);
			if (st == QAWS_STATUS_OK)
				st = map_add(m, 1, QAWS_FIELD_CONTROL_POINTS, i, c, 1, QAWS_FIELD_PARAMETER, 0, 0, dr[c]);
		}
	}
	qaws_internal_dealloc(NULL, p);
	return st;
}

qaws_status qaws_curve_split_diff(qaws_curve const* curve, qaws_scalar parameter,
	qaws_curve** out_left, qaws_curve** out_right, qaws_diff_map** out_map)
{
	qaws_curve const* in[1];
	qaws_curve* out[2];
	qaws_status st;
	if (!curve || !out_left || !out_right)
		return QAWS_STATUS_INVALID_ARGUMENT;
	in[0] = curve;
	st = run_with_map(in, 1, 1, out, 2, construct_split, &parameter, out_map);
	if (st != QAWS_STATUS_OK)
		return st;
	if (out_map)
	{
		if (curve->kind == QAWS_CURVE_KIND_BEZIER)
			st = bezier_split_parameter_column(*out_map, curve, parameter);
		else
			(*out_map)->frozen_inputs |= 1u << 1;
		if (st != QAWS_STATUS_OK)
		{
			qaws_diff_map_destroy(*out_map);
			*out_map = NULL;
			destroy_curves(out, 2);
			return st;
		}
	}
	*out_left = out[0];
	*out_right = out[1];
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_join_diff(qaws_curve const* curve_a, qaws_curve const* curve_b,
	qaws_curve** out_joined, qaws_diff_map** out_map)
{
	qaws_curve const* in[2];
	if (!curve_a || !curve_b || !out_joined)
		return QAWS_STATUS_INVALID_ARGUMENT;
	in[0] = curve_a;
	in[1] = curve_b;
	return run_with_map(in, 2, 0, out_joined, 1, construct_join, NULL, out_map);
}

qaws_status qaws_curve_convert_hermite_to_bezier_diff(qaws_curve const* curve, unsigned int span_index,
	qaws_curve** out_bezier, qaws_diff_map** out_map)
{
	if (!curve || !out_bezier)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return run_with_map(&curve, 1, 0, out_bezier, 1, construct_hermite_to_bezier, &span_index, out_map);
}

qaws_status qaws_curve_convert_bezier_to_bspline_diff(qaws_curve const* curve,
	qaws_curve** out_bspline, qaws_diff_map** out_map)
{
	if (!curve || !out_bspline)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return run_with_map(&curve, 1, 0, out_bspline, 1, construct_bezier_to_bspline, NULL, out_map);
}

qaws_status qaws_curve_convert_bspline_to_nurbs_diff(qaws_curve const* curve,
	qaws_curve** out_nurbs, qaws_diff_map** out_map)
{
	if (!curve || !out_nurbs)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return run_with_map(&curve, 1, 0, out_nurbs, 1, construct_bspline_to_nurbs, NULL, out_map);
}

qaws_status qaws_curve_elevate_degree_diff(qaws_curve const* curve,
	qaws_curve** out_elevated, qaws_diff_map** out_map)
{
	if (!curve || !out_elevated)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return run_with_map(&curve, 1, 0, out_elevated, 1, construct_elevate, NULL, out_map);
}

qaws_status qaws_curve_reduce_degree_diff(qaws_curve const* curve,
	qaws_curve** out_reduced, qaws_diff_map** out_map)
{
	if (!curve || !out_reduced)
		return QAWS_STATUS_INVALID_ARGUMENT;
	return run_with_map(&curve, 1, 0, out_reduced, 1, construct_reduce, NULL, out_map);
}

/* ================================================================== */
/*  Least-squares B-spline fit                                         */
/*                                                                    */
/*  The fit solves A Q = r with A = N^T N over interior samples and    */
/*  r = N^T (D - endpoint terms). One forward sweep per input          */
/*  direction differentiates parameters, averaged knots, the basis     */
/*  (Cox-de Boor in dual numbers) and the solve: A dQ = dr - dA Q,     */
/*  with A factored once. Spans stay at their primal values.           */
/* ================================================================== */

typedef struct fit_dual
{
	double v, d;
} fit_dual;

/* Basis values and derivatives N[0..p] on a frozen span (Piegl A2.2). */
static void fit_basis_dual(fit_dual const* knots, unsigned int p, unsigned int span, fit_dual t, fit_dual* N)
{
	fit_dual left[QAWS_DIFF_MAX_SUPPORT + 1], right[QAWS_DIFF_MAX_SUPPORT + 1];
	unsigned int j, r;
	N[0].v = 1;
	N[0].d = 0;
	for (j = 1; j <= p; j++)
	{
		fit_dual saved = { 0, 0 };
		left[j].v = t.v - knots[span + 1 - j].v;
		left[j].d = t.d - knots[span + 1 - j].d;
		right[j].v = knots[span + j].v - t.v;
		right[j].d = knots[span + j].d - t.d;
		for (r = 0; r < j; r++)
		{
			fit_dual den, tmp;
			den.v = right[r + 1].v + left[j - r].v;
			den.d = right[r + 1].d + left[j - r].d;
			tmp.v = N[r].v / den.v;
			tmp.d = (N[r].d - tmp.v * den.d) / den.v;
			N[r].v = saved.v + right[r + 1].v * tmp.v;
			N[r].d = saved.d + right[r + 1].d * tmp.v + right[r + 1].v * tmp.d;
			saved.v = left[j - r].v * tmp.v;
			saved.d = left[j - r].d * tmp.v + left[j - r].v * tmp.d;
		}
		N[j] = saved;
	}
}

typedef struct fit_state
{
	unsigned int m, n, p, dim, sys, kc;
	int chord;
	double* data;              /* m * dim */
	qaws_scalar* params;       /* m, as the fit computes them */
	qaws_scalar* knots;        /* kc */
	unsigned int* spans;       /* m */
	unsigned int* knot_src;    /* interior knot p + j averages samples knot_src[j] - 1, knot_src[j] */
	double* knot_alpha;        /* kc */
	double* chol;              /* sys * sys lower factor of A */
	double* q;                 /* n * dim fitted control points */
	/* per-direction scratch */
	double *dt, *cum, *dknots, *dr, *dq;
	fit_dual* kd;
} fit_state;

/* Directional derivative of the control points (s->dq) and knots
   (s->dknots) for a data direction dD (m * dim, or NULL) and a parameter
   direction dt_in (m, or NULL; only used with given parameters). */
static void fit_directional(fit_state* s, double const* dD, double const* dt_in)
{
	unsigned int m = s->m, n = s->n, p = s->p, dim = s->dim, sys = s->sys;
	unsigned int i, j, k, d;
	double* dt = s->dt;

	/* Parameters. */
	for (k = 0; k < m; k++)
		dt[k] = dt_in && !s->chord ? dt_in[k] : 0.0;
	if (s->chord && dD)
	{
		double total = 0, dtotal = 0;
		s->cum[0] = 0;
		for (k = 1; k < m; k++)
		{
			double len2 = 0, dlen2 = 0, len;
			for (d = 0; d < dim; d++)
			{
				double a = s->data[k * dim + d] - s->data[(k - 1) * dim + d];
				double da = dD[k * dim + d] - dD[(k - 1) * dim + d];
				len2 += a * a;
				dlen2 += a * da;
			}
			len = sqrt(len2);
			total += len;
			dtotal += len > 0 ? dlen2 / len : 0.0;
			s->cum[k] = dtotal;
		}
		if (total > 0)
			for (k = 1; k + 1 < m; k++)
				dt[k] = (s->cum[k] - (double)s->params[k] * dtotal) / total;
		dt[0] = 0;
		dt[m - 1] = 0;
	}

	/* Averaged knots. */
	for (j = 0; j < s->kc; j++)
		s->dknots[j] = 0;
	for (j = 1; j + p < n; j++)
		s->dknots[p + j] = (1.0 - s->knot_alpha[j]) * dt[s->knot_src[j] - 1] + s->knot_alpha[j] * dt[s->knot_src[j]];
	for (j = 0; j < s->kc; j++)
	{
		s->kd[j].v = (double)s->knots[j];
		s->kd[j].d = s->dknots[j];
	}

	/* dr - dA Q over interior samples: sum_k dN_j (R_k - N Q) + N_j (dR_k - dN Q). */
	for (i = 0; i < sys * dim; i++)
		s->dr[i] = 0;
	for (k = 1; k + 1 < m; k++)
	{
		fit_dual N[QAWS_DIFF_MAX_SUPPORT + 1], t;
		unsigned int span = s->spans[k];
		t.v = (double)s->params[k];
		t.d = dt[k];
		fit_basis_dual(s->kd, p, span, t, N);
		for (d = 0; d < dim; d++)
		{
			double R = s->data[k * dim + d], dR = dD ? dD[k * dim + d] : 0.0;
			double nq = 0, dnq = 0;
			for (j = 0; j <= p; j++)
			{
				int cp = (int)span - (int)p + (int)j;
				if (cp == 0)
				{
					R -= N[j].v * s->data[d];
					dR -= N[j].d * s->data[d] + N[j].v * (dD ? dD[d] : 0.0);
				}
				else if (cp == (int)(n - 1))
				{
					R -= N[j].v * s->data[(m - 1) * dim + d];
					dR -= N[j].d * s->data[(m - 1) * dim + d] + N[j].v * (dD ? dD[(m - 1) * dim + d] : 0.0);
				}
				else
				{
					nq += N[j].v * s->q[(unsigned int)cp * dim + d];
					dnq += N[j].d * s->q[(unsigned int)cp * dim + d];
				}
			}
			for (j = 0; j <= p; j++)
			{
				int cp = (int)span - (int)p + (int)j;
				if (cp < 1 || cp > (int)(n - 2))
					continue;
				s->dr[d * sys + (unsigned int)(cp - 1)] += N[j].d * (R - nq) + N[j].v * (dR - dnq);
			}
		}
	}

	/* L L^T x = dr. */
	for (d = 0; d < dim; d++)
	{
		double* b = &s->dr[d * sys];
		for (i = 0; i < sys; i++)
		{
			double v = b[i];
			for (j = 0; j < i; j++)
				v -= s->chol[i * sys + j] * b[j];
			b[i] = v / s->chol[i * sys + i];
		}
		for (i = sys; i-- > 0;)
		{
			double v = b[i];
			for (j = i + 1; j < sys; j++)
				v -= s->chol[j * sys + i] * b[j];
			b[i] = v / s->chol[i * sys + i];
		}
	}
	for (d = 0; d < dim; d++)
	{
		s->dq[d] = dD ? dD[d] : 0.0;
		s->dq[(n - 1) * dim + d] = dD ? dD[(m - 1) * dim + d] : 0.0;
		for (i = 0; i < sys; i++)
			s->dq[(i + 1) * dim + d] = s->dr[d * sys + i];
	}
}

static qaws_status fit_emit(qaws_diff_map* m, fit_state const* s, unsigned int in_obj, qaws_diff_field in_field,
	unsigned int in_elem, unsigned int in_comp)
{
	unsigned int j, d;
	qaws_status st = QAWS_STATUS_OK;
	for (j = 0; j < s->n && st == QAWS_STATUS_OK; j++)
		for (d = 0; d < s->dim && st == QAWS_STATUS_OK; d++)
			if (s->dq[j * s->dim + d] != 0.0)
				st = map_add(m, 0, QAWS_FIELD_CONTROL_POINTS, j, d, in_obj, in_field, in_elem, in_comp,
					(qaws_scalar)s->dq[j * s->dim + d]);
	for (j = 0; j < s->kc && st == QAWS_STATUS_OK; j++)
		if (s->dknots[j] != 0.0)
			st = map_add(m, 0, QAWS_FIELD_KNOTS, j, 0, in_obj, in_field, in_elem, in_comp, (qaws_scalar)s->dknots[j]);
	return st;
}

/* Parameters, knots, spans and the Cholesky factor of A at the fit. */
static qaws_status fit_prepare(fit_state* s, struct qaws_bspline_fit_desc const* desc, qaws_curve const* curve,
	qaws_scalar* cps, double* A)
{
	unsigned int i, j, k, d, got;
	qaws_scalar const* dp = (qaws_scalar const*)desc->data_points;
	qaws_status st;

	for (i = 0; i < s->m * s->dim; i++)
		s->data[i] = (double)dp[i];
	if (desc->parameters)
		memcpy(s->params, desc->parameters, sizeof(qaws_scalar) * s->m);
	else
	{
		qaws_scalar total = 0;
		s->params[0] = 0;
		for (k = 1; k < s->m; k++)
		{
			qaws_scalar l2 = 0;
			for (d = 0; d < s->dim; d++)
			{
				qaws_scalar a = dp[k * s->dim + d] - dp[(k - 1) * s->dim + d];
				l2 += a * a;
			}
			total += QAWS_SQRT(l2);
			s->params[k] = total;
		}
		if (total > 0)
			for (k = 1; k < s->m; k++)
				s->params[k] /= total;
		s->params[s->m - 1] = QAWS_ONE;
	}
	st = qaws_curve_read_field(curve, QAWS_FIELD_KNOTS, s->knots, s->kc, &got);
	if (st == QAWS_STATUS_OK)
		st = qaws_curve_read_field(curve, QAWS_FIELD_CONTROL_POINTS, cps, s->n * s->dim, &got);
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < s->n * s->dim; i++)
		s->q[i] = (double)cps[i];
	{
		qaws_scalar dfl = (qaws_scalar)(s->m - 1) / (qaws_scalar)(s->n - s->p);
		for (j = 1; j + s->p < s->n; j++)
		{
			qaws_scalar jd = (qaws_scalar)j * dfl;
			unsigned int ii = (unsigned int)jd;
			s->knot_src[j] = ii;
			s->knot_alpha[j] = (double)(jd - (qaws_scalar)ii);
		}
	}
	for (k = 0; k < s->m; k++)
		s->spans[k] = qaws_internal_find_knot_span(s->knots, s->kc, s->p, s->n, s->params[k]);

	for (j = 0; j < s->kc; j++)
	{
		s->kd[j].v = (double)s->knots[j];
		s->kd[j].d = 0;
	}
	for (i = 0; i < s->sys * s->sys; i++)
		A[i] = 0;
	for (k = 1; k + 1 < s->m; k++)
	{
		fit_dual N[QAWS_DIFF_MAX_SUPPORT + 1], t;
		unsigned int a, b;
		t.v = (double)s->params[k];
		t.d = 0;
		fit_basis_dual(s->kd, s->p, s->spans[k], t, N);
		for (a = 0; a <= s->p; a++)
		{
			int ca = (int)s->spans[k] - (int)s->p + (int)a;
			if (ca < 1 || ca > (int)(s->n - 2))
				continue;
			for (b = 0; b <= s->p; b++)
			{
				int cb = (int)s->spans[k] - (int)s->p + (int)b;
				if (cb >= 1 && cb <= (int)(s->n - 2))
					A[(unsigned int)(ca - 1) * s->sys + (unsigned int)(cb - 1)] += N[a].v * N[b].v;
			}
		}
	}
	for (i = 0; i < s->sys; i++)
		for (j = 0; j <= i; j++)
		{
			double v = A[i * s->sys + j];
			for (k = 0; k < j; k++)
				v -= s->chol[i * s->sys + k] * s->chol[j * s->sys + k];
			if (i == j)
			{
				if (v <= 0)
					return QAWS_STATUS_NUMERICAL_FAILURE;
				s->chol[i * s->sys + i] = sqrt(v);
			}
			else
				s->chol[i * s->sys + j] = v / s->chol[j * s->sys + j];
		}
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_fit_bspline_diff(struct qaws_bspline_fit_desc const* desc, qaws_curve** out_curve,
	qaws_diff_map** out_map)
{
	fit_state s;
	qaws_diff_map* map = NULL;
	double *work = NULL, *A, *dD, *dtin;
	qaws_scalar *scal = NULL, *cps;
	unsigned int i, k, d;
	size_t nd;
	qaws_status st;

	if (!desc || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (out_map)
		*out_map = NULL;
	st = qaws_curve_fit_bspline(desc, out_curve);
	if (st != QAWS_STATUS_OK || !out_map)
		return st;

	memset(&s, 0, sizeof(s));
	s.m = desc->data_point_count;
	s.n = desc->control_point_count;
	s.p = desc->degree;
	s.dim = desc->dimension == QAWS_DIMENSION_2D ? 2u : 3u;
	s.sys = s.n - 2;
	s.kc = s.n + s.p + 1;
	s.chord = desc->parameters == NULL;
	if (s.p + 1 > QAWS_DIFF_MAX_SUPPORT)
	{
		st = QAWS_STATUS_UNSUPPORTED_OPERATION;
		goto done;
	}

	/* doubles: data, dD (m*dim each), dtin, dt, cum (m each), knot_alpha,
	   dknots (kc each), chol, A (sys^2 each), q, dq (n*dim each), dr (sys*dim),
	   kd (2 * kc as fit_dual) */
	nd = (size_t)s.m * s.dim * 2 + (size_t)s.m * 3 + (size_t)s.kc * 4 + (size_t)s.sys * s.sys * 2 +
		(size_t)s.n * s.dim * 2 + (size_t)s.sys * s.dim + 8;
	work = (double*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(double) * nd));
	scal = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (s.m + s.kc + s.n * s.dim + 4) +
		sizeof(unsigned int) * (s.m + s.kc + 4)));
	if (!work || !scal)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	s.data = work;
	dD = s.data + (size_t)s.m * s.dim;
	dtin = dD + (size_t)s.m * s.dim;
	s.dt = dtin + s.m;
	s.cum = s.dt + s.m;
	s.knot_alpha = s.cum + s.m;
	s.dknots = s.knot_alpha + s.kc;
	s.kd = (fit_dual*)(s.dknots + s.kc);
	s.chol = (double*)(s.kd + s.kc);
	A = s.chol + (size_t)s.sys * s.sys;
	s.q = A + (size_t)s.sys * s.sys;
	s.dq = s.q + (size_t)s.n * s.dim;
	s.dr = s.dq + (size_t)s.n * s.dim;
	s.params = scal;
	s.knots = s.params + s.m;
	cps = s.knots + s.kc;
	s.spans = (unsigned int*)(cps + s.n * s.dim + 4);
	s.knot_src = s.spans + s.m;

	st = fit_prepare(&s, desc, *out_curve, cps, A);
	if (st != QAWS_STATUS_OK)
		goto done;

	map = map_create(QAWS_DIFF_MAP_LINEAR_SPARSE, desc->parameters ? 2u : 1u, 1u);
	if (!map)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	for (i = 0; i < s.m * s.dim; i++)
		dD[i] = 0;
	for (i = 0; i < s.m; i++)
		dtin[i] = 0;
	for (k = 0; k < s.m && st == QAWS_STATUS_OK; k++)
		for (d = 0; d < s.dim && st == QAWS_STATUS_OK; d++)
		{
			dD[k * s.dim + d] = 1;
			fit_directional(&s, dD, NULL);
			dD[k * s.dim + d] = 0;
			st = fit_emit(map, &s, 0, QAWS_FIELD_POINTS, k, d);
		}
	if (desc->parameters)
		for (k = 0; k < s.m && st == QAWS_STATUS_OK; k++)
		{
			dtin[k] = 1;
			fit_directional(&s, NULL, dtin);
			dtin[k] = 0;
			st = fit_emit(map, &s, 1, QAWS_FIELD_PARAMETER, k, 0);
		}
	if (st == QAWS_STATUS_OK)
	{
		*out_map = map;
		map = NULL;
	}
done:
	qaws_diff_map_destroy(map);
	qaws_internal_dealloc(NULL, work);
	qaws_internal_dealloc(NULL, scal);
	if (st != QAWS_STATUS_OK && *out_curve)
	{
		qaws_curve_destroy(*out_curve);
		*out_curve = NULL;
	}
	return st;
}

/* ================================================================== */
/*  3D offset                                                          */
/*                                                                    */
/*  The offset samples Q_i = C(t_i) + d n(t_i) at fixed parameters and */
/*  fits a B-spline through them. The fit is linear in the samples, so */
/*  each column of the map is the fit of the sample tangents dQ_i.     */
/* ================================================================== */

typedef struct offset_state
{
	qaws_curve const* curve;
	qaws_scalar distance;
	int mode;
	qaws_vec3 dir;             /* mode 0: unit direction */
	qaws_scalar dir_len;
	unsigned int n_samples, n_cp;
	qaws_scalar* params;       /* n_samples */
	qaws_curve_jet_3d* jets;   /* primal jets at the samples */
	qaws_vec3* normals;        /* unit offset directions */
	qaws_curve_jet_3d* tan;    /* scratch: tangent jets */
	qaws_scalar* dq;           /* scratch: n_samples * 3 */
	qaws_scalar* cps;          /* scratch: n_cp * 3 */
} offset_state;

/* Fits dq and adds the control point column for one input scalar. */
static qaws_status offset_emit(qaws_diff_map* m, offset_state* s, unsigned int in_obj, qaws_diff_field in_field,
	unsigned int in_elem, unsigned int in_comp)
{
	qaws_curve* fit = NULL;
	unsigned int j, c, got = 0;
	qaws_status st = qaws_internal_fit_bspline(QAWS_DIMENSION_3D, 3, s->params, s->dq, s->n_samples, s->n_cp, &fit);
	if (st != QAWS_STATUS_OK)
		return st;
	st = qaws_curve_read_field(fit, QAWS_FIELD_CONTROL_POINTS, s->cps, s->n_cp * 3, &got);
	qaws_curve_destroy(fit);
	if (st != QAWS_STATUS_OK)
		return st;
	for (j = 0; j < s->n_cp && st == QAWS_STATUS_OK; j++)
		for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
			if (s->cps[j * 3 + c] != QAWS_ZERO)
				st = map_add(m, 0, QAWS_FIELD_CONTROL_POINTS, j, c, in_obj, in_field, in_elem, in_comp, s->cps[j * 3 + c]);
	return st;
}

/* dq for a tangent of the input curve's fields. */
static qaws_status offset_curve_column(offset_state* s, qaws_diff_views const* views)
{
	unsigned int i;
	qaws_status st = qaws_curve_eval_batch_tangent_3d(NULL, s->curve, s->params, NULL, s->n_samples,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, views, NULL, s->tan);
	if (st != QAWS_STATUS_OK)
		return st;
	for (i = 0; i < s->n_samples; i++)
	{
		qaws_vec3 dq = s->tan[i].d[0];
		if (s->mode == 1)
		{
			qaws_curve_geometry_3d g, dg;
			qaws_diff_validity validity;
			qaws_curve_geometry_eval_3d(&s->jets[i], &s->tan[i], NULL, &g, &dg, NULL, &validity);
			dq = qaws_v3_axpy(dq, dg.normal, s->distance);
		}
		s->dq[i * 3 + 0] = dq.x;
		s->dq[i * 3 + 1] = dq.y;
		s->dq[i * 3 + 2] = dq.z;
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_curve_offset_3d_diff(
	qaws_curve const* curve,
	qaws_scalar distance,
	int direction_mode,
	qaws_vec3 const* direction,
	unsigned int sample_count,
	qaws_curve** out_curve,
	qaws_diff_map** out_map)
{
	offset_state s;
	qaws_diff_map* map = NULL;
	qaws_field_desc fields[MAP_MAX_FIELDS];
	qaws_range range;
	unsigned int nf = 0, f, e, c, i;
	qaws_status st;

	if (!curve || !out_curve)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (direction_mode != 0 && direction_mode != 1)
		return QAWS_STATUS_UNSUPPORTED_OPERATION;
	if (out_map)
		*out_map = NULL;
	st = qaws_curve_offset_3d(curve, distance, direction_mode, direction, NULL, sample_count, out_curve);
	if (st != QAWS_STATUS_OK || !out_map)
		return st;

	memset(&s, 0, sizeof(s));
	s.curve = curve;
	s.distance = distance;
	s.mode = direction_mode;
	s.n_samples = sample_count > 0 ? sample_count : 256;
	s.n_cp = s.n_samples / 4;
	if (s.n_cp < 8) s.n_cp = 8;
	if (s.n_cp > s.n_samples) s.n_cp = s.n_samples;
	if (direction_mode == 0)
	{
		s.dir = *direction;
		s.dir_len = QAWS_SQRT(qaws_v3_dot(s.dir, s.dir));
		if (s.dir_len > QAWS_LITERAL(1e-12))
			s.dir = qaws_v3_scale(s.dir, QAWS_ONE / s.dir_len);
	}
	s.params = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (s.n_samples * 4 + s.n_cp * 3 + 4)));
	s.jets = (qaws_curve_jet_3d*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_curve_jet_3d) * s.n_samples * 2));
	s.normals = (qaws_vec3*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_vec3) * s.n_samples));
	if (!s.params || !s.jets || !s.normals)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	s.dq = s.params + s.n_samples;
	s.cps = s.dq + s.n_samples * 3;
	s.tan = s.jets + s.n_samples;

	/* Samples and directions exactly as the offset computes them. */
	range = curve->parameter_range;
	for (i = 0; i < s.n_samples; i++)
		s.params[i] = range.min_value + (range.max_value - range.min_value) * (qaws_scalar)i / (qaws_scalar)(s.n_samples - 1);
	st = qaws_curve_eval_batch_tangent_3d(NULL, curve, s.params, NULL, s.n_samples,
		QAWS_EVAL_FLAG_POSITION | QAWS_EVAL_FLAG_D1 | QAWS_EVAL_FLAG_D2 | QAWS_EVAL_FLAG_D3, NULL, s.jets, s.tan);
	if (st != QAWS_STATUS_OK)
		goto done;
	for (i = 0; i < s.n_samples; i++)
	{
		if (direction_mode == 1)
		{
			qaws_curve_geometry_3d g;
			qaws_diff_validity validity;
			qaws_curve_geometry_eval_3d(&s.jets[i], NULL, NULL, &g, NULL, NULL, &validity);
			if (validity == QAWS_DIFF_INVALID || validity == QAWS_DIFF_ILL_CONDITIONED)
			{
				/* straight pieces: the Frenet normal is undefined there */
				st = QAWS_STATUS_UNSUPPORTED_OPERATION;
				goto done;
			}
			s.normals[i] = g.normal;
		}
		else
			s.normals[i] = s.dir;
	}

	map = map_create(QAWS_DIFF_MAP_LINEAR_SPARSE, 3u, 1u);
	if (!map)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}

	/* input 0: every differentiable field of the curve */
	if (qaws_curve_describe_fields(curve, fields, MAP_MAX_FIELDS, &nf) != QAWS_STATUS_OK)
		nf = 0;
	for (f = 0; f < nf && f < MAP_MAX_FIELDS && st == QAWS_STATUS_OK; f++)
	{
		unsigned int comps = fields[f].value_type == QAWS_VALUE_SCALAR ? 1u : (unsigned int)fields[f].value_type;
		qaws_scalar* seed;
		if (!(fields[f].capabilities & QAWS_CAP_TANGENT))
			continue;
		seed = (qaws_scalar*)qaws_internal_alloc(NULL, (unsigned long)(sizeof(qaws_scalar) * (fields[f].count * comps + 1)));
		if (!seed)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			break;
		}
		memset(seed, 0, sizeof(qaws_scalar) * fields[f].count * comps);
		for (e = 0; e < fields[f].count && st == QAWS_STATUS_OK; e++)
			for (c = 0; c < comps && st == QAWS_STATUS_OK; c++)
			{
				qaws_field_view fv = qaws_field_view_make(fields[f].field, seed, fields[f].count, comps);
				qaws_diff_views views;
				views.fields = &fv;
				views.field_count = 1;
				views.children = NULL;
				views.child_count = 0;
				seed[e * comps + c] = QAWS_ONE;
				st = offset_curve_column(&s, &views);
				seed[e * comps + c] = QAWS_ZERO;
				if (st == QAWS_STATUS_OK)
					st = offset_emit(map, &s, 0, fields[f].field, e, c);
			}
		qaws_internal_dealloc(NULL, seed);
	}

	/* input 1: the distance */
	if (st == QAWS_STATUS_OK)
	{
		for (i = 0; i < s.n_samples; i++)
		{
			s.dq[i * 3 + 0] = s.normals[i].x;
			s.dq[i * 3 + 1] = s.normals[i].y;
			s.dq[i * 3 + 2] = s.normals[i].z;
		}
		st = offset_emit(map, &s, 1, QAWS_FIELD_PARAMETER, 0, 0);
	}

	/* input 2: the direction (mode 0): d (I - n n^T) du / |u| */
	if (direction_mode == 0)
		for (c = 0; c < 3 && st == QAWS_STATUS_OK; c++)
		{
			qaws_vec3 du = qaws_v3_zero(), dn;
			if (c == 0) du.x = QAWS_ONE; else if (c == 1) du.y = QAWS_ONE; else du.z = QAWS_ONE;
			dn = qaws_v3_scale(qaws_v3_axpy(du, s.dir, -qaws_v3_dot(s.dir, du)), QAWS_ONE / s.dir_len);
			for (i = 0; i < s.n_samples; i++)
			{
				s.dq[i * 3 + 0] = distance * dn.x;
				s.dq[i * 3 + 1] = distance * dn.y;
				s.dq[i * 3 + 2] = distance * dn.z;
			}
			st = offset_emit(map, &s, 2, QAWS_FIELD_DIRECTION, 0, c);
		}
	else
		map->frozen_inputs |= 1u << 2;

	if (st == QAWS_STATUS_OK)
	{
		*out_map = map;
		map = NULL;
	}
done:
	qaws_diff_map_destroy(map);
	qaws_internal_dealloc(NULL, s.params);
	qaws_internal_dealloc(NULL, s.jets);
	qaws_internal_dealloc(NULL, s.normals);
	if (st != QAWS_STATUS_OK && *out_curve)
	{
		qaws_curve_destroy(*out_curve);
		*out_curve = NULL;
	}
	return st;
}
