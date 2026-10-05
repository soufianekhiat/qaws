#include "qaws_diff_map.h"
#include "qaws_diff.h"
#include "qaws_curve.h"
#include "qaws_operations.h"
#include "qaws_convert.h"
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
