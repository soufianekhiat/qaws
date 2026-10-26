/*
 * Volume meshes: cells as runs of corners, facets matched pairwise at
 * creation (sorted by their vertex sets), the border read back as a surface.
 */

#include "qaws_mesh_internal.h"
#include "../internal/qaws_internal_memory.h"
#include "../qaws_volume_mesh.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

unsigned char const vm_corners[4] = { 4, 8, 6, 5 }, vm_facets[4] = { 4, 6, 5, 5 };

static unsigned char const vm_tet_f[4][4] = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
static unsigned char const vm_hex_f[6][4] = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
static unsigned char const vm_prism_f[5][4] = { { 0, 2, 1 }, { 3, 4, 5 }, { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 2, 0, 3, 5 } };
static unsigned char const vm_pyr_f[5][4] = { { 0, 3, 2, 1 }, { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } };
static unsigned char const vm_prism_n[5] = { 3, 3, 4, 4, 4 }, vm_pyr_n[5] = { 4, 3, 3, 3, 3 };

unsigned char const* vm_facet(unsigned int type, unsigned int f, unsigned int* n)
{
	switch (type)
	{
	case QAWS_CELL_TETRAHEDRON: *n = 3; return vm_tet_f[f];
	case QAWS_CELL_HEXAHEDRON: *n = 4; return vm_hex_f[f];
	case QAWS_CELL_PRISM: *n = vm_prism_n[f]; return vm_prism_f[f];
	default: *n = vm_pyr_n[f]; return vm_pyr_f[f];
	}
}

/* each type's edges, as corner pairs */
static unsigned char const vm_edges_tet[6][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 1, 3 }, { 2, 3 } };
static unsigned char const vm_edges_hex[12][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
static unsigned char const vm_edges_prism[9][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 3, 4 }, { 4, 5 }, { 5, 3 }, { 0, 3 }, { 1, 4 }, { 2, 5 } };
static unsigned char const vm_edges_pyr[8][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 4 }, { 1, 4 }, { 2, 4 }, { 3, 4 } };

static unsigned char const (*vm_edges(unsigned int type, unsigned int* n))[2]
{
	switch (type)
	{
	case QAWS_CELL_TETRAHEDRON: *n = 6; return vm_edges_tet;
	case QAWS_CELL_HEXAHEDRON: *n = 12; return vm_edges_hex;
	case QAWS_CELL_PRISM: *n = 9; return vm_edges_prism;
	default: *n = 8; return vm_edges_pyr;
	}
}

static double vm_det(double const* a, double const* b, double const* c, double const* d)
{
	double u[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] }, v[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
	double w[3] = { d[0] - a[0], d[1] - a[1], d[2] - a[2] };
	return u[0] * (v[1] * w[2] - v[2] * w[1]) + u[1] * (v[2] * w[0] - v[0] * w[2]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
}

/* ======================================================================== */
/*  Creation                                                                */
/* ======================================================================== */

typedef struct vm_key
{
	unsigned int v[4];              /* the facet's vertices sorted (HM_NONE padding) */
	unsigned int facet;
} vm_key;

static int vm_key_cmp(void const* x, void const* y)
{
	vm_key const *a = (vm_key const*)x, *b = (vm_key const*)y;
	int k;
	for (k = 0; k < 4; k++)
		if (a->v[k] != b->v[k])
			return a->v[k] < b->v[k] ? -1 : 1;
	return a->facet < b->facet ? -1 : a->facet > b->facet ? 1 : 0;
}

/* the cell owning facet f */
static unsigned int vm_facet_cell(qaws_volume_mesh const* m, unsigned int f)
{
	unsigned int lo = 0, hi = m->nc;
	while (hi - lo > 1)
	{
		unsigned int mid = (lo + hi) / 2;
		if (m->f0[mid] <= f)
			lo = mid;
		else
			hi = mid;
	}
	return lo;
}

/* facet f's vertices, outward; their count */
static unsigned int vm_facet_vertices(qaws_volume_mesh const* m, unsigned int f, unsigned int* out)
{
	unsigned int c = vm_facet_cell(m, f), n, k;
	unsigned char const* loc = vm_facet(m->type[c], f - m->f0[c], &n);
	for (k = 0; k < n; k++)
		out[k] = m->cv[m->c0[c] + loc[k]];
	return n;
}

/* do two facets on the same vertices turn opposite ways */
static int vm_opposite(unsigned int const* a, unsigned int const* b, unsigned int n)
{
	unsigned int k, i;
	for (k = 0; k < n && b[k] != a[0]; k++) {}
	if (k == n)
		return 0;
	for (i = 0; i < n; i++)
		if (b[(k + i) % n] != a[(n - i) % n])
			return 0;
	return 1;
}

void qaws_volume_mesh_destroy(qaws_volume_mesh* m)
{
	unsigned int a;
	if (!m)
		return;
	for (a = 0; a < m->na; a++)
		qaws_internal_free(m->attr[a].v);
	qaws_internal_free(m->attr);
	qaws_internal_free(m->x); qaws_internal_free(m->type); qaws_internal_free(m->c0); qaws_internal_free(m->cv); qaws_internal_free(m->f0); qaws_internal_free(m->fo);
	qaws_internal_free(m);
}

qaws_status vm_create(double const* x, unsigned int nv, unsigned int const* cv, unsigned char const* type, unsigned int nc, qaws_volume_mesh** out)
{
	qaws_volume_mesh* m;
	vm_key* keys = NULL;
	unsigned int c, k, f, ncorner = 0, nfacet = 0, i;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	*out = NULL;
	for (c = 0; c < nc; c++)
	{
		if (type[c] > QAWS_CELL_PYRAMID)
			return QAWS_STATUS_INVALID_ARGUMENT;
		for (k = 0; k < vm_corners[type[c]]; k++)
			if (cv[ncorner + k] >= nv)
				return QAWS_STATUS_INVALID_ARGUMENT;
		ncorner += vm_corners[type[c]];
		nfacet += vm_facets[type[c]];
	}
	m = (qaws_volume_mesh*)qaws_internal_zalloc(1, sizeof *m);
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m->nv = nv;
	m->nc = nc;
	m->x = (double*)qaws_internal_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)nv + 1)));
	m->type = (unsigned char*)qaws_internal_malloc((size_t)nc + 1);
	m->c0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nc + 1)));
	m->f0 = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nc + 1)));
	m->cv = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)ncorner + 1)));
	m->fo = (unsigned int*)qaws_internal_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)nfacet + 1)));
	keys = (vm_key*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(vm_key), ((size_t)nfacet + 1)));
	if (!m->x || !m->type || !m->c0 || !m->f0 || !m->cv || !m->fo || !keys)
		goto done;
	memcpy(m->x, x, sizeof(double) * 3 * nv);
	memcpy(m->type, type, nc);
	memcpy(m->cv, cv, sizeof(unsigned int) * ncorner);
	m->c0[0] = m->f0[0] = 0;
	for (c = 0; c < nc; c++)
	{
		m->c0[c + 1] = m->c0[c] + vm_corners[type[c]];
		m->f0[c + 1] = m->f0[c] + vm_facets[type[c]];
	}
	/* facets with the same vertex set side by side */
	for (f = 0; f < nfacet; f++)
	{
		unsigned int v[4], n = vm_facet_vertices(m, f, v), j;
		for (k = n; k < 4; k++)
			v[k] = HM_NONE;
		for (k = 1; k < 4; k++)
			for (j = k; j > 0 && v[j - 1] > v[j]; j--)
			{
				unsigned int t = v[j];
				v[j] = v[j - 1];
				v[j - 1] = t;
			}
		memcpy(keys[f].v, v, sizeof v);
		keys[f].facet = f;
		m->fo[f] = HM_NONE;
	}
	qsort(keys, nfacet, sizeof(vm_key), vm_key_cmp);
	st = QAWS_STATUS_OK;
	for (i = 0; i < nfacet && st == QAWS_STATUS_OK; )
	{
		unsigned int j = i + 1;
		while (j < nfacet && !memcmp(keys[j].v, keys[i].v, sizeof keys[i].v))
			j++;
		if (j - i == 2)
		{
			unsigned int a[4], b[4], n = vm_facet_vertices(m, keys[i].facet, a);
			vm_facet_vertices(m, keys[i + 1].facet, b);
			if (!vm_opposite(a, b, n))
				st = QAWS_STATUS_INVALID_TOPOLOGY;
			m->fo[keys[i].facet] = keys[i + 1].facet;
			m->fo[keys[i + 1].facet] = keys[i].facet;
		}
		else if (j - i > 2)
			st = QAWS_STATUS_INVALID_TOPOLOGY;
		i = j;
	}
done:
	qaws_internal_free(keys);
	if (st != QAWS_STATUS_OK)
	{
		qaws_volume_mesh_destroy(m);
		return st;
	}
	*out = m;
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_create(qaws_volume_mesh** out_mesh, qaws_volume_mesh_desc const* desc)
{
	double* x;
	unsigned char* type;
	unsigned int i;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!desc || (desc->vertex_count && !desc->positions) || (desc->cell_count && !desc->cell_vertices))
		return QAWS_STATUS_INVALID_ARGUMENT;
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)desc->vertex_count + 1)));
	type = (unsigned char*)qaws_internal_scratch_malloc((size_t)desc->cell_count + 1);
	if (!x || !type)
	{
		qaws_internal_free(x);
		qaws_internal_free(type);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (i = 0; i < 3 * desc->vertex_count; i++)
		x[i] = (double)desc->positions[i];
	for (i = 0; i < desc->cell_count; i++)
		type[i] = (unsigned char)(desc->cell_types ? ((unsigned int)desc->cell_types[i] <= QAWS_CELL_PYRAMID ? desc->cell_types[i] : 255) : QAWS_CELL_TETRAHEDRON);
	st = vm_create(x, desc->vertex_count, desc->cell_vertices, type, desc->cell_count, out_mesh);
	qaws_internal_free(x);
	qaws_internal_free(type);
	return st;
}

qaws_status qaws_volume_mesh_clone(qaws_volume_mesh** out_mesh, qaws_volume_mesh const* m)
{
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!m)
		return QAWS_STATUS_INVALID_ARGUMENT;
	{
		qaws_status st = vm_create(m->x, m->nv, m->cv, m->type, m->nc, out_mesh);
		if (st == QAWS_STATUS_OK && (st = vm_attr_carry(*out_mesh, m, NULL, NULL, NULL, NULL)) != QAWS_STATUS_OK)
		{
			qaws_volume_mesh_destroy(*out_mesh);
			*out_mesh = NULL;
		}
		return st;
	}
}

/* ======================================================================== */
/*  Attributes                                                              */
/* ======================================================================== */

unsigned int vm_attr_size(qaws_volume_mesh const* m, unsigned int element)
{
	return element == QAWS_MESH_VERTEX ? m->nv : element == QAWS_MESH_CORNER ? m->c0[m->nc] : element == QAWS_MESH_CELL ? m->nc : m->f0[m->nc];
}

/* a new attribute of zeros: its index, HM_NONE out of memory */
static unsigned int vm_attr_add(qaws_volume_mesh* m, char const* name, unsigned int element, unsigned int rule, unsigned int n)
{
	vm_attr* g = (vm_attr*)qaws_internal_realloc(m->attr, sizeof(vm_attr) * ((size_t)m->na + 1));
	vm_attr* a;
	if (!g)
		return HM_NONE;
	m->attr = g;
	a = &m->attr[m->na];
	memset(a, 0, sizeof *a);
	a->v = (double*)qaws_internal_zalloc((size_t)n * vm_attr_size(m, element) + 1, sizeof(double));
	if (!a->v)
		return HM_NONE;
	if (name)
		strncpy(a->name, name, sizeof a->name - 1);
	a->element = element;
	a->rule = rule;
	a->n = n;
	return m->na++;
}

qaws_status vm_attr_carry(qaws_volume_mesh* dst, qaws_volume_mesh const* src, unsigned int const* vmap, unsigned int const* cmap,
	unsigned int const* cellmap, unsigned int const* fmap)
{
	unsigned int a, i;
	for (a = 0; a < src->na; a++)
	{
		vm_attr const* s = &src->attr[a];
		unsigned int const* map = s->element == QAWS_MESH_VERTEX ? vmap : s->element == QAWS_MESH_CORNER ? cmap : s->element == QAWS_MESH_CELL ? cellmap : fmap;
		unsigned int k = vm_attr_add(dst, s->name, s->element, s->rule, s->n), count;
		if (k == HM_NONE)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		count = vm_attr_size(dst, s->element);
		for (i = 0; i < count; i++)
		{
			unsigned int j = map ? map[i] : i;
			if (j != HM_NONE)
				memcpy(dst->attr[k].v + (size_t)s->n * i, s->v + (size_t)s->n * j, sizeof(double) * s->n);
		}
	}
	return QAWS_STATUS_OK;
}

int vm_attr_continuous(qaws_volume_mesh const* m, unsigned int f)
{
	unsigned int g = m->fo[f], c, d, a, k, i, j, v[4], nf;
	if (g == HM_NONE)
		return 1;
	c = vm_facet_cell(m, f);
	d = vm_facet_cell(m, g);
	nf = vm_facet_vertices(m, f, v);
	for (a = 0; a < m->na; a++)
	{
		vm_attr const* t = &m->attr[a];
		if (t->element == QAWS_MESH_CELL && memcmp(t->v + (size_t)t->n * c, t->v + (size_t)t->n * d, sizeof(double) * t->n))
			return 0;
		if (t->element == QAWS_MESH_FACET)
			for (k = 0; k < t->n; k++)
				if (t->v[(size_t)t->n * f + k] != 0 || t->v[(size_t)t->n * g + k] != 0)
					return 0;
		if (t->element == QAWS_MESH_CORNER)
			for (k = 0; k < nf; k++)
			{
				unsigned int ci = HM_NONE, di = HM_NONE;
				for (i = m->c0[c]; i < m->c0[c + 1]; i++)
					if (m->cv[i] == v[k])
						ci = i;
				for (j = m->c0[d]; j < m->c0[d + 1]; j++)
					if (m->cv[j] == v[k])
						di = j;
				if (memcmp(t->v + (size_t)t->n * ci, t->v + (size_t)t->n * di, sizeof(double) * t->n))
					return 0;
			}
	}
	return 1;
}

/* a facet's key: its vertices from the smallest, in their turn (the side kept) */
static void vm_facet_key(qaws_volume_mesh const* m, unsigned int f, unsigned int* key)
{
	unsigned int v[4], n = vm_facet_vertices(m, f, v), s = 0, k;
	for (k = 1; k < n; k++)
		if (v[k] < v[s])
			s = k;
	for (k = 0; k < 4; k++)
		key[k] = k < n ? v[(s + k) % n] : HM_NONE;
}

int vm_facet_map(unsigned int* out, qaws_volume_mesh const* dst, qaws_volume_mesh const* src)
{
	unsigned int nd = dst->f0[dst->nc], ns = src->f0[src->nc], f;
	vm_key* keys = (vm_key*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(vm_key), (size_t)ns + 1));
	if (!keys)
		return 0;
	for (f = 0; f < ns; f++)
	{
		vm_facet_key(src, f, keys[f].v);
		keys[f].facet = f;
	}
	qsort(keys, ns, sizeof(vm_key), vm_key_cmp);
	for (f = 0; f < nd; f++)
	{
		unsigned int key[4], lo = 0, hi = ns;
		vm_key probe;
		vm_facet_key(dst, f, key);
		memcpy(probe.v, key, sizeof key);
		probe.facet = 0;   /* (before every facet on these vertices) */
		out[f] = HM_NONE;
		while (lo < hi)
		{
			unsigned int mid = lo + (hi - lo) / 2;
			if (vm_key_cmp(&keys[mid], &probe) < 0)
				lo = mid + 1;
			else
				hi = mid;
		}
		if (lo < ns && !memcmp(keys[lo].v, key, sizeof key))
			out[f] = keys[lo].facet;
	}
	qaws_internal_free(keys);
	return 1;
}

qaws_status qaws_volume_mesh_add_attribute(unsigned int* out_attribute, qaws_volume_mesh* mesh, qaws_scalar const* values,
	qaws_attribute_desc const* desc)
{
	unsigned int a, i, count;
	if (out_attribute)
		*out_attribute = QAWS_MESH_NONE;
	if (!mesh || !desc || (desc->element != QAWS_MESH_VERTEX && desc->element != QAWS_MESH_CORNER && desc->element != QAWS_MESH_CELL
		&& desc->element != QAWS_MESH_FACET) || !hm_rule_valid((unsigned int)desc->rule, desc->components))
		return QAWS_STATUS_INVALID_ARGUMENT;
	a = vm_attr_add(mesh, desc->name, (unsigned int)desc->element, (unsigned int)desc->rule, desc->components);
	if (a == HM_NONE)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	count = desc->components * vm_attr_size(mesh, (unsigned int)desc->element);
	for (i = 0; values && i < count; i++)
		mesh->attr[a].v[i] = (double)values[i];
	if (out_attribute)
		*out_attribute = a;
	return QAWS_STATUS_OK;
}

unsigned int qaws_volume_mesh_get_attribute_count(qaws_volume_mesh const* mesh)
{
	return mesh ? mesh->na : 0;
}

qaws_status qaws_volume_mesh_find_attribute(unsigned int* out_attribute, qaws_volume_mesh const* mesh, char const* name)
{
	unsigned int a;
	if (!mesh || !name || !out_attribute)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_attribute = QAWS_MESH_NONE;
	for (a = 0; a < mesh->na; a++)
		if (!strncmp(mesh->attr[a].name, name, sizeof mesh->attr[a].name - 1))
		{
			*out_attribute = a;
			return QAWS_STATUS_OK;
		}
	return QAWS_STATUS_OUT_OF_RANGE;
}

qaws_status qaws_volume_mesh_get_attribute_desc(qaws_attribute_desc* out_desc, qaws_volume_mesh const* mesh, unsigned int attribute)
{
	vm_attr const* at;
	if (!mesh || attribute >= mesh->na || !out_desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	at = &mesh->attr[attribute];
	out_desc->name = at->name;
	out_desc->element = (qaws_mesh_element)at->element;
	out_desc->components = at->n;
	out_desc->rule = (qaws_attribute_rule)at->rule;
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_get_attribute(qaws_scalar* out_values, qaws_volume_mesh const* mesh, unsigned int attribute)
{
	unsigned int i, count;
	if (!mesh || attribute >= mesh->na || !out_values)
		return QAWS_STATUS_INVALID_ARGUMENT;
	count = mesh->attr[attribute].n * vm_attr_size(mesh, mesh->attr[attribute].element);
	for (i = 0; i < count; i++)
		out_values[i] = (qaws_scalar)mesh->attr[attribute].v[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_set_attribute(qaws_volume_mesh* mesh, unsigned int attribute, qaws_scalar const* values)
{
	unsigned int i, count;
	if (!mesh || attribute >= mesh->na || !values)
		return QAWS_STATUS_INVALID_ARGUMENT;
	count = mesh->attr[attribute].n * vm_attr_size(mesh, mesh->attr[attribute].element);
	for (i = 0; i < count; i++)
		mesh->attr[attribute].v[i] = (double)values[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_remove_attribute(qaws_volume_mesh* mesh, unsigned int attribute)
{
	if (!mesh || attribute >= mesh->na)
		return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_internal_free(mesh->attr[attribute].v);
	memmove(mesh->attr + attribute, mesh->attr + attribute + 1, sizeof(vm_attr) * (mesh->na - attribute - 1));
	mesh->na--;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Reading                                                                 */
/* ======================================================================== */

unsigned int qaws_volume_mesh_get_vertex_count(qaws_volume_mesh const* m)
{
	return m ? m->nv : 0;
}

unsigned int qaws_volume_mesh_get_cell_count(qaws_volume_mesh const* m)
{
	return m ? m->nc : 0;
}

unsigned int qaws_volume_mesh_get_corner_count(qaws_volume_mesh const* m)
{
	return m ? m->c0[m->nc] : 0;
}

unsigned int qaws_volume_mesh_get_facet_count(qaws_volume_mesh const* m)
{
	return m ? m->f0[m->nc] : 0;
}

qaws_status qaws_volume_mesh_get_positions(qaws_scalar* out, qaws_volume_mesh const* m)
{
	unsigned int i;
	if (!m || (m->nv && !out))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * m->nv; i++)
		out[i] = (qaws_scalar)m->x[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_get_cells(unsigned int* out_cell_vertices, qaws_cell_type* out_cell_types, qaws_volume_mesh const* m)
{
	unsigned int c;
	if (!m || (m->nc && !out_cell_vertices))
		return QAWS_STATUS_INVALID_ARGUMENT;
	memcpy(out_cell_vertices, m->cv, sizeof(unsigned int) * m->c0[m->nc]);
	for (c = 0; out_cell_types && c < m->nc; c++)
		out_cell_types[c] = (qaws_cell_type)m->type[c];
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_get_cell_adjacency(unsigned int* out_cells, qaws_volume_mesh const* m)
{
	unsigned int f;
	if (!m || (m->nc && !out_cells))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (f = 0; f < m->f0[m->nc]; f++)
		out_cells[f] = m->fo[f] == HM_NONE ? QAWS_MESH_NONE : vm_facet_cell(m, m->fo[f]);
	return QAWS_STATUS_OK;
}

qaws_status qaws_volume_mesh_get_boundary(qaws_mesh_3d** out_surface, qaws_volume_mesh const* m)
{
	unsigned int *map = NULL, *fv = NULL, *fs = NULL, *src = NULL, *fcell = NULL, *ffac = NULL, *fcor = NULL, nb = 0, nfv = 0, nvb = 0, f, k;
	double* xb = NULL;
	qaws_scalar* vals = NULL;
	qaws_mesh_3d* s = NULL;
	qaws_status st = QAWS_STATUS_ALLOCATION_FAILURE;
	if (!out_surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_surface = NULL;
	if (!m)
		return QAWS_STATUS_INVALID_ARGUMENT;
	map = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->nv + 1)));
	src = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->nv + 1)));
	fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)m->f0[m->nc] + 1)));
	fs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->f0[m->nc] + 1)));
	fcell = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->f0[m->nc] + 1)));
	xb = (double*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(double), (3 * (size_t)m->nv + 1)));
	ffac = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), ((size_t)m->f0[m->nc] + 1)));
	fcor = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (4 * (size_t)m->f0[m->nc] + 1)));
	if (!map || !src || !fv || !fs || !fcell || !xb || !ffac || !fcor)
		goto done;
	for (k = 0; k < m->nv; k++)
		map[k] = HM_NONE;
	for (f = 0; f < m->f0[m->nc]; f++)
	{
		unsigned int v[4], n;
		if (m->fo[f] != HM_NONE)
			continue;
		n = vm_facet_vertices(m, f, v);
		for (k = 0; k < n; k++)
		{
			if (map[v[k]] == HM_NONE)
			{
				memcpy(xb + 3 * (size_t)nvb, m->x + 3 * (size_t)v[k], sizeof(double) * 3);
				src[nvb] = v[k];
				map[v[k]] = nvb++;
			}
			fcor[nfv] = m->c0[vm_facet_cell(m, f)] + vm_facet(m->type[vm_facet_cell(m, f)], f - m->f0[vm_facet_cell(m, f)], &n)[k];
			fv[nfv++] = map[v[k]];
		}
		ffac[nb] = f;
		fcell[nb] = vm_facet_cell(m, f);
		fs[nb++] = n;
	}
	st = hm_mesh_3d_create(xb, nvb, fv, fs, nb, &s);
	if (st != QAWS_STATUS_OK)
		goto done;
	vals = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(qaws_scalar), ((size_t)s->hm.nv + nb + 1)));
	if (!vals)
	{
		st = QAWS_STATUS_ALLOCATION_FAILURE;
		goto done;
	}
	{
		qaws_attribute_desc ad;
		memset(&ad, 0, sizeof ad);
		ad.name = "vertex";
		ad.element = QAWS_MESH_VERTEX;
		ad.components = 1;
		ad.rule = QAWS_ATTRIBUTE_NEAREST;
		/* (copies split off at creation: their source's) */
		for (k = 0; k < s->hm.nv; k++)
			vals[k] = (qaws_scalar)src[s->hm.vsrc[k] < nvb ? s->hm.vsrc[k] : 0];
		st = qaws_mesh_3d_add_attribute(NULL, s, vals, &ad);
		if (st == QAWS_STATUS_OK)
		{
			ad.name = "cell";
			ad.element = QAWS_MESH_FACE;
			for (k = 0; k < nb; k++)
				vals[k] = (qaws_scalar)fcell[k];
			st = qaws_mesh_3d_add_attribute(NULL, s, vals, &ad);
		}
	}
	/* the volume's attributes (corners in face order: each face's run of nfv) */
	for (f = 0; f < m->na && st == QAWS_STATUS_OK; f++)
	{
		vm_attr const* t = &m->attr[f];
		qaws_attribute_desc ad;
		unsigned int cnt = t->element == QAWS_MESH_VERTEX ? s->hm.nv : t->element == QAWS_MESH_CORNER ? nfv : nb, i, j;
		qaws_scalar* av = (qaws_scalar*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(qaws_scalar), t->n, (size_t)cnt + 1));
		if (!av)
		{
			st = QAWS_STATUS_ALLOCATION_FAILURE;
			break;
		}
		for (i = 0; i < cnt; i++)
		{
			unsigned int e = t->element == QAWS_MESH_VERTEX ? src[s->hm.vsrc[i] < nvb ? s->hm.vsrc[i] : 0] : t->element == QAWS_MESH_CORNER ? fcor[i]
				: t->element == QAWS_MESH_CELL ? fcell[i] : ffac[i];
			for (j = 0; j < t->n; j++)
				av[(size_t)t->n * i + j] = (qaws_scalar)t->v[(size_t)t->n * e + j];
		}
		ad.name = t->name;
		ad.element = t->element == QAWS_MESH_VERTEX ? QAWS_MESH_VERTEX : t->element == QAWS_MESH_CORNER ? QAWS_MESH_CORNER : QAWS_MESH_FACE;
		ad.components = t->n;
		ad.rule = (qaws_attribute_rule)t->rule;
		st = qaws_mesh_3d_add_attribute(NULL, s, av, &ad);
		qaws_internal_free(av);
	}
done:
	qaws_internal_free(map); qaws_internal_free(src); qaws_internal_free(fv); qaws_internal_free(fs); qaws_internal_free(fcell); qaws_internal_free(xb); qaws_internal_free(vals);
	qaws_internal_free(ffac); qaws_internal_free(fcor);
	if (st != QAWS_STATUS_OK)
	{
		qaws_mesh_3d_destroy(s);
		return st;
	}
	*out_surface = s;
	return QAWS_STATUS_OK;
}

/* a cell's signed volume: its facets' fans against its centroid */
static double vm_volume(qaws_volume_mesh const* m, unsigned int c)
{
	double cen[3] = { 0, 0, 0 }, vol = 0;
	unsigned int k, f, n = vm_corners[m->type[c]];
	for (k = 0; k < n; k++)
		for (f = 0; f < 3; f++)
			cen[f] += m->x[3 * (size_t)m->cv[m->c0[c] + k] + f] / n;
	for (f = 0; f < vm_facets[m->type[c]]; f++)
	{
		unsigned int nf;
		unsigned char const* loc = vm_facet(m->type[c], f, &nf);
		for (k = 1; k + 1 < nf; k++)
			vol += vm_det(cen, m->x + 3 * (size_t)m->cv[m->c0[c] + loc[0]], m->x + 3 * (size_t)m->cv[m->c0[c] + loc[k]],
				m->x + 3 * (size_t)m->cv[m->c0[c] + loc[k + 1]]) / 6;
	}
	return vol;
}

qaws_status qaws_volume_mesh_validate(qaws_volume_mesh const* m)
{
	unsigned int c, f;
	if (!m)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (c = 0; c < m->nc; c++)
		if (!(vm_volume(m, c) > 0))
			return QAWS_STATUS_INVALID_TOPOLOGY;
	for (f = 0; f < m->f0[m->nc]; f++)
	{
		unsigned int a[4], b[4], n;
		if (m->fo[f] == HM_NONE)
			continue;
		if (m->fo[f] >= m->f0[m->nc] || m->fo[m->fo[f]] != f)
			return QAWS_STATUS_INVALID_TOPOLOGY;
		n = vm_facet_vertices(m, f, a);
		if (vm_facet_vertices(m, m->fo[f], b) != n || !vm_opposite(a, b, n))
			return QAWS_STATUS_INVALID_TOPOLOGY;
	}
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Quality                                                                 */
/* ======================================================================== */

static void vm_sub(double const* a, double const* b, double* o)
{
	o[0] = a[0] - b[0];
	o[1] = a[1] - b[1];
	o[2] = a[2] - b[2];
}

static void vm_cross(double const* a, double const* b, double* o)
{
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

static double vm_dot(double const* a, double const* b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double vm_hex_quality(double const* x, unsigned int const* v)
{
	static unsigned char const nb[8][3] = { { 1, 3, 4 }, { 2, 0, 5 }, { 3, 1, 6 }, { 0, 2, 7 }, { 7, 5, 0 }, { 4, 6, 1 }, { 5, 7, 2 }, { 6, 4, 3 } };
	double q = 1;
	unsigned int k;
	for (k = 0; k < 8; k++)
	{
		double a[3], b[3], d[3], cr[3], l;
		vm_sub(x + 3 * (size_t)v[nb[k][0]], x + 3 * (size_t)v[k], a);
		vm_sub(x + 3 * (size_t)v[nb[k][1]], x + 3 * (size_t)v[k], b);
		vm_sub(x + 3 * (size_t)v[nb[k][2]], x + 3 * (size_t)v[k], d);
		vm_cross(b, d, cr);
		l = sqrt(vm_dot(a, a) * vm_dot(b, b) * vm_dot(d, d));
		q = fmin(q, l > 0 ? vm_dot(a, cr) / l : -1);
	}
	return q;
}

qaws_status qaws_volume_mesh_measure_quality(qaws_volume_quality* out, qaws_volume_mesh const* m)
{
	double mind = 180, maxd = 0, maxre = 0, minsj = 1, mine = HUGE_VAL, maxe = 0, sume = 0, vol = 0;
	unsigned int c, k, ne = 0;
	if (!m || !out)
		return QAWS_STATUS_INVALID_ARGUMENT;
	memset(out, 0, sizeof *out);
	for (c = 0; c < m->nc; c++)
	{
		unsigned int const* v = m->cv + m->c0[c];
		unsigned int n;
		unsigned char const(*e)[2] = vm_edges(m->type[c], &n);
		double shortest = HUGE_VAL;
		vol += vm_volume(m, c);
		for (k = 0; k < n; k++)
		{
			double d[3], l;
			vm_sub(m->x + 3 * (size_t)v[e[k][1]], m->x + 3 * (size_t)v[e[k][0]], d);
			l = sqrt(vm_dot(d, d));
			mine = fmin(mine, l);
			maxe = fmax(maxe, l);
			shortest = fmin(shortest, l);
			sume += l;
			ne++;
		}
		switch (m->type[c])
		{
		case QAWS_CELL_TETRAHEDRON:
		{
			double nrm[4][3], d1[3], d2[3], d3[3], cr[3], l1, l2, l3, det, cc[3];
			unsigned int i, j;
			out->tetrahedra++;
			for (i = 0; i < 4; i++)
			{
				unsigned int nf;
				unsigned char const* loc = vm_facet(QAWS_CELL_TETRAHEDRON, i, &nf);
				double a[3], b[3];
				vm_sub(m->x + 3 * (size_t)v[loc[1]], m->x + 3 * (size_t)v[loc[0]], a);
				vm_sub(m->x + 3 * (size_t)v[loc[2]], m->x + 3 * (size_t)v[loc[0]], b);
				vm_cross(a, b, nrm[i]);
			}
			/* the dihedral angle at the edge between facets i and j (opposite vertices i and j) */
			for (i = 0; i < 4; i++)
				for (j = i + 1; j < 4; j++)
				{
					double cs = -vm_dot(nrm[i], nrm[j]) / sqrt(vm_dot(nrm[i], nrm[i]) * vm_dot(nrm[j], nrm[j])), ang;
					ang = acos(cs > 1 ? 1 : cs < -1 ? -1 : cs) * 180 / 3.14159265358979323846;
					mind = fmin(mind, ang);
					maxd = fmax(maxd, ang);
				}
			/* the circumcentre from vertex 0 */
			vm_sub(m->x + 3 * (size_t)v[1], m->x + 3 * (size_t)v[0], d1);
			vm_sub(m->x + 3 * (size_t)v[2], m->x + 3 * (size_t)v[0], d2);
			vm_sub(m->x + 3 * (size_t)v[3], m->x + 3 * (size_t)v[0], d3);
			l1 = vm_dot(d1, d1); l2 = vm_dot(d2, d2); l3 = vm_dot(d3, d3);
			vm_cross(d2, d3, cr);
			det = 2 * vm_dot(d1, cr);
			for (i = 0; i < 3; i++)
				cc[i] = 0;
			if (det != 0)
			{
				double c23[3], c31[3], c12[3];
				vm_cross(d2, d3, c23);
				vm_cross(d3, d1, c31);
				vm_cross(d1, d2, c12);
				for (i = 0; i < 3; i++)
					cc[i] = (l1 * c23[i] + l2 * c31[i] + l3 * c12[i]) / det;
				maxre = fmax(maxre, sqrt(vm_dot(cc, cc)) / shortest);
			}
			else
				maxre = HUGE_VAL;
			break;
		}
		case QAWS_CELL_HEXAHEDRON:
			out->hexahedra++;
			minsj = fmin(minsj, vm_hex_quality(m->x, v));
			break;
		case QAWS_CELL_PRISM:
			out->prisms++;
			break;
		default:
			out->pyramids++;
			break;
		}
	}
	out->min_dihedral = (qaws_scalar)(out->tetrahedra ? mind : 0);
	out->max_dihedral = (qaws_scalar)(out->tetrahedra ? maxd : 0);
	out->max_radius_edge = (qaws_scalar)maxre;
	out->min_scaled_jacobian = (qaws_scalar)(out->hexahedra ? minsj : 0);
	out->min_edge = (qaws_scalar)(ne ? mine : 0);
	out->max_edge = (qaws_scalar)maxe;
	out->mean_edge = (qaws_scalar)(ne ? sume / ne : 0);
	out->volume = (qaws_scalar)vol;
	return QAWS_STATUS_OK;
}
