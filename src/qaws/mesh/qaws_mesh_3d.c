/*
 * The public 3D mesh over the corner store: qaws_scalar arrays in and out,
 * edits checked and compacted (dense indices) one by one.
 */

#include "qaws_mesh_core.h"
#include "../internal/qaws_internal_memory.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define M3_PI 3.14159265358979323846

qaws_status hm_mesh_3d_create(double const* x, unsigned int nv, unsigned int const* fv, unsigned int const* fs, unsigned int nf,
	qaws_mesh_3d** out)
{
	qaws_mesh_3d* m = (qaws_mesh_3d*)qaws_internal_zalloc(1, sizeof *m);
	qaws_status st;
	*out = NULL;
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = hm_build(&m->hm, 3, x, nv, fv, fs, nf);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_free(m);
		return st;
	}
	*out = m;
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_create(qaws_mesh_3d** out_mesh, qaws_mesh_3d_desc const* desc)
{
	double* x;
	unsigned int i;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!desc || (desc->vertex_count && !desc->positions) || (desc->face_count && !desc->face_vertices))
		return QAWS_STATUS_INVALID_ARGUMENT;
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (desc->vertex_count ? desc->vertex_count : 1)));
	if (!x)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (i = 0; i < 3 * desc->vertex_count; i++)
		x[i] = (double)desc->positions[i];
	st = hm_mesh_3d_create(x, desc->vertex_count, desc->face_vertices, desc->face_sizes, desc->face_count, out_mesh);
	qaws_internal_free(x);
	return st;
}

qaws_status qaws_mesh_3d_clone(qaws_mesh_3d** out_mesh, qaws_mesh_3d const* mesh)
{
	qaws_mesh_3d* m;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = (qaws_mesh_3d*)qaws_internal_zalloc(1, sizeof *m);
	if (!m)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	st = hm_copy(&m->hm, &mesh->hm);
	if (st != QAWS_STATUS_OK)
	{
		qaws_internal_free(m);
		return st;
	}
	*out_mesh = m;
	return QAWS_STATUS_OK;
}

void qaws_mesh_3d_destroy(qaws_mesh_3d* mesh)
{
	if (!mesh)
		return;
	hm_free(&mesh->hm);
	qaws_internal_free(mesh);
}

/* ======================================================================== */
/*  Reading                                                                 */
/* ======================================================================== */

unsigned int qaws_mesh_3d_get_vertex_count(qaws_mesh_3d const* mesh) { return mesh ? mesh->hm.nv : 0; }
unsigned int qaws_mesh_3d_get_edge_count(qaws_mesh_3d const* mesh) { return mesh ? hm_edge_count(&mesh->hm) : 0; }
unsigned int qaws_mesh_3d_get_face_count(qaws_mesh_3d const* mesh) { return mesh ? mesh->hm.nf : 0; }
unsigned int qaws_mesh_3d_get_corner_count(qaws_mesh_3d const* mesh) { return mesh ? mesh->hm.nc - mesh->hm.cfree : 0; }

qaws_vec3 qaws_mesh_3d_get_vertex(qaws_mesh_3d const* mesh, unsigned int vertex)
{
	qaws_vec3 r = { 0, 0, 0 };
	if (mesh && vertex < mesh->hm.nv)
	{
		double const* x = mesh->hm.x + 3 * (size_t)vertex;
		r.x = (qaws_scalar)x[0];
		r.y = (qaws_scalar)x[1];
		r.z = (qaws_scalar)x[2];
	}
	return r;
}

qaws_status qaws_mesh_3d_get_positions(qaws_scalar* out_positions, qaws_mesh_3d const* mesh)
{
	unsigned int i;
	if (!mesh || (!out_positions && mesh->hm.nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * mesh->hm.nv; i++)
		out_positions[i] = (qaws_scalar)mesh->hm.x[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_set_positions(qaws_mesh_3d* mesh, qaws_scalar const* positions)
{
	unsigned int i;
	if (!mesh || (!positions && mesh->hm.nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < 3 * mesh->hm.nv; i++)
		mesh->hm.x[i] = (double)positions[i];
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_get_faces(unsigned int* out_face_vertices, unsigned int* out_face_sizes, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	unsigned int f, c, k = 0;
	if (!mesh || (!out_face_vertices && mesh->hm.nf))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	if (!out_face_sizes && m->nc - m->cfree != 3 * m->nf)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (f = 0; f < m->nf; f++)
	{
		if (out_face_sizes)
			out_face_sizes[f] = m->fn[f];
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			out_face_vertices[k++] = m->cv[c];
	}
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_get_edges(unsigned int* out_edge_vertices, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	unsigned int c, k = 0;
	if (!mesh || !out_edge_vertices)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			out_edge_vertices[k++] = m->cv[c];
			out_edge_vertices[k++] = hm_to(m, c);
		}
	return QAWS_STATUS_OK;
}

unsigned int qaws_mesh_3d_get_triangle_count(qaws_mesh_3d const* mesh)
{
	unsigned int f, n = 0;
	if (!mesh)
		return 0;
	for (f = 0; f < mesh->hm.nf; f++)
		n += mesh->hm.fn[f] - 2;
	return n;
}

void qaws_mesh_3d_get_triangle(unsigned int out_vertices[3], qaws_mesh_3d const* mesh, unsigned int triangle)
{
	qaws_hm const* m = &mesh->hm;
	unsigned int f = 0, t = triangle, c;
	if (m->nc - m->cfree == 3 * m->nf)
		f = triangle, t = 0;   /* (a triangle mesh: the faces) */
	else
		while (f < m->nf && t >= m->fn[f] - 2)
			t -= m->fn[f++] - 2;
	if (f >= m->nf)
	{
		out_vertices[0] = out_vertices[1] = out_vertices[2] = 0;
		return;
	}
	c = m->fc[f];
	out_vertices[0] = m->cv[c];
	out_vertices[1] = m->cv[c + t + 1];
	out_vertices[2] = m->cv[c + t + 2];
}

qaws_status qaws_mesh_3d_get_triangles(unsigned int* out_triangles, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	unsigned int f, t, k = 0;
	if (!mesh || !out_triangles)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	for (f = 0; f < m->nf; f++)
		for (t = 0; t + 2 < m->fn[f]; t++)
		{
			out_triangles[k++] = m->cv[m->fc[f]];
			out_triangles[k++] = m->cv[m->fc[f] + t + 1];
			out_triangles[k++] = m->cv[m->fc[f] + t + 2];
		}
	return QAWS_STATUS_OK;
}

unsigned int qaws_mesh_3d_get_vertex_source(qaws_mesh_3d const* mesh, unsigned int vertex)
{
	return mesh && vertex < mesh->hm.nv ? mesh->hm.vsrc[vertex] : QAWS_MESH_NONE;
}

qaws_status qaws_mesh_3d_get_vertex_ring(unsigned int* out_vertices, unsigned int* out_count, qaws_mesh_3d const* mesh,
	unsigned int capacity, unsigned int vertex)
{
	unsigned int n;
	if (!mesh || vertex >= mesh->hm.nv || (!out_vertices && capacity) || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	n = hm_ring(&mesh->hm, vertex, out_vertices, capacity);
	*out_count = n;
	return n > capacity ? QAWS_STATUS_BUFFER_TOO_SMALL : QAWS_STATUS_OK;
}

int qaws_mesh_3d_is_boundary_vertex(qaws_mesh_3d const* mesh, unsigned int vertex)
{
	return mesh && vertex < mesh->hm.nv && hm_is_border_vertex(&mesh->hm, vertex);
}

static double m3_norm(double* v)
{
	double l = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	if (l > 0)
	{
		v[0] /= l; v[1] /= l; v[2] /= l;
	}
	return l;
}

/* the angle at corner c of its face */
static double m3_corner_angle(qaws_hm const* m, unsigned int c)
{
	double const *p = m->x + 3 * (size_t)m->cv[c], *a = m->x + 3 * (size_t)m->cv[hm_prev(m, c)], *b = m->x + 3 * (size_t)hm_to(m, c);
	double u[3], w[3], d;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		u[k] = a[k] - p[k];
		w[k] = b[k] - p[k];
	}
	if (m3_norm(u) == 0 || m3_norm(w) == 0)
		return 0;
	d = u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
	return acos(d > 1 ? 1 : d < -1 ? -1 : d);
}

qaws_status qaws_mesh_3d_get_vertex_normals(qaws_scalar* out_normals, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	double* acc;
	unsigned int f, c, v, k;
	if (!mesh || (!out_normals && mesh->hm.nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	acc = (double*)qaws_internal_scratch_zalloc(3 * (size_t)m->nv + 1, sizeof(double));
	if (!acc)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	for (f = 0; f < m->nf; f++)
	{
		double n[3];
		hm_face_normal(m, f, n);
		if (m3_norm(n) == 0)
			continue;
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
		{
			double w = m3_corner_angle(m, c);
			for (k = 0; k < 3; k++)
				acc[3 * (size_t)m->cv[c] + k] += w * n[k];
		}
	}
	for (v = 0; v < m->nv; v++)
	{
		m3_norm(acc + 3 * (size_t)v);
		for (k = 0; k < 3; k++)
			out_normals[3 * (size_t)v + k] = (qaws_scalar)acc[3 * (size_t)v + k];
	}
	qaws_internal_free(acc);
	return QAWS_STATUS_OK;
}

static unsigned int m3_find(unsigned int* uf, unsigned int i)
{
	while (uf[i] != i)
	{
		uf[i] = uf[uf[i]];
		i = uf[i];
	}
	return i;
}

qaws_status qaws_mesh_3d_get_topology(qaws_mesh_topology* out_topology, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	unsigned int *uf, v, c, f, roots = 0, loops = 0;
	unsigned char* seen;
	long long chi;
	int* comp_chi;
	unsigned int* comp_b;
	if (!mesh || !out_topology)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	memset(out_topology, 0, sizeof *out_topology);
	uf = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (2 * (size_t)m->nv + 1)));
	seen = (unsigned char*)qaws_internal_scratch_zalloc((size_t)m->nc + 1, 1);
	comp_chi = (int*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, sizeof(int));
	comp_b = (unsigned int*)qaws_internal_scratch_zalloc((size_t)m->nv + 1, sizeof(unsigned int));
	if (!uf || !seen || !comp_chi || !comp_b)
	{
		qaws_internal_free(uf); qaws_internal_free(seen); qaws_internal_free(comp_chi); qaws_internal_free(comp_b);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (v = 0; v < m->nv; v++)
		uf[v] = v;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE)
		{
			unsigned int a = m3_find(uf, m->cv[c]), b = m3_find(uf, hm_to(m, c));
			if (a != b)
				uf[a] = b;
		}
	/* per component: chi = V - E + F, and its boundary loops */
	for (v = 0; v < m->nv; v++)
		if (m->vc[v] != HM_NONE)
			comp_chi[m3_find(uf, v)]++;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
			comp_chi[m3_find(uf, m->cv[c])]--;
	for (f = 0; f < m->nf; f++)
		comp_chi[m3_find(uf, m->cv[m->fc[f]])]++;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && m->co[c] == HM_NONE && !seen[c])
		{
			/* a border loop: from a border corner to the next one round its far vertex, clockwise */
			unsigned int d = c, guard = 0;
			do
			{
				seen[d] = 1;
				d = hm_next(m, d);
				while (m->co[d] != HM_NONE && ++guard < 4 * m->nc)
					d = hm_next(m, m->co[d]);
			} while (d != c && !seen[d] && ++guard < 4 * m->nc);
			comp_b[m3_find(uf, m->cv[c])]++;
			loops++;
		}
	chi = 0;
	for (v = 0; v < m->nv; v++)
		if (m->vc[v] != HM_NONE && m3_find(uf, v) == v)
		{
			int g2 = 2 - comp_chi[v] - (int)comp_b[v];
			roots++;
			chi += comp_chi[v];
			if (g2 > 0)
				out_topology->genus += (unsigned int)(g2 / 2);
		}
	out_topology->components = roots;
	out_topology->boundary_loops = loops;
	out_topology->euler_characteristic = (int)chi;
	out_topology->closed = loops == 0;
	qaws_internal_free(uf); qaws_internal_free(seen); qaws_internal_free(comp_chi); qaws_internal_free(comp_b);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_validate(qaws_mesh_3d const* mesh)
{
	return mesh ? hm_check(&mesh->hm) : QAWS_STATUS_INVALID_ARGUMENT;
}

qaws_status qaws_mesh_3d_get_face_adjacency(unsigned int* out_faces, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	unsigned int f, c, k = 0;
	if (!mesh || (!out_faces && mesh->hm.nf))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	for (f = 0; f < m->nf; f++)
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			out_faces[k++] = m->co[c] == HM_NONE ? QAWS_MESH_NONE : m->cf[m->co[c]];
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_create_dual(qaws_mesh_3d** out_dual, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	double* x;
	unsigned int *fv, *fs, f, c, v, k, nf = 0, nc = 0;
	qaws_status st;
	if (!out_dual)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_dual = NULL;
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	x = (double*)qaws_internal_scratch_zalloc(3 * (size_t)m->nf + 1, sizeof(double));
	fv = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nc + 1)));
	fs = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul2(sizeof(unsigned int), (m->nv + 1)));
	if (!x || !fv || !fs)
	{
		qaws_internal_free(x); qaws_internal_free(fv); qaws_internal_free(fs);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}
	for (f = 0; f < m->nf; f++)
	{
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
			for (k = 0; k < 3; k++)
				x[3 * (size_t)f + k] += m->x[3 * (size_t)m->cv[c] + k] / m->fn[f];
	}
	for (v = 0; v < m->nv; v++)
	{
		unsigned int n = 0;
		if (m->vc[v] == HM_NONE || hm_is_border_vertex(m, v))
			continue;
		c = m->vc[v];
		do
		{
			fv[nc + n++] = m->cf[c];
			c = hm_ccw(m, c);
		} while (c != m->vc[v]);
		if (n < 3)
			continue;
		fs[nf++] = n;
		nc += n;
	}
	st = hm_mesh_3d_create(x, m->nf, fv, fs, nf, out_dual);
	qaws_internal_free(x); qaws_internal_free(fv); qaws_internal_free(fs);
	return st;
}

/* ======================================================================== */
/*  Editing                                                                 */
/* ======================================================================== */

static void m3_set_foot(qaws_hm* m, unsigned int v, qaws_mesh_location const* foot)
{
	if (!m->foot || v == HM_NONE)
		return;
	if (foot)
		m->foot[v] = *foot;
	else
	{
		m->foot[v].patch = HM_NONE;
		m->foot[v].a = m->foot[v].b = 0;
	}
}

/* the parameter of p projected on the segment from vertex a to vertex b, in [0, 1] */
static double m3_param(qaws_hm const* m, unsigned int a, unsigned int b, double const* p)
{
	double const *xa = m->x + 3 * (size_t)a, *xb = m->x + 3 * (size_t)b;
	double num = 0, den = 0, t;
	unsigned int k;
	for (k = 0; k < 3; k++)
	{
		num += (p[k] - xa[k]) * (xb[k] - xa[k]);
		den += (xb[k] - xa[k]) * (xb[k] - xa[k]);
	}
	t = den > 0 ? num / den : 0.5;
	return t < 0 ? 0 : t > 1 ? 1 : t;
}

qaws_status qaws_mesh_3d_split_edge(unsigned int* out_vertex, qaws_mesh_3d* mesh, qaws_scalar const* position, qaws_mesh_location const* footpoint,
	unsigned int a, unsigned int b)
{
	qaws_hm* m;
	unsigned int c, v, k;
	double x[3], t;
	if (out_vertex)
		*out_vertex = QAWS_MESH_NONE;
	if (!mesh || a >= mesh->hm.nv || b >= mesh->hm.nv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	c = hm_find_edge(m, a, b);
	if (c == HM_NONE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (k = 0; k < 3; k++)
		x[k] = position ? (double)position[k] : 0.5 * (m->x[3 * (size_t)a + k] + m->x[3 * (size_t)b + k]);
	t = m3_param(m, m->cv[c], hm_to(m, c), x);
	v = hm_split_edge(m, c, t, x, 1);
	if (v == HM_NONE)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m3_set_foot(m, v, footpoint);
	hm_compact(m);
	if (out_vertex)
		*out_vertex = v;
	return QAWS_STATUS_OK;
}

/* why a collapse of a into b (c a corner of their edge) is refused, or 0 */
static int m3_collapse_refused(qaws_hm const* m, unsigned int c, unsigned int a)
{
	if (m->vflag[a] & HM_CORNER)
		return 1;
	if (hm_feature_degree(m, a) && !(m->cflag[c] & HM_FEATURE))
		return 1;
	return !hm_can_collapse(m, c, a);
}

int qaws_mesh_3d_can_collapse_edge(qaws_mesh_3d const* mesh, unsigned int a, unsigned int b)
{
	unsigned int c;
	if (!mesh || a >= mesh->hm.nv || b >= mesh->hm.nv)
		return 0;
	c = hm_find_edge(&mesh->hm, a, b);
	return c != HM_NONE && !m3_collapse_refused(&mesh->hm, c, a);
}

qaws_status qaws_mesh_3d_collapse_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b, qaws_scalar const* position,
	qaws_mesh_location const* footpoint)
{
	qaws_hm* m;
	unsigned int c, k;
	double x[3], t = 1;
	if (!mesh || a >= mesh->hm.nv || b >= mesh->hm.nv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	c = hm_find_edge(m, a, b);
	if (c == HM_NONE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (m3_collapse_refused(m, c, a))
		return QAWS_STATUS_INVALID_TOPOLOGY;
	for (k = 0; k < 3; k++)
		x[k] = position ? (double)position[k] : m->x[3 * (size_t)b + k];
	if (position)
		t = m3_param(m, a, b, x);
	if (!hm_collapse(m, c, a, t, x))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m3_set_foot(m, b, footpoint);
	hm_compact(m);
	return QAWS_STATUS_OK;
}

int qaws_mesh_3d_can_flip_edge(qaws_mesh_3d const* mesh, unsigned int a, unsigned int b)
{
	unsigned int c;
	if (!mesh || a >= mesh->hm.nv || b >= mesh->hm.nv)
		return 0;
	c = hm_find_edge(&mesh->hm, a, b);
	return c != HM_NONE && !(mesh->hm.cflag[c] & HM_FEATURE) && hm_can_flip(&mesh->hm, c);
}

qaws_status qaws_mesh_3d_flip_edge(qaws_mesh_3d* mesh, unsigned int a, unsigned int b)
{
	qaws_hm* m;
	unsigned int c;
	if (!mesh || a >= mesh->hm.nv || b >= mesh->hm.nv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	c = hm_find_edge(m, a, b);
	if (c == HM_NONE)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if ((m->cflag[c] & HM_FEATURE) || !hm_can_flip(m, c))
		return QAWS_STATUS_INVALID_TOPOLOGY;
	if (!hm_flip(m, c))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	hm_compact(m);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_split_face(unsigned int* out_vertex, qaws_mesh_3d* mesh, qaws_scalar const* position,
	qaws_mesh_location const* footpoint, unsigned int face)
{
	qaws_hm* m;
	unsigned int k, v;
	double x[3] = { 0, 0, 0 };
	if (out_vertex)
		*out_vertex = QAWS_MESH_NONE;
	if (!mesh || face >= mesh->hm.nf)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	for (k = 0; position && k < 3; k++)
		x[k] = (double)position[k];
	v = hm_split_face(m, face, position ? x : NULL);
	if (v == HM_NONE)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m3_set_foot(m, v, footpoint);
	hm_compact(m);
	if (out_vertex)
		*out_vertex = v;
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Features                                                                */
/* ======================================================================== */

static void m3_set_feature(qaws_hm* m, unsigned int c, int on)
{
	unsigned int o = m->co[c];
	if (on)
		m->cflag[c] |= HM_FEATURE;
	else
		m->cflag[c] &= (unsigned char)(0xFFu & ~HM_FEATURE);
	if (o != HM_NONE)
		m->cflag[o] = (unsigned char)((m->cflag[o] & ~HM_FEATURE) | (m->cflag[c] & HM_FEATURE));
}

/* corners: 1 or more than 2 feature edges, or (corner_angle > 0) a sharp turn */
static void m3_corners(qaws_hm* m, double corner_angle)
{
	unsigned int v, c;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int n = hm_feature_degree(m, v), w[2], k = 0, last = HM_NONE;
		m->vflag[v] &= (unsigned char)(0xFFu & ~HM_CORNER);
		if (n == 0)
			continue;
		if (n != 2)
		{
			m->vflag[v] |= HM_CORNER;
			continue;
		}
		if (corner_angle <= 0)
			continue;
		for (c = m->vc[v]; c != HM_NONE; c = hm_round(m, v, c))
		{
			if ((m->cflag[c] & HM_FEATURE) && k < 2)
				w[k++] = hm_to(m, c);
			last = c;
		}
		if (k < 2 && last != HM_NONE && hm_is_border_vertex(m, v) && (m->cflag[hm_prev(m, last)] & HM_FEATURE))
			w[k++] = m->cv[hm_prev(m, last)];
		if (k == 2)
		{
			double const *p = m->x + 3 * (size_t)v, *a = m->x + 3 * (size_t)w[0], *b = m->x + 3 * (size_t)w[1];
			double u[3], s[3], d;
			unsigned int q;
			for (q = 0; q < 3; q++)
			{
				u[q] = p[q] - a[q];
				s[q] = b[q] - p[q];
			}
			if (m3_norm(u) > 0 && m3_norm(s) > 0)
			{
				d = u[0] * s[0] + u[1] * s[1] + u[2] * s[2];
				if (acos(d > 1 ? 1 : d < -1 ? -1 : d) > corner_angle)
					m->vflag[v] |= HM_CORNER;
			}
		}
	}
}

qaws_status qaws_mesh_3d_set_feature_edges(qaws_mesh_3d* mesh, unsigned int const* edge_vertices, unsigned int edge_count)
{
	qaws_hm* m;
	unsigned int i, c;
	if (!mesh || (edge_count && !edge_vertices))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	for (i = 0; i < edge_count; i++)
		if (edge_vertices[2 * i] >= m->nv || edge_vertices[2 * i + 1] >= m->nv
			|| hm_find_edge(m, edge_vertices[2 * i], edge_vertices[2 * i + 1]) == HM_NONE)
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (c = 0; c < m->nc; c++)
		m->cflag[c] &= (unsigned char)(0xFFu & ~HM_FEATURE);
	for (i = 0; i < edge_count; i++)
		m3_set_feature(m, hm_find_edge(m, edge_vertices[2 * i], edge_vertices[2 * i + 1]), 1);
	m3_corners(m, 0);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_detect_features(qaws_mesh_3d* mesh, qaws_scalar dihedral_angle, qaws_scalar corner_angle,
	unsigned int flags)
{
	qaws_hm* m;
	unsigned int c;
	double cosd = cos((double)dihedral_angle);
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	for (c = 0; c < m->nc; c++)
	{
		unsigned int o = m->co[c];
		int on = 0;
		if (m->cf[c] == HM_NONE || (o != HM_NONE && o < c))
			continue;
		if (o == HM_NONE)
			on = (flags & QAWS_FEATURES_BOUNDARY) != 0;
		else if (dihedral_angle > 0)
		{
			double n0[3], n1[3];
			hm_face_normal(m, m->cf[c], n0);
			hm_face_normal(m, m->cf[o], n1);
			if (m3_norm(n0) > 0 && m3_norm(n1) > 0)
				on = n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] < cosd;
		}
		m3_set_feature(m, c, on);
	}
	m3_corners(m, (double)corner_angle);
	return QAWS_STATUS_OK;
}

/* the corners of face f whose vertex's value is the level, up to 8 */
static unsigned int m3_iso_corners(qaws_hm const* m, unsigned int f, double const* val, unsigned int stride, double level, unsigned int* out)
{
	unsigned int c, e = m->fc[f] + m->fn[f], n = 0;
	for (c = m->fc[f]; c < e && n < 8; c++)
		if (val[(size_t)stride * m->cv[c]] == level)
			out[n++] = c;
	return n;
}

qaws_status qaws_mesh_3d_cut_isolines(unsigned int* out_edge_count, qaws_mesh_3d* mesh, unsigned int attribute, qaws_scalar const* levels,
	unsigned int level_count)
{
	qaws_hm* m;
	unsigned int l, i, n, ne;
	unsigned int* pairs = NULL;
	if (out_edge_count)
		*out_edge_count = 0;
	if (!mesh || (level_count && !levels))
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	if (attribute >= m->na || m->attr[attribute].element != QAWS_MESH_VERTEX)
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (l = 0; l < level_count; l++)
		if (!isfinite((double)levels[l]))
			return QAWS_STATUS_INVALID_ARGUMENT;
	for (l = 0; l < level_count; l++)
	{
		double L = (double)levels[l];
		unsigned int stride = m->attr[attribute].n, c, f, nf0;
		/* the edges the level crosses, as their vertex pairs (each split once; splits leave the others' ends alone) */
		pairs = (unsigned int*)qaws_internal_malloc(qaws_internal_mul3(sizeof(unsigned int), 2, (size_t)m->nc + 1));
		if (!pairs)
			return QAWS_STATUS_ALLOCATION_FAILURE;
		for (c = 0, ne = 0; c < m->nc; c++)
		{
			double const* val = m->attr[attribute].v;
			double va, vb;
			if (m->cf[c] == HM_NONE || (m->co[c] != HM_NONE && m->co[c] < c))
				continue;
			va = val[(size_t)stride * m->cv[c]];
			vb = val[(size_t)stride * hm_to(m, c)];
			if ((va - L) * (vb - L) < 0)
			{
				pairs[2 * (size_t)ne] = m->cv[c];
				pairs[2 * (size_t)ne++ + 1] = hm_to(m, c);
			}
		}
		for (i = 0; i < ne; i++)
		{
			unsigned int a = pairs[2 * (size_t)i], b = pairs[2 * (size_t)i + 1], k, v;
			double const* val = m->attr[attribute].v;
			double t = (L - val[(size_t)stride * a]) / (val[(size_t)stride * b] - val[(size_t)stride * a]), x[3];
			c = hm_find(m, a, b);
			if (c == HM_NONE)
			{
				c = hm_find(m, b, a);
				t = 1 - t;
			}
			if (c == HM_NONE)
				continue;
			for (k = 0; k < 3; k++)
				x[k] = m->x[3 * (size_t)m->cv[c] + k] + t * (m->x[3 * (size_t)hm_to(m, c) + k] - m->x[3 * (size_t)m->cv[c] + k]);
			v = hm_split_edge(m, c, t, x, 0);
			if (v == HM_NONE)
			{
				qaws_internal_free(pairs);
				return QAWS_STATUS_ALLOCATION_FAILURE;
			}
			m->attr[attribute].v[(size_t)stride * v] = L;   /* (exactly on the level, whatever the rounding) */
		}
		qaws_internal_free(pairs);
		pairs = NULL;
		/* each face with two corners on the level: cut between them (or its edge between them marked); quads left in two */
		nf0 = m->nf;
		for (f = 0; f < nf0; f++)
		{
			unsigned int on[8], k;
			if (m->fflag[f] & HM_DELETED)
				continue;
			k = m3_iso_corners(m, f, m->attr[attribute].v, stride, L, on);
			if (k != 2)
				continue;
			{
				unsigned int a = m->cv[on[0]], b = m->cv[on[1]], e;
				if (hm_next(m, on[0]) != on[1] && hm_next(m, on[1]) != on[0] && !hm_cut_face(m, on[0], on[1]))
					return QAWS_STATUS_ALLOCATION_FAILURE;
				e = hm_find_edge(m, a, b);
				if (e != HM_NONE && !(m->cflag[e] & HM_FEATURE))
				{
					m3_set_feature(m, e, 1);
					if (out_edge_count)
						(*out_edge_count)++;
				}
			}
		}
		for (f = 0, n = m->nf; f < n; f++)
			if (!(m->fflag[f] & HM_DELETED) && m->fn[f] == 4)
			{
				/* (by its shorter diagonal) */
				unsigned int c0 = m->fc[f];
				double d0 = 0, d1 = 0;
				unsigned int k;
				for (k = 0; k < 3; k++)
				{
					double u = m->x[3 * (size_t)m->cv[c0] + k] - m->x[3 * (size_t)m->cv[c0 + 2] + k];
					double w = m->x[3 * (size_t)m->cv[c0 + 1] + k] - m->x[3 * (size_t)m->cv[c0 + 3] + k];
					d0 += u * u;
					d1 += w * w;
				}
				if (!hm_cut_face(m, d0 <= d1 ? c0 : c0 + 1, d0 <= d1 ? c0 + 2 : c0 + 3))
					return QAWS_STATUS_ALLOCATION_FAILURE;
			}
	}
	/* (corners added: 1 or more than 2 feature edges; those found before by angle stay) */
	for (i = 0; i < m->nv; i++)
	{
		unsigned int d = (m->vflag[i] & HM_DELETED) ? 0 : hm_feature_degree(m, i);
		if (d && d != 2)
			m->vflag[i] |= HM_CORNER;
	}
	hm_compact(m);
	return QAWS_STATUS_OK;
}


unsigned int qaws_mesh_3d_get_feature_edge_count(qaws_mesh_3d const* mesh)
{
	unsigned int c, n = 0;
	if (!mesh)
		return 0;
	for (c = 0; c < mesh->hm.nc; c++)
		if (mesh->hm.cf[c] != HM_NONE && (mesh->hm.cflag[c] & HM_FEATURE) && (mesh->hm.co[c] == HM_NONE || c < mesh->hm.co[c]))
			n++;
	return n;
}

qaws_status qaws_mesh_3d_get_feature_edges(unsigned int* out_edge_vertices, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	unsigned int c, k = 0;
	if (!mesh || !out_edge_vertices)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->cflag[c] & HM_FEATURE) && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			out_edge_vertices[k++] = m->cv[c];
			out_edge_vertices[k++] = hm_to(m, c);
		}
	return QAWS_STATUS_OK;
}

int qaws_mesh_3d_is_feature_corner(qaws_mesh_3d const* mesh, unsigned int vertex)
{
	return mesh && vertex < mesh->hm.nv && (mesh->hm.vflag[vertex] & HM_CORNER) != 0;
}

/* ======================================================================== */
/*  Attributes                                                              */
/* ======================================================================== */

/* the values of attribute a copied out (out != NULL) or in, in public order */
static void m3_attr_io(qaws_hm* m, unsigned int a, qaws_scalar* out, qaws_scalar const* in)
{
	hm_attr* at = &m->attr[a];
	unsigned int i, j, n = at->n, k = 0, count, f, c;
	if (at->element == QAWS_MESH_CORNER)
	{
		for (f = 0; f < m->nf; f++)
			for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++, k++)
				for (j = 0; j < n; j++)
				{
					if (out)
						out[(size_t)n * k + j] = (qaws_scalar)hm_attr_at(m, a, c)[j];
					else
						hm_attr_at(m, a, c)[j] = in ? (double)in[(size_t)n * k + j] : 0;
				}
		return;
	}
	count = at->element == QAWS_MESH_VERTEX ? m->nv : m->nf;
	for (i = 0; i < count; i++)
		for (j = 0; j < n; j++)
		{
			if (out)
				out[(size_t)n * i + j] = (qaws_scalar)at->v[(size_t)n * i + j];
			else
				at->v[(size_t)n * i + j] = in ? (double)in[(size_t)n * i + j] : 0;
		}
}

qaws_status qaws_mesh_3d_add_attribute(unsigned int* out_attribute, qaws_mesh_3d* mesh, qaws_scalar const* values,
	qaws_attribute_desc const* desc)
{
	unsigned int a;
	if (out_attribute)
		*out_attribute = QAWS_MESH_NONE;
	if (!mesh || !desc || desc->components < 1 || desc->components > 16 || (unsigned int)desc->element > QAWS_MESH_FACE
		|| (unsigned int)desc->rule > QAWS_ATTRIBUTE_NEAREST)
		return QAWS_STATUS_INVALID_ARGUMENT;
	a = hm_attr_add(&mesh->hm, desc->name, (unsigned int)desc->element, (unsigned int)desc->rule, desc->components);
	if (a == HM_NONE)
		return QAWS_STATUS_ALLOCATION_FAILURE;
	m3_attr_io(&mesh->hm, a, NULL, values);
	if (out_attribute)
		*out_attribute = a;
	return QAWS_STATUS_OK;
}

unsigned int qaws_mesh_3d_get_attribute_count(qaws_mesh_3d const* mesh)
{
	return mesh ? mesh->hm.na : 0;
}

qaws_status qaws_mesh_3d_find_attribute(unsigned int* out_attribute, qaws_mesh_3d const* mesh, char const* name)
{
	unsigned int a;
	if (!mesh || !name || !out_attribute)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_attribute = QAWS_MESH_NONE;
	for (a = 0; a < mesh->hm.na; a++)
		if (!strncmp(mesh->hm.attr[a].name, name, sizeof mesh->hm.attr[a].name - 1))
		{
			*out_attribute = a;
			return QAWS_STATUS_OK;
		}
	return QAWS_STATUS_OUT_OF_RANGE;
}

qaws_status qaws_mesh_3d_get_attribute_desc(qaws_attribute_desc* out_desc, qaws_mesh_3d const* mesh, unsigned int attribute)
{
	hm_attr const* at;
	if (!mesh || attribute >= mesh->hm.na || !out_desc)
		return QAWS_STATUS_INVALID_ARGUMENT;
	at = &mesh->hm.attr[attribute];
	out_desc->name = at->name;
	out_desc->element = (qaws_mesh_element)at->element;
	out_desc->components = at->n;
	out_desc->rule = (qaws_attribute_rule)at->rule;
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_get_attribute(qaws_scalar* out_values, qaws_mesh_3d const* mesh, unsigned int attribute)
{
	if (!mesh || attribute >= mesh->hm.na || !out_values)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m3_attr_io((qaws_hm*)&mesh->hm, attribute, out_values, NULL);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_set_attribute(qaws_mesh_3d* mesh, unsigned int attribute, qaws_scalar const* values)
{
	if (!mesh || attribute >= mesh->hm.na || !values)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m3_attr_io(&mesh->hm, attribute, NULL, values);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_remove_attribute(qaws_mesh_3d* mesh, unsigned int attribute)
{
	if (!mesh || attribute >= mesh->hm.na)
		return QAWS_STATUS_INVALID_ARGUMENT;
	hm_attr_remove(&mesh->hm, attribute);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Footpoints                                                              */
/* ======================================================================== */

int qaws_mesh_3d_has_footpoints(qaws_mesh_3d const* mesh)
{
	return mesh && mesh->hm.foot != NULL;
}

qaws_status qaws_mesh_3d_get_footpoints(qaws_mesh_location* out_locations, qaws_mesh_3d const* mesh)
{
	if (!mesh || !mesh->hm.foot || (!out_locations && mesh->hm.nv))
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (mesh->hm.nv)
		memcpy(out_locations, mesh->hm.foot, sizeof(qaws_mesh_location) * mesh->hm.nv);
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_set_footpoints(qaws_mesh_3d* mesh, qaws_mesh_location const* locations)
{
	if (!mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!locations)
	{
		qaws_internal_free(mesh->hm.foot);
		mesh->hm.foot = NULL;
		return QAWS_STATUS_OK;
	}
	if (!hm_foot_enable(&mesh->hm))
		return QAWS_STATUS_ALLOCATION_FAILURE;
	if (mesh->hm.nv)
		memcpy(mesh->hm.foot, locations, sizeof(qaws_mesh_location) * mesh->hm.nv);
	return QAWS_STATUS_OK;
}

/* ======================================================================== */
/*  Quality                                                                 */
/* ======================================================================== */

qaws_status qaws_mesh_3d_measure_quality(qaws_mesh_quality* out_quality, qaws_mesh_3d const* mesh)
{
	qaws_hm const* m;
	unsigned int f, c, v, ne = 0, nt = 0, regular;
	double amin = 1e300, amax = 0, emin = 1e300, emax = 0, esum = 0, qmax = 0, qsum = 0, area = 0;
	if (!mesh || !out_quality)
		return QAWS_STATUS_INVALID_ARGUMENT;
	m = &mesh->hm;
	memset(out_quality, 0, sizeof *out_quality);
	for (f = 0; f < m->nf; f++)
	{
		double n[3];
		if (m->fn[f] == 3)
			out_quality->triangles++;
		else if (m->fn[f] == 4)
			out_quality->quads++;
		else
			out_quality->polygons++;
		hm_face_normal(m, f, n);
		area += m3_norm(n);
		for (c = m->fc[f]; c < m->fc[f] + m->fn[f]; c++)
		{
			double a = m3_corner_angle(m, c) * 180 / M3_PI;
			if (a < amin)
				amin = a;
			if (a > amax)
				amax = a;
		}
		if (m->fn[f] == 3)
		{
			/* circumradius over twice the inradius: abc s / (8 A^2) */
			unsigned int c0 = m->fc[f];
			double la = hm_length(m, c0), lb = hm_length(m, c0 + 1), lc = hm_length(m, c0 + 2), s = 0.5 * (la + lb + lc);
			double A2 = s * (s - la) * (s - lb) * (s - lc), q = A2 > 0 ? la * lb * lc * s / (8 * A2) : 1e30;
			if (q > qmax)
				qmax = q;
			qsum += q;
			nt++;
		}
	}
	for (c = 0; c < m->nc; c++)
		if (m->cf[c] != HM_NONE && (m->co[c] == HM_NONE || c < m->co[c]))
		{
			double l = hm_length(m, c);
			if (l < emin)
				emin = l;
			if (l > emax)
				emax = l;
			esum += l;
			ne++;
		}
	regular = out_quality->quads > out_quality->triangles ? 4 : 6;
	for (v = 0; v < m->nv; v++)
	{
		unsigned int val;
		if (m->vc[v] == HM_NONE || hm_is_border_vertex(m, v))
			continue;
		val = hm_valence(m, v);
		out_quality->valence[val < 3 ? 0 : val > 10 ? 7 : val - 3]++;
		if (val != regular)
			out_quality->irregular++;
	}
	out_quality->min_angle = (qaws_scalar)(m->nf ? amin : 0);
	out_quality->max_angle = (qaws_scalar)amax;
	out_quality->min_edge = (qaws_scalar)(ne ? emin : 0);
	out_quality->max_edge = (qaws_scalar)emax;
	out_quality->mean_edge = (qaws_scalar)(ne ? esum / ne : 0);
	out_quality->max_aspect = (qaws_scalar)qmax;
	out_quality->mean_aspect = (qaws_scalar)(nt ? qsum / nt : 0);
	out_quality->area = (qaws_scalar)area;
	return QAWS_STATUS_OK;
}
