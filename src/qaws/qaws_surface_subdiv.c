#include "qaws_surface_subdiv.h"
#include "qaws_surface.h"
#include "qaws_eval.h"
#include "qaws_inspect.h"
#include "internal/qaws_internal_surface.h"
#include "internal/qaws_internal_curve.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* -----------------------------------------------------------------------
 * Internal data structures
 * ----------------------------------------------------------------------- */

/* Subdivision mesh for iterative refinement */
typedef struct subdiv_mesh
{
	qaws_vec3* verts;
	unsigned int vert_count;
	unsigned int* faces;       /* flat indices */
	unsigned int* face_sizes;  /* 3 or 4 per face */
	unsigned int face_count;
	unsigned int total_indices; /* sum of face_sizes */
} subdiv_mesh;

/* Half-edge representation for adjacency queries */
typedef struct subdiv_edge
{
	unsigned int v0;       /* first vertex */
	unsigned int v1;       /* second vertex */
	unsigned int face_left;  /* face on one side (or UINT_MAX if boundary) */
	unsigned int face_right; /* face on other side (or UINT_MAX if boundary) */
} subdiv_edge;

typedef struct subdiv_adjacency
{
	subdiv_edge* edges;
	unsigned int edge_count;

	/* Per-vertex: valence and list of adjacent edges */
	unsigned int* vert_valence;
	unsigned int** vert_edges;  /* vert_edges[v][k] = edge index */

	/* Per-vertex: list of adjacent faces */
	unsigned int** vert_faces;
	unsigned int* vert_face_count;
} subdiv_adjacency;

/* Impl struct stored in the surface */
typedef struct qaws_surface_subdiv_impl
{
	qaws_vec3* grid;         /* (grid_n+1)*(grid_n+1) evaluated positions */
	unsigned int grid_n;     /* grid resolution per axis */

	/* Direct mesh data for correct rendering of closed surfaces */
	qaws_vec3* mesh_verts;
	qaws_vec3* mesh_normals; /* per-vertex smooth normals */
	unsigned int* mesh_tris; /* 3 indices per triangle */
	unsigned int mesh_vert_count;
	unsigned int mesh_tri_count;
} qaws_surface_subdiv_impl;

#define SUBDIV_INVALID_INDEX 0xFFFFFFFFu
#define SUBDIV_PI QAWS_LITERAL(3.14159265358979323846)

/* -----------------------------------------------------------------------
 * Mesh helpers
 * ----------------------------------------------------------------------- */

static void subdiv_mesh_free(subdiv_mesh* m)
{
	if (m)
	{
		if (m->verts) free(m->verts);
		if (m->faces) free(m->faces);
		if (m->face_sizes) free(m->face_sizes);
	}
}

static int subdiv_mesh_init_from_desc(
	subdiv_mesh* m,
	qaws_surface_subdiv_desc const* desc)
{
	unsigned int total;
	unsigned int i;

	m->vert_count = desc->vertex_count;
	m->face_count = desc->face_count;

	/* Compute total indices */
	total = 0;
	for (i = 0; i < desc->face_count; i++)
		total += desc->face_sizes[i];
	m->total_indices = total;

	/* Allocate and copy vertices */
	m->verts = (qaws_vec3*)malloc(m->vert_count * sizeof(qaws_vec3));
	if (!m->verts) return 0;
	memcpy(m->verts, desc->vertices, m->vert_count * sizeof(qaws_vec3));

	/* Allocate and copy faces */
	m->faces = (unsigned int*)malloc(total * sizeof(unsigned int));
	if (!m->faces) return 0;
	memcpy(m->faces, desc->face_indices, total * sizeof(unsigned int));

	/* Allocate and copy face sizes */
	m->face_sizes = (unsigned int*)malloc(m->face_count * sizeof(unsigned int));
	if (!m->face_sizes) return 0;
	memcpy(m->face_sizes, desc->face_sizes, m->face_count * sizeof(unsigned int));

	return 1;
}

/* -----------------------------------------------------------------------
 * Adjacency computation
 * ----------------------------------------------------------------------- */

static void subdiv_adjacency_free(subdiv_adjacency* adj, unsigned int vert_count)
{
	unsigned int i;
	if (!adj) return;
	if (adj->edges) free(adj->edges);
	if (adj->vert_valence) free(adj->vert_valence);
	if (adj->vert_edges)
	{
		for (i = 0; i < vert_count; i++)
			if (adj->vert_edges[i]) free(adj->vert_edges[i]);
		free(adj->vert_edges);
	}
	if (adj->vert_faces)
	{
		for (i = 0; i < vert_count; i++)
			if (adj->vert_faces[i]) free(adj->vert_faces[i]);
		free(adj->vert_faces);
	}
	if (adj->vert_face_count) free(adj->vert_face_count);
}

/* Find an existing edge between v0 and v1 in the edge array, or return
   SUBDIV_INVALID_INDEX if not found. */
static unsigned int find_edge(
	subdiv_edge const* edges,
	unsigned int edge_count,
	unsigned int v0,
	unsigned int v1)
{
	unsigned int i;
	for (i = 0; i < edge_count; i++)
	{
		if ((edges[i].v0 == v0 && edges[i].v1 == v1) ||
			(edges[i].v0 == v1 && edges[i].v1 == v0))
		{
			return i;
		}
	}
	return SUBDIV_INVALID_INDEX;
}

static int subdiv_adjacency_build(subdiv_adjacency* adj, subdiv_mesh const* m)
{
	unsigned int max_edges;
	unsigned int fi, ei, vi;
	unsigned int offset;
	unsigned int* temp_valence;

	memset(adj, 0, sizeof(subdiv_adjacency));

	/* Worst case: total_indices edges (before dedup, actual is about half) */
	max_edges = m->total_indices;
	adj->edges = (subdiv_edge*)malloc(max_edges * sizeof(subdiv_edge));
	if (!adj->edges) return 0;
	adj->edge_count = 0;

	/* Build edges from faces */
	offset = 0;
	for (fi = 0; fi < m->face_count; fi++)
	{
		unsigned int n = m->face_sizes[fi];
		unsigned int k;
		for (k = 0; k < n; k++)
		{
			unsigned int v0 = m->faces[offset + k];
			unsigned int v1 = m->faces[offset + ((k + 1) % n)];
			unsigned int existing = find_edge(adj->edges, adj->edge_count, v0, v1);
			if (existing == SUBDIV_INVALID_INDEX)
			{
				/* New edge */
				subdiv_edge* e = &adj->edges[adj->edge_count];
				e->v0 = v0;
				e->v1 = v1;
				e->face_left = fi;
				e->face_right = SUBDIV_INVALID_INDEX;
				adj->edge_count++;
			}
			else
			{
				/* Second face adjacent to this edge */
				adj->edges[existing].face_right = fi;
			}
		}
		offset += n;
	}

	/* Compute vertex valence and adjacency lists */
	adj->vert_valence = (unsigned int*)calloc(m->vert_count, sizeof(unsigned int));
	if (!adj->vert_valence) return 0;

	/* Count valences */
	for (ei = 0; ei < adj->edge_count; ei++)
	{
		adj->vert_valence[adj->edges[ei].v0]++;
		adj->vert_valence[adj->edges[ei].v1]++;
	}

	/* Allocate per-vertex edge lists */
	adj->vert_edges = (unsigned int**)calloc(m->vert_count, sizeof(unsigned int*));
	if (!adj->vert_edges) return 0;
	for (vi = 0; vi < m->vert_count; vi++)
	{
		if (adj->vert_valence[vi] > 0)
		{
			adj->vert_edges[vi] = (unsigned int*)malloc(
				adj->vert_valence[vi] * sizeof(unsigned int));
			if (!adj->vert_edges[vi]) return 0;
		}
	}

	/* Fill edge lists */
	temp_valence = (unsigned int*)calloc(m->vert_count, sizeof(unsigned int));
	if (!temp_valence) return 0;

	for (ei = 0; ei < adj->edge_count; ei++)
	{
		unsigned int va = adj->edges[ei].v0;
		unsigned int vb = adj->edges[ei].v1;
		adj->vert_edges[va][temp_valence[va]++] = ei;
		adj->vert_edges[vb][temp_valence[vb]++] = ei;
	}
	free(temp_valence);

	/* Build per-vertex face lists */
	adj->vert_face_count = (unsigned int*)calloc(m->vert_count, sizeof(unsigned int));
	if (!adj->vert_face_count) return 0;

	/* Count faces per vertex */
	offset = 0;
	for (fi = 0; fi < m->face_count; fi++)
	{
		unsigned int n = m->face_sizes[fi];
		unsigned int k;
		for (k = 0; k < n; k++)
		{
			adj->vert_face_count[m->faces[offset + k]]++;
		}
		offset += n;
	}

	adj->vert_faces = (unsigned int**)calloc(m->vert_count, sizeof(unsigned int*));
	if (!adj->vert_faces) return 0;
	for (vi = 0; vi < m->vert_count; vi++)
	{
		if (adj->vert_face_count[vi] > 0)
		{
			adj->vert_faces[vi] = (unsigned int*)malloc(
				adj->vert_face_count[vi] * sizeof(unsigned int));
			if (!adj->vert_faces[vi]) return 0;
		}
	}

	/* Fill face lists (reuse vert_face_count as write index, then restore) */
	{
		unsigned int* write_idx = (unsigned int*)calloc(m->vert_count, sizeof(unsigned int));
		if (!write_idx) return 0;

		offset = 0;
		for (fi = 0; fi < m->face_count; fi++)
		{
			unsigned int n = m->face_sizes[fi];
			unsigned int k;
			for (k = 0; k < n; k++)
			{
				unsigned int v = m->faces[offset + k];
				adj->vert_faces[v][write_idx[v]++] = fi;
			}
			offset += n;
		}
		free(write_idx);
	}

	return 1;
}

/* Check if a vertex is on the boundary (has at least one boundary edge) */
static int vertex_is_boundary(
	subdiv_adjacency const* adj,
	unsigned int vi)
{
	unsigned int k;
	for (k = 0; k < adj->vert_valence[vi]; k++)
	{
		unsigned int ei = adj->vert_edges[vi][k];
		if (adj->edges[ei].face_right == SUBDIV_INVALID_INDEX)
			return 1;
	}
	return 0;
}

/* Check if an edge is on the boundary */
static int edge_is_boundary(subdiv_adjacency const* adj, unsigned int ei)
{
	return adj->edges[ei].face_right == SUBDIV_INVALID_INDEX;
}

/* -----------------------------------------------------------------------
 * Face point computation (for CC)
 * ----------------------------------------------------------------------- */

static qaws_vec3 compute_face_centroid(
	subdiv_mesh const* m,
	unsigned int face_offset,
	unsigned int face_size)
{
	qaws_vec3 c;
	unsigned int k;
	qaws_scalar inv;

	c.x = QAWS_ZERO; c.y = QAWS_ZERO; c.z = QAWS_ZERO;
	for (k = 0; k < face_size; k++)
	{
		unsigned int vi = m->faces[face_offset + k];
		c.x += m->verts[vi].x;
		c.y += m->verts[vi].y;
		c.z += m->verts[vi].z;
	}
	inv = QAWS_ONE / (qaws_scalar)face_size;
	c.x *= inv; c.y *= inv; c.z *= inv;
	return c;
}

/* -----------------------------------------------------------------------
 * Catmull-Clark subdivision step
 * ----------------------------------------------------------------------- */

static int catmull_clark_subdivide(subdiv_mesh* in_mesh, subdiv_mesh* out_mesh)
{
	subdiv_adjacency adj;
	unsigned int nv, ne, nf;
	unsigned int new_vert_count;
	unsigned int new_face_count;
	unsigned int new_total_indices;
	qaws_vec3* new_verts;
	unsigned int* new_faces;
	unsigned int* new_face_sizes;
	qaws_vec3* face_points;
	qaws_vec3* edge_points;
	unsigned int fi, ei, vi;
	unsigned int offset;
	unsigned int face_base, edge_base;
	unsigned int write_idx;
	int ok;

	memset(&adj, 0, sizeof(adj));
	ok = subdiv_adjacency_build(&adj, in_mesh);
	if (!ok)
	{
		subdiv_adjacency_free(&adj, in_mesh->vert_count);
		return 0;
	}

	nv = in_mesh->vert_count;
	ne = adj.edge_count;
	nf = in_mesh->face_count;

	/* New vertex layout:
	   [0..nf-1]       = face points (one per face)
	   [nf..nf+ne-1]   = edge points (one per edge)
	   [nf+ne..nf+ne+nv-1] = updated vertex points */
	new_vert_count = nf + ne + nv;
	face_base = 0;
	edge_base = nf;

	/* Each face of n sides produces n quads, each quad has 4 vertices.
	   Total new faces = sum of face_sizes = in_mesh->total_indices
	   Each new face has 4 indices */
	new_face_count = in_mesh->total_indices;
	new_total_indices = new_face_count * 4;

	/* Allocate output mesh */
	new_verts = (qaws_vec3*)malloc(new_vert_count * sizeof(qaws_vec3));
	new_faces = (unsigned int*)malloc(new_total_indices * sizeof(unsigned int));
	new_face_sizes = (unsigned int*)malloc(new_face_count * sizeof(unsigned int));
	face_points = (qaws_vec3*)malloc(nf * sizeof(qaws_vec3));
	edge_points = (qaws_vec3*)malloc(ne * sizeof(qaws_vec3));

	if (!new_verts || !new_faces || !new_face_sizes || !face_points || !edge_points)
	{
		if (new_verts) free(new_verts);
		if (new_faces) free(new_faces);
		if (new_face_sizes) free(new_face_sizes);
		if (face_points) free(face_points);
		if (edge_points) free(edge_points);
		subdiv_adjacency_free(&adj, in_mesh->vert_count);
		return 0;
	}

	/* Step 1: Compute face points (centroid of face vertices) */
	offset = 0;
	for (fi = 0; fi < nf; fi++)
	{
		face_points[fi] = compute_face_centroid(in_mesh, offset, in_mesh->face_sizes[fi]);
		new_verts[face_base + fi] = face_points[fi];
		offset += in_mesh->face_sizes[fi];
	}

	/* Step 2: Compute edge points */
	for (ei = 0; ei < ne; ei++)
	{
		subdiv_edge const* e = &adj.edges[ei];
		qaws_vec3 v0 = in_mesh->verts[e->v0];
		qaws_vec3 v1 = in_mesh->verts[e->v1];

		if (edge_is_boundary(&adj, ei))
		{
			/* Boundary edge: midpoint */
			edge_points[ei].x = (v0.x + v1.x) * QAWS_LITERAL(0.5);
			edge_points[ei].y = (v0.y + v1.y) * QAWS_LITERAL(0.5);
			edge_points[ei].z = (v0.z + v1.z) * QAWS_LITERAL(0.5);
		}
		else
		{
			/* Interior edge: average of 2 endpoints + 2 face points, divided by 4 */
			qaws_vec3 fl = face_points[e->face_left];
			qaws_vec3 fr = face_points[e->face_right];
			edge_points[ei].x = (v0.x + v1.x + fl.x + fr.x) * QAWS_LITERAL(0.25);
			edge_points[ei].y = (v0.y + v1.y + fl.y + fr.y) * QAWS_LITERAL(0.25);
			edge_points[ei].z = (v0.z + v1.z + fl.z + fr.z) * QAWS_LITERAL(0.25);
		}
		new_verts[edge_base + ei] = edge_points[ei];
	}

	/* Step 3: Compute updated vertex points */
	for (vi = 0; vi < nv; vi++)
	{
		qaws_vec3 p = in_mesh->verts[vi];
		unsigned int valence = adj.vert_valence[vi];

		if (valence == 0)
		{
			/* Isolated vertex, keep as is */
			new_verts[nf + ne + vi] = p;
			continue;
		}

		if (vertex_is_boundary(&adj, vi))
		{
			/* Boundary vertex: find boundary edges and apply boundary rule
			   new_pos = (1/8)*prev + (3/4)*current + (1/8)*next
			   where prev and next are the other endpoints of the two boundary edges */
			qaws_vec3 sum;
			unsigned int boundary_count;
			unsigned int k;

			sum.x = QAWS_ZERO; sum.y = QAWS_ZERO; sum.z = QAWS_ZERO;
			boundary_count = 0;

			for (k = 0; k < valence; k++)
			{
				unsigned int eidx = adj.vert_edges[vi][k];
				if (edge_is_boundary(&adj, eidx))
				{
					unsigned int other = (adj.edges[eidx].v0 == vi)
						? adj.edges[eidx].v1 : adj.edges[eidx].v0;
					sum.x += in_mesh->verts[other].x;
					sum.y += in_mesh->verts[other].y;
					sum.z += in_mesh->verts[other].z;
					boundary_count++;
				}
			}

			if (boundary_count == 2)
			{
				new_verts[nf + ne + vi].x =
					QAWS_LITERAL(0.125) * sum.x + QAWS_LITERAL(0.75) * p.x;
				new_verts[nf + ne + vi].y =
					QAWS_LITERAL(0.125) * sum.y + QAWS_LITERAL(0.75) * p.y;
				new_verts[nf + ne + vi].z =
					QAWS_LITERAL(0.125) * sum.z + QAWS_LITERAL(0.75) * p.z;
			}
			else
			{
				/* Corner or unusual boundary: keep position */
				new_verts[nf + ne + vi] = p;
			}
		}
		else
		{
			/* Interior vertex:
			   F = average of face points of adjacent faces
			   E = average of midpoints of adjacent edges
			   P = original vertex
			   new_pos = (F + 2E + (n-3)P) / n  where n = valence */
			qaws_vec3 F, E;
			unsigned int face_cnt;
			qaws_scalar inv_n;
			unsigned int k;

			F.x = QAWS_ZERO; F.y = QAWS_ZERO; F.z = QAWS_ZERO;
			face_cnt = adj.vert_face_count[vi];
			for (k = 0; k < face_cnt; k++)
			{
				unsigned int fidx = adj.vert_faces[vi][k];
				F.x += face_points[fidx].x;
				F.y += face_points[fidx].y;
				F.z += face_points[fidx].z;
			}
			if (face_cnt > 0)
			{
				qaws_scalar inv_fc = QAWS_ONE / (qaws_scalar)face_cnt;
				F.x *= inv_fc; F.y *= inv_fc; F.z *= inv_fc;
			}

			E.x = QAWS_ZERO; E.y = QAWS_ZERO; E.z = QAWS_ZERO;
			for (k = 0; k < valence; k++)
			{
				unsigned int eidx = adj.vert_edges[vi][k];
				qaws_vec3 mid;
				mid.x = (in_mesh->verts[adj.edges[eidx].v0].x +
					in_mesh->verts[adj.edges[eidx].v1].x) * QAWS_LITERAL(0.5);
				mid.y = (in_mesh->verts[adj.edges[eidx].v0].y +
					in_mesh->verts[adj.edges[eidx].v1].y) * QAWS_LITERAL(0.5);
				mid.z = (in_mesh->verts[adj.edges[eidx].v0].z +
					in_mesh->verts[adj.edges[eidx].v1].z) * QAWS_LITERAL(0.5);
				E.x += mid.x; E.y += mid.y; E.z += mid.z;
			}
			{
				qaws_scalar inv_val = QAWS_ONE / (qaws_scalar)valence;
				E.x *= inv_val; E.y *= inv_val; E.z *= inv_val;
			}

			inv_n = QAWS_ONE / (qaws_scalar)valence;
			{
				qaws_scalar nm3 = (qaws_scalar)(valence - 3);
				new_verts[nf + ne + vi].x =
					(F.x + QAWS_LITERAL(2.0) * E.x + nm3 * p.x) * inv_n;
				new_verts[nf + ne + vi].y =
					(F.y + QAWS_LITERAL(2.0) * E.y + nm3 * p.y) * inv_n;
				new_verts[nf + ne + vi].z =
					(F.z + QAWS_LITERAL(2.0) * E.z + nm3 * p.z) * inv_n;
			}
		}
	}

	/* Step 4: Build new faces.
	   For each original face of n sides, create n quads:
	   For each vertex k in the face:
	     quad = [vertex_point[k], edge_point[k-1 to k], face_point, edge_point[k to k+1]]
	   More precisely, for edge between vertex k and vertex (k+1)%n:
	     quad[k] = [updated_vertex[face[k]], edge_point(k, k+1), face_point, edge_point(k-1, k)] */
	write_idx = 0;
	offset = 0;
	{
		unsigned int face_idx = 0;
		for (fi = 0; fi < nf; fi++)
		{
			unsigned int n = in_mesh->face_sizes[fi];
			unsigned int k;
			for (k = 0; k < n; k++)
			{
				unsigned int v_curr = in_mesh->faces[offset + k];
				unsigned int v_next = in_mesh->faces[offset + ((k + 1) % n)];
				unsigned int v_prev = in_mesh->faces[offset + ((k + n - 1) % n)];
				unsigned int edge_next_idx, edge_prev_idx;

				/* Find edge from v_curr to v_next */
				edge_next_idx = find_edge(adj.edges, adj.edge_count, v_curr, v_next);
				/* Find edge from v_prev to v_curr */
				edge_prev_idx = find_edge(adj.edges, adj.edge_count, v_prev, v_curr);

				/* Quad: [updated_vertex, edge_next, face_point, edge_prev] */
				new_faces[write_idx + 0] = nf + ne + v_curr;       /* updated vertex */
				new_faces[write_idx + 1] = edge_base + edge_next_idx; /* edge point (curr->next) */
				new_faces[write_idx + 2] = face_base + fi;           /* face point */
				new_faces[write_idx + 3] = edge_base + edge_prev_idx; /* edge point (prev->curr) */
				new_face_sizes[face_idx] = 4;
				write_idx += 4;
				face_idx++;
			}
			offset += n;
		}
	}

	/* Assign output mesh */
	out_mesh->verts = new_verts;
	out_mesh->vert_count = new_vert_count;
	out_mesh->faces = new_faces;
	out_mesh->face_sizes = new_face_sizes;
	out_mesh->face_count = new_face_count;
	out_mesh->total_indices = new_total_indices;

	free(face_points);
	free(edge_points);
	subdiv_adjacency_free(&adj, in_mesh->vert_count);
	return 1;
}

/* -----------------------------------------------------------------------
 * Loop subdivision step
 * ----------------------------------------------------------------------- */

static int loop_subdivide(subdiv_mesh* in_mesh, subdiv_mesh* out_mesh)
{
	subdiv_adjacency adj;
	unsigned int nv, ne, nf;
	unsigned int new_vert_count;
	unsigned int new_face_count;
	unsigned int new_total_indices;
	qaws_vec3* new_verts;
	unsigned int* new_faces;
	unsigned int* new_face_sizes;
	unsigned int fi, ei, vi;
	unsigned int offset;
	unsigned int write_idx;
	int ok;

	memset(&adj, 0, sizeof(adj));
	ok = subdiv_adjacency_build(&adj, in_mesh);
	if (!ok)
	{
		subdiv_adjacency_free(&adj, in_mesh->vert_count);
		return 0;
	}

	nv = in_mesh->vert_count;
	ne = adj.edge_count;
	nf = in_mesh->face_count;

	/* New vertex layout:
	   [0..nv-1]     = updated original vertices
	   [nv..nv+ne-1] = new edge vertices */
	new_vert_count = nv + ne;

	/* Each triangle is split into 4 triangles */
	new_face_count = nf * 4;
	new_total_indices = new_face_count * 3;

	new_verts = (qaws_vec3*)malloc(new_vert_count * sizeof(qaws_vec3));
	new_faces = (unsigned int*)malloc(new_total_indices * sizeof(unsigned int));
	new_face_sizes = (unsigned int*)malloc(new_face_count * sizeof(unsigned int));

	if (!new_verts || !new_faces || !new_face_sizes)
	{
		if (new_verts) free(new_verts);
		if (new_faces) free(new_faces);
		if (new_face_sizes) free(new_face_sizes);
		subdiv_adjacency_free(&adj, in_mesh->vert_count);
		return 0;
	}

	/* Step 1: Compute edge vertices */
	for (ei = 0; ei < ne; ei++)
	{
		subdiv_edge const* e = &adj.edges[ei];
		qaws_vec3 v0 = in_mesh->verts[e->v0];
		qaws_vec3 v1 = in_mesh->verts[e->v1];

		if (edge_is_boundary(&adj, ei))
		{
			/* Boundary edge: midpoint */
			new_verts[nv + ei].x = (v0.x + v1.x) * QAWS_LITERAL(0.5);
			new_verts[nv + ei].y = (v0.y + v1.y) * QAWS_LITERAL(0.5);
			new_verts[nv + ei].z = (v0.z + v1.z) * QAWS_LITERAL(0.5);
		}
		else
		{
			/* Interior edge: find the two opposite vertices.
			   For edge (v0,v1) shared by two triangles, find the vertex in each
			   triangle that is neither v0 nor v1. */
			unsigned int opp0 = SUBDIV_INVALID_INDEX;
			unsigned int opp1 = SUBDIV_INVALID_INDEX;

			/* Search face_left for opposite vertex */
			{
				unsigned int f_off = 0;
				unsigned int fi2;
				for (fi2 = 0; fi2 < e->face_left; fi2++)
					f_off += in_mesh->face_sizes[fi2];
				{
					unsigned int k;
					for (k = 0; k < in_mesh->face_sizes[e->face_left]; k++)
					{
						unsigned int v = in_mesh->faces[f_off + k];
						if (v != e->v0 && v != e->v1)
						{
							opp0 = v;
							break;
						}
					}
				}
			}

			/* Search face_right for opposite vertex */
			{
				unsigned int f_off = 0;
				unsigned int fi2;
				for (fi2 = 0; fi2 < e->face_right; fi2++)
					f_off += in_mesh->face_sizes[fi2];
				{
					unsigned int k;
					for (k = 0; k < in_mesh->face_sizes[e->face_right]; k++)
					{
						unsigned int v = in_mesh->faces[f_off + k];
						if (v != e->v0 && v != e->v1)
						{
							opp1 = v;
							break;
						}
					}
				}
			}

			if (opp0 != SUBDIV_INVALID_INDEX && opp1 != SUBDIV_INVALID_INDEX)
			{
				/* New vertex = (3*v0 + 3*v1 + opp0 + opp1) / 8 */
				qaws_vec3 o0 = in_mesh->verts[opp0];
				qaws_vec3 o1 = in_mesh->verts[opp1];
				new_verts[nv + ei].x = (QAWS_LITERAL(3.0) * v0.x +
					QAWS_LITERAL(3.0) * v1.x + o0.x + o1.x) * QAWS_LITERAL(0.125);
				new_verts[nv + ei].y = (QAWS_LITERAL(3.0) * v0.y +
					QAWS_LITERAL(3.0) * v1.y + o0.y + o1.y) * QAWS_LITERAL(0.125);
				new_verts[nv + ei].z = (QAWS_LITERAL(3.0) * v0.z +
					QAWS_LITERAL(3.0) * v1.z + o0.z + o1.z) * QAWS_LITERAL(0.125);
			}
			else
			{
				/* Fallback: midpoint */
				new_verts[nv + ei].x = (v0.x + v1.x) * QAWS_LITERAL(0.5);
				new_verts[nv + ei].y = (v0.y + v1.y) * QAWS_LITERAL(0.5);
				new_verts[nv + ei].z = (v0.z + v1.z) * QAWS_LITERAL(0.5);
			}
		}
	}

	/* Step 2: Compute updated vertex positions */
	for (vi = 0; vi < nv; vi++)
	{
		qaws_vec3 p = in_mesh->verts[vi];
		unsigned int n = adj.vert_valence[vi];

		if (n == 0)
		{
			new_verts[vi] = p;
			continue;
		}

		if (vertex_is_boundary(&adj, vi))
		{
			/* Boundary vertex: (1/8)*prev + (3/4)*current + (1/8)*next */
			qaws_vec3 sum;
			unsigned int boundary_count;
			unsigned int k;

			sum.x = QAWS_ZERO; sum.y = QAWS_ZERO; sum.z = QAWS_ZERO;
			boundary_count = 0;

			for (k = 0; k < n; k++)
			{
				unsigned int eidx = adj.vert_edges[vi][k];
				if (edge_is_boundary(&adj, eidx))
				{
					unsigned int other = (adj.edges[eidx].v0 == vi)
						? adj.edges[eidx].v1 : adj.edges[eidx].v0;
					sum.x += in_mesh->verts[other].x;
					sum.y += in_mesh->verts[other].y;
					sum.z += in_mesh->verts[other].z;
					boundary_count++;
				}
			}

			if (boundary_count == 2)
			{
				new_verts[vi].x =
					QAWS_LITERAL(0.125) * sum.x + QAWS_LITERAL(0.75) * p.x;
				new_verts[vi].y =
					QAWS_LITERAL(0.125) * sum.y + QAWS_LITERAL(0.75) * p.y;
				new_verts[vi].z =
					QAWS_LITERAL(0.125) * sum.z + QAWS_LITERAL(0.75) * p.z;
			}
			else
			{
				new_verts[vi] = p;
			}
		}
		else
		{
			/* Interior vertex: Loop's formula
			   beta = (1/n) * (5/8 - (3/8 + 1/4 * cos(2*pi/n))^2)
			   new_pos = (1 - n*beta) * v + beta * sum(neighbors) */
			qaws_scalar beta;
			qaws_vec3 nbr_sum;
			unsigned int k;

			{
				qaws_scalar cos_val = QAWS_COS(
					QAWS_LITERAL(2.0) * SUBDIV_PI / (qaws_scalar)n);
				qaws_scalar tmp = QAWS_LITERAL(0.375) +
					QAWS_LITERAL(0.25) * cos_val;
				beta = (QAWS_LITERAL(0.625) - tmp * tmp) / (qaws_scalar)n;
			}

			nbr_sum.x = QAWS_ZERO;
			nbr_sum.y = QAWS_ZERO;
			nbr_sum.z = QAWS_ZERO;
			for (k = 0; k < n; k++)
			{
				unsigned int eidx = adj.vert_edges[vi][k];
				unsigned int other = (adj.edges[eidx].v0 == vi)
					? adj.edges[eidx].v1 : adj.edges[eidx].v0;
				nbr_sum.x += in_mesh->verts[other].x;
				nbr_sum.y += in_mesh->verts[other].y;
				nbr_sum.z += in_mesh->verts[other].z;
			}

			{
				qaws_scalar one_minus_nb = QAWS_ONE - (qaws_scalar)n * beta;
				new_verts[vi].x = one_minus_nb * p.x + beta * nbr_sum.x;
				new_verts[vi].y = one_minus_nb * p.y + beta * nbr_sum.y;
				new_verts[vi].z = one_minus_nb * p.z + beta * nbr_sum.z;
			}
		}
	}

	/* Step 3: Build new faces.
	   For each original triangle (a, b, c), with edge midpoints
	   m_ab, m_bc, m_ca, create four triangles:
	     (a, m_ab, m_ca)
	     (b, m_bc, m_ab)
	     (c, m_ca, m_bc)
	     (m_ab, m_bc, m_ca) */
	write_idx = 0;
	offset = 0;
	for (fi = 0; fi < nf; fi++)
	{
		unsigned int a = in_mesh->faces[offset + 0];
		unsigned int b = in_mesh->faces[offset + 1];
		unsigned int c = in_mesh->faces[offset + 2];
		unsigned int m_ab = nv + find_edge(adj.edges, adj.edge_count, a, b);
		unsigned int m_bc = nv + find_edge(adj.edges, adj.edge_count, b, c);
		unsigned int m_ca = nv + find_edge(adj.edges, adj.edge_count, c, a);

		/* Triangle 1: (a, m_ab, m_ca) */
		new_faces[write_idx + 0] = a;
		new_faces[write_idx + 1] = m_ab;
		new_faces[write_idx + 2] = m_ca;
		new_face_sizes[fi * 4 + 0] = 3;
		write_idx += 3;

		/* Triangle 2: (b, m_bc, m_ab) */
		new_faces[write_idx + 0] = b;
		new_faces[write_idx + 1] = m_bc;
		new_faces[write_idx + 2] = m_ab;
		new_face_sizes[fi * 4 + 1] = 3;
		write_idx += 3;

		/* Triangle 3: (c, m_ca, m_bc) */
		new_faces[write_idx + 0] = c;
		new_faces[write_idx + 1] = m_ca;
		new_faces[write_idx + 2] = m_bc;
		new_face_sizes[fi * 4 + 2] = 3;
		write_idx += 3;

		/* Triangle 4: (m_ab, m_bc, m_ca) */
		new_faces[write_idx + 0] = m_ab;
		new_faces[write_idx + 1] = m_bc;
		new_faces[write_idx + 2] = m_ca;
		new_face_sizes[fi * 4 + 3] = 3;
		write_idx += 3;

		offset += 3;
	}

	out_mesh->verts = new_verts;
	out_mesh->vert_count = new_vert_count;
	out_mesh->faces = new_faces;
	out_mesh->face_sizes = new_face_sizes;
	out_mesh->face_count = new_face_count;
	out_mesh->total_indices = new_total_indices;

	subdiv_adjacency_free(&adj, in_mesh->vert_count);
	return 1;
}

/* -----------------------------------------------------------------------
 * Grid sampling: map subdivided mesh onto a regular grid
 * ----------------------------------------------------------------------- */

/* Compute axis-aligned bounding box of mesh */
static void mesh_bounds(
	subdiv_mesh const* m,
	qaws_vec3* out_min,
	qaws_vec3* out_max)
{
	unsigned int i;
	*out_min = m->verts[0];
	*out_max = m->verts[0];
	for (i = 1; i < m->vert_count; i++)
	{
		if (m->verts[i].x < out_min->x) out_min->x = m->verts[i].x;
		if (m->verts[i].y < out_min->y) out_min->y = m->verts[i].y;
		if (m->verts[i].z < out_min->z) out_min->z = m->verts[i].z;
		if (m->verts[i].x > out_max->x) out_max->x = m->verts[i].x;
		if (m->verts[i].y > out_max->y) out_max->y = m->verts[i].y;
		if (m->verts[i].z > out_max->z) out_max->z = m->verts[i].z;
	}
}

/* Compute parameter coordinates (u,v) for each vertex of the subdivided mesh.
   Uses planar parameterization: project onto the plane of largest extent.
   Returns allocated arrays of u,v coordinates (caller frees). */
static int parameterize_mesh(
	subdiv_mesh const* m,
	qaws_scalar** out_u,
	qaws_scalar** out_v)
{
	qaws_vec3 mn, mx;
	qaws_scalar dx, dy, dz;
	qaws_scalar* u_coords;
	qaws_scalar* v_coords;
	unsigned int i;
	int axis0, axis1; /* indices: 0=x, 1=y, 2=z */
	qaws_scalar range0, range1, min0, min1;

	mesh_bounds(m, &mn, &mx);
	dx = mx.x - mn.x;
	dy = mx.y - mn.y;
	dz = mx.z - mn.z;

	/* Choose 2 axes with largest extents for parameterization */
	if (dx >= dy && dx >= dz)
	{
		axis0 = 0; /* x */
		if (dy >= dz)
			axis1 = 1; /* y */
		else
			axis1 = 2; /* z */
	}
	else if (dy >= dx && dy >= dz)
	{
		axis0 = 1; /* y */
		if (dx >= dz)
			axis1 = 0; /* x */
		else
			axis1 = 2; /* z */
	}
	else
	{
		axis0 = 2; /* z */
		if (dx >= dy)
			axis1 = 0; /* x */
		else
			axis1 = 1; /* y */
	}

	/* Get range for each chosen axis */
	{
		qaws_scalar mins[3], maxs[3];
		mins[0] = mn.x; mins[1] = mn.y; mins[2] = mn.z;
		maxs[0] = mx.x; maxs[1] = mx.y; maxs[2] = mx.z;
		range0 = maxs[axis0] - mins[axis0];
		range1 = maxs[axis1] - mins[axis1];
		min0 = mins[axis0];
		min1 = mins[axis1];
	}

	if (range0 < QAWS_LITERAL(1e-12)) range0 = QAWS_ONE;
	if (range1 < QAWS_LITERAL(1e-12)) range1 = QAWS_ONE;

	u_coords = (qaws_scalar*)malloc(m->vert_count * sizeof(qaws_scalar));
	v_coords = (qaws_scalar*)malloc(m->vert_count * sizeof(qaws_scalar));
	if (!u_coords || !v_coords)
	{
		if (u_coords) free(u_coords);
		if (v_coords) free(v_coords);
		return 0;
	}

	for (i = 0; i < m->vert_count; i++)
	{
		qaws_scalar coords[3];
		coords[0] = m->verts[i].x;
		coords[1] = m->verts[i].y;
		coords[2] = m->verts[i].z;
		u_coords[i] = (coords[axis0] - min0) / range0;
		v_coords[i] = (coords[axis1] - min1) / range1;

		/* Clamp to [0,1] */
		if (u_coords[i] < QAWS_ZERO) u_coords[i] = QAWS_ZERO;
		if (u_coords[i] > QAWS_ONE) u_coords[i] = QAWS_ONE;
		if (v_coords[i] < QAWS_ZERO) v_coords[i] = QAWS_ZERO;
		if (v_coords[i] > QAWS_ONE) v_coords[i] = QAWS_ONE;
	}

	*out_u = u_coords;
	*out_v = v_coords;
	return 1;
}

/* Compute barycentric coordinates of point (px,py) in triangle
   (ax,ay),(bx,by),(cx,cy). Returns 1 if inside, 0 if outside.
   Writes lambda0, lambda1, lambda2. */
static int barycentric_2d(
	qaws_scalar px, qaws_scalar py,
	qaws_scalar ax, qaws_scalar ay,
	qaws_scalar bx, qaws_scalar by,
	qaws_scalar cx, qaws_scalar cy,
	qaws_scalar* l0, qaws_scalar* l1, qaws_scalar* l2)
{
	qaws_scalar v0x = bx - ax, v0y = by - ay;
	qaws_scalar v1x = cx - ax, v1y = cy - ay;
	qaws_scalar v2x = px - ax, v2y = py - ay;
	qaws_scalar d00 = v0x * v0x + v0y * v0y;
	qaws_scalar d01 = v0x * v1x + v0y * v1y;
	qaws_scalar d11 = v1x * v1x + v1y * v1y;
	qaws_scalar d20 = v2x * v0x + v2y * v0y;
	qaws_scalar d21 = v2x * v1x + v2y * v1y;
	qaws_scalar denom = d00 * d11 - d01 * d01;
	qaws_scalar inv_denom;
	qaws_scalar eps = QAWS_LITERAL(1e-8);
	qaws_scalar margin = QAWS_LITERAL(-1e-4);

	if (QAWS_FABS(denom) < eps)
		return 0;

	inv_denom = QAWS_ONE / denom;
	*l1 = (d11 * d20 - d01 * d21) * inv_denom;
	*l2 = (d00 * d21 - d01 * d20) * inv_denom;
	*l0 = QAWS_ONE - *l1 - *l2;

	return (*l0 >= margin && *l1 >= margin && *l2 >= margin);
}

/* Sample the subdivided mesh onto a regular grid using parameterization.
   For each grid point, find the face containing it in parameter space and
   interpolate the 3D position. */
static int sample_mesh_to_grid(
	subdiv_mesh const* m,
	unsigned int grid_n,
	qaws_vec3** out_grid)
{
	unsigned int grid_size;
	qaws_vec3* grid;
	qaws_scalar* u_param;
	qaws_scalar* v_param;
	unsigned int gi, gj;
	unsigned int fi;

	if (!parameterize_mesh(m, &u_param, &v_param))
		return 0;

	grid_size = (grid_n + 1) * (grid_n + 1);
	grid = (qaws_vec3*)malloc(grid_size * sizeof(qaws_vec3));
	if (!grid)
	{
		free(u_param);
		free(v_param);
		return 0;
	}

	/* Initialize grid with zeros */
	memset(grid, 0, grid_size * sizeof(qaws_vec3));

	/* For each grid point, find containing face and interpolate */
	for (gj = 0; gj <= grid_n; gj++)
	{
		for (gi = 0; gi <= grid_n; gi++)
		{
			qaws_scalar gu = (qaws_scalar)gi / (qaws_scalar)grid_n;
			qaws_scalar gv = (qaws_scalar)gj / (qaws_scalar)grid_n;
			unsigned int idx = gj * (grid_n + 1) + gi;
			int found = 0;
			qaws_scalar best_dist = QAWS_LITERAL(1e30);
			unsigned int best_face = 0;
			unsigned int best_offset = 0;
			unsigned int offset;

			/* Search all faces */
			offset = 0;
			for (fi = 0; fi < m->face_count && !found; fi++)
			{
				unsigned int n = m->face_sizes[fi];

				if (n >= 3)
				{
					/* Triangulate face as a fan and test each triangle */
					unsigned int t;
					for (t = 1; t < n - 1 && !found; t++)
					{
						unsigned int ia = m->faces[offset + 0];
						unsigned int ib = m->faces[offset + t];
						unsigned int ic = m->faces[offset + t + 1];
						qaws_scalar l0, l1, l2;

						if (barycentric_2d(
							gu, gv,
							u_param[ia], v_param[ia],
							u_param[ib], v_param[ib],
							u_param[ic], v_param[ic],
							&l0, &l1, &l2))
						{
							grid[idx].x = l0 * m->verts[ia].x +
								l1 * m->verts[ib].x + l2 * m->verts[ic].x;
							grid[idx].y = l0 * m->verts[ia].y +
								l1 * m->verts[ib].y + l2 * m->verts[ic].y;
							grid[idx].z = l0 * m->verts[ia].z +
								l1 * m->verts[ib].z + l2 * m->verts[ic].z;
							found = 1;
						}
					}
				}

				if (!found)
				{
					/* Track closest face centroid as fallback */
					qaws_scalar cu = QAWS_ZERO, cv = QAWS_ZERO;
					unsigned int k;
					qaws_scalar dist;
					for (k = 0; k < n; k++)
					{
						cu += u_param[m->faces[offset + k]];
						cv += v_param[m->faces[offset + k]];
					}
					cu /= (qaws_scalar)n;
					cv /= (qaws_scalar)n;
					dist = (gu - cu) * (gu - cu) + (gv - cv) * (gv - cv);
					if (dist < best_dist)
					{
						best_dist = dist;
						best_face = fi;
						best_offset = offset;
					}
				}

				offset += n;
			}

			if (!found)
			{
				/* Fallback: use nearest face centroid's barycentric projection */
				unsigned int n = m->face_sizes[best_face];
				if (n >= 3)
				{
					unsigned int ia = m->faces[best_offset + 0];
					unsigned int ib = m->faces[best_offset + 1];
					unsigned int ic = m->faces[best_offset + 2];
					qaws_scalar l0, l1, l2;

					barycentric_2d(
						gu, gv,
						u_param[ia], v_param[ia],
						u_param[ib], v_param[ib],
						u_param[ic], v_param[ic],
						&l0, &l1, &l2);

					/* Clamp barycentric coords */
					if (l0 < QAWS_ZERO) l0 = QAWS_ZERO;
					if (l1 < QAWS_ZERO) l1 = QAWS_ZERO;
					if (l2 < QAWS_ZERO) l2 = QAWS_ZERO;
					{
						qaws_scalar sum = l0 + l1 + l2;
						if (sum > QAWS_LITERAL(1e-12))
						{
							qaws_scalar inv = QAWS_ONE / sum;
							l0 *= inv; l1 *= inv; l2 *= inv;
						}
						else
						{
							l0 = QAWS_LITERAL(0.333333);
							l1 = QAWS_LITERAL(0.333333);
							l2 = QAWS_LITERAL(0.333334);
						}
					}

					grid[idx].x = l0 * m->verts[ia].x +
						l1 * m->verts[ib].x + l2 * m->verts[ic].x;
					grid[idx].y = l0 * m->verts[ia].y +
						l1 * m->verts[ib].y + l2 * m->verts[ic].y;
					grid[idx].z = l0 * m->verts[ia].z +
						l1 * m->verts[ib].z + l2 * m->verts[ic].z;
				}
			}
		}
	}

	free(u_param);
	free(v_param);
	*out_grid = grid;
	return 1;
}

/* -----------------------------------------------------------------------
 * Surface evaluation via bilinear interpolation on grid
 * ----------------------------------------------------------------------- */

static void compute_normal(qaws_vec3 du, qaws_vec3 dv, qaws_vec3* out)
{
	qaws_scalar nx = du.y * dv.z - du.z * dv.y;
	qaws_scalar ny = du.z * dv.x - du.x * dv.z;
	qaws_scalar nz = du.x * dv.y - du.y * dv.x;
	qaws_scalar len = QAWS_SQRT(nx * nx + ny * ny + nz * nz);
	if (len > QAWS_LITERAL(1e-12))
	{
		out->x = nx / len; out->y = ny / len; out->z = nz / len;
	}
	else
	{
		out->x = 0; out->y = 0; out->z = 1;
	}
}

/* Bilinear patch:
   S(u,v) = (1-fu)(1-fv)P00 + fu(1-fv)P10 + (1-fu)fv*P01 + fu*fv*P11
   dS/du  = grid_n * [(1-fv)(P10-P00) + fv(P11-P01)]
   dS/dv  = grid_n * [(1-fu)(P01-P00) + fu(P11-P10)]
   d2S/du2  = 0
   d2S/dv2  = 0
   d2S/dudv = grid_n^2 * (P00 - P10 - P01 + P11) */
static qaws_status subdiv_surface_eval(
	qaws_surface const* surface,
	qaws_scalar u,
	qaws_scalar v,
	unsigned int eval_flags,
	qaws_surface_eval_result* out_result)
{
	qaws_surface_subdiv_impl const* impl =
		(qaws_surface_subdiv_impl const*)surface->impl;
	unsigned int gn = impl->grid_n;
	qaws_scalar gn_s = (qaws_scalar)gn;
	qaws_scalar gu, gv, fu, fv;
	unsigned int iu, iv;
	qaws_vec3 p00, p10, p01, p11;
	qaws_scalar one_minus_fu, one_minus_fv;

	/* Map u,v to grid coordinates */
	gu = u * gn_s;
	gv = v * gn_s;
	iu = (unsigned int)QAWS_FLOOR(gu);
	iv = (unsigned int)QAWS_FLOOR(gv);

	/* Clamp to valid grid cell range */
	if (iu >= gn) iu = gn - 1;
	if (iv >= gn) iv = gn - 1;

	fu = gu - (qaws_scalar)iu;
	fv = gv - (qaws_scalar)iv;
	one_minus_fu = QAWS_ONE - fu;
	one_minus_fv = QAWS_ONE - fv;

	/* Get four grid corners */
	p00 = impl->grid[iv * (gn + 1) + iu];
	p10 = impl->grid[iv * (gn + 1) + iu + 1];
	p01 = impl->grid[(iv + 1) * (gn + 1) + iu];
	p11 = impl->grid[(iv + 1) * (gn + 1) + iu + 1];

	/* Position: bilinear interpolation */
	if (eval_flags & QAWS_SURFACE_EVAL_POSITION)
	{
		out_result->position.x = one_minus_fu * one_minus_fv * p00.x
			+ fu * one_minus_fv * p10.x
			+ one_minus_fu * fv * p01.x
			+ fu * fv * p11.x;
		out_result->position.y = one_minus_fu * one_minus_fv * p00.y
			+ fu * one_minus_fv * p10.y
			+ one_minus_fu * fv * p01.y
			+ fu * fv * p11.y;
		out_result->position.z = one_minus_fu * one_minus_fv * p00.z
			+ fu * one_minus_fv * p10.z
			+ one_minus_fu * fv * p01.z
			+ fu * fv * p11.z;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_POSITION;
	}

	/* dS/du = grid_n * [(1-fv)(P10-P00) + fv(P11-P01)] */
	if (eval_flags & (QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_NORMAL))
	{
		out_result->du.x = gn_s * (one_minus_fv * (p10.x - p00.x)
			+ fv * (p11.x - p01.x));
		out_result->du.y = gn_s * (one_minus_fv * (p10.y - p00.y)
			+ fv * (p11.y - p01.y));
		out_result->du.z = gn_s * (one_minus_fv * (p10.z - p00.z)
			+ fv * (p11.z - p01.z));
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DU;
	}

	/* dS/dv = grid_n * [(1-fu)(P01-P00) + fu(P11-P10)] */
	if (eval_flags & (QAWS_SURFACE_EVAL_DV | QAWS_SURFACE_EVAL_NORMAL))
	{
		out_result->dv.x = gn_s * (one_minus_fu * (p01.x - p00.x)
			+ fu * (p11.x - p10.x));
		out_result->dv.y = gn_s * (one_minus_fu * (p01.y - p00.y)
			+ fu * (p11.y - p10.y));
		out_result->dv.z = gn_s * (one_minus_fu * (p01.z - p00.z)
			+ fu * (p11.z - p10.z));
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DV;
	}

	/* d2S/du2 = 0 (bilinear within each cell) */
	if (eval_flags & QAWS_SURFACE_EVAL_DUU)
	{
		out_result->duu.x = QAWS_ZERO;
		out_result->duu.y = QAWS_ZERO;
		out_result->duu.z = QAWS_ZERO;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUU;
	}

	/* d2S/dv2 = 0 */
	if (eval_flags & QAWS_SURFACE_EVAL_DVV)
	{
		out_result->dvv.x = QAWS_ZERO;
		out_result->dvv.y = QAWS_ZERO;
		out_result->dvv.z = QAWS_ZERO;
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DVV;
	}

	/* d2S/dudv = grid_n^2 * (P00 - P10 - P01 + P11) */
	if (eval_flags & QAWS_SURFACE_EVAL_DUV)
	{
		qaws_scalar gn2 = gn_s * gn_s;
		out_result->duv.x = gn2 * (p00.x - p10.x - p01.x + p11.x);
		out_result->duv.y = gn2 * (p00.y - p10.y - p01.y + p11.y);
		out_result->duv.z = gn2 * (p00.z - p10.z - p01.z + p11.z);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_DUV;
	}

	/* Normal */
	if (eval_flags & QAWS_SURFACE_EVAL_NORMAL)
	{
		compute_normal(out_result->du, out_result->dv, &out_result->normal);
		out_result->valid_flags |= QAWS_SURFACE_EVAL_NORMAL;
	}

	return QAWS_STATUS_OK;
}

/* -----------------------------------------------------------------------
 * Destroy / vtable
 * ----------------------------------------------------------------------- */

static void subdiv_surface_destroy(void* impl, qaws_allocator const* allocator)
{
	qaws_surface_subdiv_impl* si = (qaws_surface_subdiv_impl*)impl;
	if (si)
	{
		if (si->grid) free(si->grid);
		if (si->mesh_verts) free(si->mesh_verts);
		if (si->mesh_normals) free(si->mesh_normals);
		if (si->mesh_tris) free(si->mesh_tris);
	}
	qaws_internal_dealloc(allocator, impl);
}

static int subdiv_surface_is_rational(qaws_surface const* s)
{
	(void)s;
	return 0;
}

static qaws_surface_vtable const subdiv_surface_vtable = {
	subdiv_surface_eval,
	subdiv_surface_destroy,
	subdiv_surface_is_rational
};

/* -----------------------------------------------------------------------
 * Build triangle mesh + smooth normals from subdivided mesh
 * ----------------------------------------------------------------------- */

static int build_triangle_mesh(
	subdiv_mesh const* m,
	qaws_vec3** out_verts,
	qaws_vec3** out_normals,
	unsigned int** out_tris,
	unsigned int* out_vert_count,
	unsigned int* out_tri_count)
{
	unsigned int tri_count, fi, offset, ti, vi;
	qaws_vec3* verts;
	qaws_vec3* normals;
	unsigned int* tris;
	qaws_vec3* face_normals;

	/* Count triangles: each n-sided face produces n-2 triangles */
	tri_count = 0;
	for (fi = 0; fi < m->face_count; fi++)
	{
		if (m->face_sizes[fi] >= 3)
			tri_count += m->face_sizes[fi] - 2;
	}

	verts = (qaws_vec3*)malloc(m->vert_count * sizeof(qaws_vec3));
	normals = (qaws_vec3*)calloc(m->vert_count, sizeof(qaws_vec3));
	tris = (unsigned int*)malloc(tri_count * 3 * sizeof(unsigned int));
	face_normals = (qaws_vec3*)malloc(tri_count * sizeof(qaws_vec3));

	if (!verts || !normals || !tris || !face_normals)
	{
		free(verts); free(normals); free(tris); free(face_normals);
		return 0;
	}

	/* Copy vertices */
	memcpy(verts, m->verts, m->vert_count * sizeof(qaws_vec3));

	/* Triangulate faces (fan from first vertex) and compute face normals */
	ti = 0;
	offset = 0;
	for (fi = 0; fi < m->face_count; fi++)
	{
		unsigned int n = m->face_sizes[fi];
		unsigned int k;
		for (k = 1; k < n - 1; k++)
		{
			unsigned int ia = m->faces[offset + 0];
			unsigned int ib = m->faces[offset + k];
			unsigned int ic = m->faces[offset + k + 1];
			qaws_vec3 ab, ac, fn;

			tris[ti * 3 + 0] = ia;
			tris[ti * 3 + 1] = ib;
			tris[ti * 3 + 2] = ic;

			/* Face normal = (B-A) x (C-A) */
			ab.x = verts[ib].x - verts[ia].x;
			ab.y = verts[ib].y - verts[ia].y;
			ab.z = verts[ib].z - verts[ia].z;
			ac.x = verts[ic].x - verts[ia].x;
			ac.y = verts[ic].y - verts[ia].y;
			ac.z = verts[ic].z - verts[ia].z;
			fn.x = ab.y * ac.z - ab.z * ac.y;
			fn.y = ab.z * ac.x - ab.x * ac.z;
			fn.z = ab.x * ac.y - ab.y * ac.x;

			/* Area-weighted: don't normalize yet */
			face_normals[ti] = fn;

			/* Accumulate onto vertices */
			normals[ia].x += fn.x; normals[ia].y += fn.y; normals[ia].z += fn.z;
			normals[ib].x += fn.x; normals[ib].y += fn.y; normals[ib].z += fn.z;
			normals[ic].x += fn.x; normals[ic].y += fn.y; normals[ic].z += fn.z;

			ti++;
		}
		offset += n;
	}

	/* Normalize per-vertex normals */
	for (vi = 0; vi < m->vert_count; vi++)
	{
		qaws_scalar len = QAWS_SQRT(
			normals[vi].x * normals[vi].x +
			normals[vi].y * normals[vi].y +
			normals[vi].z * normals[vi].z);
		if (len > QAWS_LITERAL(1e-12))
		{
			qaws_scalar inv = QAWS_ONE / len;
			normals[vi].x *= inv;
			normals[vi].y *= inv;
			normals[vi].z *= inv;
		}
		else
		{
			normals[vi].x = QAWS_ZERO;
			normals[vi].y = QAWS_ZERO;
			normals[vi].z = QAWS_ONE;
		}
	}

	free(face_normals);

	*out_verts = verts;
	*out_normals = normals;
	*out_tris = tris;
	*out_vert_count = m->vert_count;
	*out_tri_count = tri_count;
	return 1;
}

/* -----------------------------------------------------------------------
 * Public creation function
 * ----------------------------------------------------------------------- */

qaws_status qaws_surface_create_subdiv(
	qaws_surface_subdiv_desc const* desc,
	qaws_surface** out_surface)
{
	qaws_surface* surface;
	qaws_surface_subdiv_impl* impl;
	qaws_range u_range, v_range;
	subdiv_mesh mesh;
	unsigned int level;
	unsigned int grid_n;
	qaws_vec3* grid;
	unsigned int iter;
	unsigned int i;

	if (!desc || !out_surface) return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->vertices || desc->vertex_count < 3)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (!desc->face_indices || !desc->face_sizes || desc->face_count < 1)
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate face sizes */
	for (i = 0; i < desc->face_count; i++)
	{
		if (desc->face_sizes[i] < 3)
			return QAWS_STATUS_INVALID_ARGUMENT;
	}

	/* Validate face indices are in range */
	{
		unsigned int total = 0;
		for (i = 0; i < desc->face_count; i++)
			total += desc->face_sizes[i];
		for (i = 0; i < total; i++)
		{
			if (desc->face_indices[i] >= desc->vertex_count)
				return QAWS_STATUS_INVALID_ARGUMENT;
		}
	}

	/* For Loop scheme, all faces must be triangles */
	if (desc->scheme == QAWS_SUBDIV_LOOP)
	{
		for (i = 0; i < desc->face_count; i++)
		{
			if (desc->face_sizes[i] != 3)
				return QAWS_STATUS_INVALID_ARGUMENT;
		}
	}

	/* Determine subdivision level (default 3) */
	level = desc->subdivision_level;
	if (level == 0) level = 3;

	/* Initialize mesh from descriptor */
	memset(&mesh, 0, sizeof(mesh));
	if (!subdiv_mesh_init_from_desc(&mesh, desc))
	{
		subdiv_mesh_free(&mesh);
		return QAWS_STATUS_ALLOCATION_FAILURE;
	}

	/* Perform subdivision iterations */
	for (iter = 0; iter < level; iter++)
	{
		subdiv_mesh new_mesh;
		int ok;

		memset(&new_mesh, 0, sizeof(new_mesh));

		if (desc->scheme == QAWS_SUBDIV_LOOP)
			ok = loop_subdivide(&mesh, &new_mesh);
		else
			ok = catmull_clark_subdivide(&mesh, &new_mesh);

		subdiv_mesh_free(&mesh);

		if (!ok)
		{
			subdiv_mesh_free(&new_mesh);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		mesh = new_mesh;
	}

	/* Build triangle mesh + normals from subdivided mesh (before freeing) */
	{
		qaws_vec3* mv = NULL;
		qaws_vec3* mn = NULL;
		unsigned int* mt = NULL;
		unsigned int mvc = 0, mtc = 0;

		if (!build_triangle_mesh(&mesh, &mv, &mn, &mt, &mvc, &mtc))
		{
			subdiv_mesh_free(&mesh);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		/* Determine grid resolution */
		grid_n = 1u << (level + 2); /* 2^(level+2) */
		if (grid_n < 16) grid_n = 16;
		if (grid_n > 256) grid_n = 256;

		/* Sample subdivided mesh onto regular grid (best-effort for eval) */
		grid = NULL;
		if (!sample_mesh_to_grid(&mesh, grid_n, &grid))
		{
			subdiv_mesh_free(&mesh);
			free(mv); free(mn); free(mt);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		subdiv_mesh_free(&mesh);

		/* Allocate surface */
		u_range.min_value = 0; u_range.max_value = 1;
		v_range.min_value = 0; v_range.max_value = 1;

		surface = qaws_internal_surface_alloc(
			QAWS_SURFACE_KIND_SUBDIV,
			1, 1, u_range, v_range,
			&subdiv_surface_vtable);
		if (!surface)
		{
			free(grid); free(mv); free(mn); free(mt);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		impl = (qaws_surface_subdiv_impl*)malloc(sizeof(qaws_surface_subdiv_impl));
		if (!impl)
		{
			free(grid); free(mv); free(mn); free(mt);
			qaws_internal_surface_free(surface);
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}

		impl->grid = grid;
		impl->grid_n = grid_n;
		impl->mesh_verts = mv;
		impl->mesh_normals = mn;
		impl->mesh_tris = mt;
		impl->mesh_vert_count = mvc;
		impl->mesh_tri_count = mtc;

		surface->impl = impl;
		*out_surface = surface;
		return QAWS_STATUS_OK;
	}
}

/* -----------------------------------------------------------------------
 * Public mesh accessor
 * ----------------------------------------------------------------------- */

qaws_status qaws_surface_subdiv_get_mesh(
	qaws_surface const* surface,
	qaws_vec3 const** out_vertices,
	qaws_vec3 const** out_normals,
	unsigned int const** out_triangles,
	unsigned int* out_vertex_count,
	unsigned int* out_triangle_count)
{
	qaws_surface_subdiv_impl const* impl;

	if (!surface || !out_vertices || !out_normals || !out_triangles ||
		!out_vertex_count || !out_triangle_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (qaws_surface_get_kind(surface) != QAWS_SURFACE_KIND_SUBDIV)
		return QAWS_STATUS_INVALID_ARGUMENT;

	impl = (qaws_surface_subdiv_impl const*)surface->impl;
	if (!impl || !impl->mesh_verts)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_vertices = impl->mesh_verts;
	*out_normals = impl->mesh_normals;
	*out_triangles = impl->mesh_tris;
	*out_vertex_count = impl->mesh_vert_count;
	*out_triangle_count = impl->mesh_tri_count;
	return QAWS_STATUS_OK;
}
