#include "qaws_internal_mesh.h"
#include "../qaws_platform.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* ------------------------------------------------------------------ */
/*  Hash table for edge deduplication                                  */
/* ------------------------------------------------------------------ */

typedef struct edge_key {
	unsigned int v0, v1; /* canonical: v0 < v1 */
} edge_key;

static unsigned int edge_hash(unsigned int v0, unsigned int v1, unsigned int cap)
{
	unsigned int h = v0 * 2654435761u ^ v1 * 40503u;
	return h % cap;
}

/* ------------------------------------------------------------------ */
/*  Build mesh from tessellation                                       */
/* ------------------------------------------------------------------ */

qaws_status qaws_internal_mesh_build(
	qaws_internal_mesh* mesh,
	qaws_scalar const* positions,
	qaws_scalar const* uvs,
	unsigned int vertex_count,
	unsigned int const* indices,
	unsigned int index_count)
{
	unsigned int tri_count;
	unsigned int max_edges;
	unsigned int hash_cap;
	unsigned int* hash_next = NULL;
	unsigned int* hash_table = NULL;
	edge_key* hash_keys = NULL;
	unsigned int* hash_edge_idx = NULL;
	unsigned int edge_count = 0;
	unsigned int* vert_degree = NULL;
	unsigned int total_vert_edges = 0;
	unsigned int ti, ei, vi;

	if (!mesh || !positions || !indices)
		return QAWS_STATUS_INVALID_ARGUMENT;

	memset(mesh, 0, sizeof(*mesh));
	tri_count = index_count / 3;
	if (tri_count == 0)
		return QAWS_STATUS_INVALID_ARGUMENT;

	mesh->positions = positions;
	mesh->uvs = uvs;
	mesh->vertex_count = vertex_count;
	mesh->tri_count = tri_count;

	mesh->tris = (qaws_internal_tri*)calloc(tri_count, sizeof(qaws_internal_tri));
	if (!mesh->tris) goto fail;

	/* Maximum possible edges: 3 * tri_count (before dedup) */
	max_edges = 3 * tri_count;
	mesh->edges = (qaws_internal_edge*)calloc(max_edges, sizeof(qaws_internal_edge));
	if (!mesh->edges) goto fail;

	/* Hash table for edge dedup */
	hash_cap = max_edges * 2 + 1;
	hash_table = (unsigned int*)malloc(hash_cap * sizeof(unsigned int));
	hash_next = (unsigned int*)malloc(max_edges * sizeof(unsigned int));
	hash_keys = (edge_key*)malloc(max_edges * sizeof(edge_key));
	hash_edge_idx = (unsigned int*)malloc(max_edges * sizeof(unsigned int));
	if (!hash_table || !hash_next || !hash_keys || !hash_edge_idx) goto fail;

	memset(hash_table, 0xFF, hash_cap * sizeof(unsigned int)); /* UINT_MAX = empty */

	/* Fill triangles and build edges */
	for (ti = 0; ti < tri_count; ti++)
	{
		unsigned int v0 = indices[ti * 3 + 0];
		unsigned int v1 = indices[ti * 3 + 1];
		unsigned int v2 = indices[ti * 3 + 2];
		unsigned int local_edge;

		mesh->tris[ti].v[0] = v0;
		mesh->tris[ti].v[1] = v1;
		mesh->tris[ti].v[2] = v2;

		/* Three edges: opposite vertex 0 = edge(v1,v2), opposite 1 = edge(v0,v2), opposite 2 = edge(v0,v1) */
		for (local_edge = 0; local_edge < 3; local_edge++)
		{
			unsigned int a = mesh->tris[ti].v[(local_edge + 1) % 3];
			unsigned int b = mesh->tris[ti].v[(local_edge + 2) % 3];
			unsigned int opp = mesh->tris[ti].v[local_edge];
			unsigned int ea = a < b ? a : b;
			unsigned int eb = a < b ? b : a;
			unsigned int h = edge_hash(ea, eb, hash_cap);
			unsigned int slot;
			int found = 0;

			/* Search hash chain */
			for (slot = hash_table[h]; slot != UINT_MAX; slot = hash_next[slot])
			{
				if (hash_keys[slot].v0 == ea && hash_keys[slot].v1 == eb)
				{
					/* Second adjacency */
					unsigned int eidx = hash_edge_idx[slot];
					mesh->edges[eidx].tri[1] = ti;
					mesh->edges[eidx].opposite[1] = opp;
					mesh->tris[ti].edge[local_edge] = eidx;
					found = 1;
					break;
				}
			}

			if (!found)
			{
				unsigned int eidx = edge_count++;
				mesh->edges[eidx].v[0] = ea;
				mesh->edges[eidx].v[1] = eb;
				mesh->edges[eidx].tri[0] = ti;
				mesh->edges[eidx].tri[1] = UINT_MAX;
				mesh->edges[eidx].opposite[0] = opp;
				mesh->edges[eidx].opposite[1] = UINT_MAX;
				mesh->tris[ti].edge[local_edge] = eidx;

				hash_keys[eidx].v0 = ea;
				hash_keys[eidx].v1 = eb;
				hash_edge_idx[eidx] = eidx;
				hash_next[eidx] = hash_table[h];
				hash_table[h] = eidx;
			}
		}
	}

	mesh->edge_count = edge_count;

	free(hash_table); hash_table = NULL;
	free(hash_next); hash_next = NULL;
	free(hash_keys); hash_keys = NULL;
	free(hash_edge_idx); hash_edge_idx = NULL;

	/* Build vertex-to-edge adjacency (CSR) */
	vert_degree = (unsigned int*)calloc(vertex_count, sizeof(unsigned int));
	if (!vert_degree) goto fail;

	for (ei = 0; ei < edge_count; ei++)
	{
		vert_degree[mesh->edges[ei].v[0]]++;
		vert_degree[mesh->edges[ei].v[1]]++;
	}

	mesh->vert_edge_offset = (unsigned int*)malloc((vertex_count + 1) * sizeof(unsigned int));
	if (!mesh->vert_edge_offset) { free(vert_degree); goto fail; }

	mesh->vert_edge_offset[0] = 0;
	for (vi = 0; vi < vertex_count; vi++)
	{
		mesh->vert_edge_offset[vi + 1] = mesh->vert_edge_offset[vi] + vert_degree[vi];
		total_vert_edges += vert_degree[vi];
	}

	mesh->vert_edges = (unsigned int*)malloc(total_vert_edges * sizeof(unsigned int));
	if (!mesh->vert_edges) { free(vert_degree); goto fail; }

	memset(vert_degree, 0, vertex_count * sizeof(unsigned int));
	for (ei = 0; ei < edge_count; ei++)
	{
		unsigned int va = mesh->edges[ei].v[0];
		unsigned int vb = mesh->edges[ei].v[1];
		mesh->vert_edges[mesh->vert_edge_offset[va] + vert_degree[va]++] = ei;
		mesh->vert_edges[mesh->vert_edge_offset[vb] + vert_degree[vb]++] = ei;
	}

	free(vert_degree);
	return QAWS_STATUS_OK;

fail:
	free(hash_table);
	free(hash_next);
	free(hash_keys);
	free(hash_edge_idx);
	qaws_internal_mesh_destroy(mesh);
	return QAWS_STATUS_ALLOCATION_FAILURE;
}

void qaws_internal_mesh_destroy(qaws_internal_mesh* mesh)
{
	if (!mesh) return;
	free(mesh->tris);
	free(mesh->edges);
	free(mesh->vert_edge_offset);
	free(mesh->vert_edges);
	memset(mesh, 0, sizeof(*mesh));
}

/* ------------------------------------------------------------------ */
/*  Edge length                                                        */
/* ------------------------------------------------------------------ */

qaws_scalar qaws_internal_mesh_edge_length(
	qaws_internal_mesh const* mesh,
	unsigned int v0,
	unsigned int v1)
{
	qaws_scalar dx = mesh->positions[v1 * 3 + 0] - mesh->positions[v0 * 3 + 0];
	qaws_scalar dy = mesh->positions[v1 * 3 + 1] - mesh->positions[v0 * 3 + 1];
	qaws_scalar dz = mesh->positions[v1 * 3 + 2] - mesh->positions[v0 * 3 + 2];
	return QAWS_SQRT(dx * dx + dy * dy + dz * dz);
}

/* ------------------------------------------------------------------ */
/*  Dijkstra shortest path                                             */
/* ------------------------------------------------------------------ */

/* Simple min-heap for Dijkstra */
typedef struct dijk_entry {
	qaws_scalar dist;
	unsigned int vertex;
} dijk_entry;

static void heap_swap(dijk_entry* a, dijk_entry* b)
{
	dijk_entry tmp = *a;
	*a = *b;
	*b = tmp;
}

static void heap_push(dijk_entry* heap, unsigned int* size, dijk_entry e)
{
	unsigned int i = (*size)++;
	heap[i] = e;
	while (i > 0)
	{
		unsigned int parent = (i - 1) / 2;
		if (heap[parent].dist <= heap[i].dist) break;
		heap_swap(&heap[parent], &heap[i]);
		i = parent;
	}
}

static dijk_entry heap_pop(dijk_entry* heap, unsigned int* size)
{
	dijk_entry top = heap[0];
	(*size)--;
	if (*size > 0)
	{
		unsigned int i = 0;
		heap[0] = heap[*size];
		for (;;)
		{
			unsigned int left = 2 * i + 1;
			unsigned int right = 2 * i + 2;
			unsigned int smallest = i;
			if (left < *size && heap[left].dist < heap[smallest].dist)
				smallest = left;
			if (right < *size && heap[right].dist < heap[smallest].dist)
				smallest = right;
			if (smallest == i) break;
			heap_swap(&heap[i], &heap[smallest]);
			i = smallest;
		}
	}
	return top;
}

qaws_status qaws_internal_mesh_dijkstra(
	qaws_internal_mesh const* mesh,
	unsigned int src,
	unsigned int dst,
	unsigned int* out_path,
	unsigned int* out_path_count)
{
	unsigned int n = mesh->vertex_count;
	qaws_scalar* dist = NULL;
	unsigned int* prev = NULL;
	int* visited = NULL;
	dijk_entry* heap = NULL;
	unsigned int heap_size = 0;
	unsigned int vi, count;
	qaws_status status = QAWS_STATUS_OK;

	dist = (qaws_scalar*)malloc(n * sizeof(qaws_scalar));
	prev = (unsigned int*)malloc(n * sizeof(unsigned int));
	visited = (int*)calloc(n, sizeof(int));
	heap = (dijk_entry*)malloc(n * 4 * sizeof(dijk_entry));
	if (!dist || !prev || !visited || !heap)
	{
		status = QAWS_STATUS_ALLOCATION_FAILURE;
		goto cleanup;
	}

	for (vi = 0; vi < n; vi++)
	{
		dist[vi] = QAWS_LITERAL(1e30);
		prev[vi] = UINT_MAX;
	}
	dist[src] = QAWS_ZERO;

	{
		dijk_entry e;
		e.dist = QAWS_ZERO;
		e.vertex = src;
		heap_push(heap, &heap_size, e);
	}

	while (heap_size > 0)
	{
		dijk_entry cur = heap_pop(heap, &heap_size);
		unsigned int u = cur.vertex;
		unsigned int ei_start, ei_end, ei;

		if (visited[u]) continue;
		visited[u] = 1;

		if (u == dst) break;

		ei_start = mesh->vert_edge_offset[u];
		ei_end = mesh->vert_edge_offset[u + 1];

		for (ei = ei_start; ei < ei_end; ei++)
		{
			unsigned int edge_idx = mesh->vert_edges[ei];
			unsigned int neighbor = (mesh->edges[edge_idx].v[0] == u) ?
				mesh->edges[edge_idx].v[1] : mesh->edges[edge_idx].v[0];
			qaws_scalar w = qaws_internal_mesh_edge_length(mesh,
				mesh->edges[edge_idx].v[0], mesh->edges[edge_idx].v[1]);
			qaws_scalar alt = dist[u] + w;

			if (alt < dist[neighbor])
			{
				dijk_entry e;
				dist[neighbor] = alt;
				prev[neighbor] = u;
				e.dist = alt;
				e.vertex = neighbor;
				heap_push(heap, &heap_size, e);
			}
		}
	}

	if (!visited[dst])
	{
		*out_path_count = 0;
		status = QAWS_STATUS_NUMERICAL_FAILURE;
		goto cleanup;
	}

	/* Trace path backwards */
	count = 0;
	{
		unsigned int v = dst;
		while (v != UINT_MAX)
		{
			out_path[count++] = v;
			v = prev[v];
		}
	}

	/* Reverse path to get src -> dst order */
	{
		unsigned int i;
		for (i = 0; i < count / 2; i++)
		{
			unsigned int tmp = out_path[i];
			out_path[i] = out_path[count - 1 - i];
			out_path[count - 1 - i] = tmp;
		}
	}

	*out_path_count = count;

cleanup:
	free(dist);
	free(prev);
	free(visited);
	free(heap);
	return status;
}

/* ------------------------------------------------------------------ */
/*  Cotangent weight                                                   */
/* ------------------------------------------------------------------ */

static qaws_scalar cotan_angle(
	qaws_scalar const* p,
	qaws_scalar const* a,
	qaws_scalar const* b)
{
	/* cot(angle at p between edges pa and pb) */
	qaws_scalar pax = a[0] - p[0], pay = a[1] - p[1], paz = a[2] - p[2];
	qaws_scalar pbx = b[0] - p[0], pby = b[1] - p[1], pbz = b[2] - p[2];
	/* dot = pa . pb */
	qaws_scalar dot = pax * pbx + pay * pby + paz * pbz;
	/* cross magnitude = |pa x pb| */
	qaws_scalar cx = pay * pbz - paz * pby;
	qaws_scalar cy = paz * pbx - pax * pbz;
	qaws_scalar cz = pax * pby - pay * pbx;
	qaws_scalar cross_mag = QAWS_SQRT(cx * cx + cy * cy + cz * cz);

	if (cross_mag < QAWS_LITERAL(1e-12))
		return QAWS_ZERO;

	return dot / cross_mag;
}

qaws_scalar qaws_internal_mesh_cotan_weight(
	qaws_internal_mesh const* mesh,
	unsigned int edge_idx)
{
	qaws_internal_edge const* e = &mesh->edges[edge_idx];
	qaws_scalar w = QAWS_ZERO;
	qaws_scalar const* va = mesh->positions + e->v[0] * 3;
	qaws_scalar const* vb = mesh->positions + e->v[1] * 3;

	/* Sum cotangent contributions from both adjacent triangles */
	if (e->tri[0] != UINT_MAX && e->opposite[0] != UINT_MAX)
	{
		qaws_scalar const* opp = mesh->positions + e->opposite[0] * 3;
		w += cotan_angle(opp, va, vb);
	}

	if (e->tri[1] != UINT_MAX && e->opposite[1] != UINT_MAX)
	{
		qaws_scalar const* opp = mesh->positions + e->opposite[1] * 3;
		w += cotan_angle(opp, va, vb);
	}

	w *= QAWS_LITERAL(0.5);

	/* Clamp for robustness */
	if (w < QAWS_LITERAL(-1e6)) w = QAWS_LITERAL(-1e6);
	if (w > QAWS_LITERAL(1e6))  w = QAWS_LITERAL(1e6);

	return w;
}

/* ------------------------------------------------------------------ */
/*  Edge flip                                                          */
/* ------------------------------------------------------------------ */

int qaws_internal_mesh_try_flip_edge(
	qaws_internal_mesh* mesh,
	unsigned int edge_idx)
{
	qaws_internal_edge* e = &mesh->edges[edge_idx];
	unsigned int t0, t1;
	unsigned int a, b, c, d; /* a-b is the edge; c, d are opposite vertices */
	qaws_internal_tri* tri0;
	qaws_internal_tri* tri1;
	unsigned int local0, local1;
	unsigned int li;
	qaws_scalar ax, ay, bx, by, cx_2d, cy_2d, dx_2d, dy_2d;
	qaws_scalar cross1, cross2;

	/* Cannot flip boundary edges */
	if (e->tri[1] == UINT_MAX)
		return 0;

	t0 = e->tri[0];
	t1 = e->tri[1];
	a = e->v[0];
	b = e->v[1];
	c = e->opposite[0];
	d = e->opposite[1];

	if (c == UINT_MAX || d == UINT_MAX)
		return 0;

	/* Check convexity: the quad a-c-b-d must be convex.
	   Project to 2D using largest axis of the quad normal. */
	{
		qaws_scalar const* pa = mesh->positions + a * 3;
		qaws_scalar const* pb = mesh->positions + b * 3;
		qaws_scalar const* pc = mesh->positions + c * 3;
		qaws_scalar const* pd = mesh->positions + d * 3;

		/* Use UV coords for 2D check (simpler and works for our surface meshes) */
		if (mesh->uvs)
		{
			ax = mesh->uvs[a * 2]; ay = mesh->uvs[a * 2 + 1];
			bx = mesh->uvs[b * 2]; by = mesh->uvs[b * 2 + 1];
			cx_2d = mesh->uvs[c * 2]; cy_2d = mesh->uvs[c * 2 + 1];
			dx_2d = mesh->uvs[d * 2]; dy_2d = mesh->uvs[d * 2 + 1];
		}
		else
		{
			/* Fall back to XY projection */
			ax = pa[0]; ay = pa[1];
			bx = pb[0]; by = pb[1];
			cx_2d = pc[0]; cy_2d = pc[1];
			dx_2d = pd[0]; dy_2d = pd[1];
		}

		/* Quad vertices in order: c, a, d, b (or c, b, d, a).
		   Check that diagonals c-d and a-b intersect, i.e., the quad is convex.
		   Equivalently, c and d must be on opposite sides of line a-b,
		   and a and b must be on opposite sides of line c-d. */
		cross1 = (bx - ax) * (cy_2d - ay) - (by - ay) * (cx_2d - ax);
		cross2 = (bx - ax) * (dy_2d - ay) - (by - ay) * (dx_2d - ax);
		if (cross1 * cross2 >= QAWS_ZERO)
			return 0; /* Not convex */

		cross1 = (dx_2d - cx_2d) * (ay - cy_2d) - (dy_2d - cy_2d) * (ax - cx_2d);
		cross2 = (dx_2d - cx_2d) * (by - cy_2d) - (dy_2d - cy_2d) * (bx - cx_2d);
		if (cross1 * cross2 >= QAWS_ZERO)
			return 0; /* Not convex */

		(void)pa; (void)pb; (void)pc; (void)pd;
	}

	tri0 = &mesh->tris[t0];
	tri1 = &mesh->tris[t1];

	/* Find local index of the flipped edge in each triangle */
	local0 = UINT_MAX;
	local1 = UINT_MAX;
	for (li = 0; li < 3; li++)
	{
		if (tri0->edge[li] == edge_idx) local0 = li;
		if (tri1->edge[li] == edge_idx) local1 = li;
	}
	if (local0 == UINT_MAX || local1 == UINT_MAX)
		return 0;

	/* Flip: replace edge a-b with edge c-d
	   tri0 becomes (c, d, a), tri1 becomes (d, c, b)
	   But we need to preserve the edges that connect to adjacent triangles */
	{
		/* Before flip:
		   tri0 = (v[local0]=c, next=a, next=b) with edge[local0] = edge_idx (opposite c)
		   tri1 = (v[local1]=d, next=...) with edge[local1] = edge_idx (opposite d)

		   After flip, edge becomes c-d:
		   tri0 = (a, c, d) - has edges opposite a (= edge c-d = edge_idx), opposite c, opposite d
		   tri1 = (b, d, c) - has edges opposite b (= edge d-c = edge_idx), opposite d, opposite c
		*/
		unsigned int a_local0 = (local0 + 1) % 3; /* local index of a in tri0 */
		unsigned int b_local0 = (local0 + 2) % 3; /* local index of b in tri0 */
		unsigned int d_local1_next = (local1 + 1) % 3;
		unsigned int b_local1 = (local1 + 2) % 3; /* find where b is in tri1 */

		/* tri0 has vertices v[local0]=c (opposite),
		   v[(local0+1)%3] and v[(local0+2)%3] are a and b (the edge).
		   edge[local0] is opposite c = the edge a-b being flipped.
		   edge[(local0+1)%3] is opposite v[(local0+1)%3].
		   edge[(local0+2)%3] is opposite v[(local0+2)%3]. */

		unsigned int v_a = tri0->v[a_local0]; /* a or b */
		unsigned int v_b_tri0 = tri0->v[b_local0]; /* the other one */

		/* Actually, the order is:
		   tri0->v[local0] = c (the opposite vertex)
		   tri0->v[(local0+1)%3] = one of {a, b}
		   tri0->v[(local0+2)%3] = the other of {a, b}
		   tri0->edge[local0] = edge_idx (edge a-b, opposite c)
		   tri0->edge[(local0+1)%3] = edge opposite tri0->v[(local0+1)%3]
		   tri0->edge[(local0+2)%3] = edge opposite tri0->v[(local0+2)%3]

		   Similarly for tri1 with d as opposite.
		*/

		unsigned int edge_opp_a_in_t0, edge_opp_b_in_t0;
		unsigned int v1_next, v1_next2;
		unsigned int edge_opp_next_in_t1, edge_opp_next2_in_t1;
		unsigned int edge_for_new_t0_opp_c, edge_for_new_t0_opp_d;
		unsigned int edge_for_new_t1_opp_d, edge_for_new_t1_opp_b;

		/* Determine which vertex is a and which is b in tri0 */
		if (v_a == a)
		{
			/* tri0->v = [c, a, b] at indices [local0, a_local0, b_local0] */
			edge_opp_a_in_t0 = tri0->edge[a_local0]; /* edge c-b (opposite a) */
			edge_opp_b_in_t0 = tri0->edge[b_local0]; /* edge c-a (opposite b) */
		}
		else
		{
			/* tri0->v = [c, b, a] */
			/* v_a is b, v_b_tri0 is a */
			edge_opp_a_in_t0 = tri0->edge[b_local0]; /* edge opposite a */
			edge_opp_b_in_t0 = tri0->edge[a_local0]; /* edge opposite b = edge c-a */

			(void)v_b_tri0;
		}

		/* For tri1: v[local1] = d, the other two include a and b */
		v1_next = tri1->v[d_local1_next];
		v1_next2 = tri1->v[b_local1];
		edge_opp_next_in_t1 = tri1->edge[d_local1_next];
		edge_opp_next2_in_t1 = tri1->edge[b_local1];

		/* We need: edge opposite a in tri1 and edge opposite b in tri1 */
		{
			unsigned int edge_opp_a_in_t1, edge_opp_b_in_t1;
			if (v1_next == a)
			{
				edge_opp_a_in_t1 = edge_opp_next_in_t1;
				edge_opp_b_in_t1 = edge_opp_next2_in_t1;
			}
			else
			{
				edge_opp_a_in_t1 = edge_opp_next2_in_t1;
				edge_opp_b_in_t1 = edge_opp_next_in_t1;
			}

			(void)v1_next2;

			/* New triangles after flip:
			   new_tri0 = (a, c, d) with:
			     edge opposite a = edge_idx (c-d, the flipped edge)
			     edge opposite c = edge a-d (was in tri1, opposite the vertex not a/b... opposite a? no)
			                     = edge_opp_a_in_t1 (edge opposite a in old tri1, which is edge d-b... no)

			   Hmm, let me think more carefully.

			   new_tri0 = (a, c, d):
			     edge opposite a = c-d = edge_idx
			     edge opposite c = a-d: this edge was in old tri1 (which had a, b, d).
			       In old tri1, the edge a-d is opposite b. So it's edge_opp_b_in_t1.
			     edge opposite d = a-c: this edge was in old tri0 (which had a, b, c).
			       In old tri0, the edge a-c is opposite b. So it's edge_opp_b_in_t0.

			   new_tri1 = (b, d, c):
			     edge opposite b = d-c = edge_idx
			     edge opposite d = b-c: this was in old tri0 (a, b, c).
			       In old tri0, edge b-c is opposite a. So it's edge_opp_a_in_t0.
			     edge opposite c = b-d: this was in old tri1 (a, b, d).
			       In old tri1, edge b-d is opposite a. So it's edge_opp_a_in_t1.
			*/

			/* Set new tri0 = (a, c, d) */
			tri0->v[0] = a;
			tri0->v[1] = c;
			tri0->v[2] = d;
			tri0->edge[0] = edge_idx;
			edge_for_new_t0_opp_c = edge_opp_b_in_t1;
			edge_for_new_t0_opp_d = edge_opp_b_in_t0;
			tri0->edge[1] = edge_for_new_t0_opp_c;
			tri0->edge[2] = edge_for_new_t0_opp_d;

			/* Set new tri1 = (b, d, c) */
			tri1->v[0] = b;
			tri1->v[1] = d;
			tri1->v[2] = c;
			tri1->edge[0] = edge_idx;
			edge_for_new_t1_opp_d = edge_opp_a_in_t0;
			edge_for_new_t1_opp_b = edge_opp_a_in_t1;
			tri1->edge[1] = edge_for_new_t1_opp_d;
			tri1->edge[2] = edge_for_new_t1_opp_b;

			/* Update the flipped edge itself */
			e->v[0] = c < d ? c : d;
			e->v[1] = c < d ? d : c;
			e->opposite[0] = a; /* a is opposite in t0 */
			e->opposite[1] = b; /* b is opposite in t1 */
			/* tri references stay t0, t1 */

			/* Update neighbor edges that moved between triangles.
			   edge_for_new_t0_opp_c (was in tri1) now belongs to tri0.
			   edge_for_new_t1_opp_d (was in tri0) now belongs to tri1. */
			{
				qaws_internal_edge* moved_to_t0 = &mesh->edges[edge_for_new_t0_opp_c];
				if (moved_to_t0->tri[0] == t1) moved_to_t0->tri[0] = t0;
				else if (moved_to_t0->tri[1] == t1) moved_to_t0->tri[1] = t0;
			}
			{
				qaws_internal_edge* moved_to_t1 = &mesh->edges[edge_for_new_t1_opp_d];
				if (moved_to_t1->tri[0] == t0) moved_to_t1->tri[0] = t1;
				else if (moved_to_t1->tri[1] == t0) moved_to_t1->tri[1] = t1;
			}

			/* Update opposite vertex references for the moved edges */
			{
				qaws_internal_edge* me = &mesh->edges[edge_for_new_t0_opp_c];
				if (me->tri[0] == t0) me->opposite[0] = c;
				if (me->tri[1] == t0) me->opposite[1] = c;
			}
			{
				qaws_internal_edge* me = &mesh->edges[edge_for_new_t1_opp_d];
				if (me->tri[0] == t1) me->opposite[0] = d;
				if (me->tri[1] == t1) me->opposite[1] = d;
			}
		}
	}

	/* Note: vert_edge adjacency is not updated after flip for simplicity.
	   The flip-based geodesic algorithm re-builds as needed. */

	return 1;
}

/* ------------------------------------------------------------------ */
/*  Closest vertex by UV                                               */
/* ------------------------------------------------------------------ */

unsigned int qaws_internal_mesh_closest_vertex_uv(
	qaws_internal_mesh const* mesh,
	qaws_scalar u,
	qaws_scalar v)
{
	unsigned int best = 0;
	qaws_scalar best_dist = QAWS_LITERAL(1e30);
	unsigned int vi;

	if (!mesh->uvs)
	{
		/* Fall back: use XY of position */
		for (vi = 0; vi < mesh->vertex_count; vi++)
		{
			qaws_scalar du = mesh->positions[vi * 3 + 0] - u;
			qaws_scalar dv = mesh->positions[vi * 3 + 1] - v;
			qaws_scalar d = du * du + dv * dv;
			if (d < best_dist)
			{
				best_dist = d;
				best = vi;
			}
		}
		return best;
	}

	for (vi = 0; vi < mesh->vertex_count; vi++)
	{
		qaws_scalar du = mesh->uvs[vi * 2 + 0] - u;
		qaws_scalar dv = mesh->uvs[vi * 2 + 1] - v;
		qaws_scalar d = du * du + dv * dv;
		if (d < best_dist)
		{
			best_dist = d;
			best = vi;
		}
	}
	return best;
}
