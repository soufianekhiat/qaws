#include "qaws_brep.h"
#include "qaws_surface.h"
#include "qaws_inspect.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"

/* ===================================================================
 * Internal data structures
 * =================================================================== */

#define BREP_INITIAL_CAPACITY 8
#define BREP_VOLUME_SAMPLES 20

typedef struct brep_vertex
{
	qaws_vec3 position;
} brep_vertex;

typedef struct brep_edge
{
	qaws_brep_vertex_id v_start;
	qaws_brep_vertex_id v_end;
	qaws_curve const* curve_3d;
	qaws_brep_face_id face_left;    /* QAWS_BREP_INVALID_ID if boundary */
	qaws_brep_face_id face_right;   /* QAWS_BREP_INVALID_ID if boundary */
} brep_edge;

typedef struct brep_face
{
	qaws_surface const* surface;
	int is_reversed;
	qaws_brep_edge_id* edge_ids;
	int* edge_orientations;
	unsigned int edge_count;
} brep_face;

struct qaws_brep_shell
{
	brep_vertex* vertices;
	unsigned int vertex_count;
	unsigned int vertex_capacity;

	brep_edge* edges;
	unsigned int edge_count;
	unsigned int edge_capacity;

	brep_face* faces;
	unsigned int face_count;
	unsigned int face_capacity;
};

/* ===================================================================
 * Helpers
 * =================================================================== */

static int brep_grow_vertices(qaws_brep_shell* shell)
{
	unsigned int new_cap;
	brep_vertex* new_buf;

	if (shell->vertex_count < shell->vertex_capacity)
		return 1;

	new_cap = shell->vertex_capacity == 0
		? BREP_INITIAL_CAPACITY
		: shell->vertex_capacity * 2;
	new_buf = (brep_vertex*)realloc(
		shell->vertices, new_cap * sizeof(brep_vertex));
	if (!new_buf)
		return 0;

	shell->vertices = new_buf;
	shell->vertex_capacity = new_cap;
	return 1;
}

static int brep_grow_edges(qaws_brep_shell* shell)
{
	unsigned int new_cap;
	brep_edge* new_buf;

	if (shell->edge_count < shell->edge_capacity)
		return 1;

	new_cap = shell->edge_capacity == 0
		? BREP_INITIAL_CAPACITY
		: shell->edge_capacity * 2;
	new_buf = (brep_edge*)realloc(
		shell->edges, new_cap * sizeof(brep_edge));
	if (!new_buf)
		return 0;

	shell->edges = new_buf;
	shell->edge_capacity = new_cap;
	return 1;
}

static int brep_grow_faces(qaws_brep_shell* shell)
{
	unsigned int new_cap;
	brep_face* new_buf;

	if (shell->face_count < shell->face_capacity)
		return 1;

	new_cap = shell->face_capacity == 0
		? BREP_INITIAL_CAPACITY
		: shell->face_capacity * 2;
	new_buf = (brep_face*)realloc(
		shell->faces, new_cap * sizeof(brep_face));
	if (!new_buf)
		return 0;

	shell->faces = new_buf;
	shell->face_capacity = new_cap;
	return 1;
}

/* ===================================================================
 * Creation / destruction
 * =================================================================== */

qaws_status qaws_brep_create(qaws_brep_shell** out_shell)
{
	qaws_brep_shell* shell;

	if (!out_shell)
		return QAWS_STATUS_INVALID_ARGUMENT;

	shell = (qaws_brep_shell*)calloc(1, sizeof(qaws_brep_shell));
	if (!shell)
		return QAWS_STATUS_ALLOCATION_FAILURE;

	*out_shell = shell;
	return QAWS_STATUS_OK;
}

void qaws_brep_destroy(qaws_brep_shell* shell)
{
	unsigned int i;

	if (!shell)
		return;

	for (i = 0; i < shell->face_count; i++) {
		free(shell->faces[i].edge_ids);
		free(shell->faces[i].edge_orientations);
	}

	free(shell->vertices);
	free(shell->edges);
	free(shell->faces);
	free(shell);
}

/* ===================================================================
 * Add topological entities
 * =================================================================== */

qaws_status qaws_brep_add_vertex(
	qaws_brep_shell* shell,
	qaws_vec3 position,
	qaws_brep_vertex_id* out_id)
{
	brep_vertex* v;

	if (!shell || !out_id)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (!brep_grow_vertices(shell))
		return QAWS_STATUS_ALLOCATION_FAILURE;

	v = &shell->vertices[shell->vertex_count];
	v->position = position;

	*out_id = shell->vertex_count;
	shell->vertex_count++;
	return QAWS_STATUS_OK;
}

qaws_status qaws_brep_add_edge(
	qaws_brep_shell* shell,
	qaws_brep_edge_desc const* desc,
	qaws_brep_edge_id* out_id)
{
	brep_edge* e;

	if (!shell || !desc || !out_id)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (desc->v_start >= shell->vertex_count ||
		desc->v_end >= shell->vertex_count)
		return QAWS_STATUS_INVALID_ARGUMENT;

	if (!brep_grow_edges(shell))
		return QAWS_STATUS_ALLOCATION_FAILURE;

	e = &shell->edges[shell->edge_count];
	e->v_start = desc->v_start;
	e->v_end = desc->v_end;
	e->curve_3d = desc->curve_3d;
	e->face_left = QAWS_BREP_INVALID_ID;
	e->face_right = QAWS_BREP_INVALID_ID;

	*out_id = shell->edge_count;
	shell->edge_count++;
	return QAWS_STATUS_OK;
}

qaws_status qaws_brep_add_face(
	qaws_brep_shell* shell,
	qaws_brep_face_desc const* desc,
	qaws_brep_edge_id const* edge_ids,
	int const* edge_orientations,
	unsigned int edge_count,
	qaws_brep_face_id* out_id)
{
	brep_face* f;
	qaws_brep_face_id face_id;
	unsigned int i;

	if (!shell || !desc || !out_id)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (edge_count > 0 && (!edge_ids || !edge_orientations))
		return QAWS_STATUS_INVALID_ARGUMENT;

	/* Validate edge ids */
	for (i = 0; i < edge_count; i++) {
		if (edge_ids[i] >= shell->edge_count)
			return QAWS_STATUS_INVALID_ARGUMENT;
		if (edge_orientations[i] != 1 && edge_orientations[i] != -1)
			return QAWS_STATUS_INVALID_ARGUMENT;
	}

	/* Check that assigning this face to the edges won't violate manifold */
	face_id = shell->face_count;
	for (i = 0; i < edge_count; i++) {
		brep_edge* e = &shell->edges[edge_ids[i]];
		if (edge_orientations[i] == 1) {
			if (e->face_left != QAWS_BREP_INVALID_ID)
				return QAWS_STATUS_INVALID_ARGUMENT;
		} else {
			if (e->face_right != QAWS_BREP_INVALID_ID)
				return QAWS_STATUS_INVALID_ARGUMENT;
		}
	}

	if (!brep_grow_faces(shell))
		return QAWS_STATUS_ALLOCATION_FAILURE;

	f = &shell->faces[face_id];
	f->surface = desc->surface;
	f->is_reversed = desc->is_reversed;
	f->edge_count = edge_count;
	f->edge_ids = NULL;
	f->edge_orientations = NULL;

	if (edge_count > 0) {
		f->edge_ids = (qaws_brep_edge_id*)malloc(
			edge_count * sizeof(qaws_brep_edge_id));
		f->edge_orientations = (int*)malloc(
			edge_count * sizeof(int));
		if (!f->edge_ids || !f->edge_orientations) {
			free(f->edge_ids);
			free(f->edge_orientations);
			f->edge_ids = NULL;
			f->edge_orientations = NULL;
			return QAWS_STATUS_ALLOCATION_FAILURE;
		}
		memcpy(f->edge_ids, edge_ids,
			edge_count * sizeof(qaws_brep_edge_id));
		memcpy(f->edge_orientations, edge_orientations,
			edge_count * sizeof(int));
	}

	/* Update edge face adjacency */
	for (i = 0; i < edge_count; i++) {
		brep_edge* e = &shell->edges[edge_ids[i]];
		if (edge_orientations[i] == 1)
			e->face_left = face_id;
		else
			e->face_right = face_id;
	}

	*out_id = face_id;
	shell->face_count++;
	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Query topology
 * =================================================================== */

unsigned int qaws_brep_get_vertex_count(qaws_brep_shell const* shell)
{
	return shell ? shell->vertex_count : 0;
}

unsigned int qaws_brep_get_edge_count(qaws_brep_shell const* shell)
{
	return shell ? shell->edge_count : 0;
}

unsigned int qaws_brep_get_face_count(qaws_brep_shell const* shell)
{
	return shell ? shell->face_count : 0;
}

qaws_status qaws_brep_get_vertex_position(
	qaws_brep_shell const* shell,
	qaws_brep_vertex_id id,
	qaws_vec3* out_position)
{
	if (!shell || !out_position)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (id >= shell->vertex_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	*out_position = shell->vertices[id].position;
	return QAWS_STATUS_OK;
}

qaws_status qaws_brep_get_face_surface(
	qaws_brep_shell const* shell,
	qaws_brep_face_id id,
	qaws_surface const** out_surface)
{
	if (!shell || !out_surface)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (id >= shell->face_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	*out_surface = shell->faces[id].surface;
	return QAWS_STATUS_OK;
}

qaws_status qaws_brep_get_edge_vertices(
	qaws_brep_shell const* shell,
	qaws_brep_edge_id id,
	qaws_brep_vertex_id* out_v_start,
	qaws_brep_vertex_id* out_v_end)
{
	if (!shell || !out_v_start || !out_v_end)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (id >= shell->edge_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	*out_v_start = shell->edges[id].v_start;
	*out_v_end = shell->edges[id].v_end;
	return QAWS_STATUS_OK;
}

qaws_status qaws_brep_get_edge_faces(
	qaws_brep_shell const* shell,
	qaws_brep_edge_id id,
	qaws_brep_face_id* out_face_left,
	qaws_brep_face_id* out_face_right)
{
	if (!shell || !out_face_left || !out_face_right)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (id >= shell->edge_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	*out_face_left = shell->edges[id].face_left;
	*out_face_right = shell->edges[id].face_right;
	return QAWS_STATUS_OK;
}

qaws_status qaws_brep_get_face_edges(
	qaws_brep_shell const* shell,
	qaws_brep_face_id id,
	qaws_brep_edge_id* out_edges,
	unsigned int capacity,
	unsigned int* out_count)
{
	brep_face const* f;

	if (!shell || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (id >= shell->face_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	f = &shell->faces[id];
	*out_count = f->edge_count;

	if (out_edges) {
		if (capacity < f->edge_count)
			return QAWS_STATUS_BUFFER_TOO_SMALL;
		memcpy(out_edges, f->edge_ids,
			f->edge_count * sizeof(qaws_brep_edge_id));
	}

	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Adjacency queries
 * =================================================================== */

qaws_status qaws_brep_get_vertex_edges(
	qaws_brep_shell const* shell,
	qaws_brep_vertex_id id,
	qaws_brep_edge_id* out_edges,
	unsigned int capacity,
	unsigned int* out_count)
{
	unsigned int i;
	unsigned int count;

	if (!shell || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (id >= shell->vertex_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	count = 0;
	for (i = 0; i < shell->edge_count; i++) {
		if (shell->edges[i].v_start == id ||
			shell->edges[i].v_end == id) {
			if (out_edges) {
				if (count >= capacity)
					return QAWS_STATUS_BUFFER_TOO_SMALL;
				out_edges[count] = i;
			}
			count++;
		}
	}

	*out_count = count;
	return QAWS_STATUS_OK;
}

qaws_status qaws_brep_get_vertex_faces(
	qaws_brep_shell const* shell,
	qaws_brep_vertex_id id,
	qaws_brep_face_id* out_faces,
	unsigned int capacity,
	unsigned int* out_count)
{
	unsigned int i;
	unsigned int j;
	unsigned int count;

	if (!shell || !out_count)
		return QAWS_STATUS_INVALID_ARGUMENT;
	if (id >= shell->vertex_count)
		return QAWS_STATUS_OUT_OF_RANGE;

	count = 0;

	/* Scan edges touching this vertex, collect their faces (deduplicated) */
	for (i = 0; i < shell->edge_count; i++) {
		brep_edge const* e = &shell->edges[i];
		if (e->v_start != id && e->v_end != id)
			continue;

		/* Try to add face_left */
		if (e->face_left != QAWS_BREP_INVALID_ID) {
			int found = 0;
			for (j = 0; j < count; j++) {
				if (out_faces && out_faces[j] == e->face_left) {
					found = 1;
					break;
				}
			}
			if (!found) {
				if (out_faces) {
					if (count >= capacity)
						return QAWS_STATUS_BUFFER_TOO_SMALL;
					out_faces[count] = e->face_left;
				}
				count++;
			}
		}

		/* Try to add face_right */
		if (e->face_right != QAWS_BREP_INVALID_ID) {
			int found = 0;
			for (j = 0; j < count; j++) {
				if (out_faces && out_faces[j] == e->face_right) {
					found = 1;
					break;
				}
			}
			if (!found) {
				if (out_faces) {
					if (count >= capacity)
						return QAWS_STATUS_BUFFER_TOO_SMALL;
					out_faces[count] = e->face_right;
				}
				count++;
			}
		}
	}

	*out_count = count;
	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Validation
 * =================================================================== */

qaws_status qaws_brep_validate(
	qaws_brep_shell const* shell,
	unsigned int* out_flags)
{
	unsigned int flags;
	unsigned int i;

	if (!shell || !out_flags)
		return QAWS_STATUS_INVALID_ARGUMENT;

	flags = QAWS_BREP_VALID;

	for (i = 0; i < shell->edge_count; i++) {
		brep_edge const* e = &shell->edges[i];

		/* Check for boundary edges (not closed) */
		if (e->face_left == QAWS_BREP_INVALID_ID ||
			e->face_right == QAWS_BREP_INVALID_ID) {
			flags |= QAWS_BREP_NOT_CLOSED;
		}

		/* Check non-manifold: both faces same (degenerate adjacency) */
		if (e->face_left != QAWS_BREP_INVALID_ID &&
			e->face_right != QAWS_BREP_INVALID_ID &&
			e->face_left == e->face_right) {
			flags |= QAWS_BREP_NOT_MANIFOLD;
		}

		/* Check degenerate edge: zero length */
		if (e->v_start < shell->vertex_count &&
			e->v_end < shell->vertex_count) {
			qaws_vec3 const* a = &shell->vertices[e->v_start].position;
			qaws_vec3 const* b = &shell->vertices[e->v_end].position;
			qaws_scalar dx = b->x - a->x;
			qaws_scalar dy = b->y - a->y;
			qaws_scalar dz = b->z - a->z;
			qaws_scalar len2 = dx * dx + dy * dy + dz * dz;
			if (len2 < QAWS_EPSILON * QAWS_EPSILON)
				flags |= QAWS_BREP_HAS_DEGENERATE;
		}
	}

	/* Check orientation consistency: for each edge shared by two faces,
	   the two faces should traverse the edge in opposite directions.
	   Since face_left is set when orientation == +1 and face_right when
	   orientation == -1, a properly oriented manifold edge always has
	   one face on each side by construction. We verify that no edge
	   is referenced by the same face on both sides (already checked above),
	   and additionally verify the stored orientations are consistent. */
	for (i = 0; i < shell->edge_count; i++) {
		brep_edge const* e = &shell->edges[i];

		if (e->face_left == QAWS_BREP_INVALID_ID ||
			e->face_right == QAWS_BREP_INVALID_ID)
			continue;

		/* For a properly oriented shell, the left face uses orientation +1
		   and the right face uses orientation -1. We verify this by checking
		   the face edge lists. */
		{
			brep_face const* fl = &shell->faces[e->face_left];
			brep_face const* fr = &shell->faces[e->face_right];
			int left_ori = 0;
			int right_ori = 0;
			unsigned int j;

			for (j = 0; j < fl->edge_count; j++) {
				if (fl->edge_ids[j] == i) {
					left_ori = fl->edge_orientations[j];
					break;
				}
			}
			for (j = 0; j < fr->edge_count; j++) {
				if (fr->edge_ids[j] == i) {
					right_ori = fr->edge_orientations[j];
					break;
				}
			}

			/* Consistent orientation: opposite signs */
			if (left_ori != 0 && right_ori != 0 &&
				left_ori == right_ori) {
				flags |= QAWS_BREP_NOT_ORIENTABLE;
			}
		}
	}

	/* Check degenerate faces: face with fewer than 3 edges */
	for (i = 0; i < shell->face_count; i++) {
		if (shell->faces[i].edge_count < 3)
			flags |= QAWS_BREP_HAS_DEGENERATE;
	}

	*out_flags = flags;
	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Euler characteristic
 * =================================================================== */

qaws_status qaws_brep_euler_characteristic(
	qaws_brep_shell const* shell,
	int* out_chi)
{
	if (!shell || !out_chi)
		return QAWS_STATUS_INVALID_ARGUMENT;

	*out_chi = (int)shell->vertex_count
		- (int)shell->edge_count
		+ (int)shell->face_count;
	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Volume computation (divergence theorem)
 * =================================================================== */

qaws_status qaws_brep_compute_volume(
	qaws_brep_shell const* shell,
	qaws_scalar* out_volume)
{
	unsigned int fi;
	unsigned int i;
	unsigned int j;
	unsigned int n;
	qaws_scalar total_volume;

	if (!shell || !out_volume)
		return QAWS_STATUS_INVALID_ARGUMENT;

	n = BREP_VOLUME_SAMPLES;
	total_volume = QAWS_ZERO;

	for (fi = 0; fi < shell->face_count; fi++) {
		brep_face const* face = &shell->faces[fi];
		qaws_surface const* surf = face->surface;
		qaws_range u_range;
		qaws_range v_range;
		qaws_scalar u_min, u_max, v_min, v_max;
		qaws_scalar du, dv;
		qaws_scalar face_vol;
		qaws_scalar sign;

		if (!surf)
			continue;

		u_range = qaws_surface_get_u_range(surf);
		v_range = qaws_surface_get_v_range(surf);
		u_min = u_range.min_value;
		u_max = u_range.max_value;
		v_min = v_range.min_value;
		v_max = v_range.max_value;
		du = (u_max - u_min) / (qaws_scalar)n;
		dv = (v_max - v_min) / (qaws_scalar)n;

		(void)sign;
		face_vol = QAWS_ZERO;

		/* Use the divergence theorem: V = (1/3) integral x.n dA.
		   For each grid cell, evaluate position and normal, and sum
		   x*nx + y*ny + z*nz weighted by the surface element area.
		   This correctly handles arbitrary surface orientations. */
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				qaws_scalar u_c = u_min + ((qaws_scalar)i + QAWS_LITERAL(0.5)) * du;
				qaws_scalar v_c = v_min + ((qaws_scalar)j + QAWS_LITERAL(0.5)) * dv;
				qaws_surface_eval_result rc;
				qaws_status s;
				qaws_scalar nx, ny, nz, area_elem, dot_val;

				memset(&rc, 0, sizeof(rc));
				s = qaws_surface_evaluate(surf, u_c, v_c,
					QAWS_SURFACE_EVAL_POSITION | QAWS_SURFACE_EVAL_DU | QAWS_SURFACE_EVAL_DV,
					&rc);
				if (s != QAWS_STATUS_OK) continue;

				/* Surface element normal (unnormalized = du x dv, length = area element) */
				nx = rc.du.y * rc.dv.z - rc.du.z * rc.dv.y;
				ny = rc.du.z * rc.dv.x - rc.du.x * rc.dv.z;
				nz = rc.du.x * rc.dv.y - rc.du.y * rc.dv.x;

				/* Flip for reversed faces */
				if (face->is_reversed) {
					nx = -nx; ny = -ny; nz = -nz;
				}

				/* Divergence theorem: dV = (1/3) * (x*nx + y*ny + z*nz) * du * dv */
				dot_val = rc.position.x * nx + rc.position.y * ny + rc.position.z * nz;
				area_elem = du * dv;
				face_vol += dot_val * area_elem / QAWS_LITERAL(3.0);
			}
		}

		total_volume += face_vol;
	}

	*out_volume = QAWS_FABS(total_volume);
	return QAWS_STATUS_OK;
}

/* ===================================================================
 * Surface area computation
 * =================================================================== */

qaws_status qaws_brep_compute_surface_area(
	qaws_brep_shell const* shell,
	qaws_scalar* out_area)
{
	unsigned int i;
	qaws_scalar total_area;

	if (!shell || !out_area)
		return QAWS_STATUS_INVALID_ARGUMENT;

	total_area = QAWS_ZERO;

	for (i = 0; i < shell->face_count; i++) {
		qaws_surface const* surf = shell->faces[i].surface;
		qaws_scalar face_area;
		qaws_status s;

		if (!surf)
			continue;

		s = qaws_surface_compute_area(surf, &face_area);
		if (s != QAWS_STATUS_OK)
			return s;

		total_area += face_area;
	}

	*out_area = total_area;
	return QAWS_STATUS_OK;
}
