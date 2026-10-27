/*
 * Meshes from a host's buffers (qaws_mesh_view.h): validation, an owned mesh
 * made in one pass, a frame's positions set without allocating.
 */

#include "../internal/qaws_internal_memory.h"
#include "../internal/qaws_internal_array.h"
#include "qaws_mesh_internal.h"
#include <stdint.h>
#include <string.h>
#include <math.h>

unsigned int qaws_internal_mesh_view_index(qaws_mesh_view const* v, unsigned int t, unsigned int k)
{
	size_t is = v->index_type == QAWS_INDEX_U16 ? 2 : 4, stride = v->index_stride ? v->index_stride : 3 * is;
	unsigned char const* p = (unsigned char const*)v->indices + (size_t)t * stride + k * is;
	if (is == 2)
	{
		uint16_t s;
		memcpy(&s, p, sizeof s);
		return s;
	}
	{
		uint32_t s;
		memcpy(&s, p, sizeof s);
		return s;
	}
}

int qaws_internal_mesh_view_positions(double* x, qaws_array_view const* positions)
{
	unsigned int i, k;
	for (i = 0; i < positions->count; i++)
		for (k = 0; k < 3; k++)
			if (!isfinite(x[3 * (size_t)i + k] = qaws_internal_array_get(positions, i, k)))
				return 0;
	return 1;
}

qaws_status qaws_mesh_view_validate(qaws_mesh_view const* view)
{
	unsigned int t, k, i;
	if (!view || qaws_array_view_validate(&view->positions) != QAWS_STATUS_OK || view->positions.components != 3 ||
		(view->index_type != QAWS_INDEX_U16 && view->index_type != QAWS_INDEX_U32) || (view->triangle_count && !view->indices) ||
		(view->index_stride && view->index_stride < 3 * (size_t)(view->index_type == QAWS_INDEX_U16 ? 2 : 4)))
		return QAWS_STATUS_INVALID_ARGUMENT;
	for (i = 0; i < view->positions.count; i++)
		for (k = 0; k < 3; k++)
			if (!isfinite(qaws_internal_array_get(&view->positions, i, k)))
				return QAWS_STATUS_INVALID_ARGUMENT;
	for (t = 0; t < view->triangle_count; t++)
		for (k = 0; k < 3; k++)
			if (qaws_internal_mesh_view_index(view, t, k) >= view->positions.count)
				return QAWS_STATUS_INVALID_ARGUMENT;
	return QAWS_STATUS_OK;
}

qaws_status qaws_mesh_3d_create_from_view(qaws_mesh_3d** out_mesh, qaws_mesh_view const* view)
{
	double* x;
	unsigned int* f;
	unsigned int t, k;
	qaws_status st;
	if (!out_mesh)
		return QAWS_STATUS_INVALID_ARGUMENT;
	*out_mesh = NULL;
	if ((st = qaws_mesh_view_validate(view)) != QAWS_STATUS_OK)
		return st;
	x = (double*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(double), 3, (size_t)view->positions.count + 1));
	f = (unsigned int*)qaws_internal_scratch_malloc(qaws_internal_mul3(sizeof(unsigned int), 3, (size_t)view->triangle_count + 1));
	if (!x || !f)
		st = QAWS_STATUS_ALLOCATION_FAILURE;
	else
	{
		qaws_internal_mesh_view_positions(x, &view->positions);
		for (t = 0; t < view->triangle_count; t++)
			for (k = 0; k < 3; k++)
				f[3 * (size_t)t + k] = qaws_internal_mesh_view_index(view, t, k);
		st = hm_mesh_3d_create(x, view->positions.count, f, NULL, view->triangle_count, out_mesh);
	}
	qaws_internal_free(x);
	qaws_internal_free(f);
	return st;
}

qaws_status qaws_mesh_3d_set_positions_from_view(qaws_mesh_3d* mesh, qaws_array_view const* positions)
{
	unsigned int i, k;
	if (!mesh || !positions || qaws_array_view_validate(positions) != QAWS_STATUS_OK || positions->components != 3 || positions->count != mesh->hm.nv)
		return QAWS_STATUS_INVALID_ARGUMENT;
	/* (checked first: a value not finite leaves the mesh as it was) */
	for (i = 0; i < positions->count; i++)
		for (k = 0; k < 3; k++)
			if (!isfinite(qaws_internal_array_get(positions, i, k)))
				return QAWS_STATUS_INVALID_ARGUMENT;
	qaws_internal_mesh_view_positions(mesh->hm.x, positions);
	return QAWS_STATUS_OK;
}
