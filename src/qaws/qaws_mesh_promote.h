#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_PROMOTE_H_F32
#define QAWS_MESH_PROMOTE_H_F32
#define QAWS_MESH_PROMOTE_H_BODY
#endif
#else
#ifndef QAWS_MESH_PROMOTE_H_F64
#define QAWS_MESH_PROMOTE_H_F64
#define QAWS_MESH_PROMOTE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_PROMOTE_H_BODY
#undef QAWS_MESH_PROMOTE_H_BODY
#define QAWS_MESH_PROMOTE_H

#include "qaws_mesh.h"

/*
 * Attribute promotion (as Houdini's Attribute Promote): an attribute moved
 * from one element to another, the values meeting at each target combined.
 *
 * Which source elements meet a target element:
 *   to a vertex   its corners, the faces about it, the edges at it;
 *   to a corner   its vertex, its face, its two edges on that face;
 *   to a face     its vertices, its corners, its edges;
 *   to an edge    its two vertices, the corners at its ends on its faces, its
 *                 faces.
 * With a piece attribute (one NEAREST component on the target element) the
 * target elements of equal value form a piece: the method runs over all the
 * sources of its elements and the result goes to every one of them; the
 * source and target element may then be the same.
 *
 * Methods (per component, but AVERAGE, which follows the attribute's rule:
 * a NORMAL renormalized, colours averaged in their space, quaternions on one
 * side): MAXIMUM, MINIMUM, AVERAGE, MODE (the most common value, a tie the
 * least), MEDIAN (an even count: the higher middle), SUM, SUM_OF_SQUARES,
 * RMS, FIRST and LAST (the source of least and greatest index). The index
 * attribute, when named, receives the source element that gave the value
 * (MAXIMUM and MINIMUM: the first component's; MODE, FIRST, LAST; others
 * none, -1). An element no source meets (an isolated vertex) takes zeros.
 */

#ifndef QAWS_ENUM_QAWS_PROMOTE_METHOD
#define QAWS_ENUM_QAWS_PROMOTE_METHOD
typedef enum qaws_promote_method
{
	QAWS_PROMOTE_MAXIMUM = 0,
	QAWS_PROMOTE_MINIMUM,
	QAWS_PROMOTE_AVERAGE,
	QAWS_PROMOTE_MODE,
	QAWS_PROMOTE_MEDIAN,
	QAWS_PROMOTE_SUM,
	QAWS_PROMOTE_SUM_OF_SQUARES,
	QAWS_PROMOTE_RMS,
	QAWS_PROMOTE_FIRST,
	QAWS_PROMOTE_LAST
} qaws_promote_method;
#endif

typedef struct qaws_attribute_promote_desc
{
	char const* name;                   /* the attribute promoted */
	qaws_mesh_element element;          /* the element it goes to */
	qaws_promote_method method;
	char const* new_name;               /* NULL: name (the original then replaced) */
	int delete_original;                /* nonzero: the original removed (always when the names are one) */
	char const* piece;                  /* NULL, or a 1-component attribute of the target element: its equal values pieces */
	char const* index_name;             /* NULL, or an attribute to receive the index of the source that gave each value */
} qaws_attribute_promote_desc;

/* INVALID_ARGUMENT for no such attribute, an unknown element or method, the same element without a piece, a piece that is
   no 1-component attribute of the target element, a new name another attribute holds */
qaws_status qaws_mesh_3d_promote_attribute(qaws_mesh_3d* mesh, qaws_attribute_promote_desc const* desc);

/* The whole mesh's value (Houdini's detail class): the method over every element of the attribute into out_values (its
   components), out_index (may be NULL) the element that gave it (or -1). INVALID_ARGUMENT as above. */
qaws_status qaws_mesh_3d_reduce_attribute(qaws_scalar* out_values, int* out_index, qaws_mesh_3d const* mesh, char const* name,
	qaws_promote_method method);

#endif /* QAWS_MESH_PROMOTE_H */
