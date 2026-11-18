#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_SHATTER_H_F32
#define QAWS_SHATTER_H_F32
#define QAWS_SHATTER_H_BODY
#endif
#else
#ifndef QAWS_SHATTER_H_F64
#define QAWS_SHATTER_H_F64
#define QAWS_SHATTER_H_BODY
#endif
#endif
#ifdef QAWS_SHATTER_H_BODY
#undef QAWS_SHATTER_H_BODY
#define QAWS_SHATTER_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"
#include "qaws_path.h"

/*
 * Seams and shatter of curved surfaces, as Houdini's Boolean "Seam" and
 * "Shatter".
 *
 * Two operands, A and B, are lists of surfaces (trimmed or not): the faces
 * of a body or a sheet. Each surface is laid on a grid of resolution^2
 * cells over its (u, v) domain, two triangles a cell; the triangles of
 * surfaces that may meet are paired by a grid of their boxes, each pair's
 * crossing found, and the crossings chained into curves by the triangle
 * edges they cross. Every point is then moved onto both true surfaces by
 * Newton's method on S1(u1, v1) = S2(u2, v2) (the least change of the four
 * parameters, as long as the gap shrinks), and a point of a trimmed surface
 * outside its trim is left out (the curve broken there).
 *
 * Seams: the curves of the kinds asked, QAWS_SEAM_A_B where a surface of A
 * crosses one of B, QAWS_SEAM_A_A (QAWS_SEAM_B_B) where two surfaces of A
 * (of B) cross - not along a border they share. A curve runs between ends
 * (a surface's border, three surfaces meeting) or is closed (its first point
 * not repeated). Each point carries its position and its parameters on the
 * two surfaces (surface1 the lower operand, then the lower index).
 *
 * Shatter: every surface of A cut along its A-B seams into pieces, each a
 * trimmed surface over the A surface's base (its trim loops, and the seams,
 * flattened into polygons in (u, v)). With b_is_solid (B's surfaces bound a
 * solid, their normals outward), each piece is told inside or outside B by
 * the winding number of B's triangles at a point of it (inside from 1/2).
 * A piece is valid while the shatter lives and A's surfaces do.
 *
 * Tangential contacts and crossings finer than the grid are not resolved:
 * raise the resolution. Thread-safe on immutable surfaces.
 */

#define QAWS_SEAM_A_B 1u
#define QAWS_SEAM_A_A 2u
#define QAWS_SEAM_B_B 4u

typedef struct qaws_surface_seam_desc
{
	qaws_surface const* const* a;     /* operand A's surfaces */
	unsigned int a_count;
	qaws_surface const* const* b;     /* operand B's surfaces (may be none for A_A seams) */
	unsigned int b_count;
	unsigned int kinds;               /* QAWS_SEAM_*: the seams wanted (shatter: A_B always) */
	unsigned int resolution;          /* grid cells along each side of a surface's domain; 0: 64 */
	qaws_scalar tolerance;            /* a point's gap between the two surfaces; 0: 1e-9 (f32: 1e-5) of the extent */
} qaws_surface_seam_desc;

typedef struct qaws_seam_point
{
	qaws_vec3 position;
	qaws_scalar u1, v1;               /* on surface1 */
	qaws_scalar u2, v2;               /* on surface2 */
} qaws_seam_point;

typedef struct qaws_seam_info
{
	unsigned int kind;                /* QAWS_SEAM_* */
	unsigned int operand1, surface1;  /* operand 0: A, 1: B; the surface's index in it */
	unsigned int operand2, surface2;
	unsigned int point_count;
	int closed;
} qaws_seam_info;

typedef struct qaws_surface_seams qaws_surface_seams;

/* INVALID_ARGUMENT: no surface, a NULL one, kinds 0 */
qaws_status qaws_surface_seams_create(qaws_surface_seams** out_seams, qaws_surface_seam_desc const* desc);
void qaws_surface_seams_destroy(qaws_surface_seams* seams);
unsigned int qaws_surface_seams_get_count(qaws_surface_seams const* seams);
/* seam index's description, and its points (info.point_count of them; out_points may be NULL) */
qaws_status qaws_surface_seams_get(qaws_seam_info* out_info, qaws_seam_point* out_points, qaws_surface_seams const* seams, unsigned int index);

typedef struct qaws_surface_shatter qaws_surface_shatter;

/* INVALID_ARGUMENT as qaws_surface_seams_create, or no surface in B */
qaws_status qaws_surface_shatter_create(qaws_surface_shatter** out_shatter, qaws_surface_seam_desc const* desc, int b_is_solid);
void qaws_surface_shatter_destroy(qaws_surface_shatter* shatter);
unsigned int qaws_surface_shatter_get_count(qaws_surface_shatter const* shatter);
/* piece index: the trimmed surface (owned by the shatter), the A surface it was cut from, and 1 inside B, 0 outside, -1 not
   told (B not a solid); any output may be NULL */
qaws_status qaws_surface_shatter_get(qaws_surface const** out_piece, unsigned int* out_surface, int* out_inside,
	qaws_surface_shatter const* shatter, unsigned int index);

/*
 * Seams and shatter of 2D regions of curves.
 *
 * Operands are lists of paths (qaws_path.h). A is a region (closed paths
 * read with a_rule); B is a region (b_rule) or, with b_is_region 0, curves
 * with no inside (open or closed paths), as a surface is in 3D. Every curve
 * is flattened into chords within flatness, the chords of both arranged,
 * and every crossing moved onto both curves by Newton's method.
 *
 * Seams: the points where a curve of A crosses one of B (QAWS_SEAM_A_B), or
 * two curves of one operand cross (QAWS_SEAM_A_A, QAWS_SEAM_B_B; not where a
 * path's curves join), each with its parameter on both curves (curve1 the
 * lower operand, path and curve).
 *
 * Shatter: A's region cut along B's curves into pieces - A's own curves
 * where its region lies on both sides (crossings of its paths inside it) do
 * not cut it - each a region of closed paths (its outer path
 * counter-clockwise, its holes clockwise) whose
 * curves are exact pieces of A's and B's (qaws_curve_extract, reversed where
 * walked backwards), each told inside or outside B when B is a region. The
 * paths and their curves are owned by the shatter. Inside and outside are
 * told by the winding numbers carried from face to face across the curves
 * (each operand's fill rule on its own), not by sampling.
 *
 * Tangencies and crossings closer than the flatness are not resolved.
 */
typedef struct qaws_region_seam_desc
{
	qaws_path_2d const* a;
	unsigned int a_count;
	qaws_fill_rule a_rule;
	qaws_path_2d const* b;
	unsigned int b_count;
	qaws_fill_rule b_rule;
	int b_is_region;                  /* 0: B's paths are curves with no inside */
	unsigned int kinds;               /* QAWS_SEAM_*: the seams wanted */
	qaws_scalar flatness;             /* the chords' deviation from the curves; 0: 1e-3 of the extent */
	qaws_scalar tolerance;            /* a crossing's gap between its curves; 0: 1e-12 (f32: 1e-6) of the extent */
} qaws_region_seam_desc;

typedef struct qaws_region_seam_point
{
	qaws_vec2 position;
	unsigned int kind;                /* QAWS_SEAM_* */
	unsigned int operand1, path1, curve1;
	qaws_scalar t1;
	unsigned int operand2, path2, curve2;
	qaws_scalar t2;
} qaws_region_seam_point;

/* The crossing points (at most capacity written; out_points may be NULL to count), their number in *out_count;
   BUFFER_TOO_SMALL when there are more than capacity. INVALID_ARGUMENT: no path in A (or in B), a path without curves,
   kinds 0. */
qaws_status qaws_region_seams_2d(qaws_region_seam_point* out_points, unsigned int* out_count, qaws_region_seam_desc const* desc,
	unsigned int capacity);

typedef struct qaws_region_shatter qaws_region_shatter;

/* INVALID_ARGUMENT as qaws_region_seams_2d (kinds ignored) */
qaws_status qaws_region_shatter_create_2d(qaws_region_shatter** out_shatter, qaws_region_seam_desc const* desc);
void qaws_region_shatter_destroy(qaws_region_shatter* shatter);
unsigned int qaws_region_shatter_get_count(qaws_region_shatter const* shatter);
/* piece index: its paths (outer first, then holes) and their count, and 1 inside B, 0 outside, -1 not told; any output
   may be NULL */
qaws_status qaws_region_shatter_get(qaws_path_2d const** out_paths, unsigned int* out_path_count, int* out_inside,
	qaws_region_shatter const* shatter, unsigned int index);

#endif /* QAWS_SHATTER_H */
