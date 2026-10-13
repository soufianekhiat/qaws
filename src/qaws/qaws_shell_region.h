#ifndef QAWS_SHELL_REGION_H
#define QAWS_SHELL_REGION_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_surface_types.h"
#include "qaws_clip.h"
#include "qaws_surface_region.h"
#include "qaws_brep.h"

/*
 * Curves and regions across the faces of a shell (qaws_surface_region.h
 * carried over shared edges).
 *
 * The atlas is the faces whose parameter domains meet along their sides: a
 * side of one face whose image coincides in space with a side of another
 * (in either direction) is a shared edge. A point crossing it moves to the
 * other face's side point at the same place (by inversion along that side).
 * The faces must be oriented alike (S_u x S_v all outward, or all inward):
 * a shared edge runs opposite ways in the two faces' counter-clockwise
 * boundaries. Sides along a face's closed direction are its seam, and a
 * side that is one point is a pole; neither is an edge.
 *
 * Geodesics are continued over edges as the straightest curves: the
 * direction keeps its components along the edge and across it (the
 * across part turned into the next face's tangent plane). Across a smooth
 * join this is the geodesic itself; across a crease it unfolds the two
 * faces about the edge. A geodesic reaching a face's corner (a shell
 * vertex) is not continued.
 *
 * A region on the shell is its part on every face, each a qaws_surface_region
 * of that face: the boundary's pieces on the face closed along the face's
 * sides (counter-clockwise, from where a piece leaves to where the next
 * comes in), faces without pieces inside or outside as their neighbours
 * say. Areas, perimeters and Booleans go face by face.
 */

typedef struct qaws_shell_atlas qaws_shell_atlas;
typedef struct qaws_shell_curve qaws_shell_curve;
typedef struct qaws_shell_region qaws_shell_region;

/* a point of the shell: a face and its (u, v) */
typedef struct qaws_shell_point
{
	unsigned int face;
	qaws_vec2 uv;
} qaws_shell_point;

/* The atlas of the surfaces (borrowed: they must outlive it). tolerance:
   the distance in space at which sides coincide; 0 = 1e-7 of the shell's
   size (1e-4 in float). INVALID_ARGUMENT when faces meet oriented
   oppositely. */
qaws_status qaws_shell_atlas_create(qaws_surface const* const* surfaces, unsigned int surface_count, qaws_scalar tolerance,
	qaws_shell_atlas** out_atlas);
/* The atlas of a B-rep shell's faces. */
qaws_status qaws_shell_atlas_create_from_brep(qaws_brep_shell const* shell, qaws_scalar tolerance, qaws_shell_atlas** out_atlas);
void qaws_shell_atlas_destroy(qaws_shell_atlas* atlas);

unsigned int qaws_shell_atlas_get_face_count(qaws_shell_atlas const* atlas);
qaws_surface const* qaws_shell_atlas_get_face(qaws_shell_atlas const* atlas, unsigned int face);
/* The shared edges found (each once). */
unsigned int qaws_shell_atlas_get_edge_count(qaws_shell_atlas const* atlas);

/* The point and velocity reached from p along the straightest geodesic
   leaving it in the (u, v) direction of its face, for the surface length
   `length`, across edges; the faces crossed. OUT_OF_RANGE where it leaves
   the shell (a free side) or reaches a vertex. */
qaws_status qaws_shell_exp(qaws_shell_atlas const* atlas, qaws_shell_point p, qaws_vec2 direction, qaws_scalar length,
	qaws_shell_point* out_point, qaws_vec3* out_velocity, unsigned int* out_crossings);

/* That geodesic as a shell curve: one surface curve per face run. */
qaws_status qaws_shell_curve_create_geodesic_ray(qaws_shell_atlas const* atlas, qaws_shell_point p, qaws_vec2 direction,
	qaws_scalar length, qaws_scalar tolerance, qaws_shell_curve** out_curve);

/* The geodesic circle of radius r about the centre across faces, counter-
   clockwise: its runs on each face fitted as (u, v) Hermite curves through
   points of the exponential map, the edge crossings found by bisection on
   the angle. */
qaws_status qaws_shell_curve_create_geodesic_circle(qaws_shell_atlas const* atlas, qaws_shell_point centre, qaws_scalar radius,
	qaws_scalar tolerance, qaws_shell_curve** out_curve);

void qaws_shell_curve_destroy(qaws_shell_curve* curve);

/* The runs, in order along the curve: the face and its surface curve
   (owned by the shell curve). */
unsigned int qaws_shell_curve_get_run_count(qaws_shell_curve const* curve);
qaws_status qaws_shell_curve_get_run(qaws_shell_curve const* curve, unsigned int index, unsigned int* out_face,
	qaws_surface_curve const** out_curve);
int qaws_shell_curve_is_closed(qaws_shell_curve const* curve);
qaws_status qaws_shell_curve_get_length(qaws_shell_curve const* curve, qaws_scalar* out_length);

/* The region bounded by closed shell curves (each counter-clockwise round
   what it encloses; curves inside one face may be any of that face's
   surface curves through qaws_shell_region_create_on_face). */
qaws_status qaws_shell_region_create(qaws_shell_atlas const* atlas, qaws_shell_curve const* const* loops, unsigned int loop_count,
	qaws_shell_region** out_region);

/* The Boolean of two regions of one atlas, face by face. */
qaws_status qaws_shell_region_clip(qaws_shell_region const* a, qaws_shell_region const* b, qaws_clip_type clip_type,
	qaws_shell_region** out_region);

void qaws_shell_region_destroy(qaws_shell_region* region);

/* The region's part on a face: NULL when it has none. */
qaws_surface_region const* qaws_shell_region_get_face(qaws_shell_region const* region, unsigned int face);
qaws_status qaws_shell_region_get_area(qaws_shell_region const* region, qaws_scalar* out_area);
qaws_status qaws_shell_region_get_perimeter(qaws_shell_region const* region, qaws_scalar* out_length);
/* 1 inside, 0 outside, -1 on the boundary. */
qaws_status qaws_shell_region_contains(qaws_shell_region const* region, qaws_shell_point p, int* out_inside);

#endif /* QAWS_SHELL_REGION_H */
