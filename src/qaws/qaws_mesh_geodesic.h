#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_GEODESIC_H_F32
#define QAWS_MESH_GEODESIC_H_F32
#define QAWS_MESH_GEODESIC_H_BODY
#endif
#else
#ifndef QAWS_MESH_GEODESIC_H_F64
#define QAWS_MESH_GEODESIC_H_F64
#define QAWS_MESH_GEODESIC_H_BODY
#endif
#endif
#ifdef QAWS_MESH_GEODESIC_H_BODY
#undef QAWS_MESH_GEODESIC_H_BODY
#define QAWS_MESH_GEODESIC_H

#include "qaws_mesh.h"

/*
 * Geodesic distances and paths on a mesh, by the heat method (Crane,
 * Weischedel, Wardetzky 2013): heat diffused from the sources for a short
 * time, its gradient normalized, and the distance the function whose
 * gradient best matches it (a Poisson equation). The two matrices (the
 * heat step's M + t L, the cotangent Laplacian L) are factored once at
 * creation by the sparse LDL^T; each distance query is then a few solves,
 * so many sources or many queries cost little.
 *
 * The mesh is read as its triangles (polygons fanned from their first
 * vertex, as qaws_mesh_3d_get_triangle): locations name those triangles.
 * The time is t = factor h^2, h the mean edge length (factor 1, the
 * paper's choice, by default: larger smooths the distance, smaller
 * follows the mesh more closely but grows noisier). On a mesh with a
 * border, the heat step averages its Neumann and Dirichlet solutions
 * (the paper's advice: the first bends isolines toward the border, the
 * second away from it).
 *
 * Distances are per vertex, from the nearest source, each connected part
 * shifted so its nearest source is at 0; vertices no source reaches
 * (another part) get HUGE_VAL. The accuracy is that of the triangles: about
 * the mean edge length's order in the worst places, much better on
 * fine, well-shaped meshes.
 *
 * Paths go down the distance's gradient (constant in each triangle) from a
 * location to the nearest source, across triangles, along an edge where
 * both sides lead onto it and through a vertex where the flow meets one,
 * ending straight at the source once in its triangle.
 *
 * The object copies what it needs: the mesh may change or go after
 * creation. Thread-safe for concurrent queries.
 */

#ifndef QAWS_ENUM_QAWS_GEODESIC_BOUNDARY
#define QAWS_ENUM_QAWS_GEODESIC_BOUNDARY
typedef enum qaws_geodesic_boundary
{
	QAWS_GEODESIC_BOUNDARY_MEAN = 0,    /* the mean of the two below */
	QAWS_GEODESIC_BOUNDARY_NEUMANN,     /* no heat through the border */
	QAWS_GEODESIC_BOUNDARY_DIRICHLET    /* no heat on the border */
} qaws_geodesic_boundary;
#endif

typedef struct qaws_mesh_geodesic_desc
{
	qaws_scalar time;                   /* the heat's time over the mean edge length squared; 0: 1 */
	qaws_geodesic_boundary boundary;    /* on meshes with a border */
} qaws_mesh_geodesic_desc;

typedef struct qaws_mesh_geodesic qaws_mesh_geodesic;

/* desc NULL: the defaults. INVALID_ARGUMENT for a mesh without a triangle of area */
qaws_status qaws_mesh_geodesic_create(qaws_mesh_geodesic** out_geodesic, qaws_mesh_3d const* mesh, qaws_mesh_geodesic_desc const* desc);
void qaws_mesh_geodesic_destroy(qaws_mesh_geodesic* geodesic);

/* the distance at every vertex (vertex count) to the nearest of the sources */
qaws_status qaws_mesh_geodesic_distance(qaws_scalar* out_distances, qaws_mesh_geodesic const* geodesic, qaws_mesh_location const* sources,
	unsigned int source_count);

/* the path from a location down the distances (from qaws_mesh_geodesic_distance or qaws_mesh_3d_geodesic_exact with these
   sources) to the nearest source: its points (3 each) and locations, from the start to the source; either output may be NULL.
   BUFFER_TOO_SMALL with the count when capacity is short */
qaws_status qaws_mesh_geodesic_path(qaws_scalar* out_points, qaws_mesh_location* out_locations, unsigned int* out_count,
	qaws_mesh_geodesic const* geodesic, qaws_scalar const* distances, qaws_mesh_location const* sources, unsigned int source_count,
	qaws_mesh_location from, unsigned int capacity);

/* Exact polyhedral distances at the vertices (the triangles taken as they
   are, polygons fanned): Chen and Han's windows with Xin and Wang's
   filtering (ICH, 2009). A window is an interval of an edge seen straight
   from a source unfolded into the plane of the triangle it enters; windows
   go cheapest first, each triangle's far edges taking the part its rays
   cross; a saddle vertex (its angles summing to more than 2 pi) or a border
   vertex reached becomes a pseudo-source, expanded once at its distance; a
   window splitting at a vertex no sooner than an earlier one through the
   same edge keeps one child (one angle, one split). Without a factor or a
   solve: no object needed. HUGE_VAL where no source reaches. The cost grows
   with the windows (about n sqrt n on smooth meshes, n^2 at worst).
   max_distance stops the windows once the cheapest left lies beyond it: the
   distances up to it are the same as the full run's, the farther ones
   HUGE_VAL, the cost the windows within it (HUGE_VAL: every distance).
   INVALID_ARGUMENT for a negative or NaN max_distance. */
qaws_status qaws_mesh_3d_geodesic_exact(qaws_scalar* out_distances, qaws_mesh_3d const* mesh, qaws_mesh_location const* sources,
	unsigned int source_count, qaws_scalar max_distance);

/* The exact distances at points anywhere on the surface (point_count
   locations), not only at the vertices: the windows of
   qaws_mesh_3d_geodesic_exact kept, each point's distance the least of a
   window entering its triangle that sees it through its interval, a vertex
   of its triangle then straight, a source in its triangle straight.
   INVALID_ARGUMENT for a point outside the triangles. */
qaws_status qaws_mesh_3d_geodesic_exact_at(qaws_scalar* out_distances, qaws_mesh_3d const* mesh, qaws_mesh_location const* sources,
	qaws_mesh_location const* points, unsigned int source_count, unsigned int point_count);

/* The nearest source of every vertex (an index into sources; QAWS_MESH_NONE
   where the distance is not finite), from the distances to the nearest of
   them (the heat method's or the exact ones): the geodesic Voronoi regions.
   Vertices go nearest first, each taking the source of what explains its
   distance best: a source in its triangle straight, or a neighbour already
   labelled plus their edge (Dijkstra's parent). Vertices farther than an
   edge from the regions' borders are labelled as the distances say; along
   the borders a vertex may take either side. */
qaws_status qaws_mesh_3d_geodesic_nearest(unsigned int* out_sources, qaws_mesh_3d const* mesh, qaws_scalar const* distances,
	qaws_mesh_location const* sources, unsigned int source_count);

#endif /* QAWS_MESH_GEODESIC_H */
