#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_PARAM_H_F32
#define QAWS_MESH_PARAM_H_F32
#define QAWS_MESH_PARAM_H_BODY
#endif
#else
#ifndef QAWS_MESH_PARAM_H_F64
#define QAWS_MESH_PARAM_H_F64
#define QAWS_MESH_PARAM_H_BODY
#endif
#endif
#ifdef QAWS_MESH_PARAM_H_BODY
#undef QAWS_MESH_PARAM_H_BODY
#define QAWS_MESH_PARAM_H

#include "qaws_mesh.h"

/*
 * Parameterization: a mesh flattened into the plane, its triangles' corners
 * given (u, v) coordinates. Seams (edges given by their two vertices) are
 * cut open; a vertex on a seam gets one (u, v) per side. Each connected
 * piece of the cut mesh (a chart) is flattened on its own and the charts
 * are laid side by side along u. A chart needs a border (the mesh's, or
 * seams): a closed one is INVALID_TOPOLOGY.
 *
 *   LSCM   least squares conformal maps (Levy et al. 2002), as the
 *          conformal energy E_D - A (Mullen et al. 2008) with two border
 *          vertices pinned: angles kept as well as the triangles allow,
 *          free border, one solve;
 *   ARAP   as rigid as possible (Liu et al. 2008): each triangle's map
 *          nearest a rotation, by local-global rounds from LSCM's:
 *          lengths and areas kept better, at the cost of some angle;
 *   BFF    boundary first flattening (Sawhney, Crane 2017): the border
 *          set first, the scale factors u from the curvatures (the
 *          cotangent Laplacian K u = the angle defects inside; on the
 *          border u = 0, its lengths kept, or with bff_disk its turning
 *          set to 2 pi times its share of the border's length, a circle),
 *          the border's new turning k - (K u) and lengths
 *          e^((u_i + u_j) / 2) l closed by the least change of the lengths,
 *          its curve integrated, the inside its harmonic extension. A chart
 *          must be a disk (one border loop: INVALID_TOPOLOGY else).
 *
 * auto_seams cuts the mesh into charts by itself: each triangle grouped by
 * the axis direction nearest its normal (+-x, +-y, +-z), the groups'
 * borders smoothed (a triangle joins the group most of its neighbours are
 * in when its normal faces that group's axis), the edges between groups
 * seams (with the given ones), each connected group a chart: a closed mesh
 * opens without seams given. With max_stretch > 0, a chart whose
 * triangles' area scales (in the plane over on the mesh) differ by more
 * than that factor, or with a triangle turned, is split in two across its
 * longest extent and the whole flattened again (32 splits at most).
 *
 * pack turns each chart to its smallest bounding box (wider than tall),
 * packs them by their shapes (rasterized, largest first at the lowest free
 * spot, a quarter turn allowed, a gap of 2% of the mesh's size between
 * them; a skyline of their boxes if that fails) and scales them all at once
 * into the unit square, the charts keeping their relative sizes; the
 * report's fill is the share of the square they cover.
 *
 * Polygons are fanned from their first vertex (qaws_mesh_3d_get_triangle):
 * the output is per triangle corner, in that order.
 */

#ifndef QAWS_ENUM_QAWS_PARAM_METHOD
#define QAWS_ENUM_QAWS_PARAM_METHOD
typedef enum qaws_param_method
{
	QAWS_PARAM_LSCM = 0,
	QAWS_PARAM_ARAP,
	QAWS_PARAM_BFF
} qaws_param_method;
#endif

typedef struct qaws_param_desc
{
	qaws_param_method method;
	unsigned int iterations;            /* ARAP's local-global rounds; 0: 10 */
	int bff_disk;                       /* BFF: 1 the border a circle; 0 its lengths kept */
	int auto_seams;                     /* 1: cut also between triangles facing different axis directions (6 groups) */
	int pack;                           /* 1: the charts packed into the unit square (shelves, one scale); 0: side by side along u */
	qaws_scalar max_stretch;            /* auto_seams: a chart whose triangles' area scales differ by more (largest over least) split; 0: none */
} qaws_param_desc;

typedef struct qaws_param_report
{
	unsigned int charts;
	unsigned int flipped;               /* triangles turned over in the plane */
	qaws_scalar angle_distortion;       /* the mean of |the corner's angle in the plane - on the mesh| (radians) */
	qaws_scalar area_distortion;        /* the sum over the triangles of |their share of the plane's area - of the mesh's| (0 to 2) */
	qaws_scalar fill;                   /* with pack: the charts' share of the unit square */
} qaws_param_report;

/* out_uv: 2 per triangle corner (3 per triangle); report may be NULL; seam_edges: 2 vertices each (INVALID_ARGUMENT if one is
   not an edge); desc NULL: LSCM */
qaws_status qaws_mesh_3d_parameterize(qaws_scalar* out_uv, qaws_param_report* report, qaws_mesh_3d const* mesh, unsigned int const* seam_edges,
	unsigned int seam_count, qaws_param_desc const* desc);

#endif /* QAWS_MESH_PARAM_H */
