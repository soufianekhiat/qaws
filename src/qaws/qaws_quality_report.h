#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_QUALITY_REPORT_H_F32
#define QAWS_QUALITY_REPORT_H_F32
#define QAWS_QUALITY_REPORT_H_BODY
#endif
#else
#ifndef QAWS_QUALITY_REPORT_H_F64
#define QAWS_QUALITY_REPORT_H_F64
#define QAWS_QUALITY_REPORT_H_BODY
#endif
#endif
#ifdef QAWS_QUALITY_REPORT_H_BODY
#undef QAWS_QUALITY_REPORT_H_BODY
#define QAWS_QUALITY_REPORT_H

#include "qaws_mesh.h"
#include "qaws_triangulate.h"
#include "qaws_volume_mesh.h"

/*
 * Quality reports: a measure per element of a 2D, surface or volume mesh,
 * and a summary of those values (statistics, a histogram, the worst
 * elements, those past a solver's threshold).
 *
 * Measures (1 is ideal where it is a ratio; NaN where the measure does not
 * apply to the element, which the summary skips):
 *   MIN_ANGLE, MAX_ANGLE  degrees: a face's corner angles; a cell's dihedral
 *                         angles between facets sharing an edge;
 *   SHAPE                 the mean ratio of simplices: 4 sqrt 3 A / sum l^2
 *                         for a triangle, 12 (3 V)^(2/3) / sum l^2 for a
 *                         tetrahedron (1 regular, 0 degenerate, negative
 *                         inverted); NaN for other elements;
 *   ASPECT_RATIO          simplices: circumradius over the inradius,
 *                         scaled to 1 regular (R / 2r, R / 3r); NaN else;
 *   SCALED_JACOBIAN       the smallest corner scaled Jacobian (the corner's
 *                         edges' determinant over their lengths' product,
 *                         scaled to 1 for the regular simplex, the square
 *                         or the cube; negative where the element folds;
 *                         a surface face's corners measured along its
 *                         Newell normal); a pyramid's base corners;
 *   EDGE_RATIO            longest edge over shortest;
 *   SIZE                  area or volume (signed in 2D and for cells:
 *                         negative when turned).
 *
 * A 2D mesh is measured in its plane: SHAPE, SCALED_JACOBIAN and SIZE are
 * signed by the faces' turn (counter-clockwise positive), so inverted
 * elements show as negative. Polygons of more than four sides have angles,
 * edge ratios and sizes only.
 */

#ifndef QAWS_ENUM_QAWS_ELEMENT_MEASURE
#define QAWS_ENUM_QAWS_ELEMENT_MEASURE
typedef enum qaws_element_measure
{
	QAWS_MEASURE_MIN_ANGLE = 0,
	QAWS_MEASURE_MAX_ANGLE,
	QAWS_MEASURE_SHAPE,
	QAWS_MEASURE_ASPECT_RATIO,
	QAWS_MEASURE_SCALED_JACOBIAN,
	QAWS_MEASURE_EDGE_RATIO,
	QAWS_MEASURE_SIZE
} qaws_element_measure;
#endif

/* out_values: one per face (qaws_mesh_3d_get_faces order) */
qaws_status qaws_mesh_3d_measure_elements(qaws_scalar* out_values, qaws_mesh_3d const* mesh, qaws_element_measure measure);
/* out_values: one per face, triangles then quads */
qaws_status qaws_mesh_2d_measure_elements(qaws_scalar* out_values, qaws_mesh_2d const* mesh, qaws_element_measure measure);
/* out_values: one per cell */
qaws_status qaws_volume_mesh_measure_elements(qaws_scalar* out_values, qaws_volume_mesh const* mesh, qaws_element_measure measure);

typedef struct qaws_quality_summary_desc
{
	qaws_element_measure measure;      /* the values' measure: says which way is worse (MAX_ANGLE, ASPECT_RATIO, EDGE_RATIO higher;
	                                      SIZE and the others lower) */
	int check;                         /* count the elements past threshold */
	qaws_scalar threshold;             /* worse than it: below (or above, as the measure goes) */
	unsigned int bins;                 /* histogram bins (out_histogram's size); 0: no histogram */
	qaws_scalar range[2];              /* the histogram's span; equal: the values' min and max (values outside it go to the end bins) */
	unsigned int worst_count;          /* the worst elements wanted (out_worst's size) */
} qaws_quality_summary_desc;

typedef struct qaws_quality_summary
{
	unsigned int count;                /* values measured (NaN skipped) */
	unsigned int skipped;              /* NaN values */
	qaws_scalar min, max, mean, median;
	qaws_scalar range[2];              /* the histogram's span used */
	unsigned int failing;              /* past the threshold (check) */
	unsigned int worst_found;          /* entries written to out_worst (at most worst_count, at most count) */
} qaws_quality_summary;

/* out_histogram (desc->bins) and out_worst (desc->worst_count, element indices worst first) may be NULL when their counts are 0 */
qaws_status qaws_quality_summarize(qaws_quality_summary* out_summary, unsigned int* out_histogram, unsigned int* out_worst,
	qaws_scalar const* values, unsigned int count, qaws_quality_summary_desc const* desc);

#endif /* QAWS_QUALITY_REPORT_H */
