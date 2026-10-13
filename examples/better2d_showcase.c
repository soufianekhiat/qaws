/*
 * better2d_showcase.c - 2D clipping, offsetting and 2.5D stacks on curves
 *
 * Writes SVG figures into ./showcase/:
 *   b2d1_extract.svg   the exact piece of a curve of every kind between two
 *                      parameters, as a new curve, over its source
 *   b2d2_boolean.svg   Boolean operations on polygons, NURBS and B-spline
 *                      regions, fill rules, nesting, open paths
 *   b2d3_clipper2.svg  Clipper2's own test records through qaws_clip
 *   b2d4_offset.svg    offsets with every join and end type, curves and holes
 *   b2d5_stack.svg     2.5D stacks: interpolated sections, splits, Booleans by level
 *   b2d11_region_sample.svg  points uniform or by a density over Boolean
 *                      results of lines, arcs and Beziers (qaws_region_sampler)
 *   b2d12_skeleton.svg straight skeletons: roofs, wavefronts, degenerate events
 *   b2d13_curved_skeleton.svg  skeletons of lines and arcs, mitered and medial axis
 *   b2d14_freeform_skeleton.svg  skeletons of Bezier and NURBS regions through arc splines
 *   b2d15_mean_value.svg  mean value coordinates over curved regions: a cage
 *                      deformation by the curves' own control points, interpolants
 *   b2d16_surface_geodesic.svg  geodesic polygons and circles traced on a sphere,
 *                      a bicubic patch and a torus, beside their parameter planes
 *   b2d17_surface_boolean.svg  Booleans of geodesic squares and circles on a patch
 *                      and a sphere, with their surface areas
 *   b2d18_surface_seams.svg  regions across seams and round poles on a sphere,
 *                      a cylinder and a torus, with their parameter planes
 *   b2d19_surface_offset.svg  geodesic offsets: grown, shrunk and tubes on a sphere,
 *                      a wavy patch and a torus, across seams
 *   b2d20_surface_sample.svg  sampling and meshing in the metric: Sobol points, blue
 *                      noise, a density, meshes on a sphere, a wavy patch and a torus
 *   b2d21_surface_exact.svg  exact Booleans on surfaces: shared and touching edges,
 *                      a ring of squares across the seam united into a band
 *   b2d22_surface_medial.svg  medial axes on surfaces: junctions and the largest
 *                      inscribed geodesic circles on a sphere, a patch and a torus
 *   b2d23_shell_region.svg  curves and regions across shell faces: geodesics over
 *                      a cube's edges, circles across faces, two hemispheres
 *   b2d24_shell_saddle.svg  geodesic disks through saddle vertices: five and six
 *                      faces round a vertex, flat and bent; Booleans on a bent one
 *   b2d25_shell_trimmed.svg  shells of trimmed faces: a plane split by a circle, a
 *                      plate with a hole and a tube, a cube of trimmed planes, a T-junction,
 *                      a cylinder and a torus with holes on their seams
 *
 * Each figure lives in better2d_showcase/NN_*.c, compiled as one translation
 * unit.
 */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

#include "example_svg.h"

#include "better2d_showcase/01_extract.c"
#include "better2d_showcase/02_boolean.c"
#include "better2d_showcase/03_clipper2.c"
#include "better2d_showcase/04_offset.c"
#include "better2d_showcase/05_stack.c"
#include "better2d_showcase/06_curve_offset.c"
#include "better2d_showcase/07_offset_extremes.c"
#include "better2d_showcase/08_triangulate.c"
#include "better2d_showcase/09_stack_mesh.c"
#include "better2d_showcase/10_clip_exact.c"
#include "better2d_showcase/11_region_sample.c"
#include "better2d_showcase/12_skeleton.c"
#include "better2d_showcase/13_curved_skeleton.c"
#include "better2d_showcase/14_freeform_skeleton.c"
#include "better2d_showcase/15_mean_value.c"
#include "better2d_showcase/16_surface_geodesic.c"
#include "better2d_showcase/17_surface_boolean.c"
#include "better2d_showcase/18_surface_seams.c"
#include "better2d_showcase/19_surface_offset.c"
#include "better2d_showcase/20_surface_sample.c"
#include "better2d_showcase/21_surface_exact.c"
#include "better2d_showcase/22_surface_medial.c"
#include "better2d_showcase/23_shell_region.c"
#include "better2d_showcase/24_shell_saddle.c"
#include "better2d_showcase/25_shell_trimmed.c"
#include "better2d_showcase/26_shell_paths.c"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_extract();
	demo_boolean();
	demo_clipper2();
	demo_offset();
	demo_stack();
	demo_curve_offset();
	demo_offset_extremes();
	demo_triangulate();
	demo_stack_mesh();
	demo_clip_exact();
	demo_region_sample();
	demo_skeleton();
	demo_curved_skeleton();
	demo_freeform_skeleton();
	demo_mean_value();
	demo_surface_geodesic();
	demo_surface_boolean();
	demo_surface_seams();
	demo_surface_offset();
	demo_surface_sample();
	demo_surface_exact();
	demo_surface_medial();
	demo_shell_region();
	demo_shell_saddle();
	demo_shell_trimmed();
	demo_shell_paths();
	return 0;
}
