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
#include "better2d_showcase/27_geodesic_enclosure.c"
#include "better2d_showcase/28_geodesic_query.c"
#include "better2d_showcase/29_shell_mesh.c"
#include "better2d_showcase/30_mitered_certified.c"
#include "better2d_showcase/31_repulsive_energy.c"
#include "better2d_showcase/32_repulsive_solve.c"
#include "better2d_showcase/33_repulsive_mesh.c"
#include "better2d_showcase/34_repulsive_bh.c"
#include "better2d_showcase/35_repulsive_hsolve.c"
#include "better2d_showcase/36_repulsive_surface.c"
#include "better2d_showcase/37_repulsive_packing.c"
#include "better2d_showcase/38_weighted_skeleton.c"

/* figure n wanted: every one with no arguments, else those listed */
static int b2d_want(int argc, char** argv, int n)
{
	int i;
	for (i = 1; i < argc; i++)
		if (atoi(argv[i]) == n)
			return 1;
	return argc < 2;
}

int main(int argc, char** argv)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	if (b2d_want(argc, argv, 1)) demo_extract();
	if (b2d_want(argc, argv, 2)) demo_boolean();
	if (b2d_want(argc, argv, 3)) demo_clipper2();
	if (b2d_want(argc, argv, 4)) demo_offset();
	if (b2d_want(argc, argv, 5)) demo_stack();
	if (b2d_want(argc, argv, 6)) demo_curve_offset();
	if (b2d_want(argc, argv, 7)) demo_offset_extremes();
	if (b2d_want(argc, argv, 8)) demo_triangulate();
	if (b2d_want(argc, argv, 9)) demo_stack_mesh();
	if (b2d_want(argc, argv, 10)) demo_clip_exact();
	if (b2d_want(argc, argv, 11)) demo_region_sample();
	if (b2d_want(argc, argv, 12)) demo_skeleton();
	if (b2d_want(argc, argv, 13)) demo_curved_skeleton();
	if (b2d_want(argc, argv, 14)) demo_freeform_skeleton();
	if (b2d_want(argc, argv, 15)) demo_mean_value();
	if (b2d_want(argc, argv, 16)) demo_surface_geodesic();
	if (b2d_want(argc, argv, 17)) demo_surface_boolean();
	if (b2d_want(argc, argv, 18)) demo_surface_seams();
	if (b2d_want(argc, argv, 19)) demo_surface_offset();
	if (b2d_want(argc, argv, 20)) demo_surface_sample();
	if (b2d_want(argc, argv, 21)) demo_surface_exact();
	if (b2d_want(argc, argv, 22)) demo_surface_medial();
	if (b2d_want(argc, argv, 23)) demo_shell_region();
	if (b2d_want(argc, argv, 24)) demo_shell_saddle();
	if (b2d_want(argc, argv, 25)) demo_shell_trimmed();
	if (b2d_want(argc, argv, 26)) demo_shell_paths();
	if (b2d_want(argc, argv, 27)) demo_geodesic_enclosure();
	if (b2d_want(argc, argv, 28)) demo_geodesic_query();
	if (b2d_want(argc, argv, 29)) demo_shell_mesh();
	if (b2d_want(argc, argv, 30)) demo_mitered_certified();
	if (b2d_want(argc, argv, 31)) demo_repulsive_energy();
	if (b2d_want(argc, argv, 32)) demo_repulsive_solve();
	if (b2d_want(argc, argv, 33)) demo_repulsive_mesh();
	if (b2d_want(argc, argv, 34)) demo_repulsive_bh();
	if (b2d_want(argc, argv, 35)) demo_repulsive_hsolve();
	if (b2d_want(argc, argv, 36)) demo_repulsive_surface();
	if (b2d_want(argc, argv, 37)) demo_repulsive_packing();
	if (b2d_want(argc, argv, 38)) demo_weighted_skeleton();
	return 0;
}
