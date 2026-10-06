/*
 * diff_showcase.c - Applications of qaws differentiability
 *
 * Writes SVG figures into ./showcase/:
 *   1_curve_fit.svg        B-spline fitted to data by adjoint gradient descent
 *   2_sensitivity.svg      "Affected by" (adjoint) and "Affects" (tangent) maps
 *   3_fairing.svg          curvature fairing through the geometry kernel
 *   4_nurbs_weight.svg     learning a NURBS weight: a parabola becomes a circle
 *   5_surface.svg          height-field fit of a B-spline surface (z only, by
 *                          component mask) and Gaussian curvature sensitivity
 *   6_vase.svg             surface of revolution fitted to scan points through
 *                          its profile curve (child views)
 *   7_coons.svg            Coons patch faired by editing two boundary curves
 *   8_networks.svg ... 13_knot_placement.svg: curve networks, projection
 *                          fitting, soap film, parameter correction, arch
 *                          features, knot placement
 *   14_arc_length.svg      constant-speed samples: first and second order
 *                          tangents, and a fit by Newton-CG on exact HVPs
 *   15_cdf_measures.svg    inverse-CDF samples under arc length, curvature and
 *                          a density field, with their sample tangents
 *   16_surface_cdf.svg     stratified points warped onto a patch by the area,
 *                          density and curvature inverse CDFs
 *   17_blue_noise.svg      blue noise on a patch: repulsion of the warped samples
 *                          descended through the xi adjoints of the warp
 *   18_knot_sampling.svg   knot tangents of arc-length samples, and interior
 *                          knots moved by the knot adjoint toward arc-length
 *                          parameters
 *
 * Each figure lives in diff_showcase/NN_*.c; they are compiled as one
 * translation unit, in order (later figures reuse earlier helpers).
 *
 * Only the public qaws API is used. No finite differences anywhere: every
 * gradient comes from the library's adjoint rules.
 */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

#include "example_svg.h"

#include "diff_showcase/01_curve_fit.c"
#include "diff_showcase/02_sensitivity.c"
#include "diff_showcase/03_fairing.c"
#include "diff_showcase/04_nurbs_weight.c"
#include "diff_showcase/05_surface.c"
#include "diff_showcase/06_vase.c"
#include "diff_showcase/07_coons.c"
#include "diff_showcase/08_networks.c"
#include "diff_showcase/09_projection.c"
#include "diff_showcase/10_soap_film.c"
#include "diff_showcase/11_parameter_correction.c"
#include "diff_showcase/12_arch.c"
#include "diff_showcase/13_knot_placement.c"
#include "diff_showcase/14_arc_length.c"
#include "diff_showcase/15_cdf_measures.c"
#include "diff_showcase/16_surface_cdf.c"
#include "diff_showcase/17_blue_noise.c"
#include "diff_showcase/18_knot_sampling.c"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_curve_fit();
	demo_sensitivity();
	demo_fairing();
	demo_nurbs_weight();
	demo_surface();
	demo_vase();
	demo_coons();
	demo_networks();
	demo_projection();
	demo_soap_film();
	demo_parameter_correction();
	demo_arch();
	demo_knot_placement();
	demo_arc_length();
	demo_cdf_measures();
	demo_surface_cdf();
	demo_blue_noise();
	demo_knot_sampling();
	return 0;
}
