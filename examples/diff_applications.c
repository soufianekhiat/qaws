/*
 * Applications of the differentiation API on real data.
 *
 *   1. Hair strands from a photo: orientation field (structure tensor),
 *      evenly spaced streamlines, B-spline strands optimized to follow the
 *      field through unit-tangent adjoints (qaws_diff_geometry.h).
 *   2. Photo vectorization: isocontours become interpolating splines
 *      (centripetal Catmull-Rom and Yuksel C2) whose points are optimized
 *      onto the contours; curvature combs compare the two families.
 *   3. Triangulated mesh to six bicubic patches with ADMM: exact per-patch
 *      least squares at foot points, seam control points in consensus.
 *   4. Non-rigid registration of a template curve to a scan: CMA-ES on the
 *      pose with gradient refinements inside the fitness (basin hopping).
 *   5. Hair grooming on a head: 3D B-spline strands with fixed roots under
 *      length, bending, gravity, collision and guide-alignment energies.
 *   6. Hair from a photo: strands keep their 3D priors while their
 *      projections follow the photo orientation field (photos/hair.ppm,
 *      or the ribbons photo as the flow source).
 *   7. Rational patches: NURBS weights freed after the ADMM fit (blob and
 *      sphere meshes), seams kept closed.
 *   8. Real haircuts: vector hair strands from six portrait photos with a
 *      seeded hair mask (photos/hair_*.ppm, see photos/CREDITS.txt).
 *   9. 3D hair from four frontal portraits: head placed from the face,
 *      visible strands follow the photo, hidden ones the 3D priors.
 *  10. Single-view hair modeling on plain-background portraits: Gabor
 *      orientation, a diffused 3D orientation field, strands grown from
 *      the scalp, rendered by a small rasterizer (photos/plain_*.ppm).
 *
 * Each application lives in diff_applications/NN_*.c; they are compiled as
 * one translation unit, in order (later applications reuse the helpers of
 * earlier ones).
 *
 * QAWS_APP=n in the environment runs only application n.
 *
 * Photos are read as binary PPM; examples/photo_to_ppm.ps1 converts any
 * image (and writes a PNG used as figure background): ribbons.ppm and
 * folds.ppm in the photos directory. Figures are written to showcase/ in
 * the working directory.
 *
 *   qaws_diff_applications [photos directory] [mesh.obj]
 */

#include "example_svg.h"

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

static char const* g_photos = "photos";

#include "diff_applications/00_images.c"
#include "diff_applications/01_hair_strands.c"
#include "diff_applications/02_vectorization.c"
#include "diff_applications/03_mesh_patches.c"
#include "diff_applications/04_registration.c"
#include "diff_applications/05_hair_groom.c"
#include "diff_applications/06_hair_photo.c"
#include "diff_applications/07_rational_patches.c"
#include "diff_applications/08_haircuts.c"
#include "diff_applications/09_hair3d.c"
#include "diff_applications/10_hair_volume.c"

int main(int argc, char** argv)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	if (argc > 1)
		g_photos = argv[1];
	MAKE_DIR("showcase");
	{
		/* QAWS_APP=n runs only application n */
		char const* only = getenv("QAWS_APP");
		int pick = only ? atoi(only) : 0;
		if (!pick || pick == 1) app_hair();
		if (!pick || pick == 2) app_vectorize();
		if (!pick || pick == 3) app_mesh_patches(argc > 2 ? argv[2] : NULL);
		if (!pick || pick == 4) app_registration();
		if (!pick || pick == 5) app_hair_groom();
		if (!pick || pick == 6) app_hair_photo();
		if (!pick || pick == 7) app_rational_patches();
		if (!pick || pick == 8) app_haircuts();
		if (!pick || pick == 9) app_hair3d();
		if (!pick || pick == 10) app_hair_volume();
	}
	return 0;
}
