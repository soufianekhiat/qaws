/*
 * exact_showcase.c - The exact / certified geometry kernel at work
 *
 * Writes SVG figures into ./showcase/:
 *   exact1_orient2d.svg   orientation of points a few ulps apart near a line:
 *                         naive f64 sign map, certified sign map, the cells
 *                         where f64 is wrong, and where the filter needed
 *                         the exact fallback (Kettner et al., "Classroom
 *                         examples of robustness problems")
 *   exact2_bezier.svg     a degree-16 rational Bezier evaluated exactly, the
 *                         reference for f32 and f64 errors (position, C', C'', C''')
 *   exact3_winding.svg    point-in-region near a loop of rational conics:
 *                         sampled winding against the certified one
 *
 * Each figure lives in exact_showcase/NN_*.c, compiled as one translation
 * unit.
 */

#include "qaws.h"
#include "qaws_exact.h"
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

#include "exact_showcase/01_orient2d.c"
#include "exact_showcase/02_exact_bezier.c"
#include "exact_showcase/03_winding.c"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_orient2d();
	demo_exact_bezier();
	demo_winding();
	return 0;
}
