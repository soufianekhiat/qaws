/*
 * exact_showcase.c - The exact / certified geometry kernel at work
 *
 * Writes SVG figures into ./showcase/:
 *   exact1_orient2d.svg   orientation of points a few ulps apart near a line:
 *                         naive f64 sign map, certified sign map, the cells
 *                         where f64 is wrong, and where the filter needed
 *                         the exact fallback (Kettner et al., "Classroom
 *                         examples of robustness problems")
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

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_orient2d();
	return 0;
}
