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

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_extract();
	demo_boolean();
	demo_clipper2();
	demo_offset();
	demo_stack();
	return 0;
}
