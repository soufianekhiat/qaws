/*
 * better2d_showcase.c - 2D clipping, offsetting and 2.5D stacks on curves
 *
 * Writes SVG figures into ./showcase/:
 *   b2d1_extract.svg   the exact piece of a curve of every kind between two
 *                      parameters, as a new curve, over its source
 *
 * Each figure lives in better2d_showcase/NN_*.c, compiled as one translation
 * unit.
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

#include "better2d_showcase/01_extract.c"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_extract();
	return 0;
}
