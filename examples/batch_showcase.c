/*
 * batch_showcase.c - Operations over many curves at once
 *
 * Writes SVG figures into ./showcase/:
 *   batch1_heightfield.svg   contour lines and gradient lines of a
 *                            heightfield, every crossing in one batched call,
 *                            and its time against all pairs
 *   batch2_exact_heightfield.svg  the same crossings certified by the exact
 *                            batch, against the pairwise exact call
 *   batch3_terrain.svg       the field as a B-spline terrain: contours as
 *                            intersections with planes, ballistic impacts
 *   batch4_exact_terrain.svg the terrain contours certified by the exact batch
 *   batch5_closest.svg       points snapped to the nearest contour line
 *   batch6_surface_closest.svg points snapped to the nearest terrain point
 *   batch7_raytrace.svg      the terrain ray traced: camera and shadow rays
 *   batch8_threads.svg       every batch on 1 to 16 threads of the caller
 *
 * Each figure lives in batch_showcase/NN_*.c, compiled as one translation
 * unit.
 */

#include "qaws.h"
#include "qaws_exact.h"
#include <math.h>
#include <stdint.h>
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

#include "batch_showcase/01_heightfield.c"
#include "batch_showcase/02_exact_heightfield.c"
#include "batch_showcase/03_terrain.c"
#include "batch_showcase/04_exact_terrain.c"
#include "batch_showcase/05_closest.c"
#include "batch_showcase/06_surface_closest.c"
#include "batch_showcase/07_raytrace.c"
#include "batch_showcase/08_threads.c"

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_heightfield();
	demo_exact_heightfield();
	demo_terrain();
	demo_exact_terrain();
	demo_closest();
	demo_surface_closest();
	demo_raytrace();
	demo_threads();
	return 0;
}
