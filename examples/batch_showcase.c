/*
 * batch_showcase.c - Operations over many curves at once
 *
 * Writes SVG figures into ./showcase/:
 *   batch1_heightfield.svg   contour lines and gradient lines of a
 *                            heightfield, every crossing in one batched call,
 *                            and its time against all pairs
 *
 * Each figure lives in batch_showcase/NN_*.c, compiled as one translation
 * unit.
 */

#include "qaws.h"
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

int main(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	demo_heightfield();
	return 0;
}
