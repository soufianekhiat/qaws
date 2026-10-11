#ifndef QAWS_REGION_SAMPLE_H
#define QAWS_REGION_SAMPLE_H

#include "qaws_types.h"
#include "qaws_status.h"
#include "qaws_path.h"
#include "qaws_clip.h"

/*
 * Points uniform in area over a region bounded by curves of any kind
 * (concave, with holes, under a fill rule), from the boundary alone: no
 * mesh.
 *
 * The boundary is cut into pieces along which x only grows or only falls.
 * Between two consecutive ends of pieces (the critical abscissas), the
 * same pieces cross every vertical line, so the area left of x,
 *
 *     A(x) = sum over the pieces of -sign(dx/dt) * integral y dx,
 *
 * is a smooth function, each term exact by Green's theorem on the true
 * curve. The map of the unit square onto the region
 *
 *     u -> x with A(x) = u * area,
 *     v -> y at the fraction v of the length of the vertical chords at x,
 *
 * keeps area: uniform (u, v) give points uniform over the region, and
 * stratified or low-discrepancy (u, v) stay so. Holes and concavities are
 * more chords at x.
 */

typedef struct qaws_region_sampler qaws_region_sampler;

typedef struct qaws_region_sampler_desc
{
	qaws_path_2d const* paths;      /* closed paths */
	unsigned int path_count;
	qaws_fill_rule fill_rule;       /* the region: the paths under this rule */
	qaws_scalar tolerance;          /* absolute area tolerance; 0 = 1e-13 of the box's area (1e-6 in float) */
} qaws_region_sampler_desc;

/* Sampling patterns. The first five are point sets of the unit square sent
   through the area-keeping map: uniform over the region, their structure
   sheared by the map (the map follows x, then the chords at x). Blue noise
   works on the points themselves: evenly spaced in the region. */
#define QAWS_REGION_SAMPLE_RANDOM     0u   /* independent uniform (u, v), seeded */
#define QAWS_REGION_SAMPLE_SOBOL      1u   /* the 2D Sobol sequence, digitally shifted by the seed */
#define QAWS_REGION_SAMPLE_STRATIFIED 2u   /* correlated multi-jittered (Kensler 2013): one point per row and column of an m x n grid, jittered */
#define QAWS_REGION_SAMPLE_HALTON     3u   /* the Halton sequence in bases 2 and 3, shifted by the seed (Cranley-Patterson) */
#define QAWS_REGION_SAMPLE_R2         4u   /* the R2 Kronecker sequence (plastic number), shifted by the seed */
#define QAWS_REGION_SAMPLE_BLUE_NOISE 5u   /* sample elimination (Yuksel 2015) from 5 n uniform points: Poisson-disk-like spacing */

qaws_status qaws_region_sampler_create(qaws_region_sampler_desc const* desc, qaws_region_sampler** out_sampler);
void qaws_region_sampler_destroy(qaws_region_sampler* sampler);

/* The region's area, and the area of its part left of x (0 left of the
   region, the area right of it). */
qaws_scalar qaws_region_sampler_get_area(qaws_region_sampler const* sampler);
qaws_status qaws_region_sampler_area_left(qaws_region_sampler const* sampler, qaws_scalar x, qaws_scalar* out_area);

/* The area-keeping map of the unit square onto the region; u, v in [0, 1].
   Thread-safe. */
qaws_status qaws_region_sampler_map(qaws_region_sampler const* sampler, qaws_scalar u, qaws_scalar v, qaws_vec2* out_point);
qaws_status qaws_region_sampler_map_many(qaws_region_sampler const* sampler, qaws_vec2 const* uv, unsigned int count,
	qaws_vec2* out_points);

/* count points uniform over the region: QAWS_REGION_SAMPLE_* through the
   map; the same seed gives the same points. */
qaws_status qaws_region_sampler_generate(qaws_region_sampler const* sampler, unsigned int count, unsigned int kind,
	unsigned long long seed, qaws_vec2* out_points);

#endif /* QAWS_REGION_SAMPLE_H */
