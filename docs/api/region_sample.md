# API Reference: qaws_region_sample.h

Points uniform in area over a region bounded by curves of any kind:
concave, with holes, under a fill rule. Everything comes from the boundary;
there is no mesh.

```c
qaws_region_sampler_desc d = { paths, path_count, QAWS_FILL_NON_ZERO, 0 };
qaws_region_sampler* s;
qaws_region_sampler_create(&d, &s);

qaws_region_sampler_get_area(s);
qaws_region_sampler_area_left(s, x, &a);            /* A(x): the area left of x */
qaws_region_sampler_map(s, u, v, &p);               /* the unit square onto the region */
qaws_region_sampler_map_many(s, uv, n, points);
qaws_region_sampler_generate(s, n, QAWS_REGION_SAMPLE_SOBOL, seed, points);
qaws_region_sampler_destroy(s);
```

## How it works

1. **Region.** The paths are united under the fill rule (`qaws_clip`), so
   the boundary is simple: outer loops counter-clockwise, holes clockwise.
   A difference is a union with the subtracted paths reversed and
   `QAWS_FILL_POSITIVE`.
2. **Pieces.** Every curve is cut at its span ends and wherever dx/dt
   changes sign. Along each piece x only grows (sigma = +1) or only falls
   (sigma = -1). Vertical pieces bound no area to their side and are dropped.
3. **Slabs.** The pieces' end abscissas, sorted, cut the plane into slabs.
   Inside a slab the same pieces cross every vertical line, so
   ```
   A(x) = A(x_i) - sum over the slab's pieces of sigma * (integral of y dx from x_i to x)
   ```
   Each integral is exact on the true curve, by Green's theorem:
   ∫ y dx = [x y] − ∫ x dy, the last term computed span by span to the
   tolerance. A at every slab's left side is summed once, at creation.
4. **The map.**
   - `u` locates the slab holding `u * area`. Inside it, `x` comes from
     Newton on A(x) − u * area, kept in a shrinking bracket. The derivative
     of A is the chord length at `x`, so the iteration converges
     quadratically.
   - `v` then picks `y` along the vertical chords at `x`. The pieces' y
     values are sorted, and a gap is inside when the winding below it is
     positive (a piece running to +x has the region above it). `y` sits at
     the fraction `v` of the chords' total length.

   The map keeps area, so uniform (u, v) give uniform points. Stratified or
   low-discrepancy (u, v) stay stratified: the Sobol variant gives moment
   errors around 1e-4 with 65,536 points.

## Generators

| `kind` | (u, v) |
|---|---|
| `QAWS_REGION_SAMPLE_RANDOM` | independent uniform pairs from splitmix64, seeded |
| `QAWS_REGION_SAMPLE_SOBOL` | the 2D Sobol sequence (van der Corput, then the direction numbers of x + 1), digitally shifted by the seed (seed 0: unshifted) |

The same seed gives the same points. `map` is thread-safe on one sampler.

## Precision and cost

- **Tolerance.** The area tolerance defaults to 1e-13 of the box's area,
  or 1e-6 in f32. A(x) and the area match Mathematica to about 1e-12 in f64.
- **Cost.** Creation is one union plus one integral per piece and slab.
  A point costs a few Newton steps, each one integral per piece crossing
  its slab: about 30 µs on the test regions.

## Tests

Test 98 (`tests/98_region_sample.c`) uses three regions against
Mathematica (`tests/reference/98_region_sample.wls`, Green's theorem at 30
digits):
- **Regions:** a crescent (a disc minus an overlapping disc, so two chords
  at some x), a disc with a hole, and a bean of four cubic Béziers with a
  dent.
- **Exact quantities:** the area, and A(x) at eight abscissas.
- **The map:** A(map(u, v).x) = u * area.
- **Containment:** every point lies inside the region.
- **Sobol points:** the centroid and the second moments.
- **Random points:** the distribution of x against A(x) / area, under the
  1 % Kolmogorov-Smirnov bound.

It runs in f64 and f32.
