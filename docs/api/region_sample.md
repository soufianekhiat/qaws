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
qaws_region_sampler_generate(s, n, QAWS_REGION_SAMPLE_BLUE_NOISE, seed, points);  /* or RANDOM, SOBOL, STRATIFIED, HALTON, R2 */
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

## Patterns

| `kind` | Points |
|---|---|
| `QAWS_REGION_SAMPLE_RANDOM` | independent uniform (u, v) from splitmix64, seeded |
| `QAWS_REGION_SAMPLE_SOBOL` | the 2D Sobol sequence (van der Corput, then the direction numbers of x + 1), digitally shifted by the seed (seed 0: unshifted) |
| `QAWS_REGION_SAMPLE_STRATIFIED` | correlated multi-jittered sampling (Kensler 2013): one point per row and column of a near-square m x n grid, any count |
| `QAWS_REGION_SAMPLE_HALTON` | Halton in bases 2 and 3, shifted by the seed (Cranley-Patterson) |
| `QAWS_REGION_SAMPLE_R2` | the R2 Kronecker sequence (the plastic number), shifted by the seed |
| `QAWS_REGION_SAMPLE_BLUE_NOISE` | weighted sample elimination (Yuksel 2015) from 5 n uniform points: evenly spaced in the region |

- **Unit-square patterns.** The first five are point sets of the unit
  square sent through the map. They are uniform in distribution, so their
  moments converge as each pattern should. But the map follows x, then the
  chords at x, so it is anisotropic: short chords compress y, and slowly
  growing area stretches x. A grid's rows are visibly sheared in thin parts
  of a region.
- **Blue noise.** This works on the points themselves, so they come out
  evenly spaced: the closest pair is above half the hexagonal packing radius
  of n points over the area.
- **Determinism.** The same seed gives the same points.
- **Thread safety.** `map` is thread-safe on one sampler.

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
- **Every pattern:** the points lie inside, and their centroid and second
  moments match (to 1.5e-2 of the size for random, 2e-3 for the others).
- **Blue noise:** the closest pair is above half the packing radius, ten
  times the random closest pair.
- **XOR:** two discs under even-odd: the area, and the sliver of the small
  disc outside the big one getting its share of the points.
- **Random points:** the distribution of x against A(x) / area, under the
  1 % Kolmogorov-Smirnov bound.
- **The sweep:** nine shapes against six patterns.
  - **Shapes:** an ellipse (NURBS, turned), disc minus disc, disc xor disc,
    a disc with a hole, the bean, a star polygon, the union of five
    overlapping circles, a wavy polar blob, and an ellipse with three holes.
  - **Areas:** the sampler's area matches the united region's area
    (`qaws_path_compute_area_2d`) and the closed forms to 12 digits.
  - **Area from samples:** points uniform in the box, located by winding
    numbers on the input paths (`qaws_region_locate_point_2d`), give the
    area within 4 sigma.
  - **Cells:** over a 6 x 6 grid, each cell's area is exact by clipping
    the region with the cell. Each pattern's samples estimate the same cell
    areas. The total variation stays under 0.6 % for the low-discrepancy
    and stratified patterns, about 2 % for random, and about 2 % for blue
    noise with a quarter of the points. No point falls outside.
  - **Full size:** `QAWS_RST_FULL` runs the sweep five times larger.

It runs in f64 and f32.
