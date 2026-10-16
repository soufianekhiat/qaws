# API Reference: qaws_region_sample.h

Points uniform in area over a region bounded by curves of any kind:
concave, with holes, under a fill rule. Everything comes from the boundary;
there is no mesh.

```c
qaws_region_sampler_desc d = { paths, path_count, QAWS_FILL_NON_ZERO, 0 };
d.density = rho; d.density_user = user;                 /* optional: importance sampling */
d.executor = &executor;                                 /* optional: parallel table and points */
qaws_region_sampler* s;
qaws_region_sampler_create(&s, &d);

qaws_region_sampler_get_area(s);
qaws_region_sampler_get_mass(s);                        /* the density's integral (the area without one) */
qaws_region_sampler_area_left(&a, s, x);            /* A(x): the area left of x */
qaws_region_sampler_contains(&inside, s, p);         /* the slabs' winding below p */
qaws_region_sampler_map(&p, s, u, v);               /* the unit square onto the region */
qaws_region_sampler_map_many(points, s, uv, n);
qaws_region_sampler_generate(points, s, n, QAWS_REGION_SAMPLE_BLUE_NOISE, seed);
qaws_region_sampler_generate(points, s, n, QAWS_REGION_SAMPLE_SOBOL | QAWS_REGION_SAMPLE_BOX, seed);
/* kinds: RANDOM, SOBOL, STRATIFIED, HALTON, R2 (through the map, or | BOX), BLUE_NOISE */
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
   - A(x) is the integral of the chord length l(x). At creation, l is
     tabulated in panels (the same table as a density's below, with rho =
     1): 8-point Gauss per panel, kept as the degree-7 Legendre series
     through those values. `u` locates the panel holding `u` times the
     table's total, and `x` inverts that panel's series by safeguarded
     Newton. A(x) itself, and `area_left`, stay exact by Green's theorem.
     The map agrees with them to about 1e-15 of the area, against 2e-13
     with Newton on the exact integrals, which this replaced.
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
| `QAWS_REGION_SAMPLE_SOBOL` | the 2D Sobol sequence (van der Corput, then the direction numbers of x + 1), Owen-scrambled by the seed (hashed nested uniform scrambling, Burley 2020; seed 0: plain) |
| `QAWS_REGION_SAMPLE_STRATIFIED` | correlated multi-jittered sampling (Kensler 2013): one point per row and column of a near-square m x n grid, any count |
| `QAWS_REGION_SAMPLE_HALTON` | Halton in bases 2 and 3, shifted by the seed (Cranley-Patterson) |
| `QAWS_REGION_SAMPLE_R2` | the R2 Kronecker sequence (the plastic number), shifted by the seed |
| `QAWS_REGION_SAMPLE_BLUE_NOISE` | weighted sample elimination (Yuksel 2015) from 5 n uniform points: evenly spaced in the region |

- **Unit-square patterns through the map.** By default, the first five are
  point sets of the unit square sent through the map. They are uniform in
  distribution, so their moments converge as each pattern should, and the
  count is exact by construction. But the map follows x, then the chords at
  x, so it is anisotropic: short chords compress y, and slowly growing area
  stretches x. A grid's rows are visibly sheared in thin parts of a region.
- **Over the box: `| QAWS_REGION_SAMPLE_BOX`.** OR-ed into a unit-square
  pattern, the pattern is laid over the region's bounding box and the
  points outside are dropped.
  - Restricting a uniform or low-discrepancy set to a subregion keeps it
    uniform and keeps its spacing, the same in every direction, so there
    is no shear.
  - A sequence (random, Sobol, Halton, R2) runs until n points fall inside.
  - For the stratified grid, every stratum inside is kept. Dropping strata
    leaves holes, and eliminating the crowded ones favours the boundary,
    since points there have fewer neighbours. Instead the grid's size is
    searched until the strata inside number n (or a few more, dropped at
    random).
  - The cost is box / area candidates per point kept.
  - The inside test is `qaws_region_sampler_contains`: the winding of the
    pieces below the point in its slab, with no curve winding numbers.
- **Blue noise.** This works on the points themselves, so they come out
  evenly spaced: the closest pair is above half the hexagonal packing radius
  of n points over the area.
- **Determinism.** The same seed gives the same points, with or without an
  executor: every point is a function of its index (random included: point
  i takes draws 2i + 1 and 2i + 2 of the seed's splitmix64 stream).
- **Thread safety.** `map` is thread-safe on one sampler.
- **Executor.** With `desc.executor` (`qaws_batch_executor`), creation
  tabulates the density's panels in chunks, and `generate` and `map_many`
  map points in chunks of 64. Over the box, a sequence's inside tests run
  in parallel blocks, and the points are then kept in index order. Blue
  noise maps its candidates and computes their radii in parallel; the
  elimination itself is serial.

## Density (importance sampling)

Set `desc.density` (and `density_user`), and points come in proportion to
rho(x, y) instead of uniformly. The density must be finite and >= 0, and
smooth for the map's accuracy. `get_mass` and `mass_left` return its
integral over the region and over the part left of x.

- **Through the map.** The map keeps mass.
  - `u` picks x by the mass left of it, M(x) = ∫ g, where g(x) is the
    density integrated along the chords at x.
  - `v` picks y by the mass along those chords.
  - M is tabulated at creation: 256 panels over the region's width, 8-point
    Gauss in each, and composite 8-point Gauss along the chords. The panels
    are uniform in s, with x = x_k + w (3s² − 2s³) in each slab: where the
    boundary turns vertical, chords grow like sqrt(x − x_k), and the
    smoothstep makes the integrand smooth again.
  - Each panel also keeps its integrand as a degree-7 polynomial: the
    Legendre series through its 8 Gauss values, whose integral is the
    panel's Gauss sum. x inverts this polynomial by safeguarded Newton,
    with no density call. Along the chords, y finds the Gauss segment
    holding its share, then inverts that segment's polynomial the same
    way: one pass of density calls along the chords, and one more over the
    segment.
  - `mass_left` keeps exact Gauss on the partial panel; the map agrees with
    it to about 1e-15 of the mass on smooth densities.
  - Every unit-square pattern follows the density, with its structure.
- **Over the box.** A candidate is kept when it is inside and a hash of its
  index falls under rho / rho_max. rho_max is the largest value met while
  tabulating, and on a 64 x 64 grid over the box, plus 25 %; it is an
  estimate, so a density with a sharp, narrow peak can exceed it.
- **Blue noise.** 5n candidates are drawn through the density's map. Each
  is eliminated against neighbours within its own radius, rmax ·
  sqrt(mean / rho), kept within a factor of 4 of rmax, so the spacing
  follows the density.

## Precision and cost

- **Tolerance.** The area tolerance defaults to 1e-13 of the box's area,
  or 1e-6 in f32. A(x) and the area match Mathematica to about 1e-12 in f64.
- **Cost.** Creation is one union plus one integral per piece and slab.
  The panel table adds about 5 ms. A uniform point costs one panel
  lookup, its series inverted, and the crossings at x: on the crescent,
  2.4 µs serial and 0.34 µs on 8 threads (26 µs with Newton on the exact
  integrals).
- **Density cost.** Creation is 256 panels x 8 chord integrals (a few ms).
  On the crescent, 20,000 Sobol points took 2 to 3 µs each, against 35 to
  57 µs with Newton on the integrals before the interpolants. With 8
  threads it was 0.4 to 0.6 µs. Under a narrow bump, the old Newton along
  the chords also missed the chord fraction v by up to 0.19; the
  interpolants hit it.

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
  - **Patterns:** the five unit-square patterns run through the map and over
    the box, plus blue noise. The sampler's inside test agrees with the
    winding numbers on every Monte Carlo point.
  - **Cells:** over a 6 x 6 grid, each cell's area is exact by clipping
    the region with the cell. Each pattern's samples estimate the same cell
    areas. The total variation stays under 0.6 % for the low-discrepancy
    and stratified patterns, about 2 % for random, and about 2 % for blue
    noise with a quarter of the points. No point falls outside.
  - **Density:** with rho = 2 + x² + y, the mass matches Green's theorem on
    the united region to 12 digits. Cell masses from Sobol (through the
    map) stay under 0.7 %, from blue noise about 2 %, and from stratified
    over the box about 2.5 %, because rejection by a density costs the
    stratification.
  - **Full size:** `QAWS_RST_FULL` runs the sweep five times larger.
- **Density against Mathematica:** for rho = 1 + x² + y over the three
  regions, the mass, the mass left of x at eight abscissas, the map keeping
  mass, and every pattern's weighted centroid (this needs the third-order
  moments, which the reference script also gives).
- **Executors:** chunks run backwards in uneven ranges, and on 4 threads.
  Uniform and with a density, every pattern's points and `map_many` are
  identical bit for bit to the serial run.

It runs in f64 and f32.
