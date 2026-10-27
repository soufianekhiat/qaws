# Meshing: generation, remeshing, quads, tets and hexes — design discussion

Status: draft for discussion, 2026-10-17; decisions taken 2026-10-18 (below).
It proposes the shape of the API, the algorithms behind each entry point,
and an order of work; the questions it raised are kept in section 13.

## Decisions (2026-10-18)

1. **Three mesh types**: one for 2D, one for 3D surfaces, one for volumes.
   The 2D and 3D ones are the existing `qaws_mesh_2d` (the planar mesher's
   output) and `qaws_mesh_3d` (the stack mesher's), generalized: created
   from arrays, editable, with attributes, features and footpoints. The
   volume type is new: `qaws_volume_mesh`. One internal half-edge core
   (dimension 2 or 3) serves both surface types.
2. **Editing operations are public** (split, collapse, flip, insert), with
   their validity checks.
3. **Footpoints are always kept and exposed.**
4. **PN is the default interpolant** of mesh references.
5. **Mesh cleaning is another arc.** This one lays out the modules for
   everything needed, each operation with a single method first behind a
   `method` field (the other methods to come; `UNSUPPORTED_OPERATION`
   until they do).
6. **Hexahedra include frame fields** (after the grid method).
7. **Branch `meshing`** (from `repulsive`), pushed as work lands, results
   shown as images along the way.
8. **Storage: a corner store**, as geogram's facet corners (and Atelier's
   structure): faces are runs of corners, each corner holding its vertex
   and the opposite corner across its edge (the corner table). A corner is
   a half-edge, so navigation is a half-edge view of the faces, and faces'
   adjacency, the dual graph, is stored (`qaws_mesh_3d_get_face_adjacency`,
   `qaws_mesh_3d_create_dual`). Volumes follow the same pattern (cell
   corners, cell facets' adjacency).
9. **Boundary with farq (2026-10-19).** farq (exact mesh Booleans) does
   everything exact and tolerance-free; qaws owns every tolerance and every
   quality judgement. qaws: welding near-duplicates, removing degenerate
   faces, hole filling, gap bridging, orientation repair by heuristics,
   scans, remeshing, quads, quality tetrahedra and hexahedra. farq:
   Booleans, self-union, intersection resolution, the exact defect network,
   exact orientation and winding-number solids of soups, exact Booleans on
   volume cells. So qaws's repair only has to make a surface closed: it does
   not resolve self-intersections (farq's self-union does, exactly). Neither
   library links the other; the caller composes them through arrays.
10. **Numerical Booleans in qaws (2026-10-21).** farq keeps exact Booleans
   (against Solidean); qaws adds fast, robust numerical ones (against
   trueform): filtered predicates with consistent fallbacks, no rational
   arithmetic, self-intersections resolved the same way. Decision 9 holds
   for repair's tolerances; qaws's Boolean output stays closed and
   manifold so the two compose.

## As built (2026-10-18)

On branch `meshing`, one method per operator (the others reserved behind
`method`, `UNSUPPORTED_OPERATION`); the API in [api/meshing.md](api/meshing.md):

| Step | Landed | Reserved / next |
|---|---|---|
| Core | the corner store (`qaws_mesh_core.h`), `qaws_mesh_3d` over it: creation from arrays (non-manifold input split), edits with checks, features, attributes, footpoints, quality, topology, dual; `qaws_mesh_2d` held in it at z = 0 (`qaws_mesh_2d_get_mesh`), made from arrays and edited there | the planar mesher working in the store directly (it keeps its own arrays and mirrors them) |
| References | PN (default) and LINEAR, creases, curvature, project, walk, attributes; Mathematica fixture | (all done 2026-10-27: PHONG, NAGATA, G1, SUBDIVISION, MLS) |
| Fields | size fields (uniform, curvature with tolerances and gradation, local feature size (shrinking balls, thickness and gaps, `cells_across`), vertices, callback, min); cross fields: LOCAL (Instant Meshes) and SMOOTHEST (Knoppel et al. 2013, n = 2, 4, 6, lines held, curvature_aligned) on a sparse LDL^T, indices; frame fields by a global spherical-harmonic solve (Ray et al. 2016), singular edges, combing | anisotropic metrics; meshable frame fields (singularity constraints) |
| Remeshing | LOCAL (Botsch-Kobbelt as pmp, geometry-central and Mmgs refine it: ideal valences, size-weighted relaxation, needles and caps, Delaunay flips; angle bounds after Wang, Yan et al. 2019; a distance guard as Mmgs's hausd), CVT (restricted Voronoi, Lloyd then L-BFGS, creases locked), features, attribute transfer | FIELD_ALIGNED (6-RoSy, on the smoothest field and the flow); simplification |
| Quads | INSTANT: Instant Meshes' fields or the smoothest global ones, whichever has fewer singularities; QuadriFlow's consistency by a min-cost flow, cones on lattice points; half cells paired, doublets removed; finishing by local quad operations and angle relaxation; quad-dominant or pure at the target spacing; GLOBAL: integer-grid map on the flow's integers, QEx extraction, INSTANT as fallback | PAIRING, SPLIT, LAYOUT (QuadWild); GLOBAL's own quantization |
| Generation | REMESH for surfaces and shells, snapped onto the exact faces, footpoints (face, u, v) | PARAMETRIC, RESTRICTED_DELAUNAY, STRUCTURED; trim curves followed |
| Volumes | `qaws_volume_mesh` (tetrahedra, hexahedra, prisms, pyramids; facets' adjacency stored) | attributes on volumes |
| Tetrahedra | DELAUNAY: exact lattice Bowyer-Watson, border recovery, BCC fill, refinement, flips and smoothing | ENVELOPE (fTetWild), STUFFING; Klingner-Shewchuk improvement; CDT recovery (Diazzi et al. 2023) |
| Hexahedra | GRID: staircase padded to the surface, creases and corners kept; FRAME_FIELD: integer-grid map of a field without singular edges, its lattice cubes; POLYCUBE: a polycube deformation feeding that map | FRAME_FIELD through singular edges (HexEx extraction); octree |
| Measures | `qaws_mesh_measure.h`. Shape in papers' forms: Q, angles, shares under 30 and over 90, V567, scaled Jacobian, AD, ED; ideal valence from the corners' angles (a creased quad cube 0 irregular, against 8 from `qaws_mesh_quality`). Hausdorff certified both ways by branch and bound: Bezier nets of the difference over pieces matched through unfoldings round a shared vertex, and a Lipschitz bound through the net where the closest place jumps across a concave edge, zero-area fan triangles left to their neighbours, the PN projection polished along and across patch sides; gap within 1e-4 of the diagonal. 10.7k faces: shape 1.8 ms, deviation 0.5 to 0.8 s (Release, contended). Test 130 (with a PN saddle against its own mesh, witnesses checked by dense sampling, and a straight-cornered pentagon); figure 9 and its stdout table, the baseline for 131..135 | mesh-to-PN lower bound only up to the (polished) patch projection; a triangle whose corners share no vertex star falls back to one triangle at a time (linear convergence until split); report the certified deviation from `rm_end` and quadrangulate |
| Smoothest cross fields | SMOOTHEST direction fields (Knoppel et al. 2013; n = 2, 4, 6; lines held, corners free; curvature_aligned), on a new sparse LDL^T (nested dissection, supernodal; 50k vertices in 1.3 s); INSTANT tries the global fields, relaxed at the lattice's scale (MBO, one factor for both) where lines hold them: can all quads 22 to 10 irregular (6 on the top cap, 4 on the bottom), blob 12, torus 0. Figure 4 against the prep build: blob pure 1384 to 1278 quads (irregular 12 both; mean edge 0.115 for 0.11, about 0.110 before by the quad count), blob quad-dominant 571q+5t+1p to 555q+6t+0p (deviation 2.37e-2 to 2.01e-2), can pure 758 to 764 quads with 82 to 76 crease edges, torus unchanged; figure 5's wavy quads 476 to 368; figure 3 unchanged; figure 4 3.5-4.1 s to 4.3-5.2 s (loaded machine). Indices whole on coarse faces; each connected part solved by its own lines or eigenvector; test 133, figure 12 | the holonomy-corrected connection; the smoothest field's coinciding zeros parted in the public API; relaxation at 2x the lattice gives rounder O-grids (cones at 0.6 R, not 0.8 R) but 13 irregular with the greedy consistency: retry with 134's flow |
| CVT remeshing | `QAWS_REMESH_CVT`: restricted Voronoi diagram on sub-triangles of the reference (bisector clipping to the radius of security, per region; a locked seed's cell started on every side of its creases, so the cells tile the reference), density size^-4; Lloyd 5 then L-BFGS 30 (H0 the cells' Lloyd step, scaled); crease and border vertices spread by a 1D Lloyd and locked; vertex count from the size; started from LOCAL's full 10 iterations; every move guarded (no face turned by more than 60 degrees or against the reference's normals, seeds walked along the reference, never onto a thin part's other side); a coarse level first (a quarter of the seeds, every other feature-line vertex merged, then 1-to-4 restoring LOCAL's line counts), dropped where the split would fold; LOCAL's connectivity kept extrinsic-Delaunay by flips. Icosphere 95.2% valence 6 (LOCAL 77.6%), blob 5.4% irregular (22.4%), angles 36..99 on every closed or creased case tested (slab 0.02 thick 18..97, LOCAL 18..91 since its cap flips (16..129 before), the max asserted within 10 degrees of it; ellipsoid 0.3 thick at 0.06 39..99 with no face turned over, LOCAL 34..102), edge CV half LOCAL's, mean edge within 1% of target, crease counts equal to LOCAL's; Release, under load: blob 2.4k vertices about 2 s, torus 2.5k about 2 s, can of 48 sides 1.3k 1.0 s, 8 sides 0.7 s; test 132, figure 11 | exact RDT extraction (no LOCAL connectivity); Lp-CVT and anisotropic metrics; Newton with the exact Hessian (tried: no better per second than L-BFGS with the coarse level); cells limited to sub-triangles facing the seed (tried: no measurable change once moves are guarded); cells reaching across thin parts; a smooth rim about as thin as the size (folds as LOCAL's) |
| LOCAL quality | LOCAL by default: size-weighted area-centroid relaxation (full step, half on refusal; 3 sweeps in the last third only with a size field, angle bounds or a tolerance, 1 otherwise, since 3 moved INSTANT's dense mesh: figure 4 blob pure 12 to 16 singular, can 22 to 11, test 125's torus 529 of 530 quads to 479 of 482), collapses onto the midpoint (crease midpoint along a line) with Zheng and Lv's 120-degree guard, flips toward per-sector ideal valences (twice per iteration) with a 30-degree guard, caps at 120, needles under 15, a last Delaunay queue, a convergence stop; options `min_angle`/`max_angle` (Wang et al. 2019's angle-bounded pass: flips, pattern-search relocation in 16 directions, splits with the best settled flip, a sector with too few faces given one, collapses; 2k + 16 attempts a round; exempt sectors) and `tolerance` (Mmgs-style: sizes from the curvature at a chord gap of 3/4 of it, a guard on every operation, lenient in the loop where the faces replaced strayed already, quarter points sampled after the loop, longest-edge splits of the faces still farther). Can 33.9..108.0 degrees, blob mean edge +1.5%, 82.5% valence 6, Delaunay 100%, adaptive ratio 0.386; [35, 86] reached on icosphere, blob, adaptive sphere, can, cube at 0.08, 0.15, 0.25, 0.6 and an open square, [30, 90] on the can; tolerance held at every sample (2e-3: +66% vertices; 5e-4: 13.9..141.3 degrees; 2e-4: 9.8..144.9, 0..180 before; with [35, 86] the bounds held too). Figure 4 with the gate: blob quad-dominant 7 singular (8 before), pure 12 (12), can 14 (22), torus 0; test 125's torus 525 of 525. Figure 3 1.05x its time; angle pass 1.1x on the blob. Test 131, figure 10 | angle bounds by feasible regions (Wang's exact relocation) rather than a pattern search: a strip narrower than the target and the faces next to a 20-degree corner stay out; the tolerance samples certified by 130's two-sided Hausdorff (dense sampling finds 2.38e-3 at 2e-3); the tolerance costs 5 to 10 times the default at the same density; the angle pass up to 80 times a tiny remesh where it fights to the last corners (the cube at 0.6, 1.5 s) |
| Quads by flow | INSTANT's consistency by a min-cost flow (generic `qaws_internal_flow`: successive shortest paths), cones on lattice points by parity paths and onto vertices, flip removal, positions refitted, the mesh read off half cells paired across their diagonals, `pure` at the target spacing (fallback: twice it, split); test 134, means over six inputs scaled by 1 + k 1e-12 (before: e280d24 on the same six): blob quad-dominant 0 non-quads (were 6 to 8) but 11.0 irregular (were 8.3; target 10 missed, floor 12), blob pure 9.3 at 0.11 (were 13.2; target 8, floor 10), can 11.0 (were 17.2), torus 0 with 697 quads, wavy patch 0 (were 4); figure 13 (figure 4's blob) 11 and 12 (were 8 and 12); coarse targets about a torus's tube valid (a pillow, two half cells on the same three lattice points, made a quad with a vertex twice; now dropped as a fold, the hole filled), the fallback split of the first mesh when twice the target leaves none (a disk at 0.9), open borders capped far less often (a loop all on its lines left open; a short loop holding an inside vertex still can be); 0.5 s on 50k dense faces (target 0.2) | a cone no vertex reaches moved by routed units; a global position solve (133's LDLT); odd cones by min-cost matching; GLOBAL, LAYOUT and FIELD_ALIGNED on the same flow |
| Quad finishing | INSTANT's quads finished (`qaws_quad_improve.c`): corner-store rotations, diagonal collapses, doublet removals (`qaws_mesh_core.h`); E against a feature-aware ideal valence lowered greedily, each operation judged on a relaxed copy of its patch (scaled Jacobian, normals), a rotation and a second operation where none alone helps (judged against the quads before the rotation); relaxation down a shape energy (angles against the regular polygons', lengths against the target at weight 0.1; Zhou-Shimada target or a Newton step, walked; lines along the creases) guarded by the scaled Jacobian, the shortest edge (not under min(before, L/2)) and the quads' edge ratio (not under min(before, 0.3)), worst quads searched. Test 135, figure 14: can 22 to 20 irregular, AD 5.8 to 4.7, SJ mean 0.974; blob 12 kept, AD 5.7 to 4.8; wavy patch 4 kept, AD 4.6 to 3.7; torus 12 to 11, SJ min -0.55 to 0.74; thin quads (edge ratio under 0.2) and the shortest edge at the baseline or better, ED 0.124 to 0.146 (can), 0.105 to 0.118 (blob); about 20% more time; on the integrated branch (the flow's extraction leaves quads turned inside out on the blob, SJ min -0.82 after the plain relaxation) those are collapsed: blob 15 to 12 irregular, SJ min 0.66, AD 3.0 to 2.4; can 8 kept, AD 3.8 to 2.4, ED 0.066 to 0.099; torus 0 kept, SJ min 0.98 | targets missed: can <= 16, blob <= 10, wavy <= 2 irregular (E-lowering single operations exist but each adds vertices off their ideal, a 6 into two 5s and the like, and is refused on purpose; the gains come from pairs; left are 3-5 pairs two edges apart and lone cones), can SJ min >= 0.5 (rim slivers need a collapse along the crease), blob AD 20% lower (16%), ED at the baseline (angles bought with lengths; within 20% asserted, on the integrated can 1.6 times); deviation: candidates collected once per pass, not again round each commit; dislocation moves along grid lines toward partners, crease-aware collapses |
| Quads by integer-grid maps | GLOBAL (`qaws_quad_global.c`): per dense vertex its place in its cell, every edge a quarter turn and a whole shift from INSTANT's flow (no cut graph); cones on lattice points (on a vertex, else a face split round its centroid), corners charted from their neighbours, steps between creases at a corner read again in the face's plane; one LDL^T solve of the MIQ energy with lines on whole isolines; folds undone by kernel centres; QEx-style extraction (lattice points registered round each triangle, unit isolines traced through the transitions nudged to their left, faces by left quarter turns, all quads or fallback); INSTANT's finishing; report iterations 1 or 0. Test 136, figure 15 (key f): torus AD 3.13 (INSTANT 3.60), box AD 0, blob 8 irregular SJ min 0.815 (0.772) AD 2.74 (3.27), can SJ min 0.853 (0.801), wavy patch AD 2.97 (3.22); 14 of 20 by the map over targets 0.75 to 2.2 times. Its own quantization when the flow's integers give no map: the cones (singular faces, corners all the way round) as nodes of a spanning tree whose paths cross no crease or border, whole steps along node-to-parent paths and singular faces' paths to the border as the unknowns, the map's energy over them a dense quadratic (one back substitution per step on the base factor), greedy rounding as MIQ then whole-step descent, cones then onto their nearest free vertex; corners' places carried round their triangles from the turn's fixed point; steps between creases read in the right plane (face normals normalized); a vertex wound twice is no map. Test 136's sweep: 20 of 20 by the map (can at 1.25 SJ min 0.860 against INSTANT's 0.178, wavy patch at 1.25 2 irregular against 3), three of them by default, all 20 heavy; figure 15's can at 0.175 by the quantized map | crease isovalues and torus periods quantized too (a node per line); both integer sets tried and the better kept (the box at 0.75 keeps the flow's poorer map); a corner off the lattice moved by an odd step; Bi-MDF on a motorcycle layout for a guaranteed quantization |
| Repair for farq | `qaws_mesh_repair` (`qaws_mesh_repair.h`, `qaws_mesh_repair.c`), decision 9's side of the split: arrays in, a closed outward `qaws_mesh_3d` of triangles out with a `source_face` attribute (-1 for fills). Validation (bad index, face under 3, position not finite: INVALID_ARGUMENT naming the first), polygons ear-clipped; welding within `weld_tolerance` (1e-6 D, union-find over the grid broad phase, groups at their means); degenerate faces under `degenerate_tolerance` (needles collapsed, caps removed with the faces across their long edge split, no hole); duplicates kept once, opposite pairs cancelled; gaps within `gap_tolerance` (1e-3 D: boundary vertices merged, or snapped onto boundary edges whose faces are split, until the borders meet); small parts removed; orientation by propagation over manifold edges (non-orientable parts reported and left); holes up to `hole_perimeter` filled by Liepa's programme (max dihedral, then area; 500-vertex loops), refined (centroid splits, Delaunay flips) and faired (scale-dependent bi-Laplacian, one LDL^T); parts outward by signed volume, cavities inward by winding-number depth (or removed); non-manifold edges and self-intersections left for farq, `closed` when every directed edge cancels. Test 139 (1.5 s): cube soup to 8 vertices, sphere soup off by 1e-7 rewelded, three sphere holes filled to 1.4% of the area removed (flat fills 8.6% under), fill edges +12% of the surrounding, smallest fill angle 25.9; box halves at 6 and 5 cells 3e-4 apart stitched closed without fills; turned faces turned back; cavity and island; cap, needles, zero-area face; duplicates; Moebius strip; two cubes on one edge. Figure 18 (key i): a broken torus soup (3 holes, a gap of 0.004, 3 faces in 10 turned, 46 copies) closed in 0.02 s | loops past 500 vertices by a fan (halve by the shortest diagonal); fills ignore features; gaps wider than the tolerance are filled, not bridged by strips; nesting judged at 16 samples; self-intersections are farq's |

Tests 121 to 139; showcase figures 1 to 18 (`qaws_meshing_showcase`, keys 1 to 9 then a to i).

### Next, beyond the reserved methods

User guidance (what a modelling tool's user steers with):

| Input | Today | Next |
|---|---|---|
| Edges to keep | `qaws_mesh_3d_set_feature_edges`, honoured by LOCAL, CVT, quads and generation; corners where lines meet or end | per-line options (keep straight) |
| Density paint | `qaws_size_field_create_vertices` (a size per reference vertex), callback, min; LOCAL, CVT and quads (`qaws_quadrangulate_desc.size`: the lattice at a local spacing, dislocations kept where it changes; test 137) | smooth gradings without dislocations (cones placed for them), the finishing graded too |
| Directions | `guide_directions`, `guide_weights` (pulls and holds per vertex) in the direction field (LOCAL, SMOOTHEST) and quads | strokes as curves (sampled for the user), symmetry |
| Singularities | `cone_vertices`, `cone_indices` (placed: SMOOTHEST by a connection adjusted to them, Poincare-Hopf checked on closed surfaces; quads with irregular vertices exactly on them, INSTANT and GLOBAL) and `attractors` with strengths and a reach (SMOOTHEST: a share of a cone's turn round the point; LOCAL: lighter smoothing) in the direction field and quads; test 137, figure 16 | a cone budget, cones on feature lines, chosen turns round handles |
| Edge loops | `guide_edges` in quads: segments the quads' edges run along (one lattice line, quad vertices on them) | crossings pinned, loops through cones |
| Pinned vertices | `pinned` in quads and in LOCAL and CVT remeshing: kept where they are, as corners | pins that may slide |

Landed (quality and reach): curved elements for surfaces (`qaws_mesh_generate_nodes_surface` / `_shell`, test 157, figure 20), boundary layers (`qaws_volume_mesh_grow_layers`, test 156, figure 19), metric-driven adaptation (`qaws_remesh_desc.metric`, `qaws_metric_field_create_hessian`, test 158, figure 21), sizing by local feature size (`qaws_size_field_create_feature_size`, test 138, figure 17), input repair (`qaws_mesh_repair`, test 139, figure 18), faster tests (a default tier and `QAWS_TEST_HEAVY=1`). Their details are in the As built table and in [api/meshing.md](api/meshing.md).

### Roadmap (2026-10-21)

All 19 items landed or started by 2026-10-22; what is left of them moved to the enrichments roadmap below. In order of work. Each item says what it unblocks; the order follows those dependencies first, then what users ask for first.

Foundations:

1. **Spatial queries on meshes: landed.** `qaws_mesh_bvh.h` (test 159, figure 22 key m): binned SAH, closest points and point location, watertight rays, fast winding numbers to second order, signed distance, refit, arrays or any mesh type. Box and ball overlap queries landed. Next: reinsertion optimization (as tinyBVH), TLAS instancing. Planned as: A public BVH over `qaws_mesh_3d`, `qaws_mesh_2d` and `qaws_volume_mesh`: closest point, ray casting, signed distance, inside tests, and the generalized winding number as a query (fast winding numbers, Barill et al. 2018, for soups and point clouds). The repair, watertight, Boolean, connector and offset items all build on it.
2. **Quality report: landed.** `qaws_quality_report.h` (test 160, figure 23 key n). Next: measures on curved (high-order) elements and per-solver presets. Planned as: One entry for 2D, 3D and volume meshes over the measures that exist (`qaws_mesh_measure.h`, `qaws_mesh_3d_measure_quality`, `qaws_volume_mesh_measure_quality`): histograms, the worst elements with their indices, and checks against a solver's thresholds (minimum angle, scaled Jacobian, aspect ratio). Cheap, and every later item reports through it.
3. **Repair, completed: landed** (test 139, figure 24 key o): creases through holes, long loops halved, unused vertices, volumes, non-manifold split, exact self-intersection finding. Gaps wider than the tolerance are bridged by the mesh connector (item 13). Planned as: **Repair, completed.** On top of `qaws_mesh_repair`: non-manifold edges and vertices split (or reported), isolated and floating vertices removed, small floating parts by volume as well as count, self-intersections detected and reported (resolved by item 5), holes filled along features, gaps bridged by strips past the tolerance, loops past 500 vertices.
4. **Watertight from any input: landed** (`qaws_mesh_make_watertight`, test 161, figure 25 key p): triangle input; point clouds landed (`qaws_mesh_make_watertight_points`, point dipoles); next: an adaptive (octree) grid, manifold dual contouring. Planned as: **Watertight from any input.** A closed manifold surface from soups, scans and point clouds: the generalized winding number field (item 1) or a signed distance contoured by dual contouring (sharp features kept), and an envelope-bounded variant (as fTetWild's) when the result must stay within a distance of the input.

Robust numerical Booleans (decision 10):

5. **Self-intersections resolved numerically: landed** (`qaws_mesh_resolve_intersections`, test 162, figure 26 key q); next: the face triangulation's point location and cut recovery are quadratic per face (fine for the cuts of a face, to revisit at scale, item 7). Planned as: Triangle-triangle intersections by filtered predicates with a consistent fallback (no rational arithmetic), coplanar and near-coplanar pairs handled together, the arrangement made conforming, the result snapped and checked.
6. **Mesh Booleans: landed** (`qaws_mesh_boolean`, test 163, figure 27 key r): union, intersection, difference, self-union, patches classified once by winding numbers. Next: the benchmark corpus (Thingi10K and the trueform cases, read by a harness outside the library), then speed under item 7 (the per-face setup is global and serial: 0.86 s for a 130k-triangle union). Planned as: Union, intersection, difference and self-union on closed meshes, by winding numbers on the resolved arrangement; output closed and manifold so farq can take it when an exact result is needed. A benchmark corpus (Thingi10K and the trueform cases) from the first commit, measuring time, failures and output quality.
7. **Scale: started.** Resolution and Booleans parallel through the executor, fitted pair grid (522k triangles resolved in 0.55 s on 8 threads, 0.77 s on one); BVH AVX2 with runtime CPU dispatch in progress. Next: the serial setup (weld, grid build), the Boolean's BVHs and result mesh, field solves. Planned as: Benchmarks past 10^6 elements, parallel field solves and local operations, on the executor interface the batches use. Placed here because the Booleans compete on speed.

Mesh processing:

8. **Geodesics on meshes: landed** (`qaws_mesh_geodesic`, test 164, figure 28 key s): the heat method factored once, paths shortened through their triangles, the intrinsic Delaunay triangulation for stretched and obtuse triangles. Exact polyhedral distances by ICH windows (`qaws_mesh_3d_geodesic_exact`). Exact paths: the funnel path down the exact distances. Next: the windows' cost bounded on large meshes (Xin and Wang's later variants). Planned as: The heat method (Crane et al. 2013) on the sparse LDL^T, then exact polyhedral geodesics where needed; geodesic paths traced with `qaws_mesh_reference_walk`. Needed by items 12, 13 and 14.
9. **Smoothing and fairing: landed** (`qaws_mesh_3d_smooth`, `qaws_mesh_3d_fair`, test 165, figure 29 key t). Next: the repair's own fairing moved onto it, anisotropic denoising landed (`QAWS_SMOOTH_BILATERAL`, test 165). Planned as: Laplacian, Taubin and bi-Laplacian smoothing with features held; fair patches for holes and bridges (the repair's fairing made public).
10. **Deformation and parameterization: landed** (`qaws_mesh_3d_deform`, `qaws_mesh_3d_parameterize`, test 166, figure 30 key u): ARAP, Laplacian and biharmonic editing; LSCM and ARAP flattening with seams. A prepared deformer landed (`qaws_mesh_deformer_*`). Next: BFF (reserved), seams chosen automatically, chart packing. Planned as: ARAP deformation with handles, Laplacian and biharmonic editing; ARAP, LSCM and BFF parameterization with seams cut. One local-global solver serves both.
11. **Decimation: landed** (`qaws_mesh_3d_decimate`, test 167, figure 32 key w): quadrics, a count or an error bound, borders and features held, certified deviation. Attribute seams kept as features. Next: attribute errors in the quadric (Hoppe 1999), memoryless simplification, the executor. Planned as: Quadric error metrics (Garland-Heckbert) with attributes, seams and features, and the certified deviation of item 2's report; a target count or error.
12. **Remeshing from attributes: landed** (`qaws_mesh_3d_cut_isolines`, metric intersection and gradation, test 168, figure 34 key y; metric CVT, test 158). Metric CVT measures each Voronoi cell in its seed's metric (Du and Wang 2005), density sqrt det, Delaunay flips in the metric. Next: the true anisotropic Voronoi (bisectors from both seeds' metrics, the metric's derivative in the gradient). Planned as: Isolines of an attribute kept as feature lines; metrics from several attributes combined by metric intersection; metric gradation; metric CVT (the follow-ups of metric adaptation).
13. **Mesh connector: landed** (`qaws_mesh_3d_connect`, `qaws_mesh_2d_connect`, `qaws_volume_mesh_connect`, test 169, figures 35 and 36 keys z, A): loops or chains matched by arc length, a Hermite bridge tangent to both meshes, graded, quads or triangles, attributes and features carried, the bridge measured; volume patches of the same triangulation joined by prisms and hexahedra. Next: volume patches of different triangulations (one remeshed to the other's border), attributes on planar and volume bridges. Planned as: Two meshes (2D, 3D surfaces, volumes), a boundary selected on each, joined by a bridge: boundaries matched by arc length, a strip or a tube meshed between them (quads, triangles, or prisms and hexahedra for volumes), faired (item 9), graded to the sizes on both sides, checked by item 2.
14. **Repulsive curves on surfaces with geodesic distances: landed** (`qaws_repulsive_mesh_curves.h`, test 170, figure 37 key B): curves on a mesh's triangles, each pair's distance the larger of the chord and the heat method's geodesic, steps walked across faces; patches and shells through a triangulation. Next: the parametric solver (R5) measuring geodesics on its patch directly, a hierarchical (Barnes-Hut) mode for long curves.
15. **Offsets of meshes: landed** (`qaws_mesh_3d_offset`, test 171, figure 39 key D): offset surfaces, shells and thickened sheets, contoured from the exact signed distance or moved along mitred directions with the self-union resolving overlaps. Next: an adaptive (octree) grid for fine contoured offsets, rounded joins for MOVE at convex corners, variable distances by contouring (per vertex landed for MOVE).
16. **Subdivision on `qaws_mesh_3d`: landed** (`qaws_mesh_3d_subdivide`, test 172, figure 40 key E): Loop and Catmull-Clark levels, feature edges as creases (semi-sharp by a sharpness per edge), corners, attributes by their rules. Limit positions landed (`limit`), limit normals landed (`limit_normals`, exact at smooth vertices). The modified butterfly landed (`QAWS_MESH_SUBDIVIDE_BUTTERFLY`, figure 50 key O). Next: the creases' tangent masks for exact normals on them, sqrt(3), adaptive (feature-driven) levels.

Volumes and reserved methods:

17. **Volume improvement: landed.** Border kept: `keep_border` in `qaws_tetrahedralize_desc` (test 175, figure 43 key H), the recovery's border points moved inside (Steiner point suppression rather than Diazzi et al.'s CDT; the CDT itself stays a next step for points no inside point frees). Landed: `qaws_volume_mesh_improve` (test 173, figure 41 key F): smoothing, edge removal (rings of 3 to 10) and face removal, the border kept; a shaken grid from 0.8 to 40 degrees. Landed: `qaws_volume_mesh_quadratic_nodes` (test 174, figure 42 key G): P2 tetrahedra and Q2 hexahedra on a mesh reference, Jacobians sampled and relaxed. Tetrahedra's and hexahedra's Jacobians bounded (Bernstein). Next: cubic nodes. Planned as: Sliver removal, tetrahedral smoothing and flips (Klingner-Shewchuk), CDT boundary recovery (Diazzi et al. 2023); curved tetrahedra and hexahedra at the border with a Jacobian check, nodes on a mesh reference.
18. **Hexahedra: started.** Landed: POLYCUBE (test 176, figure 44 key I): a polycube deformation (Gregson et al. 2011) feeding FRAME_FIELD's integer-grid map; spheres, rounded boxes, rings of genus one. Next: labels by graph cuts (no corner across a bent bar), FRAME_FIELD through singular edges (HexEx extraction), meshable frame fields.
19. **Quads: started.** Landed: PAIRING and SPLIT (test 177, figure 45 key J); LAYOUT as GLOBAL's quads parted by their motorcycle graph (test 178, figure 46 key K; QuadWild's own patch tracing and quantization stay next). Planned as: GLOBAL's quantization of crease isovalues and torus periods, LAYOUT (QuadWild), the cone budget and other guidance follow-ups of the table above.

### Roadmap (2026-10-22): enrichments

The features above, deepened. In order of work: bounded items with an exact check first, then larger ones, research-scale last. Each says what it adds and how it is checked.

Subdivision, geodesics, decimation, offsets (bounded, exact checks):

1. **Butterfly border rules: landed (checked).** Zorin, Schröder and Sweldens (1996, rule 4) reflect the stencils' vertices beyond the border into virtual ones, which the butterfly already did; test 172 now holds it: a sheet of valences 4 and 8 stays in its plane to 4e-16 up to the border, a cubic along the border follows the four-point curve to 1e-16. Next: the refined border rules of Zorin's thesis (a border vertex's valence in the stencil), for smoothness at irregular border vertices. Done (2026-10-22): Zorin's extraordinary crease rule (SIGGRAPH 2000 notes, figure 4.7; their typeset theta and border signs corrected by the masks' affine and linear checks), on each side of a crease or border vertex with n != 3 triangles; the characteristic map's spokes halved exactly (the reflection: 0.09 off). The regular interior-crease and crease-crease masks (figure 4.6) not taken: the reflection keeps planes and the four-point border exactly there.
2. **Geodesics up to a radius: landed** (`max_distance` of `qaws_mesh_3d_geodesic_exact`, test 164): the windows stop once the cheapest left lies beyond it; within it identical to the full run, beyond HUGE_VAL; the hill up to 0.5 in 14 ms against 525 ms.
3. **Nearest sources: landed** (`qaws_mesh_3d_geodesic_nearest`, test 164, figure 28 key s): per vertex its nearest source from either method's distances, by Dijkstra's parent in distance order; two sources on a grid part it along their bisector. Next: the regions' borders placed inside the triangles (per source distances where two meet). Done (2026-10-26): `qaws_mesh_3d_geodesic_borders`, from each source's exact distances: the squared distances' difference is linear on a triangle unfolded, so the crossings and junctions are exact within one window per source; a grid's ends equidistant to 4e-16, the junction the circumcentre, a cube's ends equal by `exact_at` to 3e-15.
4. **Distances anywhere: landed** (`qaws_mesh_3d_geodesic_exact_at`, test 164): the windows kept, a point's distance the least of the windows entering its triangle that see it, its triangle's vertices, a source in it; inside the grid's and the cube's triangles to 6e-16.
5. **Exact geodesics measured: landed** (docs/api/meshing.md): flat grids grow as n^1.6 (0.59 s at 20k), a hill as n^1.9 (98 s at 100k) through its saddles' pseudo-sources; the windows' slots reused once processed, memory 59 MB at 100k (was 460 MB at 20k). Next: fewer pseudo-source windows at saddles (a saddle's windows only past its shadow). Tried and dropped (2026-10-26): a saddle emitting only into its shadow (pi to its angle sum less pi past its path's arrival, the arrival kept per vertex) left the corners of a bumpy grid unreached: ICH's vertex filter drops windows whose region only the saddles' all-round windows then cover; and on a hill of 14,641 vertices it saved nothing (1.10 s either way). Pruning there would need the filter and the shadow decided together.
6. **Attributes in the quadric: landed** (`attribute_weight` of `qaws_decimate_desc`, Hoppe 1999, test 167, figure 51 key P): per vertex LINEAR and NORMALIZED attributes in the quadric, the position and the attributes placed by one solve; a sheet's colour step kept exactly where plain quadrics smear it (rms 0.265), a linear colour exact. Next: corner attributes (UVs across seams, Hoppe's wedges). Done (2026-10-26): wedges round each collapse (equal values at a vertex, joined across the edge's faces), each valued at the new position by interpolation on its nearest old triangle; a seamed grid's UVs exact through 512 -> 31 faces (before: 1.4 off). Not taken: the wedges' quadrics in the cost (the seams are features already; a UV-stretching collapse is not charged).
7. **Rounded offsets: landed** (`QAWS_MESH_OFFSET_ROUND`, test 171, figure 52 key Q): the points within |d| as a union of closed pieces (planar groups' slabs, convex edges' cylinder wedges, convex corners' cones over their normals' spherical hull, balls at saddles) by the Booleans; the cube's vertices within the chords of 0.1, its erosion exact. A convolution surface (faces, strips, caps) was tried first: its winding numbers go wrong at fully reflex corners (decided against). Next: the pieces' cost (a cube 0.6 s). Done in part (2026-10-26): the Boolean's patches sampled through a BVH past 2e6 solid angles (25% off the cube), `executor` in the offset's desc (the cube 0.68 s to 0.36 s on 4 threads). Left: the resolution of the pieces' 3,904 crossing pairs (0.14 ms each), fewer overlaps between the pieces.
8. **Exact limit normals on creases and borders: deferred.** A crease vertex's limit surface has a tangent plane per side, so its normals must be per corner (an API change from the vertex attribute), and our creases follow Hoppe's rules, which have no published tangent masks (Biermann's hold for Biermann's rules): it needs the left eigenvectors of the local subdivision matrix per side. Taken up with item 10.
9. **The sqrt(3) scheme: landed** (`QAWS_MESH_SUBDIVIDE_SQRT3`, Kobbelt 2000, test 172, figure 50 key O): centroids, the relaxation, the old edges flipped; counts, exact values, Euler, a linear colour exact. Next: Kobbelt's ternary border rule every second level (the border now keeps its polyline), limit positions. Ternary rule done (2026-10-22): borders and manifold creases cut in three every second level by the ternary cubic B-spline, chain ends interpolated (reflected), the triangle on a cut edge fanned in three; a sheet's 185 vertices, its border on the rule to 2e-16; the creased octahedron closed, 36 creases, its volume kept.
10. **Limit normals for the butterfly: deferred** with item 8: the same eigen-analysis, on the butterfly's two-ring stencil.

Volumes, points, smoothing (larger):

11. **Cubic nodes: landed** (`qaws_volume_mesh_cubic_nodes`, test 174, figure 53 key R): P3 tetrahedra and Q3 hexahedra from one degree-p implementation, edge nodes at equal arc lengths, inside faces by Coons and hexahedra's own nodes by Gordon and Hall, the Jacobians certified (84 and 729 Bernstein coefficients, hexahedra subdivided where loose); a coarse ball's border 0.0006 off the sphere against Q2's 0.0015. Next: simplex subdivision for the tetrahedra's bound. Done (2026-10-26): longest-edge bisection of the reference tetrahedron, 6 deep at most, only where the bound fails and the samples pass; the cubic sphere 5 relaxed cells to 3.
12. **Untangling curved elements: landed** (`untangled` in the report, test 174, figure 53): an element under `min_jacobian` has its free nodes moved by pattern search on the least sampled Jacobian round them until its certified bound passes, before any relaxing; the O-grid balls' corner cells kept on the sphere. Next: an optimizer with the bound's gradient (now coordinate steps), the border nodes sliding on the reference too. Done in part (2026-10-26): border nodes slide (tangent-plane steps, reprojected) for cells the free nodes cannot fix, and a move may not fail another cell round the node: the spheres 6 relaxed to 5, the cubic one's bound 0.167 to 0.178, more border nodes on the sphere (354 to 359). The gradient was tried and dropped: the least sampled Jacobian is nonsmooth (a min over samples), steps up its one-sided gradient lowered other samples and untangled none of the O-grid ball's cells; the O-grid ball's two relaxed cells stay.
13. **Points without normals: landed** (`normals` NULL in `qaws_mesh_make_watertight_points`, test 161, figure 48): PCA of the 10 nearest, oriented along a minimum spanning tree from each part's top, turned by the volume's sign; a sphere's, a cube's and a torus's volumes within 0.5% of the true normals'. Next: cavities (orientation by winding numbers per part). Done (2026-10-26): each part outward of its own volume, then turned when inside an odd number of others (their dipole winding numbers at its highest point): a hollow sphere's volume 3.6783 against its true normals' 3.6782.
14. **Denoised corners: landed** (test 165): BILATERAL's vertex update solved in least squares over its filtered planes (a pseudo-inverse where they span), corners to where their planes meet. Measured: the corners were not rounded before either, their error the three noisy planes' (0.0030 against faces' 0.0021 on a cube shaken by 0.01; the mean gave 0.0032): the check "within the noise's tenth" was beyond what the planes allow.
15. **Memoryless simplification: landed** (`QAWS_DECIMATE_MEMORYLESS`, test 167): the volume constraint and the swept volumes' optimization from the current faces, the border's and features' planes; the icosphere's volume kept to 1e-12, its deviation 1.4 times the quadrics'. Next: Lindstrom and Turk's boundary-area and triangle-shape terms (the regularization now pulls to the edge's middle). Border done (2026-10-26): the chain's area vector kept (two constraints by Gram-Schmidt against the volume's): a wavy flat sheet's area to 1e-15 (planes alone 6e-10). Shape tried and dropped: to the ring's centroid, deviation 0.0227 to 0.0286 and the smallest angle 26.2 to 25.7; in the cost (weights 0.001 to 0.1) the smallest angle 24.4 to 19.0.
16. **Adaptive subdivision levels: landed** (`adaptive_angle`, Loop and the butterfly, test 172, figure 54 key S): the faces at a turn and their edges split, one or two split edges cut into two or three triangles; closed, the cube's middles coarse. Next: Catmull-Clark (quads need transition patterns), curvature of the limit rather than the faces' turn. Catmull-Clark done (2026-10-22): the transition a fan of triangles from the face point (first tried: the split edges' points kept in the neighbour as a polygon; its fan folded where they lined up, dark streaks in the figure). The cube of 4 x 4 quads: three levels 4,680 faces against 6,144.
17. **Adaptive offsets: landed** (CONTOUR, `evaluated` in the report, test 171): the distance evaluated only in blocks the level may cross, blocks halved from a coarse stride (the distance 1-Lipschitz: a corner farther than the diagonal clears its block); the same surfaces as before, 54%, 38%, 22% of the grid at 40, 64, 128 cells. Next: the grid's memory (still dense: an octree's storage for very fine offsets). Done (2026-10-26), sparse rather than an octree: the dual contouring's cells by the crossings (sorted pairs, binary search; bitwise the same surfaces), the grid's positions made when needed: 10 bytes a point from 142, a 64-cell offset 51.5 MB to 25 MB. Left: the values (8 bytes a point) in a hash of the evaluated points, worth it past 256 cells only.

Parameterization:

18. **BFF: landed** (`QAWS_PARAM_BFF`, `bff_disk`, test 166, figure 55 key T): the Yamabe scale factors, the border's curve from its new turning and lengths closed by least change, the inside harmonic; lengths kept to 1e-4 on a cap (angles below LSCM's), a disc round to 0.1%, a sphere cut pole to pole flattened. Next: prescribed border angles in general (cone points), the inside by the conjugate harmonic rather than both coordinates' harmonic extension.
19. **Seams and chart packing: landed** (`auto_seams`, `pack`, `fill`, test 166, figure 56 key U): charts grouped by the axis direction nearest each normal, packed on shelves into the unit square at one scale; groups smoothed by iterated conditional modes; charts split across their longest extent while their area scales differ by more than `max_stretch`; charts turned to their least box before the shelves, the rows' width searched. A closed icosphere: 6 charts, fill 0.52 (0.37 before turning), no overlap; `max_stretch` 1.2: 12 charts, area distortion 0.090 to 0.040. Next: splits along the distortion's level lines rather than a median plane, a skyline packing, the groups by graph cuts (shared with polycube labels).

Hexahedra and quads (research-scale):

20. **Polycube labels by graph cuts: landed** (`sm_label_cut`, alpha-expansion over a Dinic max-flow, test 176, figure 57 key V): data the facets' areas times 1 − their normals along the axes, smoothness the shared edges' lengths (twice between opposite sides). The L and the bent bar unchanged (valid, 0.98 and 0.26); a box with a bump of 0.3 / 0.4 at 0.68 / 0.67 (0.36 / 0.29 before); over bumps 0.25 to 0.65 the mean 0.48 (0.37). Next: spikes (bumps past half the height) given their own patches, the atlas groups by the same cut.
21. **GLOBAL torus periods: landed** (`qg_periods`, test 136, figure 58 key W): tree-cotree cuts, a whole step along each changing one period; the map's energy over the four steps fitted from 15 solves on one factor, its least integer point kept. Spirals in 14 of 20 tori before, 2 after (radii 1 and 0.3 to 0.7, at 0.12 to 0.4), at the same time; those left follow a field that turns slowly round the tube (INSTANT keeps the smoothest on a tie; guided along the parallels 17 of 20 close). Next: higher genus (4 g steps), the curvature-aligned field preferred on closed tubes.
22. **FRAME_FIELD through singular edges: landed** (`hs_hexes`, `ff_frames` with MBO, test 128, figure 59 key X): the field relaxed by MBO (20 rounds of 48 diffusion steps on cells at 0.28 of the cell size: a can's singular loop near its axis becomes four straight lines along its side, a ball's graph a cube's), a map with transitions across a reduced cut (CubeCover) rounded greedily (singular lines, border planes, transitions), its lattice read off as HexEx does (registrations merged by the transitions, cubes from their centres and, where flat cells hid them, from their corners), the hexahedra relaxed. A can at 0.4: 16 hexahedra, a ball at 0.5: 27, all valid. Tried first without MBO (kept in the plan's history): the sweeps' field collapses the map. Next: border pits (cubes whose corner cells squeezed flat), robustness across sizes (the rounding's choices), genus beyond zero, speed (the field's MBO and one factor per rounding round).

### Roadmap (2026-10-22): follow-ups

The enrichments' "Next" items, in the order agreed: the larger ones first, then the smaller ones.

Larger:

1. **FRAME_FIELD through singular edges, robust:** no border pits, results stable across cell sizes (rounding that turns
   no cell over), genus beyond zero (a torus), speed (MBO's projections, one factor per rounding round). Check: a can,
   two balls and a torus over five cell sizes, all valid and without pits.
   Progress (2026-10-22): flip-aware rounding, the relaxation holding only creases and stepping up the gradient, a
   handle's disc kept in the cut: a can and a 32 x 16 ball valid at all five sizes (scaled Jacobian 0.26 and up), a
   torus (radii 2 and 0.7) at 0.6 valid (test 128). Left: the 16 x 8 ball at 0.3 and 0.35 (folded before rounding),
   thinner and coarser tori (a zero map: singular curves round the hole), dents, speed (a 0.3 ball about a minute).
2. **Torus periods, further:** higher genus (4 g steps, the same fit); the curvature-aligned field preferred on closed
   tubes. Check: a double torus with no spiral, the thin torus without guides.
   Done (2026-10-22): the aligned field kept on a tie for GLOBAL on closed surfaces with a handle (the thin torus clean
   at all five sizes, test 136; 3 of 20 tori still spiral, the fat one); the handles' cuts shared (`qg_generators`) and
   made unknowns of `qg_quantize` (any genus, with cones). Left: the double torus gets no map at all (a vertex beside a
   cone wound twice, a cone placement fault); a preference for an aligned period lattice when its energy is nearly as low
   (the fat torus).
3. **Polycube labels, further:** spikes (bumps past half the height) given their own patches; the atlas groups by the
   same cut. Check: the bump at 0.6 valid; atlas borders from the cut.
   Done (2026-10-22): the cut gives a spike its own wall patches, too small for hexahedra (two cell faces each), so the
   patches under 3 h^2 are merged into the label they share most border with: the 0.6 bump at 0.62 (0.025 before),
   the bent bar and the L unchanged (test 176); the atlas groups by the same cut (sm_label_cut), the blob's fill 0.52
   to 0.57, borders cleaner (test 166, figure 56).
4. **Atlas, further:** splits along the distortion's level lines, a skyline packing, groups by the graph cut. Check:
   stretch bound met with fewer charts, fill above 0.6.
   Done (2026-10-22): packing by the charts' shapes (bitmap, as xatlas; a skyline of boxes as the fallback): the blob
   0.57 to 0.60, its split atlas 0.47 to 0.51, the icosphere's 12 charts 0.46 to 0.49; its 6 equal square charts stay
   at 0.51 (three by two in a square: about 0.53 at best). Groups by the graph cut: item 3. Tried and dropped: splits
   along the scale's level lines (a two-label cut on the log scale): on a cap they cut a disc out of a ring, 19 charts
   for 12 and more distortion; the median plane stays.
5. **BFF, further:** prescribed border angles (cone points), the inside by the conjugate harmonic. Check: a square
   chart from a disc with four corners, angles kept.
   Done (2026-10-22): `bff_turning` per vertex, the rest of 2 pi spread by length (the disc's rounds): the cap to a
   square, corners and sides exact to 1e-14, none turned (an ear with three corners on a side goes flat). Tried and
   dropped: v the conjugate harmonic of u (least E_D - A, u held): the cap's angles 0.0447 for 0.0445 and the border's
   lengths off by 0.4%; the harmonic extension stays. Interior cones (cut to them) not done.
6. **Weights per edge on every skeleton** (asked 2026-10-22). Today `qaws_skeleton` weighs its lines only
   (`weight_fn`); arcs and the circles opened at reflex vertices move at speed 1, so the medial axis mode (splice angle
   pi) is unweighted; the surface regions' `qaws_surface_region_medial_axis` and `qaws_surface_region_mitered_skeleton`
   take no weights. To do: a weight for arcs (the space-time cone |p - m| = r -+ w t, its events still linear or one
   quadratic), for the splice circles (by default the larger of the reflex vertex's two edges' weights, or a weight per
   vertex given), so the weighted medial axis (each feature's distance divided by its weight) comes out of the same
   events; on surfaces, a weight per boundary piece and per corner in the labels (the least weighted geodesic distance)
   and in the crossings and junctions (Newton on the weighted equalities). Check: a rectangle with one side at weight 2
   (the roof's known faces); a square's weighted medial axis against the weighted bisectors; a planar surface region
   equal to the 2D skeleton with the same weights.
   Done (2026-10-22): circles carry a speed (|p - m| = r0 - s w t); three supports with cones of different speeds meet
   by a quartic (x, y quadratic in t; roots by bisection between the derivative's), the parallel-rows case by a
   quadratic in t; splice circles at the faster edge's speed, a slower neighbour a conic corner; a smooth joint whose
   edges part opens a circle at once; corners born touching between different speeds take the crossing where the
   edge coming in enters the other's swept side. Fixed on the way: the united region's edges named the wrong input
   curve (polylines joined, a vertex naming the curve ending there), now the candidate nearest the edge's middle.
   Checks: the rectangle (top 2/3), the disc at 2 (1/2), the stadium with arcs at 2 and 0.5 (1, symmetric), the L's
   and the bump's weighted medial axes against the weighted distance (vertices at the larger weight). Surfaces:
   `weight_fn` in `qaws_surface_medial_desc` (region path, curve, knot span), corners at the larger weight; a flat
   rectangle with its bottom at 2: radii the weighted distances to 1e-7, the largest 1.6 / 3, both the medial axis and
   the mitered skeleton. Not done: a weight per vertex given (the default rule only).

Smaller (each with its check): Zorin's refined butterfly border rules; Kobbelt's ternary sqrt(3) border; Catmull-Clark
adaptive levels; geodesic regions' borders inside triangles; fewer saddle windows; corner attributes in decimation (UV
wedges); Lindstrom and Turk's boundary and shape terms; an octree for fine offsets; rounded offset pieces cheaper; the
curved-element optimizer by the bound's gradient with border nodes sliding; simplex subdivision for the tetrahedra's
bound; point clouds with cavities.

All done or settled (2026-10-26): butterfly border (e761aae), sqrt(3) ternary (210d10e), Catmull-Clark adaptive (03018a1),
geodesic borders (08ed5c6), saddle windows tried and dropped (6af4afc), decimation wedges (4410d75), Lindstrom and Turk's
border (3da3f8d, shape dropped), sparse contouring memory (e28dd20), rounded offsets cheaper (d39dda9), untangling with
border sliding (18a0b56, gradient dropped), the tetrahedra's bound bisected (17a1f9a), point clouds with cavities
(a80a4e0). Next: the sweeps and continuous Booleans roadmap below.

### Roadmap (2026-10-22): sweeps and continuous Booleans

After the follow-ups and building on the Booleans (items 5 and 6): animation as a first-class input, so a mesh moving
between poses carves, adds or meets another one continuously, not at sampled instants. One instance, the continuous
difference (carving), generalized to every Boolean: sweep(A, motion) op B, op in union, intersection, difference,
and both moving (their relative motion).

Use cases: CNC simulation (the stock less the tool's sweep along its path, the removed volume measured), continuous
collision detection (the sweep's intersection with an obstacle: whether, when first, how much volume), virtual sculpting
by sweeps alone (a continuous union along the stroke).

Inspiration and prior work: Wang, Zhang et al., Implicit Swept Volume SDF (SIGGRAPH 2024, TOG; video
https://www.youtube.com/watch?v=XbcX-jPE89U): the swept volume's signed distance as a semi-infinite program (the least
over time of the moving shape's distance), queried implicitly, no tunnelling between samples; Sellan, Aigerman,
Jacobson, Swept Volumes via Spacetime Numerical Continuation (2021): the sweep as a 4D object's projection, its
boundary traced in space-time; the lifted surfacing of generalized sweep volumes (Adobe research).

1. **Motion as input.** A mesh and its motion: rigid (rotation and translation), affine, or per vertex (a deforming
   mesh). The path either given as a qaws curve or spline (qaws's own curves: exact, rational where they are), or built
   from two poses with their first, and optionally second, derivatives, the way PN triangles build a patch from its
   corners' positions and normals: cubic Hermite from positions and velocities, quintic from accelerations too, each
   vertex then on a polynomial (Bezier) trajectory, so the interpolation is exact and the swept triangle is a polynomial
   patch (barycentric in space, Bezier in time). Rotations as quaternions with rational Hermite (Cayley's rational
   rotations) so a rigid motion stays rational and exact. API sketch: `qaws_motion` (from a curve, from two poses and
   derivatives, evaluated with its derivatives at t), `qaws_mesh_3d_sweep` (the swept solid), `qaws_mesh_sweep_boolean`
   (out = sweep(A, motion) op B). Check: a pose evaluated at 0 and 1 exactly the given ones, derivatives matched, a
   Hermite path reproducing a cubic exactly.
   Done (2026-10-26): `qaws_motion.h` (docs/api/motion.md, test 179): two poses (Hermite translations, cubic or
   quintic; R0 Cay(c(t)) with c a Hermite polynomial in the Gibbs vector, rational and orthonormal), a 3D curve with a
   constant rotation, or per vertex; ends to 3e-16, cubics exact, the angular velocity against R's differences 5e-10.
   Not done: affine motions; a curve motion with a rotation along it (frames).
2. **The swept solid, conservative and certified.** Time cut adaptively until each slab's motion is within the
   tolerance of a straight one (bounded by the trajectory's derivatives); each triangle's slab sweep enclosed (its start,
   its end, its edges' ruled surfaces, the bulge bound); the union of the slabs' sweeps by the Booleans. The envelope's
   characteristic curves (where the velocity runs along the surface) on each slab mark the silhouette that bounds it.
   Check: a cube translated along a line is a box exactly; a box turned about an axis gives the annular sector's
   volume; a ball along a Hermite curve gives a tube (pi r^2 L + 4/3 pi r^3) within the tolerance.
   Done in part (2026-10-26): `qaws_mesh_3d_sweep` (test 180, figure 60 key Y): slabs halved until the chords are within
   tolerance (measured at the middles and quarters, not certified from the derivatives), each slab's front faces
   extruded between layers of vertices shared by consecutive slabs, the solid at 0 and all extrusions through the
   Boolean's self-union. A cube translated is exact (volume 3 to 1e-15); a box turned a quarter 2.37525 against 2.37627;
   a ball's tube within 1.7%. Not done: the certified bound (the trajectories' derivatives), the twist of rotating walls
   (bilinear quads split), the envelope's characteristic curves; the Boolean sometimes leaves a planar motion's
   coplanar walls open.
3. **Exact where the motion allows.** Translations along polylines and polynomial paths with the exact kernel:
   the sweep of a polyhedron along a segment is exact (its faces' prisms and the union), so G-code's linear moves give
   exact carving; polynomial and rational motions exact up to the slabs' certified bound, the Booleans exact on them.
   Check: a ball-end or flat-end tool along G-code lines in a block: the slots' volume exact (flat end) or within the
   bound (ball end).
   Not taken here (2026-10-26): exact Booleans are farq's (decision 9). The numerical sweep is already exact in
   geometry for translations (a cube along a segment is its Minkowski sum, volume to 1e-15; a flat-end tool's slots to
   1e-15 in test 180); exact arithmetic on them belongs with farq.
4. **Continuous Booleans.** `qaws_mesh_sweep_boolean`: difference (carving: the stock less the tool's sweep, the
   removed volume reported), union (sculpting: the stroke added), intersection (contact: the overlap's volume). Several
   motions in sequence (a tool path) carved in one pass, the stock's mesh updated incrementally. Check: carving a block
   by a sequence equals carving by the union of the sweeps; volumes conserved (stock = result + removed).
   Done (2026-10-26): `qaws_mesh_sweep_boolean` (test 180, figure 60): each move swept, all the sweeps one operand
   (their winding numbers add) against the stock in one Boolean; a cube tool's slot removing 0.5 and an L of two moves
   0.375 to 1e-15, the intersection (a separate Boolean) the same 0.375, stock = result + removed. Not done: the stock
   updated incrementally (each call carves all its moves at once).
5. **Continuous collision detection.** The first time of contact between a moving mesh and an obstacle: the slabs
   tested by the BVH (item 1), the first slab that touches bisected in time to the tolerance (exactly for translations by
   the exact kernel); the overlap's volume over the whole motion. Check: a ball falling on a plane touches at the
   analytic time; a thin wall is not tunnelled through at any speed.
   Done (2026-10-27): `qaws_mesh_3d_first_contact` (test 183, figure 71 key &): the sweep's slabs, each triangle's swept
   prism (caps, split sides, three tetrahedra) against the obstacle's triangles by a BVH and exact triangle predicates,
   the first touching slab halved until a step moves no vertex past the tolerance; a ball on a slab at 0.5 to 1.5e-5, a
   cube at 10 per unit time stopped by a wall 0.001 thick, a rod turning into a wall, a vertex motion. Then two
   moving meshes and the exact kernel for straight motions (obstacle_motion, report.exact; two cubes at 0.7
   and boxes at 2.65 / 6 exactly). Then self-collision
   (qaws_mesh_3d_first_self_contact: pairs sharing no vertex; a strip pinching at its cusp at 0.5, two parts of one
   mesh at 0.8). Then exact times for straight vertex motions (cubics isolated by Bernstein halving). Not done: folds
   within one vertex's ring.
6. **The implicit route, for comparison and sculpting.** The swept volume's signed distance as the least over time of
   the moving shape's distance (ISV-SDF; the BVH's signed distance for the shape, a 1D minimization in t per query
   bracketed by the slabs), contoured by the offsets' contouring (item 17 of the enrichments, adaptive blocks): smooth
   sweeps for sculpting, a cross-check of 2 and 4. Check: the implicit and the explicit sweeps within the tolerance of
   each other on the tube and the turned box.
   Done (2026-10-27): `qaws_mesh_3d_sweep_implicit` (test 184, figure 72 key *): time samples no vertex moving past two
   cells apart, pruned per point by bounding spheres, a rigid motion's time refined by golden sections; contoured by
   wt_contour_field. Within a cell of the explicit sweep (the cube translated, the turned box), the cutter of figure 69
   without revolution (its revolved sweep past it by 0.045 under the drift 0.05). The plan's 1D minimization is the
   golden refinement; a tolerance-bounded result needs a finer grid (the cell is the error).

Ideas noted (not ordered):

- Dynamic local remeshing for elastoplastic simulation (video: https://www.youtube.com/watch?v=moGZi1kim-8, likely Wicke
  et al. 2010): local tetrahedral remeshing during a simulation, quality-driven edge and face operations; next to
  `qaws_volume_mesh_improve`.

### Roadmap (2026-10-27): the reserved methods

Back on the meshing roadmap after the sweeps: the methods still reserved behind `method` (UNSUPPORTED_OPERATION), the
cheapest and most reused first:

1. **References:** all done 2026-10-27 (test 122, figure 2): PHONG (Boubekeur and Alexa, a quadratic patch
   degree-elevated into the PN net), NAGATA (quadratic edge curves from the normals), G1 (Walton-Meek boundary, quartic
   Gregory interior, Chiyokura-Kimura cross tangents), SUBDIVISION (exact Loop limit: Stam's box splines, local
   subdivision near extraordinary vertices and creases, against qaws_mesh_3d_subdivide's limit to 1e-15) and MLS (APSS
   with value-gradient-Hessian jets; IMLS with normals was tried first and bulges outward on curved surfaces).
2. **Remeshing:** FIELD_ALIGNED done 2026-10-27 (test 124, figure 77): the quadrangulation's INSTANT pipeline with
   mf levels of symmetry 6 (mf_compat_orientation_n, mf_compat_index_n, a hexagonal lattice), qd_tri_extract reading
   triangles off the lattice steps (votes, each directed edge once, small holes fanned; line clusters from both sides
   merged, edges skipping a line cluster or crossing a crease dropped), 3 LOCAL passes to finish.
3. **Generation:** STRUCTURED done 2026-10-27 (mg_structured: grids per untrimmed face, sides by arc length, counts by
   classes across shared edges, per-direction doubling to the tolerances; one surface's lines placed by need since
   2026-10-28; test 126, figure 78); PARAMETRIC done 2026-10-27 (mg_parametric: sides
   by their own curvature, a fan start remeshed by LOCAL in (u, v) with the first fundamental form as a metric sampled
   on a 64 x 64 reference grid, boundary pinned; figure 78 second row); RESTRICTED_DELAUNAY done 2026-10-27 (mg_delaunay: Lawson
   flips, Chew/Ruppert circumcentre insertion snapped on the faces, encroached feature edges split on their curves; 30
   degrees, Delaunay off the features; test 126, figure 78). Generation complete.
4. **Tetrahedra:** STUFFING done 2026-10-27 (tt_stuffing: BCC lattice, winding-signed distance, secant cuts,
   Labelle-Shewchuk warping, Dompierre prisms; test 127, figure 79) and ENVELOPE (tt_envelope: clean input kept or
   watertight contour refined to half the envelope, LOCAL within half the envelope on a LINEAR reference, DELAUNAY,
   improve; test 127, figure 79).
5. **Hexahedra:** an octree: deferred 2026-10-28 (an adaptive all-hex octree needs 3-refinement transition templates, whose face
   template we could not derive or find with coordinates, or a dual hybrid octree; GRID gives uniform all-hex meanwhile).

### Attributes through Booleans and sweeps (2026-10-26)

Done: `qaws_mesh_3d_boolean` carries vertex, corner and face attributes (by name; per corner in the result, each
interpolated in its source face; NORMAL negated on turned pieces), `qaws_mesh_3d_sweep` and `qaws_mesh_sweep_boolean`
carry them (directions turned with the motion; walls' face normals their own). Rules NORMAL, SRGB and OKLAB (alpha
premultiplied), honoured by edits, remeshing, subdivision, decimation, the connector and references. Test 181,
figure 61. Then (2026-10-26): edge attributes (QAWS_MESH_EDGE, on both corners of an edge as features are; an edit keeps
them on an edge and its pieces, new edges zeros; Booleans and sweeps carry them), the TANGENT rule, feature edges through
Booleans (kept, and the cuts marked); figure 62. Then the VECTOR and QUATERNION rules, vertex and corner attributes
matched in Booleans; figure 63 (four motion schemes carving, the rules compared). Then the ANGLE rule, edge attributes
through subdivision (tracked edges, each piece its input edge's value), quaternions and angles mixed sign-free there;
figure 64. Open:
- edge attributes as semi-sharp crease weights for subdivision (through the connector: done);
- volume attributes: done (vertex, corner, cell and facet ones through cloning, the border, improvement,
  tetrahedralization, hexahedralization, boundary layers and the connector; test 182, figures 65 and 66);
- features along a sweep's walls: done for feature corners' paths; envelopes (an edge sliding along itself) not marked;
- footpoints through Booleans.

## 0. Scope

In:

- Meshes as first-class, editable objects: polygonal surface meshes
  (triangles, quads, mixed; in 2D or 3D) and volume meshes (tetrahedra,
  hexahedra, mixed), with attributes and features.
- Generation from qaws geometry: 2D regions, parametric surfaces (NURBS and
  every other `qaws_surface`), trimmed surfaces and shells (the B-rep of
  qaws), implicit fields; volumes bounded by a mesh or a shell.
- Remeshing of meshes: uniform (isotropic), adaptive (sizing fields from
  curvature, error, attributes, the caller), anisotropic (metric fields),
  field-aligned; simplification; refinement.
- Element structure: triangulate (polygons), quadrangulate (pairing,
  all-quad, field-aligned, patch layouts), tetrahedralize, hexahedralize
  (grid, polycube, frame field), hex-dominant.
- The geometry between the vertices: the continuous surface a mesh with
  normals implies (section 4), used wherever a new vertex is placed.
- Quality and error measurement.

Out:

- File input and output (OBJ, STL, PLY, STEP, ...): the caller's role.
  Meshes enter and leave as arrays.
- Rendering.
- Exact Booleans and exact intersection resolution: farq's (decision 9).
  Numerical ones are in (decision 10, roadmap items 5 and 6).

## 1. What qaws has today

| Entry point | Output | What it does |
|---|---|---|
| `qaws_triangulate_region_2d` | `qaws_mesh_2d` | Constrained Delaunay (Bowyer-Watson, boundary recovered by flips), Ruppert refinement (min angle up to 33 degrees, area, size function), constraint curves, exact orientation and in-circle tests (coordinates snapped to a 2^-48 lattice, 384-bit integers); boundary vertices keep their curve and parameter |
| `qaws_mesh_2d_remesh` | in place | Botsch-Kobbelt split / collapse / flip / smooth toward a length or a size function, boundary vertices sliding on their curves |
| `qaws_mesh_2d_quadrangulate` | in place | Triangles paired into convex quads; all-quad by splitting through midpoints and centroids |
| `qaws_surface_region_triangulate` | `qaws_mesh_2d` in (u, v) | The 2D mesher in the parameter domain, isotropic size from the metric's largest eigenvalue; seams cut |
| `qaws_shell_region_triangulate` | `qaws_shell_mesh` | Faces meshed one by one, welded in space, boundary edges split so shared edges conform |
| `qaws_surface_tessellate` | caller buffers | Adaptive quadtree on a curvature threshold and a maximum edge |
| `qaws_surface_create_subdiv` | surface + mesh | Loop and Catmull-Clark refined to a level, the limit sampled on a grid |
| `qaws_stack_mesh_3d` | `qaws_mesh_3d` | Stacked sections stitched |
| `qaws_repulsive_mesh` and its figures | arrays | Repulsive descent; the figures carry a remesher (split above 1.5 L, collapse below L/2 with the link condition, Delaunay flips, tangential smoothing with a fold guard), collision checks (vertex-face, edge-edge, conservative advancement), Barnes-Hut trees, obstacles |
| `qaws_shell_atlas` | atlas | Trimmed faces bounded by (u, v) loops, shared sides found in space, seams and poles: the B-rep in practice (`qaws_brep_shell` has topology and 3D curves but no p-curves) |
| exact kernel | predicates, rationals | Filtered predicates on raw doubles, rational arithmetic |

Six mesh types, each read-only, none with attributes or features, none
editable by the caller. The 2D mesher is solid; the surface meshers are
isotropic in (u, v) with no error control; there is no volume meshing and
no quad meshing beyond pairing and splitting.

## 2. Principles

1. **Reference-driven.** Every operation that creates or moves a vertex
   asks a *reference* where the surface is: an exact surface, a shell, an
   implicit field, or a mesh with an interpolant (section 4) — never the
   mesh being edited. Each output vertex keeps its *footpoint* on the
   reference (a patch and coordinates in it). This removes drift over
   iterations, transfers attributes exactly (evaluated at the footpoint,
   not averaged from neighbours), keeps features where the reference has
   them, and makes the error measurable.
2. **Arrays in, arrays out.** Opaque handles owned by the library; caller
   arrays copied unless a field says borrowed. No I/O.
3. **The qaws conventions.** Descriptors zero-initialized mean defaults;
   `qaws_status` for the call, a report struct for what happened (counts,
   quality reached, error, reasons), as the repulsive solvers do.
4. **Deterministic.** The same result for any executor: parallel loops over
   fixed partitions or graph colourings, results gathered in a fixed order.
5. **Robust combinatorics.** Orientation, in-circle, in-sphere and
   segment-triangle decisions through the exact kernel's filtered
   predicates; constructions in double; exact constructions (rationals,
   indirect predicates) where an algorithm needs them (3D constrained
   Delaunay, envelope tetrahedralization).
6. **`qaws_scalar` at the API, double inside.** The f32 build meshes in
   double.
7. **Features explicit.** Sharp edges, boundaries, corners and attribute
   seams, detected (dihedral angle, attribute discontinuity) or given; every
   operator respects them: no flip across, sliding along, corners fixed.
8. **Additive.** The existing mesh outputs get adapters to the new type;
   nothing is removed before a major version.

## 3. The state of the art

Grouped by what qaws would offer. Each table: the work, its idea, what
qaws takes from it.

### 3.1 Surfaces through a mesh with normals

| Work | Idea | For qaws |
|---|---|---|
| Vlachos, Peters, Boyd, Mitchell 2001, *Curved PN triangles* (I3D) | Cubic Bezier triangle from 3 positions and 3 normals; each edge's control points from that edge alone, so neighbours share the edge curve (no cracks); a separate quadratic normal field | The default interpolant |
| McDonald, Kilgard 2010, *Crack-free PN triangles using adjacent edge normals* (NVIDIA) | PN with split normals at creases made crack-free from the adjacent edges' normals | Creases with per-corner normals |
| Boubekeur, Alexa 2008, *Phong tessellation* (SIGGRAPH Asia) | Quadratic: the linear point projected onto the three vertex tangent planes, blended by barycentrics, shape factor alpha (3/4) | The cheap smooth option |
| Nagata 2005, *Simple local interpolation of surfaces using normal vectors* (CAGD) | Quadratic patch; a per-edge curvature coefficient from the two normals; exact on circular arcs; singular when a normal nearly contains the edge (a limiter is needed) | Meshes tessellated from CAD (cylinders, spheres, fillets come back nearly exact) |
| Walton, Meek 1996, *A triangular G1 patch from boundary curves* (CAD); Loop, Schaefer, Ni, Castano 2009, *Gregory patches for hardware tessellation* (SIGGRAPH Asia) | Boundary curves first, cross-boundary tangents made consistent along each edge, Gregory-type interior: G1 across edges | The smooth, G1 surface "leaving" the mesh |
| Loop 1987; Stam 1998, *Exact evaluation of Catmull-Clark subdivision surfaces* (SIGGRAPH) and the Loop counterpart | Approximating subdivision, C2 away from extraordinary vertices, closed-form evaluation at any parameter | A subdivision reference (does not pass through the vertices) |
| Zorin, Schroder, Sweldens 1996, *Interpolating subdivision for meshes with arbitrary topology* (SIGGRAPH) | Modified butterfly: interpolating, C1, defined by refinement only | Optional |
| Alexa et al. 2003, *Computing and rendering point set surfaces* (TVCG); Guennebaud, Gross 2007, *Algebraic point set surfaces* (SIGGRAPH); Oztireli, Guennebaud, Gross 2009, *Feature preserving point set surfaces* (EG) | Moving-least-squares projection from points and normals; spheres fitted locally; robust kernel keeps sharp features | Point clouds and noisy scans as references |
| Macedo, Gois, Velho 2011, *Hermite radial basis functions implicits* (CGF) | A global implicit interpolating points and normals | Small and medium inputs |
| Sharp, Soliman, Crane 2019, *Navigating intrinsic triangulations* (SIGGRAPH); Gillespie, Sharp, Crane 2021, *Integer coordinates for intrinsic geometry processing* (SIGGRAPH Asia) | Points kept exactly on the input's surface through (face, barycentric) records and straightest walks | The footpoint mechanics on a linear reference |
| Lee et al. 1998, *MAPS* (SIGGRAPH) | A bijection between a remeshed surface and the original, kept through simplification | Correspondence for attribute transfer |

### 3.2 Remeshing, uniform and adaptive

| Work | Idea | For qaws |
|---|---|---|
| Botsch, Kobbelt 2004, *A remeshing approach to multiresolution modeling* (SGP) | Split above 4/3 L, collapse below 4/5 L, flip toward valence 6, tangential smoothing, projection back | The main remesher (the figures already carry one) |
| Surazhsky, Gotsman 2003, *Explicit surface remeshing* (SGP) | Local operations with area equalization in local parameterizations | Smoothing that equalizes density |
| Alliez, Meyer, Desbrun 2002, *Interactive geometry remeshing* (SIGGRAPH); Alliez et al. 2003, *Anisotropic polygonal remeshing* (SIGGRAPH) | Remeshing through a parameterization; curvature-line-aligned anisotropic meshes | Background for the field-aligned modes |
| Dunyach, Vanderhaeghe, Barthe, Botsch 2013, *Adaptive remeshing for real-time mesh deformation* (EG short) | Botsch-Kobbelt with a curvature-driven target length per vertex | The adaptive remesher |
| Hu, Yan, Bommes, Alliez, Benes 2017, *Error-bounded and feature preserving surface remeshing with minimal angle improvement* (TVCG) | Every operation checked against a Hausdorff bound; the smallest angle improved | The error-bounded mode |
| Du, Faber, Gunzburger 1999, *Centroidal Voronoi tessellations* (SIAM Review); Liu et al. 2009, *On centroidal Voronoi tessellation: energy smoothness and fast computation* (TOG); Yan, Levy, Liu, Sun, Wang 2009, *Isotropic remeshing with fast and exact computation of restricted Voronoi diagram* (SGP) | CVT by Lloyd or L-BFGS over the restricted Voronoi diagram on the surface; the dual is the mesh | The CVT remesher (most uniform, slower) |
| Valette, Chassery, Prost 2008, *Generic remeshing ... discrete Voronoi diagrams* (TVCG) | Clustering of the input's faces (ACVD): fast, approximate CVT | A fast option for dense inputs |
| Levy, Liu 2010, *Lp centroidal Voronoi tessellation and its applications* (SIGGRAPH) | CVT with an Lp norm and a metric: anisotropic, quad-dominant | Anisotropic CVT |
| Garland, Heckbert 1997, *Quadric error metrics* (SIGGRAPH); Hoppe 1999, *New quadric metric for meshes with appearance attributes* (IEEE Vis); Lindstrom, Turk 1998, *Fast and memory efficient polygonal simplification* (IEEE Vis) | Edge collapses ranked by quadrics; attributes in the quadric; volume and boundary preservation | Simplification |
| Botsch, Kobbelt, Pauly, Alliez, Levy 2010, *Polygon Mesh Processing* (book) | The reference text for the above | |

Sizing and metrics:

| Work | Idea | For qaws |
|---|---|---|
| Meyer, Desbrun, Schroder, Barr 2003, *Discrete differential-geometry operators* ; Rusinkiewicz 2004, *Estimating curvatures and their derivatives on triangle meshes* (3DPVT); Cohen-Steiner, Morvan 2003, *Restricted Delaunay triangulations and normal cycle* (SoCG) | Curvature on meshes | Curvature where no smooth reference is given (otherwise read from the interpolant's second derivatives) |
| Persson 2006, *Mesh size functions for implicit geometries and PDE-based gradient limiting* (EwC); Borouchaki, Hecht, Frey 1998, *Mesh gradation control* (IJNME) | Size functions limited in gradient: h(x) <= h(y) + beta \|x - y\| | Gradation of every size field |
| Alauzet 2010, *Size gradation control of anisotropic meshes* (FEAD); Alauzet, Loseille 2016, *A decade of progress on anisotropic mesh adaptation for CFD* (CAD) | Metric fields, their intersection and gradation, unit-mesh adaptation | The metric fields |

### 3.3 Planar meshing

| Work | Idea | For qaws |
|---|---|---|
| Chew 1989, *Constrained Delaunay triangulations* (Algorithmica); Ruppert 1995 (J. Algorithms); Shewchuk 1996, *Triangle*; Shewchuk 2002, *Delaunay refinement algorithms for triangular mesh generation* (CGTA) | Constrained Delaunay and refinement with angle guarantees (Ruppert 20.7 degrees provable, about 33 in practice) | Already in `qaws_mesh_2d` |
| Ungor 2004, *Off-centers* (LATIN) | Steiner points closer to the short edge: fewer points for the same angle | Refinement upgrade |
| Bern, Eppstein, Gilbert 1994, *Provably good mesh generation* (JCSS) | Quadtree meshing with guarantees | Graded meshes |
| Chen, Xu 2004, *Optimal Delaunay triangulations* (JCM); Cheng, Dey, Shewchuk 2012, *Delaunay Mesh Generation* (book) | ODT smoothing; the reference text | Smoothing |
| Blacker, Stephenson 1991, *Paving* (IJNME); Owen, Staten, Canann, Saigal 1999, *Q-Morph* (IJNME) | Advancing-front quads (Q-Morph from a triangle mesh) | Boundary-aligned quads |
| Remacle et al. 2012, *Blossom-Quad* (IJNME) | Triangles paired by minimum-cost perfect matching: all-quad with good shape | Upgrade of the current greedy pairing |
| Kowalski, Ledoux, Frey 2013, *A PDE based approach to multi-domain partitioning and quadrilateral meshing* (IMR) | Boundary-aligned cross fields in 2D, domains partitioned along their streamlines | Structured 2D quads |
| Hu et al. 2019, *TriWild* (SIGGRAPH); Mandad, Campen 2020, *Bezier guarding* (SIGGRAPH); Khanteimouri, Mandad, Campen 2022, *Rational Bezier guarding* (SGP) | Curved (higher-order) triangles conforming exactly to curved boundaries, optionally with exact rational arithmetic | Curved elements from qaws curves (later) |

### 3.4 Meshing parametric surfaces and B-reps

| Work | Idea | For qaws |
|---|---|---|
| Tristano, Owen, Canann 1998, *Advancing front surface mesh generation in parametric space using a Riemannian surface definition* (IMR); Borouchaki, Laug, George 2000, *Parametric surface meshing using a combined advancing-front generalized-Delaunay approach* (IJNME) | Meshing in (u, v) with the first fundamental form as a metric | Anisotropic (u, v) meshing (today qaws takes the metric's largest eigenvalue: isotropic) |
| Chew 1993, *Guaranteed-quality mesh generation for curved surfaces* (SoCG); Boissonnat, Oudot 2005, *Provably good sampling and meshing of surfaces* (Graphical Models); Cheng, Dey, Ramos 2007, *Delaunay refinement for piecewise smooth complexes* (SODA) | Restricted Delaunay refinement on the surface: facet angle, size and distance criteria, guarantees on topology and Hausdorff distance | The guaranteed mode, for shells and implicits |
| Jamin, Alliez, Yvinec, Boissonnat 2015, *CGALmesh* (TOMS) | A generic Delaunay refinement framework over an oracle of the domain | The oracle idea: the reference |
| Geuzaine, Remacle 2009, *Gmsh* (IJNME) | Curves first, then faces in (u, v), then volumes; shared discretizations make the result watertight | The order of a B-rep mesher |
| Sederberg et al. 2008, *Watertight trimmed NURBS* (SIGGRAPH) | Trims removed by conversion to T-splines | Background (out of scope) |
| Reberol, Georgiadis, Remacle 2021, *Quasi-structured quadrilateral meshing in Gmsh* (arXiv) | Cross fields on CAD faces, patches, quad meshes with few irregular vertices | Quads on shells |
| Lorensen, Cline 1987, *Marching cubes* (SIGGRAPH); Ju, Losasso, Schaefer, Warren 2002, *Dual contouring of Hermite data* (SIGGRAPH) | Implicit field to mesh; dual contouring keeps sharp features | Implicit references |

CAD practice (OpenCascade's tessellator and others) tessellates for display
by two numbers: the linear deflection (chord distance) and the angular
deflection (normal turn), with the edges discretized once and shared. qaws
should offer the same pair, on top of size and quality.

### 3.5 Quad meshing

| Work | Idea | For qaws |
|---|---|---|
| Bommes et al. 2013, *Quad-mesh generation and processing: a survey* (CGF); Campen 2017, *Partitioning surfaces into quadrilateral patches: a survey* (CGF); Vaxman et al. 2016, *Directional field synthesis, design, and processing* (CGF) | Surveys | Reading |
| Ray, Li, Levy, Sheffer, Alliez 2006, *Periodic global parameterization* (TOG); Ray, Vallet, Li, Levy 2008, *N-symmetry direction field design* (TOG); Knoppel, Crane, Pinkall, Schroder 2013, *Globally optimal direction fields* (SIGGRAPH); Diamanti et al. 2014, *Designing N-PolyVector fields with complex polynomials* (SGP); Panozzo et al. 2014, *Frame fields* (SIGGRAPH) | Cross fields (smoothest by an eigenproblem), non-orthogonal and anisotropic frame fields | The direction fields |
| Kalberer, Nieser, Polthier 2007, *QuadCover* (EG); Bommes, Zimmer, Kobbelt 2009, *Mixed-integer quadrangulation* (SIGGRAPH); Myles, Pietroni, Zorin 2014, *Robust field-aligned global parametrization* (SIGGRAPH) | A seamless parameterization aligned to the field, integers rounded | The global quad mode |
| Bommes et al. 2013, *Integer-grid maps for reliable quad meshing* (SIGGRAPH); Ebke, Bommes, Chen, Kobbelt 2013, *QEx* (SIGGRAPH Asia); Campen, Bommes, Kobbelt 2015, *Quantized global parametrization* (SIGGRAPH Asia); Lyon, Campen, Kobbelt 2021, *Quad layouts via constrained T-mesh quantization* (CGF); Heistermann, Warnett, Bommes 2023, *Min-deviation-flow in bi-directed graphs for T-mesh quantization* (SIGGRAPH) | Guaranteed-valid integer-grid maps; robust extraction; quantization as a flow problem (Bi-MDF: exact, a fraction of a percent of QuadWild's time) | The reliable path |
| Jakob, Tarini, Panozzo, Sorkine-Hornung 2015, *Instant field-aligned meshes* (SIGGRAPH Asia) | Local smoothing of an orientation field and a position field on a multiresolution hierarchy, extraction by collapsing; quads or triangles; under a second for 10^5 faces | The fast quad and triangle mode (section 7) |
| Huang, Zhou, Niessner, Shewchuk, Guibas 2018, *QuadriFlow* (SGP) | Instant Meshes' objective with consistency constraints, solved by min-cost flow and local SAT: pure quads, about 4x fewer singularities, seconds | The fast mode's consistency step |
| Schertler et al. 2017, *Field-aligned online surface reconstruction* (SIGGRAPH) | Instant Meshes on a stream of points | Point clouds |
| Pietroni, Nuvoli, Alderighi, Cignoni, Tarini 2021, *Reliable feature-line driven quad-remeshing* (QuadWild, SIGGRAPH) | Field traced into patches along feature lines, sizes by integer programming, each patch quad-meshed: pure quads with clean layouts | The layout mode |
| Dong et al. 2006, *Spectral surface quadrangulation* (SIGGRAPH) | Quads from a Laplacian eigenfunction's Morse-Smale complex | Background |
| Tarini et al. 2010, *Practical quad mesh simplification* (EG); Bommes et al. 2011, *Global structure optimization of quadrilateral meshes* (EG); Campen, Bommes, Kobbelt 2012, *Dual loops meshing* (SIGGRAPH) | Quad simplification and layout optimization | Post-processing |
| NeurCross (arXiv 2024), CrossGen (arXiv 2025) | Learned cross fields | Watch; not for the library |

### 3.6 Tetrahedral meshing

| Work | Idea | For qaws |
|---|---|---|
| Shewchuk 1998, *Tetrahedral mesh generation by Delaunay refinement* (SoCG); Si 2015, *TetGen* (TOMS) | Constrained Delaunay tetrahedralization and refinement (radius-edge ratio) | The conforming mode |
| Diazzi, Panozzo, Vaxman, Attene 2023, *Constrained Delaunay tetrahedrization: a robust and practical approach* (SIGGRAPH) | CDT of any valid PLC through indirect predicates (Steiner points implicit): succeeds on all 4408 valid Thingi10k models, fixes a gap in TetGen's theory | The robust core, with qaws's exact kernel |
| Attene 2020, *Indirect predicates for geometric constructions* (CAD) | Exact predicates on constructed points without rationals | Exact kernel extension |
| Alliez, Cohen-Steiner, Yvinec, Desbrun 2005, *Variational tetrahedral meshing* (SIGGRAPH); Cheng et al. 2000, *Sliver exudation* (JACM); Klingner, Shewchuk 2007, *Aggressive tetrahedral mesh improvement* (IMR) | ODT smoothing; sliver removal by weights; topological and smoothing improvement | Optimization |
| Labelle, Shewchuk 2007, *Isosurface stuffing* (SIGGRAPH) | Tets from an implicit with dihedral bounds, via a BCC lattice | Implicit volumes |
| Hu et al. 2018, *Tetrahedral meshing in the wild* (TetWild, SIGGRAPH); Hu et al. 2020, *Fast tetrahedral meshing in the wild* (fTetWild, SIGGRAPH) | Any triangle soup: the boundary kept within an envelope, exact or hybrid arithmetic, never fails | The envelope mode |
| Schoberl 1997, *NETGEN* (CVS) | Advancing front with Delaunay | Background |
| Jiang et al. 2021, *Bijective and coarse high-order tetrahedral meshes* (TOG); Khanteimouri, Campen 2023, *3D Bezier guarding* (TOG) | Curved tetrahedra | Later |

### 3.7 Hexahedral meshing

| Work | Idea | For qaws |
|---|---|---|
| Pietroni et al. 2022, *Hex-mesh generation and processing: a survey* (TOG) | The survey (polycube, integer-grid maps, block decomposition, dual sheets, frame fields) | Reading |
| Schneiders 1996, *A grid-based algorithm for hexahedral meshing* (EwC); Zhang, Bajaj 2006, *Adaptive and quality quadrilateral/hexahedral meshing from volumetric data* (CMAME); Marechal 2009, *Advances in octree-based all-hexahedral mesh generation* (IMR); Gao, Shen, Panozzo 2019, *Feature preserving octree-based hexahedral meshing* (SGP) | Octree, boundary projected and padded | The robust mode (any closed input; poorer elements at the boundary) |
| Gregson, Sheffer, Zhang 2011, *All-hex mesh generation via volumetric polycube deformation* (SGP); Livesu et al. 2013, *PolyCut* (SIGGRAPH Asia); Fang et al. 2016, *All-hex meshing using closed-form induced polycube* (SIGGRAPH); Guo et al. 2020, *Cut-enhanced polycube-maps* (TOG) | Map to an axis-aligned polycube, mesh the polycube | Mechanical parts |
| Nieser, Reitebuch, Polthier 2011, *CubeCover* (SGP); Huang et al. 2011, *Boundary aligned smooth 3D cross-frame field* (SIGGRAPH Asia); Ray, Sokolov, Levy 2016, *Practical 3D frame field generation* (SIGGRAPH Asia); Solomon, Vaxman, Bommes 2017, *Boundary element octahedral fields in volumes* (TOG); Liu et al. 2018, *Singularity-constrained octahedral fields* (SIGGRAPH); Palmer, Bommes, Solomon 2020, *Algebraic representations for volumetric frame fields* (TOG) | Frame fields and volumetric integer-grid maps | The research path |
| Lyon, Bommes, Kobbelt 2016, *HexEx* (SIGGRAPH); Liu, Bommes 2023, *Locally meshable frame fields* (SIGGRAPH; valid maps from 2 % to 58 % on HexMe); Kohler, Heistermann, Bommes 2025 (SIGGRAPH) | Extraction; making fields meshable | The research path, still far from guaranteed |
| Sokolov et al. 2016, *Hexahedral-dominant meshing* (TOG); Gao, Jakob, Tarini, Panozzo 2017, *Robust hex-dominant mesh generation using field-guided polyhedral agglomeration* (SIGGRAPH); Livesu et al. 2020, *LoopyCuts* (SIGGRAPH) | Hex-dominant (Gao 2017 is Instant Meshes in volume) | Hex-dominant |
| Knupp 2001, *Algebraic mesh quality metrics* (SISC); Stimpson et al. 2007, *The Verdict geometric quality library* | Scaled Jacobian and the rest | Quality reports |

### 3.8 Mesh data structures

| Work | Idea | For qaws |
|---|---|---|
| Botsch et al. 2002, *OpenMesh*; Sieger, Botsch 2011, *Design, implementation, and evaluation of the Surface_mesh data structure* | Index-based half-edge, properties as arrays, garbage collection | The surface mesh |
| Kremer, Bommes, Kobbelt 2013, *OpenVolumeMesh* | Half-faces for polyhedral volumes | The volume mesh |
| Rossignac et al. 2001, corner tables; Damiand, Lienhardt 2014, *Combinatorial maps* (book) | Compact or dimension-independent alternatives | Considered, not proposed |

Libraries worth comparing against, offline only (no code enters qaws; as
for the repulsive work, algorithms are written from the papers):
pmp-library, CGAL, OpenMesh, geometry-central, libigl, Geogram, Gmsh,
Netgen, TetGen, fTetWild, Instant Meshes, QuadriFlow, QuadWild. Several
are GPL or AGPL (TetGen is AGPL, Triangle forbids commercial use without
a licence, CGAL's meshers are GPL): they must not leak into the runtime.

## 4. The surface between the vertices

A triangle mesh with vertex normals implies a smooth surface; the mesh is
a sampling of it. Remeshing should place new vertices on that surface, not
on the flat triangles: a coarse remesh of a finely curved part otherwise
loses its curvature at every split, and repeated remeshing shrinks it.

### 4.1 Construction: edge curves first, then the interior

Every scheme worth offering fits one pattern:

1. **An edge curve per edge**, from the edge's two endpoints and their
   data only (positions, normals at that edge's corners, the feature
   tangent if the edge is sharp). Two triangles sharing an edge share its
   curve: the surface has no cracks (C0), whatever the interior does.
2. **An interior per triangle** filling the three edge curves.

| Scheme | Edge curve | Interior | Through vertices | Vertex normals honoured | Across edges | Cost | Fails when |
|---|---|---|---|---|---|---|---|
| `LINEAR` | the segment | flat | yes | no | C0 | none | never (faceted) |
| `PHONG` | quadratic | quadratic | yes | approximately | C0 | very low | normals far from the face: bulges |
| `PN` | cubic Bezier, end tangents the chord projected on each tangent plane | cubic, centre point from the edge points | yes | yes (G1 at vertices) | C0 | low | normals inconsistent with the faces: ripples |
| `PN` + creases | as PN, the crease tangent at sharp corners | as PN | yes | yes, per corner | C0 at creases, no cracks with split normals | low | |
| `NAGATA` | quadratic, curvature from the two normals | quadratic | yes | yes | C0 | very low | a normal nearly containing its edge (limiter needed); exact on circular arcs |
| `G1` (Walton-Meek / Gregory) | cubic with consistent cross-boundary tangents | Gregory blend | yes | yes | G1 | moderate | inconsistent normals at a vertex (falls back to C0 there) |
| `SUBDIVISION` (Loop limit) | none (approximating) | box spline, exact evaluation | no (shrinks) | from the limit | C2, C1 at extraordinary vertices | moderate | a fitted control mesh is needed to pass through vertices |
| `MLS` (APSS / RIMLS) | none (projection) | projection | approximately | approximately | smooth | high | kernel radius vs density; for points and noise |

PN's edge curve is the cubic Hermite curve whose end tangents are the
chord projected onto each endpoint's tangent plane. That is the same rule
a crease needs with the crease direction in place of the projection, and
the same rule a G1 scheme builds its cross tangents on. So one edge-curve
module serves PN, creases and G1; NAGATA and PHONG are separate quadratics.

Normals are interpolated separately from positions (PN's quadratic normal
field, or the patch's own geometric normal): the caller chooses which one
a new vertex carries.

### 4.2 Footpoints, evaluation, projection, walking

A *location* on a mesh reference is `(triangle, b1, b2)` (barycentric,
b0 = 1 - b1 - b2); on a surface, `(face, u, v)`; on a shell, a
`qaws_shell_point`.

- **Evaluate** a location: position, normal, first derivatives, and second
  derivatives (the curvature tensor: sizing fields read curvature from the
  interpolant, not from a discrete estimator).
- **Project** a point: the closest location. Closest point on the flat
  triangles first (a BVH), then Newton on the patch in barycentric
  coordinates, hopping to the neighbour across an edge when a coordinate
  goes negative (the edge curve is shared, so the hop is continuous).
- **Walk** a location by a tangent vector: the straightest path across
  triangles (Sharp et al. 2019), carried onto the interpolant. Tangential
  smoothing moves vertices by walking, so they never leave the reference
  and no projection is needed (as the vines move in (u, v) on the vase).
- **Drift.** The reference is the input, fixed for the whole operation;
  remeshing a remeshed mesh starts from that mesh's footpoints, never from
  a new reference built on the remeshed mesh.

### 4.3 Defaults proposed

- Mesh references: `PN` with crease handling. `NAGATA` recommended for
  meshes tessellated from CAD (it recovers circular arcs); `LINEAR` when
  the mesh is the truth (scans, engineering meshes that must not move).
- Normals when none are given: angle-weighted face normals (Thurmer and
  Wuthrich), split at feature edges.
- Surface and shell references: exact, no interpolant.

## 5. Object model and API sketch

Names are proposals. Signatures follow the qaws conventions (descriptors,
`qaws_status`, out-parameters last).

### 5.1 The surface mesh

```c
/* qaws_mesh.h: a polygonal mesh, 2D or 3D, editable, with attributes and features */
typedef struct qaws_mesh qaws_mesh;

typedef struct qaws_mesh_desc
{
	unsigned int dimension;              /* 2 or 3 */
	qaws_scalar const* positions;        /* dimension per vertex */
	unsigned int vertex_count;
	unsigned int const* face_vertices;   /* all faces' vertices, counter-clockwise, concatenated */
	unsigned int const* face_sizes;      /* NULL: all triangles */
	unsigned int face_count;
	unsigned int flags;                  /* QAWS_MESH_WELD (merge vertices within weld_tolerance),
	                                        QAWS_MESH_SPLIT_NON_MANIFOLD (a fan per manifold neighbourhood) */
	qaws_scalar weld_tolerance;          /* 0: exact coincidence */
} qaws_mesh_desc;

qaws_status qaws_mesh_create(qaws_mesh_desc const* desc, qaws_mesh** out_mesh);
qaws_status qaws_mesh_clone(qaws_mesh const* mesh, qaws_mesh** out_mesh);
void qaws_mesh_destroy(qaws_mesh* mesh);

/* bulk reads into caller arrays */
unsigned int qaws_mesh_get_vertex_count(qaws_mesh const* mesh);
unsigned int qaws_mesh_get_face_count(qaws_mesh const* mesh);
unsigned int qaws_mesh_get_face_vertex_total(qaws_mesh const* mesh);
qaws_status qaws_mesh_get_positions(qaws_mesh const* mesh, qaws_scalar* out_positions);
qaws_status qaws_mesh_get_faces(qaws_mesh const* mesh, unsigned int* out_face_vertices, unsigned int* out_face_sizes);

/* topology */
qaws_status qaws_mesh_get_vertex_ring(qaws_mesh const* mesh, unsigned int vertex, unsigned int* out_vertices,
	unsigned int capacity, unsigned int* out_count);
qaws_status qaws_mesh_get_boundary_loops(...);
qaws_status qaws_mesh_get_topology(qaws_mesh const* mesh, qaws_mesh_topology* out);   /* components, boundaries, genus, manifoldness */

/* adapters from the existing outputs */
qaws_status qaws_mesh_from_mesh_2d(struct qaws_mesh_2d const* mesh, qaws_mesh** out_mesh);
qaws_status qaws_mesh_from_shell_mesh(struct qaws_shell_mesh const* mesh, qaws_mesh** out_mesh);
```

Internally: an index-based half-edge structure (struct of arrays), as
OpenMesh and pmp. Editing operations (split, collapse with the link
condition, flip, vertex insertion) exist inside; whether to make them
public is question 2.

### 5.2 Attributes

```c
typedef enum qaws_mesh_element { QAWS_MESH_VERTEX, QAWS_MESH_CORNER, QAWS_MESH_EDGE, QAWS_MESH_FACE } qaws_mesh_element;

typedef enum qaws_attribute_rule
{
	QAWS_ATTRIBUTE_LINEAR = 0,           /* barycentric at the footpoint */
	QAWS_ATTRIBUTE_NORMALIZED,           /* linear, then normalized (normals, tangents) */
	QAWS_ATTRIBUTE_ANGLE,                /* periodic */
	QAWS_ATTRIBUTE_NEAREST,              /* categorical: material ids, labels */
	QAWS_ATTRIBUTE_CUSTOM                /* a callback */
} qaws_attribute_rule;

typedef struct qaws_attribute_desc
{
	char const* name;
	qaws_mesh_element element;
	unsigned int components;
	qaws_attribute_rule rule;
	int seams_are_features;              /* corner attributes: their discontinuities become feature edges */
} qaws_attribute_desc;

qaws_status qaws_mesh_add_attribute(qaws_mesh* mesh, qaws_attribute_desc const* desc, qaws_scalar const* values,
	unsigned int* out_attribute);
qaws_status qaws_mesh_get_attribute(qaws_mesh const* mesh, unsigned int attribute, qaws_scalar* out_values);
```

Corner attributes are needed for UVs and split normals. Every operator
carries every attribute: new vertices evaluate the input's attribute at
their footpoint (exact transfer), so remeshing does not blur colours or
UVs. Attribute-driven criteria: a size field from an attribute or its
gradient; feature edges from its discontinuities.

### 5.3 Features

```c
qaws_status qaws_mesh_set_feature_edges(qaws_mesh* mesh, unsigned int const* edge_vertices, unsigned int edge_count);
qaws_status qaws_mesh_detect_features(qaws_mesh* mesh, qaws_scalar dihedral_angle, unsigned int flags);
	/* flags: QAWS_FEATURES_BOUNDARY, QAWS_FEATURES_ATTRIBUTE_SEAMS, QAWS_FEATURES_CORNERS_BY_ANGLE */
```

The feature graph: feature edges, and corners where more than two meet or
the curve turns sharply. Operators: feature vertices slide along their
feature curve (the edge curves of section 4.1), corners stay.

### 5.4 References

```c
typedef struct qaws_mesh_reference qaws_mesh_reference;

typedef enum qaws_interpolant
{
	QAWS_INTERPOLANT_PN = 0,
	QAWS_INTERPOLANT_LINEAR,
	QAWS_INTERPOLANT_PHONG,
	QAWS_INTERPOLANT_NAGATA,
	QAWS_INTERPOLANT_G1,
	QAWS_INTERPOLANT_SUBDIVISION,
	QAWS_INTERPOLANT_MLS
} qaws_interpolant;

typedef struct qaws_mesh_reference_desc
{
	qaws_interpolant interpolant;
	qaws_scalar const* normals;          /* 3 per vertex, or per corner with normals_per_corner; NULL: computed */
	int normals_per_corner;
	qaws_scalar phong_shape;             /* 0: 3/4 */
	qaws_scalar mls_radius;              /* 0: from the density */
} qaws_mesh_reference_desc;

qaws_status qaws_mesh_reference_create(qaws_mesh_reference** out, qaws_mesh const* mesh, qaws_mesh_reference_desc const* desc);
qaws_status qaws_mesh_reference_create_surface(qaws_surface const* surface, qaws_mesh_reference** out);
qaws_status qaws_mesh_reference_create_shell(qaws_shell_atlas const* shell, qaws_mesh_reference** out);
qaws_status qaws_mesh_reference_create_implicit(qaws_implicit_desc const* desc, qaws_mesh_reference** out);   /* value and gradient callbacks, bounds */
void qaws_mesh_reference_destroy(qaws_mesh_reference* reference);

typedef struct qaws_mesh_location { unsigned int patch; qaws_scalar a, b; } qaws_mesh_location;   /* triangle and barycentrics, or face and (u, v) */

qaws_status qaws_mesh_reference_evaluate(qaws_scalar* out_positions, qaws_scalar* out_normals, qaws_scalar* out_curvatures /* k1 k2 d1 d2, may be NULL */,
	qaws_mesh_reference const* reference, qaws_mesh_location const* locations, unsigned int count);
qaws_status qaws_mesh_reference_project(qaws_mesh_location* out_locations, qaws_scalar* out_distances, qaws_mesh_reference const* reference,
	qaws_scalar const* points, qaws_mesh_location const* hints /* may be NULL */, unsigned int count);
qaws_status qaws_mesh_reference_walk(qaws_mesh_location* out_to, qaws_mesh_reference const* reference, qaws_scalar const* tangent,
	qaws_mesh_location from);
```

Every mesh produced by an operator keeps its vertices' locations
(`qaws_mesh_get_footpoints`): the correspondence with the input.

### 5.5 Fields

```c
/* sizes: the target edge length at a point */
typedef struct qaws_size_field qaws_size_field;
qaws_status qaws_size_field_create_uniform(qaws_size_field** out, qaws_scalar length);
qaws_status qaws_size_field_create_curvature(qaws_size_field** out, qaws_mesh_reference const* reference,
	qaws_curvature_size_desc const* desc);   /* chord tolerance, normal-turn tolerance, min and max length */
qaws_status qaws_size_field_create_attribute(qaws_mesh const* mesh, unsigned int attribute, qaws_scalar scale, qaws_size_field** out);
qaws_status qaws_size_field_create_callback(qaws_size_field** out, qaws_scalar (*size)(void* user, qaws_scalar const* point), void* user);
qaws_status qaws_size_field_combine(qaws_size_field const* a, qaws_size_field const* b, qaws_size_combine how, qaws_size_field** out);
qaws_status qaws_size_field_limit_gradation(qaws_size_field* field, qaws_mesh_reference const* reference, qaws_scalar beta);

/* metrics: anisotropic sizes, a symmetric tensor per point (2x2 in the tangent plane, 3x3 in volumes) */
typedef struct qaws_metric_field qaws_metric_field;
qaws_status qaws_metric_field_create_curvature(...);
qaws_status qaws_metric_field_intersect(...);
qaws_status qaws_metric_field_limit_gradation(...);

/* directions: N-symmetric fields on surfaces (2 lines, 4 crosses, 6 for triangles), octahedral frames in volumes */
typedef struct qaws_direction_field qaws_direction_field;
typedef struct qaws_direction_field_desc
{
	unsigned int symmetry;               /* 4: crosses (0 = 4) */
	unsigned int method;                 /* QAWS_DIRECTIONS_LOCAL (Instant Meshes smoothing), QAWS_DIRECTIONS_SMOOTHEST (Knoppel et al.) */
	int align_features, align_boundary;
	qaws_scalar curvature_alignment;     /* 0: none; weight toward principal directions where they are reliable */
} qaws_direction_field_desc;
qaws_status qaws_direction_field_create(qaws_mesh_reference const* reference, qaws_direction_field_desc const* desc, qaws_direction_field** out);
qaws_status qaws_direction_field_get_singularities(...);
```

Sizes from tolerances, for a surface of largest curvature k at the point:
a chord of length L deviates from an arc of radius R = 1/k by
R - sqrt(R^2 - L^2/4), so the chord tolerance e gives
L = 2 sqrt(e (2R - e)); the normal-turn tolerance t gives L = 2R sin(t/2).
The size is the smaller, clamped to [min, max], then limited in gradient.

### 5.6 Operators

```c
/* generation */
qaws_status qaws_mesh_generate_region_2d(qaws_triangulate_desc const* region, qaws_generate_desc const* desc, qaws_mesh** out);
qaws_status qaws_mesh_generate_surface(qaws_surface const* surface, qaws_generate_desc const* desc, qaws_mesh** out);
qaws_status qaws_mesh_generate_shell(qaws_shell_atlas const* shell, qaws_generate_desc const* desc, qaws_mesh** out);
qaws_status qaws_mesh_generate_implicit(qaws_implicit_desc const* implicit, qaws_generate_desc const* desc, qaws_mesh** out);

typedef struct qaws_generate_desc
{
	qaws_scalar chord_tolerance;         /* the linear deflection; 0: from the size */
	qaws_scalar angle_tolerance;         /* the angular deflection, radians; 0: none */
	qaws_scalar min_length, max_length;
	qaws_size_field const* size;         /* optional */
	qaws_metric_field const* metric;     /* optional, anisotropic */
	qaws_scalar min_angle;               /* degrees; quality */
	unsigned int elements;               /* QAWS_ELEMENTS_TRIANGLES, _QUADS (pure), _QUAD_DOMINANT */
	unsigned int method;                 /* QAWS_GENERATE_PARAMETRIC (in (u, v) with the metric), _RESTRICTED_DELAUNAY (guarantees), _STRUCTURED (iso-parametric grid where a face allows) */
	qaws_batch_executor const* executor;
} qaws_generate_desc;

/* remeshing, simplification, refinement */
qaws_status qaws_mesh_remesh(qaws_mesh** out, qaws_mesh_report* report, qaws_mesh const* mesh, qaws_remesh_desc const* desc);
qaws_status qaws_mesh_simplify(qaws_mesh const* mesh, qaws_simplify_desc const* desc, qaws_mesh** out, qaws_mesh_report* report);
qaws_status qaws_mesh_subdivide(qaws_mesh const* mesh, qaws_subdivide_desc const* desc, qaws_mesh** out);   /* Loop, Catmull-Clark, sqrt(3), midpoint; new points on the reference */

typedef struct qaws_remesh_desc
{
	qaws_mesh_reference const* reference;   /* NULL: the input mesh, LINEAR */
	unsigned int method;                    /* QAWS_REMESH_LOCAL (Botsch-Kobbelt), QAWS_REMESH_CVT, QAWS_REMESH_FIELD_ALIGNED (Instant Meshes, triangles) */
	qaws_scalar target_length;              /* 0: the input's mean edge */
	qaws_size_field const* size;            /* adaptive */
	qaws_metric_field const* metric;        /* anisotropic */
	unsigned int iterations;                /* 0: 10 */
	qaws_scalar max_deviation;              /* > 0: every operation kept within it of the reference (error-bounded) */
	qaws_scalar min_angle;                  /* > 0: improve the smallest angle toward it */
	int keep_features, keep_boundary;       /* defaults 1 */
	int intersection_free;                  /* refuse operations creating self-intersections (the repulsive collision checks) */
	qaws_batch_executor const* executor;
} qaws_remesh_desc;

/* element structure */
qaws_status qaws_mesh_triangulate(qaws_mesh const* mesh, qaws_mesh** out);   /* polygons to triangles (constrained, on the reference) */
qaws_status qaws_mesh_quadrangulate(qaws_mesh** out, qaws_mesh_report* report, qaws_mesh const* mesh, qaws_quadrangulate_desc const* desc);
	/* method: QAWS_QUADS_PAIRING (Blossom matching), QAWS_QUADS_SPLIT (all quads, halving),
	   QAWS_QUADS_INSTANT (Instant Meshes + QuadriFlow consistency), QAWS_QUADS_GLOBAL (integer-grid map),
	   QAWS_QUADS_LAYOUT (feature-line patches, QuadWild-like) */

/* quality and smoothing */
qaws_status qaws_mesh_optimize(qaws_mesh* mesh, qaws_optimize_desc const* desc);   /* tangential smoothing on the reference, ODT, quality-driven flips */

/* measures */
qaws_status qaws_mesh_measure_quality(qaws_mesh const* mesh, qaws_mesh_quality* out);   /* angles, aspect ratios, valence, singularities, scaled Jacobian (quads) */
qaws_status qaws_mesh_measure_distance(qaws_mesh const* mesh, qaws_mesh_reference const* reference, qaws_distance_desc const* desc,
	qaws_mesh_distance* out);   /* one- or two-sided Hausdorff, mean, normal deviation */
```

### 5.7 Volume meshes

```c
/* qaws_volume_mesh.h: tetrahedra, hexahedra, prisms, pyramids (and polyhedra for hex-dominant results) */
typedef struct qaws_volume_mesh qaws_volume_mesh;
qaws_status qaws_volume_mesh_create(qaws_volume_mesh** out, qaws_volume_mesh_desc const* desc);
qaws_status qaws_volume_mesh_get_boundary(qaws_mesh** out, qaws_volume_mesh const* mesh);

qaws_status qaws_volume_mesh_tetrahedralize(qaws_tetrahedralize_desc const* desc, qaws_volume_mesh** out, qaws_mesh_report* report);
	/* boundary: a closed qaws_mesh, a shell or an implicit; mode: QAWS_TETS_CONFORMING (the boundary kept exactly: CDT),
	   QAWS_TETS_ENVELOPE (any soup, the boundary within an envelope); size field; radius-edge bound; minimal dihedral angle */
qaws_status qaws_volume_mesh_hexahedralize(qaws_hexahedralize_desc const* desc, qaws_volume_mesh** out, qaws_mesh_report* report);
	/* method: QAWS_HEXES_GRID (octree, boundary projected and padded), QAWS_HEXES_POLYCUBE, QAWS_HEXES_FRAME_FIELD;
	   allow_dominant: polyhedra where hexes do not fit */
qaws_status qaws_volume_mesh_measure_quality(...);
```

### 5.8 Reports

Each operator fills a report: what it did (counts of splits, collapses,
flips, iterations), what it reached (smallest angle, aspect ratio,
deviation to the reference, singularity count, valence histogram), and why
it stopped. A target not reached is not an error: the mesh is valid, the
report says how far it got.

## 6. Remeshing in detail

**Local operations** (the default). Per iteration:

1. Split edges longer than 4/3 of the target. The new vertex is the
   midpoint *on the reference*: the edge curve's midpoint (section 4.1),
   its footpoint stored.
2. Collapse edges shorter than 4/5 of the target, unless the link
   condition fails, a face would turn over, an edge longer than 4/3 would
   appear, a feature would be crossed, or (error-bounded) the deviation
   bound would be broken. The kept position: on the reference.
3. Flip edges toward valence 6 (4 on the boundary), or to Delaunay; never
   a feature edge.
4. Smooth tangentially: each vertex toward its neighbours' area-weighted
   centroid (or the CVT centroid), moved by walking on the reference.
   Feature vertices slide along their curve; corners stay.
5. (error-bounded) Measure the deviation at sample points of every changed
   face against the reference.

Adaptive: the length of an edge is its length over the target at its
midpoint (size field), or sqrt(e^T M e) (metric); the target is 1. The
repulsive figures' remesher (link condition, fold guard, Delaunay passes)
is the seed of this one, moved into the library and given footpoints,
features and attributes.

**CVT.** Seeds on the reference, the restricted Voronoi diagram computed
exactly (the exact predicates), Lloyd or L-BFGS on its energy (with the
density from the size field, or Lp with the metric); the dual is the mesh.
More uniform than local operations, slower, and the natural route to
anisotropic and blue-noise meshes.

**Field-aligned triangles.** Instant Meshes with 6-symmetric orientations
and a 3-symmetric position lattice (section 7): edges follow a field
(curvature lines, features), near-equilateral triangles.

**Intersection-free.** The repulsive work's collision checks (vertex-face,
edge-edge, conservative advancement) refuse any operation that would make
the mesh touch itself: thin parts survive coarse remeshing.

## 7. Instant Meshes, and what to build on it

**What it is** (Jakob, Tarini, Panozzo, Sorkine-Hornung 2015):

1. A hierarchy of the input (vertices clustered level by level).
2. **Orientation field**: a tangent direction per vertex, 4-symmetric
   (quads) or 6-symmetric (triangles), smoothed locally: each vertex takes
   the average of its neighbours' directions, each rotated to the closest
   of its symmetric copies, projected onto its tangent plane. Coarse to
   fine, Gauss-Seidel on a graph colouring (parallel).
3. **Position field**: a point per vertex on a local integer lattice
   aligned with the orientation, of spacing the target length, smoothed by
   the same scheme with lattice-compatible translations.
4. **Extraction**: vertices whose lattice points coincide collapse into one
   output vertex; neighbouring lattice points one step apart give edges;
   faces from the graph's cycles. The result is quad-dominant (some
   triangles and pentagons); pure quads by one subdivision step.

Fast, simple, scalable (10^5 faces under a second), any input (meshes,
point clouds).

**Its weaknesses**: singularities placed by local smoothing, many of them;
non-quad faces and occasional holes where the position field is
inconsistent; no explicit singularity control; extraction struggles on
thin parts; no correspondence with the input, so no attribute transfer.

**Successors**: QuadriFlow (2018) keeps its objective, makes the
integer offsets across edges consistent by a global min-cost flow and
local SAT: pure quads, about 4x fewer singularities, seconds. The global
integer-grid-map line (MIQ, IGM, QEx, quantization by Bi-MDF) gives the
best quality and guarantees at a higher cost; QuadWild (2021) and Bi-MDF
(2023) give clean patch layouts along feature lines. Gao et al. 2017 is
the same local scheme in volume (hex-dominant).

**What qaws would add**:

- Output vertices on the reference (the interpolant or the exact
  surface), with footpoints: smooth results from coarse inputs, exact
  attribute transfer.
- The sizing field as the lattice spacing (adaptive density) and metric
  fields through frames (anisotropy, Panozzo et al. 2014).
- qaws features as hard constraints of both fields (Instant Meshes aligns
  to creases and boundaries already); qaws seams (UV, material) too.
- Determinism: the colouring and the order fixed, same result for any
  executor.
- QuadriFlow's consistency step by default; Blossom matching to pair the
  remaining triangles; tangential relaxation on the reference after
  extraction.
- A report: singularity count, valence histogram, alignment deviation,
  distance to the reference.

The three tiers to offer: `QAWS_QUADS_INSTANT` (seconds, quad-dominant or
pure), `QAWS_QUADS_GLOBAL` (minutes, fewest and best-placed singularities,
guarantees through integer-grid maps), `QAWS_QUADS_LAYOUT` (clean patches
along features, for CAD and animation models).

## 8. Tessellating qaws geometry

For surfaces, trimmed surfaces and shells (the shell atlas is qaws's
B-rep):

1. **Curves first.** Every shared edge (and every trim curve) discretized
   once, by the chord and angle tolerances and the size field along it;
   both faces take the same samples, so the mesh is watertight by
   construction (today the shell mesher welds and splits afterwards).
2. **Faces in (u, v) with the metric.** The 2D mesher run in the
   parameter domain with the first fundamental form scaled by the target
   sizes as an anisotropic metric (today: isotropic, the largest
   eigenvalue), Delaunay refinement with the in-circle test in the
   metric. Seams and poles as today.
3. **Error control.** The chord tolerance checked at each face's centroid
   and edge midpoints against the exact surface; refinement until met
   (restricted Delaunay refinement, Boissonnat-Oudot, as the guaranteed
   alternative).
4. **Optional 3D pass.** The local remesher with the shell as reference
   (footpoints are (face, u, v), exact): isotropic in space whatever the
   parameterization.
5. **Quads.** An iso-parametric grid on untrimmed four-sided faces
   (structured, adaptive spacing by the tolerances); transfinite (Coons)
   maps on four-sided trimmed faces; otherwise the quad tiers of section 7
   with the shell as reference and its edges as features.

## 9. Planar meshes

The existing mesher becomes `qaws_mesh_generate_region_2d` (dimension 2),
the old entry points kept as adapters. Additions: off-center refinement;
anisotropic metrics; Blossom matching instead of greedy pairing;
boundary-aligned quads (Q-Morph, or cross fields in 2D); curved
higher-order triangles from qaws curves (Bezier guarding) as a later
option.

## 10. Volumes

**Tetrahedra.**

- Conforming: the constrained Delaunay tetrahedralization of the boundary
  (Diazzi et al. 2023, on qaws's exact predicates and indirect predicates),
  then Delaunay refinement (radius-edge ratio, size field), then
  improvement (flips, smoothing, sliver removal, Klingner-Shewchuk). Needs
  a closed, non-self-intersecting boundary.
- Envelope: for triangle soups (fTetWild): the boundary kept within an
  envelope, never fails.
- Implicit: isosurface stuffing.

**Hexahedra**, by robustness:

1. Grid / octree, boundary projected and padded (Marechal 2009; Gao et
   al. 2019): any closed input, guaranteed output, poorer elements near
   the boundary.
2. Polycube (Gregson et al. 2011; PolyCut; Fang et al. 2016): good for
   mechanical parts, fails on complex topology.
3. Frame field and integer-grid map (CubeCover, Ray et al. 2016, Liu et
   al. 2018, HexEx, Liu-Bommes 2023): the best meshes, not guaranteed; a
   research track.
4. Hex-dominant: field-guided agglomeration (Gao et al. 2017) or tets
   combined into hexes, prisms and pyramids.

Quality: scaled Jacobian and the Verdict metrics in the reports.

## 11. Testing and validation

- Interpolants against Mathematica fixtures (closed forms of PN, Phong,
  Nagata, their derivatives and curvatures), as the repulsive metric is.
- Invariants on every output: manifold, oriented, Euler characteristic,
  watertight (shells), no self-intersection (exact tests), footpoints
  exact (the position equals the reference evaluated at the footpoint).
- Quality and error distributions checked against bounds.
- Determinism across executors, f32 and f64 builds.
- Shapes generated procedurally in the tests; no files read.
- Offline comparisons against the reference implementations (pmp, Instant
  Meshes, QuadriFlow, Gmsh, TetGen or fTetWild) in a scratch area, as for
  repulsive; never linked.

## 12. Order of work (proposal)

1. `qaws_mesh`: container, attributes, features, adapters, quality and
   distance measures.
2. References: mesh with LINEAR, PN (+ creases), PHONG, NAGATA; surfaces
   and shells (exact); evaluate, project, walk; Mathematica fixtures.
3. Size and metric fields: uniform, curvature with tolerances, attribute,
   callback, combination, gradation.
4. Surface remeshing: local operations (isotropic, adaptive, anisotropic),
   features, attributes, error-bounded, intersection-free.
5. Simplification: quadrics with attributes, error bound.
6. Tessellation of surfaces and shells: curves first, anisotropic (u, v),
   tolerances; structured quads where a face allows.
7. Quads: Blossom pairing, all-quad split, direction fields, Instant
   (+ QuadriFlow consistency); then global integer-grid maps with
   quantization; then layouts.
8. Volumes: `qaws_volume_mesh`, conforming tetrahedralization and
   improvement; envelope mode; grid / octree hexes; polycubes; hex-dominant;
   frame fields last.
9. 2D unification and upgrades (off-centers, metric, Blossom, curved).
10. Later: G1 and subdivision interpolants, MLS references, CVT remeshing,
    curved elements.

## 13. Decisions to take

1. **Containers.** One surface mesh type (2D and 3D, any polygon) and one
   volume mesh type (proposed), or a single cell complex for all?
2. **Editing operations.** Public (split, collapse, flip, insert, with
   validity checks), so callers build their own remeshers; or internal
   only?
3. **Outputs.** Always a new mesh (input untouched, proposed), or in-place
   variants too?
4. **Footpoints.** Always kept and exposed (proposed; some memory), or on
   request?
5. **Default interpolant** for mesh references: PN (proposed), NAGATA,
   LINEAR?
6. **Attribute types.** `qaws_scalar` only (proposed), or integers and
   categorical types natively?
7. **Non-manifold and dirty input.** Weld and split non-manifold vertices
   (proposed); hole filling in, self-intersection resolution out (farq,
   decision 9).
8. **Exactness for volumes.** Exact predicates and indirect predicates
   (proposed), and the envelope mode with exact constructions (the
   rationals of the exact kernel)?
9. **Hex ambition.** Grid and polycube as the supported modes, frame
   fields as research; or frame fields as a goal?
10. **Implementation policy.** From the papers only, reference
    implementations for offline comparison (proposed, as repulsive).
11. **Branch.** A `meshing` branch, from `repulsive` (which has the
    remesher and collision checks to build on) or from main?
12. **Names.** `qaws_mesh` / `qaws_volume_mesh` / `qaws_mesh_reference`;
    and the fate of `qaws_mesh_2d`, `qaws_shell_mesh`, `qaws_mesh_3d`
    (adapters now, removal at a major version?).
13. **Curved elements** (higher-order triangles and tetrahedra from exact
    curves and surfaces): in scope, later?
