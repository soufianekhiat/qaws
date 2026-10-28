# Texture attributes

`qaws_mesh_texture.h`: an attribute held on a grid of texels over a mesh's (u, v), at a resolution of its own rather than
one value per vertex, corner or face. A mesh's attributes can only be as fine as its triangles; a texture keeps colours,
normals, times or labels sharp on a coarse or decimated mesh. The texels are the caller's array: qaws reads and writes no
files.

```c
qaws_scalar texels[256 * 256 * 3];
qaws_texture t = { texels, 256, 256, 3 };          /* width, height, components; rule set by baking */
qaws_texture_desc d = { 0 };                        /* "uv", no gutter, every face */
d.gutter = 4;
qaws_mesh_3d_bake_texture(coverage, &t, mesh, "color", &d);            /* a vertex, corner or face attribute */
qaws_mesh_3d_bake_texture_from(coverage, &t, coarse, detailed, "color", &d);   /* from another mesh */
qaws_texture_sample(values, &t, uv, count);                              /* any (u, v) */
qaws_mesh_3d_sample_texture(mesh, &t, "color", &d);                      /* onto the corners */
```

**Layout.** Texel (i, j) is at (j width + i) components, its centre at ((i + 0.5) / width, (j + 0.5) / height). Each axis
is clamped or repeated (`repeat_u`, `repeat_v`).

**Rules.** A texture keeps its attribute's rule and is baked and read by it: colours mixed in their space (SRGB, OKLAB),
normals unit, quaternions on one side, angles unwrapped, NEAREST labels never mixed (the nearest texel).

**Baking.** Each face (fanned from its first corner) is rasterized in (u, v): a texel whose centre lies in a triangle takes
the triangle's corner values blended at the centre's weights. A triangle across a repeated axis' edge is unwrapped to its
first corner's copy. Then `gutter` rings grow past the charts, each texel the mean of its written neighbours, so bilinear
reads and mipmaps near a chart's border do not pull in texels of no chart. Other texels are left as they were.
`out_coverage` (may be NULL) gets 1 for a texel in a triangle, 2 in the gutter, 0 untouched.

**From another mesh.** `qaws_mesh_3d_bake_texture_from` takes each texel's point on the mesh to the closest point of a
source mesh (a BVH, in one batch) and the source's attribute there by its rule: a detailed mesh's colours or normals onto
a decimated one's texture. `max_distance` leaves texels farther from the source untouched.

**Reading.** Bilinear between the four nearest centres, by the rule. `qaws_mesh_3d_sample_texture` writes a corner attribute
(new, or one of that name with as many components), only on the faces the desc selects.

**Which (u, v).** A corner (or vertex) attribute of two components, `"uv"` unless named; `qaws_mesh_3d_parameterize` makes
one (charts, seams, packing). It is carried through Booleans, remeshing and decimation like any corner attribute, so the
texture stays valid on the meshes made from it.
- `texture_id`: a one-component face attribute; only the faces whose value is `id` are baked or read. A Boolean of two
  textured meshes keeps each face's id, so each operand's texture is read onto its own faces.
- `footpoints`: a texture over a surface patch's own (u, v) instead. Each vertex's footpoint (generated meshes, and the
  meshes Booleans make of them) is mapped from the patch's `domain` (u0, u1, v0, v1) to the unit square, and only faces
  on `patch` count. Such a texture belongs to the surface: every mesh generated from it reads it alike. A surface closed
  in u or v takes `repeat_u` or `repeat_v`. A vertex on a border between patches has one footpoint; given the patch's
  `surface`, a face with a corner on the patch takes its other corners' (u, v) as their closest points on the surface,
  so the faces along the border count too.

**Sweeps per texel.** `qaws_mesh_3d_sweep_texture` (`qaws_mesh_sweep.h`) evaluates a sweep's generated attribute (`time`,
`index`, `angle`, `orientation`, `pose_position`, `velocity`, `distance`, `source_face`, `source_uv`, `passes`, `dwell`) or
a carried body attribute at each texel's point of the swept solid. Each point is taken to the moving body as the implicit
sweep takes its vertices (first touch, every pass when wanted, tracks and combines), within twice the sweep's tolerance.
The swept mesh may be any mesh of the solid with a (u, v).

**Tests.**
- Test 187: a linear field baked exactly (0 at the centres, 9e-16 between them, exact onto a finer mesh); labels never
  mixed; normals unit (2e-16); the gutter's rings, a repeated axis, a texture id, the refusals.
- Test 187: from a 32 × 32 grid onto two triangles' 32 × 32 texture, each texel the fine mesh's own interpolation (1e-16)
  where the two triangles' vertices are 1.39 off; `max_distance`.
- Test 187: a union of two textured boxes read by id (A's linear texture exact on its corners, cut ones too); a generated
  cylinder's patch texture read through its footpoints within the bilinear error (1.2e-3), vertices without footpoints
  read through the surface.
- Test 184: a cube slowing to a stop along x, its swept box of 12 triangles: the time per texel within 0.0005 of
  1 − (2 − x)^(1/3) on the top, the vertices' exact times interpolated 0.54 off.

**Figures.**
- Figure 81 (key `.`): an 8064-face ball's checker on a 400-face decimation, 87% of points off from the vertices, 23%
  (cell borders) from a 256 × 256 texture baked from the detailed ball.
- Figure 82 (key `,`): the slowing cube's swept box coloured by time from its 8 vertices (0.54 off) and from 128 × 128
  texels (0.0011).
