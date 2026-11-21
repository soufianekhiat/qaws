# Qaws API Reference

One page per public header area. Every precision-bearing function and type exists as `name_f32` and
`name_f64`; the pages document the unsuffixed names, which mean the precision selected by
`QAWS_SCALAR_IS_FLOAT` (see [precision.md](precision.md)). Include `qaws.h` for the whole API, or
`qaws_both.h` for both precisions in one translation unit.

For concepts and design, start at [../index.md](../index.md).

## Core

| Page | Headers | Contents |
|---|---|---|
| [types.md](types.md) | `qaws_types.h` | Scalars, vectors, enums, shared structs |
| [status.md](status.md) | `qaws_status.h` | Status codes and `qaws_status_to_string` |
| [precision.md](precision.md) | `qaws_platform.h`, `qaws_precision.h`, `qaws_names.h`, `qaws_both.h` | Backend and precision selection, the `_f32` / `_f64` naming |
| [alloc.md](alloc.md) | `qaws_alloc.h`, `qaws_allocator.h` | Allocators, settings, thread allocators, arenas |
| [batch_executor.md](batch_executor.md) | `qaws_batch_executor.h` | Host-supplied parallel executor |
| [curve.md](curve.md) | `qaws_curve.h` | Curve lifecycle |

## Curve families

| Page | Header |
|---|---|
| [bezier.md](bezier.md) | `qaws_bezier.h` |
| [rational_bezier.md](rational_bezier.md) | `qaws_rational_bezier.h` |
| [hermite.md](hermite.md) | `qaws_hermite.h` |
| [catmull_rom.md](catmull_rom.md) | `qaws_catmull_rom.h` |
| [bspline.md](bspline.md) | `qaws_bspline.h` |
| [nurbs.md](nurbs.md) | `qaws_nurbs.h` |
| [trajectory.md](trajectory.md) | `qaws_trajectory.h` |
| [yuksel.md](yuksel.md) | `qaws_yuksel.h` |
| [arc.md](arc.md) | `qaws_arc.h` |
| [polynomial.md](polynomial.md) | `qaws_polynomial.h` |
| [clothoid.md](clothoid.md) | `qaws_clothoid.h` |
| [subdivision.md](subdivision.md) | `qaws_subdivision.h` |
| [composite.md](composite.md) | `qaws_composite.h` |
| [interpolating.md](interpolating.md) | `qaws_interpolating.h` |
| [inline.md](inline.md) | `qaws_inline.h` (heap-free curves) |

## Evaluation, sampling, traversal, inspection

| Page | Header |
|---|---|
| [eval.md](eval.md) | `qaws_eval.h` |
| [prepare.md](prepare.md) | `qaws_prepare.h` |
| [sampling.md](sampling.md) | `qaws_sampling.h` |
| [traversal.md](traversal.md) | `qaws_traversal.h` |
| [motion.md](motion.md) | `qaws_motion.h` |
| [inspect.md](inspect.md) | `qaws_inspect.h` |

## Views and batches

| Page | Header |
|---|---|
| [view.md](view.md) | `qaws_view.h` (borrowed curves and surfaces over host arrays) |
| [curve_batch.md](curve_batch.md) | `qaws_curve_batch.h` |
| [surface_batch.md](surface_batch.md) | `qaws_surface_batch.h` |

## Operations, conversion, import and export

| Page | Header |
|---|---|
| [operations.md](operations.md) | `qaws_operations.h` |
| [convert.md](convert.md) | `qaws_convert.h` |
| [export.md](export.md) | `qaws_export.h` |
| [import.md](import.md) | `qaws_import.h` |

## Surfaces

| Page | Headers |
|---|---|
| [surface.md](surface.md) | `qaws_surface.h`, `qaws_surface_types.h` |
| [surface_families.md](surface_families.md) | `qaws_surface_bezier.h`, `_bspline`, `_nurbs`, `_bilinear`, `_biquadratic`, `_swept`, `_ruled`, `_coons`, `_extrusion`, `_revolution`, `_pipe`, `_loft`, `_gordon`, `_offset` |
| [surface_trim.md](surface_trim.md) | `qaws_surface_trim.h` |
| [surface_tspline.md](surface_tspline.md) | `qaws_surface_tspline.h` |
| [surface_subdiv.md](surface_subdiv.md) | `qaws_surface_subdiv.h` |
| [brep.md](brep.md) | `qaws_brep.h` |
| [shatter.md](shatter.md) | `qaws_shatter.h` |

## 2D regions

| Page | Header |
|---|---|
| [path.md](path.md) | `qaws_path.h` |
| [clip.md](clip.md) | `qaws_clip.h` |
| [boolean_2d.md](boolean_2d.md) | `qaws_boolean_2d.h` |

## Exact kernel

| Page | Header |
|---|---|
| [exact.md](exact.md) | `qaws_exact.h` |

## Differentiation

| Page | Headers |
|---|---|
| [diff.md](diff.md) | `qaws_diff.h`, `qaws_diff_types.h`, `qaws_diff_geometry.h`, `qaws_diff_ops.h`, `qaws_diff_map.h`, `qaws_diff_sampling.h` |

## Portable core and shader backends

| Page | Headers |
|---|---|
| [core.md](core.md) | `core/qaws_*_core.h`, `qaws_core_types.h`, `qaws_hlsl.h`, `qaws_glsl.h`, `qaws_halide.h` |
