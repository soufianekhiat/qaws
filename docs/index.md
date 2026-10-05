# Qaws Documentation

**Qaws** (قوس) -- Dependency-free C11 curve evaluation and traversal library for 2D and 3D.

## Guides

- [Overview](overview.md) -- What Qaws is, design principles, quick start
- [Building](building.md) -- CMake and Sharpmake build instructions
- [Integration](integration.md) -- Using Qaws without CMake: dropping the sources in, unity builds, the single-file amalgamation
- [Curve Families](curve_families.md) -- All curve types with examples
- [Surfaces](surfaces.md) -- Bezier, B-spline, NURBS, swept, ruled surfaces
- [Evaluation](evaluation.md) -- Position and derivative evaluation
- [Sampling](sampling.md) -- Generating point arrays from curves
- [Traversal](traversal.md) -- Motion profiles, constant-speed, arc-length
- [Inspection](inspection.md) -- Querying curve properties and analysis
- [Differentiation](differentiation.md) -- Exact tangents, adjoints and second derivatives of evaluation, geometry, solves, operations and functionals
- [Backends](backends.md) -- Multi-backend architecture (C, HLSL, GLSL, Halide)
- [Error Handling](error_handling.md) -- Status codes and error patterns

## Feature matrix

- [Feature per Backend](feature_per_backends.md) -- Supported operations for each curve type

## API reference

Per-header documentation with all function signatures and parameters:

### Core

- [qaws_types.h](api/types.md) -- Scalar, vector, enum, and struct definitions
- [qaws_status.h](api/status.md) -- Error codes
- [qaws_curve.h](api/curve.md) -- Curve lifecycle (destroy)

### Evaluation, sampling, traversal, inspection

- [qaws_eval.h](api/eval.md) -- Evaluation functions
- [qaws_sampling.h](api/sampling.md) -- Sampling functions
- [qaws_traversal.h](api/traversal.md) -- Traversal functions
- [qaws_inspect.h](api/inspect.md) -- Inspection and analysis functions

### Curve creation

- [qaws_bezier.h](api/bezier.md) -- Bezier creation
- [qaws_hermite.h](api/hermite.md) -- Hermite creation
- [qaws_catmull_rom.h](api/catmull_rom.md) -- Catmull-Rom creation
- [qaws_bspline.h](api/bspline.md) -- B-Spline creation
- [qaws_nurbs.h](api/nurbs.md) -- NURBS creation
- [qaws_trajectory.h](api/trajectory.md) -- Trajectory creation
- [qaws_yuksel.h](api/yuksel.md) -- Yuksel C2 spline creation
- [qaws_rational_bezier.h](api/rational_bezier.md) -- Rational Bezier creation
- [qaws_arc.h](api/arc.md) -- Circular and elliptical arc creation
- [qaws_polynomial.h](api/polynomial.md) -- Polynomial curve creation
- [qaws_clothoid.h](api/clothoid.md) -- Clothoid (Euler spiral) creation
- [qaws_subdivision.h](api/subdivision.md) -- Subdivision curve creation
- [qaws_composite.h](api/composite.md) -- Composite / piecewise curve creation

### Operations, conversion, import/export

- [qaws_operations.h](api/operations.md) -- Curve operations (split, join, offset, reparameterization)
- [qaws_convert.h](api/convert.md) -- Conversion between curve types and degree manipulation
- [qaws_export.h](api/export.md) -- Serialization, export, and curve fitting
- [qaws_import.h](api/import.md) -- Deserialization / import from raw point data

### Utilities

- [qaws_alloc.h](api/alloc.md) -- Custom allocator support
- [qaws_inline.h](api/inline.md) -- Stack-allocated inline curve initialization
- [qaws_prepare.h](api/prepare.md) -- Precomputation / preparation utilities

### Differentiation

- [qaws_diff*.h](api/diff.md) -- Fields, views, jets, tangents, adjoints, geometry, implicit operations, differential maps, integral functionals

### Headers without dedicated doc pages yet

The following APIs are implemented and available via their headers but do not have standalone documentation pages. See the header files in `src/qaws/` for signatures and usage.

- `qaws_surface.h` -- Surface evaluation and types (see also `qaws_surface_types.h`, `qaws_surface_bezier.h`, `qaws_surface_bspline.h`, `qaws_surface_nurbs.h`, `qaws_surface_swept.h`, `qaws_surface_ruled.h`)
