# Qaws Documentation

**Qaws** (Arabic for "arc") is a dependency-free C11 geometry library for 2D and 3D: curves and surfaces,
their evaluation, sampling, traversal and inspection, 2D regions and clipping, trimmed and subdivided
surfaces, B-rep and shatter, an exact rational kernel, and exact differentiation of all of it.

Two documentation sets:

- Concepts and design (this folder): what the library models, how the parts fit, conventions and the
  reasons behind them.
- [API reference](api/index.md): one page per public header area, every function, type and flag.

## Getting started

- [Overview](overview.md) -- Scope, design principles, quick start, repository layout
- [Building](building.md) -- Sharpmake (primary) and CMake builds, options
- [Integration](integration.md) -- Using Qaws without a build system: sources, unity builds, the amalgamation

## Geometry

- [Curve families](curve_families.md) -- Every curve type, its parameters and continuity
- [Surfaces](surfaces.md) -- Patches, swept and lofted surfaces, trims, T-splines, subdivision, offsets
- [2D regions](regions_2d.md) -- Paths, fill rules, clipping and Boolean operations
- [B-rep](brep.md) -- The boundary representation and how trimmed surfaces relate to it
- [Shatter](shatter.md) -- Shatter and seams

## Using curves and surfaces

- [Evaluation](evaluation.md) -- Position and derivatives, flags, prepared evaluation
- [Sampling](sampling.md) -- Point arrays from curves and surfaces
- [Traversal](traversal.md) -- Arc length, constant speed, motion profiles
- [Inspection](inspection.md) -- Curvature, frames, bounds, closest points, intersections
- [Operations](operations.md) -- Split, join, offset, conversion, fitting, import and export
- [Views and batches](views_and_batches.md) -- Borrowed host arrays, their lifetime, batched evaluation

## Robustness and exactness

- [Precision](precision.md) -- f32 and f64, the naming scheme, tolerances, the float robustness model
- [Exact kernel](exact_kernel.md) -- Rational arithmetic, integer budgets, certified signs
- [Differentiation](differentiation.md) -- Tangents, adjoints and second derivatives of evaluation, geometry, solves and operations

## Runtime model

- [Memory](memory.md) -- Allocators, arenas, ownership and destroy rules, heap-free curves
- [Threading](threading.md) -- No threads of its own, the batch executor, what is safe to share, cancellation
- [Error handling](error_handling.md) -- Status codes and error patterns

## Backends and coverage

- [Backends](backends.md) -- C runtime, portable core headers, HLSL, GLSL, Halide, SIMD
- [Capabilities](capabilities.md) -- What runs on which backend, and which test backs it
