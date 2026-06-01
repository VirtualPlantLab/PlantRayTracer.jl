# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**PlantRayTracer.jl** is a ray tracer for simulating light interception by plants, part of the [VirtualPlantLab](https://virtualplantlab.com) (VPL) ecosystem. It models light as rays interacting with triangulated meshes through physically-based materials (Lambertian, Phong, sensors), acceleration structures (BVH, Naive), and light sources (point, line, area, directional).

- **Julia ≥ 1.12** required
- **Package version:** 0.1.1
- **Code style:** SciML (configured in `.JuliaFormatter.toml`)

### Key Design Patterns

**Multi-wavelength via parametric types:** Materials and sources are parameterized by `nw` (number of wavelengths), stored as `SVector{nw, Float64}` from StaticArrays. This ensures type stability across wavelengths without dynamic dispatch.

**Thread-safe material accumulation:** During parallel ray tracing, `Atomix.@atomic` is used when updating absorbed/intercepted power in materials — do not replace with plain assignment.

**Ray struct pre-computation:** `Ray` stores `inv_direction` (1/direction component-wise) computed at construction. The AABB slab intersection test depends on this; never recompute inline.

**Acceleration via `accelerate`:** Calling `accelerate(mesh, rule)` converts a `PlantGeomPrimitives` mesh into an `AccMesh` with the chosen acceleration structure (`Naive()` or `BVH(SAH())` / `BVH(AvgSplit())`).

**Grid cloning:** `GridCloner` (in `Cloner.jl`) tiles meshes horizontally (x/y) and optionally vertically (z) at trace time to simulate infinite/periodic canopies. This is configured via `RTSettings` (`nx`, `ny`, `nz`, `dx`, `dy`, `dz`).

**Russian roulette:** After `RTSettings.maxiter` bounces, rays are terminated probabilistically with survival probability `RTSettings.pkill` to avoid infinite scattering loops.


### VPL Ecosystem Integration

This package consumes geometry from **PlantGeomPrimitives.jl** (triangulated meshes) and is typically used via the **VirtualPlantLab.jl** umbrella package. Tests import both packages; changes to primitive types in PlantGeomPrimitives may require updates here.


## CI

GitHub Actions (`.github/workflows/CI.yml`) runs tests on Julia 1.12 / ubuntu-latest on pushes to `master` and all tags. Coverage is uploaded to Codecov.
