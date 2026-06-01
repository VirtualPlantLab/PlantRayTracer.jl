# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**PlantRayTracer.jl** is a ray tracer for simulating light interception by plants, part of the [VirtualPlantLab](https://virtualplantlab.com) (VPL) ecosystem. It models light as rays interacting with triangulated meshes through physically-based materials (Lambertian, Phong, sensors), acceleration structures (BVH, Naive), and light sources (point, line, area, directional).

- **Julia ≥ 1.12** required
- **Package version:** 0.1.1
- **Code style:** SciML (configured in `.JuliaFormatter.toml`)

## Commands

### Testing
```julia
# Run full test suite
julia --project=. -e "using Pkg; Pkg.test()"

# Run a single test file
julia --project=. test/elements_raytracer.jl

# Format code
julia -e "using JuliaFormatter; format(\".\")"
```

### Development
```julia
# Activate environment in REPL
julia --project=.

# Update dependencies
julia --project=. -e "using Pkg; Pkg.update()"
```

## Architecture

### Module Structure (`src/`)

```
PlantRayTracer.jl          # Root: exports public API, includes all submodules
utils.jl                   # tau/rho SVector constructors, polar_to_cartesian, GVector type
├── Geometry/
│   ├── Ray.jl             # Ray type with precomputed inverse direction
│   ├── Triangle.jl        # Triangle + Möller–Trumbore intersection
│   └── AABB.jl            # Axis-aligned bounding box (slab method)
├── Materials/
│   ├── Material.jl        # Abstract Material type
│   ├── Lambertian.jl      # Diffuse reflection/transmission
│   ├── Phong.jl           # Phong specular model
│   ├── Sensor.jl          # One-sided irradiance sensor
│   ├── TwoSidedSensor.jl  # Two-sided sensor tracking front/back power separately
│   └── Black.jl           # Perfect absorber
├── Sources/
│   ├── SourceGeometry.jl  # PointSource, LineSource, AreaSource geometries
│   ├── SourceAngle.jl     # SourceAngle distributions (LambertianSource, etc.)
│   ├── Directional.jl     # Directional geometry + DirectionalSource factory
│   └── Source.jl          # Source = geometry + angle + power + nrays; FixedSource
├── Acceleration/
│   ├── Acceleration.jl    # Intersection struct + abstract Acceleration type
│   ├── Naive.jl           # Brute-force O(n) ray-mesh intersection
│   ├── BVH.jl             # Bounding Volume Hierarchy (SAH and AvgSplit rules)
│   └── Cloner.jl          # GridCloner: infinite-canopy via horizontal/vertical tiling
└── RayTracer/
    └── RayTracer.jl       # RTSettings, AccMesh, RayTracer, accelerate, trace!, get_nw
```

### Key Design Patterns

**Multi-wavelength via parametric types:** Materials and sources are parameterized by `nw` (number of wavelengths), stored as `SVector{nw, Float64}` from StaticArrays. This ensures type stability across wavelengths without dynamic dispatch.

**Thread-safe material accumulation:** During parallel ray tracing, `Atomix.@atomic` is used when updating absorbed/intercepted power in materials — do not replace with plain assignment.

**Ray struct pre-computation:** `Ray` stores `inv_direction` (1/direction component-wise) computed at construction. The AABB slab intersection test depends on this; never recompute inline.

**Acceleration via `accelerate`:** Calling `accelerate(mesh, rule)` converts a `PlantGeomPrimitives` mesh into an `AccMesh` with the chosen acceleration structure (`Naive()` or `BVH(SAH())` / `BVH(AvgSplit())`).

**Grid cloning:** `GridCloner` (in `Cloner.jl`) tiles meshes horizontally (x/y) and optionally vertically (z) at trace time to simulate infinite/periodic canopies. This is configured via `RTSettings` (`nx`, `ny`, `nz`, `dx`, `dy`, `dz`).

**Russian roulette:** After `RTSettings.maxiter` bounces, rays are terminated probabilistically with survival probability `RTSettings.pkill` to avoid infinite scattering loops.

### Public API (from `src/PlantRayTracer.jl` exports)

```julia
RayTracer, RTSettings, trace!, accelerate, reset!, materials, has_materials,
Naive, BVH, SAH, AvgSplit,
Material, Lambertian, Phong, Sensor, TwoSidedSensor, Black,
Source, FixedSource, LambertianSource, DirectionalSource,
PointSource, LineSource, AreaSource, Directional,
power, tau, rho, get_nw
```

### VPL Ecosystem Integration

This package consumes geometry from **PlantGeomPrimitives.jl** (triangulated meshes) and is typically used via the **VirtualPlantLab.jl** umbrella package. Tests import both packages; changes to primitive types in PlantGeomPrimitives may require updates here.

## Tests

| File | Coverage |
|------|----------|
| `test/elements_raytracer.jl` | Ray, Triangle, AABB primitives and intersections |
| `test/raytracer_tiles.jl` | Tiled geometry ray tracing |
| `test/test_graphs.jl` | Full scene ray tracing |
| `test/test_shade_net.jl` | Shade cloth and radiation profiles |

Tests use `Test.jl`, `Aqua.jl` (code quality), and `Documenter.jl` (doctests).

## CI

GitHub Actions (`.github/workflows/CI.yml`) runs tests on Julia 1.12 / ubuntu-latest on pushes to `master` and all tags. Coverage is uploaded to Codecov.
