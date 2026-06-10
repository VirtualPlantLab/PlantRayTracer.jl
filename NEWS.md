# PlantRayTracer release notes

We started keeping track of changes in the `NEWS.md` file after version 0.1.0.

# PlantRayTracer 0.2.0 (2026-06-10)

* `DirectionalSource` and `FixedSource` now accept an optional `α` keyword argument
  (azimuth of the X axis, default `π`) to support crop rows not aligned North-South.
  The internal `rotate_coordinates` function was updated accordingly.

* `DirectionalSource` and `FixedSource` now accept optional `alpha_soil` (slope
  inclination, default `0`) and `beta_soil` (azimuth of slope normal, default `π`)
  keyword arguments to simulate light interception on tilted, oriented soil surfaces.

* Fixed a type error in `FixedSource` where a free type variable `FT` was incorrectly
  referenced after the signature was refactored.

* Fixed numerical issues with mixed floating-point precisions and `Irrational` values
  of `π` by using explicit type promotion throughout `rotate_coordinates`, `FixedSource`,
  and `DirectionalSource`.

# PlantRayTracer 0.1.1 (2026-01-14)

* Update dependencies and make sure it works on Julia 1.12
