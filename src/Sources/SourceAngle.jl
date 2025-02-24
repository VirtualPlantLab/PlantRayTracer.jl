### This file contains public API ###
# FixedSource
# LambertianSource

"""
    FixedSource(dir)
    FixedSource(θ, Φ)

Create a fixed irradiance source by given a vector with the direction of the
rays (dir) or zenith (θ) and azimuth (Φ) angles.

## Examples
```jldoctest
julia> source_dir = FixedSource(0.0, 0.0);

julia> import PlantGeomPrimitives as PG

julia> source_dir = FixedSource(PG.Vec(0.0, 0.0, -1.0));
```
"""
struct FixedSource{FT} <: SourceAngle
    dir::PGP.Vec{FT}
end
function FixedSource(θ, Φ)
    FixedSource(rotate_coordinates(θ, Φ).z)
end
generate_direction(a::FixedSource, rng) = a.dir

# Emission of diffuser (to be used with general light sources and for thermal radiation)
"""
    LambertianSource(x, y, z)
    LambertianSource(axes)

Create a Lambertian irradiance source angle by given a local coordinate system as three
separate `Vec` objects representing the axes (x, y, z) or as tuple containing the three
axes. Rays will be generated towards the hemisphere defined by the `z` direction. See VPL
documentation for details on irradiance sources.

## Examples
```jldoctest
julia> import PlantGeomPrimitives as PG

julia> source_dir = LambertianSource(PG.X(), PG.Y(), PG.Z());

julia> source_dir = LambertianSource((PG.X(), PG.Y(), PG.Z()));
```
"""
struct LambertianSource{FT} <: SourceAngle
    x::PGP.Vec{FT}
    y::PGP.Vec{FT}
    z::PGP.Vec{FT}
end
LambertianSource(axes) = LambertianSource(axes...)

function generate_direction(a::LambertianSource{FT}, rng) where {FT}
    Φ = FT(2) * FT(π) * rand(rng, FT)
    θ = acos(sqrt(rand(rng, FT)))
    polar_to_cartesian((e1 = a.x, e2 = a.y, n = a.z), θ, Φ)
end
