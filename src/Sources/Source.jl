### This file contains public API ###
# Source
# get_nw

#= Define the following methods for a SourceGeometry:
    generate_point(g::SourceGeometry, rng) -> Return a Vec representing origin of the ray
=#
abstract type SourceGeometry end
include("SourceGeometry.jl")
include("Directional.jl")

#= Define the following methods for a SourceAngle:
    generate_direction(a::SourceAngle, rng) -> Return a Vec representing the unit vector of the direction of the ray
=#
abstract type SourceAngle end
include("SourceAngle.jl")

#= A source is a combination of:
    - Geometry: Samples the origin of the ray
    - Angle: Samples the direction of the ray
    - power: Power of each ray (irradiance × area/nrays)
    - Number of rays shot by the source
=#
struct Source{G, A, nw}
    geom::G
    angle::A
    power::SA.SVector{nw, Float64}
    nrays::Int
end

"""
    Source(geom, angle, power::Number, nrays)
    Source(geom, angle, power::Tuple, nrays)

Create an irradiance source given a source geometry, a source angle, the power per ray and
the total number of rays to be generated from this source. When simulating more than one
wavelength simultaneously, a tuple of power values should be given, of the same length as
in the materials used in the scene. See VPL documentation for details on source geometries
and source angles.

## Examples
```jldoctest
julia> import PlantGeomPrimitives as PG

julia> source_geom = PointSource(PG.O());

julia> source_dir = FixedSource(0.0, 0.0);

julia> Source(source_geom, source_dir, 1.0, 1_000);
```
"""
function Source(geom, angle, power::Number, nrays::Integer)
    Source(geom, angle, SA.SVector{1, Float64}(power), nrays)
end
function Source(geom, angle, power::Tuple, nrays::Integer)
    Source(geom, angle, SA.SVector{length(power), Float64}(power...), nrays)
end

"""
    get_nw(s::Source)

Retrieve the number of wavelengths that rays from a source will contain.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> mesh = Ellipse();

julia> source = DirectionalSource(mesh, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);

julia> get_nw(source);
```
"""
get_nw(s::Source{G, A, nw}) where {G, A, nw} = nw

# Forward methods to the components
generate_point(source::Source, rng) = generate_point(source.geom, rng)
generate_direction(source::Source, rng) = generate_direction(source.angle, rng)

# Sample ray from a light source using Monte Carlo
function shoot_ray!(source::Source, power, rng)
    # Transfer power from the source to the ray payload
    power .= source.power
    # Generate the origin and direction of the ray
    origin = generate_point(source, rng)
    dir = generate_direction(source, rng)
    return Ray(origin, dir)
end
