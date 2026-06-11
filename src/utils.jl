### This file contains public API ###
# tau
# rho

# Interface to generate SVectors more easily (for the materials)
"""
    tau(vals...)

Generate values of transmisivity to be used in material object. `vals...` is a
list of one or more comma-separated values, corresponding to the different
wavelengths/wavebands to be simulated in a ray tracer.

## Examples
```jldoctest
julia> tau(1.0, 0.0, 2.0);
```
"""
tau(vals...) = SA.SVector(vals)

"""
    rho(vals...)

Generate values of reflectivity to be used in material object. `vals...` is a
list of one or more comma-separated values, corresponding to the different
wavelengths/wavebands to be simulated in a ray tracer.

## Examples
```jldoctest
julia> rho(1.0, 0.0, 2.0);
```
"""
rho(vals...) = SA.SVector(vals)

###############################################################################
################################## Geometry ###################################
###############################################################################

# Calculate direction unit vector on Cartesian system axes from the
# angles θ and Φ (in degrees)
function polar_to_cartesian(axes, θ, Φ)
    e1, e2, n = axes
    dir1 = @. n * cosd(θ)
    dir2 = @. e1 * cosd(Φ) * sind(θ)
    dir3 = @. e2 * sind(Φ) * sind(θ)
    dir = @. dir1 + dir2 + dir3
    PGP.Vec(normalize(dir))
end

# Project a point (p) onto a plane (defined by point pp and normal pn)
function project(p, pp, pn)
    v = p .- pp
    dist = v ⋅ pn
    p .- dist .* pn
end

# Help write coordinate transformations (CT = CoordinateTransformations.jl)
translate(x, y, z) = CT.Translation(x, y, z)
rotatex(x) = CT.LinearMap(Rotations.RotX(x))
rotatey(x) = CT.LinearMap(Rotations.RotY(x))
rotatez(x) = CT.LinearMap(Rotations.RotZ(x))

#=
 Given angles θ (solar zenith), Φ (solar azimuth), α (azimuth of X axis),
 alpha_soil (slope inclination) and beta_soil (azimuth of slope normal),
 calculate the flipped coordinate system of the plane (i.e. a coordinate system
 for rays generated from the sky towards the Earth).
 All angles are in degrees. By default we assume:
 X axis points towards South (azimuth of 180 degrees)
 Y axis points towards East (azimuth of 90 degrees)
 Z axis points towards zenith (away from ground)
 Φ clockwise looking against Z - East is positive
 θ counterclockwise - 0 = Zenith, Sunrise = 90
 α is the geocentric azimuth of the X axis (e.g. α = 90 means X points East)
 alpha_soil is the inclination of the slope (0 = horizontal, 90 = vertical)
 beta_soil is the geocentric azimuth of the slope normal (e.g. 180 = south-facing slope)
 PGP = PlantGeomPrimitives.jl package (just generates unit vectors of a Cartesian system)
=#
function rotate_coordinates(θ, Φ, α = 180.0, alpha_soil = 0.0, beta_soil = 180.0)
    # To deal with different precision levels
    FT = promote_type(typeof(θ), typeof(Φ), typeof(α), typeof(alpha_soil), typeof(beta_soil))
    rot = rotatez(deg2rad(α - beta_soil)) ∘ rotatey(deg2rad(-alpha_soil)) ∘ rotatez(deg2rad(beta_soil - FT(180) - Φ)) ∘ rotatey(deg2rad(-θ))
    (x = .-rot(PGP.X(FT)), y = .-rot(PGP.Y(FT)), z = .-rot(PGP.Z(FT)))
end

###############################################################################
################################## GVector ####################################
###############################################################################

# Vector wrapper that automatically grows the internal vector when indexing
# out-of bounds. Useful when inserting elements out-of-order as in the flat
# dense representation of a n-arity tree
const chunk = 15

struct GVector{T}
    data::Vector{T}
end

Base.IndexStyle(::Type{<:GVector}) = IndexLinear()

@inline function Base.getindex(v::GVector, i::Int)
    @boundscheck checkbounds(v.data, i)
    @inbounds return v.data[i]
end

@inline function Base.setindex!(v::GVector, x, i::Int)
    if i > length(v.data)
        resize!(v.data, i + chunk)
        @inbounds v.data[i] = x
    else
        @inbounds v.data[i] = x
    end
end

# Some methods for convenience
Base.size(v::GVector) = size(v.data)
Base.length(v::GVector) = length(v.data)
