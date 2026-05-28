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
# angles θ and  Φ
function polar_to_cartesian(axes, θ, Φ)
    e1, e2, n = axes
    dir1 = @. n * cos(θ)
    dir2 = @. e1 * cos(Φ) * sin(θ)
    dir3 = @. e2 * sin(Φ) * sin(θ)
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
 Given angles θ (solar zenith), Φ (solar azimuth), α (azimuth of X axis)
 calculate the flipped coordinate system of the plane (i.e. a coordinate system
 for rays generated from the sky towards the Earth)
 By default we assume:
 X axis points towards South (azimuth of pi radians)
 Y axis points towards East (azimuth of pi/2 radians)
 Z axis points towards zenith (away from ground)
 Φ clockwise looking against Z - East is positive
 θ counterclock wise - 0 = Zenith, Sunrise = pi/2
 α is the angle between the X axis and the North-South axis (i.e. the azimuth of X axis)
 For example: α = pi/2 means the X axis points towards East and Y axis points towards North
 PGP = PlantGeomPrimtives.jl package (just generates unit vectors of a Cartesian system)
=#
function rotate_coordinates(θ::FT, Φ::FT, α = FT(π)) where {FT}
    rot = rotatez(α - FT(π) - Φ) ∘ rotatey(-θ)
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
