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

# Help write coordinate transformations
translate(x, y, z) = CT.Translation(x, y, z)
rotatex(x) = CT.LinearMap(Rotations.RotX(x))
rotatey(x) = CT.LinearMap(Rotations.RotY(x))
rotatez(x) = CT.LinearMap(Rotations.RotZ(x))

# Given angles θ and Φ, calculate the flipped coordinate system of the plane
# Φ clockwise looking against Z - East is positive
# θ counterclock wise looking against Y - Sunrise is positive
function rotate_coordinates(θ::FT, Φ::FT) where {FT}
    rot = rotatez(-Φ) ∘ rotatex(-θ)
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
