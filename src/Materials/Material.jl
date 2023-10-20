### This file contains public API ###
# reset!
# power

# Methods that need to be implemented for Material types:
# choose_interaction(material, power, ray, intersection, rng) -> Choose type of interaction and calculate variable optical properties
# absorb_power!(material, power, interaction)                 -> Transfer from power to material based on type of interaction
# generate_ray(material, ray, intersection, interaction, rng) -> Generate a new ray based on type of interaction

# Intersection should contain the local reference system as unit vectors (e1, e2, n)
# To avoid type instability propagating throughout the code, all these functions
# should return tuples of the same type

# Implementations of different types of materials
include("Lambertian.jl")
include("Phong.jl")
include("Sensor.jl")
include("Black.jl")

"""
    reset!(material::Material)

Reset the power stored inside a material back to zero

## Examples
```jldoctest
julia> l = Lambertian(τ = 0.1, ρ = 0.2);

julia> l.power[1].value = 10.0;

julia> reset!(l);

julia> l;
```
"""
function reset!(material::Material)
    for i in eachindex(material.power)
        @inbounds material.power[i].value = 0.0
    end
    return nothing
end
function reset!(materials::Vector{<:Material})
    @threads for mat in materials
        reset!(mat)
    end
end

"""
    power(material::Material)

Extract the power stored inside a material.

## Examples
```jldoctest
julia> l = Lambertian(τ = 0.1, ρ = 0.2);

julia> power(l);
```
"""
function power(mat::Material)
    @inbounds nw = typeof(mat).parameters[1]
    @inbounds typ = eltype(mat.power.data).parameters[1]
    SVector{nw, typ}(pow.value for pow in mat.power)
end
