### This file contains public API ###
# Black

#=
  Material that has no effect on the ray and its associated power but accumulates
all the power associated to the ray
=#
struct Black{nw} <: Material
    power::MVector{nw, Float64}
end

"""
    Black(nw::Int)

Create a black material object to store power for `nw` wavelengths. See VPL documentation for
details.

## Examples
```jldoctest
julia> b = Black(1);

julia> b = Black(3);
```
"""
Black(nw::Int = 1) = Black(@MVector zeros(nw))

###############################################################################
################################## API ########################################
###############################################################################

#=
    Does not really matter as no ray will come out
=#
function calculate_interaction(material::Black, power, ray, intersection, rng)
    return (mode = :black, coef = 1.0, θ = 1.0, Φ = 1.0)
end

#=
    Add all the power to the material and set power to 0
=#
function absorb_power!(material::Black, power, interaction)
    @inbounds for i in eachindex(power)
        @atomic material.power[i] += power[i]
    end
    power .= 0.0
    return nothing
end

#=
    Just return a ray (we will never get here due to Russian roulette if it is black surface)
=#
function generate_ray(material::Black,
    ray::Vec{FT},
    disp::Vec{FT},
    intersection,
    interaction,
    rng) where {FT}
    Ray(O(FT), O(FT))
end
