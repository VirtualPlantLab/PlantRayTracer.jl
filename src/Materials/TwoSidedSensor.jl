### This file contains public API ###
# TwoSidedSensor

#=
  Material that has no effect on the ray and its associated power but accumulates
the power associated to rays hitting on the front and back side
=#
struct TwoSidedSensor{nw} <: Material
    power_front::MVector{nw, Float64}
    power_back::MVector{nw, Float64}
end

"""
    TwoSidedSensor(nw::Int)

Create a sensor material object to store power for `nw` wavelengths. A two-sided sensor
material will let rays pass through without altering the direction or irradiance.
They will also not count for the total number of ray iterations. The accumulated power
is stored separately for rays hitting on the front and back side of the sensor.

## Examples
```jldoctest
julia> s = TwoSidedSensor(1);

julia> s = TwoSidedSensor(3);
```
"""
function TwoSidedSensor(nw::Int = 1)
    power_front = @MVector zeros(nw)
    power_back = @MVector zeros(nw)
    TwoSidedSensor(power_front, power_back)
end



###############################################################################
################################## API ########################################
###############################################################################

#=
    There is actually no interaction with the sensor but we do keep track of the side the
    ray is hitting (we reuse the coef variable for this purpose to help the compiler by
    ensuring that we always return the same type of tuple, although it is not strictly necessary)
=#
function calculate_interaction(material::TwoSidedSensor, power, ray, intersection, rng)
    return (mode = :sensor, coef = 1.0*intersection.front, θ = 1.0, Φ = 1.0)
end

#=
    Add all the power to the material but do not affect the power of the ray. Distinguish
    between front and back side (see calculate_interaction above)
=#
function absorb_power!(material::TwoSidedSensor, power, interaction)
    if interaction.coef > 0.0
        @inbounds for i in eachindex(power)
            @atomic material.power_front[i] += power[i]
        end
    else
        @inbounds for i in eachindex(power)
            @atomic material.power_back[i] += power[i]
        end
    end
    return nothing
end

#=
    Just return the same ray but displaced a certain distance
=#
function generate_ray(material::TwoSidedSensor,
    ray::Ray{FT},
    disp::Vec{FT},
    intersection,
    interaction,
    rng) where {FT}
    Ray(intersection.pint .+ ray.dir .* sqrt.(eps.(FT)), ray.dir, ray.idir, ray.extra)
end


"""
    power(material::TwoSidedSensor; front = true)

Extract the power stored inside a two sided material.
```
"""
function power(mat::TwoSidedSensor; front = true)
    if front
        return mat.power_front
    else
        return mat.power_back
    end
end
