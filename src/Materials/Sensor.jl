### This file contains public API ###
# Sensor

#=
  Material that has no effect on the ray and its associated power but accumulates
all the power associated to the ray
=#
struct Sensor{nw} <: Material
    power::MVector{nw, Float64}
end

"""
    Sensor(nw::Int)

Create a sensor material object to store power for `nw` wavelengths. A sensor
material will let rays pass through without altering the direction or irradiance.
They will also not count for the total number of ray iterations.

## Examples
```jldoctest
julia> s = Sensor(1);

julia> s = Sensor(3);
```
"""
Sensor(nw::Int = 1) = Sensor(@MVector zeros(nw))

###############################################################################
################################## API ########################################
###############################################################################

#=
    There is actually no interaction with the sensor
=#
function calculate_interaction(material::Sensor, power, ray, intersection, rng)
    return (mode = :sensor, coef = 1.0, θ = 1.0, Φ = 1.0)
end

# TODO: A sensor should always distinguish front and back
#=
    Add all the power to the material but do not affect the power of the ray
=#
function absorb_power!(material::Sensor, power, interaction)
    @inbounds for i in eachindex(power)
        @atomic material.power[i] += power[i]
    end
    return nothing
end

#=
    Just return the same ray but displaced a certain distance
=#
function generate_ray(material::Sensor,
    ray::Ray{FT},
    disp::Vec{FT},
    intersection,
    interaction,
    rng) where {FT}
    Ray(intersection.pint .+ ray.dir .* eps(FT) .* FT(2), ray.dir, ray.idir, ray.extra)
end
