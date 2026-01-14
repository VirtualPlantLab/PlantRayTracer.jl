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

# Abstract type for materials
abstract type Material end

"""
    materials(mesh::PGP.Mesh)

Extract the materials properties from a given mesh.

# Arguments
- `mesh::PGP.Mesh`: The mesh object from which to extract material properties.

# Returns
- A dictionary containing the material properties of the mesh.

# Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> r = Rectangle();

julia> add_property!(r, :materials, Lambertian(τ = 0.1, ρ = 0.2));

julia> materials(r);
```
"""
function materials(mesh::PGP.Mesh)
    return PGP.properties(mesh)[:materials]
end


"""
    has_materials(mesh::Mesh)

Check if a mesh contains materials among its properties.

# Arguments
- `mesh`: The mesh object being tested,

# Returns
- A boolean (`true` or `false`)

# Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> r = Rectangle();

julia> add_property!(r, :materials, Lambertian(τ = 0.1, ρ = 0.2));

julia> has_materials(r);

julia> r2 = Rectangle();

julia> has_materials(r2);
```
"""
function has_materials(mesh::PGP.Mesh)
    :materials in keys(PGP.properties(mesh))
end

# Implementations of different types of materials
include("Lambertian.jl")
include("Phong.jl")
include("Sensor.jl")
include("TwoSidedSensor.jl")
include("Black.jl")

"""
    reset!(material::Material)

Reset the power stored inside a material back to zero

## Examples
```jldoctest
julia> l = Lambertian(τ = 0.1, ρ = 0.2);

julia> l.power[1] = 10.0;

julia> reset!(l);

julia> l;
```
"""
function reset!(material::Material)
    material.power .= 0.0
    return nothing
end
function reset!(materials::Vector{<:Material})
    T.@threads for mat in materials
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
    mat.power
end
