### This file contains public API ###
# Naive

"""
    Naive

Allow to run the ray tracer without an acceleration structure. This should be
assigned to the argument `acceleration` in the `RayTracer` function.
"""
struct Naive{FT} <: Acceleration{FT}
    gbox::AABB{FT}
    triangles::Vector{Triangle{FT}}
    id::Vector{Int}
end

"""
    Naive(tris::Vector{Triangle{FT}}, ids::Vector{Int}, rule) where {FT}

Wrap the scene into a global bounding box (no acceleration), given the triangles in a scene,
the ids linking each triangle to the corresponding material objects. The argument `rule` is
left for compatibility with `BVH`, it does not do anything.

## Examples
```jldoctest
julia> using PlantGeomPrimitives

julia> tris = PlantRayTracer.Triangle(Ellipse());

julia> ids  = repeat([1], length(tris));

julia> Naive(tris, ids);
```
"""
function Naive(triangles, ids, rule = nothing)
    gbox = AABB(triangles)
    Naive(gbox, triangles, ids)
end

# Return closest hit (if any)
function intersect(ray::Ray{FT}, acc::Naive, nodestack, dstack, dmin) where {FT}
    @inbounds begin
        #dmin = Inf
        frontmin = true
        posmin = -1
        for i in eachindex(acc.triangles)
            triangle = acc.triangles[i]
            hit, d, front = intersect(ray, triangle)
            if hit && d <= dmin
                dmin = d
                frontmin = front
                posmin = i
            end
        end
        if posmin == -1
            return false, Intersection(FT), dmin
        else
            triangle = acc.triangles[posmin]
            intersection = Intersection(ray.o .+ dmin .* ray.dir, # pint
                axes(triangle),         # axes
                frontmin,               # front
                acc.id[posmin])         # material
            return true, intersection, dmin
        end
    end
end
