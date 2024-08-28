### This file contains public API ###
# Triangle

# Triangle
# The triangle is defined in barycentric coordinates to speed up ray-triangle
# intersections.
# The normal vector is also included as well as the index of the material object
# p = A vertex on the triangle
# e1 = An edge of the triangle having p as origin (length = side of triangle)
# e2 = Like e1, but for the other side of the triangle
# n = Normal unit vector (defined by cross product of e1 and e2)
struct Triangle{FT}
    p::Vec{FT}
    e1::Vec{FT}
    e2::Vec{FT}
    n::Vec{FT}
end

# Compute edges and normal vector of a triangle for barycentric coordinates
function barycentric(p1, p2, p3)
    e1 = p2 .- p1
    e2 = p3 .- p1
    n = normalize(cross(e2, e1))
    return e1, e2, n
end

"""
    Triangle()

Create a ray tracing `Triangle` object with default vertices (unit vectors in each axis).

## Examples
```jldoctest
julia> t = PlantRayTracer.Triangle();
```
"""
function Triangle()
    Triangle(X(), Y(), Z())
end

"""
    Triangle(p1::Vec, p2::Vec, p3::Vec)

Create a ray tracing `Triangle` object given the three vertices `p1`, `p2` and `p3`.

## Examples
```jldoctest
julia> using PlantGeomPrimitives

julia> t = PlantRayTracer.Triangle(Vec(1.0, 0.0, 1.0), Vec(0.0, 1.0, .0), Vec(1.0, 1.0, 1.0));
```
"""
function Triangle(p1::Vec, p2::Vec, p3::Vec)
    e1, e2, n = barycentric(p1, p2, p3)
    Triangle(p1, e1, e2, n)
end

"""
    Triangle(mesh::Mesh)

Create a vector of ray tracing `Triangle` objects from a `Mesh` object.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> e = Ellipse();

julia> t = PlantRayTracer.Triangle(e);
```
"""
function Triangle(mesh::Mesh{VT}) where {VT}
    update_normals!(mesh)
    FT = eltype(VT)
    output = Triangle{FT}[]
    sizehint!(output, ntriangles(mesh))
    for i in 1:ntriangles(mesh)
        j  = (i - 1)*3
        p1 = vertices(mesh)[j + 1]
        p2 = vertices(mesh)[j + 2]
        p3 = vertices(mesh)[j + 3]
        e1 = p2 .- p1
        e2 = p3 .- p1
        push!(output, Triangle(p1, e1, e2, normals(mesh)[i]))
    end
    return output
end

"""
    Triangle(scene::Scene)

Create a vector of ray tracing `Triangle` objects from a `Scene` object.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> sc = Scene(mesh = Ellipse());

julia> t = PlantRayTracer.Triangle(sc);
```
"""
function Triangle(scene::Scene)
    Triangle(mesh(scene))
end

# Moller-Trumbore intersection test with early exits
# Returns intersection_test (T/F), intersection_distance, front (T/F)
# The front of the triangle is simply the side the normal uVec points to
# Note: This will fail when the ray hits exactly on the border of the triangle
function intersect(ray::Ray{FT}, t::Triangle{FT}) where {FT}
    # Check if the ray intercepts the plane containing the triangle
    pvec = ray.dir × t.e2
    det = t.e1 ⋅ pvec
    idet = one(FT) / det
    # Calculate coordinates of point in barycentric coordinates
    tvec = ray.o - t.p
    u = idet * (tvec ⋅ pvec)
    (u <= zero(FT) || u >= one(FT)) && (return (false, zero(FT), true))
    qvec = tvec × t.e1
    v = idet * (ray.dir ⋅ qvec)
    (v <= zero(FT) || u + v >= one(FT)) && (return (false, zero(FT), true))
    # Calculate distance to hit
    d = idet * (t.e2 ⋅ qvec)
    (d <= zero(FT)) && (return (false, zero(FT), true))
    return true, d, det < zero(FT)
end

# Return a local system of coordinates defined by two vectors (third axes is the cross product of these two)
function axes(t::Triangle)
    (normalize(t.e1), normalize(t.n × t.e1), t.n) # e1, e2, n
end

# Sample a point from the parallelogram and then apply reflection
function generate_point(t::Triangle{FT}, rng) where {FT}
    u1 = rand(rng, FT)
    u2 = rand(rng, FT)
    if u1 + u2 > FT(1)
        u1 = FT(1) - u1
        u2 = FT(1) - u2
    end
    @. t.p + t.e1 * u1 + t.e2 * u2
end

# Calculate area of a triangle (the norm of cross product is the area of the parallelogram defined by those vectors)
area(t::Triangle{FT}) where {FT} = FT(0.5) * norm(t.e1 × t.e2)
