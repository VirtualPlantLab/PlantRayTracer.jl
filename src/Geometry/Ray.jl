### This file DOES NOT contain public API ###

# Ray
# o = Origin of the Ray
# dir = Direction unit vector
# idir = Inverse of dir (precalculation for ray - triangle intersection)
# extra = - o * idir (precalculation for ray - aabb intersection)
struct Ray{FT}
    o::PGP.Vec{FT}
    dir::PGP.Vec{FT}
    idir::PGP.Vec{FT}
    extra::PGP.Vec{FT}
end

# Construct Ray from point of origin and direction
function Ray(o::PGP.Vec{FT}, dir::PGP.Vec{FT}) where {FT}
    idir = one(FT) ./ dir
    Ray(o, dir, idir, .-o .* idir)
end
