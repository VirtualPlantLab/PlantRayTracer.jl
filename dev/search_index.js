var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PlantRayTracer","category":"page"},{"location":"#PlantRayTracer","page":"Home","title":"PlantRayTracer","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PlantRayTracer.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PlantRayTracer]","category":"page"},{"location":"#PlantRayTracer.AreaSource-Tuple{PlantGeomPrimitives.Mesh}","page":"Home","title":"PlantRayTracer.AreaSource","text":"AreaSource(mesh)\n\nCreate an area irradiance source geometry given a triangular mesh.\n\nExamples\n\njulia> using PlantGeomPrimitives\n\njulia> e = Ellipse();\n\njulia> source_geom = AreaSource(e);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.AvgSplit","page":"Home","title":"PlantRayTracer.AvgSplit","text":"AvgSplit(minN, maxL)\n\nRule to be used in RayTracer when acceleration = BVH. It will divide each node along the longest axis through the mean coordinate value. The rule is parameterized by the minimum number of triangles in a leaf node (minN) and the maximum depth of the tree (maxL).\n\nExamples\n\njulia> rule = AvgSplit(1,5);\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.BVH","page":"Home","title":"PlantRayTracer.BVH","text":"BVH(gbox, nodes, tris, ids)\n\nConstruct a Bounding Volume Hierarchy around the triangular meshes to accelerate the ray tracer. This should be assigned to the argument acceleration in the RayTracer function, in combination with a corresponding rule.\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.Black","page":"Home","title":"PlantRayTracer.Black","text":"Black(nw::Int)\n\nCreate a black material object to store power for nw wavelengths. See VPL documentation for details.\n\nExamples\n\njulia> b = Black(1);\n\njulia> b = Black(3);\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.FixedSource","page":"Home","title":"PlantRayTracer.FixedSource","text":"FixedSource(dir)\nFixedSource(θ, Φ)\n\nCreate a fixed irradiance source by given a vector with the direction of the rays (dir) or zenith (θ) and azimuth (Φ) angles.\n\nExamples\n\njulia> source_dir = FixedSource(0.0, 0.0);\n\njulia> import PlantGeomPrimitives as PG\n\njulia> source_dir = FixedSource(PG.Vec(0.0, 0.0, -1.0));\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.Lambertian-Tuple{}","page":"Home","title":"PlantRayTracer.Lambertian","text":"Lambertian(;τ = 0.0, ρ = 0.0)\n\nCreate a Lambertian material object from the values of transmittance (τ) and reflectance (ρ). When more than one wavelength is being simulated, a tuple of values should be passed for each optical property (as in τ = (0.1,0.2)).\n\nExamples\n\njulia> l = Lambertian(τ = 0.1, ρ = 0.2);\n\njulia> l = Lambertian(τ = (0.1, 0.45), ρ = (0.2, 0.45));\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.LambertianSource","page":"Home","title":"PlantRayTracer.LambertianSource","text":"LambertianSource(x, y, z)\nLambertianSource(axes)\n\nCreate a Lambertian irradiance source angle by given a local coordinate system as three separate Vec objects representing the axes (x, y, z) or as tuple containing the three axes. Rays will be generated towards the hemisphere defined by the z direction. See VPL documentation for details on irradiance sources.\n\nExamples\n\njulia> import PlantGeomPrimitives as PG\n\njulia> source_dir = LambertianSource(PG.X(), PG.Y(), PG.Z());\n\njulia> source_dir = LambertianSource((PG.X(), PG.Y(), PG.Z()));\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.LineSource","page":"Home","title":"PlantRayTracer.LineSource","text":"LineSource(p, line)\n\nCreate a line irradiance source geometry given an origin (p) and a segment (line) both specified as vector of Cartesian coordinates (Vec(x, y, z)). This will create a line source between the points p and p .+ line.\n\nExamples\n\njulia> import PlantGeomPrimitives as PG\n\njulia> source_geom = LineSource(PG.Vec(1.0, 1.0, 1.0), PG.Y());\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.Naive","page":"Home","title":"PlantRayTracer.Naive","text":"Naive(tris::Vector{Triangle{FT}}, ids::Vector{Int}, rule) where {FT}\n\nWrap the scene into a global bounding box (no acceleration), given the triangles in a scene, the ids linking each triangle to the corresponding material objects. The argument rule is left for compatibility with BVH, it does not do anything.\n\nExamples\n\njulia> using PlantGeomPrimitives\n\njulia> tris = PlantRayTracer.Triangle(Ellipse());\n\njulia> ids  = repeat([1], length(tris));\n\njulia> Naive(tris, ids);\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.Naive-2","page":"Home","title":"PlantRayTracer.Naive","text":"Naive\n\nAllow to run the ray tracer without an acceleration structure. This should be assigned to the argument acceleration in the RayTracer function.\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.Phong-Tuple{}","page":"Home","title":"PlantRayTracer.Phong","text":"Phong(;τ = 0.0, ρd = 0.0, ρsmax = 0.0, n = 2)\n\nCreate a Phong material object from the values of transmittance (τ) diffuse reflectance (ρd), maximum Phong specular reflectance (ρsmax) and n is the specular exponent that controls the \"Phong reflectance lobe\". When more than one wavelength is being simulated, a tuple of values should be passed for each optical property (as in τ = (0.1, 0.3)).\n\nExamples\n\njulia> p = Phong(τ = 0.1, ρd = 0.2, ρsmax = 0.5);\n\njulia> p = Phong(τ = (0.1, 0.45), ρd = (0.2, 0.45), ρsmax = (0.5, 0.5));\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.PointSource","page":"Home","title":"PlantRayTracer.PointSource","text":"PointSource(vec)\n\nCreate a point irradiance source geometry at given 3D location vec, defined as vector of Cartesian coordinates (Vec(x, y, z)).\n\nExamples\n\njulia> import PlantGeomPrimitives as PG\n\njulia> source_geom = PointSource(PG.Vec(1.0, 1.0, 1.0));\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.RTSettings","page":"Home","title":"PlantRayTracer.RTSettings","text":"RTSettings(;parallel = false, pkill = 0.2, maxiter = 2, sampler = Random.Xoshiro(123456789),\n            nx = 3, ny = 3, nz = 0, dx = 0.0, dy = 0.0, dz = 0.0)\n\nSettings for the ray tracer: parallel indicates if the raytracer will run on a single core or make use of multiple cores in the machine based on Julia's multithreading support. pkill is the probably that a ray is terminated by the Russian roulette after it has been scattered a maxiter number of times. sampler is the pseudo-random number generator to be used by the ray tracer. nx and ny are the number of times the scene will be clone by the grid cloner in each direction along the x and y axis (e.g., setting nx = 1 and ny = 1 will generate a grid of 3 x 3 clones of the original scene), whereas dx and dy will be distance at which each new clone will be generated (along the axis). nz is the number of times the scene will be cloned in the vertical direction. Unlike horizontal cloning, the vertical cloning is always done in the positive direction of z axis and the number of clones will be exactly nz. See VPL documentation for more details on the ray  tracer.\n\nExamples\n\njulia> RTSettings(parallel = true, maxiter = 3);\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.RayTracer","page":"Home","title":"PlantRayTracer.RayTracer","text":"RayTracer(scene, materials, sources, settings)\n\nCreate a ray tracer object from an acceleration structure built around a 3D mesh, a grid cloner structure around the acceleration structure (scene), a vector of materials associated to the mesh (materials), a vector of sources of irradiance (sources) and settings. (as generated by RTSettings()). See VPL documentation for more details on the ray tracer.\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.RayTracer-Tuple{PlantGeomPrimitives.Scene, Any}","page":"Home","title":"PlantRayTracer.RayTracer","text":"RayTracer(scene, sources; settings = RTSettings(), acceleration = Naive, rule = nothing)\n\nCreate a RayTracer object from a scene (as generated by Scene), a tuple of sources (objects that inherit from Source), or a single source, settings and acceleration function (choose from Naive or  BVH). The argument rule is only required for the accelerator BVH and it must be an object of type SAH or AvgSplit (it is ignored for the Naive accelerator).\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> sc = Scene(mesh = Ellipse());\n\njulia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);\n\njulia> rt = RayTracer(sc, source);\n\njulia> sources = (DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000),\n                  DirectionalSource(sc, θ = 45.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000));\n\njulia> rt = RayTracer(sc, sources);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.SAH","page":"Home","title":"PlantRayTracer.SAH","text":"SAH{K}(minN, maxL)\n\nRule to be used in RayTracer when acceleration = BVH. It will divide each node at the axis and location using the Surface Area Heuristics. To speed up the construction, only K cuts will be tested per axis. These cuts will correspond to the quantiles along each axis. The rule is parameterized by the minimum number of triangles in a leaf node (minN) and the maximum depth of the tree (maxL).\n\nExamples\n\njulia> rule = SAH{3}(1,5);\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.Sensor","page":"Home","title":"PlantRayTracer.Sensor","text":"Sensor(nw::Int)\n\nCreate a sensor material object to store power for nw wavelengths. A sensor material will let rays pass through without altering the direction or irradiance. They will also not count for the total number of ray iterations.\n\nExamples\n\njulia> s = Sensor(1);\n\njulia> s = Sensor(3);\n\n\n\n\n\n","category":"type"},{"location":"#PlantRayTracer.Source-Tuple{Any, Any, Number, Integer}","page":"Home","title":"PlantRayTracer.Source","text":"Source(geom, angle, power::Number, nrays)\nSource(geom, angle, power::Tuple, nrays)\n\nCreate an irradiance source given a source geometry, a source angle, the power per ray and the total number of rays to be generated from this source. When simulating more than one wavelength simultaneously, a tuple of power values should be given, of the same length as in the materials used in the scene. See VPL documentation for details on source geometries and source angles.\n\nExamples\n\njulia> import PlantGeomPrimitives as PG\n\njulia> source_geom = PointSource(PG.O());\n\njulia> source_dir = FixedSource(0.0, 0.0);\n\njulia> Source(source_geom, source_dir, 1.0, 1_000);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.Triangle-Tuple{PlantGeomPrimitives.Mesh}","page":"Home","title":"PlantRayTracer.Triangle","text":"Triangle(mesh::Mesh)\n\nCreate a vector of ray tracing Triangle objects from a Mesh object.\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> e = Ellipse();\n\njulia> t = PlantRayTracer.Triangle(e);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.Triangle-Tuple{PlantGeomPrimitives.Scene}","page":"Home","title":"PlantRayTracer.Triangle","text":"Triangle(scene::Scene)\n\nCreate a vector of ray tracing Triangle objects from a Scene object.\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> sc = Scene(mesh = Ellipse());\n\njulia> t = PlantRayTracer.Triangle(sc);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.Triangle-Tuple{StaticArraysCore.SVector{3}, StaticArraysCore.SVector{3}, StaticArraysCore.SVector{3}}","page":"Home","title":"PlantRayTracer.Triangle","text":"Triangle(p1::Vec, p2::Vec, p3::Vec)\n\nCreate a ray tracing Triangle object given the three vertices p1, p2 and p3.\n\nExamples\n\njulia> using PlantGeomPrimitives\n\njulia> t = PlantRayTracer.Triangle(Vec(1.0, 0.0, 1.0), Vec(0.0, 1.0, .0), Vec(1.0, 1.0, 1.0));\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.Triangle-Tuple{}","page":"Home","title":"PlantRayTracer.Triangle","text":"Triangle()\n\nCreate a ray tracing Triangle object with default vertices (unit vectors in each axis).\n\nExamples\n\njulia> t = PlantRayTracer.Triangle();\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.DirectionalSource-Tuple{PlantRayTracer.AABB}","page":"Home","title":"PlantRayTracer.DirectionalSource","text":"DirectionalSource(box::AABB; θ, Φ, radiosity, nrays)\nDirectionalSource(scene::Scene; θ, Φ, radiosity, nrays)\n\nCreate a Directional source (including geometry and angle components) by providing an axis-aligned bounding box (box) or an Scene object (scene) as well as the zenith (θ) and azimuth (Φ) angles, the radiosity of the source and the number of rays to be generated. Directional sources may generate incorrect results in the absence of a grid cloner that extendes the scenes. This is because the rays are generated from the upper face of the scene's bounding box. See VPL documentation for details on light sources.\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> sc = Scene(mesh = Ellipse());\n\njulia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.accelerate-Tuple{PlantGeomPrimitives.Scene}","page":"Home","title":"PlantRayTracer.accelerate","text":"accelerate(scene::Scene; settings = RTSettings(), acceleration = Naive, rule = nothing)\n\nCreate an AccScene object from a scene, settings and acceleration function (choose from Naive or  BVH). The argument rule is only required for the accelerator BVH and it must be an object of type SAH or AvgSplit (it is ignored for the Naive accelerator). The AccScene object contains the acceleration structure and the grid cloner structure built on top of the original 3D meshes in scene. See VPL documentation for details.\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> sc = Scene(mesh = Ellipse());\n\njulia> acc = accelerate(sc);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.get_nw-Tuple{RayTracer}","page":"Home","title":"PlantRayTracer.get_nw","text":"get_nw(rt::RayTracer)\n\nRetrieve the number of wavelengths being simulated by the ray tracer. See VPL documentation for more details on the ray tracer.\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> sc = Scene(mesh = Ellipse());\n\njulia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);\n\njulia> rt = RayTracer(sc, source);\n\njulia> get_nw(rt)\n1\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.get_nw-Union{Tuple{Source{G, A, nw}}, Tuple{nw}, Tuple{A}, Tuple{G}} where {G, A, nw}","page":"Home","title":"PlantRayTracer.get_nw","text":"get_nw(s::Source)\n\nRetrieve the number of wavelengths that rays from a source will contain.\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> sc = Scene(mesh = Ellipse());\n\njulia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);\n\njulia> get_nw(source);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.power-Tuple{PlantGeomPrimitives.Material}","page":"Home","title":"PlantRayTracer.power","text":"power(material::Material)\n\nExtract the power stored inside a material.\n\nExamples\n\njulia> l = Lambertian(τ = 0.1, ρ = 0.2);\n\njulia> power(l);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.reset!-Tuple{PlantGeomPrimitives.Material}","page":"Home","title":"PlantRayTracer.reset!","text":"reset!(material::Material)\n\nReset the power stored inside a material back to zero\n\nExamples\n\njulia> l = Lambertian(τ = 0.1, ρ = 0.2);\n\njulia> l.power[1] = 10.0;\n\njulia> reset!(l);\n\njulia> l;\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.rho-Tuple","page":"Home","title":"PlantRayTracer.rho","text":"rho(vals...)\n\nGenerate values of reflectivity to be used in material object. vals... is a list of one or more comma separted values, corresponding to the different wavelengths/wavebands to be simulated in a ray tracer.\n\nExamples\n\njulia> rho(1.0, 0.0, 2.0);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.tau-Tuple","page":"Home","title":"PlantRayTracer.tau","text":"tau(vals...)\n\nGenerate values of transmisivity to be used in material object. vals... is a list of one or more comma separted values, corresponding to the different wavelengths/wavebands to be simulated in a ray tracer.\n\nExamples\n\njulia> tau(1.0, 0.0, 2.0);\n\n\n\n\n\n","category":"method"},{"location":"#PlantRayTracer.trace!-Tuple{RayTracer}","page":"Home","title":"PlantRayTracer.trace!","text":"trace!(rt)\n\nRun the ray tracing simulations. This function will overwrite the power component of any material object that is included in the scene. It returns the total number of rays being traced (primary and secondary). See VPL documentation for more details on the ray tracer.\n\nExamples\n\njulia> using PlantGeomPrimitives;\n\njulia> sc = Scene(mesh = Ellipse());\n\njulia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);\n\njulia> rt = RayTracer(sc, source);\n\njulia> trace!(rt);\n\n\n\n\n\n","category":"method"}]
}
