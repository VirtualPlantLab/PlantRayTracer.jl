### This file DOES NOT contain public API ###

###############################################################################
################### Grid cloner with transformation nodes #####################
###############################################################################

# A node in the grid cloner is defined by an AABB, a leaf flag and, if a leaf, a displacement vector
struct TNODE{FT}
    box::AABB{FT}
    leaf::Bool
    disp::PGP.Vec{FT}
end

# Shortcut to create intermediate TNODES in the binary tree of the grid cloner
TNODE(box::AABB{FT}, leaf) where {FT} = TNODE(box, leaf, PGP.Vec{FT}(0, 0, 0))

# Structure that contains the binary tree required to traverse the grid cloner
# We store dx, dy and dz so that it can be used by other algorithms like sky domes
struct GridCloner{FT}
    nodes::GVector{TNODE{FT}}
    nleaves::Int
    dx::Float64
    dy::Float64
    dz::Float64
end

# Create a grid cloner around an acceleration structure
# Create (2*n + 1) clones in each direction -> symmetry enforced by odd numbers
# Create n in z direction (always go up)
function GridCloner(acc::Acceleration{FT}; nx = 3, ny = 3, nz = 0,
    dx = nothing, dy = nothing, dz = nothing) where {FT}
    # Note that nothing is a special case to signal that dx, dy and dz should
    # be derived from the bounding box
    dx == nothing && (dx = acc.gbox.max[1] - acc.gbox.min[1])
    dy == nothing && (dy = acc.gbox.max[2] - acc.gbox.min[2])
    dz == nothing && (dz = acc.gbox.max[3] - acc.gbox.min[3])
    # Special case when the global box is empty
    gbox = acc.gbox
    isempty(acc.gbox) ||
        nx == ny == nz == 0 &&
            (return GridCloner(GVector([TNODE(gbox, true, PGP.Vec{FT}(0, 0, 0))]),
                               1, 0.0, 0.0, 0.0))
    # Create all the leaf TNODES, their centers and indices
    scene = create_tnodes(nx, ny, nz, dx, dy, dz, gbox)
    indices = collect(1:length(scene.nodes))
    # Create empty grid structure with the two arrays
    grid = GridCloner(GVector(TNODE{FT}[]), length(scene.nodes), dx, dy, dz)
    # Call the recusive function to add (packets of) nodes and triangles
    parentbox = AABB(scene.boxes, indices)
    addNode!(grid, scene, parentbox, 0, indices)
    # Return the update grid object
    return grid
end

# Create all the transformation nodes in the grid
function create_tnodes(nx, ny, nz, dx::FT, dy::FT, dz::FT, box) where {FT}
    # Total number of tnodes in the grid
    nnodes = (2 * nx + 1) * (2 * ny + 1) * (nz + 1)

    # Special case when only one node in the grid cloner
    if nnodes == 1
        leaves = [TNODE(box, true, PGP.Vec{FT}(0, 0, 0))]
        centers = hcat(FT(0), FT(0), FT(0))
        return leaves, centers
    end

    # Create all leaf nodes in the trees
    nodes = Vector{TNODE}(undef, nnodes)
    boxes = Vector{AABB{FT}}(undef, nnodes)
    centers = zeros(FT, nnodes, 3)
    p0 = center(box)
    counter = 0
    for i in (-nx):1:nx
        for j in (-ny):1:ny
            for k in 0:1:nz
                counter += 1
                disp = PGP.Vec{FT}(i * dx, j * dy, k * dz)
                nbox = AABB(box.min .+ disp, box.max .+ disp)
                nodes[counter] = TNODE(nbox, true, disp)
                boxes[counter] = nbox
                centers[counter, :] = p0 .+ disp
            end
        end
    end
    return (nodes = nodes, centers = centers, boxes = boxes)
end

###############################################################################
######################### Grid cloner construction ############################
###############################################################################

# Recursive function to create a binary tree of transformation nodes for the grid cloner
function addNode!(grid::GridCloner, scene, parentbox, parentid, indices)
    # Split the current node and distribute children to each side of the split
    newindices, childrenboxes, leaves = splitnode(grid, parentbox, indices, scene)
    # Ids of the children in the flattened version of the tree.
    # First child of a node with order i is given by 2i + 1
    nodesid = Tuple(2 * parentid + i for i in 1:2)
    # Add children to nodes
    pushnodes!(grid, scene, nodesid, childrenboxes, newindices, leaves)
    # Recursion trigerred for children that are inner nodes (depth-first construction)
    c = 0
    for id in nodesid
        c += 1
        node = grid.nodes[id]
        if !node.leaf
            addNode!(grid, scene, node.box, id, newindices[c])
        end
    end
    return nothing
end

# Determine the split of a inner grid cloner node and distribute nodes to each side
# The rule is simple: Always split along the longest axis and do it in the middle
function splitnode(grid::GridCloner, parentbox, indices, scene)
    # Find the longest distance among centers along axes
    centers = scene.centers[indices, :]
    distances = (maximum(centers[:, 1]) - minimum(centers[:, 1]),
        maximum(centers[:, 2]) - minimum(centers[:, 2]),
        maximum(centers[:, 3]) - minimum(centers[:, 3]))
    index = findmax(distances)[2]
    # Split along the chosen axis (arithmetic mean is good enough for TNODEs)
    pos = Statistics.mean(centers[:, index])
    # Distribute the nodes according to position of centers relative to pos
    leftindices = indices[centers[:, index] .< pos]
    rightindices = setdiff(indices, leftindices)
    # Create bounding boxes for inner nodes
    leftbox = AABB(scene.boxes, leftindices)
    rightbox = AABB(scene.boxes, rightindices)
    # Determine if the children should be leaf nodes
    leaves = (length(leftindices) == 1, length(rightindices) == 1)
    # Return the new indices and bounding boxes
    (newindices = (leftindices, rightindices),
        childrenboxes = (leftbox, rightbox),
        leaves = leaves)
end

# Add TNODEs to the grid cloner (either pre-constructed leaf nodes or create new inner node)
function pushnodes!(grid::GridCloner, scene, nodesid, childrenboxes, newindices, leaves)
    #@inbounds
    for i in 1:2
        if leaves[i]
            grid.nodes[nodesid[i]] = scene.nodes[newindices[i]...]
        else
            grid.nodes[nodesid[i]] = TNODE(childrenboxes[i], false)
        end
    end
    return nothing
end

###############################################################################
########################### Grid cloner traversal #############################
###############################################################################

# Traverse the grid cloner (if a leaf TNODE is hit, update the ray and traverse the scene)
function Base.intersect(ray::Ray{FT}, grid::GridCloner, acc::Acceleration{FT}, tnodestack,
    tdstack, nodestack, dstack) where {FT}
    #@inbounds begin
    begin
        # Initialize solution state
        dmin = Inf
        disp = PGP.Vec{FT}(0, 0, 0) # Do we still need to return disp?
        hit = false
        intersection = Intersection(FT)
        # Special case when we only have one node (go directly to acceleration structure)
        if grid.nleaves == 1
            hit, intersection, dmin = intersect(ray, acc, nodestack, dstack, dmin)
            return hit, intersection, disp
        end
        # Check the two children of the global box to start the process
        update!(tnodestack, tdstack, 1, ray, grid, dmin)
        # Depth-first traversal using LIFO stack
        while (length(tnodestack) > 0)
            # Pop the next node and hit distance from the stack
            nodecur = pop!(tnodestack)
            node = grid.nodes[nodecur]
            dnode = pop!(tdstack)
            # If the node is a leaf, update the ray and traverse the accelerated scene
            if node.leaf && dnode < dmin
                ray_disp = Ray(ray.o .- node.disp, ray.dir)
                l_hit, l_intersection, l_dmin = intersect(ray_disp,
                    acc,
                    nodestack,
                    dstack,
                    dmin)
                # If it hits a triangle inside the node, update everything
                if l_hit && l_dmin < dmin
                    dmin = l_dmin
                    hit = true
                    intersection = l_intersection
                    disp = node.disp
                end
                # If the node is inner, compute the index of the first child and update nodestack
            else
                childid = 2 * nodecur + 1
                update!(tnodestack, tdstack, childid, ray, grid, dmin)
            end
        end
        return hit, intersection, disp
    end
end

# Check intersection with two children and store them in the nodestack
function update!(tnodestack, tdstack, i, ray, grid::GridCloner, dmin)
    # Intersect the two nodes
    hit1, d1 = intersect(ray, grid.nodes[i], dmin)
    hit2, d2 = intersect(ray, grid.nodes[i + 1], dmin)
    # If both nodes are hit, store them in the order of hit
    if hit1 && hit2
        first, second, d1, d2 = ifelse(d1 < d2, (i, i + 1, d1, d2), (i + 1, i, d2, d1))
        push!(tnodestack, first)
        push!(tnodestack, second)
        # TODO: Avoid redundant operation on d1 and d2 when sorting?
        push!(tdstack, d1)
        push!(tdstack, d2)
        # If only one (or none) of the nodes are hit
    else
        hit1 && push!(tnodestack, i)
        hit2 && push!(tnodestack, i + 1)
        hit1 && push!(tdstack, d1)
        hit2 && push!(tdstack, d2)
    end
    return nothing
end

# Intersect a node in the tree and compare against dmin
function Base.intersect(ray::Ray, node::TNODE, dmin)
    hit, tmin = intersect(ray, node.box)
    hit && tmin < dmin, tmin
end
