using LinearAlgebra: det
using PyCall: pyimport, PyObject

# Convert a GDS file to an OpenCASCADE geo file using the python package gdspy

# The main complication is that GDS files can have polygons that contain holes,
# and in this case there are seams, i.e. repeated segments.
# These seams are not allowed in OpenCASCADE, so we must remove them
# and replace them with subtraction commands in OpenCASCADE.
# The functions below find these seams.

# To illustrate the following functions, here is an example.
# If points = [0,1,2,3,4,5,6,5,7,3,2,8,9] then the correct answer
# for remove_seams(points) is [[0,1,2,8,9], [3,4,5,7], [5,6]].
# The idea for the algorithm is that if you find a seam then the problem splits
# into two independent pieces which can be treated recursively.

# first find all indices that are at the second half of a seam, in linear time.
# in our example, this is seam_end_inds = [8,10,11]
# corresponding to the points [5, 3, 2].
function identify_seam_end_inds(points::AbstractVector{T}) where T
    point_set = Set{T}()
    seam_end_inds = Int[]
    for (i, p) in enumerate(points)
        if p in point_set
            push!(seam_end_inds, i)
        else
            push!(point_set, p)
        end
    end
    return seam_end_inds
end

# identify the latest ending seam. in our example, the first time
# this is begin_inds = [3,4], end_inds = [11,10], ind = 1.
function find_latest_ending_seam(points::AbstractVector{T}, seam_end_inds::AbstractVector{Int}) where T
    begin_inds, end_inds = Int[], Int[]
    for ind in reverse(1:length(seam_end_inds))
        j = seam_end_inds[ind]
        if (length(end_inds) != 0) && (j != end_inds[end]-1)
            return begin_inds, end_inds, ind
        end
        i = findfirst(isequal(points[j]), points)
        if (length(begin_inds) != 0) && (i != begin_inds[end]+1)
            return begin_inds, end_inds, ind
        end
        push!(begin_inds, i); push!(end_inds, j)
    end
    return begin_inds, end_inds, 0
end

# this method is a helper that performs the recursion
function remove_seams(points::AbstractVector{T}, seam_end_inds::AbstractVector{Int}) where T
    if length(seam_end_inds) == 0
        return [points]
    else
        begin_inds, end_inds, ind = find_latest_ending_seam(points, seam_end_inds)
        outside = [points[1:begin_inds[1]]; points[(end_inds[1]+1):end]]
        inside = points[begin_inds[end]:(end_inds[end]-1)]
        outside_seam_end_inds = [filter(x -> x <= begin_inds[1], seam_end_inds[1:ind]);
                                 filter(x -> end_inds[1] < x, seam_end_inds[1:ind]) .- (end_inds[1] - begin_inds[1])]
        inside_seam_end_inds = filter(x -> begin_inds[end] <= x < end_inds[end], seam_end_inds[1:ind]) .- (begin_inds[end]-1)
        return [remove_seams(outside, outside_seam_end_inds); remove_seams(inside, inside_seam_end_inds)]
    end
end

# this method puts all the pieces together
remove_seams(points::AbstractVector{T}) where T = remove_seams(points, identify_seam_end_inds(points))

# one special case is that when a polygon intersects itself at a vertex that
# vertex is not necessarily repeated, and this causes means that remove_seams
# will not correctly identify the seams. in this function, we add these missing points.
function demarcate_self_intersections(points_array::AbstractVector{<:Tuple{<:Real, <:Real}})
    output = typeof(points_array)()
    for i1 in 1:length(points_array)
        i2 = mod1(i1+1, length(points_array))
        p1, p2 = points_array[i1], points_array[i2]
        # find all points on the line from p1 to p2 using determinant
        collinears = filter(p -> det([[p1..., 1] [p2..., 1] [p..., 1]]) == 0, points_array)
        # p = t * p1 + (1-t) * p2 implies t = (p - p2) ⋅ (p1 - p2)/|p1 - p2|^2
        p_diff = collect(p1) - collect(p2)
        p_norm_square = transpose(p_diff) * p_diff
        fractional_dists = [transpose(p_diff) * (collect(p) - collect(p2))/p_norm_square for p in collinears]
        # remove anything on the line but not strictly between p1 and p2
        between_inds = 0 .< fractional_dists .< 1
        collinears = collinears[between_inds]
        fractional_dists = fractional_dists[between_inds]
        # add all the points in order from p1 to p2, not including p2
        append!(output, [p1; collinears[sortperm(fractional_dists)]])
    end
    return output
end

# find the outermost polygon, which is the one containing an extremal point,
# and move it to the end of the list
function identify_outer_polygon!(points_array::AbstractVector{<:AbstractVector{<:Tuple{<:Real, <:Real}}})
    all_points = vcat(points_array...)
    extremal_point = all_points[argmax([p[1] for p in all_points])]
    ind = findfirst(arr -> extremal_point in arr, points_array)
    extremal_arr = points_array[ind]
    return push!(deleteat!(points_array, ind), extremal_arr)
end

function GeoCode(cell::PyObject;
    layer_translation::Dict{Int, <:Real}=Dict{Int, Float64}(),
    layer_extrusion::Dict{Int, <:Real}=Dict{Int, Float64}(),
    layer_volume_difference::AbstractVector{Tuple{Int, Int, Bool}}=Tuple{Int, Int, Bool}[],
    layer_surface_groups::Dict{String, <:AbstractVector{Int}}=Dict{String, Vector{Int}}(),
    layer_volume_groups::Dict{String, <:AbstractVector{Int}}=Dict{String, Vector{Int}}(),
    layer_mesh_size::Dict{Int, <:Real}=Dict{Int, Float64}())
    elements = cell.copy("", exclude_from_current=true).flatten().elements # flattening ensures all are PolygonSet objects
    geo_code = GeoCode()
    surface_group_layers = unique(vcat(values(layer_surface_groups)...))
    volume_group_layers = unique(vcat(values(layer_volume_groups)...))
    @assert issubset(volume_group_layers, keys(layer_extrusion))
    layer_to_planes = Dict(l => GeoPlane[] for l in surface_group_layers)
    layer_to_extrusions = Dict(l => GeoSurfaceExtrusion[] for l in keys(layer_extrusion))
    for polygon_set in elements
        for (layer, points_array) in zip(polygon_set.layers, polygon_set.polygons)
            layer = convert(Int, layer) # convert from PyObject
            point_tuples = [(points_array[i,1], points_array[i,2]) for i in 1:size(points_array,1)]
            # find seams and remove them, then create inner polygons and outer polygon
            seamless_polygons = filter(p -> length(p) > 0, remove_seams(demarcate_self_intersections(point_tuples)))
            identify_outer_polygon!(seamless_polygons)
            planes = GeoPlane[] # will need to subtract all planes 1:end-1 from the last
            z = get(layer_translation, layer, 0.0)
            for points in seamless_polygons
                mesh_size = get(layer_mesh_size, layer, 0.0)
                g = GeoCode([(1.0 * p[1], 1.0 * p[2], 1.0 * z, 1.0 * mesh_size)
                    for p in points]) # p is just an x and y coordinate
                append!(geo_code, g)
                push!(planes, g.planes[1])
            end
            # subtract inner polygons from outer polygon
            for p in planes[1:end-1]
                push!(geo_code.surface_differences,
                    GeoSurfaceDifference(planes[end], p, false))
            end
            # add last plane to a physical group
            if layer in keys(layer_to_planes)
                push!(layer_to_planes[layer], planes[end])
            end
            # put extrusion code in ending_code since it creates points
            if layer in keys(layer_extrusion)
                extrusion = GeoSurfaceExtrusion(0, 0, layer_extrusion[layer], planes[end])
                push!(geo_code.extrusions, extrusion)
                push!(layer_to_extrusions[layer], extrusion)
            end
        end
    end
    # add any volume differences
    for (l1, l2, keep_tool) in layer_volume_difference
        for v1 in layer_to_extrusions[l1]
            for v2 in layer_to_extrusions[l2]
                push!(geo_code.volume_differences,
                    GeoVolumeDifference(v1, v2, keep_tool))
            end
        end
    end
    # add any physical groups
    for name in keys(layer_surface_groups)
        planes = vcat([layer_to_planes[l] for l in layer_surface_groups[name]]...)
        push!(geo_code.surface_groups, GeoSurfaceGroup(name, planes))
    end
    for name in keys(layer_volume_groups)
        extrusions = vcat([layer_to_extrusions[l] for l in layer_volume_groups[name]]...)
        push!(geo_code.volume_groups, GeoVolumeGroup(name, extrusions))
    end
    return geo_code
end

function gdspy_open(file_name::AbstractString, cell_name::AbstractString)
    gds = pyimport("gdspy").GdsLibrary().read_gds(file_name)
    return gds.cell_dict[cell_name]
end
