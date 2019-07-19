using Test, DiscretePDEs
const DPE = DiscretePDEs
using PyCall: pyimport

@testset "remove_seams" begin
    points = [0,1,2,3,4,5,6,5,7,3,2,8,9]
    seam_end_inds = DPE.identify_seam_end_inds(points)
    @test seam_end_inds == [8,10,11]
    begin_inds, end_inds, recur_ind = DPE.find_latest_ending_seam(points, seam_end_inds)
    @test begin_inds == [3,4]
    @test end_inds == [11,10]
    @test recur_ind == 1
    outside = [points[1:begin_inds[1]]; points[(end_inds[1]+1):end]]
    inside = points[begin_inds[end]:(end_inds[end]-1)]
    seam_end_inds[1:recur_ind] .- (begin_inds[end]-1)
    @test DPE.remove_seams(points) == [[0,1,2,8,9], [3,4,5,7], [5,6]]
end


@testset "demarcate_self_intersections" begin
    points = [(1.65, -0.34), (0.79, -0.34), (0.79, -0.53), (1.38, -0.53),
             (1.38, -0.67), (1.22, -0.67), (1.22, -0.71), (1.48, -0.71), (1.48, -0.76), (1.22, -0.76), (1.22, -0.67),
             (1.1, -0.67), (1.1, -0.8), (1.38, -0.8), (1.38, -0.92), (0.79, -0.92), (0.79, -0.34), (0.55, -0.34), (0.55, -0.71),
             (0.67, -0.71), (0.67, -0.88), (0.55, -0.88), (0.55, -0.34), (0.48, -0.34), (0.48, -1.11), (1.65, -1.11)]
    demarcated_points = [points[1:10]; points[7]; points[11:16]; points[3]; points[17:22]; points[19]; points[23:end]]
    @test DPE.demarcate_self_intersections(points) == demarcated_points
end

@testset "GeoCode" begin
    gdspy = pyimport("gdspy")
    cell = gdspy.Cell("test", exclude_from_current=true)
    cell.add(gdspy.Rectangle((0,0), (4,8)))
    geo_code = GeoCode(cell, layer_extrusion=Dict(0 => 1),
        layer_surface_groups = Dict("surface" => [0]),
        layer_volume_groups = Dict("volume" => [0]))
    @test length(geo_code.points) == 4 # subtract 1 for leading comment in each section
    @test length(geo_code.lines) == 4
    @test length(geo_code.curves) == 1
    @test length(geo_code.planes) == 1
    @test length(geo_code.surface_differences) == 0
    @test length(geo_code.extrusions) == 1
    @test length(geo_code.volume_differences) == 0
    @test length(geo_code.surface_groups) == 1
    @test length(geo_code.volume_groups) == 1
end
