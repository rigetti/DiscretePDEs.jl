using Test, DiscreteExteriorCalculus, DiscretePDEs
using UniqueVectors: UniqueVector

radius, height, offset1, offset2, height_fraction = 3, 2, 1, 1.5, .1
file_name = joinpath(@__DIR__, "current_carrying_wire_cylinder.geo")
geo_write!(file_name, characteristic_length_factor=.3,
    footer="""
    Cylinder(1) = {0, 0, 0, 0, 0, $height, $radius, 2*Pi};
    Point(3) = {0, 0, 0, 1.0};
    Point(4) = {0, 0, $height, 1.0};
    Point(5) = {$offset1, 0, $(height_fraction * height), 1.0};
    Point(6) = {$offset1, 0, $((1-height_fraction) * height), 1.0};
    Point(7) = {$offset2, 0, $(height_fraction * height), 1.0};
    Point(8) = {$offset2, 0, $((1-height_fraction) * height), 1.0};
    Point(9) = {$(-offset1), 0, 0, 1.0};
    Point(10) = {$(-offset1), 0, $height, 1.0};
    Line(4) = {3, 4};
    Line(5) = {5, 6};
    Line(6) = {6, 8};
    Line(7) = {8, 7};
    Line(8) = {7, 5};
    Line(9) = {9, 10};
    Point{3, 9} In Surface{3};
    Point{4, 10} In Surface{2};
    Point{6, 8, 7, 5} In Volume{1};
    Curve{4, 5, 6, 7, 8, 9} In Volume{1};
    Physical Curve("Wire center") = {4};
    Physical Curve("Loop") = {6, 7, 8, 5};
    Physical Curve("Wire offset") = {9};
    """)
initialize!()
gmsh_open!(file_name)
N, K = 3, 4
mesh!(K)

node_tags, points, tcomp = get_triangulated_complex(N, K)
group_dict = get_physical_groups(node_tags, points)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
orient!(comp)
sources = [get_current_source(comp, group_dict[k]) for k in
    ["Wire center", "Loop", "Wire offset"]]
m = Metric(N)
mesh = Mesh(tcomp, circumcenter(m))

μ⁻, Λ⁻ = 2, 0
μ⁻_form = get_material(comp, μ⁻, 3)
Λ⁻_form = get_material(comp, Λ⁻, 2)

bbox, null_basis = magnetostatics_blackbox(m, mesh, sources, boundary(comp),
    μ⁻_form, Λ⁻_form)
inductance = impedance_matrix(bbox, null_basis)
@testset "mutual inductance of wires and loops" begin
    @test isapprox(abs(inductance[1,3]),
        abs((log(radius) - log(offset1)) * height/(2π*μ⁻)), rtol=3e-2)
    @test isapprox(abs(inductance[1,2]),
        abs((log(offset1) - log(offset2)) * height * (1-2*height_fraction)/(2π*μ⁻)), rtol=3e-2)
end

if false # set to true to display plots with Gmsh
    A, _ = solve_statics(bbox, null_basis, [0.0, 1.0, 0.0])
    vec_A = sharp(m, comp, A)
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    add_field!("Vector potential", node_tags, vec_A[ordering])
    gui!()
end
rm(file_name)
