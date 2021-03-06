using Test, DiscreteExteriorCalculus, DiscretePDEs
using UniqueVectors: UniqueVector

geo_code = GeoCode(gdspy_open(joinpath(@__DIR__, "two_pads.gds"), "main"),
    layer_translation=Dict(24 => -50),
    layer_extrusion=Dict(22 => -50, 23 => -50, 24 => 250),
    layer_volume_difference=[(23, 22, true), (24, 23, true), (24, 22, true)],
    layer_volume_groups=Dict("Substrate" => [22, 23], "Vacuum" => [24]))
file_name = joinpath(@__DIR__, "two_pads.geo")
geo_write!(file_name, geo_code, characteristic_length_factor=.35,
    footer="""
    Physical Surface("Ground") = {4};
    Physical Surface("Pad 1") = {5};
    Physical Surface("Pad 2") = {6};
    """)

initialize!()
gmsh_open!(file_name)
N, K = 3, 4
mesh!(K)

node_tags, points, tcomp = get_triangulated_complex(N, K, 1e-3)
group_dict = get_physical_groups(node_tags, points)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
orient!(comp)
sources = [get_charge_source(comp, group_dict[k])
    for k in ["Pad 1", "Pad 2"]]
exterior = boundary(comp)
m = Metric(N)
mesh = Mesh(tcomp, circumcenter(m))

ϵr = 2
ϵ_form = DiscretePDEs.ϵ₀ * (get_material(comp, 1, 2) +
    get_material(comp, group_dict["Substrate"], ϵr-1, 2))

bbox, null_basis = electrostatics_blackbox(m, mesh, sources,
    append!(subcomplex(comp, group_dict["Ground"]), boundary(comp)), ϵ_form)

capacitance = admittance_matrix(bbox, null_basis)

q3d_capacitance = [[80.995, -1.2219] [-1.2219, 80.994]] * 1e-12
@testset "capacitance two pads" begin
    @test isapprox(capacitance, q3d_capacitance, rtol=.25)
end

if false # set to true to display plots with Gmsh
    φ, q = solve_statics(bbox, null_basis, [1.0, -1.0])
    ρ = source_density(m, mesh, q, 1)
    vec_E = electric_field(m, mesh, φ)
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    add_field!("Electric field", node_tags, vec_E[ordering])
    sub_points = UniqueVector([c.points[1] for c in subcomplex(comp,
        vcat(group_dict["Pad 1"], group_dict["Pad 2"], group_dict["Ground"])).cells[1]])
    sub_inds = findall([p in sub_points for p in points])
    add_field!("Charge density", node_tags[sub_inds], ρ[ordering[sub_inds]])
    add_field!("Voltage", node_tags[sub_inds], φ[ordering[sub_inds]])
    gui!()
end
rm(file_name)
