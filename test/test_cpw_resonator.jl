using Test, DiscreteExteriorCalculus, DiscretePDEs
using UniqueVectors: UniqueVector
using AdmittanceModels: lossless_modes_dense, apply_transform

geo_code = GeoCode(gdspy_open(joinpath(@__DIR__, "cpw_resonator.gds"), "main"),
    layer_translation=Dict(23 => 50),
    layer_extrusion=Dict(21 => 10, 22 => 10, 23 => -100, 24 => 50),
    layer_volume_difference=[(24, 22, true), (23, 24, true), (23, 22, true)],
    layer_volume_groups=Dict("Superconductor" => [22], "Substrate" => [23], "Vacuum" => [24]),
    layer_mesh_size=Dict(21 => 10))
file_name = joinpath(@__DIR__, "cpw_resonator.geo")
geo_write!(file_name, geo_code, characteristic_length_factor=1.5, random_factor=1e-8,
    random_factor_3D=1e-8)

initialize!()
gmsh_open!(file_name)
N, K = 3, 4
mesh!(K)

node_tags, points, tcomp = get_triangulated_complex(N, K, 1e-6)
group_dict = get_physical_groups(node_tags, points)
comp = tcomp.complex
orient!(comp)
m = Metric(N)
mesh = Mesh(tcomp, circumcenter(m))

μ⁻_form = get_material(comp, 1/DiscretePDEs.μ₀, 3)
Λ⁻_form = get_material(comp, 0, 2)
σ_form = get_material(comp, 0, 2)
ϵr = 11.9 # silicon
ϵ_form = DiscretePDEs.ϵ₀ * (get_material(comp, 1, 2) +
    get_material(comp, group_dict["Substrate"], ϵr-1, 2))

pso, null_basis = electrodynamics_pso(m, mesh, Vector{Cell{N}}[],
    append!(subcomplex(comp, group_dict["Superconductor"]), boundary(comp)),
    μ⁻_form, Λ⁻_form, σ_form, ϵ_form)
constrained_pso = apply_transform(pso, null_basis)

λs, vs = lossless_modes_dense(constrained_pso, min_freq=1e9)
freqs = imag.(λs)/(2π)

hfss_answers = [47.0946, 140.140] * 1e9
@testset "cpw resonator" begin
    @test isapprox(freqs[1:2], hfss_answers, rtol=1.5e-1)
end

if false # set to true to display plots with Gmsh
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    for i in 1:2
        v = null_basis * vs[:,i]
        v /= maximum(abs.(v))
        vec_A = sharp(m, comp, v)
        add_field!("Vector potential mode $i", node_tags, vec_A[ordering])
    end
    gui!()
end
rm(file_name)
