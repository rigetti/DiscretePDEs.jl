using Test, DiscreteExteriorCalculus, DiscretePDEs
using UniqueVectors: UniqueVector
using AdmittanceModels: lossless_modes_dense, apply_transform

geo_code = GeoCode(gdspy_open(joinpath(@__DIR__, "boxes.gds"), "main"),
    layer_translation=Dict(12 => -50, 14 => -20),
    layer_extrusion=Dict(12 => 100, 14 => 40),
    layer_volume_difference=[(12, 14, false)])
file_name = joinpath(@__DIR__, "boxes.geo")
geo_write!(file_name, geo_code, characteristic_length_factor=.9)

initialize!()
gmsh_open!(file_name)
N, K = 3, 4
mesh!(K)
node_tags, points, tcomp = get_triangulated_complex(N, K)
group_dict = get_physical_groups(node_tags, points)
comp = tcomp.complex
orient!(comp)
m = Metric(N)
mesh = Mesh(tcomp, circumcenter(m))

μ⁻, Λ⁻, σ, ϵ = 1/DiscretePDEs.μ₀, 0, 0, DiscretePDEs.ϵ₀
μ⁻_form = get_material(comp, μ⁻, 3)
Λ⁻_form = get_material(comp, Λ⁻, 2)
σ_form = get_material(comp, σ, 2)
ϵ_form = get_material(comp, ϵ, 2)

pso, null_basis = electrodynamics_pso(m, mesh, Vector{Cell{N}}[], boundary(comp),
    μ⁻_form, Λ⁻_form, σ_form, ϵ_form)
constrained_pso = apply_transform(pso, null_basis)
λs, vs = lossless_modes_dense(constrained_pso, min_freq=1e6)
freq = imag.(λs[1])/(2π)

hfss_answer = 1.97194e6
@testset "modes of box with superconductors" begin
    @test isapprox(freq, hfss_answer, rtol=6e-2)
end

if false # set to true to display plots with Gmsh
    v = null_basis * vs[:,1]
    v /= maximum(abs.(v))

    vec_A = sharp(m, comp, v)
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    add_field!("Vector potential", node_tags, vec_A[ordering])
    gui!()
end
rm(file_name)
