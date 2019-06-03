using Test, DiscreteExteriorCalculus, DiscretePDEs
const DEC = DiscreteExteriorCalculus
const DPE = DiscretePDEs
using UniqueVectors: UniqueVector
using AdmittanceModels: lossless_modes_dense, apply_transform

geo_code = DPE.GeoCode(DPE.gdspy_open(joinpath(@__DIR__, "boxes.gds"), "main"),
    layer_translation=Dict(12 => -50, 14 => -20),
    layer_extrusion=Dict(12 => 100, 14 => 40),
    layer_volume_difference=[(12, 14, false)])
file_name = joinpath(@__DIR__, "boxes.geo")
DPE.geo_write!(file_name, geo_code, characteristic_length_factor=.9)

DPE.initialize!()
DPE.gmsh_open!(file_name)
N, K = 3, 4
DPE.mesh!(K)
node_tags, points, tcomp = DPE.get_triangulated_complex(N, K)
group_dict = DPE.get_physical_groups(node_tags, points)
comp = tcomp.complex
DEC.orient!(comp)
boundary = DEC.boundary(comp)
m = Metric(N)
mesh = Mesh(tcomp, DEC.circumcenter(m))

μ⁻, Λ⁻, σ, ϵ = 1/DPE.μ₀, 0, 0, DPE.ϵ₀
μ⁻_form = DPE.get_material(comp, μ⁻, 3)
Λ⁻_form = DPE.get_material(comp, Λ⁻, 2)
σ_form = DPE.get_material(comp, σ, 2)
ϵ_form = DPE.get_material(comp, ϵ, 2)

pso, null_basis = DPE.coulomb_pso(m, mesh, Vector{Cell{N}}[], boundary,
    μ⁻_form, Λ⁻_form, σ_form, ϵ_form)
constrained_pso = apply_transform(pso, null_basis)
λs, vs = lossless_modes_dense(constrained_pso, min_freq=1e6)
freq = imag.(λs[1])/(2π)

hfss_answer = 1.97194e6
@testset "box with boxes" begin
    @test isapprox(freq, hfss_answer, rtol=6e-2)
end

if false
    v = null_basis * vs[:,1]
    v /= maximum(abs.(v))

    vec_A = DEC.sharp(m, comp, v)
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    DPE.add_field!("Vector potential", node_tags, vec_A[ordering])
    DPE.gui!()
end
rm(file_name)
