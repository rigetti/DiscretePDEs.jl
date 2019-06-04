using Test, DiscreteExteriorCalculus, DiscretePDEs
const DPE = DiscretePDEs
using LinearAlgebra: norm, normalize
using UniqueVectors: UniqueVector
using Statistics: median

function disks!(r1::Real, r2::Real, factor::Real, N::Int)
    @assert (N == 2) || (N == 3)
    file_name = joinpath(@__DIR__, "disks.geo")
    if N == 2
        DPE.geo_write!(file_name, characteristic_length_factor=factor,
            footer="""
            Circle(1) = {0, 0, 0, $r1, 0, 2*Pi};
            Circle(2) = {0, 0, 0, $r2, 0, 2*Pi};
            Curve Loop(1) = {2};
            Curve Loop(2) = {1};
            Plane Surface(1) = {1, 2};
            Physical Curve("Conductor") = {1};
            Physical Surface("Dielectric") = {1};
            """)
    elseif N == 3
        DPE.geo_write!(file_name, characteristic_length_factor=factor,
            footer="""
            Sphere(1) = {0, 0, 0, $r1, -Pi/2, Pi/2, 2*Pi};
            Sphere(2) = {0, 0, 0, $r2, -Pi/2, Pi/2, 2*Pi};
            BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
            Physical Surface("Conductor") = {1};
            Physical Volume("Dielectric") = {2};
            """)
    end
    return file_name
end

function voltage(R, outer_radius, ϵ, N)
    @assert (N == 2) || (N == 3)
    if N == 2
        return (log(outer_radius) - log(norm(R)))/(2π * ϵ)
    elseif N == 3
        return (1/norm(R) - 1/outer_radius)/(4π * ϵ)
    end
end

DPE.initialize!()
for (N, factor) in [(2, .3), (3, .45)]
    r1, r2 = 1.5, 3
    file_name = disks!(r1, r2, factor, N)
    DPE.gmsh_open!(file_name)
    K = N+1
    DPE.mesh!(K)

    node_tags, points, tcomp = DPE.get_triangulated_complex(N, K)
    group_dict = DPE.get_physical_groups(node_tags, points)
    @test typeof(tcomp) <: TriangulatedComplex{N, K}
    comp = tcomp.complex
    orient!(comp)
    sources = [DPE.get_charge_source(comp, group_dict["Conductor"])]
    _, exterior = boundary_components_connected(comp)
    m = Metric(N)
    mesh = Mesh(tcomp, circumcenter(m))

    ϵ = 2
    ϵ_form = DPE.get_material(comp, group_dict["Dielectric"], ϵ, 2)
    bbox, null_basis = DPE.electrostatics_blackbox(m, mesh, sources, exterior, ϵ_form)

    elastance = DPE.impedance_matrix(bbox, null_basis)[1,1]
    charge = 2
    φ, q = DPE.solve_statics(bbox, null_basis, [charge])
    @testset "concentric $N-disks" begin
        @test isapprox(elastance[1,1], voltage(r1, r2, ϵ, N), rtol= N == 2 ? 1e-4 : 4e-2)
        @test bbox.Y[1] * φ == q
        inner_charges = q[[findfirst(isequal(c), comp.cells[1]) for c in sources[1]]]
        @test sum(inner_charges) ≈ charge
        @test all(abs.(diff(inner_charges)) .< charge * (N == 2 ? .0003 : .005))
        outer_charges = q[[findfirst(isequal(c), comp.cells[1]) for c in exterior.cells[1]]]
        @test sum(outer_charges) ≈ -charge
        @test all(abs.(diff(outer_charges)) .< charge * (N == 2 ? .0003 : .005))
        @test isapprox(φ, charge * voltage.([c.points[1].coords
            for c in comp.cells[1]], r2, ϵ, N), rtol= N == 2 ? 7e-4 : 4e-2)
    end
    if false # make plots with Gmsh
        ρ = DPE.source_density(m, mesh, q, 1)
        vec_E = [N == 2 ? [p..., 0] : p for p in DPE.φ_to_vec_E(m, mesh, φ)]
        comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
        ordering = [findfirst(isequal(p), comp_points) for p in points]
        DPE.add_field!("Electric field", node_tags, vec_E[ordering])
        inds = N == 2 ? (1:length(points)) : findall([p.coords[1] > 0 for p in points])
        DPE.add_field!("Electric field magnitude", node_tags[inds], norm.(vec_E)[ordering][inds])
        DPE.add_field!("Voltage", node_tags[inds], φ[ordering][inds])
        DPE.add_field!("Charge density", node_tags[inds], ρ[ordering][inds])
        DPE.gui!()
    end
    rm(file_name)
end
