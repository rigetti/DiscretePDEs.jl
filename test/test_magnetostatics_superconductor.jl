using Test, DiscreteExteriorCalculus, DiscretePDEs
const DEC = DiscreteExteriorCalculus
const DPE = DiscretePDEs
using LinearAlgebra: norm, normalize
using UniqueVectors: UniqueVector

file_name = joinpath(@__DIR__, "superconductor.geo")
DPE.geo_write!(file_name, characteristic_length_factor=.7,
    footer="""
    Box(1) = {-5, -5, -1, 10, 10, 2};
    Box(2) = {-10, -10, -5, 20, 20, 10};
    BooleanDifference{Volume{2}; Delete;}{Volume{1};}
    Point(17) = {10, 0, 0, 1.0};
    Point(18) = {5, 0, 0, 1.0};
    Point(19) = {-5, 0, 0, 1.0};
    Point(20) = {-10, 0, 0, 1.0};
    Line(25) = {17, 18};
    Line(26) = {19, 20};
    Point{17} In Surface{12};
    Point{18} In Surface{2};
    Point{19} In Surface{1};
    Point{20} In Surface{7};
    Curve{25} In Volume{2};
    Curve{26} In Volume{2};
    Physical Curve("1") = {25};
    Physical Curve("2") = {26};
    Physical Volume("Superconductor") = {1};
    """)
DPE.initialize!()
DPE.gmsh_open!(file_name)
N, K = 3, 4
DPE.mesh!(K)

node_tags, points, tcomp = DPE.get_triangulated_complex(N, K)
group_dict = DPE.get_physical_groups(node_tags, points)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
DEC.orient!(comp)
sources = [DPE.get_current_source(comp, group_dict[k])
    for k in ["1", "2"]]
boundary = DEC.boundary(comp)
m = Metric(N)
mesh = Mesh(tcomp, DEC.circumcenter(m))

μ⁻, Λ⁻ = 2, 0
μ⁻_form = DPE.get_material(comp, μ⁻, 3)
Λ⁻_form = DPE.get_material(comp, Λ⁻, 2)

b = CellComplex{N, K}()
append!(b, boundary)
append!(b, DEC.submanifold(comp, group_dict["Superconductor"]))
bbox, null_basis = DPE.magnetostatics_blackbox(m, mesh,
    sources, b, μ⁻_form, Λ⁻_form)

if false
    A, I = DPE.solve_statics(bbox, null_basis, [1.0, 1.0])
    vec_A = DEC.sharp(m, comp, A)
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    DPE.add_field!("Vector potential", node_tags, vec_A[ordering])
    J = DPE.source_density(m, mesh, I, 2)
    vec_J = DEC.sharp(m, comp, J)
    # current on superconductor
    sub_points = UniqueVector([c.points[1] for c in DEC.submanifold(comp, group_dict["Superconductor"]).cells[1]])
    sub_inds = findall([p in sub_points for p in points])
    DPE.add_field!("Current density superconductor", node_tags[sub_inds], vec_J[ordering[sub_inds]])
    # current on boundary
    sub_points = UniqueVector([c.points[1] for c in boundary.cells[1]])
    sub_inds = findall([p in sub_points for p in points])
    DPE.add_field!("Current density boundary", node_tags[sub_inds], vec_J[ordering[sub_inds]])
    DPE.gui!()
end
rm(file_name)
