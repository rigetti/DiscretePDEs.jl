module DiscretePDEs
    include("gmsh_interface.jl")
    include("geo.jl")
    include("gds.jl")
    include("constraints.jl")
    include("electromagnetics.jl")
end
