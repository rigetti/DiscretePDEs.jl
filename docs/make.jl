using Documenter, DiscretePDEs

makedocs(modules = [DiscretePDEs],
    sitename="DiscretePDEs.jl",
    pages = ["Main" => "index.md"])

deploydocs(
    repo = "github.com/rigetti/DiscretePDEs.jl.git",
)
