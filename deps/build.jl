using DotEnv, Pkg

# Ensure ENV["PYTHON"] is set for PyCall. Should match the python env created
# using conda-env.yml.
DotEnv.config(path=joinpath(@__DIR__, "../.env"))
Pkg.build("PyCall")
