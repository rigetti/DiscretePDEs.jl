using DotEnv, Pkg

"""
See README.md for build notes.
"""

env = joinpath(@__DIR__, ".env")
run(`touch $env`) # make sure .env file exists
DotEnv.config(path=env)
env_keys = ["gmshjl", "PYTHON"]

# configure gmsh
find_gmshjl() =
    try # search for gmsh.jl in /usr/local/
        path = read(pipeline(`find /usr/local/ -name "gmsh.jl"`, `head -n 1`), String)
        path[1:end-1]
    catch
        ""
    end
key = env_keys[1]
if !(key in keys(ENV))
    value = find_gmshjl()
    if value == ""
        # install gmsh
        install_gmsh_sh = joinpath(@__DIR__, "install_gmsh.sh")
        run(`bash $install_gmsh_sh`)
        value = find_gmshjl()
        @assert value != ""
    end
    ENV[key] = value
    println("$key=$value")
end

# configure python
find_conda() =
    try # use conda in $PATH
        path = read(`which conda`, String)
        path[1:end-1]
    catch
        try # check if miniconda3 is already in home directory
            path = joinpath(homedir(), "miniconda3", "bin", "conda")
            run(`test -f $path`) # check if file exists
            path
        catch
            ""
        end
    end
find_python(conda::String) =
    try
        path = read(pipeline(`$conda env list`, `grep DiscretePDEs`, `awk '{print $2}'`),
            String)
        joinpath(path[1:end-1], "bin", "python")
    catch
        ""
    end
key = env_keys[2]
if !(key in keys(ENV))
    conda = find_conda()
    if conda == ""
        # install conda
        install_conda_sh = joinpath(@__DIR__, "install_conda.sh")
        run(`bash $install_conda_sh`)
        conda = find_conda()
        @assert conda != ""
    end
    value = find_python(conda)
    if value == ""
        # create a virtual environment
        conda_env_yml = joinpath(@__DIR__, "conda_env.yml")
        run(`$conda env create -f $conda_env_yml`)
        value = find_python(conda)
        @assert value != ""
    end
    ENV[key] = value
    println("$key=$value")
end
Pkg.build("PyCall")

# save to .env
env_contents = reduce(*, ["$k=$(ENV[k])\n" for k in env_keys])
open(env, "w") do io
    write(io, env_contents)
end
