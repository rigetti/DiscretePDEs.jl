# DiscretePDEs.jl

## Installation & setup

Run [setup.sh](./setup/setup.sh):

```
cd setup && bash setup.sh
```

This will setup both gmsh as well as a conda environment for gdspy (accessed via
pycall). To learn more about this, read the respective scripts that are run by
setup.sh.

Note: Before you run `setup.sh`, confirm that your machine is supported by
referring to the section below.

### Supported systems

DiscretePDEs currently supports

- Mac OSX Mojave and
- Ubuntu 18.10.

Ask the maintainers about support for another system if you need it.

### Installation on Linux (Ubuntu 18.10)

At this time (2019/04/20) the only known reliable installation method for gmsh
on Linux is to build the binaries from source. Building from source
takes time, so the setup procedure on Linux is to

1. Do this only occasionally for any new machine, or
2. Use the `DiscretePDEs.jl.base` docker image created by the `build_gmsh.sh` file.

### Installation on Mac

The setup for Mac OSX is quite a bit simpler than Linux, since gmsh is installed
via Homebrew. Please confirm directly once installed, that the command
`gmsh` works from the command line; we have seen anecdotally that installation
will run fine, but that there may be a missing dependency due to a brew linking
error.

### Building gmsh binaries from source

See [build_gmsh.sh](./setup/build_gmsh.sh).
