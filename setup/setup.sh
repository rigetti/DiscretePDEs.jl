#!/bin/bash
: '
This script runs a series of other scripts necessary for setting up your
system for development with DiscretePDEs.jl. Specifically, it

1. Sets up gmsh using setup_gmsh.sh (needed for gmsh),
2. Sets up a conda environment using setup_conda_env.sh (needed for gdspy and
   pycall).

See these respective files for more information on what they are doing.

TODO: Failures at any stage in this setup process do not terminate the setup.
Adding this timely exiting is captured in this ticket:
https://gitlab.com/rigetti/cms/coulomb.jl/issues/3
'
set -e
pushd "$(dirname "$0")"

echo "Setting up gmsh..."
source setup_gmsh.sh
echo "Setting up Julia..."
source setup_julia.sh
echo "Setting up conda..."
source setup_conda_env.sh
