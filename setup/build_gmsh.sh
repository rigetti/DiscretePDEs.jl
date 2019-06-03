#!/bin/bash
: '
This script configures a machine (expected to be Ubuntu:18.10) to have gmsh. The configuration
is largely taken from these instructions:
https://gitlab.onelab.info/gmsh/gmsh/wikis/Gmsh-compilation#opencascade.
'

pushd "$(dirname "$0")"

echo "*****************************************************"
echo "Installing gmsh from source. This may take a while..."
echo "*****************************************************"
echo ""

set -e

# Install dependencies
apt-get update
# Need noninteractive here so that tcl/tk tzdata dependency doesn't prompt user
# for timezone input
# TODO: Determine what from this list can be removed.
DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential cmake curl gfortran git libxi-dev libxmu-dev libxmu-headers wget tk-dev tcl-dev libgpm2 libmagic-mgc libmagic1 vim

# Install OpenCascade
apt-get install -y libocct-foundation-dev libocct-data-exchange-dev

# Make gmsh from source
apt-get install -y libfltk1.3-dev
git clone http://gitlab.onelab.info/gmsh/gmsh.git
pushd gmsh

# pin to v4.1.1 tag
git checkout gmsh_4_3_0

mkdir build
pushd build
cmake -DENABLE_BUILD_DYNAMIC=1 -DENABLE_FLTK=1 -DENABLE_MPI=0 -DCMAKE_PREFIX_PATH=/usr/local/include/opencascade ..
make -j5
make install
popd
popd
