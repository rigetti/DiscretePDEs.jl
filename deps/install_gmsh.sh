#!/bin/bash
: '
This script installs Gmsh for either Mac OSX or Linux. On Mac, Gmsh is installed with brew.
On Linux, it is made from source. The Linux build script expects Ubuntu:18.10.
'
pushd "$(dirname "$0")" >&2

echo "Installing gmsh."
machine=$(bash determine_os.sh)
if [ "$machine" = "Mac" ]; then
  brew install gmsh >&2
elif [ "$machine" = "Linux" ]; then
  source build_gmsh.sh >&2
fi

popd >&2
