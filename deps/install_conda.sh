#!/bin/bash
: '
This script installs conda via the downloadable miniconda tarball.
'
pushd "$(dirname "$0")" >&2

machine=$(bash determine_os.sh)
[ "$machine" = "Linux" ] && dist="Linux" || dist="MacOSX"
miniconda_file="Miniconda3-latest-$dist-x86_64.sh"
miniconda_url="https://repo.anaconda.com/miniconda/$miniconda_file"
curl -O $miniconda_url >&2
bash $miniconda_file -b  >&2

popd
