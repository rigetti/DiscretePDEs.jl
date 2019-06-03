#!/bin/bash
: '
TODO: Docstring.
'

install_julia_linux () {
  apt-get update
  apt-get install --assume-yes curl

  julia_folder=julia-1.1.0
  julia_file=$julia_folder-linux-x86_64.tar.gz
  julia_url=https://julialang-s3.julialang.org/bin/linux/x64/1.1/$julia_file

  pushd ~
  curl -O -sSL $julia_url
  tar -xvf $julia_file
  pushd $julia_folder
  ln -s $(pwd)/bin/julia /usr/local/bin/julia
  popd
  rm $julia_file
  popd
}

install_julia_mac () {
  brew cask install julia # TODO: Hardcode julia version
}

if [ ! -z $(which julia) ]; then
  echo ""
  echo "*******************************************"
  echo "Julia already installed at $(which julia), skipping installation."
  echo "*******************************************"
  echo ""
  exit 0
fi

echo ""
echo "*******************************************"
echo "No julia command found, installing julia..."
echo "*******************************************"
echo ""


if [ $(bash determine_machine.sh) = "Mac" ]; then
  install_julia_mac
else
  install_julia_linux
fi

echo ""
echo "**********************************"
echo "Julia installed at $(which julia)."
echo "**********************************"
echo ""
