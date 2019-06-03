#!/bin/bash
: '
This script sets up the gmsh dependency for DiscretePDEs.jl for either Mac or Ubuntu
Linux. It does this via the following steps:

1. Install gmsh (if not installed)
  a. If Mac, install via brew.
  b. If Linux, make from source.
2. Create a .env file with a value for the path to the gmsh.jl file, included
  via the installation above.

## Undoing the changes made by this script

Undo the changes made by setup_gmsh.sh by performing the following steps:

On Mac:
1. Run `brew uninstall gmsh`
2. Remove the entry in the .env file in DiscretePDEs.jl for gmshjlpath=...

On Linux:
1. Since gmsh was made from source, many files were added to
   /usr/local/(lib|bin|include). Care should be taken when attempting to undo
   this process.
2. Remove the entry in the .env file in DiscretePDEs.jl for gmshjlpath=...
3. *Note*: Numerous dependencies were installed via apt-get/yum, as seen in
   build_gmsh, but since some of these may already have been installed,
   it is recommended that you do not uninstall any of them.
'

: '
  Helper functions
'
install_gmsh () {
  : '
  This function takes an operating system machine type as input (e.g. Mac or
  Linux), installs gmsh, and then echos the path to the installed gmsh as
  output.
  '
  machine=$1
  if [ "$machine" = "Mac" ]; then
    brew install gmsh >&2 # Direct brew output to something other than stdout
  elif [ "$machine" = "Linux" ]; then
    install_gmsh_linux >&2 # Direct output to something other than stdout
  fi
  echo $(which gmsh)
}

install_gmsh_linux () {
  : '
  This function installs gmsh for ubuntu machines. The SDK is unreliable, so
  gmsh must be made from source.
  TODO: Fix SDK installation, allow as option.
  '
  source build_gmsh.sh >&2
}

install_gmsh_linux_from_sdk () {
  : '
  NOTE: Not currently in use, since this does not work.
  '
  # Install gmsh sdk
  gmsh_sdk_folder=gmsh-4.2.3-Linux64-sdk
  gmsh_sdk_tarfile=$gmsh_sdk_folder.tgz
  gmsh_sdk_url=http://gmsh.info/bin/Linux/$gmsh_sdk_tarfile
  pushd ~
  curl -O $gmsh_sdk_url
  tar -xvf $gmsh_sdk_tarfile
  rm $gmsh_sdk_tarfile

  # Add content to bashrc file for running
  pushd $gmsh_sdk_folder/bin
  cat << EOF >> ~/.bashrc
# Added by setup_gmsh.sh in DiscretePDEs.jl. Adds the /usr/lib/ folder to the
# LD_LIBRARY_PATH, and sets an alias for gmsh to the executable in the
# dowloaded SDK folder.
# For more information, see https://gitlab.com/rigetti/swe/DiscretePDEs.jl
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu"
EOF
  ln -s $(pwd)/gmsh /usr/local/bin/gmsh
  source ~/.bashrc

  popd
  popd
}

install_gmsh_linux_from_build () {
  : '
  NOTE: Not currently in use. Only linux option is to build gmsh from source.
  '
  # TEMPORARY RESTRICTION: Exit installation if /usr/local/(bin|lib|include)
  # have any files in them. TODO: Test if installation is robust against
  # existing files in these folders.
  if [ $(bash folder_is_empty.sh "/usr/local/lib") = "false" ]; then
    echo "ERROR: /usr/local/lib is not empty, exiting setup."
    exit 1
  fi
  if [ $(bash folder_is_empty.sh "/usr/local/include") = "false" ]; then
    echo "ERROR: /usr/local/include is not empty, exiting setup."
    exit 1
  fi
  # apt-get install --assume-yes git-lfs
  # git lfs install
  build_file=ubuntu-18.10-gmsh-4.2.3-usrlocal.tar.gz
  tar -xvf $build_file
  cp -r ./local/bin/* /usr/local/bin/
  cp -r ./local/include/* /usr/local/include/
  cp -r ./local/lib/* /usr/local/lib/
  rm -rf ./local
}

get_gmshjl_path () {
  : '
  This function takes the output from install_gmsh function (a supposed path to
  the gmsh command, and echoes as output the path to gmsh.jl, via any symlinks.
  '
  gmshpath=$1 >&2
  gmshdir=${gmshpath%/*}

  # 2b. Point to actual gmsh, if symlinked from above
  if test -h $gmshpath; then
    pushd $gmshdir >&2
    relgmshpath=$(readlink -n $gmshpath)
    relgmshdir=${relgmshpath%/*}
    pushd $relgmshdir >&2
    gmshdir="$(pwd)"
    gmshpath="$gmshdir/gmsh"
    popd >&2
    popd >&2
  fi

  # 2c. Point to gmsh.jl
  pushd $gmshdir/../lib/ >&2
  gmshjlpath=$(pwd)/gmsh.jl

  if test -f $gmshjlpath; then
    echo "Path to gmsh.jl: $gmshjlpath." >&2
  else
    echo "File $gmshjlpath not found, exiting setup."
    exit 1
  fi

  popd >&2
  echo $gmshjlpath
}

: '
  Begin setup script.
'
# 1. Determine operating system, from https://stackoverflow.com/questions/3466166/how-to-check-if-running-in-cygwin-mac-or-linux
echo "Determining operating system..."
machine=$(bash determine_machine.sh)
echo "System is $machine."

# 2. Install gmsh
echo "Determining gmsh installation..."
gmshpath=$(which gmsh)

# 2a. If no gmsh found, install it
if [ -z "$gmshpath" ]; then
  echo ""
  echo "**********************************************"
  echo "No gmsh installation found. Installing gmsh..."
  echo "**********************************************"
  echo ""
  gmshpath=$(install_gmsh $machine)
fi

gmshjlpath=$(get_gmshjl_path $gmshpath)

# 3. Set correct value for gmsh.jl path in a .env file
echo ""
echo "****************************************************************************"
echo "gmsh installed at $gmshpath, using this installation of gmsh for DiscretePDEs.jl."
echo "****************************************************************************"
echo ""

bash replace_envvar.sh "gmshjlpath=$gmshjlpath" ".."
