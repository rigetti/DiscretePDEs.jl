#!/bin/bash
: '
Setup a conda environment, named DiscretePDEs, for DiscretePDEs.jl, since it depends on
PyCall and gdspy (see build.jl). This script

1. Checks for a conda installation, and if it does not find one, installs one
   via miniconda, placing the installation folder in your home directory. Note:
   this will add content to your ~/.bash_profile to make the conda command
   available.
2. Checks for a conda virtual environment called DiscretePDEs, and if it does not
   find one, creates it.
3. Adds PYTHON=... to the .env file in DiscretePDEs.jl.

## Undoing the changes made by this script.

Undo the changes made by setup_conda_env.sh by performing the following steps:

1. Remove the content from your ~/.bash_profile or ~/.bash_aliases for conda.
2. Remove the miniconda folder from your home directory.
3. Remove the entry in the .env file in DiscretePDEs.jl for PYTHON=...

TODO: Investigate using conda configuration in .julia folder.
'
# Helper functions


install_miniconda () {
  : '
  This function installs conda via the downloadable miniconda tarball. When
  finished, it sources ~/.bash_profile so that the conda command is available
  for use.
  '
  machine=$(bash determine_machine.sh)
  pushd ~ >&2
  echo "Machine: $machine" >&2
  [ "$machine" = "Linux" ] && dist="Linux" || dist="MacOSX"
  miniconda_file="Miniconda3-latest-$dist-x86_64.sh"
  miniconda_url="https://repo.anaconda.com/miniconda/$miniconda_file"
  curl -O $miniconda_url >&2
  bash $miniconda_file -b  >&2

  # NOTE: The following cat to .bash_profile/.bash_aliases was lifted from the
  # interactive miniconda installation, since this step is ignored when
  # installation is run silently with -b.
  [ "$machine" = "Linux" ] && bashfile="$HOME/.bash_aliases" || bashfile="$HOME/.bash_profile"
  touch $bashfile
  cat << EOF >> $bashfile
# added by Miniconda3 installer
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$(CONDA_REPORT_ERRORS=false '$HOME/miniconda3/bin/conda' shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
    \eval "$__conda_setup"
else
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        . "$HOME/miniconda3/etc/profile.d/conda.sh"
        CONDA_CHANGEPS1=false conda activate base
    else
        \export PATH="$HOME/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda init <<<
EOF
  source $bashfile >&2
  popd >&2
}

extract_env_path () {
  : '
  This function checks for the path to the DiscretePDEs conda virtual environment
  via the `conda env list` command, and echoes the result as output.
  '
  echo $(conda env list | grep DiscretePDEs | awk '{print $2}')
}

: '
  Begin setup script.
'
set +e # Wrap this command in +e as it can give status code 1
condapath=$(which conda) >/dev/null
set -e

if [ -z "$condapath" ]; then
  echo ""
  echo "********************************************************"
  echo "No conda installation found. Installing via miniconda..."
  echo "********************************************************"
  echo ""

  install_miniconda
else
  echo ""
  echo "********************************************************"
  echo "Conda installation found at $(which conda), using this one."
  echo "********************************************************"
  echo ""
fi

DiscretePDEs_env_path=$(extract_env_path)
if [ -z "$DiscretePDEs_env_path" ]; then
  echo ""
  echo "***************************************************"
  echo "No DiscretePDEs conda environment found, creating one..."
  echo "***************************************************"
  echo ""
  conda env create -f ../deps/conda-env.yml
  DiscretePDEs_env_path=$(extract_env_path)
  echo ""
  echo "***************************************************"
  echo "Conda environment DiscretePDEs created"
  echo "***************************************************"
  echo ""
else
  echo ""
  echo "***************************************************"
  echo "Found DiscretePDEs conda env at $DiscretePDEs_env_path, using this one."
  echo "***************************************************"
  echo ""
fi

bash replace_envvar.sh "PYTHON=$DiscretePDEs_env_path/bin/python" ".."
