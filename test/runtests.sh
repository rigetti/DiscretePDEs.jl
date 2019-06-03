#!/bin/bash

# Currently needed to clone AdmittanceModels.jl and DiscreteExteriorCalculus.jl
# which are private github repos.
eval $(ssh-agent -s) && ssh-add <(echo "$SSH_PRIVATE_KEY_GITHUB")

pushd "$(dirname "$0")"

pushd ..
source ~/.bashrc
julia --project=@. -e "import Pkg; Pkg.build(); Pkg.test(; coverage=true)"
popd
