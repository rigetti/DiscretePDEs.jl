: '
This script checks the machine via uname and echoes as output either
Mac, Linux, or UKNOWN:...
'
unameOut="$(uname -s)"
case "${unameOut}" in
  Linux*)     machine=Linux;;
  Darwin*)    machine=Mac;;
  *)          machine="UNKNOWN:${unameOut}"
esac
echo $machine
