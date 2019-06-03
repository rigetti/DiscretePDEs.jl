#!/bin/bash
: '
This script is used to replace values in a .env file. It
is executed by calling bash replace_envvar.sh "<key>=<val>" "</path/to/.env>".
If <key> already exists in .env, its value will be replaced with <val>. If <key>
does not exist, it will be added.
'
inkey=$(cut -d'=' -f1 <<<$1)
inval=$(cut -d'=' -f2 <<<$1)

directory=$2

pushd $2 >&2
infile=".env"
# Create .env if it doesn't exist.
if [ ! -f .env ]; then
  touch .env
fi

# If inkey doesn't exist in .env, we can simply add it and quit.
if [ -z $(egrep -v '^#' .env | grep $inkey= | cut -d'=' -f2) ]; then
  echo "$inkey=$inval" >> .env
  exit 0
fi

# inkey does exist in .env, so we must replace it.
outfile=""
while IFS= read line; do
  if [[ $line = '#'* ]]; then # Leave comments as is
    outfile="$outfile$line\n"
  else
    key=$(cut -d'=' -f1 <<<$line)
    val=$(cut -d'=' -f2 <<<$line)
    if [[ "$key" = "$inkey" ]]; then
      val=$inval
    fi
    outfile="$outfile$key=$val\n"
  fi
done < "$infile"
echo -e $outfile > .env
popd >&2
