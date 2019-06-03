#!/bin/bash
: '
This file takes the absolute path of a directory as input, and echoes "true" if
that directory is empty, or "false" if not.

This solution was taken from
https://www.cyberciti.biz/faq/linux-unix-shell-check-if-directory-empty/.
'
directory=$1
[ "$(ls -A $directory)" ] && echo "false" || echo "true"
