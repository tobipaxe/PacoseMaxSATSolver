#!/bin/bash

directory=$(dirname $0)
rand=$RANDOM
echo "c Random prefix: $rand"
oldwcnf=/tmp/"$rand".wcnf

echo "c Converting to old format:"
#"$directory"/../TobiasPaxian/maxsat_benchmarks_code_base/bin/to_old_fmt "$1" > "$oldwcnf"
"$directory"/wcnfverifier -c "$oldwcnf" "$1"

eval "$directory"/../bin/Pacose "$oldwcnf"
ret=$?

rm -f "$oldwcnf"
exit $ret
