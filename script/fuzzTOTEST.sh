#! /bin/bash

wcnf="/tmp/$RANDOM$RANDOM$RANDOM.wcnf"

rv=1
while [[ $rv -eq "1" ]]; do
    ./wcnfuzz --wcnf > "$wcnf"
    /usr/bin/timeout --kill-after=2s 2s ../bin/Pacose --TrimMaxSAT "$wcnf"
    grep "^* TOTEST!" maxsat_proof.pbp
    rv=$?
done
cp "$wcnf" .
rm "$wcnf"
