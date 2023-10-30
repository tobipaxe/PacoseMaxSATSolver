#! /bin/bash

wcnf="/tmp/$RANDOM$RANDOM$RANDOM.wcnf"

(cd ../build && make -j)
(cd wcnffuzz && make)
cp wcnffuzz/wcnfuzz .

rv=1
while [[ $rv -eq "1" ]]; do
    ./wcnfuzz --wcnf > "$wcnf"
    /usr/bin/timeout --kill-after=2s 2s ../bin/Pacose --TrimMaxSAT "$wcnf"
    grep -E "TOTEST|SHOULDNEVERHAPPEN" maxsat_proof.pbp
    rv=$?
done
cp "$wcnf" .
rm "$wcnf"
