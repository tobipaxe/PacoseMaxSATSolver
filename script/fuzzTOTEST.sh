#! /bin/bash

ToTest=$1
if [ -z "$1" ]; then
  ToTest="TOTEST"
fi

wcnf="/tmp/$RANDOM$RANDOM$RANDOM.wcnf"

(cd ../build && make -j)
(cd wcnffuzz && make)
#(cd wcnfddmin && make)
mv wcnffuzz/wcnfuzz .
#mv wcnfddmin/wcnfddmin .

rv=1
counter=0
while [[ $rv -eq "1" ]]; do
    ((counter++))
    ./wcnfuzz --wcnf > "$wcnf"
    /usr/bin/timeout --kill-after=2s 2s ../bin/Pacose --TrimMaxSAT "$wcnf"
    grep -E "$ToTest|SHOULDNEVERHAPPEN" maxsat_proof.pbp
    rv=$?
done

echo "Found $ToTest after $counter loops"
cp "$wcnf" .
# TODO Tobias: Minimize wcnf
rm "$wcnf"
