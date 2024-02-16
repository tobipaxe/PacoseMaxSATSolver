#! /bin/bash

ToTest=$1
if [ -z "$1" ]; then
  ToTest="TOTEST"
fi

echo "Try to find $ToTest or SHOULDNEVERHAPPEN in the proof file!"

wcnf="/tmp/$RANDOM$RANDOM$RANDOM$RANDOM.wcnf"

(cd ../build && make -j)
(cd wcnfdddmin && make -j)
mv wcnfdddmin/wcnfddmin .

rv=0
counter=0
while [[ $rv -eq "0" ]]; do
    ((counter++))
    echo ""
    echo -E "ROUND $counter"
    ./wcnfuzz.py -u 4611686018427387904 > "$wcnf"
    # ./wcnfuzz.py --gbmo -u 4611686018427387904 > "$wcnf"
    ./callPacose.sh "$ToTest" "$wcnf"
    rv=$?
done

echo "Grep found either $ToTest or SHOULDNEVERHAPPEN, reduce instance!"
baseWcnf=$(basename "$wcnf")
reducedWcnf="red-$baseWcnf"
var="./wcnfddmin -e 0 -r $reducedWcnf -s ./callPacose.sh_$ToTest $wcnf"
echo "$var"
eval "$var"

echo ""
echo "Found $ToTest after $counter loops. It is reduced to $reducedWcnf"

cp "$wcnf" .
rm "$wcnf"
