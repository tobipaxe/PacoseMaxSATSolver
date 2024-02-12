#! /bin/bash

rm maxsat_proof.pbp
echo "./callPacose.sh $1 $2"
/usr/bin/timeout --kill-after=2s 3600s ../bin/Pacose --verbose=20 "$2"
rv=$?
echo "c ReturnValue: $rv"
exit $rv

/usr/bin/timeout --kill-after=2s 2s ../bin/Pacose --GBMO --TrimMaxSAT "$2"
#echo "grep -E \"$1|SHOULDNEVERHAPPEN\" maxsat_proof.pbp"
grep -E "$1|SHOULDNEVERHAPPEN" maxsat_proof.pbp
if [[ $? -eq "0" ]]; then
  exit 1
fi
exit 0
