#! /bin/bash

echo "./callPacose.sh $1 $2"
/usr/bin/timeout --kill-after=2s 60s ./certified-pacose.sh "$2" "$1"
result=$?

pacose_rv=0
if [ $result -eq 0 ]; then
  echo "certified-pacose.sh executed successfully within time limit."
elif [ $result -eq 124 ]; then
  echo "certified-pacose.sh timed out."
elif [ $result -eq 137 ]; then
  echo "certified-pacose.sh killed after additional grace period."
else
  echo "certified-pacose.sh exited with code: $result"
  pacose_rv=$result
fi

echo "c ReturnValue of certified pacose script: $pacose_rv"
exit $pacose_rv

# /usr/bin/timeout --kill-after=2s 2s ../bin/Pacose --GBMO --TrimMaxSAT "$2"
# #echo "grep -E \"$1|SHOULDNEVERHAPPEN\" maxsat_proof.pbp"
# grep -E "$1|SHOULDNEVERHAPPEN" maxsat_proof.pbp
# if [[ $? -eq "0" ]]; then
#   exit 1
# fi
# exit 0
