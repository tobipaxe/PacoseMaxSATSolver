while read p; do 
  echo "$p";
  bn=$(basename "$p")
  dn=$(dirname "$p")

  echo '/usr/local/scratch/paxiant/MaxSATFuzzer/TobiasPaxian/runfuzz/wcnfddmin -e 0 -r "$dn"/red-"$bn" -s /usr/local/scratch/paxiant/MaxSATFuzzer/MaxSATSolver/compare.py_--solvers_PacoseNew_--logAll_--saveWCNF -p 0 "$p"'
  echo "$dn/red-$bn $p"
  /usr/local/scratch/paxiant/MaxSATFuzzer/TobiasPaxian/runfuzz/wcnfddmin -e 0 -r "$dn"/red-"$bn" -s /usr/local/scratch/paxiant/MaxSATFuzzer/MaxSATSolver/compare.py_--solvers_PacoseNew_--logAll_--saveWCNF -p 0 "$p"
  echo ""
done < $1
