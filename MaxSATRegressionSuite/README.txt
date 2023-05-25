All unique not empty files from the reduction process.

First version of a test script running all WCNFs in these three subfolders.

./run_tests.sh SOLVERBINARY

If the binaries of kissat or / and wcnfverifier do not work build them by yourself! 
You can find the wcnfverifier in the wcnfverifier subfolder and kissat in: https://github.com/arminbiere/kissat
 

TODO: Call WCNFddmin for faulty WCNFs from the reduction process of faultyWCNFs to produce more interesting WCNFs from other categories
TODO: Improve script with command line options to choose which set to check, improve output, ...
