# MaxPre2 MaxSAT preprocessor


## Basic use and flags

The first argument is the instance file, the second is preprocess, reconstruct or solve.

An example of using the preprocessor:  
	./preprocessor test.wcnf preprocess -techniques=[bu]#[buvsrg] -mapfile=test.map > preprocessed.wcnf  
	./solver < preprocessed.wcnf > sol0.sol  
	./preprocessor sol0.sol reconstruct -mapfile=test.map > solution.sol  

Another way to do the same thing:  
	./preprocessor test.wcnf solve -solver=./solver -techniques=[bu]#[buvsrg] > solution.sol  

-techniques (default: [bu]#[buvsrgc])  
	string:  
	This string defines the preprocessing techniques to use and the order of them.  
	Each letter corresponds to a preprocessing technique. Each preprocessing technique is applied until its fixpoint.  
	Techniques inside brackets are applied until all of them are in fixpoint. The brackets work recursively.  
	If # character is given, all techniques before it are applied before group detection and adding labels (techniques available before labeling are BCE, UP and SE).  
	Techniques:  
	b = blocked clause elimination  
	u = unit propagation  
	v = bounded variable elimination  
	s = subsumption elimination  
	r = self subsuming resolution  
	l = subsumed label elimination  
	c = binary core removal  
	a = bounded variable addition  
	g = Group subsumed label elimination  
	e = equivalence elimination  
	h = unhiding techniques (on binary implication graph, failed literals, hidden tautology elimination, hidden literal elimination)  
	t = structure labeling  
	G = intrinsic atmost1 constraints  
	T = TrimMaxSat  
	V = TrimMaxSat -algorithm on all variables to detect backbones  
	H = hardening  
	R = failed literal elimination  

-solver (default: disabled)
	string:
	The solver to use to solve the preprocessed instance

-solverflags (default: disabled)
	string:
	The flags to use with the solver
	For example -solver=./LMHS -solverflags="--infile-assumps --no-preprocess" results in using the command ./LMHS preprocessed.wcnf --infile-assumps --no-preprocess > sol0.sol

-mapfile (default: disabled)
	string:
	The file to write the solution reconstruction map

-problemtype (default: maxsat)
	string: {maxsat, sat}
	Should the problem be preprocessed as a MaxSAT or SAT instance

-outputformat (default: wpms)
	string: {original}
	By default the preprocessor always gives the output in weighted partial MaxSAT format
	Output in SAT format by setting this to original when preprocessing SAT instances

-timelimit (default: inf)
	double: [0, 500000000]
	Limit for preprocessing time in seconds

-skiptechnique (default: disabled)
	int: [1, 1000000000]
	Skip a preprocessing technique if it seems to be not effective in x tries (x is given in this flag)
	Recommended values for this could be something between 10 and 1000

-matchlabels (default: 0)
	bool: {0, 1}
	Use label matching technique to reduce the number of labels

-bvegate (default: 0)
	bool: {0, 1}
	Use BVE gate extraction to extend BVE
	Note: applying BCE will destroy all recognizable gate structures

-verb (default: 1)
	int: [0, 1]
	If verb is 0 the preprocessor will output less stuff to the standard error


## More detailed information

## Actual order of simplifyings

* remove tautologies
* remove empty clauses
* remove duplicate clauses
* apply preprocessing techniques specified before the #-character
* remove empty clauses
* group detection
* label creation
* label matching (if enabled)
* remove duplicate clauses
* apply preprocessing techniques
* remove empty clauses
* remove duplicate clauses

## The techniques string

Each technique is applied somewhat modularly indenpendent of each other. The user
can specify the exact order of used preprocessing techniques with the techniques string.
Per default, each technique is always applied until "fixpoint" - the point after which that 
technique is unable to simplify the instance further. The exception to this is 
The techniques inside brackets are applied
until none of them change the instance when applying all of them in the given order.
The brackets work recursively, for example [[[[vu]b]sr]ea] is valid syntax.

## Timelimit

You can set internal time limit for the preprocessor running time with the timelimit flag.
The preprocessor will try to preprocess the instance in less than the given time,
however it is hard to estimate I/O times with large files, so when preprocessing very
large instances you should take that into account. The preprocessor tries to limit
the time used by each of the preprocessing techniques somewhat independly of each other.
Each technique will be allocated its own time limit (the proportions of time given
to each preprocessing technique are hardcoded in log.cpp file in Log::timePlan
function. When some techniques do not use all of their allocated time, it will be
given to other techniques with some heuristics. Note that by specifying the timelimit,
the preprocessor might work less efficiently overall. For example it could preprocess
everything to fixpoint in 60 seconds without the timelimit flag, but with -timelimit=60
it could use only 30 seconds not get as much preprocessing done. It is not recommended
to try to optimize the time used by the preprocessor by using the timelimit flag, but
rather to use it for an upper bound for the time used by the preprocessor.

## API

Maxpre offers an API for integration with MaxSAT solvers. Use make lib to make
the static library file and include preprocessorinterface.hpp to use it. The API
is documented in preprocessorinterface.hpp.