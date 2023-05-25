#include <iostream>
#include <fstream>
#include <chrono>
#include <cassert>
#include <sstream>
#include <map>
#include <iomanip>

#include "preprocessorinterface.hpp"
#include "inputreader.hpp"
#include "outputreader.hpp"
#include "timer.hpp"
#include "utility.hpp"
using namespace std;

map<string, string> getFlags(int argc, char* argv[]) {
	map<string, string> ret;
	for (int i = 3; i < argc; i++) {
		string s(argv[i]);
		if (s.size() == 0 || s[0] != '-') {
			cout<<"c Invalid arg "<<s<<endl;
			cerr<<"Invalid arg "<<s<<endl;
		}
		else {
			int p = -1;
			int prd = 0;
			for (int j = 0; j < (int)s.size(); j++) {
				if (s[j] == '=') {
					p = j;
					break;
				}
				if (j == prd && s[j] == '-') prd++;
			}
			if (p == -1) {
				ret[s.substr(prd)] = "";
			}
			else {
				ret[s.substr(prd, p-prd)] = s.substr(p+1);
			}
		}
	}
	return ret;
}

int parseInt(string& s) {
	if (s=="inf") return ~(1<<31);
	if (s=="-inf") return (1<<31);
	stringstream ss;
	ss<<s;
	int r;
	ss>>r;
	return r;
}
unsigned parseUnsigned(string& s) {
	if (s=="inf") return ~(0u);
	stringstream ss;
	ss<<s;
	unsigned r;
	ss>>r;
	return r;
}
bool parseBool(string& s) {
	if (s=="1") return 1;
	if (s=="0" || s=="false" || s=="False" || s=="FALSE") return 0;
	return 1;
}
double parseDouble(string& s) {
	if (s=="inf") return 1e9;
	if (s=="-inf") return -(1e9);
	stringstream ss;
	ss<<s;
	int r;
	ss>>r;
	return r;
}
double parseUint64(string& s) {
	if (s=="inf") return ~(0ull);
	stringstream ss;
	ss<<s;
	uint64_t r;
	ss>>r;
	return r;
}

string parseTechniques(map<string, string>& flags) {
	string techniques = "[bu]#[buvsrgc]";
	if (flags.count("techniques")) {
		techniques = flags["techniques"];
		cout<<"c Techniques "<<techniques<<endl;
		cerr<<"Techniques "<<techniques<<endl;
	}
	else {
		cout<<"c No -techniques= given, defaulting to "<<techniques<<endl;
		cerr<<"No -techniques= given, defaulting to "<<techniques<<endl;
	}
	return techniques;
}

double parseTimeLimit(map<string, string>& flags) {
	double timeLimit = 1e9;
	if (flags.count("timelimit")) {
		stringstream ss;
		ss<<flags["timelimit"];
		ss>>timeLimit;
		cout<<"c Preprocess time limit "<<timeLimit<<endl;
		cerr<<"Preprocess time limit "<<timeLimit<<endl;
	}
	else {
		cout<<"c No -timelimit= given, defaulting to inf"<<endl;
		cerr<<"No -timelimit= given, defaulting to inf"<<endl;
	}
	return timeLimit;
}

bool parseBVEgate(map<string, string>& flags) {
	bool BVEgate = false;
	if (flags.count("bvegate")) {
		if (flags["bvegate"] == "1") {
			BVEgate = true;
			cout<<"c BVE gate extraction enabled"<<endl;
			cerr<<"BVE gate extraction enabled"<<endl;
		}
		else if (flags["bvegate"] == "0") {
			BVEgate = false;
			cout<<"c BVE gate extraction disabled"<<endl;
			cerr<<"BVE gate extraction disabled"<<endl;
		}
		else {
			cout<<"Invalid bvegate flag"<<endl;
			cerr<<"Invalid bvegate flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -bvegate= given, defaulting to disabled"<<endl;
		cerr<<"No -bvegate= given, defaulting to disabled"<<endl;
	}
	return BVEgate;
}

bool parseLabelMatching(map<string, string>& flags) {
	bool labelMatching = false;
	if (flags.count("matchlabels")) {
		if (flags["matchlabels"] == "1") {
			labelMatching = true;
			cout<<"c Label matching enabled"<<endl;
			cerr<<"Label matching enabled"<<endl;
		}
		else if (flags["matchlabels"] == "0") {
			labelMatching = false;
			cout<<"c Label matching disabled"<<endl;
			cerr<<"Label matching disabled"<<endl;
		}
		else {
			cout<<"Invalid matchlabels flag"<<endl;
			cerr<<"Invalid matchlabels flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -matchlabels given, defaulting to disabled"<<endl;
		cerr<<"No -matchlabels given, defaulting to disabled"<<endl;
	}
	return labelMatching;
}

bool parseProblemType(map<string, string>& flags) {
	bool maxSat = true;
	if (flags.count("problemtype")) {
		for (char& c : flags["problemtype"]) {
			c = tolower(c);
		}
		if (flags["problemtype"] == "sat") {
			maxSat = false;
			cout<<"c Problem type is SAT"<<endl;
			cerr<<"Problem type is SAT"<<endl;
		}
		else if (flags["problemtype"] == "maxsat" || flags["problemtype"] == "max-sat") {
			maxSat = true;
			cout<<"c Problem type is Max-SAT"<<endl;
			cerr<<"Problem type is Max-SAT"<<endl;
		}
		else {
			cout<<"Invalid problemtype flag"<<endl;
			cerr<<"Invalid problemtype flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -problemtype given, defaulting to Max-SAT"<<endl;
		cerr<<"No -problemtype given, defaulting to Max-SAT"<<endl;
	}
	return maxSat;
}

int parseSkipTechnique(map<string, string>& flags) {
	int skipTechnique = 0;
	if (flags.count("skiptechnique")) {
		stringstream ss;
		ss<<flags["skiptechnique"];
		ss>>skipTechnique;
		cout<<"c Skiptechnique "<<skipTechnique<<endl;
		cerr<<"Skiptechnique "<<skipTechnique<<endl;

		if (skipTechnique <= 0 || skipTechnique > 1000000000) {
			cout<<"Invalid skiptechnique flag"<<endl;
			cerr<<"Invalid skiptechinque flag"<<endl;
			exit(0);
		}
	}
	else {
		cout<<"c No -skiptechnique given, defaulting to disabled"<<endl;
		cerr<<"No -skiptechnique given, defaulting to disabled"<<endl;
	}
	return skipTechnique;
}

bool parseBVEsortMaxFirst(map<string, string>& flags) {
  bool BVEsortMaxFirst = false;
  if (flags.count("bvesortmaxfirst")) {
    if (flags["bvesortmaxfirst"] == "1") {
      BVEsortMaxFirst = true;
      cout<<"c BVEsortMaxFirst enabled"<<endl;
      cerr<<"BVEsortMaxFirst enabled"<<endl;
    }
    else if (flags["bvesortmaxfirst"] == "0") {
      BVEsortMaxFirst = false;
      cout<<"c BVEsortMaxFirst disabled"<<endl;
      cerr<<"BVEsortMaxFirst disabled"<<endl;
    }
    else {
      cout<<"Invalid bvesortmaxfirst flag"<<endl;
      cerr<<"Invalid bvesortmaxfirst flag"<<endl;
      exit(0);
    }
  }
  else {
    cout<<"c No -bvesortmaxfirst given, defaulting to disabled"<<endl;
    cerr<<"No -bvesortmaxfirst given, defaulting to disabled"<<endl;
  }
  return BVEsortMaxFirst;
}

int parseBVElocalGrow(map<string, string>& flags) {
  int BVElocalGrow = 0;
  if (flags.count("bvelocalgrow")) {
    stringstream ss;
    ss<<flags["bvelocalgrow"];
    ss>>BVElocalGrow;
    cout<<"c BVElocalgrow "<<BVElocalGrow<<endl;
    cerr<<"BVElocalgrow "<<BVElocalGrow<<endl;

    if (BVElocalGrow <= 0 || BVElocalGrow > 1000000000) {
      cout<<"Invalid bvelocalgrow flag"<<endl;
      cerr<<"Invalid bvelocalgrow flag"<<endl;
      exit(0);
    }
  }
  else {
    cout<<"c No -bvelocalgrow given, defaulting to disabled"<<endl;
    cerr<<"No -bvelocalgrow given, defaulting to disabled"<<endl;
  }
  return BVElocalGrow;
}

int parseBVEglobalGrow(map<string, string>& flags) {
  int BVEglobalGrow = 0;
  if (flags.count("bveglobalgrow")) {
    stringstream ss;
    ss<<flags["bveglobalgrow"];
    ss>>BVEglobalGrow;
    cout<<"c BVEglobalgrow "<<BVEglobalGrow<<endl;
    cerr<<"BVEglobalgrow "<<BVEglobalGrow<<endl;

    if (BVEglobalGrow <= 0 || BVEglobalGrow > 1000000000) {
      cout<<"Invalid bveglobalgrow flag"<<endl;
      cerr<<"Invalid bveglobalgrow flag"<<endl;
      exit(0);
    }
  }
  else {
    cout<<"c No -bveglobalgrow given, defaulting to disabled"<<endl;
    cerr<<"No -bveglobalgrow given, defaulting to disabled"<<endl;
  }
  return BVEglobalGrow;
}



// common variable parsing
void parseIntVars(map<string, string>& flags, map<string, int>& intVars) {
	for (auto& v : intVars) {
		const string& key=v.first;
		if (flags.count(key)) {
			intVars[key]=parseInt(flags[key]);
			cout << "c " << key << " " << intVars[key] << endl;
			cerr << key << " " << intVars[key] << endl;
		}
	}
}
void parseBoolVars(map<string, string>& flags, map<string, bool>& boolVars) {
	for (auto& v : boolVars) {
		const string& key=v.first;
		if (flags.count(key)) {
			boolVars[key]=parseBool(flags[key]);
			cout << "c " << key << " " << boolVars[key] << endl;
			cerr << key << " " << boolVars[key] << endl;
		}
	}
}
void parseDoubleVars(map<string, string>& flags, map<string, double>& doubleVars) {
	for (auto& v : doubleVars) {
		const string& key=v.first;
		if (flags.count(key)) {
			doubleVars[key]=parseDouble(flags[key]);
			cout << "c " << key << " " << doubleVars[key] << endl;
			cerr << key << " " << doubleVars[key] << endl;
		}
	}
}
void parseUint64Vars(map<string, string>& flags, map<string, uint64_t>& uint64Vars) {
	for (auto& v : uint64Vars) {
		const string& key=v.first;
		if (flags.count(key)) {
			uint64Vars[key]=parseUint64(flags[key]);
			cout << "c " << key << " " << uint64Vars[key] << endl;
			cerr << key << " " << uint64Vars[key] << endl;
		}
	}
}



pair<unsigned, unsigned> parseSizeLimit(map<string, string>& flags) {
	pair<unsigned, unsigned> sizeLimit={0,~0u};
	if (flags.count("sizelimit")) {
		stringstream ss;
		ss<<flags["sizelimit"];
		string t1;
		string t2;
		getline(ss, t1, '-');
		getline(ss, t2, '-');
		sizeLimit.first=parseUnsigned(t1);
		sizeLimit.second=parseUnsigned(t2);
		cout<<"c sizeLimit: "<< (sizeLimit.first==(~0u)?"inf":to_string(sizeLimit.first)) << "-" << (sizeLimit.second==(~0u)?"inf":to_string(sizeLimit.second)) <<endl;
		cerr<<"sizeLimit: "<< (sizeLimit.first==(~0u)?"inf":to_string(sizeLimit.first)) << "-" << (sizeLimit.second==(~0u)?"inf":to_string(sizeLimit.second)) <<endl;
	}
	return sizeLimit;
}

string parsePrepFilename(map<string, string>& flags) {
	if (flags.count("prepfile")) {
		return flags["prepfile"];
	}
	return "preprocessed.wcnf";
}

string parseSolFilename(map<string, string>& flags) {
	if (flags.count("solfile")) {
		return flags["solfile"];
	}
	return "sol0.sol";
}


int parseVerb(map<string, string>& flags) {
	int verb = 1;
	if (flags.count("verb")) {
		if (flags["verb"] == "0") {
			verb = 0;
		}
		else if (flags["verb"] == "1") {
			verb = 1;
		}
		else if (flags["verb"] == "2") {
			verb = 2;
		}
		else {
			cout<<"Invalid verb flag"<<endl;
			cerr<<"Invalid verb flag"<<endl;
			exit(0);
		}
		cout<<"c Verb "<<verb<<endl;
		cerr<<"Verb "<<verb<<endl;
	}
	else {
		cout<<"c No -verb given, defaulting to 1"<<endl;
		cerr<<"No -verb given, defaulting to 1"<<endl;
	}
	return verb;
}

void printHelp(ostream& out, map<string, int>& intVars, map<string, bool>& boolVars, map<string, double>& doubleVars, map<string, uint64_t> uint64Vars, bool shrt) {
	out<<"The first argument is the instance file, the second is preprocess, reconstruct or solve."<<endl;
	out<<endl;

	out<<"An example of using the preprocessor:"<<endl;
	out<<"\t./preprocessor test.wcnf preprocess -techniques=[bu]#[buvsrg] -mapfile=test.map > preprocessed.wcnf"<<endl;
	out<<"\t./solver < preprocessed.wcnf > sol0.sol"<<endl;
	out<<"\t./preprocessor sol0.sol reconstruct -mapfile=test.map > solution.sol"<<endl;
	out<<endl;
	if (!shrt) {
		out<<"Another way to do the same thing:"<<endl;
		out<<"\t./preprocessor test.wcnf solve -solver=./solver -techniques=[bu]#[buvsrg] > solution.sol"<<endl;
		out<<endl;
	}
	out<<"Parameters:"<<endl;
	out<<endl;

	out<<"-techniques (default: [bu]#[buvsrgc])"<<endl;
	if (!shrt) {
		out<<"\tstring:"<<endl;
		out<<"\tThis string defines the preprocessing techniques to use and the order of them."<<endl;
		out<<"\tEach letter corresponds to a preprocessing technique. Each preprocessing technique is applied until its fixpoint."<<endl;
		out<<"\tTechniques inside brackets are applied until all of them are in fixpoint. The brackets work recursively. "<<endl;
		out<<"\tIf # character is given, all techniques before it are applied before group detection and adding labels (techniques available before labeling are BCE, UP and SE)."<<endl;
		out<<"\tTechniques:"<<endl;
	}
	out<<"\tb = blocked clause elimintation"<<endl;
	out<<"\tu = unit propagation"<<endl;
	out<<"\tv = bounded variable elimination"<<endl;
	out<<"\ts = subsumption elimination"<<endl;
	out<<"\tr = self subsuming resolution"<<endl;
	out<<"\tl = subsumed label elimintion"<<endl;
	out<<"\tc = binary core removal"<<endl;
	out<<"\ta = bounded variable addition"<<endl;
	out<<"\tg = generalized subsumed label elimination"<<endl;
	out<<"\te = equivalence elimination"<<endl;
	out<<"\th = unhiding techniques on binary implication graph (failed literals, hidden tautology elimination, hidden literal elimination)"<<endl;
	out<<"\tt = structure labeling"<<endl;
	out<<"\tG = intrinsic atmost1 constraints"<<endl;
	out<<"\tT = TrimMaxSat"<<endl;
	out<<"\tV = TrimMaxsat based backbone detection"<<endl;
	out<<"\tH = hardening"<<endl;
	out<<"\tR = failed literal elimination + unhiding (extended with redundancy detection)"<<endl;
	if (!shrt) out<<endl;

	out<<"-solver (default: disabled)"<<endl;
	if (!shrt) {
		out<<"\tstring:"<<endl;
		out<<"\tThe solver to use to solve the preprocessed instance"<<endl;
		out<<endl;
	}
	out<<"-solverflags (default: disabled)"<<endl;
	if (!shrt) {
		out<<"\tstring:"<<endl;
		out<<"\tThe flags to use with the solver"<<endl;
		out<<"\tFor example -solver=./LMHS -solverflags=\"--infile-assumps --no-preprocess\" results in using the command ./LMHS preprocessed.wcnf --infile-assumps --no-preprocess > sol0.sol"<<endl;
		out<<endl;
	}
	out<<"-mapfile (default: disabled)"<<endl;
	if (!shrt) {
		out<<"\tstring:"<<endl;
		out<<"\tThe file to write the solution reconstruction map"<<endl;
		out<<endl;
	}
	out<<"-problemtype (default: maxsat)"<<endl;
	if (!shrt) {
		out<<"\tstring: {maxsat, sat}"<<endl;
		out<<"\tShould the problem be preprocessed as a MaxSAT or SAT instance"<<endl;
		out<<endl;
	}
	out<<"-outputformat (default: wpms22)"<<endl;
	if (!shrt) {
		out<<"\tstring: {original, wpms, wpms22, sat}"<<endl;
		out<<"\tBy default the preprocessor always gives the output in weighted partial MaxSAT format"<<endl;
		out<<"\tOutput in SAT format by setting this to original when preprocessing SAT instances"<<endl;
		out<<endl;
	}
	out<<"-timelimit (default: inf)"<<endl;
	if (!shrt) {
		out<<"\tdouble: [0, 500000000]"<<endl;
		out<<"\tLimit for preprocessing time in seconds"<<endl;
		out<<endl;
	}
	out<<"-skiptechnique (default: disabled)"<<endl;
	if (!shrt) {
		out<<"\tint: [1, 1000000000]"<<endl;
		out<<"\tSkip a preprocessing technique if it seems to be not effective in x tries (x is given in this flag)"<<endl;
		out<<"\tRecommended values for this could be something between 10 and 1000"<<endl;
		out<<endl;
	}
	out<<"-matchlabels (default: 0)"<<endl;
	if (!shrt) {
		out<<"\tbool: {0, 1}"<<endl;
		out<<"\tUse label matching technique to reduce the number of labels"<<endl;
		out<<endl;
	}
	out<<"-bvegate (default: 0)"<<endl;
	if (!shrt) {
		out<<"\tbool: {0, 1}"<<endl;
		out<<"\tUse BVE gate extraction to extend BVE"<<endl;
		out<<"\tNote: applying BCE will destroy all recognizable gate structures"<<endl;
		out<<endl;
	}
	out<<"-sizelimit (default: 0-inf)"<<endl;
	if (!shrt) {
		out<<"\tstring: int-int, range for ints [0, 2^32-1], constant inf=2^32-1"<<endl;
		out<<"\tUse preprocessing only if the number of clauses in the instance is on the given range"<<endl;
		out<<"\tFor example -sizelimit=1000000-inf skips preprocessing when there are less than a million clauses on the instance"<<endl;
		out<<endl;
	}
	out<<"-ignore-exit-code (default: not set)"<<endl;
	if (!shrt) {
		out<<"\tBy default if MaxSAT solver exits with a nonzero exit code, maxpre will halt and any solution of the solver is ignored."<<endl;
		out<<"\tUse flag -ignore-exit-code to ignore the exit value of the MaxSAT solver and try to parse and analyze the result anyways."<<endl;
		out<<endl;
	}
	out<<"-prepfile (default: preprocessed.wcnf)"<<endl;
	if (!shrt) {
		out<<"\tstring"<<endl;
		out<<"\tSet the auxiliary file into which the preprocessed instance is saved when type solve is used."<<endl;
		out<<endl;
	}
	out<<"-solfile (default: sol0.sol)"<<endl;
	if (!shrt) {
		out<<"\tstring"<<endl;
		out<<"\tSet the auxiliary file into which the output of the MaxSAT solver is piped when type solve is used."<<endl;
		out<<endl;
	}

	if (intVars.size()) {
		if (!shrt) out<<"Following integer values"<<endl;
		for (auto& s : intVars) out<<(shrt?"":"\t")<<"-"<<s.first<<" (default: "<<s.second<<")"<<endl;
		if (!shrt) out << endl;
	}
	if (boolVars.size()) {
		if (!shrt) out<<"Following boolean values"<<endl;
		for (auto& s : boolVars) out<<(shrt?"":"\t")<<"-"<<s.first<<" (default: "<<(s.second?"true":"false")<<")"<<endl;
		if (!shrt) out << endl;
	}
	if (doubleVars.size()) {
		if (!shrt) out<<"Following double values"<<endl;
		out << std::fixed << setprecision(1);
		for (auto& s : doubleVars) out<<(shrt?"":"\t")<<"-"<<s.first<<" (default: "<<s.second<<")"<<endl;
		if (!shrt) out << endl;
	}
	if (uint64Vars.size()) {
		if (!shrt) out<<"Following uint64 values"<<endl;
		for (auto& s : uint64Vars) out<<(shrt?"":"\t")<<"-"<<s.first<<" (default: "<<s.second<<")"<<endl;
		if (!shrt) out << endl;
	}
	out<<endl;

	out<<"-verb (default: 1)"<<endl;
	if (!shrt) {
		out<<"\tint: [0, 1, 2]"<<endl;
		out<<"\tIf verb is 0 the preprocessor will output less stuff to the standard error."<<endl;
	}
	out<<endl;
}

int isHelp(char* arg) {
	return string(arg) == "-help" || string(arg) == "--help" || string(arg) == "-h" || string(arg) == "--h";
}

int main(int argc, char* argv[]){
	ios_base::sync_with_stdio(0);
	cin.tie(0);

	map<string, int> intVars;
	map<string, bool> boolVars;
	map<string, double> doubleVars;
	map<string, uint64_t> uint64Vars;

	intVars["BBTMS_maxVars"]=200;
	boolVars["hardenInModelSearch"]=false;
	intVars["modelSearchAsearchType"]=1;
	intVars["modelSearchIterLimit"]=-1;
	intVars["MRED_minLitsInClause"]=3;
	intVars["MRED_maxLitsInClause"]=100;
	boolVars["MRED_trimBefore"]=false;
	intVars["MRED_randomizedTries"]=0;
	intVars["LRED_minPartitions"]=1;


	if ((argc > 1 && isHelp(argv[1])) || (argc > 2 && isHelp(argv[2])) || (argc > 3 && isHelp(argv[3]))) {
		printHelp(cout, intVars, boolVars, doubleVars, uint64Vars, false);
		return 0;
	}
	if (argc < 3) {
		printHelp(cout, intVars, boolVars, doubleVars, uint64Vars, true);
		cout<<"Use -help for more detailed information"<<endl;
		return 0;
	}
	auto flags = getFlags(argc, argv);

	if (flags.count("h") || flags.count("help")) {
		printHelp(cout, intVars, boolVars, doubleVars, uint64Vars, true);
	}

	string type(argv[2]);
	assert(type == "solve" || type == "preprocess" || type == "reconstruct");
	string file(argv[1]);
	if (type == "reconstruct") {
		if (!flags.count("mapfile")) {
			cout<<"Give mapfile with -mapfile= flag"<<endl;
			cerr<<"Give mapfile with -mapfile= flag"<<endl;
			return 0;
		}
		string mapFile = flags["mapfile"];
		maxPreprocessor::OutputReader opr;
		ifstream in(file);
		int readStatus = opr.readSolution(in);
		in.close();
		if (readStatus > 0) {
			cout<<"Failed to parse solution"<<endl;
			cerr<<"Failed to parse solution"<<endl;
			return 0;
		}
		if (opr.status == 2) {
			cout<<"s UNSATISFIABLE"<<endl;
		}
		else {
			ifstream mapF(mapFile);
			int vars, ppVars, originalVars;
			mapF>>vars>>ppVars>>originalVars;
			vector<int> varMap(vars);
			for (int i = 0; i < vars; i++) mapF>>varMap[i];
			vector<int> trueLits;
			for (int lit : opr.trueLits) {
				if (abs(lit) > vars) continue;
				if (lit > 0) lit = varMap[abs(lit)-1];
				else lit = -varMap[abs(lit)-1];
				trueLits.push_back(lit);
			}
			maxPreprocessor::Trace trace;
			int traceLines;
			mapF>>traceLines;
			trace.operations.resize(traceLines);
			trace.data.resize(traceLines);
			for (int i = 0; i < traceLines; i++) {
				mapF>>trace.operations[i];
				int sz;
				mapF>>sz;
				trace.data[i].resize(sz);
				for (int j = 0; j < sz; j++) {
					mapF>>trace.data[i][j];
				}
			}
			trace.printSolution(cout, trueLits, opr.ansW, ppVars, originalVars);
		}
		return 0;
	}

	if (type == "solve") {
		// just check that -solver flag is used
		if (!(flags.count("solver") && flags["solver"].size() > 0)) {
			cout<<"Please specify the solver"<<endl;
			cerr<<"Please specify the solver"<<endl;
			return 0;
		}
	}

	string techniques = parseTechniques(flags);
	double timeLimit = parseTimeLimit(flags);
	bool BVEgate = parseBVEgate(flags);
	bool labelMatching = parseLabelMatching(flags);
	bool maxSat = parseProblemType(flags);
	int skipTechnique = parseSkipTechnique(flags);
	bool BVEsortMaxFirst = parseBVEsortMaxFirst(flags);
	int BVElocalGrow = parseBVElocalGrow(flags);
	int BVEglobalGrow = parseBVEglobalGrow(flags);



	pair<unsigned, unsigned> sizeLimit = parseSizeLimit(flags);
	bool ignoreExitCode = flags.count("ignore-exit-code");
	string prepFile = parsePrepFilename(flags);
	string solFile = parseSolFilename(flags);

	parseIntVars(flags, intVars);
	parseBoolVars(flags, boolVars);
	parseDoubleVars(flags, doubleVars);
	parseUint64Vars(flags, uint64Vars);

	ifstream instanceFile(file);
	if (instanceFile.fail()) {
		cout<<"Failed to read the input file"<<endl;
		cerr<<"Failed to read the input file"<<endl;
		return 0;
	}
	maxPreprocessor::InputReader inputReader;
	int readStatus = inputReader.readClauses(instanceFile, maxSat);
	instanceFile.close();

	if (readStatus > 0) {
		cout<<"Failed to parse input instance: "<<inputReader.readError<<endl;
		cerr<<"Failed to parse input instance: "<<inputReader.readError<<endl;
		return 0;
	}

	int outputFormat = maxPreprocessor::INPUT_FORMAT_WPMS22;
	if (flags.count("outputformat")) {
		if (flags["outputformat"] == "original")	outputFormat = inputReader.inputFormat;
		else if (flags["outputformat"] == "wpms")	outputFormat = maxPreprocessor::INPUT_FORMAT_WPMS;
		else if (flags["outputformat"] == "wpms22")	outputFormat = maxPreprocessor::INPUT_FORMAT_WPMS22;
		else if (flags["outputformat"] == "sat")	outputFormat = maxPreprocessor::	INPUT_FORMAT_SAT;
		else {
			cout << "c invalid outputformat value " << flags["outputformat"] << endl;
			cerr << "Invalid outputformat value " << flags["outputformat"] << endl;
		}
		if (outputFormat == maxPreprocessor::INPUT_FORMAT_MS) {
			// preprocessor works in labeled cnf so it cannot output pure maxsat
			outputFormat = maxPreprocessor::INPUT_FORMAT_WPMS;
		}
		string outf;
		if (outputFormat == maxPreprocessor::INPUT_FORMAT_WPMS) {
			outf = "weighted partial Max-SAT (pre 2022)";
		}
		else if (outputFormat == maxPreprocessor::INPUT_FORMAT_SAT) {
			outf = "SAT";
		}
		else if (outputFormat == maxPreprocessor::INPUT_FORMAT_WPMS22) {
			outf = "weighted partial Max-SAT (2022 ->)";
		}
		else {
			return 0;
		}
		cout<<"c Outputformat "<<outf<<endl;
		cerr<<"Outputformat "<<outf<<endl;
	}

	int verb = parseVerb(flags);

	maxPreprocessor::PreprocessorInterface pif(inputReader.clauses, inputReader.weights, inputReader.top);
	pif.setBVEGateExtraction(BVEgate);
	pif.setLabelMatching(labelMatching);
	pif.setSkipTechnique(skipTechnique);
	pif.setBVEsortMaxFirst(BVEsortMaxFirst);
	pif.setBVElocalGrowLimit(BVElocalGrow);
	pif.setBVEglobalGrowLimit(BVEglobalGrow);
	if (intVars["BBTMS_maxVars"]<0) {
		intVars["BBTMS_maxVars"] = (int) ((long long) pif.getOriginalVariables() * (long long)(-intVars["BBTMS_maxVars"]) / 1000000);
	}
	pif.setOptionVariables(intVars, boolVars, doubleVars, uint64Vars);

	maxPreprocessor::Timer preprocessTimer;
	preprocessTimer.start();

	// TODO: better size limiting...
	if (inputReader.clauses.size() < sizeLimit.first || sizeLimit.second < inputReader.clauses.size()) timeLimit=0;
	pif.preprocess(techniques, verb, timeLimit);

	preprocessTimer.stop();
	pif.printPreprocessorStats(cerr);

	cerr<<"Preprocess time: "<<preprocessTimer.getTime().count()<<endl;
	if (verb > 0) pif.printTimeLog(cerr);
	if (type == "preprocess") {
		if (verb > 0) pif.printTechniqueLog(cerr);
		pif.printTechniqueLog(cout);
		pif.printTimeLog(cout);
		pif.printInfoLog(cout);
		pif.printInstance(cout, outputFormat);

		if (flags.count("mapfile")) {
			string mapFile = flags["mapfile"];
			cout<<"c Outputting reconstruction map to "<<mapFile<<endl;
			cerr<<"Outputting reconstruction map to "<<mapFile<<endl;
			ofstream out(mapFile);
			pif.printMap(out);
			out.close();
		}
		else {
			cout<<"c No -mapfile= given, will not ouput reconstruction map"<<endl;
			cerr<<"No -mapfile= given, will not ouput reconstruction map"<<endl;
		}
	}
	if (type == "solve") {
		string solver;
		if (flags.count("solver") && flags["solver"].size() > 0) {
			solver = flags["solver"];
			cout<<"c Using solver "<<solver<<endl;
			cerr<<"Using solver "<<solver<<endl;
			cout<<"c Solver flags "<<flags["solverflags"]<<endl;
			cerr<<"Solver flags "<<flags["solverflags"]<<endl;
		}
		else {
			cout<<"Please specify the solver"<<endl;
			cerr<<"Please specify the solver"<<endl;
			return 0;
		}

		cout << "c Saving preprocessed instance into file " << prepFile << endl;
		cerr << "Saving preprocessed instance into file " << prepFile << endl;
		ofstream out(prepFile);
		pif.printInstance(out, outputFormat);
		out.close();

		string command = (solver + " " + prepFile +" " + flags["solverflags"] + " > " + solFile);
		cout << "c Invoking solver... command: "  << command << endl;
		cerr << "Invoking solver... command: " << command << endl;

		maxPreprocessor::Timer solveTimer;
		solveTimer.start();
		int rv = system(command.c_str());
		int exit_status = (WEXITSTATUS(rv));
		if (exit_status && !ignoreExitCode) {
			cout << "Solver error, returned nonzero value (" << exit_status << ")" << endl;
			cerr << "Solver error, returned nonzero value (" << exit_status << ")" << endl;
			return 0;
		}
		solveTimer.stop();

		maxPreprocessor::OutputReader opr;
		ifstream in(solFile);
		readStatus = opr.readSolution(in);
		in.close();
		if (readStatus > 0) {
			cout<<"Failed to parse solution"<<endl;
			cerr<<"Failed to parse solution"<<endl;
			return 0;
		}

		if (opr.status == 2) {
			cout<<"s UNSATISFIABLE"<<endl;
		}
		else {
			pif.printSolution(opr.trueLits, cout, opr.ansW);
		}
		cerr<<"Preprocess time: "<<preprocessTimer.getTime().count()<<", Solve time: "<<solveTimer.getTime().count()<<endl;
		if (verb > 0) pif.printTimeLog(cerr);
	}
}
