#ifndef MAXPP_PREPROCESSOR_HPP
#define MAXPP_PREPROCESSOR_HPP

#include <vector>
#include <queue>
#include <set>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "global.hpp"
#include "preprocessedinstance.hpp"
#include "probleminstance.hpp"
#include "trace.hpp"
#include "timer.hpp"
#include "clause.hpp"
#include "log.hpp"
#include "AMSLEX.hpp"
#include "satsolver/satsolverinterface.hpp"
#include "satsolver/glucose3.cpp"

namespace maxPreprocessor {
// Variables are indexed 0..n-1
// Positive literal is v*2, negative is v*2+1
class Preprocessor {
public:
	// parameters, magic constants etc..
	struct Options {
		int skipTechnique;
		int skipSkipConstant;
		bool BVEgate;
		bool BVEsortMaxFirst;

		bool hardenInModelSearch;
		int modelSearchIterLimit;
		int modelSearchAsearchType; // 0: linear seatch, 1: binary search, -1: vary

		int BIG_tries;
		int BVE_sizelimit;
		int GSLE_tries;
		int MRED_minLitsInClause, MRED_maxLitsInClause;
		bool MRED_trimBefore;
		int MRED_randomizedTries;
		int SE_Lim, SE_HashLim, SE_AmsLexLim, SE_HashAmsLexLim2;
		uint64_t SSR_Lim; int SSR_HashLim, SSR_AmsLexLim, SSR_HashAmsLexLim2;
		int BBTMS_maxVars;
		int HARD_asearchIterLimit; double HARD_improveTimeLimit; int HARD_satLikeTries; double HARD_satLikeTimeLimit;
		int MRED_asearchIterLimit; double MRED_improveTimeLimit; int MRED_satLikeTries; double MRED_satLikeTimeLimit;
		int LRED_minPartitions, LRED_maxUNSATReductChecksPerLabel;
		Options() :
				skipTechnique(0), skipSkipConstant(4), BVEgate(true), BVEsortMaxFirst(false),
				hardenInModelSearch(false), modelSearchIterLimit(-1), modelSearchAsearchType(1),
				BIG_tries(9), BVE_sizelimit(6), GSLE_tries(10),
				MRED_minLitsInClause(3), MRED_maxLitsInClause(100), MRED_trimBefore(false), MRED_randomizedTries(0),
				SE_Lim(64), SE_HashLim(10), SE_AmsLexLim(30), SE_HashAmsLexLim2(10000),
				SSR_Lim(50000), SSR_HashLim(4), SSR_AmsLexLim(8), SSR_HashAmsLexLim2(10000),
				BBTMS_maxVars(200),
				HARD_asearchIterLimit(-1), HARD_improveTimeLimit(1), HARD_satLikeTries(0), HARD_satLikeTimeLimit(0),
				MRED_asearchIterLimit(1), MRED_improveTimeLimit(0), MRED_satLikeTries(0), MRED_satLikeTimeLimit(0),
				LRED_minPartitions(1), LRED_maxUNSATReductChecksPerLabel(20)
		{}

		void parseValues(map<string, int>& intVars, map<string, bool>& boolVars, map<string, double>& doubleVars, map<string, uint64_t>& uint64Vars) {
			map<string, int>::iterator intIt;
			map<string, bool>::iterator boolIt;
			map<string, double>::iterator doubleIt;
			map<string, uint64_t>::iterator uint64It;

			if ((intIt=intVars.find("skipTechnique")) != intVars.end()) { skipTechnique = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("skipSkipConstant")) != intVars.end()) { skipSkipConstant = intIt->second; intVars.erase(intIt); }
			if ((boolIt=boolVars.find("BVEgate")) != boolVars.end()) { BVEgate = boolIt->second; boolVars.erase(boolIt); }
			if ((boolIt=boolVars.find("BVEsortMaxFirst")) != boolVars.end()) { BVEsortMaxFirst = boolIt->second; boolVars.erase(boolIt); }
			if ((boolIt=boolVars.find("hardenInModelSearch")) != boolVars.end()) { hardenInModelSearch = boolIt->second; boolVars.erase(boolIt); }
			if ((intIt=intVars.find("modelSearchIterLimit")) != intVars.end()) { modelSearchIterLimit = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("modelSearchAsearchType")) != intVars.end()) { modelSearchAsearchType = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("BIG_tries")) != intVars.end()) { BIG_tries = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("BVE_sizelimit")) != intVars.end()) { BVE_sizelimit = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("GSLE_tries")) != intVars.end()) { GSLE_tries = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("MRED_minLitsInClause")) != intVars.end()) { MRED_minLitsInClause = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("MRED_maxLitsInClause")) != intVars.end()) { MRED_maxLitsInClause = intIt->second; intVars.erase(intIt); }
			if ((boolIt=boolVars.find("MRED_trimBefore")) != boolVars.end()) { MRED_trimBefore = boolIt->second; boolVars.erase(boolIt); }
			if ((intIt=intVars.find("MRED_randomizedTries")) != intVars.end()) { MRED_randomizedTries = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("SE_Lim")) != intVars.end()) { SE_Lim = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("SE_HashLim")) != intVars.end()) { SE_HashLim = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("SE_AmsLexLim")) != intVars.end()) { SE_AmsLexLim = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("SE_HashAmsLexLim2")) != intVars.end()) { SE_HashAmsLexLim2 = intIt->second; intVars.erase(intIt); }
			if ((uint64It=uint64Vars.find("SSR_Lim")) != uint64Vars.end()) { SSR_Lim = uint64It->second; uint64Vars.erase(uint64It); }
			if ((intIt=intVars.find("SSR_HashLim")) != intVars.end()) { SSR_HashLim = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("SSR_AmsLexLim")) != intVars.end()) { SSR_AmsLexLim = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("SSR_HashAmsLexLim2")) != intVars.end()) { SSR_HashAmsLexLim2 = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("BBTMS_maxVars")) != intVars.end()) { BBTMS_maxVars = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("HARD_asearchIterLimit")) != intVars.end()) { HARD_asearchIterLimit = intIt->second; intVars.erase(intIt); }
			if ((doubleIt=doubleVars.find("HARD_improveTimeLimit")) != doubleVars.end()) { HARD_improveTimeLimit = doubleIt->second; doubleVars.erase(doubleIt); }
			if ((intIt=intVars.find("HARD_satLikeTries")) != intVars.end()) { HARD_satLikeTries = intIt->second; intVars.erase(intIt); }
			if ((doubleIt=doubleVars.find("HARD_satLikeTimeLimit")) != doubleVars.end()) { HARD_satLikeTimeLimit = doubleIt->second; doubleVars.erase(doubleIt); }
			if ((intIt=intVars.find("MRED_asearchIterLimit")) != intVars.end()) { MRED_asearchIterLimit = intIt->second; intVars.erase(intIt); }
			if ((doubleIt=doubleVars.find("MRED_improveTimeLimit")) != doubleVars.end()) { MRED_improveTimeLimit = doubleIt->second; doubleVars.erase(doubleIt); }
			if ((intIt=intVars.find("MRED_satLikeTries")) != intVars.end()) { MRED_satLikeTries = intIt->second; intVars.erase(intIt); }
			if ((doubleIt=doubleVars.find("MRED_satLikeTimeLimit")) != doubleVars.end()) { MRED_satLikeTimeLimit = doubleIt->second; doubleVars.erase(doubleIt); }
			if ((intIt=intVars.find("LRED_minPartitions")) != intVars.end()) { LRED_minPartitions = intIt->second; intVars.erase(intIt); }
			if ((intIt=intVars.find("LRED_maxUNSATReductChecksPerLabel")) != intVars.end()) { LRED_maxUNSATReductChecksPerLabel = intIt->second; intVars.erase(intIt); }

		}
	};



	Options opt;

	ProblemInstance pi;
	int originalClauses;
	int originalVars;

	Log rLog;

	int logLevel = 1;
	bool printComments = true;

	Trace trace;

	SATSolver* satSolver;
	void prepareSatSolver();


	Preprocessor(const std::vector<std::vector<int> >& clauses_, const std::vector<uint64_t>& weights_, uint64_t topWeight_);

	bool isTautology(const Clause& clause) const;

	// Returns number of clauses removed
	int setVariable(int var, bool value);

	// This is called only in the beginning since no tautologies are added
	void removeTautologies();

	int eliminateReduntantLabels();
	// This is called only in the beginning
	void identifyLabels();

	// This is called only in the beginning
	void createLabels();

	const uint64_t polyHashMul = 1000000007; // magic constant

	int tryUP(int lit);
	int tryUPAll();

	int doUP();
	void doUP2();

	int removeDuplicateClauses();

	// Is clause a subsumed by clause b?
	// Supposes that a and b are sorted
	bool isSubsumed(const std::vector<int>& a, const std::vector<int>& b) const;

	AMSLEX amsLex;
	void trySEHash(std::vector<int>& clauses, int tLit, std::vector<int>& toRemove);
	void trySEAmsLex(std::vector<int>& clauses, std::vector<int>& toRemove);
	int trySESlow(int lit);
	void trySE(std::vector<int>& clauses, std::vector<int>& toRemove);
	void trySEgen(int lit, std::vector<int>& toRemove);

	int doSE();
	void doSE2();

	int BVElocalGrow;
	int BVEglobalGrow;

	void printC(int c) const;
	std::pair<std::vector<int>, int> searchXor(int var) const;
	std::pair<std::vector<int>, int> searchITE(int var) const;
	std::pair<std::vector<int>, int> searchAndOr(int dLit, int defVars, const std::set<int>& binaryClauseLits) const;
	// Dont give labels to this
	int tryBVE(int var);
	int tryBVE2(int var);
	std::vector<int> tryBVEGE(int var);

	std::vector<uint64_t> getBVEHash(const std::vector<int>& cs, int var, int sw) const;

	std::unordered_map<std::string, double> stats;

	int doBVE();
	void doBVE2();


	// Supposes that the clauses are in sorted order
	// Can we do SSR such that var is removed from c2?
	bool canSSR(int var, const Clause& c1, const Clause& c2);


	bool SSRC(int c1, int c2, int var);
	int trySSRAmsLex(int var);
	int trySSRHash(int var);
	int trySSR2(int var);
	int trySSR(int var);
	int trySSRgen(int var);

	int doSSR();
	void doSSR2();

	int tryBCE(int lit);

	int doBCE();
	void doBCE2();

	bool vSubsumed(std::vector<int>& v1, std::vector<int>& v2);
	int trySLESlow(int lb1, int lb2);
	int doSLE();
	void doSLE2();

	bool tryBCR(int c, int l11);
	int doBCR();
	void doBCR2();

	bool SIErndCheck(int litX, int litY);
	int try2SIE(int litX, int litY);
	int trySIE(int lit);
	int doSIE();
	void doSIE2();



	void genIndex(std::vector<std::vector<int> >& g, std::vector<std::vector<int> >& rg, int x, int u1, int u2, int& stamp, std::vector<int>& le, std::vector<int>& ri, std::vector<int>& up1, std::vector<int>& up2, int order);
	int BIGIt;
	std::vector<int> BIGu, BIGu2, BIGid;
	void BIGdfs1(int x, std::vector<int>& ns);
	void BIGdfs2(std::vector<std::vector<int> >& g, int x, std::vector<int>& ns);
	void BIGdfs3(std::vector<std::vector<int> >& rg, int x, std::vector<int>& scc);
	bool BIGisPath(int x, int to, std::vector<int>& leIndex, std::vector<int>& riIndex, std::vector<int>& up1, std::vector<int>& up2);
	int tryBIG(int lit, bool doTC);
	int doBIG(bool doTC);
	void doBIG2(bool doTC);
	bool doneUnhiding;

	std::vector<uint64_t> sfH;
	std::vector<uint64_t> tMul;
	std::unordered_map<uint64_t, int> BVAHashTable;
	void addBVAHash(std::vector<int>& lits, std::unordered_map<uint64_t, int>& hashes);
	int canBVA(int c, int d, int lit);
	int tryBVA(int lit, std::unordered_map<uint64_t, int>& hashes);
	int doBVA();
	void doBVA2();

	void GSLEBT(int i, uint64_t w, std::vector<int>& sel, std::vector<uint64_t>& weights, std::vector<std::vector<int> >& hs, bool& found, uint64_t& itLim);
	bool GSLEtryBackTrack(std::vector<std::vector<int> >& hs, std::vector<uint64_t>& weights, uint64_t w, uint64_t itLim);
	int tryGSLE(int lb);
	int doGSLE();
	void doGSLE2();

	int tryFLP(std::vector<int> fLit, int clause);
	int doFLP();

	void tryLFF(int lb);
	void findLabeledFormula();

	int tryLS(int lbl);
	void tryLSBCE(int lit, std::unordered_set<int>& deletedClauses, std::unordered_set<int>& touchedList, std::vector<std::pair<int, int> >& blockedClauses);
	int doLS();

	int tryAM1(vector<int>& vars, bool weight_aware, bool stratification, bool greedy_cost);
	int doAM1(bool weight_aware, bool stratification, bool greedy_cost);

	unordered_set<int> canSatLits;
	int tryTMS(vector<int>& vars);
	int TMSMaxVars; // tmp TODO: test and remove
	int doBBTMS(int maxVars);
	int doTMS();

	int redSatSolverCalls; // statistics
	bool checkPositiveReduct(const vector<int>& literals, const vector<int>& assumptions);
	bool checkExtendedPositiveReduct(const vector<int>& literals, const vector<int>& assumptions);
	bool checkFilteredPositiveReduct(vector<int>& literals, const vector<int>& assumptions, bool f2=false);
	bool checkTrimmedFilteredPositiveReduct(const vector<int>& literals, const vector<int>& assumptions, bool f2=false);
	int trimReductClause(vector<int>& clause);
	int findREDPartitionForLit(int lit, unsigned n, vector<pair<vector<pair<uint64_t, int> >, vector<int> > >& clauses, uint64_t mcost);
	int tryREDOnLit(int lit);
	int tryREDOnClause(int c);
	int tryModelCuttingRED();
	int tryUPLitRED(int lit);
	int doLabelRED();
	int doClauseRED();
	int modelCuttingNewLit(vector<int>& clause, vector<int>& bclause, vector<int>& up, int litsInClause, int lit);
	int doModelCuttingRED();
	int doUPLitRED();


	vector<bool> bestModel;
	uint64_t bestCost;
	pair<uint64_t, uint64_t> modelCostCheck(vector<bool>& model);
	int updateBestModel(vector<bool>& best_model, uint64_t& best_cost, vector<bool>& model);
	void findGoodModelSub(vector<pair<uint64_t, int> >& labels, int m, vector<int>& idx, vector<int> base_assumps, uint64_t& best_cost, vector<bool>& best_model, int iters_left, double leaveTime, bool binarySearch);
	uint64_t findGoodModelA(vector<pair<uint64_t, int> >& labels, vector<bool>& model, int iterLimit, double leaveTime);
	uint64_t findGoodModel(vector<bool>& model, int assumpsSearchIterLimit, double timeImprove, int satLikeTries, double timeSatLike);


	int tryHARD();
	int doHARD();

	int flePos;
	void replaceLit(int lit1, int lit2);
	void handleEqLits(vector<int>& lits);
	int tryFLE(int lit, vector<int>& up, bool doRLE);
	int tryFLE(int var, bool doRLE, bool findEqs, bool findRedEqs, bool findForced, bool findRedForced);
	int doFLE(bool doRLE, bool findEqs, bool findRedEqs, bool findForced, bool findRedForced);

	void CBIGdfs1(int x, std::vector<int>& ns, std::vector<std::pair<int, std::pair<int, int> > >& condEdges);
	int findConditionalGraph(int lit, std::vector<std::pair<int, std::pair<int, int> > >& condEdges);
	std::vector<std::pair<int, std::pair<int, int> > > findConditionalComponents();

	int removeEmptyClauses();

	bool validTechniques(std::string techniques) const;
	bool validPreTechniques(std::string techniques) const;

	PreprocessedInstance getPreprocessedInstance();

	int doPreprocess(const std::string& techniques, int l, int r, bool debug, bool topLevel);
	void preprocess(std::string techniques, double timeLimit, bool debug, bool BVEgate, bool initialCall, bool matchLabels);


	std::mt19937 randGen;

	template<typename T>
	void log(T t) {
		if (logLevel < 2) return;
		std::cerr << t << std::endl;
	}

	template<typename T, typename... Args>
	void log(T t, Args... args) {
		if (logLevel < 2) return;
		std::cerr << t;
		log(args...);
	}

	template<typename T>
	void print(T t) {
		if (!printComments) return;
		std::cout << t << std::endl;
	}

	template<typename T, typename... Args>
	void print(T t, Args... args) {
		if (!printComments) return;
		std::cout << t;
		print(args...);
	}

	template<typename T>
	T getRand(T lo, T hi) {
		return std::uniform_int_distribution<T>(lo, hi)(randGen);
	}

	void printStats(std::ostream& o, std::string b = "c MAXPRE-STATS ") {
		for (auto& t : stats) {
			o << b << t.first << " = " << t.second << "\n";
		}
		o.flush();
	}

	static std::string version(int l = 0);
};
}
#endif
