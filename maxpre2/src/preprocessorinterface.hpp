#ifndef MAXPP_PREPROCESSORINSTANCE_HPP
#define MAXPP_PREPROCESSORINSTANCE_HPP

#include <cstdint>
#include <iostream>

#include "preprocessor.hpp"

namespace maxPreprocessor {
class PreprocessorInterface {
private:
	int variables;
	int originalVariables;
	Preprocessor preprocessor;
	PreprocessedInstance preprocessedInstance;
	uint64_t topWeight;
	bool inProcessMode;
	bool preprocessed;
	bool useBVEGateExtraction;
	bool useLabelMatching;
	int BVElocalGrow;
	int BVEglobalGrow;

	Preprocessor::Options opt;

	std::vector<int> solverVarToPPVar;
	std::vector<int> PPVarToSolverVar;
	int litToSolver(int lit);
	int litToPP(int lit);
public:
	/* clauses should contain the clauses of the maxsat instance. The variables
	 * are indexed with positive integers, x is the positive literal of variable
	 * x and -x is the negative literal of variable x.
	 * weights are the weights of the clauses. The sizes of clauses and weights
	 * vectors shoudl be equal. topWeight paremeter gives the topweigth.
	 * Setting inProcessMode true enables the inprocessing mode. It disables
	 * compression of the variable indexes and enables the functionality to
	 * alter the instance.
	 */
	PreprocessorInterface(const std::vector<std::vector<int> >& clauses, const std::vector<uint64_t>& weights, uint64_t topWeight_, bool inProcessMode_ = false);

	/* Preprocesses the current maxsat instance with the given techniques
	 * string, loglevel and timelimit.
	 */
	void preprocess(std::string techniques, int logLevel = 0, double timeLimit = 1e9);

	// Returns the topweight. This should be the same as given in the constructor.
	uint64_t getTopWeight();


	/* These functions work only if inprocessing mode is enabled. Use them only
	 * if you know what you are doing. Correctness of some preprocessing
	 * techniques is not preseved when adding arbitrary clauses.
	 */
	/* Adds a variable with index var to the instance. If var = 0, the added
	 * variable will have the next free index. Returns the index of the variable
	 * if successful and 0 if not.
	 */
	int addVar(int var = 0);
	/* Adds a clause to the instace. Unknown variables will be added with addVar
	 * function. Return true if successful.
	 */
	bool addClause(std::vector<int> clause);
	/* Adds a variable that is a label with the given weight to the instance.
	 * The negative literal will be in an unit clause with the given weight.
	 * Works similarly to the addVar.
	 */
	int addLabel(int lbl, uint64_t weight);
	/* Changes the weight of a label. lbl is the index of the label variable
	 * i.e. its always positive. If weight >= topweight, the label will become
	 * a normal variable. Normal variables can never be changed into labels.
	 */
	bool alterWeight(int lbl, uint64_t weight);
	/* Changes a label into a variable and deletes the soft clause. lbl is the
	 * index of the label variable i.e. its always positive.
	 */
	bool labelToVar(int lbl);
	// resets removed weight
	bool resetRemovedWeight();
	uint64_t getRemovedWeight();

	// Functions for enabling/disabling some functionality
	void setBVEGateExtraction(bool use);
	void setLabelMatching(bool use);
	// Set value=0 to disable skiptechnique
	void setSkipTechnique(int value);

	// Functions for controlling BVE for some experiments
	void setBVEsortMaxFirst(bool use);
	void setBVElocalGrowLimit(int limit);
	// This limit is per one call of preprocess
	void setBVEglobalGrowLimit(int limit);

	// functions for controlling back bone computing tms for some experiments
	void setMaxBBTMSVars(int limit);

	void setHardenInModelSearch(bool harden);
	void setModelSearchIterLimit(int limit);
	// lazy way of setting option variables...
	void setOptionVariables(map<string, int>& intVars, map<string, bool>& boolVars, map<string, double>& doubleVars, map<string, uint64_t>& uint64Vars);

	int getOriginalVariables() { return originalVariables;}

	uint64_t getUpperBound() { return preprocessor.bestCost; }

	// Returns the current instance
	void _getInstance(std::vector<std::vector<int> >& retClauses, std::vector<uint64_t> & retWeights); // doesn't convert lits
	void getInstance(std::vector<std::vector<int> >& retClauses, std::vector<uint64_t>& retWeights, std::vector<int>& retLabels);
	/* Returns the assignment of the original variables given the assignment of
	 * variable in the solution of the preprocessed istance
	 */
	std::vector<int> reconstruct(const std::vector<int>& trueLiterals, bool convertLits = 1);
	std::vector<std::pair<int, std::pair<int, int> > > getCondEdges();

	void printInstance(std::ostream& output, int outputFormat = 0);
	void printSolution(const std::vector<int>& trueLiterals, std::ostream& output, uint64_t ansWeight);
	void printMap(std::ostream& output);
	void printTechniqueLog(std::ostream& output);
	void printTimeLog(std::ostream& output);
	void printInfoLog(std::ostream& output);
	void printPreprocessorStats(std::ostream& output);

	static std::string version(int l = 0);
};
}
#endif
