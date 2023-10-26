
/********************************************************************************************
dgpw.h -- Copyright (c) 2014-2016, Sven Reimer

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************/

#ifndef DGPW_H_
#define DGPW_H_

#include <vector>
// Include standard headers.
#include "Settings.h"
#include "timemeasurement.h"
//#include "Softclause.h"


class SATSolverProxy;
namespace Pacose {
class Pacose;
struct SoftClause;

namespace DGPW {

// Some definitions.
#define UNKNOWN 0
#define SATISFIABLE 10
#define UNSAT 20

class Sorter;
class Cascade;
class MultipleCascade;
struct TimeVariables;

// The "MaxAntom" class.
class DGPW {
 public:
  // Constructor.
  explicit DGPW(Pacose *pacose);

  // Destructor.
  ~DGPW();

  // Creates and returns a fresh variable index not used so far.
  uint32_t NewVariable(void);

  bool AddClause(std::vector<uint32_t> &clause, uint32_t lbd = 1);
  bool AddUnit(uint32_t lit);

  uint32_t Solve(void);
  uint32_t Solve(std::vector<uint32_t> &assumptions);
  std::vector<uint32_t> Model(void) const;
  uint32_t Model(uint32_t var) const;

  uint32_t Variables(void) const;
  // Returns the current number of clauses within the clause database.
  uint32_t Clauses(void) const;
  uint32_t StaticClauses(void) const;

  uint32_t CurrentBinaryClauses() const;
  uint32_t CurrentTernaryClauses() const;

  //    void EncodeOR(uint32_t output, uint32_t input1, uint32_t input2);
  //    void EncodeAND(uint32_t output, uint32_t input1, uint32_t input2);

  Settings *GetSetting();
  bool HasEncoding();
  uint64_t GetMinWeight();
  uint64_t GetMaxWeight();
  uint64_t GetSATWeight();
  uint32_t GetEncodingClauses() { return _addedClauses; };
  uint32_t GetEncodingVariables() { return _addedVariables; };
  uint32_t GetNoUnsatisfiableSCs() { return _unsatisfiableSCs; };
  uint32_t GetGreedyPrepro() { return _greedyPrepro; };
  uint32_t GetApproxSolverCalls() { return _approxSolverCalls; };
  uint32_t GetSolverCalls() { return _solverCalls; };
  std::vector<std::pair<uint64_t, uint32_t>> GetTareVector(uint64_t weightDiff);
  std::vector<std::pair<uint64_t, uint32_t>> GetWatchdogs(uint64_t weightDiff);
  bool GetHasHardClauses();
  bool GetHasMoreThanTwoWeights();
  uint32_t GetLastResult();
  std::vector<uint32_t> GetLastSatisfiableAssignment();

  // Add weighted SoftClauses to _softClauseVector
  bool AddWeightedSoftClause(std::vector<uint32_t> &clause,
                             uint64_t weight = 1);

  // Add whole _softclause Vector to clause database
  //    void AddSCfromGlucose(uint64_t weight, const std::vector<u_int32_t>&
  //    literals, uint32_t antBlock);

  // Solve weighted partial maxSAT Problem
  uint32_t MaxSolveWeightedPartial(int64_t &optimum);
  uint32_t MaxSolveWeightedPartial(std::vector<uint32_t> &externalAssumptions,
                                   int64_t &optimum);

  /* Some maxSAT related interface functions */
  void SetSoftClauseVector(std::vector<SoftClause *> *softClauses);

  // Weighted MaxSat Stuff
  StructureInfo AnalyzeandConvertStructure();

  void SetGreatestCommonDivisor(uint64_t val);
  void SetMoreThanTwoWeights(bool val);
  void SetTopWeight(uint64_t val);
  void SetMinWeight(uint64_t val);
  void SetMaxWeight(uint64_t val);
  void SetHasHardClauses(bool val);
  void SetSatWeight(uint64_t val);

  void SetInitialAssumptions(std::vector<uint32_t> externalAssumptions);

  friend class Cascade;
  friend class MultipleCascade;
  friend class Bucket;
  friend class Sorter;
  friend class Counter;

  uint64_t _sumOfSoftWeights;

  uint64_t CalcGlobalOpt();
  std::vector<uint32_t> GetLastAssumptions();
  void SetFixedAssumptions(std::vector<uint32_t>);
  void RemoveFixedAssumptions();
  uint32_t GreedyMaxInitSATWeight(int greedyPrepro = 1, uint32_t maxTime = 5,
                                  uint32_t maxRounds = 500);

  TimeVariables *_timeVariables;

  uint32_t MaxSolveIncremental();

 private:
  // Copy constructor.
  DGPW(const DGPW &) = default;

  uint64_t CountSatisfiedSoftClauses(Sorter *sorter,
                                     const std::vector<uint32_t> &model);
  uint64_t CountSatisfiedSoftClauses(std::vector<SoftClause *> softclauses,
                                     const std::vector<uint32_t> &model);

  // void mergePartTrigger ( Sorter& trigger1, Sorter& trigger2 );
  void CheckAllConflictingSoftclauses(void);

  // Cout the model of the buckets and which position has the first set bit.
  void DumpBucketModel(const std::vector<uint32_t> &model);

  // Called from Cascade and Bucket - they only calculate their own optimum.
  uint64_t CalculateOverallOptimum(uint64_t SatWeight, bool countAgain);

  std::vector<SoftClause *> _softClauses;

  std::vector<std::vector<Sorter *>> _sorterTree;

  SATSolverProxy *_solver;
  Settings *_dgpwSetting;
  Pacose *_pacose;

  uint32_t _maxSorterDepth;

  uint32_t _maxsatResult;

  // needed for updating optimum value from Weighted Mode - Cascade and Bucket
  int64_t *_optimum;
  // Weighted MaxSAT stuff

  Cascade *_mainCascade;
  MultipleCascade *_mainMultipleCascade;
  int64_t _greatestCommonDivisor;
  StructureInfo AnalyzeStructure();
  // StructureInformation AnalyzeStructure();
  void ConvertFormulaToMaxSAT(uint64_t maxWeight);
  void DivideAllSoftClausesByFactor(uint64_t factor);
  uint64_t GreatestCommonDivisor(uint64_t a, uint64_t b);

  uint32_t _currentBucketForTare;
  uint64_t _satWeight;
  bool _moreThanTwoWeights;
  uint64_t _topWeight;
  uint64_t _minWeight;
  uint64_t _maxWeight;
  bool _hasEncoding;
  bool _hasHardClauses;

  // if last solver call returned ANTOMUNKNOWN - this variable is set to true
  bool _resultUnknown;

  // maybe later on as local vars again - just for testing!!
  uint32_t _clausesBefore;
  uint32_t _variablesBefore;
  uint32_t _binaryClausesBefore;
  uint32_t _ternaryClausesBefore;
  uint32_t _addedClauses;
  uint32_t _addedVariables;
  uint32_t _addedBinaryClauses;
  uint32_t _addedTernaryClauses;

  // for coding clause estimation in standard mode
  double _binaryTopClauses;
  uint64_t _ternaryTopClauses;
  uint64_t _binaryBottomClauses;
  uint64_t _ternaryBottomClauses;

  uint32_t _solverCalls;
  uint32_t _approxSolverCalls;

  // partial status variables
  uint64_t _highestBucketMultiplicator;
  uint64_t _currentBucketMultiplicator;
  // the value calculated highestBucketMultiplicator * bucketInfo[1];
  uint64_t _estimatedSatWeight;
  uint64_t _diffToSatWeight;

  //    std::vector<uint32_t> _collectedAssumptions;
  std::vector<uint32_t> _externalAssumptions;
  std::vector<uint32_t> _fixedAssumptions;
  uint32_t _lastResult;

  // Store last model
  std::vector<uint32_t> _lastModel;
  uint64_t _preproPropagationLimit = 0;

  ///
  /// \brief SatisfiedSCsInfo
  /// \param sortedSCIndices
  /// \return tuple < vector of all relaxLiterals which are satisfied, can be
  /// used directly for assumptions
  ///                 vector of all SoftCaluses which are not satisfied
  ///                 actual satisfied weight >
  std::tuple<std::vector<uint32_t>, std::vector<SoftClause *>, uint64_t>
  SatisfiedSCsInfo(std::vector<uint32_t> *sortedSCIndices);

  uint32_t BinarySearchSatisfySCs(std::vector<uint32_t> &nextAssumptions,
                                  std::vector<SoftClause *> *UNSATSCs);
  double _previousPosition;
  double _previousLowerBound;
  double _noClauses;  // number of clauses added, each clause forcing at least
                      // one SC to be true
  double _noVars;  // number of vars in one clause, one SC out of noVars has to
                   // be satisfied.
  uint32_t _unsatisfiableSCs;
  int _greedyPrepro;

  // reset for incremental MaxSAT Solving
  void IncrementalReset();
};
}  // namespace DGPW
} // Namespace Pacose
#endif  // DGPW_H
