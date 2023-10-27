/********************************************************************************************
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************************************/

#ifndef PACOSE_H
#define PACOSE_H

#include <iostream>
#include <vector>
#include <tuple>
#include "Settings.h"
#include "../VeriPB_Prooflogger/VeriPBProoflogger.h"
#include "../VeriPB_Prooflogger/MaxSATProoflogger.h"
#include "../VeriPB_Prooflogger/cadicalprooftracer.hpp"




class SATSolverProxy;
enum class SATSolverType;

namespace Pacose {
// #define SaveCNF

struct ClauseDB;
struct SoftClause;
class Encodings;
namespace DGPW
{
  class DGPW;
}

// Pacose is an extension to the glucose solver
// containing
/**
 * @brief The Pacose class
 *      1.  (Partial) (Weighted) MaxSAT Solver
 *          Initially derived from QMaxSAT (2017 Competition Version)
 */
class Pacose
{
public:
  Pacose();

  ~Pacose();

  uint32_t SolveProcedure(ClauseDB& clauseDB);
  bool ExternalPreprocessing(ClauseDB& clauseDB, VeriPbProofLogger &vPL, MaxSATProoflogger &mPL);
  // void CallMaxPre2(ClauseDB &clauseDB);

  uint32_t CalculateNextResult();
  uint32_t CalculateNextSoftclauseCombination();
  
  uint32_t SolveQMax(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                     EncodingType *encodingType = nullptr);
  
  void Preprocess();
  void ExcludeCurrentSoftclauseCombination();
  void ExcludeCurrentResult();
  

  // Variables
  Settings _settings;
  //    Encodings _encodings;

  SATSolverProxy *_satSolver;
  uint32_t _nbVars;
  uint32_t _nbClauses;
  uint32_t _nbOriginalClauses;
  uint64_t _top; // hard clause weight
  std::vector<SoftClause *> _originalSoftClauses;
  std::vector<SoftClause *> *_actualSoftClauses;
  std::vector<std::vector<SoftClause *>> _sClauses;

  uint64_t CalculateSATWeight();
  long long int GetSATWeight()
  {
    if (_satWeight < 0 || _satWeight > _sumOfSoftWeights)
    {
      _satWeight = _sumOfSoftWeights - _unSatWeight;
    }
    return _satWeight;
  };
  long long int GetOValue()
  {
    return _unSatWeight;
  };

  /**
   * @brief InitSatSolver
   * @param solver
   *          0 == GLUCOSE421
   *          1 == GLUCOSE3
   *          2 == CADICAL
   *          3 == CRYPTOMINISAT
   */
  void InitSatSolver(int solver = 0);
  void InitSatSolver(SATSolverType solverType);
  void AddSoftClause(std::vector<uint32_t> &clause, MaxSATProoflogger &mPL, VeriPbProofLogger& vPL, std::vector<std::tuple<uint64_t, uint32_t, uint32_t>>& unitsoftclauses,  uint64_t weight = 1);

  /**
   * @brief AddNextCNF from the clause vector - for incremental CNF with main
   * @param number    0 means main CNF including Soft Clauses
   *                  i stands for CNF i
   */
  void AddNextCNF(int number);

  /**
   * @brief AddUpcomingClausesAsAssumptions
   * @return bool    true  - upcoming clauses are added temporarily
   *                 false - upcoming clauses are added as hard clauses
   */
  bool AddUpcomingClausesAsAssumptions();


  bool _hasHardClauses;
  void SetSumOfSoftWeights(uint64_t softWeights);

  // literal given as 2 x int (+1)
  void AddAssumption(uint32_t literal);
  void ClearAssumptions();
  void DeactivateLastCNF();

  bool AddClause(std::vector<uint32_t> &clause);

  uint32_t NewVariable();
  void NewVariables(uint32_t noVars);
  uint32_t GetModel(int var);

  uint32_t Solve();

  bool DumpWCNF();

private:
  struct partitionInformation
  {
    partitionInformation() : dgpw(nullptr) {}

    DGPW::DGPW *dgpw;
    uint32_t Points;
    uint64_t weightsTillPoint;
    uint64_t ggtTillPoint;
    bool allWeightsAreEqual;
  };

  // maxPreprocessor::PreprocessorInterface* maxpre;

  std::vector<partitionInformation> _cascCandidates;
  EncodingType _encoding;

  // Variables
  std::vector<std::vector<uint32_t>> _CNF;


  // std::vector<std::vector<std::vector<uint32_t>>> _incrementalCNF;
  int _cpuLimit;
  int _memLimit;
  int _nbOfOrigVars;
  uint64_t  _sumOfSoftWeights;
  uint64_t _overallSoftWeights;
  // for incremental MaxSAT
  // bool _withAssumption;
  // std::vector<uint32_t> _assumptions;
  // uint32_t _relaxVar;
  // fulfilled softclauses
  uint64_t _satWeight;
  // o-value
  uint64_t _unSatWeight;
  uint64_t _lastCalculatedUnsatWeight;
  uint64_t _localUnSatWeight;
  uint64_t _localSatWeight;
  uint64_t _minWeight;
  uint64_t _maxWeight;
  uint64_t _GCD;
  // uint32_t _incrementalMaxSATCalls;

  

  // statistics
  uint32_t _alwaysSATSCs;
  uint32_t _alwaysUNSATSCs;
  uint64_t _alwaysSATWeight;
  uint64_t _alwaysUNSATWeight;
  double _trimSATTime;
  uint32_t _noTrimSAT;
  uint32_t _noTrimSATSolverCalls;
  uint32_t _GBMOPartitions;
  double _GBMOTime;
  uint32_t _variablesOfEncoding;
  uint32_t _clausesOfEncoding;
  uint32_t _noSolverCalls;

  Encodings *_encodings;

  // relaxation literals or unit clause values!
  std::vector<uint32_t> _blockings; 
  uint32_t _negRelaxLit;
  // weights of relaxation literals
  std::vector<long long int> _weights;

  uint32_t SignedTouint32_tLit(int literal);
  void CalcGCDAndDivideIfPossible();
  uint64_t GreatestCommonDivisor(uint64_t a, uint64_t b);
  void AnalyzeSCsAndConvertIfPossible();
  uint64_t DivideSCsIfPossible();
  void genCardinals(long long tmpUnSATWeight, long long &divisor,
                    std::vector<uint32_t> &lits,
                    std::vector<uint32_t> &linkingVar,
                    std::vector<long long> &linkingWeight,
                    std::vector<long long> &divisors,
                    std::vector<std::vector<uint32_t>> &linkingVars,
                    std::vector<std::vector<long long>> &linkingWeights,
                    int compression);
  void HeuristicQMaxSAT(long long sum, long long k);
  void wbSortAndFilter(uint64_t UnSATWeight);
  void DumpCascCandidates();
  void AddSoftClauseTo(std::vector<SoftClause *> *softClauseVector,
                       std::vector<uint32_t> &clause, uint64_t weight);
  void GreedyMaximizeInitialSATWeight(uint32_t maxTime = 10,
                                      uint32_t maxSolves = 1000);
  std::vector<uint32_t> GetBestSCAssignment();
  void DivideSCs(std::vector<uint32_t> &sortedSCs, int acceptedMode = 0);
  void CalculateOverallTimes();
  void RemoveCascCand(uint32_t i);
  bool CheckMinWeightDist(std::vector<uint32_t> &sortedSCs, uint32_t firstPoint,
                          uint64_t biggerThan, uint32_t index);

  bool AddEncoding(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                   EncodingType *encodingType = nullptr);

  /**
   * @brief SolveMaxSAT straight MaxSAT solving of SoftClauseVector plus
   * encoding type without any preprocessing
   * @param tmpSoftClauses
   * @param encodingType
   */
  uint32_t SolveMaxSAT(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                       EncodingType *encodingType = nullptr);
  uint64_t CalculateLocalSATWeight(
      std::vector<SoftClause *> *tmpSatClauses = nullptr);
  void PrintResult();
  void DumpSolvingInformation();
  bool TreatBorderCases();
  void ChooseEncoding();
};
} // Namespace Pacose

#endif // PACOSE_H
