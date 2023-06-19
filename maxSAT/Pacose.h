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
#include "Settings.h"
// #include "../maxpre2/src/preprocessorinterface.hpp"



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

  unsigned SolveProcedure(ClauseDB& clauseDB);
  bool ExternalPreprocessing(ClauseDB& clauseDB);
  // void CallMaxPre2(ClauseDB &clauseDB);

  unsigned CalculateNextResult();
  unsigned CalculateNextSoftclauseCombination();
  
  unsigned SolveQMax(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                     EncodingType *encodingType = nullptr);
  
  void Preprocess();
  void ExcludeCurrentSoftclauseCombination();
  void ExcludeCurrentResult();
  

  // Variables
  Settings _settings;
  //    Encodings _encodings;

  SATSolverProxy *_satSolver;
  unsigned _nbVars;
  unsigned _nbClauses;
  unsigned _nbOriginalClauses;
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
  void AddSoftClause(std::vector<uint32_t> &clause, uint64_t weight = 1);

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
  void SetSumOfSoftWeights(unsigned long softWeights);

  // literal given as 2 x int (+1)
  void AddAssumption(unsigned literal);
  void ClearAssumptions();
  void DeactivateLastCNF();

  bool AddClause(std::vector<unsigned> &clause);

  unsigned NewVariable();
  void NewVariables(unsigned noVars);
  unsigned GetModel(int var);

  unsigned Solve();

  bool DumpWCNF();

private:
  struct partitionInformation
  {
    partitionInformation() : dgpw(nullptr) {}

    DGPW::DGPW *dgpw;
    unsigned Points;
    uint64_t weightsTillPoint;
    uint64_t ggtTillPoint;
    bool allWeightsAreEqual;
  };

  // maxPreprocessor::PreprocessorInterface* maxpre;

  std::vector<partitionInformation> _cascCandidates;
  EncodingType _encoding;

  // Variables
  std::vector<std::vector<unsigned>> _CNF;


  // std::vector<std::vector<std::vector<unsigned>>> _incrementalCNF;
  int _cpuLimit;
  int _memLimit;
  int _nbOfOrigVars;
  uint64_t  _sumOfSoftWeights;
  uint64_t _overallSoftWeights;
  // for incremental MaxSAT
  // bool _withAssumption;
  // std::vector<unsigned> _assumptions;
  // unsigned _relaxVar;
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
  // unsigned _incrementalMaxSATCalls;

  

  // statistics
  unsigned _alwaysSATSCs;
  unsigned _alwaysUNSATSCs;
  uint64_t _alwaysSATWeight;
  uint64_t _alwaysUNSATWeight;
  double _trimSATTime;
  unsigned _noTrimSAT;
  unsigned _noTrimSATSolverCalls;
  unsigned _GBMOPartitions;
  double _GBMOTime;
  unsigned _variablesOfEncoding;
  unsigned _clausesOfEncoding;
  unsigned _noSolverCalls;

  Encodings *_encodings;

  // relaxation literals or unit clause values!
  std::vector<unsigned> _blockings; 
  unsigned _negRelaxLit;
  // weights of relaxation literals
  std::vector<long long int> _weights;

  unsigned SignedToUnsignedLit(int literal);
  void CalcGCDAndDivideIfPossible();
  uint64_t GreatestCommonDivisor(uint64_t a, uint64_t b);
  void AnalyzeSCsAndConvertIfPossible();
  uint64_t DivideSCsIfPossible();
  void genCardinals(long long tmpUnSATWeight, long long &divisor,
                    std::vector<unsigned> &lits,
                    std::vector<unsigned> &linkingVar,
                    std::vector<long long> &linkingWeight,
                    std::vector<long long> &divisors,
                    std::vector<std::vector<unsigned>> &linkingVars,
                    std::vector<std::vector<long long>> &linkingWeights,
                    int compression);
  void HeuristicQMaxSAT(long long sum, long long k);
  void wbSortAndFilter(uint64_t UnSATWeight);
  void DumpCascCandidates();
  void AddSoftClauseTo(std::vector<SoftClause *> *softClauseVector,
                       std::vector<unsigned> &clause, uint64_t weight);
  void GreedyMaximizeInitialSATWeight(unsigned maxTime = 10,
                                      unsigned maxSolves = 1000);
  std::vector<unsigned> GetBestSCAssignment();
  void DivideSCs(std::vector<unsigned> &sortedSCs, int acceptedMode = 0);
  void CalculateOverallTimes();
  void RemoveCascCand(unsigned i);
  bool CheckMinWeightDist(std::vector<unsigned> &sortedSCs, unsigned firstPoint,
                          unsigned long biggerThan, unsigned index);

  bool AddEncoding(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
                   EncodingType *encodingType = nullptr);

  /**
   * @brief SolveMaxSAT straight MaxSAT solving of SoftClauseVector plus
   * encoding type without any preprocessing
   * @param tmpSoftClauses
   * @param encodingType
   */
  unsigned SolveMaxSAT(std::vector<SoftClause *> *tmpSoftClauses = nullptr,
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
