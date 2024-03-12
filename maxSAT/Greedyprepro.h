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

#ifndef GREEDYPREPRO_H
#define GREEDYPREPRO_H

#include <cstdint>
#include <vector>
#include "../solver-proxy/SATSolverProxy.h"
#include "Settings.h"
#include "Softclause.h"
#include "DGPW/timemeasurement.h"
#include "DGPW/timevariables.h"

namespace Pacose {

class Pacose;
// struct TimeVariables;

class GreedyPrepro {
 public:
  GreedyPrepro(std::vector<SoftClause *> &softClauses, Settings *settings,
               SATSolverProxy *satSolver, Pacose *pacose, int fixSCs);
  ~GreedyPrepro();

  uint64_t StartPrepro();

  uint64_t GetAlwaysUNSATWeight() { return _unsatisfiableSCWeight; };
  uint64_t GetAlwaysSATWeight() { return _satisfiableSCWeight; };
  uint64_t GetAlwaysUNSATSCs() { return _unsatisfiableSCs; };
  uint64_t GetAlwaysSATSCs() { return _satisfiableSCs; };

  void CreateSimplifiedWCNF(std::string filename,
                            std::vector<std::vector<uint32_t>> &hardClauses);

  void DumpPreProInformation();
  
  uint64_t _satWeight;

 private:
  uint32_t Solve(std::vector<uint32_t> &assumptions);
  uint32_t SolveLimited(std::vector<uint32_t> &assumptions);
  uint32_t Solve();
  bool AddClause(std::vector<uint32_t> &clause, uint32_t lbd = 1);
  uint32_t NewVariable();

  SATSolverProxy *_solver;
  Pacose *_pacose;

  Settings *_settings;
  DGPW::TimeVariables *_timeVariables;
  DGPW::TimeMeasurement *_timeSolvedFirst;
  std::vector<SoftClause *> &_softClauses;
  std::vector<std::vector<uint32_t>> _newHardClauses;
  uint64_t _unsatisfiableSCWeight;
  uint64_t _satisfiableSCWeight;

  uint64_t _sumOfSoftWeights;
  uint64_t _opti;
  int64_t _preproPropagationLimit;
  int _unknownSolverCalls;
  int _fixSoftClauses;
  double _timeLimit;
  uint32_t _addedVariables;
  uint32_t _addedClauses;
  uint32_t _solverCalls;
  uint32_t _unsatisfiableSCs;
  uint32_t _satisfiableSCs;
  double _previousPosition;
  double _previousLowerBound;
  double _noClauses;
  bool _allWeightsSat;
  bool _pos1Unsat;

  uint32_t GreedyMaxInitSATWeightV2(int greedyPrepro, uint32_t maxRounds);
  std::tuple<std::vector<uint32_t>, std::vector<SoftClause *>, uint64_t>
  SatisfiedSCsInfo(std::vector<uint32_t> *sortedSCIndices);
  uint32_t Model(uint32_t var) const;
  std::vector<SoftClause *> BuildOrderedIntersection(
      std::vector<SoftClause *> *UNSATSCs,
      std::vector<SoftClause *> *neverSATSCs);
  uint32_t BinarySearchSatisfySCs(std::vector<uint32_t> &nextAssumptions,
                                  std::vector<SoftClause *> *unsatSCs);
  void SaveHighestAssignment();
  void RemoveAlwaysSatisfiedSoftClauses(std::vector<uint32_t> &sortedSCIndices);
};
} // Namespace Pacose

#endif  // GREEDYPREPRO_H
