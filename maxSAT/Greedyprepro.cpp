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

#include <assert.h>
#include <cmath>
#include <iomanip> 
#include <iostream>
#include <tuple>

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>
#include "Greedyprepro.h"
#include "Pacose.h"
#include "../solver-proxy/SATSolverProxy.h"
#include "Settings.h"
#include "Softclause.h"
//#include "dgpw.h"  // for SAT UNSATISFIABLE
#include <chrono>

#include "Greedyprepro.h"
#include "timemeasurement.h"
#include "timevariables.h"

#include <iostream>

namespace Pacose {

#ifndef UNKNOW
#define UNKNOW 0
#endif
#ifndef SAT
#define SAT 10
#endif
#ifndef UNSAT
#define UNSAT 20
#endif


GreedyPrepro::GreedyPrepro(std::vector<SoftClause *> &softClauses,
                           Settings *settings, SATSolverProxy *solverProxy,
                           Pacose *pacose, int fixSCs)
    : _solver(solverProxy),
      _pacose(pacose),
      _settings(settings),
      _softClauses(softClauses),
      _satWeight(0),
      _unsatisfiableSCWeight(0),
      _satisfiableSCWeight(0),
      _sumOfSoftWeights(0),
      _unknownSolverCalls(0),
      _fixSoftClauses(fixSCs),
      _addedVariables(0),
      _addedClauses(0),
      _solverCalls(0),
      _unsatisfiableSCs(0),
      _satisfiableSCs(0),
      _previousPosition(1 - (_settings->greedyPPSATPercentage * 0.01)),
      _previousLowerBound(0.0),
      _noClauses(-1),
      _allWeightsSat(false),
      _pos1Unsat(false) {
  _timeVariables = new DGPW::TimeVariables();
  //  _solver = SATSolverProxy::InitSATSolver(SATSolverType::GLUCOSE421);
}

GreedyPrepro::~GreedyPrepro() {
  // SAT solver will be destroyed by Pacose
  // Pacose will be destroyed by application or main
  // Settings will be destroyed by Pacose
  // _softClauses will be destroyed by Pacose
  //  delete _timeVariables;
  //  delete _timeSolvedFirst;
}

uint32_t GreedyPrepro::NewVariable(void) {
  _addedVariables++;
  return static_cast<uint32_t>(_solver->NewVariable());
}

// TODO: optimize -> change interface to Glucose data structures
bool GreedyPrepro::AddClause(std::vector<uint32_t> &clause, uint32_t lbd) {
  _addedClauses++;
  return _solver->AddClause(clause);
}

uint32_t GreedyPrepro::Solve(void) {
  _solver->ClearAssumption();

  _solverCalls++;
  return _solver->Solve();
}

uint32_t GreedyPrepro::Solve(std::vector<uint32_t> &assumptions) {
  _solver->ClearAssumption();
  _solver->AddAssumptions(assumptions);
  _solverCalls++;
  return _solver->Solve();
}

uint32_t GreedyPrepro::SolveLimited(std::vector<uint32_t> &assumptions) {
  _solver->ClearAssumption();
  _solver->AddAssumptions(assumptions);
  _solverCalls++;
  return _solver->SolveLimited();
}

uint32_t GreedyPrepro::Model(uint32_t var) const {
  uint32_t value = _solver->GetModel(static_cast<int>(var));

  return value;
}

uint64_t GreedyPrepro::StartPrepro() {
  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

  //      static_cast<uint32_t>(_settings->greedyPPTimeLimit);
  uint32_t maxRounds =
      static_cast<uint32_t>(ceil(pow(_softClauses.size(), 0.3)) * 100);
  _timeLimit =
      (_settings->greedyPPTimeoutFactor * pow(_softClauses.size(), 0.35)) / 10;
  _preproPropagationLimit = static_cast<uint64_t>(
      _timeLimit * _settings->greedyPPPropagationPerSecond);
  if (_settings->verbosity > 0) {
    std::cout << "time Limit: " << _timeLimit << std::endl;
    std::cout << "preproPropagationLimit: " << _preproPropagationLimit
              << std::endl;
  }

  GreedyMaxInitSATWeightV2(1, maxRounds);

  return _opti;
}

void GreedyPrepro::DumpPreProInformation() {
  if (_settings->verbosity > 0) {
    if (_satisfiableSCs > 0) {
      std::cout << "c in any case SAT SCs....: " << _satisfiableSCs
                << std::endl;
    }

    if (_unsatisfiableSCs > 0) {
      std::cout << "c unsatisfiable SCs......: " << _unsatisfiableSCs
                << std::endl;
      std::cout << "c remaining SCs..........: " << _softClauses.size()
                << std::endl;
      std::cout << "c weight of unsat SCs....: " << _unsatisfiableSCWeight
                << std::endl;
    }

    //    if (onlyOneSCSatByChance > 0)
    //      std::cout << "c only 1 SC SAT by chance: " << onlyOneSCSatByChance
    //                << std::endl;
    //  std::cout << "c max SatWeight Prepro...: " << actualSATWeight <<
    //  std::endl;

    std::cout << "c max SatWeight Prepro...: " << _satWeight << std::endl;

    std::cout << "c sumOfSoftWeights.......: " << _sumOfSoftWeights
              << std::endl;

    std::cout << "c localPrePro o..........: " << _opti << std::endl;

    //            << std::endl;
    std::cout << "c time greedyPrepro......: "
              << _timeSolvedFirst->CurrentTimeDiff() << std::endl;
  }
}

void GreedyPrepro::RemoveAlwaysSatisfiedSoftClauses(
    std::vector<uint32_t> &sortedSCIndices) {
  
  // Same as WBSortAndFilter
  while (!_softClauses.empty() &&
         _opti < _softClauses[sortedSCIndices.back()]->weight) {
    _satisfiableSCs++;
    _satisfiableSCWeight += _softClauses[sortedSCIndices.back()]->originalWeight;

    std::vector<uint32_t> unitclause;
    unitclause.push_back(_softClauses[sortedSCIndices.back()]->relaxationLit ^
                         1);

    AddClause(unitclause);
    
    for (uint32_t i = 0; i < sortedSCIndices.size(); i++) {
      if (sortedSCIndices[i] > sortedSCIndices.back()) sortedSCIndices[i]--;
    }

    if (_settings->verbosity > 0) {
      std::cout << "c Softclause with weight "
                << _softClauses[sortedSCIndices.back()]->weight
                << " is always satisfied! - Remove this softclause!"
                << std::endl;
    }
    if (_settings->createSimplifiedWCNF) {
      // add this SC to the set of HardClauses!
      _newHardClauses.push_back(_softClauses[sortedSCIndices.back()]->clause);
    }
    assert(_satWeight >= _softClauses[sortedSCIndices.back()]->weight);
    // _satWeight -= _softClauses[sortedSCIndices.back()]->weight;
    assert(_pacose->_localSatWeight >= _softClauses[sortedSCIndices.back()]->weight);
    _pacose->_localSatWeight -= _softClauses[sortedSCIndices.back()]->weight;
    _pacose->_sumOfActualSoftWeights -= _softClauses[sortedSCIndices.back()]->weight;
    _softClauses.erase(_softClauses.begin() + sortedSCIndices.back());
    sortedSCIndices.pop_back();
  }
}

uint32_t GreedyPrepro::GreedyMaxInitSATWeightV2(int greedyPrepro,
                                                uint32_t maxRounds) {
  uint64_t softWeights = 0;
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { softWeights += SC->weight; });
  _sumOfSoftWeights = softWeights;

  if (_settings->verbosity > 0)
    std::cout << "c _sumOfSoftWeights : " << _sumOfSoftWeights << std::endl;

  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_settings->verbosity > 0) {
    std::cout << "c maxTime: " << _timeLimit << std::endl;
    std::cout << "c maxRounds: " << maxRounds << std::endl;
  }

  //  DGPW::TimeMeasurement timeSolvedFirst(&_timeVariables->solvedFirst, true);
  _timeSolvedFirst =
      new DGPW::TimeMeasurement(&_timeVariables->solvedFirst, true);
  uint32_t currentresult = 0;
  std::vector<uint32_t> nextAssumptions;
  std::vector<uint32_t> highestSATAssignment;
  std::vector<SoftClause *> UNSATSCs = {};
  std::vector<SoftClause *> neverSATSCs = {};
  uint64_t actualSATWeight = 0;
  std::vector<uint32_t> sortedSCIndices(_softClauses.size());
  uint32_t localMax = 0;
  uint32_t previousMax = 0;

  //  uint64_t unsatisfiableSCWeight = 0;

  uint32_t lastNewMaxRound = 0;

  _opti = static_cast<uint64_t>(-1);

  // to get the right order of weights to test them!
  std::size_t n(0);
  std::generate(
      std::begin(sortedSCIndices),
      std::begin(sortedSCIndices) + static_cast<uint32_t>(_softClauses.size()),
      [&] { return n++; });
  std::stable_sort(
      std::begin(sortedSCIndices), std::end(sortedSCIndices),
      [&](std::size_t i1, std::size_t i2) {
        return (_softClauses[i2]->weight > _softClauses[i1]->weight);
      });

  if (_softClauses[sortedSCIndices[0]]->weight ==
      _softClauses[sortedSCIndices.back()]->weight) {
    if (_settings->verbosity > 0)
      std::cout << "c it is unweighted MaxSAT!!!" << std::endl;
  }

  // first solve, with external assumptions, if given
  // could be the last satisfying model of another DGPW
  // for DivideDGPW mode!
  currentresult = Solve();
  assert(currentresult == SAT);
   
  if (greedyPrepro == 0) {
    return currentresult;
  }

  if (_settings->verbosity > 0) {
    std::cout << "c PrePro _timeLimit: " << _timeLimit << std::endl;
    std::cout << "c PrePro propagation limit: " << _preproPropagationLimit
              << std::endl;
  }

  if (_settings->verbosity > 1)
    std::cout << "c start solving model!" << std::endl;
  for (uint32_t rounds = 0; rounds < maxRounds; ++rounds) {
    // informations...
    if (_settings->verbosity > 0) {
      std::cout << "c time of greedy prepro: "
                << _timeSolvedFirst->CurrentTimeDiff() << std::endl;
      if (_settings->verbosity > 4)
        std::cout << "c local Max: " << localMax
                  << "  rounds / 2: " << rounds / 2
                  << "  rounds - lastNewMaxRound: " << rounds - lastNewMaxRound
                  << "  ceil(maxRounds/40): " << ceil(maxRounds / 40)
                  << std::endl;
      if (_settings->verbosity > 1)
        std::cout << "c rounds: " << rounds << std::endl;
    }

    //        std::vector<uint32_t> indices;
    //        for (uint32_t i = 0; i < UNSATSCs.size(); i++) {
    //            indices.push_back(i);
    //        }

    // get Info about the actual SC Model, the return values are sorted
    // highest weight first!
    tie(nextAssumptions, UNSATSCs, actualSATWeight) =
        SatisfiedSCsInfo(&sortedSCIndices);

    if (rounds == 0) {
      //      std::cout << "_SATWeight: " << _satWeight << std::endl;
      
      RemoveAlwaysSatisfiedSoftClauses(sortedSCIndices);

    }

    if (actualSATWeight > _satWeight) {
      _satWeight = actualSATWeight;
      //      std::cout << "ASW: " << actualSATWeight << std::endl;
      SaveHighestAssignment();

      _pacose->CalculateSATWeight();
      _satWeight = _pacose->_localSatWeight;
      _opti = _pacose->_localUnSatWeight;
      _sumOfSoftWeights = _opti + _satWeight;
      //      highestSATAssignment = GetLastSatisfiableAssignment();
      if (_settings->verbosity > 0) {
        std::cout << "c                                        NEW MAX FOUND: "
                  << _satWeight << std::endl;
        std::cout << "c local o " << _opti << std::endl;
      }
      lastNewMaxRound = rounds;
      if (_satWeight == _sumOfSoftWeights) {
        if (_settings->verbosity > 0)
          std::cout << "c All weights could be satisfied!!!" << std::endl;

        break;
      }
    }

    if (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit) {
      if (_settings->verbosity > 0)
        std::cout << "c timeout in greedy preprocessor!" << std::endl;
      break;
    } else if (_fixSoftClauses != 0 &&
               (_timeSolvedFirst->CurrentTimeDiff() >
                (_timeLimit * _settings->greedyPPSSCSwitchFactor * 0.01))) {
      // after 2/3 of the time, change strategy to only satisfy all SC at least
      // once!
      if (_settings->verbosity > 0) {
        std::cout << "_timeSolvedFirst->CurrentTimeDiff() > (_timeLimit * "
                     "_settings->greedyPPSSCSwitchFactor * 0.01))"
                  << _timeSolvedFirst->CurrentTimeDiff() << " < " << _timeLimit
                  << " * " << _settings->greedyPPSSCSwitchFactor << " * 0.01 "
                  << std::endl;
      }
      _fixSoftClauses = 0;
    }

    //    std::cout << std::setw(100) << "previousMax: " << previousMax <<
    //    std::endl; std::cout << std::setw(100) << "localMax: " << localMax <<
    //    std::endl;
    if (previousMax != localMax) {
      _previousPosition = 1 - (_settings->greedyPPSATPercentage * 0.01);
      _previousLowerBound = 0.0;
    }
    if (_fixSoftClauses == 0 || (_fixSoftClauses == 1 && localMax > 0)) {
      nextAssumptions.clear();
      neverSATSCs = BuildOrderedIntersection(&UNSATSCs, &neverSATSCs);
      currentresult = BinarySearchSatisfySCs(nextAssumptions, &neverSATSCs);
    } else {
      currentresult = BinarySearchSatisfySCs(nextAssumptions, &UNSATSCs);
    }
    previousMax = localMax;

    if (currentresult != SAT) {
      localMax++;
      if (_settings->verbosity > 1 && currentresult == UNSAT)
        std::cout << "c                                              The "
                  << localMax << " th local Max is found: " << actualSATWeight
                  << std::endl;

      neverSATSCs = BuildOrderedIntersection(&UNSATSCs, &neverSATSCs);
      if (neverSATSCs.empty()) {
        if (_settings->verbosity > 1)
          std::cout << "c All clauses were at least once satisfied!"
                    << std::endl;
        break;
      }

      if (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit) {
        if (_settings->verbosity > 0)
          std::cout << "c timeout in greedy preprocessor!" << std::endl;
        break;
      }
      if (_settings->verbosity > 1)
        std::cout << "c neverSATSCs.size(): " << neverSATSCs.size()
                  << std::endl;

      //      if (_noClauses != 1 || localMax == 1) {
      if ((currentresult == UNKNOW && _fixSoftClauses == 0) ||
          (_fixSoftClauses == 0 && _noClauses != 1) ||
          (_fixSoftClauses == 1 && (localMax == 1 || _noClauses != 1))) {
        nextAssumptions.clear();
        currentresult = Solve(nextAssumptions);
        continue;
      }

      bool rewrite_objective = false;
      for (uint32_t iter = 0; iter < neverSATSCs.size(); ++iter) {
        if (_fixSoftClauses == 2) {
          // no additional SCs can be satisfied!
          // all remaining SCs are unsat!

          nextAssumptions.clear();
          nextAssumptions.push_back(neverSATSCs[iter]->relaxationLit ^ 1);
          //                  _solver->DeactivateLimits();
          currentresult = Solve(nextAssumptions);
          if (currentresult == SAT) {
            if (_settings->verbosity > 0) {
              std::cout << "c                   SAT" << std::endl;
              std::cout << "c weight: " << neverSATSCs[iter]->weight
                        << "  relaxLit: " << neverSATSCs[iter]->relaxationLit
                        << std::endl;
            }
            break;
          }
        }
        // remove Softclause from set
        if (currentresult == UNSAT) {
          //                      std::cout << "c                   SAT" <<
          //                      std::endl;
          _unsatisfiableSCs++;
          _unsatisfiableSCWeight += neverSATSCs[iter]->originalWeight;

          if (_settings->verbosity > 0) {
            std::cout << "c Softclause with weight "
                      << neverSATSCs[iter]->weight
                      << " couldn't be satisfied! - Remove this softclause! - "
                         "new SC.size: "
                      << _softClauses.size() - 1 << "  neverSATSCs.size() "
                      << neverSATSCs.size() << "  UNSATSCs.size() "
                      << UNSATSCs.size()
                      << "   sortedSize: " << sortedSCIndices.size()
                      << std::endl;
          }
                    
          _opti -= neverSATSCs[iter]->weight;
          _pacose->_localUnSatWeight -= neverSATSCs[iter]->weight;
          _pacose->_sumOfActualSoftWeights -= neverSATSCs[iter]->weight;
          // TODO-Test Dieter: Rewrite objective to remove objective literal from the objective.
          // There was a solver-call where the negation of this call was an assumption, hence it was found as a core and is therefore implied by RUP.
          
          std::vector<uint32_t> unitclause; 
          unitclause.push_back(neverSATSCs[iter]->relaxationLit); 
          
          AddClause(unitclause);

          for (uint32_t j = 0; j < _softClauses.size(); ++j) {
            if (_softClauses[j]->relaxationLit ==
                neverSATSCs[iter]->relaxationLit) {
              _sumOfSoftWeights -= _softClauses[j]->weight;
              if (_settings->createSimplifiedWCNF) {
                // add this SC to the set of HardClauses!
                for (auto lit : _softClauses[j]->clause) {
                  std::vector<uint32_t> tmpUnitClause = {lit ^ 1};
                  _newHardClauses.push_back(tmpUnitClause);
                }
              }
              _softClauses.erase(_softClauses.begin() + j);
              break;
            }
          }

          neverSATSCs.erase(neverSATSCs.begin() + iter);
          _pos1Unsat = false;

          // to get the right order of weights to test them!
          sortedSCIndices.pop_back();
          std::size_t n(0);
          std::generate(std::begin(sortedSCIndices),
                        std::begin(sortedSCIndices) +
                            static_cast<uint32_t>(_softClauses.size()),
                        [&] { return n++; });
          std::stable_sort(
              std::begin(sortedSCIndices), std::end(sortedSCIndices),
              [&](std::size_t i1, std::size_t i2) {
                return (_softClauses[i2]->weight > _softClauses[i1]->weight);
              });

          iter--;
        }
        if (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit) {
          if (_settings->verbosity > 0)
            std::cout << "c timeout in greedy preprocessor!" << std::endl;
          break;
        }
      }

      if (currentresult != SAT) {
        if (_settings->verbosity > 0)
          std::cout << "c no more softclauses can be satisfieda!" << std::endl;
        break;
      }
    }

    //    if (neverSATSCs.empty() && localMax > 0) {

    if (neverSATSCs.empty() &&
        (_fixSoftClauses == 0 || (_fixSoftClauses != 1 && localMax > 0))) {
      if (_settings->verbosity > 0)
        std::cout << "c no more softclauses are satisfiable!" << std::endl;

      break;
    }
  }

  RemoveAlwaysSatisfiedSoftClauses(sortedSCIndices);

  DumpPreProInformation();

  return SAT;
}

void GreedyPrepro::SaveHighestAssignment() {
  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

  for (auto SC : _softClauses) {
    if (Model(SC->relaxationLit >> 1) == SC->relaxationLit) {
      SC->lastassignment = SC->relaxationLit;
      //            std::cout << "RelaxLit is NOT satisfiable" << std::endl;
      for (auto lit : SC->clause) {
        if (Model(lit >> 1) == lit) {
          // BUT clause is satisfied without relaxLit --> relaxlit can be
          // satisfied!
          SC->lastassignment = SC->relaxationLit ^ 1;
          break;
        }
      }
    } else {
      //            std::cout << "RelaxLit IS satisfiable" << std::endl;
      SC->lastassignment = SC->relaxationLit ^ 1;
    }
  }
}

std::tuple<std::vector<uint32_t>, std::vector<SoftClause *>, uint64_t>
GreedyPrepro::SatisfiedSCsInfo(std::vector<uint32_t> *sortedSCIndices) {
  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<uint32_t> nextAssumptions;
  std::vector<SoftClause *> UNSATSCs = {};
  uint64_t actualWeight = 0;
  uint64_t unsatWeight = 0;

  for (auto i : *sortedSCIndices) {
    // highest not yet looked at weight!
    SoftClause *SC = _softClauses[i];
    //          std::cout << "SCWEIGHT: " << SC->weight << std::endl;
    //          std::cout << "NEXT RELAXLIT: " << SC->relaxationLit <<
    //          std::endl; std::cout << "Number of vars: " << Variables() <<
    //          std::endl; std::cout << "NEXT MODEL OF RELAXLIT: " <<
    //          Model(SC->relaxationLit>>1) << std::endl;

    if (Model(SC->relaxationLit >> 1) == SC->relaxationLit) {
      //            std::cout << "RelaxLit is NOT satisfiable" << std::endl;
      for (auto lit : SC->clause) {
        if (Model(lit >> 1) == lit) {
          // clause is satisfied without relaxLit
          nextAssumptions.push_back(SC->relaxationLit ^ 1);
          actualWeight += SC->weight;
          break;
        }
      }
    } else {
      //            std::cout << "RelaxLit IS satisfiable" << std::endl;
      nextAssumptions.push_back(SC->relaxationLit ^ 1);
      actualWeight += SC->weight;
      //            std::cout << "RelaxLit IS satisfiable" << std::endl;
    }
    if (nextAssumptions.empty() ||
        nextAssumptions.back() != (SC->relaxationLit ^ 1)) {
      UNSATSCs.push_back(SC);
      unsatWeight += SC->weight;
    }
  }
  //    std::cout << "DONE!" << std::endl;
  if (_settings->verbosity > 0)
    std::cout << "c local o: " << unsatWeight + _unsatisfiableSCWeight
              << std::endl;
  //    _sumOfSoftWeights << std::endl;

  _sumOfSoftWeights = unsatWeight + actualWeight;

  if (_opti < unsatWeight || _opti == static_cast<uint64_t>(-1)) {
    _opti = unsatWeight;
    _satWeight = actualWeight;
    _pacose->CalculateSATWeight();
  }
  assert(_sumOfSoftWeights = _satWeight + _opti);
  //    std::cout << "DONE!" << std::endl;
  return std::make_tuple(nextAssumptions, UNSATSCs, actualWeight);
}

std::vector<SoftClause *> GreedyPrepro::BuildOrderedIntersection(
    std::vector<SoftClause *> *SCsetA, std::vector<SoftClause *> *SCsetB) {
  if (_settings->verbosity > 2) std::cout << __PRETTY_FUNCTION__ << std::endl;

  if ((*SCsetB).size() == 0) {
    return *SCsetA;
  }

  std::vector<SoftClause *> v_intersection;

  // sort both sets of soft clauses
  std::sort((*SCsetB).begin(), (*SCsetB).end(), SoftClause::ids);
  std::sort((*SCsetA).begin(), (*SCsetA).end(), SoftClause::ids);

  //              std::cout << "Build intersection" << std::endl;
  std::set_intersection((*SCsetB).begin(), (*SCsetB).end(), (*SCsetA).begin(),
                        (*SCsetA).end(), std::back_inserter(v_intersection),
                        SoftClause::ids);

  std::sort(v_intersection.begin(), v_intersection.end(), SoftClause::bigger);

  return v_intersection;
}

uint32_t GreedyPrepro::BinarySearchSatisfySCs(
    std::vector<uint32_t> &nextAssumptions,
    std::vector<SoftClause *> *unsatSCs) {
  if (_settings->verbosity > 2) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    std::cout << "c greedy strategy: kind of binary search to force optimal "
                 "number of SC to be true!"
              << std::endl;
  }
  double noVars(-1);

  if (unsatSCs->empty()) {
    if (_settings->verbosity > 1) {
      std::cout << "No more SCs to satisfy!" << std::endl;
    }
    return SAT;
  }

  //  unsatSCs = &_softClauses;

  if (_settings->verbosity > 0) {
    std::cout << std::setw(30) << "_previousPosition: " << std::setw(10)
              << _previousPosition << std::endl;
    std::cout << std::setw(30) << "_previousLowerBound: " << std::setw(10)
              << _previousLowerBound << std::endl;
  }
  assert(_previousPosition <= 1);
  assert(_previousPosition > 0);

  uint32_t noUnsatSCs = unsatSCs->size();
  //  noUnsatSCs = 25;
  //  _previousPosition = 0.18;

  double sqrt = std::sqrt(noUnsatSCs);
  uint32_t numberOfPositions = (floor(sqrt * 2)) - 1;
  if (noUnsatSCs == 2) {
    numberOfPositions++;
  }

  uint32_t newPosition = 0;

  if (_settings->verbosity > 5) {
    std::cout << std::setw(30) << "newPosition: " << std::setw(10)
              << newPosition << std::endl;

    std::cout << std::setw(30) << "noUnsatSCs: " << std::setw(10) << noUnsatSCs
              << std::endl;
    std::cout << std::setw(30) << "numberOfPositions: " << std::setw(10)
              << numberOfPositions << std::endl;
    std::cout << std::setw(30) << "sqrt: " << std::setw(10) << sqrt
              << std::endl;
    std::cout << "If we have n clauses, at least n SC have to be satisfied!"
              << std::endl;
    std::cout
        << "If we have n vars, at least one out of n SC have to be satisfied!"
        << std::endl;
    for (uint32_t i = 1; i < numberOfPositions + 1; i++) {
      newPosition = i;
      if (newPosition < sqrt) {
        //        noClauses = ceil(noUnsatSCs / newPosition);
        _noClauses = noUnsatSCs / newPosition;
      } else {
        _noClauses = numberOfPositions - newPosition + 1;
      }
      noVars = noUnsatSCs / _noClauses;
      std::cout << "Position, noClauses, noVars: " << std::setw(5) << i
                << std::setw(8) << _noClauses << std::setw(8)
                << std::round(noVars * 100) / 100 << std::endl;
    }
  }

  newPosition = ceil(_previousPosition * numberOfPositions);
  if (_pos1Unsat && newPosition == 1 && numberOfPositions > 1) newPosition++;

  if (newPosition < sqrt) {
    //    noClauses = ceil(noUnsatSCs / newPosition);
    assert(newPosition != 0);
    //    if (newPosition == 0) newPosition = 1;
    _noClauses = noUnsatSCs / newPosition;
  } else {
    _noClauses = numberOfPositions - newPosition + 1;
  }
  noVars = (noUnsatSCs / _noClauses);
  if (_settings->verbosity > 2) {
    std::cout << std::endl
              << "Position, noClauses, noVars: " << std::setw(5) << newPosition
              << std::setw(8) << _noClauses << std::setw(8)
              << std::round(noVars * 100) / 100 << std::endl;
  }

  uint32_t upperBound = noUnsatSCs % static_cast<int>(floor(noVars));
  if (_settings->verbosity > 3) {
    std::cout << "Number of clauses  X  Clause Size: " << std::setw(10)
              << upperBound << " X " << ceil(noVars) << std::endl;
    std::cout << "Number of clauses  X  Clause Size: " << std::setw(10)
              << _noClauses - upperBound << "  X  " << floor(noVars)
              << std::endl;
  }

  //  std::vector<uint32_t> &nextAssumptions;

  //  uint32_t round;
  //  exit(1);


  // initialing clause vector with relaxlit!
  uint32_t relaxlit = NewVariable() << 1;
  std::vector<uint32_t> clause{relaxlit};
  std::vector<std::vector<uint32_t>> clauses(_noClauses, clause);

  // generate randomly sorted vector of indices
  std::vector<uint32_t> indices(unsatSCs->size());
  std::size_t n(0);
  std::generate(std::begin(indices),
                std::begin(indices) + static_cast<uint32_t>(unsatSCs->size()),
                [&] { return n++; });

  // uint32_t seed =
  // std::chrono::system_clock::now().time_since_epoch().count();
  uint32_t seed = 1234567890;
  std::shuffle(indices.begin(), indices.end(),
               std::default_random_engine(seed));

  int varsPerClause = ceil(noVars);
  //  std::cout << "unsatSCs->size(): " << unsatSCs->size()
  //            << " clauses.size(): " << clauses.size() << std::endl;

  size_t i = 0;
  for (uint32_t j = 0; j < clauses.size(); j++) {
    if (j == upperBound) varsPerClause = floor(noVars);
    for (int k = 0; k < varsPerClause; k++) {
      //      std::cout << "i: " << i << " varsPerClause: " << varsPerClause
      //                << " j: " << j << " UB: " << upperBound << "; ";
      if (_settings->verbosity > 4)
        std::cout << ((*unsatSCs)[indices[i]]->relaxationLit ^ 1) << " ";
      //                << std::endl;
      clauses[j].push_back((*unsatSCs)[indices[i]]->relaxationLit ^ 1);
      i++;
    } 
    //  TODO-Test Dieter: 
    // Here, clauses of the form 
    //             a + b+ c + ... + relaxLit >= 1
    //             x + y + z + ... + relaxLit >= 1
    //             ...
    // are created for a new literal relaxLit. 
    // These clauses can be proven by RBS with witness relaxLit -> 1. 
    // Next, the SAT solver is called with assumption ~relaxLit.  Next time a new set of such clauses is created. 
    // In the special case where we only have one clause of such form, i.e., only the clause a + b+ c + ... + relaxLit >= 1 is created,
    // we receive a core relaxLit >= 1 from the solver. Since there is no assignment that assigns relaxLit false, we know that we can derive the clauses 
    // ~a >= 1, ~b >= 1, ... To derive these clauses, we first need to derive the other direction, namely 20 ~relaxationLit 1 ~x1 + 1 ~x2 + ... + 1 ~x20 >= 20 by RBS (witness relaxationLit -> 0).
    // This constraint propagates to 1 ~x1 + 1 ~x2 + ... + 1 x20 >= 20 after propagating core relaxationLit >= 1, to which clause ~a >= 1 is RUP. 
    // This constraint for the other direction only needs to be derived whenever _noClauses == 1.
    // This is only the case if nextAssumptions is empty!!!
    //
    // If this is not the case, the clauses  a + b+ c + ... + relaxLit >= 1 are removed again by introducing the unit clause relaxLit >= 1. Can be proven by RBS using witness relaxLit -> 1. Only proof obligations on clauses containing r, which are all satisfied by witness.
      
    AddClause(clauses[j]);
    
    if (_settings->verbosity > 4) std::cout << std::endl;
  }

  assert(i == unsatSCs->size());

  //  //      _solver->SetMaxPropagations(_preproPropagationLimit);

  nextAssumptions.push_back(relaxlit ^ 1);

  double solvingTime = _timeSolvedFirst->CurrentTimeDiff();
  _solver->SetPropagationBudget(_preproPropagationLimit);

  // Solver call with time / propagation Limit
  uint32_t currentresult = SolveLimited(nextAssumptions);

  //  _solver->SetPropagationBudget(-1);
  solvingTime = _timeSolvedFirst->CurrentTimeDiff() - solvingTime;
  if (_settings->verbosity > 0)
    std::cout << "c Solver Call Time: " << solvingTime << std::endl;

  nextAssumptions.pop_back();
  std::vector<uint32_t> UC = {relaxlit};
  AddClause(UC); // deactivate added clauses again
  _solver->ClearAssumption();

  //  double noClauses;  // number of clauses added, each clause forcing at
  //  least
  // one SC to be true
  //  double noVars;  // number of vars in one clause, one SC out of noVars has
  //  to
  // be satisfied.

  //  greedyPPSATPercentage;
  //  greedyPPUNSATPercentage;
  //  greedyMinSizeOfSet;

  bool firstTimeLimitReached =
      (_timeSolvedFirst->CurrentTimeDiff() >
       (_timeLimit * _settings->greedyPPSSCSwitchFactor * 0.01));
  //  bool firstTimeLimitReached =
  //      (_timeSolvedFirst->CurrentTimeDiff() > (_timeLimit / 3) * 2);
  bool timeLimitReached = (_timeSolvedFirst->CurrentTimeDiff() > _timeLimit);

  if (currentresult == SAT) {
    
    if (_settings->verbosity > 0) {
      std::cout << std::setw(100) << "SAT" << std::endl;
      if (noVars == 1) {
        std::cout << std::setw(100) << "c ALL REMAINING SCs COULD BE SATISFIED AT LEAST ONCE!"
                  << std::endl;
      }
    }
    _previousPosition -= (_previousPosition - _previousLowerBound) /
                         (1 / (_settings->greedyPPSATPercentage * 0.01));
    if (_previousPosition == 0) _previousPosition = 0.000001;

    return currentresult;
  } else if (currentresult == UNSAT) {
    if (_settings->verbosity > 0)
      std::cout << std::setw(100) << "UNSAT" << std::endl;
    _previousLowerBound = _previousPosition;
    _previousPosition += (1 - _previousLowerBound) /
                         (1 / (_settings->greedyPPUNSATPercentage * 0.01));
    if (_previousPosition > 1) _previousPosition = 1;

    // no more soft clauses can be satisfied
    if (_noClauses == 1) {
      if (_settings->verbosity > 0)
        std::cout
            << std::setw(100)
            << "c NO ADDITIONAL SOFTCLAUSES FROM THIS SET ARE SAT!"
            << std::endl;

      return currentresult;
    }

    if (newPosition == 1) {
      // we do not try again to satisfy all remaining softclauses at once
      _pos1Unsat = true;
    }

    if (timeLimitReached || (firstTimeLimitReached && _fixSoftClauses == 2)) {
      return UNKNOW;
    }
    return BinarySearchSatisfySCs(nextAssumptions, unsatSCs);
  } else {
    _unknownSolverCalls++;
    if (_settings->verbosity > 0)
      std::cout << std::setw(100)
                << "PROPAGATION LIMIT REACHED FOR SINGLE SOLVER CALL"
                << std::endl;
    _previousLowerBound = _previousPosition;
    _previousPosition += (1 - _previousLowerBound) /
                         (1 / (_settings->greedyPPUNSATPercentage * 0.01));
    _previousLowerBound = _previousPosition;
    _previousPosition += (1 - _previousLowerBound) /
                         (1 / (_settings->greedyPPUNSATPercentage * 0.01));
    if (_previousPosition > 1) _previousPosition = 1;

    // no more soft clauses can be satisfied
    if (_noClauses == 1 && _fixSoftClauses == 0) {
      // the last remaining SC could not be satisfied, reduce timelimit to kill
      // prepro!
      _timeLimit = 0;
      return currentresult;
    } else if (_noClauses == 1) {
      _fixSoftClauses = 0;
      return currentresult;
    }

    if (_unknownSolverCalls < 2 &&
        ((!firstTimeLimitReached && _fixSoftClauses == 2) ||
         !timeLimitReached)) {
      return BinarySearchSatisfySCs(nextAssumptions, unsatSCs);
    } else if (_unknownSolverCalls == 2 || firstTimeLimitReached) {
      if (_fixSoftClauses != 0 || firstTimeLimitReached) {
        _fixSoftClauses = 0;
        return currentresult;
      } else {
        return BinarySearchSatisfySCs(nextAssumptions, unsatSCs);
      }
    } else if (_unknownSolverCalls > 2) {
      _timeLimit = 0;
      return currentresult;
    }
  }
  return 0;
}

void GreedyPrepro::CreateSimplifiedWCNF(
    std::string filename, std::vector<std::vector<uint32_t>> &hardClauses) {
  if (_newHardClauses.empty()) {
    std::cout << "c There are no soft clauses removed - exit!" << std::endl;
    exit(1);
  }

  size_t fileSep = filename.find_last_of("/");
  filename.erase(filename.end() - 5, filename.end());
  filename = filename.substr(fileSep + 1) + "_" +
             std::to_string(_unsatisfiableSCs) + "unsatSCs_" +
             std::to_string(_satisfiableSCs) + "satSCs_" + "removed_FS" +
             std::to_string(_settings->greedyPPFixSCs) + ".wcnf";

  //  std::cout << "New Filename: " << filename << std::endl;

  uint64_t topWeight = 1;

  for (auto sc : _softClauses) {
    topWeight += sc->originalWeight;
  }
  //  std::cout << "NoSoftClauses: " << _softClauses.size() << std::endl;
  //  std::cout << "TopWeight: " << topWeight << std::endl;

  uint32_t maxVar = 0;
  for (auto hc : hardClauses) {
    for (auto lit : hc) {
      if (lit > maxVar) maxVar = lit;
    }
  }

  for (auto sc : _softClauses) {
    for (auto lit : sc->clause) {
      if (lit > maxVar) maxVar = lit;
    }
  }

  for (auto hc : _newHardClauses) {
    for (auto lit : hc) {
      if (lit > maxVar) maxVar = lit;
    }
  }

  size_t SCsize = _softClauses.size();

  if (SCsize == 0)
    SCsize = 1;
  else
    SCsize = _softClauses.size();
  //  std::cout << "p wcnf " << maxVar << " "
  //            << hardClauses.size() + SCsize + _newHardClauses.size() << " "
  //            << topWeight << std::endl;
  std::ofstream out(filename);
  //  out.open(filename, std::ofstream::out | std::ofstream::app);
  std::streambuf *coutbuf = std::cout.rdbuf();  // save old buf
  std::cout.rdbuf(out.rdbuf());  // redirect std::cout to filename!
  std::cout << "c " << _unsatisfiableSCs << " soft clauses with weight "
            << _unsatisfiableSCWeight << " are in the unsatisfiable backbone!"
            << std::endl;
  std::cout << "c " << _satisfiableSCs << " soft clauses with weight "
            << _satisfiableSCWeight << " are in the satisfiable backbone!"
            << std::endl;

  if (_softClauses.size() == 0) {
    std::cout << "c " << _unsatisfiableSCs
              << " all SoftClauses could be removed. Add again first HC as "
                 "SC - such that all Solver can solve the instances!"
              << std::endl;
    topWeight = 2;
  }

  std::cout << "p wcnf " << maxVar << " "
            << hardClauses.size() + SCsize + _newHardClauses.size() << " "
            << topWeight << std::endl;

  for (auto hc : hardClauses) {
    std::cout << topWeight << " ";
    for (auto lit : hc) {
      if (lit % 2 == 0)
        std::cout << (lit >> 1) << " ";
      else
        std::cout << "-" << (lit >> 1) << " ";
    }
    std::cout << "0" << std::endl;
  }

  // Removed soft clauses have to be added as new hard clauses!
  for (auto newHC : _newHardClauses) {
    std::cout << topWeight << " ";
    for (auto lit : newHC) {
      if (lit % 2 == 0)
        std::cout << (lit >> 1) << " ";
      else
        std::cout << "-" << (lit >> 1) << " ";
    }
    std::cout << "0" << std::endl;
  }

  if (_softClauses.size() == 0) {
    // if all soft clauses are removed add first HC again as SoftClause
    // some MaxSAT solver do not work without this feature!
    std::cout << "1 ";
    for (auto lit : hardClauses[0]) {
      if (lit % 2 == 0)
        std::cout << (lit >> 1) << " ";
      else
        std::cout << "-" << (lit >> 1) << " ";
    }
    std::cout << "0" << std::endl;
  } else {
    for (auto sc : _softClauses) {
      std::cout << sc->weight << " ";
      for (auto lit : sc->clause) {
        if (lit % 2 == 0)
          std::cout << (lit >> 1) << " ";
        else
          std::cout << "-" << (lit >> 1) << " ";
      }
      std::cout << "0" << std::endl;
    }
  }

  std::cout.rdbuf(coutbuf);  // reset to standard output again

  std::cout << "c Clauses are written to CNF file: " << filename << std::endl;
  exit(1);
}
} // Namespace Pacose
