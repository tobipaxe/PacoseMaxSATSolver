/********************************************************************************************
Copyright (c) 2014-2016, Sven Reimer,
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
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

#include <assert.h>
#include <cmath> 
#include <iomanip>

#include "../Pacose.h"
#include "../../solver-proxy/SATSolverProxy.h"
#include "../Softclause.h"
#include "cascade.h"
#include "dgpw.h"
#include "multiplecascade.h"
#include "sorter.h"
#include "timemeasurement.h"
#include "timevariables.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

namespace Pacose {
namespace DGPW {
DGPW::DGPW(Pacose *pacose)
    : _sumOfSoftWeights(),
      _timeVariables(nullptr),
      _softClauses(),
      _sorterTree(),
      _solver(pacose->_satSolver),
      _dgpwSetting(&pacose->_settings),
      _pacose(pacose),
      _maxSorterDepth(0),
      _maxsatResult(UNKNOWN),
      _mainCascade(nullptr),
      _mainMultipleCascade(nullptr),
      _greatestCommonDivisor(1),
      _currentBucketForTare(0),
      _satWeight(0),
      _moreThanTwoWeights(false),
      //    _topWeight(pacose->_top),
      _minWeight(1),
      _maxWeight(0),
      _hasEncoding(false),
      _resultUnknown(false),
      _clausesBefore(pacose->_nbClauses),
      _variablesBefore(pacose->_nbVars),
      _binaryClausesBefore(0),
      _ternaryClausesBefore(0),
      _addedClauses(0),
      _addedVariables(0),
      _addedBinaryClauses(0),
      _addedTernaryClauses(0),
      _binaryTopClauses(0),
      _ternaryTopClauses(0),
      _binaryBottomClauses(0),
      _ternaryBottomClauses(0),
      _solverCalls(0),
      _approxSolverCalls(0),
      _highestBucketMultiplicator(0),
      _currentBucketMultiplicator(0),
      _estimatedSatWeight(0),
      _diffToSatWeight(0),
      //    _collectedAssumptions(),
      _externalAssumptions({}),
      _fixedAssumptions({}),
      _lastResult(0) {
  _sorterTree.resize(1);
}

void DGPW::IncrementalReset() {
  *_optimum = -1;
  _maxsatResult = UNKNOWN;
  _currentBucketForTare = 0;
  _satWeight = 0;
  _moreThanTwoWeights = true;
  _minWeight = 1;
  _maxWeight = 0;
  _resultUnknown = false;
  _clausesBefore = _pacose->_nbClauses;
  _variablesBefore = _pacose->_nbVars;
  _solverCalls = 0;
  _approxSolverCalls = 0;
  _currentBucketMultiplicator = 1;
  _estimatedSatWeight = 0, _diffToSatWeight = 0;
  _lastResult = 0;
  _mainCascade->IncrementalReset();
}

DGPW::~DGPW() {
  for (uint32_t i = 0; i != _sorterTree.size(); ++i) {
    for (uint32_t j = 0; j != _sorterTree[i].size(); ++j) {
      delete _sorterTree[i][j];
    }
  }
  delete _timeVariables;
  delete _mainCascade;
  delete _mainMultipleCascade;

  // the following are destroyed elsewhere!
  //  SATSolverProxy *_solver;
  //  Settings *_dgpwSetting;
  //  Pacose *_pacose;
}

/********** Interface to Glucose and reimplemented functions from AntomBase
 * **********/

// Creates and returns a fresh variable index not used so far.
uint32_t DGPW::NewVariable(void) {
  _addedVariables++;
  return static_cast<uint32_t>(_solver->NewVariable());
  //    Glucose::Var newVar = _glucose->newVar(true, true);
  //    return static_cast<uint32_t>(newVar) + 1;
}

// TODO: optimize -> change interface to Glucose data structures
bool DGPW::AddClause(std::vector<uint32_t> &clause, uint32_t lbd) {
  //  std::cout << "Clause: ";
  //  for (auto literal : clause) {
  //    std::cout << literal << " ";
  //  }
  //  std::cout << std::endl;
  _addedClauses++;
  return _solver->AddClause(clause);
  //    Glucose::vec<Glucose::Lit> newClause;
  //    for (uint32_t lit : clause )
  //      {
  //        assert(lit > 1);
  //        Glucose::Lit newLit = Glucose::toLit(static_cast<int>(lit - 2));
  //        newClause.push(newLit);
  //      }
  //    return _glucose->addClause_(newClause);
}

//  bool DGPW::AddClause(Glucose::vec<Glucose::Lit>& clause)
//  {
//    return _glucose->addClause_(clause);
//  }

uint32_t DGPW::GetLastResult() { return _lastResult; }

bool DGPW::AddUnit(uint32_t lit) {
  std::vector<uint32_t> unitclause;
  unitclause.push_back(lit);
  //  std::cout << " "
  //               "                          UNIT CLAUSE: "
  //            << lit << std::endl;

  return AddClause(unitclause);
}

uint32_t DGPW::Solve(void) {
  _solver->ClearAssumption();

  if (!_fixedAssumptions.empty()) {
    _solver->AddAssumptions(_fixedAssumptions);
  }
  _lastResult = _solver->Solve();
  _solverCalls++;
  return _lastResult;
}

uint32_t DGPW::Solve(std::vector<uint32_t> &assumptions) {
  //  std::cout << "SOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOLVE(assumptions)!"
  //            << std::endl;
  _solver->ClearAssumption();
  if (!_fixedAssumptions.empty()) {
    _solver->AddAssumptions(_fixedAssumptions);
  }
  _solver->AddAssumptions(assumptions);
  _lastResult = _solver->Solve();
  _solverCalls++;
  return _lastResult;
  //    std::cout << "result: " << rv << std::endl;
  //    return rv;
}

// TODO: optimize!
// copies whole model vector
std::vector<uint32_t> DGPW::Model(void) const {
  //    std::cout << __PRETTY_FUNCTION__ << "  x  " << std::endl;

  //    _solver->ClearAssumptions();
  //    std::cout << _solver->Solve() << std::endl;

  std::vector<uint32_t> model;

  for (int32_t i = 0; i < static_cast<int32_t>(_solver->GetNumberOfVariables());
       ++i) {
    model.push_back(_solver->GetModel(i));
  }
  return model;
}

// TODO: Check this for correctness!
// TODO: Try to use this interface whenever possible
uint32_t DGPW::Model(uint32_t var) const {
  //      std::cout << __PRETTY_FUNCTION__ << "  y  " << std::endl;
  //    assert(var>0);
  //      std::cout << "_solver->GetLastResult(): " << _solver->GetLastResult()
  //      << std::endl;
  if (_lastResult != 10) return 0;

  assert(var <= Variables());
  uint32_t value = _solver->GetModel(static_cast<int>(var));
  assert(value <= 2 * Variables() + 1);
  assert(value <= 2 * Variables() + 1);
  return value;
}

uint32_t DGPW::Variables() const { return _solver->GetNumberOfVariables(); }

uint32_t DGPW::Clauses() const { return _solver->GetNumberOfClauses(); }

uint32_t DGPW::StaticClauses() const {
  //    return static_cast<uint32_t>(_glucose->nClauses());
  return Clauses();
}

uint32_t DGPW::CurrentBinaryClauses() const {
  // TODO: implement me!
  return 0;
}
uint32_t DGPW::CurrentTernaryClauses() const {
  // TODO: implement me!
  return 0;
}

Settings *DGPW::GetSetting() { return _dgpwSetting; }

bool DGPW::HasEncoding() { return _hasEncoding; }

/********** End: Interface to Glucose and reimplemented functions from AntomBase
 * **********/

bool DGPW::AddWeightedSoftClause(std::vector<uint32_t> &clause,
                                 uint64_t weight) {
  SoftClause *sclaus = new SoftClause(0, clause, weight);
  _softClauses.push_back(sclaus);

  return true;
}

void DGPW::SetSoftClauseVector(std::vector<SoftClause *> *softClauses) {
  _softClauses = *softClauses;
}

// If sorter == nullptr, we want to count for all soft clauses
uint64_t DGPW::CountSatisfiedSoftClauses(Sorter *sorter,
                                         const std::vector<uint32_t> &model) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<SoftClause *> softclauses;
  if (sorter == nullptr) {
    softclauses = _softClauses;
  } else {
    softclauses = sorter->GetSoftClauses();
  }

  uint64_t result = CountSatisfiedSoftClauses(softclauses, model);

  if (sorter != nullptr) {
    sorter->SetMinSatisfied(result);
  }
  return result;
}

uint64_t DGPW::CountSatisfiedSoftClauses(std::vector<SoftClause *> softclauses,
                                         const std::vector<uint32_t> &model) {
  if (_dgpwSetting->verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "SC.size: " << softclauses.size() << std::endl;
    //        std::cout << "Clauses.size: " << Clauses() << std::endl;
  }
  uint64_t result(0);
  // Proceed all soft clauses
  for (uint32_t i = 0; i != softclauses.size(); ++i) {
    uint32_t relaxlit = softclauses[i]->relaxationLit;

    //    if( model[relaxlit>>1] == relaxlit )
    if (Model(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause(softclauses[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (Model(clause[pos] >> 1) == clause[pos]) {
          result += softclauses[i]->weight;
          softclauses[i]->lastassignment = relaxlit ^ 1;
          //            std::cout << "++" << result;
          break;
        }
      }
      if (pos == clause.size()) {
        softclauses[i]->lastassignment = relaxlit;
      }
    } else if (Model(relaxlit >> 1) != 0) {
      assert(Model(relaxlit >> 1) == (relaxlit ^ 1));
      softclauses[i]->lastassignment = relaxlit ^ 1;
      result += softclauses[i]->weight;
      //        std::cout << " ++relax sat " << result;
    }
  }
  return result;
}

// Weigthed maxsat
// ____________________________________________________________________________________________________________________
uint32_t DGPW::MaxSolveIncremental() {
  IncrementalReset();

  // _satWeight = 0;
  // _minWeight = 1;
  // _maxWeight = 0;

  uint32_t currentresult = Solve();

  if (currentresult != SATISFIABLE) return currentresult;

  CalculateOverallOptimum(_satWeight, true);

  currentresult = _mainCascade->Solve(true, true);

  return currentresult;
}

// Weigthed maxsat
// ____________________________________________________________________________________________________________________
uint32_t DGPW::MaxSolveWeightedPartial(int64_t &optimum) {
  std::vector<uint32_t> externalAssumptions = {};
  uint32_t result = MaxSolveWeightedPartial(externalAssumptions, optimum);

  return result;
}

// ____________________________________________________________________________________________________________________
uint32_t DGPW::MaxSolveWeightedPartial(
    std::vector<uint32_t> &externalAssumptions, int64_t &optimum) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  //    _dgpwSetting->Print();

  // to change optimum value with another function called from cascade and
  // bucket
  _optimum = &optimum;
  _hasEncoding = false;
  //    _dgpwSetting->application = WEIGHTEDMAXSAT;

  _timeVariables = new TimeVariables();
  TimeMeasurement totalTime(&_timeVariables->total, true);

  uint32_t currentresult(UNKNOWN);
  if (_dgpwSetting->checkIfSolutionIsUnique &&
      _dgpwSetting->currentCascade.iteration > 0) {
    if (_dgpwSetting->currentCascade._onlyWithAssumptions) {
      std::cout << "c BAD PARAMETER COMBINATION!!!" << std::endl;
      exit(1);
    }
    _dgpwSetting->currentCascade._onlyWithAssumptions = true;
  }

  // uint32_t averageBucketEntries(0);

  //#ifndef NDEBUG
  uint64_t softWeights = 0;
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { softWeights += SC->weight; });
  _sumOfSoftWeights = softWeights;

  //	std::cout << "soft weights: " << softWeights << " sumofweights: " <<
  //_sumOfSoftWeights << std::endl;
  //    assert(softWeights == _sumOfSoftWeights);
  //#endif
  if (_sumOfSoftWeights == 1) {
  }

  //  std::cout << "c clauses before.........: " << _clausesBefore << std::endl;
  if (_dgpwSetting->verbosity > 0) {
    std::cout << "c #softclauses...........: " << _softClauses.size()
              << std::endl;
    std::cout << "c #sum of SoftWeights....: "
              << _sumOfSoftWeights * _greatestCommonDivisor << std::endl;
  }
  //    if (!_dgpwSetting->solveAtFirst)
  //      std::cout << "c #sum of SoftWeights....: " << _sumOfSoftWeights *
  //      _greatestCommonDivisor << std::endl;

  if (_dgpwSetting->onlyByTares ||
      _dgpwSetting->mcDivideStrategy != SOLVEINNORMALCASCADEMODE) {
    _dgpwSetting->encodeStrategy = ENCODEONLYIFNEEDED;
  }

  if (_dgpwSetting->solveAtFirst) {
    if (false) {
      _approxSolverCalls = _solverCalls;
      TimeMeasurement solvedFirst(&_timeVariables->solvedFirst);

      if (_greedyPrepro == 0) {
        // old variant
        currentresult = _solver->Solve();
      }
      if (!_fixedAssumptions.empty()) {
        return SATISFIABLE;
      }
      CalculateOverallOptimum(_satWeight, true);
      // now done in pacose
    } else {
      currentresult = SATISFIABLE;
    }

    if (currentresult == UNSAT) {
      _satWeight = 0;
      std::cout << "c first Solver call is UNSAT!" << std::endl;
      std::cout << "UNSAT" << std::endl;
      if (!_dgpwSetting->currentCascade._onlyWithAssumptions) {
        return UNSAT;
      }
    } else if (currentresult == UNKNOWN) {
      std::cout << "c UNKNOWN!!!!!! shouldn't be!" << std::endl;
      return UNKNOWN;
    }

    if (_softClauses.size() == 1 && _satWeight == 0 &&
        _sumOfSoftWeights == _softClauses[0]->weight) {
      // try to solve the only remaining softclause!
      std::vector<uint32_t> assumptions;
      assumptions.push_back(_softClauses[0]->relaxationLit ^ 1);
      currentresult = Solve(assumptions);
      if (currentresult == UNSAT) {
        std::cout << "c The only SC is NOT satisfiable!" << std::endl;
        AddUnit(_softClauses[0]->relaxationLit);
        if (Solve() != 10) {
          std::cout << "Strange solver result, shouldn't be possible!"
                    << std::endl;
        }
        CalculateOverallOptimum(_satWeight, true);
        return SATISFIABLE;
      } else if (currentresult == SATISFIABLE) {
        std::cout << "c The only SC is SATISFIABLE!" << std::endl;
        AddUnit(_softClauses[0]->relaxationLit ^ 1);
        CalculateOverallOptimum(_satWeight, true);
        return SATISFIABLE;
      }
    }

    //    std::cout << "c SATWeight solved first.: " << _satWeight << std::endl;
    if (_softClauses.size() == 0) {
      return SATISFIABLE;
    }
  }

  if (_satWeight == _sumOfSoftWeights) {
    std::cout << "o 0";
    std::cout << "c All SoftClauses are Satisfiable!" << std::endl;
    if (!_dgpwSetting->currentCascade._onlyWithAssumptions ||
        _dgpwSetting->checkIfSolutionIsUnique) {
      // add relaxation literals as UnitClauses!
      for (auto sc : _softClauses) {
        AddUnit(sc->relaxationLit ^ 1);
      }
    }
    return SATISFIABLE;
  }

  // calc and set all trivially conflicting softclauses
  // TOASK: maybe better seperately for each bucket - to know boundaries of
  // bucket! together with mode of solving bucket parts to get max
  // CheckAllWeightedConflictingSoftclauses();
  if (_dgpwSetting->mcDivideStrategy != SOLVEINNORMALCASCADEMODE) {
    _dgpwSetting->encodeStrategy = ENCODEONLYIFNEEDED; // cone of influence encoding
    _mainMultipleCascade = new MultipleCascade(
        this, _dgpwSetting->onlyByTares, _dgpwSetting->tareCascadeOnlyByTares,
        _dgpwSetting->cascadeDivider, _dgpwSetting->maxBucketSize,
        _dgpwSetting->nOfCasc, _dgpwSetting->interimResult, _sumOfSoftWeights);
    //        std::cout << "Cascade Divider: " << _cascadeDivider << std::endl;
    if (_mainMultipleCascade->DivideAndConnectCascades(
            _dgpwSetting->mcDivideStrategy, &_softClauses)) {
      currentresult = _mainMultipleCascade->Solve();

      if (currentresult != UNKNOWN) currentresult = SATISFIABLE;
    }
  } else {
    _dgpwSetting->cascadeDivider = 0;
  }

  // no else if, because if _cascadeDivider returns only one cascade then
  // _cascadeDivider == 0
  if (_dgpwSetting->cascadeDivider == 0) {
    _clausesBefore = StaticClauses();
    _mainMultipleCascade = nullptr;
    //        std::cout << "_dgpwSetting->onlyByTares: " <<
    //        _dgpwSetting->onlyByTares << std::endl;
    _mainCascade = new Cascade(this, nullptr, _dgpwSetting->onlyByTares);

    //        std::cout << "PartitionStrategy: " <<
    //        _dgpwSetting->partitionStrategy << std::endl;

    _mainCascade->Fill(&_softClauses, _dgpwSetting->partitionStrategy,
                       _dgpwSetting->encodeStrategy);
    _hasEncoding = true;

    // if (_antomSetting->verbosity > 2 && _antomSetting->base > 2)
    //{
    // std::cout << "How often is Softclause in Bucket, because of base " <<
    // _baseMode << std::endl; DumpBucketStructure();
    //}

    if (_dgpwSetting->cascadeDivider == 0) {
      if (_dgpwSetting->encodeStrategy == ENCODEALL) {
        _clausesBefore = StaticClauses();
      }
      if (!_mainCascade->Encode()) {
        return UNKNOWN;
      }
      //      if (_dgpwSetting->encodeStrategy == ENCODEALL) {
      //        std::cout << "c #clauses of coding.....: "
      //                  << StaticClauses() - _clausesBefore << std::endl;
      //        std::cout << "c #variables of coding...: "
      //                  << Variables() - _variablesBefore << std::endl;
      //      }
    }

    // TOBI: SetDecisionStrategiesForCascadeModel!!!!
    // SetDecisionStrategiesForMaxSAT();
    currentresult =
        _mainCascade->Solve(_dgpwSetting->currentCascade._onlyWithAssumptions,
                            _dgpwSetting->currentCascade._solveTares);
  }

  //  if (_dgpwSetting->encodeStrategy != ENCODEALL &&
  //  _dgpwSetting->onlyByTares) {
  //    std::cout << "c #clauses of coding.....: " << _addedClauses <<
  //    std::endl; std::cout << "c #variables of coding...: " << _addedVariables
  //    << std::endl;
  //  }
  if (currentresult == UNKNOWN) {
    return UNKNOWN;
  }

  // if there is more than one cascade, we have to calculate
  // the OverallOptimum again to sum up the weights.
  //    CalculateOverallOptimum(_satWeight, true);
  optimum = optimum * _greatestCommonDivisor;

  if (_dgpwSetting->checkIfSolutionIsUnique &&
      _dgpwSetting->currentCascade.iteration > 0) {
    if (optimum == 0) {
      // special case 1 - all SCs are satisfied!
      for (auto SC : _softClauses) {
        AddUnit(SC->relaxationLit ^ 1);
      }
      //            currentresult = Solve();
      //            if (currentresult != SATISFIABLE) {
      //                std::cout << "c Strange result - shouldn't be possible!"
      //                << std::endl; exit(1);
      //            }
    } else if (_satWeight == 0) {
      // special case 2 - no SC could be satisfied
      for (auto SC : _softClauses) {
        AddUnit(SC->relaxationLit);
      }
      //            currentresult = Solve();
      //            if (currentresult != SATISFIABLE) {
      //                std::cout << "c Strange result - shouldn't be possible!"
      //                << std::endl; exit(1);
      //            }
    } else {
      std::vector<uint32_t> clause;
      for (auto SC : _softClauses) {
        clause.push_back(Model(SC->relaxationLit >> 1) ^ 1);
      }
      uint32_t relaxLit = NewVariable() << 1;
      clause.push_back(relaxLit);
      AddClause(clause);

      std::vector<uint32_t> assumptions = GetLastAssumptions();

      if (assumptions.size() == 0) {
        std::cout << "c NO ASSUMPTIONS - ERROR!!!" << std::endl;
        exit(1);
      }
      assumptions.push_back(relaxLit ^ 1);

      uint32_t currentresult = Solve(assumptions);

      if (currentresult != SATISFIABLE) {
        std::cout << "c SOLUTION IS UNIQUE! Add this solution as hard clauses "
                     "to the problem!"
                  << std::endl;
        AddUnit(relaxLit);
        clause.pop_back();
        // add solution as unit clauses
        for (auto lit : clause) {
          AddUnit(lit ^ 1);
        }
        // we don't need the encoding anymore - throw it away!!
        for (auto ass : assumptions) {
          AddUnit(ass ^ 1);
        }
        currentresult = Solve();
        if (currentresult != SATISFIABLE) {
          std::cout << "c Strange result - shouldn't be possible!" << std::endl;
          exit(1);
        }
      } else {
        std::cout << "c THERE ARE MULTIPLE SOLUTIONS!" << std::endl;
        assumptions.pop_back();
        AddUnit(relaxLit);
        for (auto ass : assumptions) {
          std::cout << ass << ", ";
          AddUnit(ass);
        }
      }
    }
  }

  if (_dgpwSetting->verbosity > 0) {
    if (!_dgpwSetting->formulaIsDivided) {
      std::cout << "s SATISFIABLE" << std::endl;
      //    std::cout << "o " << optimum << std::endl;
    } else {
      if (currentresult == SATISFIABLE) {
        std::cout << "c currently SAT" << std::endl;
        //      std::cout << "c local o " << optimum << std::endl;
      } else if (currentresult == 20) {
        std::cout << "c currently UNSAT" << std::endl;
        //      std::cout << "c local o " << optimum << std::endl;
      }
    }
  }

  if (_satWeight == _sumOfSoftWeights &&
      !_dgpwSetting->currentCascade._onlyWithAssumptions) {
    // add relaxation literals as UnitClauses!
    for (auto sc : _softClauses) {
      AddUnit(sc->relaxationLit ^ 1);
    }
    if (_dgpwSetting->verbosity > 0)
      std::cout << "c All SoftClauses are Satisfiable!" << std::endl;
    return SATISFIABLE;
  }
  if (_dgpwSetting->verbosity > 0)
    std::cout << "c number of variables: " << Variables() << std::endl;
  return currentresult;

}  // namespace DGPW

void DGPW::SetInitialAssumptions(std::vector<uint32_t> assumptions) {
  _externalAssumptions = assumptions;
}

StructureInfo DGPW::AnalyzeStructure() {
  // assert(_sumOfSoftWeights > _topWeight);
  bool maxIsMin = (_maxWeight == _minWeight);
  bool withHC = Clauses() != 0;
  bool minIsSet = _minWeight != ((uint64_t)-1);
  bool maxIsSet = _maxWeight != 0;
  // is SAT formula
  bool onlyHC = !minIsSet && !maxIsSet && withHC;
  // is MaxSAT formula -- NOT PARTIAL
  bool onlySCwithOneWeight = !_moreThanTwoWeights && maxIsMin && !withHC;
  // is partial MaxSAT formula - Set _minWeight to one
  bool oneHCWeightAndSCWeight = maxIsMin && withHC;
  // if sum (_minWeight) < _maxWeight --> is partial MaxSAT formula
  bool onlySCwithTwoWeights =
      !_moreThanTwoWeights && !maxIsMin && minIsSet && maxIsSet && !withHC;
  // Calculate the greatest common divisor of all softclauseweights.
  uint64_t greatestCommonDivisor = _minWeight;

  for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
    if (greatestCommonDivisor == 1) break;
    greatestCommonDivisor =
        GreatestCommonDivisor(greatestCommonDivisor, _softClauses[ind]->weight);
  }
  _greatestCommonDivisor = greatestCommonDivisor;

  if (_dgpwSetting->verbosity > 2) {
    std::cout << std::setw(30) << "_topWeight " << _topWeight << std::endl;
    std::cout << std::setw(30) << "_minWeight " << _minWeight << std::endl;
    std::cout << std::setw(30) << "_maxWeight " << _maxWeight << std::endl;
    std::cout << std::setw(30) << "_moreThanTwoWeights " << _moreThanTwoWeights
              << std::endl;
    std::cout << std::setw(30) << "maxIsMin " << maxIsMin << std::endl;
    std::cout << std::setw(30) << "withHC " << withHC << std::endl;
    std::cout << std::setw(30) << "minIsSet " << minIsSet << std::endl;
    std::cout << std::setw(30) << "maxIsSet " << maxIsSet << std::endl;
    std::cout << std::setw(30) << "onlyHC " << onlyHC << std::endl;
    std::cout << std::setw(30) << "onlySCwithOneWeight " << onlySCwithOneWeight
              << std::endl;
    std::cout << std::setw(30) << "oneHCWeightAndSCWeight "
              << oneHCWeightAndSCWeight << std::endl;
    std::cout << std::setw(30) << "onlySCwithTwoWeights "
              << onlySCwithTwoWeights << std::endl;
    std::cout << std::setw(30) << "_sumOfSoftWeights " << _sumOfSoftWeights
              << std::endl;
  }

  if (onlyHC) {
    return ISSAT;
  }
  if (onlySCwithTwoWeights) {
    uint32_t sumOfMinWeights(0);

    // the sum of all _minWeights has to be smaller than _maxWeight
    // then it can be solved like MaxSAT formula.
    for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
      if (_softClauses[ind]->weight != _minWeight) continue;

      sumOfMinWeights += _softClauses[ind]->weight;
    }

    if (sumOfMinWeights < _maxWeight) {
      _sumOfSoftWeights = sumOfMinWeights;
      return CONVERTTOMAXSAT;
    }
  }
  if (oneHCWeightAndSCWeight || onlySCwithOneWeight || oneHCWeightAndSCWeight) {
    return ISMAXSAT;
  }
  //    if (greatestCommonDivisor > 1)
  //      {
  //        _greatestCommonDivisor = greatestCommonDivisor;
  //        std::cout << "c greatest common divisor: " <<
  //        greatestCommonDivisor
  //        << std::endl; return dgpw::DIVWEIGHTSBYDIVISOR;
  //      }
  return ISWEIGHTEDMAXSAT;
}

uint64_t DGPW::GreatestCommonDivisor(uint64_t a, uint64_t b) {
  uint64_t temp;
  while (b > 0) {
    temp = b;
    b = a % b;
    a = temp;
  }
  return a;
}

void DGPW::ConvertFormulaToMaxSAT(uint64_t maxWeight) {
  // configure all _maxWeights as hardclauses
  std::vector<SoftClause *> newSoftClauseVector;
  for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
    if (_softClauses[ind]->weight != maxWeight) {
      newSoftClauseVector.push_back(_softClauses[ind]);
      continue;
    }

    AddClause(_softClauses[ind]->clause);
  }

  _softClauses = newSoftClauseVector;
}

void DGPW::DivideAllSoftClausesByFactor(uint64_t factor) {
  for (uint32_t ind = 0; ind < _softClauses.size(); ++ind) {
    _softClauses[ind]->weight /= factor;
  }
  _sumOfSoftWeights /= factor;
}

uint64_t DGPW::CalcGlobalOpt() { return _pacose->CalculateSATWeight(); }

uint64_t DGPW::CalculateOverallOptimum(uint64_t satWeight, bool countAgain) {
  if (_dgpwSetting->verbosity > 3)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  // struct rusage resources;
  // Cascade and Bucket counts the satisfied SC only for their SC's
  // for an overall result count again.
  if (countAgain) {
    // maybe _maxSatWeight - because it could be a local optima of this Solver
    // call
    satWeight = CountSatisfiedSoftClauses(nullptr, {});
  }

  // std::cout << std::setw(80) <<std::right<< "OPTIMUM: " << *_optimum << std::endl;

  if (*_optimum > static_cast<int64_t>(_sumOfSoftWeights - satWeight) ||
      *_optimum == -1) {
    /*
      if (_antomSetting->verbosity > 0)
      {
      getrusage(RUSAGE_SELF, &resources);
      double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double)
      resources.ru_utime.tv_usec; std::cout << "c " << (timeC -
      _control->GetStartTime()) << "s" << std::endl;
      }
    */
    // better actualize satWeight once at the end.
    _satWeight = satWeight;
    //    std::cout << "_sumOfSoftWeights: " << _sumOfSoftWeights
    //              << "  satWeight: " << satWeight << std::endl;
    *_optimum = static_cast<int64_t>(_sumOfSoftWeights - satWeight);


    // std::cout <<std::setw(80) << std::right <<  "FORMULA IS DIVIDED? " << _dgpwSetting->formulaIsDivided << std::endl;
    if (!_dgpwSetting->formulaIsDivided) {
      std::cout << "o " << *_optimum * _greatestCommonDivisor << std::endl;
    } else {
      _pacose->CalculateSATWeight();
    }

    if (_dgpwSetting->verbosity > 0)
      std::cout << std::setw(50) << "Calculated Global SATWeight: " << satWeight
                << std::endl;

  } else if (_dgpwSetting->verbosity > 2) {
    /*
      if (_antomSetting->verbosity > 3)
      {
      getrusage(RUSAGE_SELF, &resources);
      double timeC = (double) resources.ru_utime.tv_sec + 1.e-6 * (double)
      resources.ru_utime.tv_usec;
      //std::cout << "c " << (timeC - _control->GetStartTime()) << "s" <<
      std::endl;
      }
    */
    if (!_dgpwSetting->formulaIsDivided)
      std::cout << "o " << *_optimum * _greatestCommonDivisor << std::endl;
    std::cout << std::setw(50) << "Calculated Global SATWeight: " << satWeight
              << std::endl;
  }

  return _satWeight;
}

void DGPW::SetMoreThanTwoWeights(bool val) { _moreThanTwoWeights = val; }

bool DGPW::GetHasMoreThanTwoWeights() { return _moreThanTwoWeights; }

void DGPW::SetTopWeight(uint64_t val) {
  assert(val > 0);
  _topWeight = val;
}

void DGPW::SetMinWeight(uint64_t val) {
  // assert( val < (uint64_t) - 1 );
  _minWeight = val;
}

uint64_t DGPW::GetMinWeight() { return _minWeight; }

void DGPW::SetMaxWeight(uint64_t val) {
  // assert( val > 0 );
  _maxWeight = val;
}

uint64_t DGPW::GetMaxWeight() { return _maxWeight; }

uint64_t DGPW::GetSATWeight() { return _satWeight; }

std::vector<std::pair<uint64_t, uint32_t>> DGPW::GetTareVector(
    uint64_t weightDiff) {
  if (_mainCascade == nullptr) return {};

  return _mainCascade->GetTareVector(weightDiff);
}

std::vector<std::pair<uint64_t, uint32_t>> DGPW::GetWatchdogs(
    uint64_t weightDiff) {
  if (_mainCascade == nullptr) {
    return {};
  }

  return _mainCascade->GetWatchdogs(weightDiff);
}

std::vector<uint32_t> DGPW::GetLastAssumptions() {
  if (_mainCascade == nullptr) {
    return {};
  }
  return _mainCascade->GetLastAssumptions();
}

void DGPW::SetFixedAssumptions(std::vector<uint32_t> fixedAssumptions) {
  _fixedAssumptions = fixedAssumptions;
  //  if (_dgpwSetting->verbosity > 0) std::cout << "c fixed added assumptions:
  //  "; for (auto fa : fixedAssumptions) {
  //    std::cout << fa << " ";
  //  }
  //  std::cout << std::endl;
}

void DGPW::RemoveFixedAssumptions() { _fixedAssumptions.clear(); }

void DGPW::SetHasHardClauses(bool val) { _hasHardClauses = val; }

void DGPW::SetSatWeight(uint64_t val) {
  _satWeight = val;
  //  _minWeight = val;
}

bool DGPW::GetHasHardClauses() { return _hasHardClauses; }

void DGPW::SetGreatestCommonDivisor(uint64_t val) {
  _greatestCommonDivisor = static_cast<int64_t>(val);
}
}  // namespace DGPW
} // Namespace Pacose
