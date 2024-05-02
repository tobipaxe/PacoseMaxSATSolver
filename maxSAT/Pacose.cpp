/********************************************************************************************
SATSolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian
    Parts are taken from the QMaxSAT 2017 MaxSAT evaluation version.

Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
********************************************************************************************/

#include "Pacose.h"
#include "Settings.h"

// to use original antom code
#include "../Helper/ClauseDB.h"
// maxpre2
// #include "../maxpre2/src/preprocessorinterface.hpp"
#include "../solver-proxy/SATSolverProxy.h"
#include "DGPW/dgpw.h"
#include "Encodings.h"
#include "Greedyprepro.h"
#include "Pacose.h"
#include "Softclause.h"
#include "DGPW/timevariables.h"

#include <iostream> // std::cout
#include <math.h>   // std::cout
// #include <sstream>      // std::setw
#include <algorithm>      // std::sort
#include <cassert>        // assertions
#include <cstring>        // fast parser
#include <fstream>        // parseWcnfFile std::ifstream
#include <iomanip>        // std::setw
#include <set>            // comparison of weight difference
#include <sstream>        // parseWcnfFile std::istringstream
#include <stdlib.h>       // strtoull
#include <sys/resource.h> // timing
#include <sys/time.h>     // timing
#include <numeric>        // std::accumulate

namespace Pacose {

Pacose::Pacose()
    : _settings(), _satSolver(nullptr), _nbVars(0), _nbClauses(0),
      _nbOriginalClauses(0), _nbUnitSoftClausesAdded(0), _top(0), _originalSoftClauses(),
      _actualSoftClauses(nullptr), _hasHardClauses(false),
#ifdef SaveCNF
      _CNF({}),
#endif // SaveCNF
      _cpuLimit(INT32_MAX), _memLimit(INT32_MAX), _nbOfOrigVars(0),
      _sumOfSoftWeights(0), _sumOfActualSoftWeights(0), _overallSoftWeights(0), _satWeight(0),
      _unSatWeight(INT64_MAX), _lastCalculatedUnsatWeight(-1),
      _currentLocalUnSatWeight(INT64_MAX),_currentLocalSatWeight(0),
      _localUnSatWeight(INT64_MAX), _minWeight(0), _maxWeight(0), _GCD(1),
      _alwaysSATSCs(0), _alwaysUNSATSCs(0), _alwaysSATWeight(0),
      _alwaysUNSATWeight(0), _trimSATTime(0), _noTrimSAT(0),
      _noTrimSATSolverCalls(0), _GBMOPartitions(0), _GBMOTime(0),
      _variablesOfEncoding(0), _clausesOfEncoding(0), _noSolverCalls(0),
      _encodings(), _negRelaxLit(0) {
}

Pacose::~Pacose() {
  if (_satSolver != nullptr) {
    delete _satSolver;
  }
  //  delete _actualSoftClauses;
  for (auto sc : _originalSoftClauses)
    delete sc;
  for (auto cand : _cascCandidates) {
    if (cand.dgpw != nullptr)
      delete cand.dgpw;
  }

  delete _encodings;
}

void Pacose::InitSatSolver(int solver) {
  if (solver == 0) {
    InitSatSolver(SATSolverType::GLUCOSE421);
  } else if (solver == 1) {
    InitSatSolver(SATSolverType::GLUCOSE3);
  } else if (solver == 2) {
    InitSatSolver(SATSolverType::CADICAL);
  } else if (solver == 3) {
    InitSatSolver(SATSolverType::CRYPTOMINISAT);
  }
}

void Pacose::InitSatSolver(SATSolverType solverType) {
  //  std::cout << __PRETTY_FUNCTION__ << std::endl;
  _satSolver = SATSolverProxy::InitSATSolver(solverType, 1);
  //    _satSolver = SATSolverProxy();
  //            SATSolverProxy::InitSATSolver(solverType, 1);
  if (_settings.reconf != 99) {
    _satSolver->SetReconf(_settings.reconf);
  }
}

// Wrapper Functions for SolverProxy
bool Pacose::AddClause(std::vector<uint32_t> &clause) {
#ifdef SaveCNF
  _CNF.push_back(clause);
#endif // SaveCNF
  bool rv = _satSolver->AddClause(clause);
  return rv;
}

bool Pacose::DumpWCNF() {
#ifdef SaveCNF
  struct timeval tp;
  gettimeofday(&tp, NULL);
  long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;

  std::string filename = std::to_string(ms) + ".wcnf";

  uint32_t maxValue = 0;
  uint64_t weight = 0;
  for (auto clause : _CNF) {
    for (auto lit : clause) {
      if ((lit >> 1) > maxValue)
        maxValue = (lit >> 1);
      assert((lit >> 1) != 0);
    }
  }
  for (auto sc : _originalSoftClauses) {
    for (auto lit : sc->clause) {
      if ((lit >> 1) > maxValue)
        maxValue = (lit >> 1);
      assert((lit >> 1) != 0);
    }
    weight += sc->weight;
  }
  _nbVars = maxValue;
  _top = weight + 1;

  std::ofstream out(filename);
  //  out.open(filename, std::ofstream::out | std::ofstream::app);
  std::streambuf *coutbuf = std::cout.rdbuf(); // save old buf
  std::cout.rdbuf(out.rdbuf()); // redirect std::cout to filename!

  std::cout << "p wcnf " << _nbVars << " " << _CNF.size() << " " << _top
            << std::endl;

  for (auto hc : _CNF) {
    assert(hc.size() > 0);
    std::cout << _top << " ";
    for (auto lit : hc) {
      assert((lit >> 1) != 0);
      if (lit % 2 == 0)
        std::cout << (lit >> 1) << " ";
      else
        std::cout << "-" << (lit >> 1) << " ";
    }

    std::cout << "0" << std::endl;
  }

  for (auto sc : _originalSoftClauses) {
    assert(sc->clause.size() > 0);
    assert(sc->weight > 0);
    std::cout << sc->weight << " ";
    for (auto lit : sc->clause) {
      assert((lit >> 1) != 0);
      if (lit % 2 == 0)
        std::cout << (lit >> 1) << " ";
      else
        std::cout << "-" << (lit >> 1) << " ";
    }
    std::cout << "0" << std::endl;
  }

  std::cout.rdbuf(coutbuf); // reset to standard output again

  std::cout << "c Clauses are written to WCNF file: " << filename << std::endl;
  return true;
#endif // SaveCNF
  return false;
}

uint32_t Pacose::SignedTouint32_tLit(int literal) {
  // if literal < 0 then lit=(2*literal)+1 else lit=2*literal

  return (static_cast<uint32_t>(abs(literal)) << 1) ^ (literal < 0);
}

void Pacose::AddSoftClause(std::vector<uint32_t> &clause, uint64_t weight) {
  uint32_t relaxLit;
  if (clause.size() == 1) {
    relaxLit = clause[0] ^ 1;
  } else {
    assert(clause.size() > 1);
    relaxLit = static_cast<uint32_t>(_satSolver->NewVariable() << 1);
  }

  SoftClause *SC = new SoftClause(relaxLit, clause, weight);
  _originalSoftClauses.push_back(SC);
  if (clause.size() != 1) {
    clause.push_back(relaxLit);
    _satSolver->AddClause(clause);
  }
}

void Pacose::AddSoftClauseTo(std::vector<SoftClause *> *softClauseVector,
                             std::vector<uint32_t> &clause, uint64_t weight) {
  uint32_t relaxLit = static_cast<uint32_t>(_satSolver->NewVariable() << 1);

  SoftClause *SC = new SoftClause(relaxLit, clause, weight);
  softClauseVector->push_back(SC);
  clause.push_back(relaxLit);

  if (!_satSolver->AddClause(clause)) {
    assert(false);
    //        exit(1);
  }
}

// koshi 2013.04.05, 2013.05.21, 2013.06.28, 2013.07.01, 2013.10.04
// koshi 20140121
void Pacose::HeuristicQMaxSAT(long long int sum, long long int k) {
  printf("c QMaxSAT17 auto-mode for generating cardinality constraints\n");
  int logk = 0;
  int logsum = 0;
  // calculate 2nd logarithm of k
  for (long long int ok = k; ok > 0; ok = ok >> 1)
    logk++;
  // calculate 2nd logarithm of sum of weights
  for (long long int osum = sum; osum > 0; osum = osum >> 1)
    logsum++;
  printf("c logk = %d, logsum = %d\n", logk, logsum);
  if (logk + logsum < 15) {
    // Bailleux
    _settings.SetEncoding(BAILLEUX);
    _settings.SetCompression(0);
    printf("c Bailleux's encoding (comp=0)\n");
  } else if (k < 3) {
    // Warners
    _settings.SetEncoding(WARNERS);
    _settings.SetCompression(1);
    printf("c Warners' encoding (comp=1)\n");
  } else if (logsum < 17) {
    // Ogawa
    _settings.SetEncoding(OGAWA);
    _settings.SetCompression(0);
    printf("c Ogawa's encoding (comp=0)\n");
  } else {
    _settings.SetEncoding(WARNERS);
    _settings.SetCompression(1);
    printf("c Warners' encoding (comp=1)\n");
  }
}

void Pacose::wbSortAndFilter() {
  if (_settings.verbosity > 5)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  assert(std::cout << "c assertion Solver call in Pacose, wbsortandfilter1" << std::endl && _satSolver->Solve() == SAT);
  
  if (_settings.createSimplifiedWCNF)
    return;

  CalculateSATWeight();
  uint32_t softClausesBefore = 0;
  uint32_t softClausesAfter = 0;

  for (uint32_t i = static_cast<uint32_t>(_sClauses.size()); i > 0; i--) {
    uint32_t clausesBefore = _sClauses[i - 1].size();
    softClausesBefore += _sClauses[i - 1].size();

    // std::cout << "_sClauses[i - 1].size() " << static_cast<uint32_t>(_sClauses[i - 1].size()) << "  _sClauses.size(): " << static_cast<uint32_t>(_sClauses.size()) << std::endl;
    if (_sClauses[i - 1].size() > 0)
      wbSortAndFilter(_sClauses[i - 1]);
    softClausesAfter += _sClauses[i - 1].size();
    if (_settings.verbosity > 1)
      std::cout << "Removed SCs Round " << i << ": " << clausesBefore - softClausesAfter << std::endl;
  }

  if (softClausesAfter < softClausesBefore || _settings.verbosity > 10) {
    if (_settings.verbosity > 0) {
      std::cout << "c removed SCs by wb......: "
                << softClausesBefore - softClausesAfter << std::endl;
      std::cout << "c remaining SCs..........: "
                << softClausesAfter << std::endl;
      // std::cout << _originalSoftClauses[0]->weight << std::endl;
      // std::cout << _originalSoftClauses[0]->relaxationLit << std::endl;
    }
  }
  assert(std::cout << "c assertion Solver call in Pacose, wbsortandfilter2" << std::endl && _satSolver->Solve() == SAT);

  #ifdef DEBUG
    if (_satSolver->Solve() != SAT) {
      std::cout << "ERROR: Solver call should've been satisfiable!"
                << std::endl;
    }

    CalculateSATWeight(); 
  #endif

}

void Pacose::wbSortAndFilter(std::vector<SoftClause *> & softClauseVector) {
  if (_settings.verbosity > 5)
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  bool currentSCs = false;

  if (softClauseVector.size() == (*_actualSoftClauses).size() && 
      softClauseVector[0]->relaxationLit == (*_actualSoftClauses)[0]->relaxationLit && 
      softClauseVector[0]->originalWeight == (*_actualSoftClauses)[0]->originalWeight) {
    std::cout << "c processing actualSoftClause vector" << std::endl;
    currentSCs = true;
  }
    
  for (uint32_t i = 0; i < softClauseVector.size(); i++) {
    uint64_t currentWeight = softClauseVector[i]->originalWeight;;
    uint64_t currentUnsatweight = _unSatWeight;
    if (currentSCs) {
      currentWeight = softClauseVector[i]->weight;
      currentUnsatweight = _localUnSatWeight;
    }
    if (_settings.verbosity > 2) {
      std::cout << "i: " << i << " size: " << softClauseVector.size() << std::endl;
      std::cout << "currentWeight: " << currentWeight << "  currentUnsatweight: " << currentUnsatweight << std::endl;
    }
    if (!(currentWeight <= currentUnsatweight)) {
      if (currentSCs) {
        assert(_localSatWeight >= currentWeight);
        _localSatWeight -= currentWeight;
        _sumOfActualSoftWeights -= currentWeight;
      }

      // SC has to be satisfied in any case!
      _alwaysSATSCs++;
      _alwaysSATWeight += softClauseVector[i]->originalWeight;
      // TODO Dieter:  This clause is added because of the fact that the objective literal as higher weight than the current optimal solution (minimizing). 
      // Same reasoning applies. What is to do is the rewriting of the objective function.
      _satSolver->ResetClause();
      _satSolver->NewClause();
      uint32_t ulit = softClauseVector[i]->relaxationLit ^ 1;
      // std::cout << "OrigWeight: " << softClauseVector[i]-> originalWeight << " Add Ulit: " << ulit << std::endl;

      _satSolver->AddLiteral(&ulit);
      _satSolver->CommitClause();
      softClauseVector.erase(softClauseVector.begin() + i);
      i--;
    }
  }

}

void Pacose::genCardinals(
    long long int tmpUnSATWeight,
    long long int &divisor, // koshi 2013.10.04
    std::vector<uint32_t> &lits, std::vector<uint32_t> &linkingVar,
    std::vector<long long int> &linkingWeight, // uemura 20161202
    std::vector<long long int> &divisors,      // uemura 20161128
    std::vector<std::vector<uint32_t>> &linkingVars,
    std::vector<std::vector<long long int>> &linkingWeights, // uemura 20161128
    int compression) {
  // simple inprocessor, after first solve call
  // filter out weights which cannot be fulfilled anymore!

   

  if (_settings.verbosity > 0)
    std::cout << "c Sum of weights = " << _sumOfActualSoftWeights
              << " _localUnSatWeight: " << _localUnSatWeight << std::endl;
  
  assert(_localUnSatWeight == tmpUnSATWeight);
  
  if (_actualSoftClauses->size() == 0) {
    linkingVar.clear();
  } else {
    // original QMaxSAT Encodings need weights and blockings for their
    // recursive calls
    //      std::cout << "NOT DGPW" << std::endl;
    _blockings.clear();
    _weights.clear();
    for (size_t i = 0; i < _actualSoftClauses->size(); i++) {
      _blockings.push_back((*_actualSoftClauses)[i]->relaxationLit);
      _weights.push_back((*_actualSoftClauses)[i]->weight);
      //      std::cout << _weights.back() << ", " << _blockings.back() <<
      //      std::endl;
    }
    if (_settings.verbosity > 10) {
      std::cout << std::setw(30) << "Weights and Blockings filled." << std::endl;
      std::cout << std::setw(30) << "blockings size: " << _blockings.size()
                << std::endl;
      std::cout << std::setw(30) << "weights size: " << _weights.size()
                << std::endl;
    }

    //    std::cout << "ASC: " << _actualSoftClauses->size() << std::endl;
    //    std::cout << "BLO: " << _blockings.size() << std::endl;
    //    std::cout << "WEI: " << _weights.size() << std::endl;
    //    std::cout << "SUM: " << sum << std::endl;
    //    std::cout << "USW: " << tmpUnSATWeight << std::endl;
    //    std::cout << "COM: " << _settings.GetCompression() << std::endl;
    //    std::cout << "LTS: " << lits.size() << std::endl;
    //    std::cout << "LVS: " << linkingVar.size() << std::endl;
    //    std::cout << "LWS: " << linkingWeight.size() << std::endl;
    //    std::cout << "LWSS: " << linkingWeights.size() << std::endl;

    // koshi 20140124 20140129
    // koshi 2013.06.28

    //    std::cout << "tmpUnSATWeight: " << tmpUnSATWeight << std::endl;

    // why using old sum of Softweights and not the newly calculated sum?
    switch (_encoding) {
    case WARNERS:
      _encodings->genWarners0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight,
                              compression, *_satSolver, lits, linkingVar);
      break;
    case BAILLEUX:
      _encodings->genBailleux0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight + 1,
                               compression, *_satSolver, lits, linkingVar);
      break;
    case ASIN:
      _encodings->genAsin(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight,
                          compression, *_satSolver, lits, linkingVar);
      break;
    case OGAWA:
      _encodings->genOgawa0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight + 1, divisor,
                            compression, *_satSolver, lits, linkingVar);
      break;
    case BAILLEUXW2:
      _encodings->genBailleuxW20(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight + 1,
                                 compression, *_satSolver, lits, linkingVar,
                                 linkingWeight);
      break;
    case WMTO:
      _encodings->genKWMTO0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight, divisors,
                            *_satSolver, lits, linkingVars, linkingWeights);
      break;
    case MRWTO:
      _encodings->genMRWTO0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight, divisors,
                            *_satSolver, lits, linkingVars, linkingWeights,
                            _encoding);
      break;
    case MRWTO2:
      _encodings->genMRWTO0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight, divisors,
                            *_satSolver, lits, linkingVars, linkingWeights,
                            _encoding);
      break;
    case MRWTO19:
      _encodings->genMRWTO19_0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight,
                               divisors, *_satSolver, lits, linkingVars,
                               linkingWeights, _encoding);
      break;
    case MRWTO19_2:
      _encodings->genMRWTO19_0(_weights, _blockings, _sumOfActualSoftWeights, _localUnSatWeight,
                               divisors, *_satSolver, lits, linkingVars,
                               linkingWeights, _encoding);
      break;
    case DGPW18:
      std::cout << "DGPW18 NOT SUPPORTED!" << std::endl;
      break;
    default:
      std::cout << "Encoding not yet defined!";
      exit(EXIT_FAILURE);
    }
  }
}

uint32_t Pacose::SolveQMax(EncodingType *encodingType) {
  if (_settings.verbosity > 5)
    std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

  if (encodingType == nullptr)
    encodingType = &_encoding;

  _nbOfOrigVars = _satSolver->GetNumberOfVariables();

  _encodings = new Encodings(&_settings);

  // Heuristic to choose which compression of warners!
  int compression = _settings._compression;
  if (compression == -1 && *encodingType == WARNERS) {
    //    assert(*encodingType == WARNERS);
    if (_actualSoftClauses->size() > 1500) {
      compression = 1;
    } else {
      compression = 99;
    }
    if (_settings.verbosity > 0)
      std::cout << "c compressionRate set....: " << compression << std::endl;
  }
  //  std::cout << std::endl;
  uint64_t answer;
  // in BAILLEUXW2 a table of length maxWeight has to be initialized
  // therefore set a n upper bound on _sumOfActualWeights
  if (_encoding == BAILLEUXW2 && _sumOfActualSoftWeights > 1073741824) {
    _settings.SetEncoding(WARNERS);
    _encoding = WARNERS;
  }
  // Set reasonable upper bound for those encodings on _sumOfActualWeights
  if ((_encoding == BAILLEUX || _encoding == ASIN || _encoding == OGAWA) && _sumOfActualSoftWeights > 2048) {
    _settings.SetEncoding(WARNERS);
    _encoding = WARNERS;
  }
  if (_encoding == MRWTO || _encoding == MRWTO2 ||
            _encoding == MRWTO19 || _encoding == MRWTO19_2) {
    answer = _sumOfActualSoftWeights + 1;
  } else {
    answer = _sumOfActualSoftWeights;
  }
  uint64_t oldanswer = _sumOfActualSoftWeights;
  if (_settings.verbosity > 0) {
    std::cout << "c _sumOfActualSoftWeights: " << _sumOfActualSoftWeights
              << std::endl;
    std::cout << "c Local No Soft Clauses..: " << _actualSoftClauses->size()
              << std::endl;
  }

  // this vector is used to create clauses!
  std::vector<uint32_t> lits;
  // loop count
  int lcnt = 0;
  std::vector<uint32_t> linkingVar;
  std::vector<long long int> linkingWeight; // uemura 20161202
  //    bool *mmodel = new bool[static_cast<uint64_t>(_nbOfOrigVars)];
  //    //uemura 20161128
  lits.clear();
  linkingVar.clear();
  linkingWeight.clear();

  long long int divisor = 1; // koshi 2013.10.04

  std::vector<long long int>
      ndivisor; // mrwto用の複数の基数を保存する変数 uemura 2016.12.05
  ndivisor.clear();

  std::vector<std::vector<uint32_t>>
      linkingVarMR; // uemura 2016.12.05 for mrwto
  linkingVarMR.clear();
  std::vector<std::vector<long long int>>
      linkingWeightMR; // uemura 2016.12.05 for mrwto
  linkingWeightMR.clear();

  std::vector<long long int> cc; // cardinality constraints
  cc.clear();
  size_t ccSizeOld = SIZE_MAX;

  // koshi 20140701        lbool ret = S.solveLimited(dummy);
  uint32_t ret = UNKNOW;
  _encodings->_relaxLit = 0;
  uint64_t answerNew = _localUnSatWeight;

  // if (answerNew != 0) {
  assert(answerNew != 0);
  uint32_t ncls = _satSolver->GetNumberOfClauses();

  //      printf("c linkingVar.size() = %zu\n", linkingVar.size());
  // uemura 20161128
  genCardinals(answerNew, divisor, lits, linkingVar, linkingWeight,
                ndivisor, linkingVarMR, linkingWeightMR, compression);

  //      printf("c linkingVar.size() = %zu\n", linkingVar.size());
  //      std::cout << "Clauses before: " << ncls << std::endl;
  //      std::cout << "Clauses after: " <<
  //      _satSolver->GetNumberOfClauses()
  //                << std::endl;


  if (_settings.verbosity > 0) {
    std::cout << "c Cardinals are generated!" << std::endl;
    std::cout << "c Clauses of encoding: " <<
      _satSolver->GetNumberOfClauses()
      - ncls
                << std::endl;
    // uemura 20161129
    if (_settings.verbosity > 1) {
    if (_encoding == WMTO || _encoding == MRWTO || _encoding == MRWTO2 ||
        _encoding == MRWTO19 || _encoding == MRWTO19_2) {
      std::cout << "linkingVarMR.size(): " << linkingVarMR.size()
                << std::endl;
      printf("c ");
      for (size_t i = 0; i < linkingVarMR.size(); i++) {
        printf("c linkingVar[%lu].size = %zu, ", i, linkingVarMR[i].size());
      }
      printf("\n");
    } else {
      printf("c linkingVar.size() = %zu\n", linkingVar.size());
    }
    printf("c Cardinality Constraints: %d variables and %d clauses\n",
            _satSolver->GetNumberOfVariables() - _nbOfOrigVars,
            _satSolver->GetNumberOfClauses() - ncls);
    }
  }
  _encoding = _settings.GetEncoding();

  while ((ret = _satSolver->Solve()) == SAT) { // koshi 20140107
    CalculateSATWeight();
    answerNew = _localUnSatWeight;
    //    std::cout << "c noOfClauses AfterSolve: " <<
    //    _satSolver->GetNumberOfClauses()
    //              << std::endl;
    if (_encodings->_relaxLit != 0) {
      _satSolver->ClearAssumption();
      _satSolver->ResetClause();
      _satSolver->NewClause();
      uint32_t rl = _encodings->_relaxLit ^ 1;
      _satSolver->AddLiteral(&rl);
      _satSolver->CommitClause();
      //      std::cout << "c RELAXLIT added: " << (_encodings->_relaxLit ^ 1)
      //                << std::endl;
      _encodings->_relaxLit = 0;
    }
    assert(std::cout << "c assertion Solver call in Pacose, SolveQMax0" << std::endl && _satSolver->Solve() == 10);
    //        std::cout << "ret: " << ret << std::endl;
    lcnt++;

    if (_settings.verbosity > 10) {
      std::cout << std::setw(30) << "blockings size: " << _blockings.size()
                << std::endl;
      std::cout << std::setw(30) << "weights size: " << _weights.size()
                << std::endl;
      std::cout << std::setw(30) << "lcnt: " << lcnt << std::endl;
      std::cout << std::setw(30) << "answer: " << answer << std::endl;
      std::cout << std::setw(30) << "answerNew: " << answerNew << std::endl;
    }

    if (lcnt > 2) {
      assert(answerNew < answer);
    }

    //    std::cout << "answer + answernew: " << answer << "  " << answerNew
    //              << std::endl;
    if (answerNew > 0) {
      uint32_t nofCl = _satSolver->GetNumberOfClauses();
      if (_encoding == WMTO || _encoding == MRWTO || _encoding == MRWTO2 ||
          _encoding == MRWTO19 || _encoding == MRWTO19_2) {
        _encodings->lessthanMR(linkingVarMR, linkingWeightMR, answer, answerNew,
                               ndivisor, cc, *_satSolver, lits, _encoding);
      } else {
        if (_encoding == BAILLEUX && lcnt == 1)
          answer = linkingVar.size();

        //        std::cout << "answer + answernew: " << answer << "  " <<
        //        answerNew
        //                  << std::endl;
        ccSizeOld = cc.size();
        //        std::cout << "ccSizeOld: " << cc.size() << std::endl;

        _encodings->lessthan(linkingVar, linkingWeight, answer, answerNew,
                             divisor, cc, *_satSolver, _encoding);
      }
      if (_settings.verbosity > 0)
        std::cout << "c Additional Encoding Clauses: "
                  << _satSolver->GetNumberOfClauses() - nofCl << std::endl;
      //      if (_satSolver->GetNumberOfClauses() - nofCl == 0) {
      //        exit(1);
      //      }
      //      std::cout << "ccSizeNew: " << cc.size() << std::endl;
      oldanswer = answer;
      answer = answerNew;
    } else {
      answer = answerNew;
      ret = UNSAT; // koshi 20140124
      break;
    }
    //    std::cout << "noOfClauses BeforeNextSolve: "
    //              << _satSolver->GetNumberOfClauses() << std::endl;
    //    std::cout << "oldAnswer: " << oldanswer << "  answer: " << answer
    //              << std::endl;
  }

  //  std::cout << ret << std::endl;
  //  std::cout << "answer + anwernew: " << answer << std::endl;

  // koshi 20140124
  if (ret == UNSAT) {
    if (lcnt > 0) {
      if (_settings.verbosity > 0)
        printf("c local opt found\n");

      if (answer != 0 && oldanswer == answer + 1 && lcnt > 1) {
        assert(std::cout << "c assertion Solver call in Pacose, SolveQMax1" << std::endl && _satSolver->Solve() == 20);
        _satSolver->ClearAssumption();
        assert(std::cout << "c assertion Solver call in Pacose, SolveQMax2" << std::endl && _satSolver->Solve() == 10);
        //        uint32_t lastResult = _satSolver->Solve();
        //        std::cout << "SolveBefore: " << lastResult <<
        //        std::endl;
      } else if (answer != 0) {
        assert(std::cout << "c assertion Solver call in Pacose, SolveQMax3" << std::endl && _satSolver->Solve() == 20);
        _satSolver->ClearAssumption();
        assert(std::cout << "c assertion Solver call in Pacose, SolveQMax4" << std::endl && _satSolver->Solve() == 10);
        //        std::cout << "SOLVE RESULT: " << _satSolver->Solve() <<
        //        std::endl;
        // deactivate last assumption!
        //        _satSolver->ClearAssumption();
        //        uint32_t lastResult = _satSolver->Solve();
        //        std::cout << "SolveBefore: " << lastResult << std::endl;
        //        assert(lastResult == SAT);
        if (_encodings->_relaxLit != 0) {
          _satSolver->ResetClause();
          _satSolver->NewClause();
          uint32_t rl = _encodings->_relaxLit;
          //        std::cout << "AddRL: " << _encodings->_relaxLit <<
          //        std::endl;
          _satSolver->AddLiteral(&rl);
          _satSolver->CommitClause();
        }
        _encodings->_relaxLit = 0;
        //        std::cout << "SolveBefore2: " << _satSolver->Solve() <<
        //        std::endl; printf("local o %lld\n",
        //               _sumOfActualSoftWeights - CalculateLocalSATWeight());
        //        std::cout << "noOfClauses: " <<
        //        _satSolver->GetNumberOfClauses()
        //                  << std::endl;
        if (ccSizeOld == SIZE_MAX || cc.size() > ccSizeOld) {
          if (_settings.verbosity > 0)
            std::cout << "c ccSize is going to be decreased Old/New: "
                      << ccSizeOld << ", " << cc.size() << std::endl;
          while (cc.size() > ccSizeOld) {
            cc.pop_back();
          }
        }
        //        std::cout << "ccSizeOld: " << cc.size() << std::endl;
        _satSolver->ClearAssumption();
        // uint32_t nofCl = _satSolver->GetNumberOfClauses();
        if (_encoding == WMTO || _encoding == MRWTO || _encoding == MRWTO2 ||
            _encoding == MRWTO19 || _encoding == MRWTO19_2) {
          _encodings->lessthanMR(linkingVarMR, linkingWeightMR, oldanswer,
                                 answer + 1, ndivisor, cc, *_satSolver, lits,
                                 _encoding);
        } else {
          //          if (_encoding == BAILLEUX && lcnt == 1) answer =
          //          linkingVar.size();

          _encodings->lessthan(linkingVar, linkingWeight, oldanswer, answer + 1,
                               divisor, cc, *_satSolver, _encoding);
        }
        //        std::cout << "Newly Added clauses: "
        //                  << _satSolver->GetNumberOfClauses() - nofCl <<
        //                  std::endl;
        //        //        std::cout << "ccSizeNew: " << cc.size() <<
        //        std::endl; uint32_t lastResult = _satSolver->Solve();

        //        std::cout << "SolveAfterWithNewAssumptions: " << lastResult
        //                  << std::endl;
        assert(std::cout << "c assertion Solver call in Pacose, SolveQMax5" << std::endl && _satSolver->Solve() == SAT);
        if (_encodings->_relaxLit != 0) {
          _satSolver->ResetClause();
          _satSolver->NewClause();
          uint32_t rl = _encodings->_relaxLit ^ 1;

          _satSolver->AddLiteral(&rl);
          //        std::cout << "AddRL: " << _encodings->_relaxLit <<
          //        std::endl;
          _satSolver->CommitClause();
        }
        _encodings->_relaxLit = 0;
        //        CalculateLocalSATWeight();
        _satSolver->ClearAssumption();
        //        lastResult = _satSolver->Solve();
        //        std::cout << "SolveAfter2: " << lastResult << std::endl;
        assert(std::cout << "c assertion Solver call in Pacose, SolveQMax6" << std::endl && _satSolver->Solve() == SAT);
        //        CalculateLocalSATWeight();
        //        CalculateSATWeight();
      } else {
        if (_settings.verbosity > 0)
          std::cout << "c ANSWER IS 0!, Add all SCs" << std::endl;
        
        wbSortAndFilter();
        //        uint32_t lastResult = _satSolver->Solve();
        //        std::cout << "SolveAfter2: " << lastResult << std::endl;
        //        assert(lastResult == SAT);

        //        case answer == 0 -- all SCs are SAT
        // for (auto SC : *_actualSoftClauses) {
        //   _satSolver->ResetClause();
        //   _satSolver->NewClause();
        //   uint32_t rlit = SC->relaxationLit ^ 1;
        //   _satSolver->AddLiteral(&rlit);
        //   _satSolver->CommitClause();
        // }
        // _satSolver->ClearAssumption();
        //        lastResult = _satSolver->Solve();
        //        std::cout << "SolveAfterAddingAllRelaxLits: " << lastResult
        //                  << std::endl;
        //        assert(lastResult == SAT);
      }
    } else {
      if (_settings.verbosity > 0)
        printf("s Hard clauses are UNSATISFIABLE\n");
    }
  } else if (ret == SAT) {
    std::cout << "c ERROR: SHOULD NEVER OCCUR!" << std::endl;
  } else if (_settings.verbosity > 0) {
    printf("s UNKNOWN\n");
    printf("c Search is stopped by a limit (maybe time-limit)\n");
  }

  if (lcnt > 0 && _settings.verbosity > 0)
    std::cout << "c Latest Answer = " << answer << "/" << _sumOfActualSoftWeights
              << "  loops " << lcnt << std::endl;

  if (_settings.verbosity > 5)
    std::cout << "END QMAXSAT" << std::endl << std::endl;

  //    delete[] mmodel;
  return ret;
}

bool Pacose::TreatBorderCases() {

  wbSortAndFilter();

  if (_actualSoftClauses->size() == 0) {
    CalculateSATWeight();
    // We cannot conclue overall optimality 
    // only for the _actualSoftClauses

    // Test this corner case 
    // 1. Eithter TrimMaxSAT removed all SoftClauses OR
    // 2. WBSortAndFilter removed already all SoftClauses
    return true;
  } else if (_localUnSatWeight == 0) {
    // case all SCs are already SAT
    // PROBABLY NEVER TO HAPPEN!!
    // Not reproducible with short fuzzer test!
    // as TrimMaxSAT or WBSortAndFilter already
    // removed all SCs
    std::cout << "c SHOULDNEVERHAPPEN!" << std::endl;

    for (auto SC : *_actualSoftClauses) {
      uint32_t relaxLit = SC->relaxationLit ^ 1;
      _satSolver->ResetClause();
      _satSolver->NewClause();
      _satSolver->AddLiteral(&relaxLit);
      _satSolver->CommitClause();
      _satSolver->ClearAssumption();
    }
    assert(std::cout << "c assertion Solver call in Pacose, TreatBorderCases1" << std::endl && _satSolver->Solve() == SAT);
    return true;
  } else if (_actualSoftClauses->size() == 1) {
    std::cout << "c Border Case ONLY ONE SOFT CLAUSE" << std::endl;
    // border case - only one SC
    // which is not yet SAT otherwise localUnsatWeight would be 0
    // TODO Dieter: Test this corner case
    _satSolver->ClearAssumption();
    uint32_t relaxLit = (*_actualSoftClauses)[0]->relaxationLit ^ 1;
    uint32_t negRelaxLit = (*_actualSoftClauses)[0]->relaxationLit;
    uint64_t wghtRelaxLit = (*_actualSoftClauses)[0]->originalWeight;
    _satSolver->AddAssumption(&relaxLit);
    _satSolver->ResetClause();
    _satSolver->NewClause();
    uint32_t rv = _satSolver->Solve();
    assert(rv == SAT or rv == UNSAT);
    
    if (rv == SAT) {
      _satSolver->AddLiteral(&relaxLit);
      _satSolver->CommitClause();
      _satSolver->ClearAssumption();
      
    } else if (rv == UNSAT) {

      _satSolver->AddLiteral(&negRelaxLit);
      _satSolver->CommitClause();
      _satSolver->ClearAssumption();
      
      // TODO TOBI: IS THIS SOLVER CALL NECESSARY?
      rv = _satSolver->Solve();
      assert(rv == SAT);
    }
    CalculateSATWeight();
    return true;
  }
  return false;
}

void Pacose::ChooseEncoding() {
  if (_encoding == HEURISTIC20) {
    std::cout << "c Use HEURISTIC Pacose 2020" << std::endl;
    if (_settings.verbosity > 1)
      std::cout << "c HEURISTIC20" << std::endl;
    if (((_overallSoftWeights / _originalSoftClauses.size()) > 850) &&
        (_overallSoftWeights < 80000000000) &&
        (_originalSoftClauses.size() < 50000)) {
      _settings.SetEncoding(WARNERS);
      _encoding = WARNERS;
      std::cout << "c Use Warners encoding!" << std::endl;
    } else {
      _settings.SetEncoding(DGPW18);
      _encoding = DGPW18;
      std::cout << "c Use DGPW encoding!" << std::endl;
    }
  } else if (_encoding == HEURISTIC1819) {
    //    (answer > 20000000000 || answer < 400000))
    if (_settings.verbosity > 0)
      std::cout << "c Use HEURISTIC Pacose 2018/2019" << std::endl;
    //    std::cout << "_overallSoftWeights: " << _overallSoftWeights <<
    //    std::endl;

    if (_overallSoftWeights > 20000000000 || _overallSoftWeights < 400000) {
      _settings.SetEncoding(DGPW18);
      _encoding = DGPW18;
      if (_settings.verbosity > 0)
        std::cout << "c Use DGPW encoding!" << std::endl;
    } else {
      _settings.SetEncoding(WARNERS);
      _encoding = WARNERS;
      _settings.SetCompression(1);
      if (_settings.verbosity > 0)
        std::cout << "c Use Warners encoding!" << std::endl;
    }
  } else if (_encoding == HEURISTICQMAXSAT17) {
    //    (answer > 20000000000 || answer < 400000))
    if (_settings.verbosity > 0)
      std::cout << "c Use HEURISTIC QMaxSAT 2017" << std::endl;

    int k = _originalSoftClauses.size();
    int logk = ceil(log2(_originalSoftClauses.size()));
    int logsum = ceil(log2(_sumOfSoftWeights));
    // calculate 2nd logarithm of k
    //    for (long long int ok = k; ok > 0; ok = ok >> 1) logk++;
    //    // calculate 2nd logarithm of sum of weights
    //    for (long long int osum = _sumOfSoftWeights; osum > 0; osum = osum >>
    //    1)
    //      logsum++;

    if (_settings.verbosity > 1)
      printf("c logk = %d, logsum = %d\n", logk, logsum);

    //    std::cout << "logk: " << ceil(log2(_originalSoftClauses.size()))
    //              << "  logsum: " << ceil(log2(_sumOfSoftWeights)) <<
    //              std::endl;
    if (logk + logsum < 15) {
      // Bailleux
      _encoding = BAILLEUX;
      _settings.SetEncoding(BAILLEUX);
      _settings.SetCompression(0);
      if (_settings.verbosity > 0)
        printf("c Bailleux's encoding (comp=0)\n");
    } else if (k < 3) {
      // Warners
      _encoding = WARNERS;
      _settings.SetEncoding(WARNERS);
      _settings.SetCompression(1);
      if (_settings.verbosity > 0)
        printf("c Warners' encoding (comp=1)\n");
    } else if (logsum < 17) {
      // Ogawa
      _encoding = OGAWA;
      _settings.SetEncoding(OGAWA);
      _settings.SetCompression(0);
      if (_settings.verbosity > 0)
        printf("c Ogawa's encoding (comp=0)\n");
    } else {
      _encoding = WARNERS;
      _settings.SetEncoding(WARNERS);
      _settings.SetCompression(1);
      if (_settings.verbosity > 0)
        printf("c Warners' encoding (comp=1).\n");
    }
  } else if (_encoding == QMAXSAT19) {
    _encoding = MRWTO19;
    _settings.SetEncoding(MRWTO19);
    _settings.SetCompression(0);
  }

  std::cout << "c encoding...............: "
            << _settings.ReturnEncodingString(_encoding) << std::endl;
}

// void Pacose::CallMaxPre2(ClauseDB &clauseDB) {
//   if (!_settings.useMaxPre2)
//     return;
//   std::cout << "c Use MaxPre2............: " << _settings.useMaxPre2
//             << std::endl;
//   if (_settings.verbosity > 1) {
//     uint32_t nbSoftClbefore =
//         std::count_if(clauseDB.weights.begin(), clauseDB.weights.end(),
//                       [&clauseDB](uint64_t weight) {
//                         return weight <= clauseDB.sumOfSoftWeights;
//                       });

//     std::cout << "c SoftClauses before: " << nbSoftClbefore << std::endl;
//     std::cout << "c Initialize maxpre2" << std::endl;
//   }
//   maxpre = new maxPreprocessor::PreprocessorInterface(
//       clauseDB.clauses, clauseDB.weights, clauseDB.sumOfSoftWeights + 1);

//   maxpre->preprocess(_settings.maxPre2Techniques, _settings.verbosity,
//                      _settings.maxPre2TimeOut);

//   std::vector<int> retLabels;
//   clauseDB.clauses.clear();
//   clauseDB.weights.clear();
//   maxpre->getInstance(clauseDB.clauses, clauseDB.weights, retLabels);
//   _unSatWeight = maxpre->getUpperBound();

//   if (_settings.verbosity > 1) {
//     std::cout << "c removedWeight: " << maxpre->getRemovedWeight() << std::endl;
//     std::cout << "c new Top Weight: " << maxpre->getTopWeight() << std::endl;
//     std::cout << "c new UpperBound: " << maxpre->getUpperBound() << std::endl;
//   }

//   // maxpre.printInstance(std::cout);
//   // std::cout << "MAP:" << std::endl;
//   // maxpre.printMap(std::cout);
//   // std::cout << "TechniqueLog:" << std::endl;
//   // maxpre.printTechniqueLog(std::cout);
//   // std::cout << "TimeLog:" << std::endl;
//   // maxpre.printTimeLog(std::cout);
//   // std::cout << "InfoLog:" << std::endl;
//   // maxpre.printInfoLog(std::cout);
//   // std::cout << "PreproStats:" << std::endl;
//   // maxpre.printPreprocessorStats(std::cout);
//   // for (size_t i = 0; i < clauseDB.clauses.size(); ++i) {
//   //   if (clauseDB.weights[i] > clauseDB.sumOfSoftWeights)
//   //     std::cout << "h ";
//   //   else
//   //     std::cout << clauseDB.weights[i] << " ";
//   //   for (auto lit : clauseDB.clauses[i]) {
//   //     std::cout << lit << " ";
//   //   }
//   //   // std::cout << "0 " << retLabels[i] << std::endl;
//   //   std::cout << "0 " << std::endl;
//   // }
// }

bool Pacose::ExternalPreprocessing(ClauseDB &clauseDB) {

  // CallMaxPre2(clauseDB);
  // count soft clauses after
  // uint32_t nbSoftCl =
  //     std::count_if(clauseDB.weights.begin(), clauseDB.weights.end(),
  //                   [&clauseDB](uint64_t weight) {
  //                     return weight <= clauseDB.sumOfSoftWeights;
  //                   });
  // uint64_t sumOfWeightsAfter = std::accumulate(
  //     clauseDB.weights.begin(), clauseDB.weights.end(), 0ull,
  //     [&clauseDB](uint64_t sum, uint64_t weight) {
  //       return sum + (weight <= clauseDB.sumOfSoftWeights ? weight : 0);
  //     });
  uint32_t nbSoftCl = clauseDB.nbSoftClauses;
  uint64_t sumOfWeightsAfter = clauseDB.sumOfSoftWeights;

  if (_settings.verbosity > 1) {
    std::cout << "c clauses size: " << clauseDB.clauses.size() << std::endl;
    std::cout << "c sumOfWeightsBefore: " << clauseDB.sumOfSoftWeights
              << std::endl;
    std::cout << "c SoftClauses: " << nbSoftCl << std::endl;
    std::cout << "c sum of Weights: " << clauseDB.sumOfSoftWeights << "/"
              << sumOfWeightsAfter << std::endl;
  }

  _nbClauses =  clauseDB.clauses.size();
  _nbOfOrigVars = 0;
  for (auto clause : clauseDB.clauses) {
    for (auto lit : clause) {
      if (abs(lit) > _nbOfOrigVars) {
        _nbOfOrigVars = abs(lit);
      }
    }
  }
  _nbOfOrigPlusSCRelaxVars = _nbVars = clauseDB.nbVars = _nbOfOrigVars;
  
  // uint64_t maxSumOfWeights = 9223372036854775808;
  uint64_t maxSumOfWeights = 4611686018427387904ULL; // 2^62
  if (sumOfWeightsAfter >= maxSumOfWeights) {
    std::cout << "c ATTENTION: weights bigger than 2^62 are currently not "
                 "supported. This might produce wrong results!"
              << std::endl;
  }
  std::cout << "c sum of weights.........: " << sumOfWeightsAfter << std::endl;
  _maxWeight = clauseDB.maxWeight;

  _satSolver->NewVariables(_nbVars + 1);
  // std::cout << "_nbVars = " << _nbVars << std::endl;

  assert(clauseDB.clauses.size() == clauseDB.weights.size());
  uint64_t emptyWeight = 0; 

  for (size_t i = 0; i < clauseDB.clauses.size(); ++i) {
    // hard clause
    // for (auto clause : clauseDB.clauses) {
    if (clauseDB.weights[i] > clauseDB.sumOfSoftWeights) {
      if (clauseDB.clauses[i].empty()) {
        // Border Case, empty hard clause cannot be satisfied, thus the
        // instance is UNSATISFIABLE!
        std::cout << "s UNSATISFIABLE" << std::endl;
        return false;
      }
      _hasHardClauses = true;
      std::vector<uint32_t> *clause = new std::vector<uint32_t>;
      for (auto lit : clauseDB.clauses[i]) {
        clause->push_back(clauseDB.SignedTouint32_tLit(lit));
      }
      // _satSolver->AddClause(clauseDB.clauses[i]);
      AddClause(*clause);
    } else {
      // soft clause
      if ((clauseDB.clauses[i].empty())) {
        emptyWeight += clauseDB.weights[i];
        continue;
      }
      std::vector<uint32_t> *sclause = new std::vector<uint32_t>;
      for (auto lit : clauseDB.clauses[i]) {
        sclause->push_back(clauseDB.SignedTouint32_tLit(lit));
      }
      AddSoftClause(*sclause, clauseDB.weights[i]);
    }
    _nbOfOrigPlusSCRelaxVars = _satSolver->GetNumberOfVariables();
  }


  _sumOfSoftWeights = clauseDB.sumOfSoftWeights = sumOfWeightsAfter;
  // // std::cout << "HardClauses: " << clauseDB.clauses.size() << std::endl;
  // uint64_t emptyWeight = 0;
  // for (auto sclause : clauseDB.sclauses) {
  //   if (sclause.first->empty()) {
  //     emptyWeight += sclause.second;
  //     continue;
  //   }
  //   AddSoftClause(*sclause.first, sclause.second);
  // }

  // std::cout << "SoftClauses: " << clauseDB.sclauses.size() << std::endl;

  // CASE NO SOFT CLAUSES:
  if (_originalSoftClauses.empty()) {
    _unSatWeight = emptyWeight;
    if (!_hasHardClauses) {
      std::cout << "s OPTIMUM FOUND" << std::endl;
      std::cout << "o " << emptyWeight << std::endl;
      PrintResult();
      // TODO-Dieter: Check the way empty soft clauses are treated. They should be rewritten by objective update rule. 
      // std::cout << "v " << std::endl;
      return false;
    }
    uint32_t rv = _satSolver->Solve();
    if (rv == 10) {
      std::cout << "s OPTIMUM FOUND" << std::endl;
      std::cout << "o " << emptyWeight << std::endl;
      PrintResult();
    } else if (rv == 20) {
      std::cout << "s UNSATISFIABLE" << std::endl;
    } else {
      std::cout << "s UNKNOWN" << std::endl;
    }
    return false;
  }

  if (emptyWeight > 0) {
    // Border Case, empty soft clauses were added which are unsatisfiable.
    // Trick the solver and add one unsatisfiable soft clause.
    // TODO Dieter: Only one clause for all empty soft clauses. 

    uint32_t lit = (_satSolver->NewVariable() << 1);
    std::vector<uint32_t> *sclause = new std::vector<uint32_t>;
    sclause->push_back(lit);
    
    AddClause(*sclause);
    (*sclause)[0] = lit ^ 1;
    AddSoftClause(*sclause, emptyWeight);
    _nbOfOrigPlusSCRelaxVars = _satSolver->GetNumberOfVariables();
  }

  clauseDB.clauses.clear();
  clauseDB.weights.clear();

  return true;
}


uint32_t Pacose::SolveProcedure(ClauseDB &clauseDB) {

  if (!ExternalPreprocessing(clauseDB)) {
    return 0;
  };

  _settings.formulaIsDivided = true;
  double timeStart;
  struct rusage resources;
  getrusage(RUSAGE_SELF, &resources);
  timeStart =
      resources.ru_utime.tv_sec + 1.e-6 * (double)resources.ru_utime.tv_usec;
  std::vector<uint32_t> lastSatisfiableAssignment = {};

  // tests if the instances are possible to split with GBMO
  DivideSCsIfPossible();
  getrusage(RUSAGE_SELF, &resources);
  _GBMOTime =
      (resources.ru_utime.tv_sec + 1.e-6 * (double)resources.ru_utime.tv_usec) -
      timeStart;
  _GBMOPartitions = _sClauses.size();
  _originalSClauses = {};
  for ( auto clauses : _sClauses ) {
    std::vector<SoftClause *> softClauses;
    for (auto clause : clauses) {
      softClauses.push_back(clause);
    }
    _originalSClauses.push_back(softClauses);
  }
  

  _settings.Print();

  if (_satSolver->Solve() == 20) {
    std::cout << "c Hard Clauses are not Satisfiable!" << std::endl;
    std::cout << "s UNSATISFIABLE" << std::endl;
    return 20;
  }

  CalculateSATWeight();

  _encoding = _settings._encoding;
  ChooseEncoding();
  
  std::cout << "c #SoftClauses...........: " << _originalSoftClauses.size()
            << std::endl;
  std::cout << "c #HardClauses...........: " << _satSolver->GetNumberOfClauses()
            << std::endl;
  std::cout << "c #Variables.............: "
            << _satSolver->GetNumberOfVariables() << std::endl;

  // std::cout << "Sclauses.size() " << _sClauses.size() << std::endl;

  for (uint32_t i = static_cast<uint32_t>(_sClauses.size()); i > 0; i--) {
    if (_settings.verbosity > 0)
      std::cout << std::endl << "--- NEW GBMO LEVEL " << i << "/" << _sClauses.size() << " -------------- soft clauses in this level: " << _sClauses[i - 1].size() << std::endl;
    
    // GBMO starts
    _localUnSatWeight = UINT64_MAX;
    _localSatWeight = 0;
    _settings.currentCascade.iteration = i - 1;
    _actualSoftClauses = &_sClauses[i - 1];
    _originalActualSoftClauses = &_originalSClauses[i-1];
    _cascCandidates[i - 1].dgpw = nullptr;

    // calc greatest common divisor
    // convert into unweighted MaxSAT if possible!
    Preprocess();
    if (i != _sClauses.size() && _satSolver->Solve() == 20) {
      std::cout << "c At the beginning of a new level, the clauses are not Satisfiable!" << std::endl;
      assert(false);
    }

    CalculateSATWeight();
    // If this is here it throws sometimes errors!

    if (_localSatWeight == _sumOfActualSoftWeights) {
      if (_actualSoftClauses->size() > 0)
        wbSortAndFilter();
      // std::cout << "localSatWeight == sumOfActualWeights: " << localSatWeight  << " == " << sumOfActualWeights << std::endl;
      // std::cout << "_satWeight ?? _overallSoftWeights: " << _satWeight  << " ?? " << _overallSoftWeights << std::endl;

      if (_satWeight == _overallSoftWeights) {
        wbSortAndFilter();
        break;
      }
      continue;
    }

    if (_settings.verbosity > 0) {
      std::cout << "c local o: " << _localUnSatWeight << std::endl;
      if (_settings.verbosity > 3) {
        std::cout << "c number of SCs: " << _actualSoftClauses->size()
                  << std::endl;
        std::cout << "c sumOfActualWeights: " << _sumOfActualSoftWeights
                  << std::endl;
        std::cout << "c locSatWeight: " << _localSatWeight << std::endl;
      }
    }

    int fixSCs = _settings.greedyPPFixSCs;
    size_t minSizeSCs = _settings.greedyMinSizeOfSet;

    // heuristic to choose when to use TrimSAT
    if (fixSCs == -1) {
      if (_encoding == DGPW18) {
        minSizeSCs = 650;
        fixSCs = 2;
      } else if (_encoding == WARNERS) {
        minSizeSCs = 12000;
        fixSCs = 2;
      } else {
        minSizeSCs = 12000;
        fixSCs = 2;
        //        std::cout << "ERROR: no heuristic set yet!" << std::endl;
        //        assert(true);
        //        exit(1);
      }
      if (_settings.verbosity > 0)
        std::cout << "c MinSize for PP chose...: " << minSizeSCs << std::endl;
    }

    //    std::cout << "sumOfActualWeights != _satWeight: "
    //              << (sumOfActualWeights != _satWeight) << std::endl;
    //    std::cout << "_actualSoftClauses->size(): " <<
    //    _actualSoftClauses->size()
    //              << std::endl;
    
    assert(std::cout << "c assertion Solver call in Pacose, SolveProcedure1" << std::endl && _satSolver->Solve() == SAT);

    // TRIMMaxSAT
    if (_actualSoftClauses->size() != 0 && _sumOfActualSoftWeights != _satWeight &&
        ((_settings.greedyPrepro != 0 && _settings.greedyPPFixSCs != -1) ||
         (_settings.greedyPrepro != 0 &&
          _actualSoftClauses->size() > minSizeSCs))) {
      _noTrimSAT++;

      double tmpTimeTrimming;
      struct rusage resources;
      getrusage(RUSAGE_SELF, &resources);
      tmpTimeTrimming = resources.ru_utime.tv_sec +
                        1.e-6 * (double)resources.ru_utime.tv_usec;
      uint32_t tmpNoSolverCalls = _satSolver->_noSolverCalls;

      GreedyPrepro greedyPrePro = GreedyPrepro(*_actualSoftClauses, &_settings,
                                               _satSolver, this, fixSCs);

      greedyPrePro.StartPrepro();
      // uint64_t gPPUnSatWeight = greedyPrePro.StartPrepro();
      // if (gPPUnSatWeight < _localUnSatWeight) {
      //   _localUnSatWeight = gPPUnSatWeight;
      //   _localSatWeight = greedyPrePro._satWeight;
      // }

      _noTrimSATSolverCalls += (_satSolver->_noSolverCalls - tmpNoSolverCalls);
      _alwaysUNSATSCs += greedyPrePro.GetAlwaysUNSATSCs();
      _alwaysSATSCs += greedyPrePro.GetAlwaysSATSCs();
      _alwaysUNSATWeight += greedyPrePro.GetAlwaysUNSATWeight();
      _alwaysSATWeight += greedyPrePro.GetAlwaysSATWeight();

      uint32_t rv = _satSolver->Solve();
      assert(rv == SAT);
      CalculateSATWeight();

      getrusage(RUSAGE_SELF, &resources);
      tmpTimeTrimming = resources.ru_utime.tv_sec +
                        1.e-6 * (double)resources.ru_utime.tv_usec -
                        tmpTimeTrimming;
      _trimSATTime += tmpTimeTrimming;


      if (_settings.createSimplifiedWCNF) {
        // Write out the Simplified WCNF and EXIT
        DumpSolvingInformation();

        greedyPrePro.CreateSimplifiedWCNF(_settings.maxCnfFile,
                                          _satSolver->CNF);
      }

      if (_settings.verbosity > 2)
        std::cout << "c local PrePro o: " << _localUnSatWeight << std::endl;
    }

    if (TreatBorderCases()) {
      // it is a border case, then continue with next GBMO MaxSAT problem
      continue;
    }

    // calc greatest common divisor
    // convert into unweighted MaxSAT if possible!
    Preprocess();

    int64_t optimum = -1;
    if (_encoding == DGPW18 || _actualSoftClauses->size() <= 2) {
      _cascCandidates[i - 1].dgpw = new DGPW::DGPW(this);
      _cascCandidates[i - 1].dgpw->SetSoftClauseVector(_actualSoftClauses);
      if (_settings.verbosity > 0)
        std::cout << "c greatest Common Divisor: " << _GCD << std::endl;
      _cascCandidates[i - 1].dgpw->SetGreatestCommonDivisor(_GCD);
      _settings.currentCascade._solveTares = true;
      if (_settings.onlyWithAssumptions)
        _settings.currentCascade._onlyWithAssumptions = true;
      else
        _settings.currentCascade._onlyWithAssumptions = false;
      CalculateLocalSATWeight();
      _cascCandidates[i - 1].dgpw->SetSatWeight(_localSatWeight);

      // DGPW solve procedure
      _cascCandidates[i - 1].dgpw->MaxSolveWeightedPartial(optimum);

      _clausesOfEncoding += _cascCandidates[i - 1].dgpw->GetEncodingClauses();
      _variablesOfEncoding +=
          _cascCandidates[i - 1].dgpw->GetEncodingVariables();
      
    } else {
      // QMaxSAT Encodings WARNERS ENCODING
      _satSolver->ClearAssumption();
      uint32_t tmpNoClauses = _satSolver->GetNumberOfClauses();
      uint32_t tmpNoVariables = _satSolver->GetNumberOfVariables();
      SolveQMax(&_encoding);
      _clausesOfEncoding += _satSolver->GetNumberOfClauses() - tmpNoClauses;
      _variablesOfEncoding +=
          _satSolver->GetNumberOfVariables() - tmpNoVariables;
    }

    _satSolver->ClearAssumption();

    if (_lastCalculatedUnsatWeight != _unSatWeight) {
      uint32_t cr = _satSolver->Solve();
      if (cr != 10) {
        std::cerr << "c ERROR: STRANGE SOLVING RESULT " << cr << " -> QUIT!"
                  << std::endl;
        exit(0);
      }
      CalculateSATWeight();
    }
  }

  // if (_alwaysUNSATWeight == _overallSoftWeights) {
  //   _unSatWeight = _overallSoftWeights;
  // } else if (_alwaysSATWeight == _overallSoftWeights) {
  //   _unSatWeight = 0;
  // }

  DumpSolvingInformation();

  uint32_t solutionCount = 1;
  if (_settings.calculateAllSolutions) {
    PrintResult();
    while (CalculateNextResult() == SAT) {
      solutionCount++;
      std::cout << "c " << solutionCount << " result:" << std::endl;
      PrintResult();
    }
  } else if (_settings.calculateAllSoftClauseCombinations) {
    PrintResult();
    while (CalculateNextSoftclauseCombination() == SAT) {
      solutionCount++;
      std::cout << std::endl
                << "c " << solutionCount << " result:" << std::endl;
      PrintResult();
    }
    std::cout << "c NoOfSolutionsFound.....: " << solutionCount << std::endl;
  } else {
    PrintResult(true);
  }
  return 10;
}

uint32_t Pacose::CalculateNextSoftclauseCombination() {
  ExcludeCurrentSoftclauseCombination();
  return _satSolver->Solve();
}

uint32_t Pacose::CalculateNextResult() {
  ExcludeCurrentResult();
  return _satSolver->Solve();
}

void Pacose::DumpSolvingInformation() {
  std::cout << "c #SolverCalls...........: " << _satSolver->_noSolverCalls
            << std::endl;
  if (_settings.testIfDividable > 0) {
    std::cout << "c GBMO time..............: " << _GBMOTime << std::endl;
    std::cout << "c GBMO #partitions.......: " << _GBMOPartitions << std::endl;
  }
  if (_noTrimSATSolverCalls > 0) {
    std::cout << "c TrimSAT time...........: " << _trimSATTime << std::endl;
    std::cout << "c TrimSAT #used..........: " << _noTrimSAT << std::endl;
    std::cout << "c TrimSAT #solverCalls...: " << _noTrimSATSolverCalls
              << std::endl;
    if (_alwaysUNSATSCs > 0) {
      std::cout << "c #always UNSAT SCs......: " << _alwaysUNSATSCs
                << std::endl;
      std::cout << "c #always UNSAT weight...: " << _alwaysUNSATWeight
                << std::endl;
    }
  }

  if (_alwaysSATSCs > 0) {
    std::cout << "c #always SAT SCs........: " << _alwaysSATSCs << std::endl;
    std::cout << "c #always SAT weight.....: " << _alwaysSATWeight << std::endl;
  }

  std::cout << "c encoding #clauses......: " << _clausesOfEncoding << std::endl;
  std::cout << "c encoding #variables....: " << _variablesOfEncoding
            << std::endl;

  //  }

  std::cout << "o " << _unSatWeight << std::endl;
  std::cout << "s OPTIMUM FOUND" << std::endl;
}

void Pacose::PrintResult(bool savedModel) {
  if (!_settings._printModel)
    return;

  // if (!_settings.useMaxPre2) {
    
    // TO CHECK
    // THIS SOLVER CALL SHOULDN'T BE NECESSARY!!!
    // Assertion with Paxiant fuzzed instance with seed: 2319112478511009383 or 12384149468690201067
    // _satSolver->Solve();
    std::cout << "v ";
    if (savedModel) {
      for (uint32_t i = 1; i <= _nbVars; i++) {
        std::cout << ((_bestModel[i] ^ 1) % 2);
      }
    } else {
      for (uint32_t i = 1; i <= _nbVars; i++) {
      std::cout << ((_satSolver->GetModel(i) ^ 1) % 2);
    }
    }

    std::cout << std::endl;
    return;
  // }
  // MaxPreReconstructResult();

  // MAXPRE RECONSTRUCTION
  // std::vector<int> model;
  // for (uint32_t i = 1; i <= _nbVars; i++) {
  //   int var = (_satSolver->GetModel(i) >> 1);
  //   if ((_satSolver->GetModel(i) ^ 1) % 2 == 0)
  //     var = -var;
  //   model.push_back(var);
  // }
  // // std::cout << "v ";
  // // for (auto var : model)
  // //   std::cout << var << " ";
  // // std::cout << "<-- old Model" << std::endl;
  // std::vector<int> newModel = maxpre->reconstruct(model);
  // std::cout << "v ";
  // for (size_t i = 0; i < newModel.size(); ++i) {
  //   assert(abs(newModel[i]) == i + 1);
  //   std::cout << ((newModel[i] > 0) ? 1 : 0);
  // }
  // std::cout << std::endl;
  // // maxpre->printSolution(model, std::cout, _unSatWeight)

  // //    old print model version
  // //  std::cout << "v ";
  // //  for (uint32_t i = 1; i < _nbVars; i++) {
  // //    if ((_satSolver->GetModel(i) ^ 1) % 2 == 0) std::cout << "-";
  // //    std::cout << (_satSolver->GetModel(i) >> 1) << " ";
  // //  }
  // //  std::cout << std::endl;
}

void Pacose::CalculateOverallTimes() {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }
  DGPW::TimeVariables overallTime;
  int iter = 0;
  for (auto casc : _cascCandidates) {
    iter++;
    overallTime.AddTimeStruct(casc.dgpw->_timeVariables);
    casc.dgpw->_timeVariables->DumpVariables(iter);
  }
  overallTime.DumpVariables();
}

std::vector<uint32_t> Pacose::GetBestSCAssignment() {
  std::vector<uint32_t> SCModel;
  for (auto SC : _originalSoftClauses) {
    SCModel.push_back(SC->lastassignment);
  }
  return SCModel;
}

uint64_t
Pacose::CalculateLocalSATWeight(std::vector<SoftClause *> *tmpSoftClauses) {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }
  bool setLocalSatWeight = false;

  if (tmpSoftClauses == nullptr) {
    setLocalSatWeight = true;
    tmpSoftClauses = _actualSoftClauses;
  }

  uint64_t unSatWeight = 0;
  uint64_t satWeight = 0;

  // Process all soft clauses
  for (uint32_t i = 0; i != tmpSoftClauses->size(); ++i) {
    uint32_t relaxlit = (*tmpSoftClauses)[i]->relaxationLit;
    if (_satSolver->GetModel(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause((*tmpSoftClauses)[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (_satSolver->GetModel(clause[pos] >> 1) == clause[pos]) {
          satWeight += (*tmpSoftClauses)[i]->weight;
          // std::cout << "SATWEIGHT/CL: " << (*tmpSoftClauses)[i]->weight << std::endl;
          break;
        }
      }
      if (pos == clause.size()) {
        //        (*tmpSoftClauses)[i]->lastassignment = relaxlit;
        unSatWeight += (*tmpSoftClauses)[i]->weight;
        // std::cout << "UNSATWEIGHT: " << (*tmpSoftClauses)[i]->weight << std::endl;
      }
    } else if (_satSolver->GetModel(relaxlit >> 1) != 0) {
      assert(_satSolver->GetModel(relaxlit >> 1) == (relaxlit ^ 1));
      //      (*tmpSoftClauses)[i]->lastassignment = relaxlit ^ 1;
      satWeight += (*tmpSoftClauses)[i]->weight;
      // std::cout << "SATWEIGHT/RL: " << (*tmpSoftClauses)[i]->weight << std::endl;
    }
  }

  
  if (setLocalSatWeight) {
    if (unSatWeight < _localUnSatWeight) {
      _localUnSatWeight = unSatWeight;
      _localSatWeight = satWeight;
    }
    _currentLocalSatWeight = satWeight; 
    _currentLocalUnSatWeight = unSatWeight;
    _sumOfActualSoftWeights = satWeight + unSatWeight;
  }
  if (_settings.verbosity > 0) {
    std::cout << "c actual soft clause satWeight   / best: " << _currentLocalSatWeight << " / " << _localSatWeight << std::endl;
    std::cout << "c actual soft clause unSatWeight / best: " << _currentLocalUnSatWeight << " / " << _localUnSatWeight << std::endl;
    std::cout << "c actual soft claus sum of weights.....: " << _sumOfActualSoftWeights << std::endl;
  }

  return satWeight;
}

uint64_t Pacose::CalculateSATWeight() {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "c SC.size: " << _originalSoftClauses.size() << std::endl;
    std::cout << "c Clauses.size: " << _satSolver->GetNumberOfClauses()
              << std::endl;
  }

  if (_actualSoftClauses != nullptr) {
    CalculateLocalSATWeight();
  }

  uint64_t unSatWeight = 0;
  uint64_t satWeight = 0;

  // Process all soft clauses
  for (uint32_t i = 0; i != _originalSoftClauses.size(); ++i) {
    uint32_t relaxlit = _originalSoftClauses[i]->relaxationLit;
    if (_satSolver->GetModel(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause(_originalSoftClauses[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (_satSolver->GetModel(clause[pos] >> 1) == clause[pos]) {
          satWeight += _originalSoftClauses[i]->originalWeight;
          break;
        }
      }
      if (pos == clause.size()) {
        //        _originalSoftClauses[i]->lastassignment = relaxlit;
        unSatWeight += _originalSoftClauses[i]->originalWeight;
      }
    } else if (_satSolver->GetModel(relaxlit >> 1) != 0) {
      assert(_satSolver->GetModel(relaxlit >> 1) == (relaxlit ^ 1));
      //      _originalSoftClauses[i]->lastassignment = relaxlit ^ 1;
      satWeight += _originalSoftClauses[i]->originalWeight;
    }
  }

  //  _actualSoftClauses->size()
  if (satWeight > _satWeight || unSatWeight < _unSatWeight) {
    if (_settings.verbosity > 3)
      std::cout << "c A BETTER SOLUTION IS FOUND!!!! Overwrite _unsatweight = " << _unSatWeight << " with " << unSatWeight << " and _satWeight = " << _satWeight << " with " << satWeight << std::endl;
    // save current model!
    assert(_unSatWeight >= unSatWeight);
    SaveModel();
    _satWeight = satWeight;
    _unSatWeight = unSatWeight;
    std::cout << "o " << _unSatWeight << std::endl;
  }
  if (_settings.verbosity > 0) {
    std::cout << "c current #sat / best: " << satWeight << "/" << _satWeight << std::endl;
    std::cout << "c current o    / best: " << unSatWeight <<  "/" << _unSatWeight << std::endl;
  }
  _lastCalculatedUnsatWeight = unSatWeight;
  assert(_lastCalculatedUnsatWeight >= _unSatWeight);

  return _satWeight;
}

void Pacose::SaveModel() {
  _bestModel.clear();
  _bestModel.push_back(0);
  for (uint32_t i = 1; i <= _nbVars; i++) {
    _bestModel.push_back(_satSolver->GetModel(i));
  }
}

void Pacose::ExcludeCurrentResult() {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }
  _satSolver->NewClause();
  for (uint32_t i = 1; i < _nbVars; i++) {
    _satSolver->AddLiteral(_satSolver->GetModel(i) ^ 1);
  }
  _satSolver->CommitClause();
  _satSolver->ResetClause();
}

void Pacose::ExcludeCurrentSoftclauseCombination() {
  if (_settings.verbosity > 4) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }

  uint64_t unSatWeight = 0;
  uint64_t satWeight = 0;

  std::vector<uint32_t> resultClause;

  // Process all soft clauses
  for (uint32_t i = 0; i != _originalSoftClauses.size(); ++i) {
    uint32_t relaxlit = _originalSoftClauses[i]->relaxationLit;
    if (_satSolver->GetModel(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause(_originalSoftClauses[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (_satSolver->GetModel(clause[pos] >> 1) == clause[pos]) {
          satWeight += _originalSoftClauses[i]->originalWeight;
          if (!_originalSoftClauses[i]->positiveForcePossible) {
            for (auto lit : clause) {
              _satSolver->NewClause();
              _satSolver->AddLiteral(lit ^ 1);
              _satSolver->AddLiteral(relaxlit ^ 1);
              _satSolver->CommitClause();
              _satSolver->ResetClause();
            }
            _originalSoftClauses[i]->positiveForcePossible = true;
          }
          assert(_originalSoftClauses[i]->positiveForcePossible);
          resultClause.push_back(relaxlit);
          break;
        }
      }
      if (pos == clause.size()) {
        //        _originalSoftClauses[i]->lastassignment = relaxlit;
        unSatWeight += _originalSoftClauses[i]->originalWeight;
      }
    } else if (_satSolver->GetModel(relaxlit >> 1) != 0) {
      assert(_satSolver->GetModel(relaxlit >> 1) == (relaxlit ^ 1));
      //      _originalSoftClauses[i]->lastassignment = relaxlit ^ 1;
      satWeight += _originalSoftClauses[i]->originalWeight;
      std::vector<uint32_t> clause(_originalSoftClauses[i]->clause);
      if (!_originalSoftClauses[i]->positiveForcePossible) {
        for (auto lit : clause) {
          _satSolver->NewClause();
          _satSolver->AddLiteral(lit ^ 1);
          _satSolver->AddLiteral(relaxlit ^ 1);
          _satSolver->CommitClause();
          _satSolver->ResetClause();
        }
        _originalSoftClauses[i]->positiveForcePossible = true;
      }
      assert(_originalSoftClauses[i]->positiveForcePossible);
      resultClause.push_back(relaxlit);
    }
  }
  _satSolver->AddClause(resultClause);

  assert(satWeight == _satWeight);
  assert(unSatWeight == _unSatWeight);
}

void Pacose::DumpCascCandidates() {
  if (_settings.verbosity == 0)
    return;
  //    int i = 0;
  //    for (auto cascade : _cascCandidates) {
  //        std::cout << "c Points[" << i << "]..............: " <<
  //        cascade.Points << std::endl; std::cout << "c weightsTillPoint[" <<
  //        i
  //        << "]....: " << cascade.weightsTillPoint
  //                  << std::endl;
  //        std::cout << "c ggtTillPoint[" << i << "]........: " <<
  //        cascade.ggtTillPoint << std::endl; std::cout << "c
  //        allWeightsAreEqual[" << i << "]..: " << cascade.allWeightsAreEqual
  //                  << std::endl;

  //        i++;
  //    }

  int i = 0;
  for (auto cascade : _cascCandidates) {
    if (i == 0) {
      std::cout << "c Points[" << i << "]..............: " << cascade.Points
                << std::endl;
      std::cout << "c weightsTillPoint[" << i
                << "]....: " << cascade.weightsTillPoint << std::endl;
    } else {
      std::cout << "c Points[" << i << "]..............: "
                << _cascCandidates[i].Points - _cascCandidates[i - 1].Points
                << std::endl;
      std::cout << "c weightsTillPoint[" << i << "]....: "
                << _cascCandidates[i].weightsTillPoint -
                       _cascCandidates[i - 1].weightsTillPoint
                << std::endl;
    }
    std::cout << "c ggtTillPoint[" << i << "]........: " << cascade.ggtTillPoint
              << std::endl;
    std::cout << "c allWeightsAreEqual[" << i
              << "]..: " << cascade.allWeightsAreEqual << std::endl;

    i++;
  }
}

void Pacose::DivideSCs(std::vector<uint32_t> &sortedSCs, int acceptedMode) {
  partitionInformation nextCandidate;

  if (_cascCandidates.size() == 0) {
    _cascCandidates.push_back(nextCandidate);
  } else {
    assert(_cascCandidates.size() == 1);
    _cascCandidates[0] = nextCandidate;
  }
  _cascCandidates.back().ggtTillPoint =
      _originalSoftClauses[sortedSCs[0]]->weight;
  uint32_t j = 0;

  //    std::cout << "partitionPoints.size(): " << SCs.Points.size() <<
  //    std::endl;

  _overallSoftWeights = 0;
  //    uint64_t weightRange =
  //    _originalSoftClauses[sortedSCs.back()]->weight -
  //    _originalSoftClauses[sortedSCs[0]]->weight;

  //    uint64_t weightRange = 5;
  //    std::cout << "c weightRange: " << weightRange << std::endl;
  

  for (uint32_t i = 0; i < sortedSCs.size(); ++i) {
    j++;

    int mode = -1;
    if (acceptedMode == 0 &&
        _overallSoftWeights <= _originalSoftClauses[sortedSCs[i]]->weight &&
        i > 0 && j > 2) {
      //      std::cout << "se" << _originalSoftClauses[sortedSCs[i]]->weight
      //                << std::endl;
      if (!(_settings.testIfDividable != 2 &&
            _overallSoftWeights ==
                _originalSoftClauses[sortedSCs[i]]->weight)) {
        //        std::cout << "mod0" << std::endl;
        mode = 0;
        // if it is only the first element, then undo the mode decision unless
        // it is a strict smaller!
        if (i == 1 &&
            _overallSoftWeights == _originalSoftClauses[sortedSCs[i]]->weight) {
          //          std::cout << "mod-1" << std::endl;
          mode = -1;
        }
      }
    }
    if (acceptedMode == 1 && j > _originalSoftClauses.size() / 15 &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 15) &&
        (10 * _originalSoftClauses[sortedSCs[i - 1]]->weight) <
            _originalSoftClauses[sortedSCs[i]]->weight)
      mode = 1;
    if (acceptedMode == 2 && j > _originalSoftClauses.size() / 4 &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 3) &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 4) &&
        (2 * _overallSoftWeights < _sumOfSoftWeights - _overallSoftWeights))
      mode = 2;
    if (acceptedMode == 3 && j > _originalSoftClauses.size() / 6 &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 6) &&
        _originalSoftClauses[sortedSCs[i - 1]]->weight <
            _originalSoftClauses[sortedSCs[i]]->weight &&
        _originalSoftClauses[sortedSCs[i - 1]]->weight ==
            _cascCandidates.back().ggtTillPoint)
      mode = 3;
    if (acceptedMode == 4 && j > _originalSoftClauses.size() / 4 &&
        _originalSoftClauses.size() > 150 &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 1.5) &&
        i < _originalSoftClauses.size() - (_originalSoftClauses.size() / 4) &&
        (5 * _overallSoftWeights < _sumOfSoftWeights - _overallSoftWeights))
      mode = 4;
    if (acceptedMode == 5 && j > 20 &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 20) &&
        i < _originalSoftClauses.size() - 20 &&
        _overallSoftWeights < 2 * (_sumOfSoftWeights - _overallSoftWeights))
      mode = 5;
    if (acceptedMode == 6 && j > 20 &&
        _cascCandidates.back().allWeightsAreEqual &&
        (_originalSoftClauses[sortedSCs[i]]->weight >
         _originalSoftClauses[sortedSCs[i - 1]]->weight * 10) &&
        i < _originalSoftClauses.size() - 20 &&
        _overallSoftWeights < 5 * (_sumOfSoftWeights - _overallSoftWeights))
      mode = 6;

    // minimal size fo rone Softclause
    //    if (j < _settings.minSize) mode = -1;
    if (mode == acceptedMode) {
      j = 0;

      _cascCandidates.back().allWeightsAreEqual =
          _cascCandidates.back().ggtTillPoint ==
          _originalSoftClauses[sortedSCs[i - 1]]->weight;
      _cascCandidates.back().Points = i;
      _cascCandidates.back().weightsTillPoint = _overallSoftWeights;
      // already next cascade!
      _cascCandidates.push_back(nextCandidate);
      _cascCandidates.back().ggtTillPoint =
          _originalSoftClauses[sortedSCs[i]]->weight;
    }

    if (_cascCandidates.back().ggtTillPoint != 1) {
      _cascCandidates.back().ggtTillPoint = GreatestCommonDivisor(
          _cascCandidates.back().ggtTillPoint,
          _originalSoftClauses[sortedSCs[i]]->weight);
    }

    _overallSoftWeights += _originalSoftClauses[sortedSCs[i]]->weight;
  }

  _cascCandidates.back().allWeightsAreEqual =
      (_cascCandidates.back().ggtTillPoint ==
       _originalSoftClauses[sortedSCs.back()]->weight);

  _cascCandidates.back().Points = _originalSoftClauses.size();
  _cascCandidates.back().weightsTillPoint = _overallSoftWeights;

  //    std::cout << "_cascCandidates.size() " << _cascCandidates.size() <<
  //    std::endl;

  //    if ((_cascCandidates[0].Points < _originalSoftClauses.size() / 15)
  //        && (_cascCandidates[0].weightsTillPoint >
  //        _cascCandidates[1].ggtTillPoint)) { _cascCandidates.clear();
  //        _cascCandidates.push_back(nextCandidate);
  //    }
  //    for (uint32_t i = 1; i < _cascCandidates.size() - 1; ++i) {
  //        DumpCascCandidates();
  //        if ((_cascCandidates[i].weightsTillPoint > _cascCandidates[i +
  //        1].ggtTillPoint)) {
  //            &&(_cascCandidates[i].Points - _cascCandidates[i - 1].Points
  //               < _originalSoftClauses.size() / 15)
  //               _cascCandidates.clear();
  //            _cascCandidates.push_back(nextCandidate);
  //        }
  //        std::cout << "" << std::endl;
  //    }
}

void Pacose::RemoveCascCand(uint32_t i) {
  _cascCandidates[i - 1].Points = _cascCandidates[i].Points;
  _cascCandidates[i - 1].weightsTillPoint = _cascCandidates[i].weightsTillPoint;
  _cascCandidates[i - 1].ggtTillPoint = GreatestCommonDivisor(
      _cascCandidates[i].ggtTillPoint, _cascCandidates[i - 1].ggtTillPoint);
  _cascCandidates[i - 1].allWeightsAreEqual = false;
  _cascCandidates.erase(_cascCandidates.begin() + i);
}

bool Pacose::CheckMinWeightDist(std::vector<uint32_t> &sortedSCs,
                                uint32_t firstPoint, uint64_t biggerThan,
                                uint32_t index) {
  if (_settings.verbosity > 3) {
    //    std::cout << "c first element: "
    //              << _originalSoftClauses[sortedSCs[firstPoint]]->weight
    //              << std::endl;
    std::cout << "c biggerThan: " << biggerThan << std::endl;
  }

  // at first easy check simple distance between the weights is already
  // smaller than the needed distance
  for (uint32_t i = firstPoint; i < sortedSCs.size() - 1; i++) {
    // test if original weights have at least a distance of biggerThan
    if (_originalSoftClauses[sortedSCs[i + 1]]->weight !=
            _originalSoftClauses[sortedSCs[i]]->weight &&
        _originalSoftClauses[sortedSCs[i + 1]]->weight - biggerThan <
            _originalSoftClauses[sortedSCs[i]]->weight) {
      //      std::cout << "i: " << i << " "
      //                << _originalSoftClauses[sortedSCs[i + 1]]->weight << "
      //                "
      //                << _originalSoftClauses[sortedSCs[i]]->weight <<
      //                std::endl;
      std::cout << "c valid weight distance..: false" << std::endl;
      return false;
    }
  }

  // check all possible weight combinations
  // ATTENTION: NP hard problem!!!
  std::set<uint64_t> allWeightCombinations = {
      0, static_cast<uint64_t>(-1)};
  std::set<uint64_t>::iterator it, it2;
  std::pair<std::set<uint64_t>::iterator, bool> ret;
  for (uint32_t i = firstPoint; i < sortedSCs.size(); i++) {
    uint64_t actualWeight = _originalSoftClauses[sortedSCs[i]]->weight;
    std::set<uint64_t> allCurrentCombinations = allWeightCombinations;
    for (it = ++allCurrentCombinations.begin();
         it != --allCurrentCombinations.end(); ++it) {
      ret = allWeightCombinations.insert(actualWeight + *it);
      if (ret.second == false) {
        // element was already in, continue with next element
        continue;
      }

      std::set<uint64_t>::iterator before = ret.first;
      std::set<uint64_t>::iterator after = ret.first;

      if (ret.second && (!(*ret.first - *--before > biggerThan) ||
                         !(*++after - *ret.first > biggerThan))) {
        if (_settings.verbosity > 0)
          std::cout << "c valid weight distance..: false2" << std::endl;
        return false;
      }
    }

    ret = allWeightCombinations.insert(actualWeight);
    if (allWeightCombinations.size() > 10000000) {
      if (_settings.verbosity > 0) {
        std::cout << "c too many weight combinations [" << index << ": "
                  << allWeightCombinations.size() << std::endl;
        std::cout << "c valid weight distance..: unknown" << std::endl;
      }
      return false;
    }
    if (ret.second == false) {
      // element was already in, continue with next element
      //      std::cout << "Element was already in!" << std::endl;
      continue;
    }
    // Check weight distance for actual weight
    std::set<uint64_t>::iterator before = ret.first;
    std::set<uint64_t>::iterator after = ret.first;
    if (ret.second && (!(*ret.first - *--before > biggerThan) ||
                       !(*++after - *ret.first > biggerThan))) {
      if (_settings.verbosity > 0)
        std::cout << "c valid weight distance..: false3" << std::endl;
      return false;
    }
  }

  if (_settings.verbosity > 0)
    std::cout << "c valid weight distance..: true" << std::endl;

  std::cout << "c weight combinations[" << index
            << "].: " << allWeightCombinations.size() << std::endl;
  return true;
}

uint64_t Pacose::DivideSCsIfPossible() {
  _settings.divideDGPW = USEONLYGCD;

  if (_settings.divideDGPW == NODIVISION) {

    _sClauses.push_back(_originalSoftClauses);
    partitionInformation nextCandidate;
    _cascCandidates.push_back(nextCandidate);

    _overallSoftWeights = 0;
    for (auto sc : _originalSoftClauses) {
      _overallSoftWeights += sc->weight;
    }
    return _sClauses.size();
  }

  // if (_settings.minSize == 0)
  //   _settings.minSize =
  //     std::floor(_settings.minSize * (0.0001 * _originalSoftClauses.size()));
  if (_settings.minSize != 0)
    std::cout << "c min size of DGPW.......: " << _settings.minSize
              << std::endl;

  // get indices of sorted Bucket entries
  // generate Indice Vector to sort that vector!
  std::vector<uint32_t> sortedSCIndices(_originalSoftClauses.size());
  std::size_t n(0);
  std::generate(std::begin(sortedSCIndices),
                std::begin(sortedSCIndices) +
                    static_cast<uint32_t>(_originalSoftClauses.size()),
                [&] { return n++; });

  // stable sort - not changing order of SC's important for some of the
  // instances! especially spot5!
  std::stable_sort(std::begin(sortedSCIndices), std::end(sortedSCIndices),
                   [&](std::size_t i1, std::size_t i2) {
                     return (_originalSoftClauses[i2]->weight >
                             _originalSoftClauses[i1]->weight);
                   });
  if (_settings.divisionMode == -1) {
    for (int mode = 0; mode <= 6; mode++) {
      DivideSCs(sortedSCIndices, mode);
      if (_cascCandidates.size() > 1 ||
          (mode == 0 && _settings.divideDGPW == USEONLYGCD)) {
        if (_settings.verbosity > 0)
          std::cout << "c Accepted Mode..........: " << mode << std::endl;
        break;
      }
    }
  } else if (_settings.divisionMode >= 0 && _settings.divisionMode <= 6) {
    DivideSCs(sortedSCIndices, _settings.divisionMode);
  } else {
    std::cout << "c ERROR: no valid division!!" << std::endl;
    exit(0);
  }
  if (_cascCandidates.size() > 1) {
    //        std::cout << "c partitionFactor........: " <<
    //        _settings.partitionFactor << std::endl;

    // test if gcd of the bigger potential cascades i greater than the sum of
    // weights of the lower cascades
    for (uint32_t i = static_cast<uint32_t>(_cascCandidates.size() - 1); i > 0;
         --i) {
      uint64_t minWeightDistance = _cascCandidates[i - 1].weightsTillPoint;
      //      std::cout << "c wTP[" << i - 1 << "]!"
      //                << _cascCandidates[i - 1].weightsTillPoint <<
      //                std::endl;

      //      std::cout << "ggtTillPoint: " << _cascCandidates[i].ggtTillPoint
      //                << std::endl;
      //      std::cout << "minSize: " << _settings.minSize << std::endl;
      //      std::cout << "Points: "
      //                << _cascCandidates[i].Points - _cascCandidates[i -
      //                1].Points
      //                << std::endl;

      if (_settings.testIfDividable == 2) {
        // test for bigger equal
        minWeightDistance = _cascCandidates[i - 1].weightsTillPoint - 1;
      }

      bool isBigger = false;
      if ((_cascCandidates[i].Points - _cascCandidates[i - 1].Points >=
           _settings.minSize) &&
          (i != 1 || _cascCandidates[0].Points >= _settings.minSize))
        isBigger = true;

          //  std::cout << "i " << i << " isbigger: " << isBigger <<
          //  std::endl; std::cout << "i " << i << " minWeightDistance: " <<
          //  minWeightDistance
          //            << std::endl;
          //  std::cout << "i " << i << " _cascCandidates[i].ggtTillPoint: "
          //            << _cascCandidates[i].ggtTillPoint << std::endl;

      if ((_cascCandidates[i].ggtTillPoint > minWeightDistance) && isBigger) {
        if (_settings.verbosity > 0) {
          std::cout << "c VALID subCascade[" << i << "]!" << std::endl;
          std::cout << "c SPLIT due to GCD!" << std::endl;
        }
      } else if (_settings.exhaustiveWeightDistCheck) {
        if (!CheckMinWeightDist(sortedSCIndices, _cascCandidates[i - 1].Points,
                                minWeightDistance, i) ||
            !isBigger) {
          if (_settings.verbosity > 0)
            std::cout << "c INVALID subCascade[" << i << "]!!" << std::endl;
          RemoveCascCand(i);
        } else {
          if (_settings.verbosity > 0) {
            std::cout << "c VALID subCascade[" << i << "]!!" << std::endl;
            std::cout << "c SPLIT due to exhaustive weight distance check!" << std::endl;
          }
        }
      } else {
        if (_settings.verbosity > 0)
            std::cout << "c INVALID subCascade[" << i << "]!!" << std::endl;
          RemoveCascCand(i);
      }
    }
    if (_cascCandidates.size() > 1) {
      if (_settings.verbosity > 0)
        std::cout << "c number of subproblems..: " << _cascCandidates.size()
                  << std::endl;
    }

    _settings.formulaIsDivided = true;
    uint32_t nextPartition = 0;
    std::vector<SoftClause *> SC;
    _sClauses.push_back(SC);
    for (uint32_t i = 0; i < sortedSCIndices.size(); ++i) {
      if (i == _cascCandidates[nextPartition].Points) {
        nextPartition++;
        _sClauses.push_back(SC);
      }
      _sClauses.back().push_back(_originalSoftClauses[sortedSCIndices[i]]);
    }
  }

  if (_cascCandidates.size() <= 1) {
    // COPY CONSTRUCTOR??
    if (_sClauses.size() == 0) {
      std::vector<SoftClause *> SC;
      _sClauses.push_back(SC);
    }
    _sClauses.back() = _originalSoftClauses;
  } else {
    DumpCascCandidates();
  }

  if (_settings.verbosity > 3) {
    std::cout << std::endl;
    for (auto i : _sClauses) {
      std::cout << "c ";
      for (auto j : i) {
        std::cout << j->weight << " ";
      }
      std::cout << std::endl << std::endl;
    }
  }

  if (_settings.divCheck)
    exit(0);

  return _sClauses.size();
}

void Pacose::Preprocess() {
  //    GreedyMaximizeInitialSATWeight();

  // not compatible with dividing SC's into subformulas!
  if (_settings.divideDGPW == NODIVISION && !_settings.createSimplifiedWCNF) {
    AnalyzeSCsAndConvertIfPossible();
  }
  uint64_t gcdBefore = _GCD;
  CalcGCDAndDivideIfPossible();
  if (gcdBefore != _GCD && _localUnSatWeight!= UINT64_MAX) {
    if (_settings.verbosity > 2) 
      std::cout << "c GCD changed, recalculate (UN)SAT weight:" << _GCD << std::endl;
    assert(_GCD % gcdBefore == 0);
    _localSatWeight /= (_GCD / gcdBefore);
    _localUnSatWeight /= (_GCD / gcdBefore);
  }

}

void Pacose::AnalyzeSCsAndConvertIfPossible() {

  if (_settings.verbosity > 8)
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  if (!_settings.GetAnalyzeFormula())
    return;

  if (_actualSoftClauses->size() == 0) {
    _settings.SetFormulaType(FormulaType::NORMALSAT);
    std::cout << "c is pure SAT formula......: true" << std::endl;
    return;
  }

  uint64_t minWeight = (*_actualSoftClauses)[0]->weight;
  uint64_t maxWeight = (*_actualSoftClauses)[0]->weight;
  uint64_t prevMaxWeight = (*_actualSoftClauses)[0]->weight;
  bool hasMoreThanTwoWeights = false;
  bool recalc = false;
  uint64_t sumOfMinWeights = 0;

  for (size_t i = 0; i < _actualSoftClauses->size(); i++) {
    //        std::cout << (*_actualSoftClauses)[i]->weight << " ";

    if ((*_actualSoftClauses)[i]->weight < minWeight) {
      minWeight = (*_actualSoftClauses)[i]->weight;
      recalc = true;
    } else if ((*_actualSoftClauses)[i]->weight > maxWeight) {
      prevMaxWeight = maxWeight;
      maxWeight = (*_actualSoftClauses)[i]->weight;
    }

    if (minWeight != maxWeight && minWeight != prevMaxWeight &&
        maxWeight != prevMaxWeight) {
      hasMoreThanTwoWeights = true;
    }

    if ((*_actualSoftClauses)[i]->weight == minWeight) {
      if (recalc) {
        sumOfMinWeights = 0;
        recalc = false;
      }
      sumOfMinWeights += minWeight;
    }
  }

  //  recognize maxSAT instances
  if (minWeight == maxWeight) {
    //    if (!_hasHardClauses) {
    if (false) {
      for (int i = static_cast<int>(_actualSoftClauses->size() - 1); i < 0;
           i--) {
        if ((*_actualSoftClauses)[static_cast<size_t>(i)]->clause.size() == 1) {
          // wrongly interpreted as unit soft clause, has to be added as unit
          // hard clause
          _satSolver->NewClause();
          uint32_t ulit =
              (*_actualSoftClauses)[static_cast<size_t>(i)]->clause[0];
          _satSolver->AddLiteral(&ulit);
          _satSolver->CommitClause();
          _satSolver->ResetClause();

          //                    _satSolver->AddClause(_softClauses[i]->clause);
        } else {
          // deactivate trigger literal by adding it as negated unit clause.
          _satSolver->NewClause();
          uint32_t ulit =
              (*_actualSoftClauses)[static_cast<size_t>(i)]->relaxationLit ^ 1;
          _satSolver->AddLiteral(&ulit);
          _satSolver->CommitClause();
          _satSolver->ResetClause();
        }
      }
      _hasHardClauses = true;
      _settings.SetFormulaType(FormulaType::NORMALSAT);
      std::cout << "c converted to SAT.......: true" << std::endl;
      // CAN BE CONTINUED WITH SAT SOLVER CALL, NOT YET IMPLEMENTED!
    } else {
      _settings.SetFormulaType(FormulaType::MAXSAT);
      std::cout << "c is MaxSAT formula......: true" << std::endl;
    }
  }

  //    // has exactly two weights and no hard clauses - it can be converted
  //    to pure unweighted maxSAT std::cout << "!hasMoreThanTwoWeights: " <<
  //    !hasMoreThanTwoWeights
  //              << "  !_hasHardClauses: " << _hasHardClauses
  //              << "  minWeight: " << minWeight
  //              << "  maxWeight: " << maxWeight
  //              << "  sumOfMinWeights: " << sumOfMinWeights
  //              << "  maxWeight: " << maxWeight
  //              << std::endl;

  if (!hasMoreThanTwoWeights && !_hasHardClauses && minWeight != maxWeight &&
      sumOfMinWeights < maxWeight) {
    std::vector<SoftClause *> newSoftClauses;
    //        std::vector< SoftClause* > newSoftClauses;
    //        std::vector<std::vector<u_int32_t>> newLiterals;

    //        std::cout << "maxWeight " << maxWeight << std::endl;
    for (int ind = 0; ind < static_cast<int>(_actualSoftClauses->size());
         ++ind) {
      if ((*_actualSoftClauses)[static_cast<size_t>(ind)]->weight !=
          maxWeight) {
        newSoftClauses.push_back(
            (*_actualSoftClauses)[static_cast<size_t>(ind)]);
      }
      //            else if
      //            ((*_actualSoftClauses)[static_cast<size_t>(ind)]->clause.size()
      //            == 1)
      //            {
      //                std::cout << "in here" << std::endl;
      //                // wrongly interpreted as unit soft clause, has to be
      //                added as unit hard clause _satSolver->NewClause();
      //                uint32_t ulit =
      //                (*_actualSoftClauses)[static_cast<size_t>(ind)]->clause[0];
      //                _satSolver->AddLiteral(&ulit);
      //                _satSolver->CommitClause();
      //                _satSolver->ResetClause();
      //            }
      else {
        //                std::cout << (*_actualSoftClauses)[ind]->weight << "
        //                "
        //                << (*_actualSoftClauses)[ind]->relaxationLit << " ";
        // deactivate trigger literal by adding it as negated unit clause.
        _satSolver->NewClause();
        uint32_t ulit =
            (*_actualSoftClauses)[static_cast<size_t>(ind)]->relaxationLit ^ 1;
        //                std::cout << ulit << ", " << std::endl;
        _satSolver->AddLiteral(&ulit);
        _satSolver->CommitClause();
        _satSolver->ResetClause();
      }
    }
    _hasHardClauses = true;
    //        _sClauses[1] = newSoftClauses;
    (*_actualSoftClauses) = newSoftClauses;

    _settings.SetFormulaType(FormulaType::MAXSAT);
    std::cout << "c converted to MaxSAT....: true" << std::endl;
    std::cout << "c Remaining SoftClauses..: " << _actualSoftClauses->size()
              << std::endl;

    //        for( auto SC : *_actualSoftClauses ) {
    //            std::cout << SC->weight << " ";
    //        }
    //        std::cout << std::endl;

    _sumOfSoftWeights = sumOfMinWeights;
  }
}

uint64_t Pacose::GreatestCommonDivisor(uint64_t a, uint64_t b) {
  if (a == b) {
    return a;
  }
  //    std::cout << " a: " << a << "  b: " << b << std::endl;
  uint64_t temp;
  while (b > 0) {
    temp = b;
    b = a % b;
    a = temp;
  }
  //    std::cout << "xa: " << a << "  b: " << b << std::endl;
  return a;
}

void Pacose::CalcGCDAndDivideIfPossible() {
  if (!_settings.GetAnalyzeFormula())
    return;

  _GCD = (*_actualSoftClauses)[0]->originalWeight;

  for (size_t ind = 1; ind < _actualSoftClauses->size(); ++ind) {
    if (_GCD == 1)
      break;
    _GCD = GreatestCommonDivisor(_GCD, (*_actualSoftClauses)[ind]->originalWeight);
  }

  if (_GCD > 1) {
    //        std::cout << "c greatest common divisor: " << _GGT << std::endl;
    for (size_t ind = 0; ind < _actualSoftClauses->size(); ++ind) {
      (*_actualSoftClauses)[ind]->weight = (*_actualSoftClauses)[ind]->originalWeight / _GCD;
    }
  }
}

void Pacose::SetSumOfSoftWeights(uint64_t softWeights) {
  _sumOfSoftWeights = softWeights;
}

} // Namespace Pacose
