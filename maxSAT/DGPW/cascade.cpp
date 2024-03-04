/********************************************************************************************
Copyright (c) 2014-2016, Sven Reimer,
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
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
#include <algorithm> //std::sort
#include <cmath>
#include <iomanip>
#include <numeric>
#include <unistd.h>

// Include dgpw related headers.
#include "bucket.h"
#include "cascade.h"
#include "dgpw.h"
#include "multiplecascade.h"
#include "softclausenodes.h"
#include "timemeasurement.h"
#include "timevariables.h"
#include "totalizerencodetree.h"
#include "Pacose.h"

namespace Pacose {
namespace DGPW {

// Constructor
Cascade::Cascade(DGPW *dgpw, MultipleCascade *multipleCascade, bool onlyByTares)
    : _dgpw(dgpw),
      //_control(dgpw->_control),
      _setting(dgpw->_dgpwSetting), _multipleCascade(multipleCascade),
      _base(_setting->base), _onlyByTares(onlyByTares), _maxPos(-1),
      _satWeight(0), _tareWeight(0), _weightToSubstract(0),
      _sumOfSoftWeights(0), _softClauseTreeCreated(false),
      _collectedCascadeAssumptions(), _structure(), _numberOfBuckets(0),
      _totalBucketEntriesperWeight(0), _totalBucketOccurrences(0),
      _totalBucketEntries(), _maxSorterDepth(0), _highestBucketMultiplicator(0),
      _upperWeightBoundAllLowerCascades(0), _howManyPositionsCuttedAtBottom(0),
      _softClauses(), _softClauseTree(), _processingSoftClauseTree(),
      _processingPercentOffTree(),
      _howOftenReinsertedFromProcessingPercentOffTree(0), _tareAssumptions(),
      _fixedTareAssumption() {
  assert(dgpw != nullptr);
}

Cascade::~Cascade() {
  for (auto SCT : _softClauseTree) {
    delete SCT;
  }
  for (auto bucket : _structure) {
    delete bucket;
  }
}

void Cascade::IncrementalReset() {
  //  _maxPos = -1;
  _satWeight = 0;
  //  _tareWeight(0),
  //  _weightToSubstract(0),
  //   _sumOfSoftWeights(0),
  _softClauseTreeCreated = true;

  _upperWeightBoundAllLowerCascades = 0;
  //  _howManyPositionsCuttedAtBottom(0),
  _tareAssumptions = {};
  _fixedTareAssumption = {};
  for (auto bucket : _structure) {
    bucket->IncrementalReset();
  }
}

void Cascade::Fill(std::vector<SoftClause *> *softClauses,
                   PartitionStrategy partitionStrategy,
                   EncodeStrategy encodeStrategy) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_onlyByTares)
    encodeStrategy = ENCODEONLYIFNEEDED;

  TimeMeasurement timeFillingBuckets(&_dgpw->_timeVariables->fillingBuckets,
                                     true);

  CountSumOfSoftWeights(softClauses);

  FillStructure(partitionStrategy, encodeStrategy);

  if (_setting->verbosity > 3)
    std::cout << std::endl
              << "Buckets are filled - structure is dumped" << std::endl
              << std::endl;
}

void Cascade::CountSumOfSoftWeights(std::vector<SoftClause *> *softClauses) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  _softClauses = *softClauses;
  _sumOfSoftWeights = 0;

  // sum over all SC weights of CASCADE!
  std::for_each(_softClauses.begin(), _softClauses.end(),
                [&](SoftClause *SC) { _sumOfSoftWeights += SC->weight; });
  _satWeight = CountSatisfiedSoftClauses(nullptr, _dgpw->_lastModel);

  if (_setting->verbosity > 4) {
    std::cout << "c SAT Weight of Cascade..: " << _satWeight << std::endl;
    std::cout << "c SAT Weight of dgpw.....: " << _dgpw->_satWeight
              << std::endl;
    std::cout << "c softWeights of Cascade.: " << _sumOfSoftWeights
              << std::endl;
    std::cout << "c softWeights of dgpw....: " << _dgpw->_sumOfSoftWeights
              << std::endl;
  }
}

uint32_t Cascade::CutMaxPos(bool solve) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  _structure.back()->_isLastBucket = true;
  return _structure.back()->CutMaxPos(solve);
}

uint32_t Cascade::CutMinPos(bool solve) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  _structure.back()->_isLastBucket = true;
  return _structure.back()->CutMinPos(solve);
}

void Cascade::FillStructure(PartitionStrategy partitionStrategy,
                            EncodeStrategy encodeStrategy) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_onlyByTares)
    encodeStrategy = ENCODEONLYIFNEEDED;

  PartitionSoftClauses(partitionStrategy);
  FillBuckets();

  AddTaresToBuckets();

  // Cone of influence encoding == ENCODEONLYIFNEEDED
  if (encodeStrategy == ENCODEONLYIFNEEDED) {
    UnionBucketsIntoLast();
    if (_onlyByTares)
      AddAsManyBucketsAsPossible();
    DumpBucketStructure(true, 3);

  } else {
    DumpBucketStructure(false, 3);
  }
}

void Cascade::AddAsManyBucketsAsPossible() {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  //    std::cout << "size: " << _structure.back()->size() << std::endl;
  //    std::cout << "_structure.back()->_tares[0]: " <<
  //    _structure.back()->_tares[0] << std::endl;

  if (_structure.back()->_tares.empty())
    AddTare(_structure.size() - 1);

  if (_setting->interimResult == CUTATTOP) {
    CutMaxPos();
  }
  _structure.back()->_encodeTreeGenerated = false;

  /** ATTENTION NOT ALWAYS GIVEN
   *
   * If e.g. highest number is close to uint64_t...
   */
  if (AddNewBucketsTillMultiplicatorMatches(static_cast<uint64_t>(-1), true)) {
    std::cout << "EXCEPTION - cannot be solved only by tares because of "
                 "multiplicator limit of "
                 "64 Bit!";
    _onlyByTares = false;
    return;
  }
  DumpBucketStructure(true);

  // _estimatedWeightBoundaries[0] = -_highestBucketMultiplicator + 1;
  _estimatedWeightBoundaries[0] = 0;
  _estimatedWeightBoundaries[1] = _highestBucketMultiplicator;
}

bool Cascade::Encode() {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  switch (_setting->encodeStrategy) {
  case ENCODEALL:
    EncodeTopBuckets();
    /*
    if ( _control->ReachedLimits() )
        return false;
            */

    EncodeBottomBuckets();
    CalculateBucketEntries();
    DumpBucketStructure(false, 4);
    break;
  case ENCODEONLYIFNEEDED:
    if (_onlyByTares)
      return true;
    CreateTotalizerEncodeTree();
    CalculateBucketEntries();

    DumpBucketStructure(true, 4);
    break;
  }
  /*
  if ( _control->ReachedLimits() )
      return false;
      */

  if (_setting->verbosity < 4)
    return true;

  std::cout << std::endl
            << "Buckets are encoded - structure is dumped" << std::endl
            << std::endl;
  return true;
}

uint32_t Cascade::Solve(bool onlyWithAssumptions, bool solveTares) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  uint32_t currentresult(1);

  if (_dgpw->_satWeight == _dgpw->_sumOfSoftWeights) {
    if (_setting->verbosity > 10) {
      std::cout << "_dgpw->_satWeight == _dgpw->_sumOfSoftWeights" << std::endl;
    }
    _dgpw->FixAllSoftClauses();
    _maxPos = _structure.back()->size() - 1;
    return SAT;
  }

  //    TimeMeasurement
  //    TimeSolvingLastBucket(&_dgpw->_timeVariables->solvingLastBucket);
  // STANDARD solving "LAST BUCKET"
  if (!_onlyByTares) {
    _maxPos = static_cast<int>(_structure.back()->SolveBucketReturnMaxPosition(
        onlyWithAssumptions, false));
    
    if (_maxPos != -1 && onlyWithAssumptions)
      _fixedTareAssumption.push_back(
          (_structure.back()->_sorter->GetOrEncodeOutput(_maxPos) << 1) ^ 1);
    //        std::cout << "literal of mP == " << litmP << std::endl <<
    //        std::endl;
    //        _dgpw->AddUnit((_structure.back()->_sorter->GetOrEncodeOutput(_maxPos)
    //        << 1) ^ 1); exit(0);

    if (_dgpw->_resultUnknown) {
      vPL->write_comment("SHOULDNEVERHAPPEN: result = unknown");
      vPL->write_fail();
      return UNKNOW;
    }
    if (_setting->encodeStrategy == ENCODEONLYIFNEEDED &&
        _setting->createGraphFile != "")
      _structure.back()->_sorter->_outputTree->DumpOutputTree(
          _setting->createGraphFile + "_withOutputs.tgf", true);

    //    if (_setting->encodeStrategy != ENCODEALL) {
    //      std::cout << "c #clauses of coding.....: " << _dgpw->_addedClauses
    //                << std::endl;
    //      std::cout << "c #variables of coding...: "
    //                << _dgpw->Variables() - _dgpw->_variablesBefore <<
    //                std::endl;
    //    }
  }
  if (_dgpw->_satWeight == _dgpw->_sumOfSoftWeights) {

    if (_setting->verbosity > 10) {
      std::cout << "_dgpw->_satWeight == _dgpw->_sumOfSoftWeights2" << std::endl;
    }
    _dgpw->FixAllSoftClauses();
    return currentresult;
  }
  //        std::cout << "_dgpw->_satWeight != _dgpw->_sumOfSoftWeights" <<
  //        std::endl;

  //    _dgpw->_glucose->printIncrementalStats(1);
  //      TimeSolvingLastBucket.~TimeMeasurement();
  if (solveTares) {

    // standard solving tares
    currentresult = SolveTares(onlyWithAssumptions);
    //        std::cout << "c CURRENTRESULT: AFTER SOLVE TARES! " <<
    //        currentresult << std::endl;
  }
  //    currentresult = SolveAllTares();

  return currentresult;
}

void Cascade::CreateSoftClauseTree(std::vector<SoftClause *> *softClauses,
                                   bool split) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  _processingSoftClauseTree.clear();
  _softClauseTree.clear();

  // sort without changing order of elements
  //    std::stable_sort(softClauses->begin(), softClauses->end(),
  //    SoftClause::bigger ); _numberOfBuckets =
  //    static_cast<uint32_t>(floor(log2(_softClauses[0]->weight)/log2(_base)));

  //    if (_setting->featureTest)
  //    {
  //        //sort with changing order of elements (better worst case runtime)
  //        std::sort(softClauses->begin(), softClauses->end(),
  //        SoftClause::bigger ); _numberOfBuckets =
  //        static_cast<uint32_t>(floor(log2(_softClauses[0]->weight)/log2(_base)));
  //    } else
  //    {
  uint64_t maxValue = 0;
  for (auto softclause : *softClauses) {
    maxValue = softclause->weight > maxValue ? softclause->weight : maxValue;
  }
  _numberOfBuckets = static_cast<uint32_t>(floor(log2(maxValue) / log2(_base)));

  // SoftClauseNode Structure is created. Extract function!
  for (uint32_t i = 0; i != softClauses->size(); ++i) {
    SoftClauseNodes *softClauseNode =
        new SoftClauseNodes((*softClauses)[i], _base);
    // softclause is only in one bucket
    if (split && softClauseNode->inHowManyBuckets != 1) {
      _processingSoftClauseTree.push_back(softClauseNode);
    } else {
      _softClauseTree.push_back(softClauseNode);
    }
  }
  if (split)
    _softClauseTreeCreated = true;
}

void Cascade::PartitionSoftClauseTree(
    std::vector<SoftClauseNodes *> *tmpSoftClausesTree) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  _processingSoftClauseTree.clear();
  _softClauseTree.clear();

  // SoftClauseNode Structure is created. Extract function!
  for (uint32_t i = 0; i != tmpSoftClausesTree->size(); ++i) {
    // SoftClauseNodes* softClauseNode = new SoftClauseNodes(_softClauses[i],
    // _base);
    // softclause is only in one bucket

    if ((*tmpSoftClausesTree)[i]->inHowManyBuckets != 1) {
      _processingSoftClauseTree.push_back((*tmpSoftClausesTree)[i]);
    } else {
      _softClauseTree.push_back((*tmpSoftClausesTree)[i]);
    }
  }
  _softClauseTreeCreated = true;
}

std::vector<std::vector<uint32_t>> Cascade::GetTareVectors() {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<std::vector<uint32_t>> tares;

  uint32_t addNoTareToLastBucket = (_setting->cascadeDivider > 0) ? 0 : 1;
  if (_setting->verbosity > 3)
    std::cout << "addNoTareToLastBucket: " << addNoTareToLastBucket
              << std::endl;
  for (uint32_t ind = 0; ind < _structure.size() - addNoTareToLastBucket;
       ind++) {
    tares.push_back(_structure[ind]->_tares);
  }
  if (_setting->verbosity > 3)
    std::cout << "returning tares" << std::endl;
  return tares;
}

std::vector<std::pair<uint64_t, uint32_t>>
Cascade::GetTareVector(uint64_t weightDiff) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  int64_t wd = static_cast<int64_t>(weightDiff);
  std::vector<std::pair<uint64_t, uint32_t>> tares = {};
  if (_softClauses.size() == 0) {
    return tares;
  }

  uint32_t addNoTareToLastBucket = (_setting->cascadeDivider > 0) ? 0 : 1;
  if (_setting->verbosity > 3)
    std::cout << "addNoTareToLastBucket: " << addNoTareToLastBucket
              << std::endl;
  for (uint32_t ind = 0; ind < _structure.size() - addNoTareToLastBucket;
       ind++) {
    if (wd > 0) {
      tares.push_back(std::pair<uint64_t, uint32_t>(
          _structure[ind]->_multiplicator * _dgpw->_greatestCommonDivisor,
          _structure[ind]->_tares[0] << 1));
      // if I only solve the watchdogs I have to return all tares!
      // else test if tare is used to minimize weight, then we can lift this
      // weight!
      if (!(_dgpw->_dgpwSetting->divideDGPW == DIVIDEALLSOLVEONLYWATCHDOGS) &&
          _dgpw->Model(_structure[ind]->_tares[0]) ==
              ((_structure[ind]->_tares[0] << 1) ^ 0)) {
        wd -= static_cast<int64_t>(_structure[ind]->_multiplicator) *
              _dgpw->_greatestCommonDivisor;
      }
      if (_setting->verbosity > 4)
        std::cout << "WD: " << wd
                  << "  multi: " << _structure[ind]->_multiplicator
                  << "  tare: " << _structure[ind]->_tares[0] * 2 << std::endl;
    } else {
      // tare doesn't need to be added anymore, can be set to the model value!
      _dgpw->AddUnit(_dgpw->Model(_structure[ind]->_tares[0]));
    }
  }

  if (_setting->verbosity > 3)
    std::cout << "c returning tares" << std::endl;
  return tares;
}

std::vector<std::pair<uint64_t, uint32_t>>
Cascade::GetWatchdogs(uint64_t weightDiff) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  std::vector<std::pair<uint64_t, uint32_t>> watchdogs = {};
  if (_softClauses.size() == 0) {
    return watchdogs;
  }

  if (_setting->verbosity > 1) {
    std::cout
        << "(_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - weightDiff: "
        << (_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - weightDiff
        << std::endl;
    std::cout << "_dgpw->_greatestCommonDivisor * (_dgpw->_satWeight - "
                 "_estimatedWeightBoundaries[0]): "
              << _dgpw->_greatestCommonDivisor *
                     (_dgpw->_satWeight - _estimatedWeightBoundaries[0])
              << std::endl;

    std::cout << "((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - "
                 "weightDiff) / "
                 "(_structure.back()->_multiplicator * "
                 "_dgpw->_greatestCommonDivisor): "
              << ((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) -
                  weightDiff) /
                     (_structure.back()->_multiplicator *
                      _dgpw->_greatestCommonDivisor)
              << std::endl;
    std::cout << "floor(((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - "
                 "weightDiff) / "
                 "(_structure.back()->_multiplicator * "
                 "_dgpw->_greatestCommonDivisor)): "
              << floor(((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) -
                        weightDiff) /
                       (_structure.back()->_multiplicator *
                        _dgpw->_greatestCommonDivisor))
              << std::endl;
  }

  uint32_t firstWatchdog = 0;
  uint32_t lastWatchdog =
      static_cast<uint32_t>(_maxPos) < _structure.back()->size()
          ? static_cast<uint32_t>(_maxPos)
          : _structure.back()->size() - 1;

  if (static_cast<uint64_t>(_dgpw->_greatestCommonDivisor) * _dgpw->_satWeight >
      weightDiff) {
    firstWatchdog = static_cast<uint32_t>(floor(
        ((_dgpw->_greatestCommonDivisor * _dgpw->_satWeight) - weightDiff) /
        (_structure.back()->_multiplicator * _dgpw->_greatestCommonDivisor)));
  } /*else {
      firstWatchdog = 0;
  }*/

  if (_setting->verbosity > 0) {
    std::cout << "firstWatchdog: " << firstWatchdog << std::endl;
    std::cout << "lastWatchdog: " << lastWatchdog << std::endl;
    std::cout << "(_dgpw->_greatestCommonDivisor * _dgpw->_satWeight): "
              << _dgpw->_greatestCommonDivisor * _dgpw->_satWeight << std::endl;
    std::cout << std::endl;
    std::cout << "weightDiff: " << weightDiff << std::endl;
    std::cout << "_structure.back()->_multiplicator: "
              << _structure.back()->_multiplicator << std::endl;
    std::cout << "_structure.back()->_localMaxPos: "
              << _structure.back()->_localMaxPos << std::endl;
    std::cout << "weight boundaries: ( " << _estimatedWeightBoundaries[0]
              << " / " << _estimatedWeightBoundaries[1] << " )" << std::endl;
    std::cout << "_maxPos: " << _maxPos << std::endl;
    std::cout << "dgpwSatWeight: " << _dgpw->_satWeight << std::endl;
    std::cout << "_GGT: " << _dgpw->_greatestCommonDivisor << std::endl
              << std::endl;
  }
  assert(firstWatchdog <= lastWatchdog);

  if (_maxPos >= 0 || (firstWatchdog) > 0) {
    std::cout << "c watchdog of position " << firstWatchdog
              << " added as unit clause! " << std::endl;
    _dgpw->AddUnit(
        (_structure.back()->_sorter->GetOrEncodeOutput(firstWatchdog) << 1) ^
        1);
  }

  firstWatchdog++;
  if (firstWatchdog <= lastWatchdog) {
    for (uint32_t i = firstWatchdog; i <= lastWatchdog; ++i) {
      watchdogs.push_back(std::pair<uint64_t, uint32_t>(
          _structure.back()->_multiplicator * _dgpw->_greatestCommonDivisor,
          (_structure.back()->_sorter->GetOrEncodeOutput(i) << 1) ^ 1));
    }
    //    std::cout << "watchdogs.back(): " << watchdogs.back().second <<
    //    std::endl;
  }

  if (_setting->verbosity > 0) {
  }
  return watchdogs;
}

std::vector<uint32_t> Cascade::GetLastAssumptions() {
  //  std::cout << "_fixedTareAssumption: " << _fixedTareAssumption.size()
  //            << std::endl;

  if (_fixedTareAssumption.empty() && _maxPos < (int)_structure.size()) {
    _fixedTareAssumption = _structure.back()->GetAssumptions(_maxPos);
  }

  //  std::cout << "_fixedTareAssumption: " << _fixedTareAssumption.size() << ":
  //  "; for (auto i : _fixedTareAssumption) {
  //    std::cout << i << " ";
  //  }
  //  std::cout << std::endl;
  return _fixedTareAssumption;
}

void Cascade::PartitionSoftClauses(PartitionStrategy partitionStrategy) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (!_softClauseTreeCreated)
    CreateSoftClauseTree(&_softClauses, true);

  //    std::cout << "start sorting!" << std::endl;
  if (partitionStrategy == GROUPBYWEIGHT &&
      _setting->atLeastnEqualWeights > 1) {
    // test if there are at least 8 equal elements - otherwise don't sort!
    // get indices of sorted Bucket entries
    // generate Indice Vector to sort that vector!
    std::vector<uint32_t> sortedSCIndices(_softClauses.size());
    std::size_t n(0);
    std::generate(std::begin(sortedSCIndices),
                  std::begin(sortedSCIndices) +
                      static_cast<uint32_t>(_softClauses.size()),
                  [&] { return n++; });

    // stable sort - not changing order of SC's important for some of the
    // instances! especially spot5!
    std::stable_sort(std::begin(sortedSCIndices), std::end(sortedSCIndices),
                     [&](std::size_t i1, std::size_t i2) {
                       return (_softClauses[i2]->weight >
                               _softClauses[i1]->weight);
                     });

    //        std::cout << "before for!" << std::endl;
    //        std::cout << "softclauses.size(): " << _softClauses.size() <<
    //        std::endl;
    partitionStrategy = NOPARTITION;
    for (uint32_t i = _setting->atLeastnEqualWeights; i < _softClauses.size();
         ++i) {
      //            std::cout << "i: " << i << "  weight: " <<
      //            _softClauses[sortedSCIndices[i]]->weight << "   i-n: " <<
      //            i-_setting->atLeastnEqualWeights << "  weight: " <<
      //            _softClauses[sortedSCIndices[i-_setting->atLeastnEqualWeights]]->weight
      //            << "  ((_softClauses[sortedSCIndices[i]]->weight &
      //            (_softClauses[sortedSCIndices[i]]->weight - 1)): " <<
      //            (_softClauses[sortedSCIndices[i]]->weight &
      //            (_softClauses[sortedSCIndices[i]]->weight - 1)) <<
      //            std::endl;
      if (_softClauses[sortedSCIndices[i]]->weight ==
              _softClauses[sortedSCIndices[i - _setting->atLeastnEqualWeights]]
                  ->weight &&
          // weight is not power of two!
          !((_softClauses[sortedSCIndices[i]]->weight &
             (_softClauses[sortedSCIndices[i]]->weight - 1)) == 0)) {
        partitionStrategy = GROUPBYWEIGHT;
        break;
      }
    }
    //        std::cout << "PartitionStrategy changed to: " << partitionStrategy
    //        << std::endl;
  }

  switch (partitionStrategy) {
  case NOPARTITION:
    _softClauseTree.insert(_softClauseTree.end(),
                           _processingSoftClauseTree.begin(),
                           _processingSoftClauseTree.end());
    _processingSoftClauseTree.clear();
    break;
  // both cases have the same grouping algorithm, but differ in the way they
  // connect the weights.
  case GROUPBYWEIGHTADDATLAST:
  case GROUPBYWEIGHT:
    GroupByWeight();
    _softClauseTree.insert(_softClauseTree.end(),
                           _processingSoftClauseTree.begin(),
                           _processingSoftClauseTree.end());
    _processingSoftClauseTree.clear();
    break;
  case GROUPBYBIGGESTREPEATINGENTRY:
    //        _setting->verbosity = 7;
    std::cout << "c GROUPBYBIGGESTREPEATINGENTRY" << std::endl;
    GroupByWeight();

    // actualize values as sum of both trees.
    CalculateTotalBucketEntries(&_processingSoftClauseTree, false);
    CalculateTotalBucketEntries(&_softClauseTree, true);

    DumpSCNodeStructure(&_processingSoftClauseTree, 5);

    GroupByBiggestRepeatingEntry();

    assert(_processingSoftClauseTree.empty());
    //        _setting->verbosity = 0;
    break;
  case GROUPBYCOLUMNS:

    std::cout << "c GROUPBYCOLUMNS" << std::endl;
    //        _setting->verbosity = 7;
    if (_dgpw->_dgpwSetting->featureTest == 1) {
      GroupByWeight();
    }

    // actualize values as sum of both trees.
    CalculateTotalBucketEntries(&_processingSoftClauseTree, false);
    CalculateTotalBucketEntries(&_softClauseTree, true);

    //        DumpSCNodeStructure( &_processingSoftClauseTree, 2 );

    GroupByColumns();

    assert(_processingSoftClauseTree.empty());
    //        _setting->verbosity = 0;
    break;
  }

  CalculateTotalBucketEntries(&_softClauseTree, false);

  DumpSCNodeStructure(&_softClauseTree, 5);
}

void Cascade::GroupByBiggestRepeatingEntry() {
  if (_processingSoftClauseTree.empty()) {
    return;
  }

  // std::cout << std::endl << __func__ << std::endl;

  // get indices of sorted Bucket entries
  // generate Indice Vector to sort that vector!
  std::vector<uint16_t> sortedBucketIndices(
      _processingSoftClauseTree[0]->highestBucket + 1);
  std::size_t n(0);
  std::generate(std::begin(sortedBucketIndices),
                std::begin(sortedBucketIndices) +
                    _processingSoftClauseTree[0]->highestBucket + 1,
                [&] { return n++; });

  // sort Buckets by total entries.
  std::sort(std::begin(sortedBucketIndices), std::end(sortedBucketIndices),
            [&](std::size_t i1, std::size_t i2) {
              return (_totalBucketEntries[i1] > _totalBucketEntries[i2]);
            });

  uint32_t ind(0);
  while (true) {
    if (_processingSoftClauseTree.size() == 1) {
      _softClauseTree.push_back(_processingSoftClauseTree.back());
      _processingSoftClauseTree.clear();
      std::cout << "c grouping iterations....: " << ind << std::endl;
      break;
    }

    uint16_t maxNodeIndex = GetMaxNodeIndex();

    assert(maxNodeIndex < _processingSoftClauseTree.size());

    // form a new vector with all indices of buckets containing a max Bucket
    // entry
    std::vector<uint16_t> tmpBucketIndices;
    for (auto v : sortedBucketIndices) {
      if (v <= _processingSoftClauseTree[maxNodeIndex]->highestBucket &&
          _processingSoftClauseTree[maxNodeIndex]->occursHowOftenInBucket[v] >
              0)
        tmpBucketIndices.push_back(v);
    }

    if ((ind % 500 == 0 && _setting->verbosity > 2) ||
        (ind % 250 == 0 && _setting->verbosity > 3) || _setting->verbosity > 4)
      DumpMaxNodeOverlappingsAndHeuristicValues(maxNodeIndex,
                                                &tmpBucketIndices);

    // TOBI: Is there another node with exactly the same overlapping? -> then
    // merge more nodes!
    // Maybe there is even a bigger Set in the same subset to merge first with
    // --> see bwt3cc
    int32_t nodeIndexToMergeWith =
        CalculateNodeIndexToMergeWith(maxNodeIndex, &tmpBucketIndices);

    MergeNodes(&tmpBucketIndices, nodeIndexToMergeWith, maxNodeIndex);

    if (_setting->equalWeight > 0 && (ind % _setting->equalWeight == 0))
      GroupByWeight();

    if ((ind % 1000 == 0 && _setting->verbosity > 3) ||
        _setting->verbosity > 4) {
      std::cout << "ProcessingSoftClauseTree" << std::endl;
      CalculateTotalBucketEntries(&_processingSoftClauseTree, false);
      DumpSCNodeStructure(&_processingSoftClauseTree, 6);

      std::cout << "FinalSoftClauseTree" << std::endl;
      CalculateTotalBucketEntries(&_softClauseTree, false);
      DumpSCNodeStructure(&_softClauseTree, 6);

      std::cout << "_processingSoftClauseTree.size(): "
                << _processingSoftClauseTree.size() << std::endl;
    }

    ind++;
  }
  for (std::map<uint32_t, SoftClauseNodes *>::iterator it =
           _processingPercentOffTree.begin();
       it != _processingPercentOffTree.end(); ++it) {
    _softClauseTree.push_back(it->second);
  }
  DumpSCNodeStructure(&_processingSoftClauseTree, 6);
  std::cout << "c percentOff reinsertions: "
            << _howOftenReinsertedFromProcessingPercentOffTree << std::endl;
}

Cascade::BucketOverlaps
Cascade::OverlappingBuckets(std::vector<uint32_t> *bucketIndices,
                            std::vector<uint32_t> *pSCTIndices,
                            bool testOnlyLastBucket) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  BucketOverlaps overlaps;
  if (pSCTIndices == nullptr) {
    pSCTIndices = new std::vector<uint32_t>;
    for (uint32_t i = 0; i < _processingSoftClauseTree.size(); i++) {
      pSCTIndices->push_back(i);
    }
  }
  //    std::cout << "pSCTIndices: " << std::endl;
  //    std::cout << std::setw(4) << "s" << std::setw(8) << "s" << std::setw(15)
  //    << "Weight" << std::setw(3) << "|"; for (uint32_t i = 0; i <=
  //    _numberOfBuckets; ++i)
  //    {
  //        std::cout << std::setw(4) << i;
  //    }
  //    std::cout << std::endl;
  for (auto j : (*pSCTIndices)) {
    if (testOnlyLastBucket) {
      uint32_t k = (*bucketIndices).back();
      if ((_processingSoftClauseTree[j]->weight & 1UL << k)) {
        assert(_processingSoftClauseTree[j]->occursHowOftenInBucket[k] > 0);
        assert(_processingSoftClauseTree[j]->weight >
               static_cast<uint64_t>(pow(2, k)));
        overlaps.noOverlaps++;
        //                std::cout << "weight: " <<
        //                _processingSoftClauseTree[j]->weight << "  PSCTindex:
        //                " << j << " position: " << k << "  pow^position: " <<
        //                std::fixed << static_cast<uint64_t>(pow(2, k)) <<
        //                std::endl;
        //                _processingSoftClauseTree[j]->dumpStructure(true, j);

        overlaps.overlappingSCTIndices.push_back(j);
      }
      continue;
    }
    for (auto i : (*bucketIndices)) {
      if (!(_processingSoftClauseTree[j]->weight & (1UL << i))) {
        break;
      } else if (i == (*bucketIndices).back()) {
        overlaps.noOverlaps++;
        //                std::cout << j << " ";
        overlaps.overlappingSCTIndices.push_back(j);
      }
    }
  }
  //    std::cout << std::endl;

  //    if (overlaps.noOverlaps > 0) {
  //        std::cout << "Compare Buckets: ";
  //        for (auto i : (*bucketIndices)) {
  //            std::cout << i << "  ";
  //        }
  //        std::cout << "OverlSCTInd: ";
  //        for (auto i : overlaps.overlappingSCTIndices) {
  //            std::cout << i << "  ";
  //        }

  //        std::cout << "   - noOverlaps: " << overlaps.noOverlaps <<
  //        std::endl;
  //    }

  overlaps.overlappingBucketIndices = (*bucketIndices);
  return overlaps;
}

uint32_t Cascade::CalcExactTernaryClauseCosts(uint32_t size) {
  if (size == 0)
    return 0;
  //    std::cout << "size: " << size << std::endl;
  uint32_t a = size / 2;
  uint32_t b = size - a;
  uint32_t c = a * b;
  //    std::cout << "a: " << a << std::endl;
  //    std::cout << "b: " << b << std::endl;
  //    std::cout << "c: " << c << std::endl;
  uint32_t d = 0;
  if (c == 0) {
    d = c;
    //        return c;
    //    } else if (a == b) {
    //        return c + 2 * CalcExactTernaryClauseCosts(a);
  } else if (a == 0) {
    d = c + CalcExactTernaryClauseCosts(b);
    //        return c + CalcExactTernaryClauseCosts(b);
  } else if (b == 0) {
    d = c + CalcExactTernaryClauseCosts(a);
    //        return c + CalcExactTernaryClauseCosts(a);
  } else {
    d = c + CalcExactTernaryClauseCosts(a) + CalcExactTernaryClauseCosts(b);
  }
  //    std::cout << "d: " << d << std::endl;
  return d;
  //    return c + CalcExactTernaryClauseCosts(a) +
  //    CalcExactTernaryClauseCosts(b);
}

void Cascade::GroupByColumns() {
  struct GreaterThan {
    bool operator()(const uint32_t &left, const uint32_t &right) const {
      return (left > right);
    }
  };

  std::vector<std::multimap<uint32_t, BucketOverlaps, GreaterThan>>
      sortedOverlapBuckets;
  std::multimap<uint32_t, BucketOverlaps, GreaterThan> sortedOverlaps;
  uint32_t ind(0);

  if (_processingSoftClauseTree.empty()) {
    return;
  }

  // std::cout << std::endl << __func__ << std::endl;

  // get indices of sorted Bucket entries
  // generate Indice Vector to sort that vector!
  //    std::cout << "_numberOfBuckets: " << _numberOfBuckets << std::endl;
  std::vector<uint32_t> sortedBucketIndices(_numberOfBuckets + 1);
  std::size_t n(0);
  std::generate(std::begin(sortedBucketIndices),
                std::begin(sortedBucketIndices) + _numberOfBuckets + 1,
                [&] { return n++; });
  //    std::cout << "sortedBucketIndices.size(): " <<
  //    sortedBucketIndices.size() << std::endl;

  // sort Buckets by total entries.
  std::sort(std::begin(sortedBucketIndices), std::end(sortedBucketIndices),
            [&](std::size_t i1, std::size_t i2) {
              return (_totalBucketEntries[i1] > _totalBucketEntries[i2]);
            });

  std::cout << std::endl;
  //    for (uint32_t j = 0; j < sortedBucketIndices.size(); j++) {
  //        std::cout << "sortedBucketIndices[" << j << "]: " <<
  //        sortedBucketIndices[j] << std::endl;
  //    }
  //    std::cout << std::endl;

  uint32_t totalSavedTernaryClauses = 0;
  uint32_t totalSavedBinaryClauses = 0;
  while (true) {
    ind++;
    // calculate in which buckets the first bucket is in!
    std::vector<uint32_t> tmpBucketIndices;
    for (uint32_t i = 0; i <= _numberOfBuckets; i++) {
      tmpBucketIndices.push_back(i);
      BucketOverlaps bOverlaps = OverlappingBuckets(&tmpBucketIndices);
      sortedOverlaps.insert(
          std::pair<uint32_t, BucketOverlaps>(bOverlaps.noOverlaps, bOverlaps));

      tmpBucketIndices.pop_back();
    }

    sortedOverlapBuckets.clear();
    sortedOverlapBuckets.push_back(sortedOverlaps);

    while (true) {
      //        sleep(1);

      sortedOverlaps.clear();
      uint32_t counter = 0;
      for (auto overlapMM : sortedOverlapBuckets.back()) {
        counter++;
        if (counter > 40) {
          break;
        }
        for (uint32_t bucketIndex = 0; bucketIndex <= _numberOfBuckets;
             bucketIndex++) {
          if (bucketIndex <= overlapMM.second.overlappingBucketIndices.back() ||
              find(overlapMM.second.overlappingBucketIndices.begin(),
                   overlapMM.second.overlappingBucketIndices.end(),
                   bucketIndex) !=
                  overlapMM.second.overlappingBucketIndices.end()) {
            continue;
          }
          counter++;
          tmpBucketIndices = overlapMM.second.overlappingBucketIndices;
          tmpBucketIndices.push_back(bucketIndex);
          BucketOverlaps bOverlaps = OverlappingBuckets(
              &tmpBucketIndices, &overlapMM.second.overlappingSCTIndices, true);

          uint32_t minOverlaps = 3;
          if (_dgpw->_dgpwSetting->featureTest == 3) {
            minOverlaps = 2;
          } else if (_dgpw->_dgpwSetting->featureTest == 4) {
            minOverlaps = 4;
          }

          if (bOverlaps.noOverlaps < minOverlaps) {
            continue;
          }
          sortedOverlaps.insert(std::pair<uint32_t, BucketOverlaps>(
              bOverlaps.noOverlaps, bOverlaps));
        }
      }

      if (sortedOverlaps.empty()) {
        break;
      }
      sortedOverlapBuckets.push_back(sortedOverlaps);
      //        for (std::multimap<uint32_t, BucketOverlaps>::iterator it =
      //        sortedOverlaps.begin(); it != sortedOverlaps.end(); ++it) {
      //            std::cout << "  [" << it->first << ", (";
      //            for (auto j : it->second.overlappingBucketIndices) {
      //                std::cout << j << " ";
      //            }
      //            std::cout << "), (";
      //            for (auto j : it->second.overlappingSCTIndices) {
      //                std::cout << j << " ";
      //            }
      //            std::cout << ")]" << std::endl;
      //        }
      //        if (counter == 8) {
      //            break;
      //        }
    }
    //    std::cout << std::endl << "biggest Values: " << std::endl;
    uint32_t maxTernarySavings = 0;
    uint32_t maxBinarySavings = 0;
    BucketOverlaps chosenBucketOverlaps;
    for (auto overlapMMs : sortedOverlapBuckets) {
      //        std::cout << "Size: " << overlapMMs.size() << std::endl;
      //        std::cout << "  ["  << std::setw(3) << overlapMMs.begin()->first
      //        << ", "; std::cout <<
      //        overlapMMs.begin()->second.overlappingBucketIndices.size() << ":
      //        (" << std::setw(3); for (auto j :
      //        overlapMMs.begin()->second.overlappingBucketIndices) {
      //            std::cout << j << " " << std::setw(3);
      //        }
      //        std::cout << "), " << "(" << std::setw(3);
      //        for (auto j : overlapMMs.begin()->second.overlappingSCTIndices)
      //        {
      //            std::cout << j << " " << std::setw(3);
      //        }
      //        std::cout << ")]" << std::endl;
      uint32_t savedTernaryClauses = 0;
      if (_dgpw->_dgpwSetting->featureTest == 2) {
        savedTernaryClauses =
            ((static_cast<uint32_t>(pow(overlapMMs.begin()->first, 2)) -
              overlapMMs.begin()->first) /
             2) *
            (2 *
             (static_cast<uint32_t>(
                 overlapMMs.begin()->second.overlappingBucketIndices.size() -
                 1)));
      } else if (_dgpw->_dgpwSetting->featureTest == 5) {
        savedTernaryClauses =
            ((static_cast<uint32_t>(pow(overlapMMs.begin()->first, 2)) -
              overlapMMs.begin()->first) /
             2) *
            static_cast<uint32_t>(pow(
                overlapMMs.begin()->second.overlappingBucketIndices.size() - 1,
                2));
      } else if (_dgpw->_dgpwSetting->featureTest == 6) {
        savedTernaryClauses =
            ((static_cast<uint32_t>(pow(overlapMMs.begin()->first, 2)) -
              overlapMMs.begin()->first) /
             2) *
            static_cast<uint32_t>(pow(
                overlapMMs.begin()->second.overlappingBucketIndices.size() - 1,
                3));
      } else {
        savedTernaryClauses =
            (static_cast<uint32_t>(pow(overlapMMs.begin()->first, 2)) -
             overlapMMs.begin()->first) /
            2 *
            static_cast<uint32_t>(
                overlapMMs.begin()->second.overlappingBucketIndices.size() - 1);
      }
      uint32_t savedBinaryClauses =
          (static_cast<uint32_t>(log2(overlapMMs.begin()->first)) *
           overlapMMs.begin()->first *
           static_cast<uint32_t>(
               overlapMMs.begin()->second.overlappingBucketIndices.size() - 1));
      //        uint32_t savedVariables = savedBinaryClauses;
      if (savedTernaryClauses > maxTernarySavings) {
        //            std::cout << "This bucket is chosen" << std::endl;
        //            for (auto j :
        //            overlapMMs.begin()->second.overlappingBucketIndices) {
        ////                std::cout << j << " " << std::setw(3);
        //            }
        maxTernarySavings = savedTernaryClauses;
        maxBinarySavings = savedBinaryClauses;
        chosenBucketOverlaps = overlapMMs.begin()->second;
      }
      //        std::cout << "Saved ternary clauses: " << std::setw(5) <<
      //        savedTernaryClauses; std::cout << "  Saved binary clauses: " <<
      //        std::setw(5) << savedBinaryClauses << std::endl; std::cout << "
      //        Saved variables: " << std::setw(5) << savedVariables <<
      //        std::endl << std::endl; savedTernaryClauses = 0;
      //        savedTernaryClauses =
      //        CalcExactTernaryClauseCosts(overlapMMs.begin()->first) *
      //        static_cast<uint32_t>(overlapMMs.begin()->second.overlappingBucketIndices.size()
      //        - 1); savedBinaryClauses =
      //        (static_cast<uint32_t>(log2(overlapMMs.begin()->first)) *
      //        overlapMMs.begin()->first *
      //        static_cast<uint32_t>(overlapMMs.begin()->second.overlappingBucketIndices.size()
      //        - 1)); savedVariables = savedBinaryClauses; std::cout << "Saved
      //        ternary clauses: " << std::setw(5) << savedTernaryClauses;
      //        std::cout << "   Saved binary clauses: " << std::setw(5) <<
      //        savedBinaryClauses; std::cout << "   Saved variables: " <<
      //        std::setw(5) << savedVariables << std::endl << std::endl;
    }
    totalSavedTernaryClauses += maxTernarySavings;
    totalSavedBinaryClauses += maxBinarySavings;
    //    std::cout << "maxSavings: " << maxTernarySavings << std::endl;
    if (maxTernarySavings == 0) {
      break;
    }

    if (_processingSoftClauseTree.size() == 1) {
      _softClauseTree.push_back(_processingSoftClauseTree.back());
      _processingSoftClauseTree.clear();
      break;
    }

    //    std::cout << "Chosen Bucket overlaps!!!" << std::endl;
    //    std::cout << "  ["  << std::setw(3) << chosenBucketOverlaps.noOverlaps
    //    << ", "; std::cout <<
    //    chosenBucketOverlaps.overlappingBucketIndices.size() << ": (" <<
    //    std::setw(3); for (auto j :
    //    chosenBucketOverlaps.overlappingBucketIndices) {
    //        std::cout << j << " " << std::setw(3);
    //    }
    //    std::cout << "), " << "(" << std::setw(3);
    //    for (auto j : chosenBucketOverlaps.overlappingSCTIndices) {
    //        std::cout << j << " " << std::setw(3);
    //    }
    //    std::cout << ")]" << std::endl;

    MergeMultipleNodes(&chosenBucketOverlaps);

    //    if ((ind % 1000 == 0 && _setting->verbosity > 3)  ||
    //    _setting->verbosity > 4)
    //    {
    //        std::cout << "ProcessingSoftClauseTree" << std::endl;
    //        CalculateTotalBucketEntries( &_processingSoftClauseTree, false);
    //        DumpSCNodeStructure(&_processingSoftClauseTree, 3);

    //        std::cout << "FinalSoftClauseTree" << std::endl;
    //        CalculateTotalBucketEntries( &_softClauseTree, false);
    //        DumpSCNodeStructure(&_softClauseTree, 3);
    //    }

    //    if (ind == 7) {
    //        break;
    //    }
  }

  std::cout << "c adder caching on all rows of the soft clause tree in " << ind
            << " rounds!" << std::endl;

  // append all elements of the _processingSoftClauseTree to _softClauseTree
  _softClauseTree.insert(std::end(_softClauseTree),
                         std::begin(_processingSoftClauseTree),
                         std::end(_processingSoftClauseTree));
  _processingSoftClauseTree.clear();

  if ((ind % 1000 == 0 && _setting->verbosity > 3) || _setting->verbosity > 4) {
    std::cout << "FinalSoftClauseTree" << std::endl;
    DumpSCNodeStructure(&_softClauseTree, 6);
  }

  std::cout << "c TotalSavedTernaryClause: " << totalSavedTernaryClauses
            << std::endl;
  std::cout << "c TotalSavedBinaryClauses: " << totalSavedBinaryClauses
            << std::endl;
  std::cout << "c grouping iterations....: " << ind << std::endl;
  //    exit(0);
}

void Cascade::MergeMultipleNodes(BucketOverlaps *overlaps) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  assert(overlaps->overlappingBucketIndices.size() > 1);
  assert(overlaps->overlappingSCTIndices.size() > 1);

  // calc new weight for the two subweights to merge
  uint64_t newOverlappingWeight(0);
  for (auto v : overlaps->overlappingBucketIndices) {
    newOverlappingWeight += pow(_base, v);
    //        std::cout << "v: " << v << "   pow(base, v): " << std::fixed <<
    //        static_cast<uint64_t>(pow(_base,v)) << "
    //        newOverlappingWeight: " << newOverlappingWeight << std::endl;
  }
  //    std::cout << "newOverlappingWeight: " << newOverlappingWeight <<
  //    std::endl;
  Bucket *tmpBucket = new Bucket(_dgpw, this);
  for (auto w : overlaps->overlappingSCTIndices) {
    tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[w]);
  }
  SoftClauseNodes *sCNode =
      new SoftClauseNodes(tmpBucket, newOverlappingWeight, _base);
  //    sCNode->dumpStructure(true);
  // is at least in two buckets - otherwise no merge
  _processingSoftClauseTree.push_back(sCNode);

  std::vector<uint32_t> eraseList;

  // erase the old SoftClauseNodes.
  for (auto w : overlaps->overlappingSCTIndices) {
    //        std::cout << "_processingSoftClauseTree[w]->weight: " <<
    //        _processingSoftClauseTree[w]->weight << "  index: " << w <<
    //        std::endl; std::cout << "newOverlappingWeight: " <<
    //        newOverlappingWeight << "  index: " << w << std::endl; std::cout
    //        << "newWeight: " << newWeight << "  index: " << w << std::endl <<
    //        std::endl; std::cout << "pSCT[w]->weight: " << std::setw(20) <<
    //        _processingSoftClauseTree[w]->weight; std::cout << "
    //        newOverlappingWeight: "  << std::setw(20) << newOverlappingWeight
    //        << "  index: "  << std::setw(5) << w << std::endl;
    assert(_processingSoftClauseTree[w]->weight >= newOverlappingWeight);
    uint64_t newWeight =
        _processingSoftClauseTree[w]->weight - newOverlappingWeight;
    //        std::cout << "   newWeight: "  << std::setw(20) << newWeight <<
    //        std::endl;
    if (newWeight == 0) {
      eraseList.push_back(w);
    } else {
      _processingSoftClauseTree[w]->setWeight(newWeight);
      if (_processingSoftClauseTree[w]->inHowManyBuckets == 1) {
        _softClauseTree.push_back(_processingSoftClauseTree[w]);
        eraseList.push_back(w);
      }
    }
  }

  for (uint32_t eraseInd = static_cast<uint32_t>(eraseList.size());
       eraseInd > 0; eraseInd--) {
    //        std::cout << "EraseIndex: " << eraseList[eraseInd - 1] <<
    //        std::endl; _processingSoftClauseTree[eraseList[eraseInd -
    //        1]]->dumpStructure(true);
    _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                    eraseList[eraseInd - 1]);
  }

  if (_setting->verbosity < 4)
    return;

  //    std::cout << std::setw(30) << "new Weight: " << newOverlappingWeight <<
  //    std::endl;
}

void Cascade::GroupByWeight() {
  std::vector<SoftClauseNodes *> tmpSCTree;
  // sort Buckets by total entries.
  // use stable sort only if there are identical weights!
  std::stable_sort(std::begin(_processingSoftClauseTree),
                   std::end(_processingSoftClauseTree),
                   [&](SoftClauseNodes *sCN1, SoftClauseNodes *sCN2) {
                     return (sCN1->weight > sCN2->weight);
                   });

  for (uint32_t i = 0; i < _processingSoftClauseTree.size(); ++i) {
    if (i + 1 == _processingSoftClauseTree.size() ||
        _processingSoftClauseTree[i + 1]->weight !=
            _processingSoftClauseTree[i]->weight) {
      tmpSCTree.push_back(_processingSoftClauseTree[i]);
      continue;
    }
    Bucket *tmpBucket = new Bucket(_dgpw, this);
    while (i + 1 < _processingSoftClauseTree.size() &&
           _processingSoftClauseTree[i + 1]->weight ==
               _processingSoftClauseTree[i]->weight) {
      tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[i]);
      i++;
    }
    tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[i]);
    SoftClauseNodes *sCNode = new SoftClauseNodes(
        tmpBucket, _processingSoftClauseTree[i]->weight, _base);
    tmpSCTree.push_back(sCNode);
  }
  _processingSoftClauseTree = tmpSCTree;
}

uint16_t Cascade::GetMaxNodeIndex() {
  // generate Indice Vector to get highest index from!
  std::vector<uint16_t> maxSCTIndices(_processingSoftClauseTree.size());
  std::size_t m(0);
  std::generate(std::begin(maxSCTIndices),
                std::begin(maxSCTIndices) + _processingSoftClauseTree.size(),
                [&] { return m++; });

  // Get index of SoftClauseTree with most Bucket entries!
  // If 2 SCN have same # entries, take the one with more occurrences.
  return *std::max_element(
      std::begin(maxSCTIndices), std::end(maxSCTIndices),
      [&](std::size_t i1, std::size_t i2) {
        return ((_processingSoftClauseTree[i1]->inHowManyBuckets <
                 _processingSoftClauseTree[i2]->inHowManyBuckets) ||
                ((_processingSoftClauseTree[i1]->inHowManyBuckets ==
                  _processingSoftClauseTree[i2]->inHowManyBuckets) &&
                 (_processingSoftClauseTree[i1]->GetOccurrences() *
                      _processingSoftClauseTree[i1]->size() <
                  _processingSoftClauseTree[i2]->GetOccurrences() *
                      _processingSoftClauseTree[i2]->size())));
      });

  // First Idea - sort whole tree...
  // sort whole SoftClauseTree by entries and if they are equal by occurences.
  // maybe to look later on only to the top elements. Otherwise it is sufficient
  // to find in O(n) the biggest element.
  // std::sort(_processingSoftClauseTree.begin(),
  // _processingSoftClauseTree.end(),
  // SoftClauseNodes::SortByEntriesOccurrences);
}

void Cascade::DumpMaxNodeOverlappingsAndHeuristicValues(
    uint16_t maxNodeIndex, std::vector<uint16_t> *tmpBucketIndices) {
  if (_setting->verbosity < 1)
    return;

  std::cout << "_processingSoftClauseTree.size: "
            << _processingSoftClauseTree.size() << std::endl;
  std::cout << "maxBucketEntryIndex: " << maxNodeIndex << "  maxBucketEntries: "
            << _processingSoftClauseTree[maxNodeIndex]->inHowManyBuckets
            << std::endl
            << std::endl;
  for (auto v : *tmpBucketIndices)
    std::cout << std::setw(4) << v;
  std::cout << std::endl;
  for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
    std::cout << std::setw(4) << "----";
  std::cout << "--------------------------------------------------"
            << std::endl;

  for (auto v : *tmpBucketIndices) {
    if (_processingSoftClauseTree[maxNodeIndex]->occursHowOftenInBucket[v] >
            0 &&
        v < _processingSoftClauseTree[maxNodeIndex]->highestBucket + 1)
      std::cout << std::setw(4)
                << _processingSoftClauseTree[maxNodeIndex]->size();
    else
      std::cout << std::setw(4) << "";
  }
  std::cout << "   |" << std::setw(7) << "0"
            << " |" << std::setw(7) << "1"
            << " |" << std::setw(7) << "2"
            << " |" << std::setw(7) << "3"
            << " |" << std::setw(7) << "4"
            << " |" << std::setw(7) << "5"
            << " |" << std::setw(7) << "6"
            << " |  <-- heuristics" << std::endl;

  for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
    std::cout << std::setw(4) << "----";
  std::cout << "--------------------------------------------------"
            << std::endl;

  int32_t maxCosts(0);
  int32_t estimatedMaxCostIndex(0);

  for (int32_t j = 0; j < (int)_processingSoftClauseTree.size(); j++) {
    if (j == maxNodeIndex)
      continue;

    std::vector<SoftClauseNodes *> NodesToMerge = {
        _processingSoftClauseTree[maxNodeIndex], _processingSoftClauseTree[j]};
    int32_t usedCosts = GetBenefitOfMergingNodes(NodesToMerge, tmpBucketIndices,
                                                 true, _setting->groupHeuristic,
                                                 _setting->percentOff);

    if (usedCosts == 0)
      continue;

    if (maxCosts < usedCosts) {
      maxCosts = usedCosts;
      estimatedMaxCostIndex = j;
    }
  }

  for (uint64_t ind = 0; ind < (*tmpBucketIndices).size(); ind++)
    std::cout << std::setw(4) << "----";
  std::cout << "--------------------------------------------------"
            << std::endl;

  for (auto v : *tmpBucketIndices)
    std::cout << std::setw(4) << _totalBucketEntries[v];
  std::cout << std::endl;
  std::cout << std::endl
            << "MaxCosts: " << maxCosts << "  Index: " << estimatedMaxCostIndex
            << std::endl;
}

void Cascade::DumpModelOfTares(uint16_t verbosity) {
  if (_setting->verbosity < verbosity || _dgpw->GetLastResult() != 10)
    return;

  std::cout << __PRETTY_FUNCTION__ << std::endl;

  // if there wasn't a SAT call making the result bigger after constructing the
  // whole cascade, then the tares are 0.
  //    if (_dgpw->_lastModel[_structure[0]->_tares[0]] == 0)
  //        return;
  // output the variables of the Tare T
  //    _dgpw->Solve();
  std::cout << "Model of Tares: (n...0): ";
  for (int bucketInd = (_structure.size() - 1); bucketInd >= 0; --bucketInd) {
    std::cout << "(";
    for (int tareInd =
             static_cast<int>(_structure[bucketInd]->_tares.size() - 1);
         tareInd >= 0; --tareInd) {
      std::cout << _dgpw->Model(_structure[bucketInd]->_tares[tareInd]);
      if (tareInd != 0) {
        std::cout << ", ";
      }
    }
    if (bucketInd != 0) {
      std::cout << "), ";
    } else {
      std::cout << ")";
    }
  }
  std::cout << std::endl;
}

void Cascade::DumpBucketSolveInformation(uint32_t actualPos, bool _isLastBucket,
                                         uint16_t verbosity) {
  if (_setting->verbosity < verbosity)
    return;

  if (_isLastBucket) {
    _estimatedWeightBoundaries[0] = _highestBucketMultiplicator * actualPos;
    _estimatedWeightBoundaries[1] =
        _highestBucketMultiplicator * _structure.back()->size();
    std::cout << std::setw(52) << "weight boundaries: ( "
              << _estimatedWeightBoundaries[0] << " / "
              << _estimatedWeightBoundaries[1] << " )" << std::endl;
  }
  std::cout << std::setw(52) << "satisfied weights: ( " << _satWeight << " / "
            << _sumOfSoftWeights << " )" << std::endl;
  std::cout << std::setw(50) << "actualPosition: " << actualPos << std::endl;
  std::cout << "---------------------------------------------------------------"
               "------------------"
               "----------"
            << std::endl;
}

int32_t Cascade::CalculateNodeIndexToMergeWith(
    uint16_t maxNodeIndex, std::vector<uint16_t> *tmpBucketIndices) {
  int32_t maxCosts(0);
  // highest possible number to indicate that there is no clause node to merge
  // with!
  int32_t estimatedMaxCostIndex(-1);

  int32_t sumOfSizesOfPercentageFails(0);
  for (uint32_t j = 0; j < _processingSoftClauseTree.size(); j++) {
    if (j == maxNodeIndex)
      continue;

    std::vector<SoftClauseNodes *> NodesToMerge = {
        _processingSoftClauseTree[maxNodeIndex], _processingSoftClauseTree[j]};
    int32_t usedCosts = GetBenefitOfMergingNodes(
        NodesToMerge, tmpBucketIndices, false, _setting->groupHeuristic,
        _setting->percentOff);
    if (usedCosts == 0) {
      continue;
    } else if (usedCosts < 0) {
      sumOfSizesOfPercentageFails -= usedCosts;
    }

    if (maxCosts < usedCosts) {
      maxCosts = usedCosts;
      estimatedMaxCostIndex = j;
    }
  }
  // If no merge is possible
  // check if there is the possibility that the node can be merged later on
  if (estimatedMaxCostIndex == -1 && sumOfSizesOfPercentageFails > 0 &&
      _setting->percentOffReinsert) {
    if ((double)_processingSoftClauseTree[maxNodeIndex]->size() *
            (double)((100 - _setting->percentOff) / 100) <=
        (double)sumOfSizesOfPercentageFails) {
      return -2;
    }
  }

  // std::cout << "sumOfSizesOfPercentageFails: " << sumOfSizesOfPercentageFails
  // << std::endl;
  return estimatedMaxCostIndex;
}

int32_t
Cascade::GetBenefitOfMergingNodes(std::vector<SoftClauseNodes *> NodesToMerge,
                                  std::vector<uint16_t> *tmpBucketIndices,
                                  bool dump, uint32_t heuristic,
                                  int32_t minPercentOff) {
  int32_t maxSize = (NodesToMerge[0]->size() > NodesToMerge[1]->size())
                        ? NodesToMerge[0]->size()
                        : NodesToMerge[1]->size();
  int32_t minSize = (NodesToMerge[0]->size() > NodesToMerge[1]->size())
                        ? NodesToMerge[1]->size()
                        : NodesToMerge[0]->size();
  //    std::cout << "maxSize: " << maxSize << "  minSize: " << minSize <<
  //    std::endl;
  int calcPercentOff = 100 - (int)(((double)minSize / (double)maxSize) * 100);
  int32_t binaryClauses = NodesToMerge[1]->size() + NodesToMerge[0]->size();
  int32_t ternaryClauses = NodesToMerge[1]->size() * NodesToMerge[0]->size();
  int32_t clauseCosts = 2 * ternaryClauses + binaryClauses;
  int32_t usedCosts(0);
  int32_t howManyBucketsOverlapping(0);
  int32_t sumOfBucketSizes(1);
  int32_t howOftenBucketsOverlapping(0);
  int32_t bucketSizesFactor(1);
  for (auto v : *tmpBucketIndices) {
    if (!(NodesToMerge[1]->occursHowOftenInBucket[v] > 0 &&
          v < NodesToMerge[1]->highestBucket + 1))
      continue;

    bucketSizesFactor +=
        static_cast<uint32_t>(pow(_totalBucketEntries[v], 1.7) * 0.1);

    howManyBucketsOverlapping += NodesToMerge[1]->size();
    sumOfBucketSizes += _totalBucketEntries[v];
    howOftenBucketsOverlapping++;
  }
  if (howOftenBucketsOverlapping <= 1)
    return 0;

  switch (heuristic) {
  // Standard combination of other heuristics
  case 0:
    usedCosts =
        bucketSizesFactor *
        static_cast<int32_t>(pow((howOftenBucketsOverlapping - 1), 1.7) *
                             pow((clauseCosts), 0.5));
    break;
  // number of reduced (ternary clauses * 2 + binary clauses)
  case 1:
    usedCosts = (howOftenBucketsOverlapping - 1) * clauseCosts;
    //        std::cout << "howOftenBucketsOverlapping: " <<
    //        howOftenBucketsOverlapping << "   clauseCosts: " << clauseCosts
    //        << "   ternaryClauses : " << ternaryClauses << "  binaryClauses:
    //        " << binaryClauses << std::endl;
    break;
  // sum of bucket sizes the SC occurs
  case 2:
    usedCosts = sumOfBucketSizes;
    break;
  // closest size in percentage
  case 3:
    usedCosts = 100 - calcPercentOff;
    break;
  // the one with least occurences in other buckets
  case 4:
    usedCosts =
        (int)_totalBucketEntries.size() -
        ((int)NodesToMerge[1]->inHowManyBuckets - howOftenBucketsOverlapping);
    break;
  // how many merges are possible after this merge
  case 5:
    usedCosts =
        CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 0);
    break;
  // greatest depth of submerges (till depth 3, then adds possible merges of
  // depth 4) can be made more effective - but for the ones with big numbers
  // it is way to much work!
  case 6:
    usedCosts =
        CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);
    break;
  //
  case 7:
    // usedCosts = CalculateNumberOfPossibleSubmerges(tmpBucketIndices,
    // NodesToMerge[1]);
    break;
  }

  if (minPercentOff != 100) {
    // at least one difference - should be always possible, or the calculated
    // percentage.
    if (maxSize - minSize != 1 && calcPercentOff > minPercentOff) {
      // if NodesToMerge[0]->size() is bigger, then give back the negative value
      // to calculate if it is theoretically possible to reach that value again
      // if no merge at all is possible.
      if (NodesToMerge[0]->size() > NodesToMerge[1]->size()) {
        usedCosts = -NodesToMerge[1]->size();
      } else {
        usedCosts = 0;
      }
    }
  }

  if (!dump)
    return usedCosts;

  int32_t heuristic0 =
      bucketSizesFactor *
      static_cast<int32_t>(pow((howOftenBucketsOverlapping - 1), 1.7) *
                           pow((clauseCosts), 0.5));
  int32_t heuristic1 = (howOftenBucketsOverlapping - 1) * clauseCosts;
  //            CalculateHowManySCNHaveTheSameOverlap(tmpBucketIndices,
  //            NodesToMerge);
  int32_t heuristic2 = sumOfBucketSizes;
  int32_t heuristic3 = calcPercentOff;
  int32_t heuristic4 =
      (int)_totalBucketEntries.size() -
      ((int)NodesToMerge[1]->inHowManyBuckets - howOftenBucketsOverlapping);
  ;
  int32_t heuristic5 =
      CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 0);
  int32_t heuristic6 =
      CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);
  int32_t heuristic7 =
      CalculateHowManySCNHaveTheSameOverlap(tmpBucketIndices, NodesToMerge);
  //    int32_t heuristic7 =
  //    CalculateNumberOfPossibleSubmerges(tmpBucketIndices, NodesToMerge, 1);

  for (auto v : *tmpBucketIndices) {
    if (NodesToMerge[1]->occursHowOftenInBucket[v] > 0 &&
        v < NodesToMerge[1]->highestBucket + 1)
      std::cout << std::setw(4) << NodesToMerge[1]->size();
    else
      std::cout << std::setw(4) << "";
  }
  std::cout << "  ";
  std::cout << " |" << std::setw(7) << heuristic0;
  std::cout << " |" << std::setw(7) << heuristic1;
  std::cout << " |" << std::setw(7) << heuristic2;
  std::cout << " |" << std::setw(7) << heuristic3;
  std::cout << " |" << std::setw(7) << heuristic4;
  std::cout << " |" << std::setw(7) << heuristic5;
  std::cout << " |" << std::setw(7) << heuristic6;
  std::cout << " |" << std::setw(7) << heuristic7;

  if (minPercentOff != 100) {
    std::cout << " |" << std::setw(7) << calcPercentOff;

    if (usedCosts < 0)
      std::cout << "  <   " << minPercentOff << " <- too big difference";
    else
      std::cout << "  >=  " << minPercentOff;

    std::cout << " |  usedCosts: " << usedCosts << std::endl;
  } else {
    std::cout << std::endl;
  }

  return usedCosts;
}

int32_t Cascade::CalculateHowManySCNHaveTheSameOverlap(
    std::vector<uint16_t> *tmpBucketIndices,
    std::vector<SoftClauseNodes *> NodesToMerge) {
  //    std::cout << std::endl << __PRETTY_FUNCTION__ << std::endl;

  std::vector<uint16_t> overlappingIndices;
  //    std::cout << "Overlapping indices are: ";
  for (auto v : *tmpBucketIndices) {
    if (NodesToMerge[1]->occursHowOftenInBucket[v] >= 1) {
      //            std::cout << v << " ";
      overlappingIndices.push_back(v);
    }
  }
  //    std::cout << std::endl;
  //    std::cout << "They have " << overlappingIndices.size() << " overlapping
  //    indices." << std::endl;

  std::vector<SoftClauseNodes *> overlappingSCN;
  for (auto softClauseNodes : _processingSoftClauseTree) {
    size_t overlaps = 0;
    //        std::cout << std::endl;
    for (auto v : overlappingIndices) {
      if (softClauseNodes->occursHowOftenInBucket[v] >= 1 &&
          softClauseNodes->highestBucket >= v) {
        overlaps++;
        //                std::cout << "v: " << v << "  howOften: " <<
        //                softClauseNodes->occursHowOftenInBucket[v] << "
        //                overlaps: " << overlaps << std::endl;
      } else {
        break;
      }
    }
    if (overlaps == overlappingIndices.size()) {
      overlappingSCN.push_back(softClauseNodes);
      //            softClauseNodes->dumpStructure(0);
    }
  }
  //    std::cout << "Overlapping SoftClauseNodes: " << overlappingSCN.size() <<
  //    std::endl;

  return static_cast<int32_t>(overlappingSCN.size());

  //    std::cout << std::endl;
  for (auto w : NodesToMerge) {
    std::cout << "DumpStructure: " << std::endl;
    w->dumpStructure(1);
    std::cout << "Done!" << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
  exit(1);
}

int32_t Cascade::CalculateNumberOfPossibleSubmerges(
    std::vector<uint16_t> *tmpBucketIndices,
    std::vector<SoftClauseNodes *> NodesToMerge, int32_t depth) {
  int32_t incomingDepth(depth);
  int32_t possibleSubmerges(0);
  int32_t highestDepth(0);
  std::vector<uint16_t> lastNodeIndices;
  // std::cout << "SNI: ";

  for (auto v : *tmpBucketIndices) {
    if (!(NodesToMerge.back()->occursHowOftenInBucket[v] > 0 &&
          v < NodesToMerge.back()->highestBucket + 1))
      continue;
    lastNodeIndices.push_back(v);
    // std::cout << v << ", ";
  }

  for (uint32_t j = 0; j < _processingSoftClauseTree.size(); j++) {
    // checking if j'th element is one of the yet merged ones.
    if (std::find(NodesToMerge.begin(), NodesToMerge.end(),
                  _processingSoftClauseTree[j]) != NodesToMerge.end()) {
      // std::cout << "sizeofNodesToMerge: " << NodesToMerge.size() << "  Index:
      // " << j << std::endl;
      continue;
    }
    // std::cout << "j: " << j << std::endl;
    // if (_processingSoftClauseTree[j] == NodesToMerge.back())
    //    continue;
    int32_t occurencesInNode(0);
    for (auto v : lastNodeIndices) {
      if (_processingSoftClauseTree[j]->highestBucket < v)
        continue;
      // std::cout << "j: " << j << " v: " << v << std::endl;
      if (_processingSoftClauseTree[j]->occursHowOftenInBucket[v] > 0)
        occurencesInNode++;
    }
    if (depth > 0 && occurencesInNode > 1) {
      NodesToMerge.push_back(_processingSoftClauseTree[j]);
      // here we are at depth 3! - go one depth deeper and Calculate then the
      // number of possible submerges!
      if (incomingDepth == 2)
        return (incomingDepth + 1 +
                CalculateNumberOfPossibleSubmerges(&lastNodeIndices,
                                                   NodesToMerge, 0));

      int32_t newDepth = CalculateNumberOfPossibleSubmerges(
          &lastNodeIndices, NodesToMerge, incomingDepth + 1);
      // std::cout << "newDepth: " << newDepth << std::endl;
      if (newDepth > highestDepth)
        highestDepth = newDepth;
    } else if (occurencesInNode > 1) {
      // std::cout << "Index: " << j << "  oIN: " << occurencesInNode <<
      // std::endl;
      possibleSubmerges++;
    }
  }
  if (highestDepth == 0)
    highestDepth = incomingDepth;
  if (depth > 0)
    return highestDepth;
  else
    return possibleSubmerges;
}

void Cascade::MergeNodes(std::vector<uint16_t> *tmpBucketIndices,
                         int32_t nodeIndexToMergeWith, int32_t maxNodeIndex) {
  // Exception - there is no node to merge with!!
  if (nodeIndexToMergeWith == -1) {
    _softClauseTree.push_back(_processingSoftClauseTree[maxNodeIndex]);
    _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                    maxNodeIndex);
    return;
  }

  // Exception - no node to merge with, but with other merges it can become
  // possible again.
  else if (nodeIndexToMergeWith == -2) {
    _processingPercentOffTree.insert(
        std::make_pair(_processingSoftClauseTree[maxNodeIndex]->size(),
                       _processingSoftClauseTree[maxNodeIndex]));
    _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                    maxNodeIndex);
    return;
  }

  // calc new weight for the two subweights to merge
  uint64_t newWeightforBoth(0);
  for (auto v : *tmpBucketIndices) {
    if (_processingSoftClauseTree[nodeIndexToMergeWith]
            ->occursHowOftenInBucket[v] &&
        v < _processingSoftClauseTree[nodeIndexToMergeWith]->highestBucket +
                1) {
      newWeightforBoth += pow(_base, v);
    }
  }

  Bucket *tmpBucket = new Bucket(_dgpw, this);
  tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[maxNodeIndex]);
  tmpBucket->AddSoftClauseNode(_processingSoftClauseTree[nodeIndexToMergeWith]);
  SoftClauseNodes *sCNode =
      new SoftClauseNodes(tmpBucket, newWeightforBoth, _base);
  // at least in two buckets - otherwise no merge
  _processingSoftClauseTree.push_back(sCNode);

  // erase the old SoftClauseNodes.
  uint16_t nodeIndex =
      maxNodeIndex > nodeIndexToMergeWith ? maxNodeIndex : nodeIndexToMergeWith;

  for (uint16_t i = 0; i < 2; i++) {
    if (i == 1)
      nodeIndex = maxNodeIndex > nodeIndexToMergeWith ? nodeIndexToMergeWith
                                                      : maxNodeIndex;

    if (_processingSoftClauseTree[nodeIndex]->weight - newWeightforBoth == 0) {
      _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                      nodeIndex);
    } else {
      _processingSoftClauseTree[nodeIndex]->setWeight(
          _processingSoftClauseTree[nodeIndex]->weight - newWeightforBoth);
      // only in one bucket - then move it to the from processingSCT to the SCT
      if (_processingSoftClauseTree[nodeIndex]->inHowManyBuckets == 1) {
        _softClauseTree.push_back(_processingSoftClauseTree[nodeIndex]);
        _processingSoftClauseTree.erase(_processingSoftClauseTree.begin() +
                                        nodeIndex);
      }
    }
  }

  // if the new weight is bigger than a weight from the _laterToMergeMultimap
  // move the corresponding SCN to the _processing SCT.
  for (std::multimap<uint32_t, SoftClauseNodes *>::iterator it =
           _processingPercentOffTree.begin();
       it != _processingPercentOffTree.end(); ++it) {
    if (!AtLeastTwoBucketsInCommon(it->second,
                                   _processingSoftClauseTree.back()))
      continue;

    int32_t maxSize = (it->first > _processingSoftClauseTree.back()->size())
                          ? it->first
                          : _processingSoftClauseTree.back()->size();
    int32_t minSize = (it->first > _processingSoftClauseTree.back()->size())
                          ? _processingSoftClauseTree.back()->size()
                          : it->first;
    uint32_t calcPercentOff =
        100 - (uint32_t)(((double)minSize / (double)maxSize) * 100);

    // Is there a chance of merging again?
    if (calcPercentOff <= _setting->percentOff) {
      if (_setting->verbosity > 3) {
        std::cout << "_processingPercentOffTree.size(); "
                  << _processingPercentOffTree.size() << std::endl;
        std::cout << "percentOff: " << calcPercentOff
                  << "   _setting->_percentOff: " << _setting->percentOff
                  << "   minSize: " << minSize << "  maxSize: " << maxSize
                  << std::endl;
      }

      _processingSoftClauseTree.push_back(it->second);
      _processingPercentOffTree.erase(it);
      _howOftenReinsertedFromProcessingPercentOffTree++;
    }
    // because of sorted Multimap - there is no chance that calcPercentOff gets
    // biggerEqual than _dgpw->groupPercentOff again.
    else if (it->first > _processingSoftClauseTree.back()->size()) {
      break;
    }
  }

  if (_setting->verbosity < 4)
    return;

  std::cout << std::setw(30) << "new Weight: " << newWeightforBoth << std::endl;
}

void Cascade::FillBuckets() {
  // create the bucket structure.
  for (uint16_t ind = 0; ind <= _numberOfBuckets; ind++) {
    Bucket *tmpBucket = new Bucket(_dgpw, this, ind);
    _structure.push_back(tmpBucket);
  }

  // fill each Bucket according to the _softClauseTree structure
  for (auto sCNode : _softClauseTree) {
    for (uint16_t ind = 0; ind <= sCNode->highestBucket; ind++) {
      //            if (sCNode->occursHowOftenInBucket[ind] == 0)
      //                continue;
      //
      for (uint32_t howOften = 0;
           howOften < sCNode->occursHowOftenInBucket[ind]; howOften++) {
        _structure[ind]->AddSoftClauseNode(sCNode);
      }
    }
  }

  _highestBucketMultiplicator =
      static_cast<uint64_t>(pow(_base, _numberOfBuckets));

  if (_setting->verbosity < 4)
    return;

  std::cout << "Buckets are filled with SoftClauseNodes!" << std::endl;
  DumpBucketStructure(false, 5);
}

void Cascade::AddTaresToBuckets() {
  assert(_structure.size() > 0);

  //    uint32_t addNoTareToLastBucket = ( _dgpw->_cascadeDivider > 0 ||
  //    _onlyByTares ) ? 0 : 1;
  uint32_t addNoTareToLastBucket = 1;

  if (_setting->verbosity > 3) {
    std::cout << "caddNoTareToLastBucket: " << addNoTareToLastBucket
              << std::endl;
    std::cout << "cTare(s): ";
  }
  // add at the last position _base many variables T to each trigger except the
  // last one.
  for (uint32_t i = 0; i < _structure.size() - addNoTareToLastBucket; i++) {
    // assert( _dgpw->_sorterTree[i].size() <= 1 );
    if (!_structure.empty()) {
      for (uint32_t j = 0; j < _setting->base - 1; j++)
        AddTare(i);
    }
  }
  vPL->write_comment("All tares are added!");

  if (_setting->verbosity < 2)
    return;

  std::cout << std::endl << "c Tares are added to Structure!" << std::endl;
}

void Cascade::EncodeTopBuckets() {
  TimeMeasurement timeEncoding(&_dgpw->_timeVariables->encoding, true);
  for (auto bucket : _structure) {
    bucket->CalculateNumberOfClauses(true, true, false);

    if (_setting->partitionStrategy == GROUPBYWEIGHTADDATLAST) {
      bucket->EncodeTopAddAtLast();
    } else {
      bucket->EncodeTop();
    }

    bucket->CalculateNumberOfClauses(true, false, true);
  }

  if (_setting->verbosity < 1)
    return;
  if (_setting->verbosity > 3)
    std::cout << std::endl << "Top Buckets are encoded!" << std::endl;
}

void Cascade::EncodeBottomBuckets() {
  TimeMeasurement timeEncoding(&_dgpw->_timeVariables->encoding, true);
  for (uint16_t ind = 1; ind < _structure.size(); ind++) {
    _structure[ind]->CalculateNumberOfClauses(false, true, false);

    _structure[ind]->MergeSorterWith(
        _structure[ind - 1]->GetEveryNthOutput(_base));

    _structure[ind]->CalculateNumberOfClauses(false, false, true);
  }

  _structure.back()->_isLastBucket = true;

  if (_setting->verbosity < 1)
    return;

  DumpNumberOfBucketsAndClauses();

  if (_setting->verbosity > 3)
    std::cout << "Bottom Buckets are encoded!" << std::endl << std::endl;
}

void Cascade::UnionBucketsIntoLast() {
  // Push all Buckets to one final Bucket == _structure.back()
  for (uint16_t ind = 1; ind < _structure.size(); ind++) {
    _structure[ind - 1]->_nthOutputTaken = _base;
    // std::cout << "_base " << _base << std::endl;
    _structure[ind]->_subBuckets.push_back(_structure[ind - 1]);
  }

  if (_setting->verbosity > 3)
    std::cout << "All Buckets in the last bucket!" << std::endl;
}

void Cascade::CreateTotalizerEncodeTree() {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  TimeMeasurement timeEncodeTree(&_dgpw->_timeVariables->createTree, true);

  _structure.back()->_isLastBucket = true;
  _structure.back()->CreateTotalizerEncodeTree(true);
  
  uint32_t lowTare = UINT32_MAX;
  if (_structure.size() > 1)
    lowTare = _structure[0]->_tares[0];
  
  _structure.back()->_sorter->_outputTree->AddBookkeepingForPL(true, lowTare); //TODO-Tobias: Note that this makes that if we have only one bucket, the last bucket will contain a vector of ones in the _leavesWeights. This is overhead. Should we change that?
  if (_setting->createGraphFile != "")
    _structure.back()->_sorter->_outputTree->DumpOutputTree(
        _setting->createGraphFile + std::to_string(_structure.back()->size()) +
            ".tgf",
        false);

  // test if the actualized bottom bucket values are correct by checking the root
  assert(_structure.size() == 1 || _structure.back()->_sorter->_outputTree->_leaves.size() == _structure.back()->_sorter->_outputTree->_leavesWeights.size());
  if (_setting->verbosity > 5) {
    std::cout << UINT64_MAX << std::endl;
    std::cout << "c root leaf values / soft clause values: " << std::endl;
    std::vector<std::pair<uint32_t, uint64_t>> sortedSCs;
    std::transform(_dgpw->_softClauses.begin(), _dgpw->_softClauses.end(), std::back_inserter(sortedSCs), [](const auto& sc) {
      return std::make_pair(sc->relaxationLit ^ 1, sc->weight);
    });
      std::sort(sortedSCs.begin(), sortedSCs.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
          return a.first < b.first;
      });

    for (size_t i = 0; i < _structure.back()->_sorter->_outputTree->_leaves.size(); i++) {

      if (_structure.size() != 1)
        std::cout << _structure.back()->_sorter->_outputTree->_leaves[i] << "(" << _structure.back()->_sorter->_outputTree->_leavesWeights[i] << "), ";
      else
        std::cout << _structure.back()->_sorter->_outputTree->_leaves[i] << ", ";
    }
    std::cout << std::endl;

    for (auto sc : sortedSCs) {
      std::cout << sc.first << "(" << sc.second << "), ";
    }
    std::cout << std::endl;
    for (size_t i = 0; i < _structure.back()->_sorter->_outputTree->_leaves.size(); i++) {
      assert(_structure.back()->_sorter->_outputTree->_leaves[i] == sortedSCs[i].first);
      assert(_structure.size() == 1 || _structure.back()->_sorter->_outputTree->_leavesWeights[i] == sortedSCs[i].second);
    }
  }
  assert(_structure.back()->_sorter->_outputTree->_leaves.size() == _dgpw->_softClauses.size());



  if (_setting->verbosity > 0)
    std::cout << "c #max sorter depth......: "
              << _structure.back()->_sorter->_outputTree->_depth + 1
              << std::endl;
  //    std::cout << "SIZE OF TREE: " <<
  //    _structure.back()->_sorter->_outputTree->_size << std::endl;
  if (_setting->verbosity < 1)
    return;

  if (_setting->verbosity > 3)
    std::cout << std::endl << "Totalizer Tree encoded!" << std::endl;
}

void Cascade::DumpNumberOfBucketsAndClauses() {
  std::cout << std::endl;
  // Dump BucketID
  DumpNumberOfBucketEntriesOrClauses(false, false, false, false, false, false);
  std::cout << std::endl;

  // Dump TopBucketEntries
  DumpNumberOfBucketEntriesOrClauses(true, false, false, false, false, false);
  // Dump BottomBucketEntries
  DumpNumberOfBucketEntriesOrClauses(false, true, false, false, false, false);
  std::cout << std::endl;

  // Dump Estimated TopBinaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, true, false, true, false);
  // Dump Calculated TopBinaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, false, true, true, false);
  std::cout << std::endl;

  // Dump Estimated TopTernaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, true, false, false, true);
  // Dump Calculated TopTernaryClauses
  DumpNumberOfBucketEntriesOrClauses(true, false, false, true, false, true);
  std::cout << std::endl;

  // Dump Estimated Bottom Binary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, true, false, true, false);
  // Dump Calculated Bottom Binary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, false, true, true, false);
  std::cout << std::endl;

  // Dump Estimated Bottom Ternary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, true, false, false, true);
  // Dump Calculated Bottom Ternary Clauses
  DumpNumberOfBucketEntriesOrClauses(false, true, false, true, false, true);
  std::cout << std::endl;
}

void Cascade::DumpNumberOfBucketEntriesOrClauses(bool top, bool bottom,
                                                 bool estimated,
                                                 bool calculated, bool binary,
                                                 bool ternary) {
  if (!top && !bottom && !estimated && !calculated)
    std::cout << std::setw(35) << "Bucket ID: ( ";
  else if (top && !estimated && !calculated)
    std::cout << std::setw(35) << "TopBucket Entries: ( ";
  else if (bottom && !estimated && !calculated)
    std::cout << std::setw(35) << "BottomBucket Entries: ( ";
  else if (estimated && top && binary)
    std::cout << std::setw(35) << "Estimated BinaryTopClauses: ( ";
  else if (estimated && top && ternary)
    std::cout << std::setw(35) << "Estimated TernaryTopClauses: ( ";
  else if (calculated && top && binary)
    std::cout << std::setw(35) << "Calculated BinaryTopClauses: ( ";
  else if (calculated && top && ternary)
    std::cout << std::setw(35) << "Calculated TernaryTopClauses: ( ";
  else if (estimated && bottom && binary)
    std::cout << std::setw(35) << "Estimated BinaryBottomClauses: ( ";
  else if (estimated && bottom && ternary)
    std::cout << std::setw(35) << "Estimated TernaryBottomClauses: ( ";
  else if (calculated && bottom && binary)
    std::cout << std::setw(35) << "Calculated BinaryBottomClauses: ( ";
  else if (calculated && bottom && ternary)
    std::cout << std::setw(35) << "Calculated TernaryBottomClauses: ( ";

  uint32_t actualValue(0);
  uint32_t sum(0);

  // for (auto bucket : stuint16_t ind = 0; ind < _structure.size(); ind++)
  for (uint16_t ind = 0; ind < _structure.size(); ind++) {
    if (!top && !bottom && !estimated && !calculated)
      actualValue = ind;
    else if (top && !estimated && !calculated)
      actualValue = _structure[ind]->_topEntries;
    else if (bottom && !estimated && !calculated)
      actualValue = static_cast<uint32_t>(_structure[ind]->_outputs->size());
    else if (estimated && top && binary)
      actualValue = _structure[ind]->_binaryTopClEstimated;
    else if (estimated && top && ternary)
      actualValue = _structure[ind]->_ternaryTopClEstimated;
    else if (calculated && top && binary)
      actualValue = _structure[ind]->_binaryTopCl;
    else if (calculated && top && ternary)
      actualValue = _structure[ind]->_ternaryTopCl;
    else if (estimated && bottom && binary)
      actualValue = _structure[ind]->_binaryBottomClEstimated;
    else if (estimated && bottom && ternary)
      actualValue = _structure[ind]->_ternaryBottomClEstimated;
    else if (calculated && bottom && binary)
      actualValue = _structure[ind]->_binaryBottomCl;
    else if (calculated && bottom && ternary)
      actualValue = _structure[ind]->_ternaryBottomCl;

    std::cout << std::setw(5) << actualValue;
    if (ind < _structure.size() - 1)
      std::cout << ", ";
    sum += actualValue;
  }
  if (!top && !bottom && !estimated && !calculated)
    std::cout << " )" << std::endl;
  else
    std::cout << " )   sum =" << std::setw(6) << sum << std::endl;
}

uint32_t Cascade::SolveAllTares() {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  TimeMeasurement timeSolvingTares(&_dgpw->_timeVariables->solvingTares, true);
  uint32_t currentresult(UNKNOW);
  std::vector<uint32_t> collectedAssumptions;

  if (_setting->verbosity > 0) {
    std::cout << std::endl
              << std::setw(90) << "Minimize all Tares from now on!"
              << std::endl;
    if (_setting->verbosity > 1) {
      std::cout << "----------------temporary solve "
                   "structure--------------------------------------------------"
                << std::endl;
    }
  }

  uint16_t sizeMinus = 2;
  if (_onlyByTares)
    sizeMinus = 1;

  // start with second last bucket!
  for (int32_t ind = static_cast<int32_t>(_structure.size() - sizeMinus);
       ind >= 0; ind--) {
    // assert(_estimatedWeightBoundaries[1] >= _satWeight);
    collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1);
    currentresult = _dgpw->Solve(collectedAssumptions);

    //        if (currentresult == ANTOM_SAT && _multipleCascade != nullptrptr)
    //            _multipleCascade->CalculateWeightBoundaries(_structure[ind]->_multiplicator);
    //        else if (currentresult == ANTOM_UNSAT && _multipleCascade !=
    //        nullptrptr)
    //            _multipleCascade->CalculateWeightBoundaries(-_structure[ind]->_multiplicator);

    if (currentresult == SAT) {
      if (_setting->verbosity > 3)
        std::cout << "SAT" << std::endl;
//            _dgpw->_lastModel = _dgpw->Model();
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit(collectedAssumptions.back());
      assert(rst);
#else
      _dgpw->AddUnit(collectedAssumptions.back());
#endif
      assert(collectedAssumptions.size() == 1);
      _dgpw->CalculateOverallOptimum(_satWeight, true);

    } else if (currentresult == UNSAT) {
      if (_setting->verbosity > 3) {
        std::cout << "UNSAT" << std::endl;
      }
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit(collectedAssumptions.back() ^ 1);
      assert(rst);
#else
      _dgpw->AddUnit(collectedAssumptions.back() ^ 1);
#endif
    }
    if (_setting->verbosity > 2) {
      std::cout << std::setw(50) << "Current Bucket Multiplicator: "
                << _structure[ind]->_multiplicator << std::endl;
      if (currentresult == SAT) {
        std::cout << std::setw(34) << "All Tares of Bucket " << ind
                  << " could be set! " << std::setw(40) << "ANTOM_SAT"
                  << std::endl;
        DumpModelOfTares(5);
      } else if (currentresult == UNSAT)
        std::cout << std::setw(31) << "At least one Tare of Bucket " << ind
                  << " couldn't be set! " << std::setw(40) << "ANTOM_UNSAT"
                  << std::endl;
      else
        std::cout << std::setw(88) << "TIMEOUT: " << currentresult << std::endl;
      currentresult = SAT;
    }

    if (_setting->verbosity > 3)
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;
    collectedAssumptions.pop_back();
  }

  return currentresult;
}

uint32_t Cascade::SolveTares(bool onlyWithAssumptions,
                             bool solveTareInLastBucketToo) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  //    std::cout << std::setw(50) << "Weight boundaries: " <<
  //    _estimatedWeightBoundaries[0] << " / " << _estimatedWeightBoundaries[1]
  //    << std::endl;
  if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] == 0 ||
      _structure.size() == 1) {
    return SAT;
  }

  TimeMeasurement timeSolvingTares(&_dgpw->_timeVariables->solvingTares, true);
  uint32_t currentresult(UNKNOW);

  // EXACT BOUND ENCODING -- STANDARD!!
  if (_setting->weightPlusOne || onlyWithAssumptions)
    return SolveTareWeightPlusOne(onlyWithAssumptions);

  if (onlyWithAssumptions) {
    std::cout
        << "c onlyWithAssumptions not yet implemented, to use that use the "
           "solveTareWeightPlusOne option!!!";
  }
  assert(!onlyWithAssumptions);

  if (_setting->verbosity > 0) {
    if (_setting->verbosity > 4)
      std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << std::endl
              << std::setw(90) << "Minimize Tares from now on!" << std::endl;
    if (_setting->verbosity > 1) {
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;
    }
  }

  uint16_t sizeMinus = 2;
  if (solveTareInLastBucketToo || _onlyByTares)
    sizeMinus = 1;

  // start with second last bucket!
  for (int16_t ind = static_cast<int16_t>(_structure.size() - sizeMinus);
       ind >= 0; ind--) {
    uint64_t diffEstimatedCurWeight;

    diffEstimatedCurWeight = _estimatedWeightBoundaries[1] - _dgpw->_satWeight;

    currentresult = _structure[ind]->SolveTares(diffEstimatedCurWeight);

    if (_setting->verbosity > 1) {
      std::cout << std::setw(50)
                << "Weight boundaries: " << _estimatedWeightBoundaries[0]
                << " / " << _estimatedWeightBoundaries[1] << std::endl;
      if (_setting->verbosity > 2) {
        std::cout << std::setw(50)
                  << "Current SAT Weight dgpw: " << _dgpw->_satWeight
                  << std::endl;
        std::cout << std::setw(50) << "Current Bucket Multiplicator: "
                  << _structure[ind]->_multiplicator << std::endl;
        if (currentresult == SAT) {
          std::cout << std::setw(34) << "All Tares of Bucket " << ind
                    << " could be set! " << std::setw(40) << "ANTOM_SAT"
                    << std::endl;
          DumpModelOfTares(5);
        } else if (currentresult == UNSAT)
          std::cout << std::setw(31) << "At least one Tare of Bucket " << ind
                    << " couldn't be set! " << std::setw(40) << "ANTOM_UNSAT"
                    << std::endl;
        else
          std::cout << std::setw(88) << "TIMEOUT: " << currentresult
                    << std::endl;
      }

      if (ind != 0)
        std::cout << "---------------------------------------------------------"
                     "------------"
                     "----------------------"
                  << std::endl;
      else
        std::cout << "========================================================="
                     "============"
                     "======================"
                  << std::endl;
    }

    // CASE TIMEOUT
    if (currentresult == UNKNOW /*|| _control->ReachedLimits()*/) {
      return UNKNOW;
    }
    //        std::cout << _estimatedWeightBoundaries[0] << std::endl;
    //        std::cout << _dgpw->_satWeight << std::endl;
    assert(_estimatedWeightBoundaries[0] <=
           (_dgpw->_satWeight));
    // assert(static_cast<int64_t>(_dgpw->_satWeight) <=
    // _estimatedWeightBoundaries[1]);
  }
  if (_setting->verbosity < 3)
    return SAT;
  //    DumpModelOfTares(3);
  if (_setting->verbosity > 0)
    std::cout << "All Tares are solved!" << std::endl << std::endl;
  return SAT;
}

void Cascade::CalculateBucketEntries() {
  uint32_t sumOfSizes(0);
  uint32_t maxBucketEntries(0);

  for (auto bucket : _structure) {
    sumOfSizes += bucket->size();
    if (bucket->size() > maxBucketEntries)
      maxBucketEntries = bucket->size();
  }
  if (_setting->verbosity < 1)
    return;
  std::cout << "c #buckets...............: " << _structure.size() << std::endl;
  std::cout << "c max Bucket entries.....: " << GetMaxBucketSize() << std::endl;
  std::cout << "c average Bucket entries.: " << sumOfSizes / _structure.size()
            << std::endl;
}

uint64_t
Cascade::CountSatisfiedSoftClauses(std::vector<SoftClause *> softclauses,
                                   const std::vector<uint32_t> &model,
                                   bool addWeight) {
  uint64_t result(0);
  // Proceed all soft clauses
  for (uint32_t i = 0; i != softclauses.size(); ++i) {
    uint32_t relaxlit = softclauses[i]->relaxationLit;

    // Just proceed satisfied triggers
    if (_dgpw->Model(relaxlit >> 1) == relaxlit) {
      std::vector<uint32_t> clause(softclauses[i]->clause);
      uint32_t pos = 0;
      for (; pos != clause.size(); ++pos) {
        // clause satisfied without trigger?
        if (_dgpw->Model(clause[pos] >> 1) == clause[pos]) {
          if (addWeight)
            result += softclauses[i]->weight;
          else
            result += 1;
          break;
        }
      }
    } else if (_dgpw->Model(relaxlit >> 1) != 0) {
      assert(_dgpw->Model(relaxlit >> 1) == (relaxlit ^ 1));
      if (addWeight)
        result += softclauses[i]->weight;
      else
        result += 1;
    }
  }
  return result;
}

uint64_t
Cascade::CountSatisfiedSoftClauses(Bucket *bucket,
                                   const std::vector<uint32_t> &model) {
  bool addWeight;
  std::vector<SoftClause *> softclauses;

  if (bucket == nullptr || bucket->_isLastBucket) {
    addWeight = true;
    softclauses = _softClauses;
  } else {
    addWeight = false;
    softclauses = *bucket->_softClauses;
    //        std::cout << "How many Softclauses: " << softclauses.size() <<
    //        std::endl;
  }

  uint64_t result = CountSatisfiedSoftClauses(softclauses, model, addWeight);

  //    std::cout << "Fulfilled Softclauses: " << result << std::endl;
  if (addWeight)
    _satWeight = result > _satWeight ? result : _satWeight;

  return result;
}

void Cascade::CalculateTotalBucketEntries(
    std::vector<SoftClauseNodes *> *SCTree, bool add) {
  _totalBucketEntriesperWeight = 0;
  _totalBucketOccurrences = 0;
  if (!add) {
    for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
      if (i >= _totalBucketEntries.size())
        _totalBucketEntries.push_back(0);
      else
        _totalBucketEntries[i] = 0;
    }
  }
  for (uint32_t i = 0; i < SCTree->size(); ++i) {
    _totalBucketEntriesperWeight += (*SCTree)[i]->inHowManyBuckets;
    _totalBucketOccurrences +=
        (*SCTree)[i]->GetOccurrences() * (*SCTree)[i]->size();
    for (uint32_t j = 0; j <= (*SCTree)[i]->highestBucket; ++j) {
      if ((*SCTree)[i]->hasSubBuckets &&
          (*SCTree)[i]->occursHowOftenInBucket[j])
        _totalBucketEntries[j] +=
            (*SCTree)[i]->occursHowOftenInBucket[j] * (*SCTree)[i]->size();
      else
        _totalBucketEntries[j] += (*SCTree)[i]->occursHowOftenInBucket[j];

      //            if ((*SCTree)[i]->occursHowOftenInBucket[j] >= 1)
      //                _totalBucketEntries[j]++;
      //            std::cout << "_totalBucketEntries[" << j << "]: " <<
      //            _totalBucketEntries[j] << std::endl;
    }
  }
  //    std::cout << "_totalBucketEntriesperWeight: " <<
  //    _totalBucketEntriesperWeight << std::endl; std::cout <<
  //    "_totalBucketOccurrences: " << _totalBucketOccurrences << std::endl;
}

bool Cascade::AtLeastTwoBucketsInCommon(SoftClauseNodes *SCN1,
                                        SoftClauseNodes *SCN2) {
  uint32_t minSize = (SCN1->occursHowOftenInBucket.size() <
                      SCN2->occursHowOftenInBucket.size())
                         ? SCN1->occursHowOftenInBucket.size()
                         : SCN2->occursHowOftenInBucket.size();
  uint32_t overlappings(0);
  for (uint32_t ind = 0; ind < minSize; ind++) {
    if (SCN1->occursHowOftenInBucket[ind] > 0 &&
        SCN2->occursHowOftenInBucket[ind] > 0) {
      overlappings++;
      if (overlappings > 1)
        return true;
    }
  }
  return false;
}

void Cascade::DumpSCNodeStructure(std::vector<SoftClauseNodes *> *dumpingSCTree,
                                  uint16_t verbosity) {
  if (_setting->verbosity < verbosity || dumpingSCTree->size() > 1000)
    return;
  std::cout << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "#" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "O" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "c" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "" << std::setw(8) << "c" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "#" << std::setw(8) << "u" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "B" << std::setw(8) << "r" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "u" << std::setw(8) << "r" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "c" << std::setw(8) << "e" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "k" << std::setw(8) << "n" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "e" << std::setw(8) << "c" << std::setw(15) << ""
            << std::setw(3) << "|"
            << "   " << _base << " up to the power of:" << std::endl;
  std::cout << std::setw(4) << "t" << std::setw(8) << "e" << std::setw(15) << ""
            << std::setw(3) << "|" << std::endl;
  std::cout << std::setw(4) << "s" << std::setw(8) << "s" << std::setw(15)
            << "Weight" << std::setw(3) << "|";
  for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
    std::cout << std::setw(4) << i;
  }
  std::cout << std::endl;
  std::cout << "--#occurrences---------------|";

  for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
    std::cout << "----";
  }
  std::cout << "----" << std::endl;

  for (uint32_t i = 0; i != (*dumpingSCTree).size(); ++i) {
    (*dumpingSCTree)[i]->dumpStructure(true, i);
  }
  std::cout << "-----------------------------|";
  for (uint32_t i = 0; i <= _numberOfBuckets; ++i) {
    std::cout << "----";
  }
  std::cout << "----" << std::endl;

  std::cout << std::setw(4) << _totalBucketEntriesperWeight << std::setw(8)
            << _totalBucketOccurrences << std::setw(15) << _sumOfSoftWeights
            << std::setw(3) << "|";
  for (uint32_t i = 0; i < _totalBucketEntries.size(); ++i) {
    std::cout << std::setw(4) << _totalBucketEntries[i];
  }
  std::cout << std::endl;
}

void Cascade::DumpBucketStructure(bool onlyLastBucket, uint16_t verbosity) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity < verbosity)
    return;
  uint16_t depth(0);
  uint16_t maxDepth(0);
  std::cout << "(SoftClauseEntries, SorterEntries, Multiplicator, Tares.empty(), nthOutputTaken, [outputs])" << std::endl;
  if (onlyLastBucket) {
    std::cout << "Last Bucket: " << std::endl;
    depth = _structure.back()->DumpAndGetMaxDepth(0);
    std::cout << "LastBucket has a depth of: " << depth << std::endl;
    std::cout << std::endl;
    return;
  }
  for (uint32_t i = 0; i < _structure.size(); i++) {
    std::cout << "Bucket " << i << ":" << std::endl;
    depth = _structure[i]->DumpAndGetMaxDepth(0);
    if (depth > maxDepth)
      maxDepth = depth;
    std::cout << "Bucket " << i << " has a depth of: " << depth << std::endl;
    std::cout << std::endl;
  }
  std::cout << "Max depth of all buckets is: " << maxDepth << std::endl
            << std::endl;
}

uint64_t Cascade::GetHighestMultiplicator() {
  //  assert(_highestBucketMultiplicator != 0);
  return _highestBucketMultiplicator;
}

uint32_t Cascade::GetMaxBucketSize() {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  uint32_t maxSize = 0;
  for (auto bucket : _structure)
    maxSize = bucket->size(true) * bucket->_nthOutputTaken > maxSize
                  ? bucket->size(true) * bucket->_nthOutputTaken
                  : maxSize;
  //    for (uint32_t ind=0; ind < _structure.size(); ind++)
  //    {
  //        std::cout << "Bucket " << ind << " has a size of " <<
  //        _structure[ind]->size(true) * _structure[ind]->_nthOutputTaken << "
  //        and every " << _structure[ind]->_nthOutputTaken << " output is
  //        taken." << std::endl; maxSize = ((_structure[ind]->size(true) *
  //        _structure[ind]->_nthOutputTaken) > maxSize) ?
  //        (_structure[ind]->size(true) * _structure[ind]->_nthOutputTaken) :
  //        maxSize;
  //    }
  return maxSize;
}

bool Cascade::AddNewBucketsTillSizeBoundary(uint32_t maxSize,
                                            bool onlySolveWithTares,
                                            bool addTareToLastBucket) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity > 3)
    std::cout << "Add another bucket because of size boundary: " << maxSize
              << std::endl;

  while (_structure.back()->size(true) > maxSize) {
    if (!AddNewBucketsTillMultiplicatorMatches(_highestBucketMultiplicator * 2,
                                               onlySolveWithTares,
                                               addTareToLastBucket))
      return false;
  }
  return true;
}

bool Cascade::AddNewBucketsTillMultiplicatorMatches(uint64_t maxMultiplicator,
                                                    bool onlySolveWithTares,
                                                    bool addTareToLastBucket) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity > 3)
    std::cout << "_highestBucketMultiplicator: " << _highestBucketMultiplicator
              << std::endl;
  assert(_structure.back()->_multiplicator <= maxMultiplicator);

  while (_structure.back()->_multiplicator < maxMultiplicator) {
    if (_setting->verbosity > 3) {
      std::cout << "_structure.back()->size(true): "
                << _structure.back()->size(true) << std::endl;
      std::cout << "onlySolveWithTares: " << onlySolveWithTares << std::endl;
    }

    // This bucket is too small to connect with another bucket.
    // It contains only tare and one result from the bucket before!
    // Calculate the results only by the tares!
    // TOBI: VERIFY!!! - Should be possible to calculate all results only by
    // tares
    if (_structure.back()->size(true) == 1 && onlySolveWithTares) {
      //            if (_setting->verbosity > 3)
      //                std::cout << "c MORE THAN TWO ENTRIES..: TRUE" <<
      //                std::endl;
      if (_setting->verbosity > 3)
        std::cout << "Add one tare to last bucket with size one!" << std::endl;
      AddTare(_structure.size() - 1);

      // here is the last position the bucket can be dumped
      // otherwise the subbuckets are encoded and not dumpable anymore :-)
      DumpBucketStructure(true, 4);

      CreateTotalizerEncodeTree();
      CalculateBucketEntries();

      // at least the tare should be 0, then by solving max tares to 1 is asked
      // for!
      uint32_t unitClauseVar =
          (_structure.back()->_sorter->GetOrEncodeOutput(0) << 1) ^ 1;

      if (_setting->verbosity > 3)
        std::cout << "Unit Clause for first entry in last Bucket Added: "
                  << unitClauseVar << std::endl;
      _dgpw->AddUnit(unitClauseVar);
      return false;
    } else {
      // if (_structure.back()->size(true) >= 1)
      AddAdditionalBucket();
    }

    // at least one element + one tare!!
    assert(_structure.back()->size(true) > 0);
  }

  if (addTareToLastBucket && _structure.back()->_tares.empty()) {
    if (_setting->verbosity > 3)
      std::cout << std::endl
                << "There is no tare at the last bucket! Add one!" << std::endl;
    AddTare(_structure.size() - 1);
  }

  return true;
}

void Cascade::AddTare(uint64_t position) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  uint32_t tare(_dgpw->NewVariable());
  _structure[position]->AddTare(tare);
  _tareWeight += _structure[position]->_multiplicator;

  vPL->write_comment("Tare added: " + vPL->var_name(tare));

  if (_setting->verbosity < 3)
    return;

  std::cout << "   Tare added: " << tare << std::endl;
}

void Cascade::AddAdditionalBucket() {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;

  if (_setting->verbosity > 3) {
    std::cout << "_highestBucketMultiplicator: " << _highestBucketMultiplicator
              << std::endl;
    std::cout << "Add another Bucket with multiplicator: "
              << _highestBucketMultiplicator;
  }
  _structure.back()->_isLastBucket = false;
  if (_structure.back()->_tares.empty()) {
    if (_setting->verbosity > 3)
      std::cout << std::endl
                << "There is no tare at the last bucket! Add one!" << std::endl;
    AddTare(_structure.size() - 1);
  }

  _structure.push_back(new Bucket(_dgpw, this));
  _structure.back()->_isLastBucket = true;
  _highestBucketMultiplicator =
      static_cast<uint64_t>(pow(_base, _structure.size() - 1));
  _structure.back()->_multiplicator = _highestBucketMultiplicator;
  _structure[_structure.size() - 2]->_nthOutputTaken = _base;
  _structure.back()->_subBuckets.push_back(_structure[_structure.size() - 2]);
  if (_structure.back()->_subBuckets.back()->_localMaxPos !=
      static_cast<uint32_t>(-1))
    _structure.back()->_localMaxPos = static_cast<uint32_t>(
        floor(_structure.back()->_subBuckets.back()->_localMaxPos / 2));
}

std::vector<uint32_t> Cascade::CalculateAssumptionsFor(int64_t weight,
                                                       int32_t startingPos, uint64_t &assumedTareWeights) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  
  assumedTareWeights = 0;

  if (startingPos == -1) {
    return {};
  }
  // std::cout << "weight: " << weight << std::endl;
  // std::cout << "_estimatedWeightBoundaries[0]: " <<
  // _estimatedWeightBoundaries[0] << std::endl; std::cout <<
  // "_estimatedWeightBoundaries[1]: " << _estimatedWeightBoundaries[1] <<
  // std::endl;

  assert(weight <= _estimatedWeightBoundaries[1]);
  assert(weight >= _estimatedWeightBoundaries[0]);
  assert(weight >= (_dgpw->_satWeight));

  //    std::cout << "c Calculate assumptions for weight: " << weight <<
  //    std::endl;

  std::vector<uint32_t> collectedAssumptions;

  int64_t upperDiff = _estimatedWeightBoundaries[1] - weight;
  int64_t lowerDiff = weight - _estimatedWeightBoundaries[0];

  // start with second last bucket!
  for (int32_t ind = startingPos; ind >= 0; --ind) {
    int64_t actualMult = static_cast<int64_t>(_structure[ind]->_multiplicator);

    if (_setting->verbosity > 5) {
      std::cout << std::setw(50) << "c Actual Diffs: " << lowerDiff << " / "
                << upperDiff << std::endl;
      std::cout << std::setw(50) << "c actualMult: " << actualMult << std::endl;
    }

    assert(upperDiff + lowerDiff == 2 * actualMult - 1);

    // Exact Bound Encoding
    if (lowerDiff >= actualMult) {
      assert(upperDiff < actualMult);
      collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1);
      assumedTareWeights += actualMult;
      lowerDiff -= actualMult;
    } else if (upperDiff >= actualMult) {
      assert(lowerDiff < actualMult);
      collectedAssumptions.push_back(_structure[ind]->_tares[0] << 1 ^ 1);
      upperDiff -= actualMult;
    } else {
      assert(false);
    }
  }

  for (auto assumption : _fixedTareAssumption) {
    collectedAssumptions.push_back(assumption);
  }

  if (_setting->verbosity < 3)
    return collectedAssumptions;

  std::cout << std::setw(51) << "Assumptions: (";
  for (auto assumption : collectedAssumptions) {
    std::cout << assumption << ", ";
  }
  std::cout << ")" << std::endl;
  return collectedAssumptions;
}

// Function to set the unit clauses in fine convergence. 
// This function adds the unit clauses of dominant tare variables that will not be changed anymore + at the end of the fine convergence.
// It is called before every solver call. 
int32_t Cascade::SetUnitClauses(int32_t startingPos, uint64_t &fixedTareValues) {
  if (_setting->verbosity > 6) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
  }

  
  TotalizerEncodeTree* tree = _structure.back()->_sorter->_outputTree  ;
  witnessT = vPL->get_new_substitution(); 
  
  uint64_t p = _structure.size() - 1;
  uint64_t s = _dgpw->_satWeight - (_structure.back()->kopt - 1) * (1 << p);
  vPL->write_comment("p = " + std::to_string(p) + " kopt = " + std::to_string(_structure.back()->kopt) + " _dgpw->_satWeight = " + std::to_string(_dgpw->_satWeight));

  // start with second last bucket!
  for (int32_t ind = startingPos; ind >= 0; ind--) {
    int64_t actualMult = static_cast<int64_t>(_structure[ind]->_multiplicator);

    if (_setting->verbosity > 5) {
      std::cout << std::setw(50)
                << "Weight boundaries: " << _estimatedWeightBoundaries[0]
                << " / " << _estimatedWeightBoundaries[1] << std::endl;
      std::cout << std::setw(50) << "actualMult: " << actualMult << std::endl;
    }

    assert(_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] ==
           2 * actualMult - 1);
    // The actual result implies that the currentTare and maybe more tares can
    // be set directly to TRUE
    if (_estimatedWeightBoundaries[1] -
            static_cast<int64_t>(_dgpw->_satWeight) <
        actualMult) {
      fixedTareValues += actualMult;
      vPL->write_comment("_estimatedWeightBoundaries[1] - static_cast<int64_t>(_dgpw->_satWeight) = " + std::to_string(_estimatedWeightBoundaries[1] -static_cast<int64_t>(_dgpw->_satWeight))+ " actualMult = " + std::to_string(actualMult));
//            std::cout << _structure[ind]->_tares[0]<< std::endl;
      if(lbTderivedfor < s-1){
        CreateShadowCircuitPL(s-1, witnessT, false);
        vPL->write_comment("Derive T >= s-1 for setting unit clauses in fine convergence");
        derivelbT(s-1, tree, witnessT);
        lbTderivedfor = s-1;
      }
      // Add unit clause that fixes the most dominant tare value that can be set.
      vPL->write_comment("Add unit clause that fixes the most dominant tare value that can be set.");
      vPL->rup_unit_clause(_structure[ind]->_tares[0] << 1);
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit(_structure[ind]->_tares[0] << 1);
      assert(rst);
#else
      _dgpw->AddUnit(_structure[ind]->_tares[0] << 1);
#endif
      _estimatedWeightBoundaries[0] += actualMult;
      continue;
    }
    // Corner Case!
    // If new result is already larger as maximal possible weight, we do not
    // need to try,
    // -> the corresponding tare can be directly set to FALSE
    else if (_estimatedWeightBoundaries[1] - actualMult >=
             static_cast<int64_t>(_dgpw->_sumOfSoftWeights)) {
      // Set values for T if  UB - actual tare value greater than the actual value, we can set an upper bound on T.
      vPL->write_comment("Set UB on T if value of objective (max) is close to upper bound.");
      vPL->unchecked_assumption_unit_clause((_structure[ind]->_tares[0] << 1) ^ 1);
#ifndef NDEBUG
      bool rst = _dgpw->AddUnit((_structure[ind]->_tares[0] << 1) ^ 1);
      assert(rst);
#else
      _dgpw->AddUnit((_structure[ind]->_tares[0] << 1) ^ 1);
#endif

      _estimatedWeightBoundaries[1] -= actualMult;
      continue;
    }
    // set assumption vector.
    else {
      if (_setting->verbosity > 2) {
        std::cout << std::setw(50)
                  << "Weight boundaries: " << _estimatedWeightBoundaries[0]
                  << " / " << _estimatedWeightBoundaries[1] << std::endl;
        std::cout << std::setw(50) << "actualMult: " << actualMult << std::endl;
      }
      return ind;
    }
  }

  if (_estimatedWeightBoundaries[1] - _dgpw->_satWeight > 0)
    return 0;
  else
    // result is found.
    return -1;
}

uint32_t Cascade::SolveTareWeightPlusOne(bool onlyWithAssumptions) {
  TimeMeasurement timeSolvingTares(&_dgpw->_timeVariables->solvingTares, true);
  uint32_t currentresult(SAT);
  std::vector<uint32_t> collectedAssumptions;
  std::vector<uint32_t> collectedAssumptionsMinusOne;
  std::vector<uint32_t> lastCollectedAssumptions;

  if (_setting->verbosity > 0) {
    if (_setting->verbosity > 6)
      std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << std::endl
              << std::setw(90) << "Try to set tares for actual weight + 1!"
              << std::endl;
    if (_setting->verbosity > 1)
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;
  }

  int32_t startingPos = _structure.size() - 2;

  // Probably solving ONLY WITH TARES -- special case!!! - not relevant for
  // standard solving procedure
  if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] >
      static_cast<int64_t>(_highestBucketMultiplicator))
    startingPos++;
  //    else if (_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0]
  //    == 1 && _highestBucketMultiplicator <= 2)
  //        //border case, only 2 buckets... est weight diffs = 1
  //        startingPos++;
  else
    assert(_estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] >=
               static_cast<int64_t>(_highestBucketMultiplicator / 2) ||
           _estimatedWeightBoundaries[1] - _estimatedWeightBoundaries[0] == 1);

  //    std::cout << "startingPos: " << startingPos << std::endl;
  uint64_t fixedTareWeights = 0;
  uint64_t assumedTareWeights = 0;
  
  
  if(currentresult == SAT){
    vPL->write_comment("test currentresult = SAT");
  }
  else{
    vPL->write_comment("test currentresult = UNSAT");
  }

  while (currentresult == SAT) {
    if (!onlyWithAssumptions) {
      startingPos = SetUnitClauses(startingPos, fixedTareWeights);
      //            std::cout << "SP: " << startingPos << std::endl;
    }

    // some corner case
    if (startingPos == -1 || static_cast<int64_t>(_dgpw->_satWeight) ==
                                 _estimatedWeightBoundaries[1]) {
      // Setting all tare variables for current satisfied weight
      // In this case, all tares are satisfiable, hence objective is already optimal.
      // Or, might also be that all tares have been set by set unit clauses (in the case we are in the last position of the coarse convergence).
      collectedAssumptions = CalculateAssumptionsFor(
          static_cast<int64_t>(_dgpw->_satWeight), startingPos, assumedTareWeights);
      break;
    }

    //        collectedAssumptionsMinusOne =
    //        CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight),
    //        startingPos); lastCollectedAssumptions = collectedAssumptions;
    collectedAssumptions = CalculateAssumptionsFor(
        static_cast<int64_t>(_dgpw->_satWeight) + 1, startingPos, assumedTareWeights);

    // PROOF: The proof for this SAT solver call is required. Should be handled
    // directly by the SAT solver.
    
    std::cout << "c assumedTareWeights: " << assumedTareWeights << std::endl;
    std::cout << "c fixedTareWeights: " << fixedTareWeights << std::endl;
    if (assumedTareWeights > 0){ 
      std::cout << "c T = " << assumedTareWeights + fixedTareWeights - 1 << std::endl;
    }  
    else{
      std::cout << "c T = " << fixedTareWeights << std::endl;
    }

    currentresult = _dgpw->Solve(collectedAssumptions);
    
    if(currentresult == SAT){
      _dgpw->_pacose->SendVPBModel(_structure[_structure.size()-1]->_sorter->_outputTree->_tares);
    }

    //        std::cout << "tried SATWeight: " << _dgpw->_satWeight + 1 <<
    //        std::endl;
    if (_setting->verbosity > 1)
      std::cout << "-----------------------------------------------------------"
                   "--------------"
                   "------------------"
                << std::endl;

    //        std::cout << "CURRENT RESULT: " << currentresult << std::endl;
    if (currentresult == SAT) {
      _dgpw->CalculateOverallOptimum(0, true);
      //            std::cout << "real SATWeight:  " << _dgpw->_satWeight <<
      //            std::endl;
    }
  }
  //    collectedAssumptions =
  //    CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight) - 1,
  //    startingPos);

  if (currentresult == SAT) {
    if (_setting->verbosity > 0)
      std::cout << "c SAT AFTER SOLVING TARES!" << std::endl;
    //        for (auto unitClause : collectedAssumptions) {
    //            _dgpw->AddUnit(unitClause);
    //        }
  } else {
    //        std::cout << "collect assumptions for: " << _dgpw->_satWeight <<
    //        std::endl;
    if (_setting->verbosity > 0)
      std::cout
          << "c UNSAT AFTER SOLVING TARES!, solve again with right assumptions!"
          << _dgpw->_satWeight << std::endl;
    collectedAssumptions = CalculateAssumptionsFor(
        static_cast<int64_t>(_dgpw->_satWeight), startingPos, assumedTareWeights);

    //        assert(_dgpw->Solve(collectedAssumptions)==SAT);
    //        if (onlyWithAssumptions && _dgpw->Solve(collectedAssumptions) !=
    //        10)
    //        {
    //            assert(false);
    //        }
    //        std::cout << "calc!1" << std::endl;
    //        _dgpw->CalculateOverallOptimum(0,true);
    //        collectedAssumptions =
    //        CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight) +
    //        1, startingPos);
    //        assert(_dgpw->Solve(collectedAssumptions)!=SAT);
  }
  vPL->write_comment("Fine convergence has finished. We now set the tare variables as they are for the optimal solution.");

  

  std::string cmnt = "Collected assumptions = "; for(auto a : collectedAssumptions){cmnt += vPL->to_string(a) + " ";}; vPL->write_comment(cmnt);
  for (auto unitClause : collectedAssumptions) {
    //        _fixedTareAssumption.clear();
    if (onlyWithAssumptions) {
      _fixedTareAssumption.push_back(unitClause);
      //            std::cout << "TareAssumptions: " << unitClause << std::endl;
    } else {
      // PROOF: Derive that the tare can be fixed to right values.
      // Fine convergence has finished. We will now set the tare variables as they are for the optimal solution as unit clauses.
      // vPL->unchecked_assumption_unit_clause(unitClause);
      if(!ubTderived){
        TotalizerEncodeTree* tree = _structure.back()->_sorter->_outputTree  ;
        uint64_t p = _structure.size() - 1;
        uint64_t s = _dgpw->_satWeight - (_structure.back()->kopt - 1) * (1 << p);
        if(lbTderivedfor < s-1){
          CreateShadowCircuitPL(s-1, witnessT, false);
          vPL->write_comment("Derive T >= s-1 for setting unit clauses in fine convergence");
          derivelbT(s-1, tree, witnessT);
          lbTderivedfor = s-1;
        }
        deriveubT(s-1, tree, witnessT);
        ubTderived = true;
      }

      vPL->rup_unit_clause(unitClause);
      _dgpw->AddUnit(unitClause);
      //            std::cout << "AddUnit: " << unitClause << std::endl;
    }
  }
  if (onlyWithAssumptions &&
      _dgpw->Solve(_fixedTareAssumption) != SAT) {
    std::cout << "c ERROR: Wrong result!" << std::endl;
    assert(false);
  } else if (!onlyWithAssumptions && _dgpw->Solve() != SAT) {
    std::cout << "c ERROR: Wrong Solve Result!" << std::endl;
    assert(false);
  }
  _dgpw->_pacose->SendVPBModel(_structure[_structure.size()-1]->_sorter->_outputTree->_tares);
  _dgpw->CalculateOverallOptimum(0, true);
  //    std::cout << "SolveValue: " << _dgpw->Solve() << std::endl;
  //    std::cout << "SATWeight: " << _dgpw->_satWeight << std::endl;
  //    _dgpw->CalculateOverallOptimum(0,true);
  //    std::cout << "SATWeight: " << _dgpw->_satWeight << std::endl;
  //    collectedAssumptions =
  //    CalculateAssumptionsFor(static_cast<int64_t>(_dgpw->_satWeight) + 1,
  //    startingPos); currentresult = _dgpw->Solve(collectedAssumptions);
  //     assert(_dgpw->Solve() == SAT);

  //    std::cout << "SOLVE TARE WEIGHT + 1 - currentresult: " << currentresult
  //    << std::endl; assert(currentresult == SAT);

  if (_setting->verbosity > 1)
    std::cout << "============================================================="
                 "================"
                 "=============="
              << std::endl;

  if (currentresult == UNKNOW /*|| _control->ReachedLimits()*/) {
    _dgpw->_resultUnknown = true;
    std::cout << std::setw(88) << "TIMEOUT: " << currentresult << std::endl;
    return UNKNOW;
  }

  if (_setting->verbosity > 1)
    std::cout << "All Tares are solved!" << std::endl << std::endl;
  return SAT;
}


// PROOF LOGGING: Creation of shadow circuit in the proof
void Cascade::CreateShadowCircuitPL(uint64_t s, substitution& w, bool check_for_already_shadowed_lits){
  vPL->write_comment("Creation of shadow circuit for T = " + std::to_string(s));
  std::unordered_map<uint32_t, uint64_t> valuesTareVariables; // tare[var] = 2^i if variable needs to be assigned 1 such that T = s or 0 otherwise.
  valuesTareVariables.reserve(_structure.size()-1); 

  // Adder caching: We should visit the same node only once. 
  // TODO-Dieter: Might be even better to implement it for variables or to take into account that only nodes that are used multiple times have to be added.
  std::unordered_set<uintptr_t> nodesAlreadyVisited; 

  size_t witnessSize = vPL->get_substitution_size(w);

  for(int i = _structure.size()-2; i>=0; i--){
    wght m = (1 << i);

    if(s >= m){
      valuesTareVariables[_structure[i]->_tares[0]] = m;
      s -= m;
      if(witnessSize == 0) 
        vPL->add_boolean_assignment(w, _structure[i]->_tares[0], 1);
    }
    else{
      valuesTareVariables[_structure[i]->_tares[0]]  = 0;
      if(witnessSize == 0) 
        vPL->add_boolean_assignment(w, _structure[i]->_tares[0], 0);
    }
  }
  

  // Print comment for debugging reasons:
  std::string contentofvaluesTareVariables = "Tares: "; 
  for(const auto& pair : valuesTareVariables){
    if(pair.second > 0)
      contentofvaluesTareVariables += vPL->var_name(pair.first) + " -> 1 ";
    else 
      contentofvaluesTareVariables += vPL->var_name(pair.first) + " -> 0 ";
  }
  vPL->write_comment(contentofvaluesTareVariables);
  contentofvaluesTareVariables = "Values for Tares: ";
  for(const auto& pair : valuesTareVariables){
    contentofvaluesTareVariables += vPL->var_name(pair.first) + " (" + std::to_string(pair.first) + ") -> " + std::to_string(pair.second) + " ";
  }
  vPL->write_comment(contentofvaluesTareVariables);

  vPL->write_comment("Shadow Circuit creation - Start traversing the tree");
  // Create the shadow circuit recursively by traversing the encoding tree
  CreateShadowCircuitPL_rec(w, _structure.back()->_sorter->_outputTree,valuesTareVariables, nodesAlreadyVisited, true, check_for_already_shadowed_lits); 
  vPL->write_comment("Shadow Circuit creation - End traversing the tree");
  // vPL->write_comment("Done creation of shadow circuit. Created witness: " + w);

  // TODO-Dieter: Probably not necessary to do it again for all proof goals, especially during coarse convergence. 
  vPL->write_comment("Derive proofgoals for satisfied output literals in Coarse convergence.");
  constraintid cxnLBcurrentGBMO = _dgpw->_pacose->derive_LBcxn_currentGBMO();
  std::vector<uint32_t>* outputs = &_structure.back()->_sorter->_outputTree->_encodedOutputs;
  // std::vector<constraintid>* cxnoutputs = &_structure.back()->cxn_sat_outputlit;


  for(uint32_t i = 0; i < outputs->size(); i++){
    if(outputs->at(i) == 0 || i == _structure.back()->kopt) continue;
    
    vPL->write_comment("Derive that shadow-variable (in shadow circuit) for variable assigned false in coarse convergence is false as well.");
    cuttingplanes_derivation cpder = vPL->CP_constraintid(vPL->getReifiedConstraintLeftImpl(variable(vPL->get_literal_assignment(w, toVeriPbVar(outputs->at(i))))));
    cpder = vPL->CP_multiplication(cpder, _dgpw->_pacose->_GCD);    
    cpder = vPL->CP_saturation( vPL->CP_division(vPL->CP_addition(cpder, vPL->CP_constraintid(cxnLBcurrentGBMO)), _dgpw->_pacose->_GCD) ); 
    vPL->write_CP_derivation(cpder);     
  }
  
}

//TODO-Dieter: Write function to print valuesTareVariables and call it multiple times to see why it's not working.

void Cascade::CreateShadowCircuitPL_rec(substitution& w, const TotalizerEncodeTree* tree, const std::unordered_map<uint32_t, uint64_t>& valuesTareVariables, std::unordered_set<uintptr_t>& nodesAlreadyVisited, bool is_root, bool check_for_already_shadowed_lits){
  uintptr_t node_id = reinterpret_cast<uintptr_t>(tree);
  if(nodesAlreadyVisited.find(node_id) != nodesAlreadyVisited.end()) return; // Adder Caching, only visit the same node once.
  nodesAlreadyVisited.emplace(node_id);

  std::string leavesstr = std::to_string(node_id) + " - Leaves: ";
  for(int i = 0; i < tree->_leaves.size(); i++){
    if(tree->_leavesWeights.size() > 0)
      leavesstr += "(" + vPL->to_string(tree->_leaves[i]) + "," + std::to_string(tree->_leavesWeights[i]) +  ") ";
    else
      leavesstr += vPL->to_string(tree->_leaves[i]) + " ";
    
  }
  std::string encodedoutputsstr = "EncodedOutputs: ";
  for(auto output : tree->_encodedOutputs){
    encodedoutputsstr += vPL->var_name(output) + " ";
  }
  std::string taresstr = "Tares: ";
  for(auto tare : tree->_tares){
    taresstr += vPL->var_name(tare) + " ";
  }
  vPL->write_comment(leavesstr + encodedoutputsstr + taresstr);

  for(uint32_t k = 0; k < tree->_encodedOutputs.size(); k++){
    if(tree->_encodedOutputs.size() == 1) {
      vPL->write_comment("Leaf " + (tree->_leaves.size() > 0 ? vPL->to_string(tree->_leaves[0]) : vPL->to_string(tree->_tares[0])) + " represented by " + (tree->_leaves.size() > 0 ? std::to_string(variable(tree->_leaves[0])) : std::to_string(variable(tree->_tares[0]))));
      continue; // We are in a leaf.
    }
    if(tree->_encodedOutputs.size() > 0 && tree->_encodedOutputs[k] == 0){ 
      vPL->write_comment("Encoded outputs == 0");
      continue; // Current variable was not derived yet
    }

    VeriPB::Var encodedVar = toVeriPbVar(tree->_encodedOutputs[k]);

    if(check_for_already_shadowed_lits && vPL->has_literal_assignment(w, encodedVar)){
      vPL->write_comment("Variable " + vPL->var_name(encodedVar) + " is already assigned");
      continue;
    }

    VeriPB::Var shadowvar = vPL->new_variable_only_in_proof();
    VeriPB::Lit shadowlit = create_literal(shadowvar, false);
    VeriPB::Lit shadowlitneg = create_literal(shadowvar, true);

    // std::string leavesstr = (tree->_leavesWeights.size() > 0) ? vPL->sequence_to_string(tree->_leaves, tree->_leavesWeights) : vPL->sequence_to_string(tree->_leaves, tree->_leavesWeights);
    // std::string leavesstr = vPL->sequence_to_string(tree->_leaves);
    
    vPL->write_comment("Encoding shadow literal  " + vPL->to_string(shadowlitneg) + " for " + vPL->to_string(create_literal(encodedVar, true)));

    if(tree->_isBottomBucket){
      vPL->write_comment("isBottomBucket");
      // ~z_k <-> O' + ~T >= (k+1) * 2^exp <-> O'  >= (k+1) * 2^exp - ~T
      wght T = 0;
      for(auto tare : tree->_tares){
        assert(valuesTareVariables.find(tree->_tares[0]) != valuesTareVariables.end());
        T += valuesTareVariables.at(tare);
      }

      wght negT = (1 << tree->_tares.size()) - 1 - T;

      wght rhs = (k+1)*(1 << tree->_exponent);
      // Note that the VeriPB proof checker always translates a constraint with negative left hand side to a constraint with rhs = 0;
      vPL->write_comment("negT = " + std::to_string(negT) + " and rhs = " + std::to_string(rhs));
      if(negT > rhs)
        rhs = 0;
      else
        rhs = rhs - negT;

      vPL->reificationLiteralLeftImpl(shadowlitneg, tree->_leaves, tree->_leavesWeights, rhs, is_root);
      vPL->reificationLiteralRightImpl(shadowlitneg, tree->_leaves, tree->_leavesWeights, rhs, is_root);
    }
    else{
      // ~y_k <-> O' + ~t >= k+1 <-> O' >= k + t.
      wght negt = 0;

      assert(tree->_tares.size() < 2);
      if(tree->_tares.size() == 1){
        assert(valuesTareVariables.find(tree->_tares[0]) != valuesTareVariables.end());
        if(valuesTareVariables.at(tree->_tares[0]) == 0)
          negt = 1;
      }

      vPL->reificationLiteralLeftImpl(shadowlitneg, tree->_leaves, k+1-negt, is_root);
      vPL->reificationLiteralRightImpl(shadowlitneg, tree->_leaves, k+1-negt, is_root);
    }
    
    // add to substitution
    vPL->add_literal_assignment(w, encodedVar, shadowlit);
  }

  vPL->write_comment("Done creating shadow circuit for node.");

  if(tree->_child1 != nullptr) CreateShadowCircuitPL_rec(w, tree->_child1, valuesTareVariables, nodesAlreadyVisited, false, check_for_already_shadowed_lits);
  if(tree->_child2 != nullptr) CreateShadowCircuitPL_rec(w, tree->_child2, valuesTareVariables, nodesAlreadyVisited, false, check_for_already_shadowed_lits);
}

// Functions for Proof Logging
constraintid Cascade::derivelbT(uint64_t lb, TotalizerEncodeTree* tree, substitution witnessT){
  std::vector<uint32_t> Clits; std::vector<uint64_t> Cwghts;
  for(int i = 0; i < tree->_tares.size(); i++){
      Clits.push_back(create_literal(tree->_tares[i], false));
      vPL->write_comment("Variable " + vPL->var_name(tree->_tares[i]) + " with weight " + std::to_string(1 << (tree->_tares.size() - 1 -  i )));
      Cwghts.push_back(1 << (tree->_tares.size() - 1 -  i ));
  }
  
  // constraintid VeriPbProofLogger::redundanceBasedStrengthening(&lits, &weights, const wght RHS, const substitution &witness)
  return vPL->redundanceBasedStrengthening(Clits, Cwghts, lb, witnessT);
}

constraintid Cascade::deriveubT(uint64_t ub, TotalizerEncodeTree* tree, substitution witnessT){
  std::vector<uint32_t> Clits; std::vector<uint64_t> Cwghts;
  for(int i = 0; i < tree->_tares.size(); i++){
      Clits.push_back(create_literal(tree->_tares[i], true));
      Cwghts.push_back(1 << (tree->_tares.size() - 1 -  i ));
  }
  vPL->write_comment("Derive T =< s-1 for setting unit clauses in fine convergence");
  uint64_t p = _structure.size() - 1;
  return vPL->redundanceBasedStrengthening(Clits, Cwghts, (1 << p) + 1 -  ub , witnessT);
}

} // namespace DGPW
} // Namespace Pacose