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

#include <cmath>
#include <iomanip>
#include <algorithm>

// Include dgpw related headers.
#include "bucket.h"
#include "cascade.h"
#include "dgpw.h"
#include "sorter.h"
#include "../Softclause.h"
#include "timemeasurement.h"
#include "timevariables.h"
#include "totalizerencodetree.h"

namespace Pacose {
namespace DGPW {
// Constructor
Sorter::Sorter(uint32_t size, DGPW *dgpw)
    : _dgpw(dgpw), _setting(dgpw->_dgpwSetting), _outputs(),
      _outputTree(nullptr), _preprocessingOutputs(), _sortedVectors(),
      _numberOfOutput(), _softClauses(), _tare(), _minContra(size), _depth(0),
      _minUnsatisfied(0), _minSatisfied(0), _fakeLiterals(0), _sumOfWeights(0),
      _proceeded(false), _proceedNext(false), _tarePosition(0),
      _sorterType(_setting->networkType) {
  assert(dgpw != NULL);
}

Sorter::~Sorter(void) { delete _outputTree; }

// PROOF: If new parts are encoded then we need a proof for that.
uint32_t Sorter::GetOrEncodeOutput(uint32_t position, bool encodeOnlyOnes) {
  if (_setting->verbosity > 6)
    std::cout << __PRETTY_FUNCTION__ << ", position, encodeOnlyOnes: " << position << ", " << encodeOnlyOnes << std::endl;

  if (_setting->encodeStrategy == ENCODEONLYIFNEEDED) {
    //        std::cout << "ENCODEONLYIFNEEDED" << std::endl;
    // Get Number of clauses
    uint32_t tmpCl = _dgpw->Clauses();
    //        uint32_t tmpBinCl = _dgpw->CurrentBinaryClauses();
    //        uint32_t tmpTerCl = _dgpw->CurrentTernaryClauses();
    TimeMeasurement timeEncoding(&_dgpw->_timeVariables->encoding, true);
    uint32_t var = _outputTree->ReturnOutputEncodeIfNecessary(position, this,
                                                              encodeOnlyOnes);
    //    _dgpw->_addedClauses += (_dgpw->Clauses() - tmpCl);
    //        _dgpw->_addedBinaryClauses += (_dgpw->CurrentBinaryClauses() -
    //        tmpBinCl); _dgpw->_addedTernaryClauses +=
    //        (_dgpw->CurrentTernaryClauses() - tmpTerCl);
    if (_setting->verbosity >
        4 /*&& (_dgpw->CurrentBinaryClauses() - tmpBinCl) > 0*/) {
      std::cout << std::setw(50)
                << "Added Clauses: " << _dgpw->StaticClauses() - tmpCl << " / "
                << _dgpw->StaticClauses() - _dgpw->_clausesBefore << std::endl;
      //            std::cout << std::setw(50) << "Added Binary Clauses: " <<
      //            _dgpw->CurrentBinaryClauses() - tmpBinCl << " / " <<
      //            _dgpw->_addedBinaryClauses  << std::endl; std::cout <<
      //            std::setw(50) << "Added Ternary Clauses: "<<
      //            _dgpw->CurrentTernaryClauses() - tmpTerCl << " / " <<
      //            _dgpw->_addedTernaryClauses << std::endl;
    }
    //        std::cout << "c OutputVariable: " << var << std::endl;
    return var;
  } else {
    //        std::cout << "ENCODEALWAYS!" << std::endl;
    return _outputs[position];
  }
}

uint32_t Sorter::GetCurrentTare(void) const { return _tare[_tarePosition]; }

void Sorter::SimpleMergeWithSorter(Sorter &sorter) {
  uint32_t thissortersize = size();

  // Merge current triggers
  _outputs.insert(_outputs.end(), sorter._outputs.begin(),
                  sorter._outputs.end());

  // Merge original trigger indices
  _softClauses.insert(_softClauses.end(), sorter._softClauses.begin(),
                      sorter._softClauses.end());

  if (_sorterType == TOTALIZER) {
    //        if (_dgpw->_featureTest)
    EncodingUVDiagonal(0, thissortersize, size());
    //        else
    //            EncodingUV(0, thissortersize, size());
  }
}

// Set clause bypasses in sorter, given by found trivial cases in
// "findSorterBypasses()" Assumes that the triggervariables in "m_inputsParts"
// are up to date
// void Sorter::SetSorterCSC(void)
//{
//    assert(_dgpw != NULL);

//    // Currently does not work with weighted maxsat
//    if( (_setting->setCSC == 0 ) || (_setting->application == WEIGHTEDMAXSAT )
//    ) { return; }

//    // No of soft clauses in this part
//    uint32_t noofsoftclauses(size());
//    uint32_t softclausessize(_softClauses.size());

//    uint32_t diff = noofsoftclauses - softclausessize;

//    if ( _minContra < noofsoftclauses && _minContra > 0 )
//    {
//        uint32_t var( _outputs[noofsoftclauses - _minContra - diff] );
//#ifndef NDEBUG
//        bool rst = _dgpw->AddUnit( var<<1 ); assert( rst );
//#else
//        _dgpw->AddUnit( var<<1 );
//#endif
//    }

//    for (uint32_t i = 0; i != softclausessize; ++i)
//    {
//        uint32_t softClauseContra = _softClauses[i]->contra;
//        if (softClauseContra > _minContra )
//        {
//            uint32_t var( _outputs[noofsoftclauses - softClauseContra - diff]
//            );
////            AddBypassClause(_softClauses[i]->relaxationLit, var<<1);
//        }
//    }
//}

// Pushes new part in front of "inputsParts"
void Sorter::EncodeSorter(uint64_t weight) {
  if (size() <= 1) {
    return;
  }

  // FindCSC();

  /*
  std::cout << "inputs: " << std::endl;
  for(uint32_t i = 0; i != size(); ++i )
    {
      std::cout << outputs[i] << " ";
    }
  std::cout << std::endl;
  */

  // clear model
  //_dgpw->TrivialAssignment();

  // Creates and adds sorter network clauses
  if (_sorterType == TOTALIZER) {
    TotalizerPhi(0, size());
  } else {
    // TODO
    assert(false);
  }

  /*
  std::cout << "outputs: " << std::endl;
  for(uint32_t i = 0; i != size(); ++i )
    {
      std::cout << outputs[i] << " ";
    }
  std::cout << std::endl;
  */

  //    SetVerticalBypasses();
  _depth = weight;
  //    SetSorterCSC();
}

void Sorter::MergeTotalizer(std::vector<uint32_t> &X) {
  /*
  std::cout << "inputs: " << std::endl;
  for(uint32_t i = 0; i != size(); ++i )
    {
      std::cout << outputs[i] << " ";
    }
  std::cout << std::endl << "and" << std::endl;
  for(uint32_t i = 0; i != X.size(); ++i )
    {
      std::cout << X[i] << " ";
    }
  std::cout << std::endl;
  */

  uint32_t outputsize = _outputs.size();

  _outputs.insert(_outputs.end(), X.begin(), X.end());
  uint32_t newoutputsize = _outputs.size();

  // Skip, if X is empty
  if (outputsize != newoutputsize) {
    EncodingUV(0, outputsize, newoutputsize);
    //        SetVerticalBypasses();
  }

  /*
  std::cout << "outputs: " << std::endl;
  for(uint32_t i = 0; i != size(); ++i )
    {
      std::cout << outputs[i] << " ";
    }
  std::cout << std::endl;
  */
}

void Sorter::TotalizerPhi(uint32_t lo, uint32_t hi) {
  if ((hi - lo) > 1) {
    assert(lo < hi);
    uint32_t m = ((hi - lo) >> 1);

    // TOBI::: here function to get to encode list for lo and list for high!

    // Calculate recursively the unary representation.
    TotalizerPhi(lo, lo + m);
    TotalizerPhi(lo + m, hi);
    //        if (_dgpw->_featureTest)
    EncodingUVDiagonal(lo, lo + m, hi);
    //        else
    //            EncodingUV(lo, lo+m, hi);
  }
}

// void Sorter::TotalizerPhiChosenOutputs(uint32_t lo, uint32_t hi,
// std::vector<uint32_t>* outputList)
void Sorter::TotalizerPhiChosenOutputs(uint32_t lo, uint32_t hi,
                                       std::vector<uint32_t> *outputList) {
  // TotalizerEncodedOutputsNode
  if ((hi - lo) > 1) {
    assert(lo < hi);
    uint32_t m = ((hi - lo) >> 1);

    std::vector<uint32_t> *outputListLo = nullptr;
    std::vector<uint32_t> *outputListHi = nullptr;

    // TOBI::: here function to get to encode list for lo and list for high!
    //        if (outputList != NULL)
    //            CalculateNextOutputLists(lo, lo+m, hi, outputListLo,
    //            outputListHi)

    // Calculate recursively the unary representation.
    TotalizerPhiChosenOutputs(lo, lo + m, outputListLo);
    TotalizerPhiChosenOutputs(lo + m, hi, outputListHi);
    EncodingUVDiagonal(lo, lo + m, hi, outputList);
  }
}

void Sorter::EncodingUV(uint32_t beginA, uint32_t endA, uint32_t endB) {
  assert(_dgpw != nullptr);
  // std::cout << __func__ << " " << beginA << " " << endA << " " << endB <<
  // std::endl; Variables----------------------------------------------------
  std::vector<uint32_t> clause;
  // if direction is true 0's at the beginning, 1's at the end
  // if direction is false the other way around.
  bool direction(true);
  // std::cout << "beginA: " << beginA << "   endA: " << endA << "   endB: " <<
  // endB << "   outputs.size(): " << _outputs.size() << std::endl;
  assert(beginA <= endA);
  assert(endA < endB);
  //    if( endA == endB ) {
  //        std::cout << "A==B" << std::endl;
  //        return;
  //    }
  assert(endB <= _outputs.size());
  // -------------------------------------------------------------

  std::vector<uint32_t> W;
  for (uint32_t a = beginA; a < endB; ++a) {
    W.push_back(_dgpw->NewVariable());
  }

  uint32_t beginB = endA;
  // std::cout << std::endl;
  for (uint32_t a = beginA; a <= endA; a++) {
    uint32_t an = a - beginA;
    for (uint32_t b = beginB; b <= endB; b++) {
      uint32_t bn = b - beginB;
      // To assure to get 1's at the end of A and B.
      /*
      if (a > beginA && b > beginB)
          std::cout << "( " << a - 1 << ", " << b - 1 << ", " << an + bn - 1 <<
      " ), "; else if (a == beginA && b > beginB) std::cout << "(  , " << b - 1
      << ", " << an + bn - 1 << " ), "; else if (b == beginB && a > beginA)
          std::cout << "( " << a - 1 << ",  , " << an + bn - 1 << " ), ";
      std::cout << std::endl;
      */
      if (_setting->encode01 && ((a > beginA) || (b > beginB))) {
        if (a > beginA && b > beginB) {
          clause.push_back((_outputs[a - 1] << 1) ^ !direction);
          clause.push_back((_outputs[b - 1] << 1) ^ !direction);
          clause.push_back((W[an + bn - 1] << 1) ^ direction);
        } else if (a == beginA && b > beginB) {
          clause.push_back((_outputs[b - 1] << 1) ^ !direction);
          clause.push_back((W[bn - 1] << 1) ^ direction);
        } else if (b == beginB && a > beginA) {
          clause.push_back((_outputs[a - 1] << 1) ^ !direction);
          clause.push_back((W[an - 1] << 1) ^ direction);
        }
        _dgpw->AddClause(clause);
        clause.clear();
      }

      if (a == endA && b == endB) {
        continue;
      }
      //            //To assure to get 0's at the beginning of A and B.
      //            if (a < endA && b < endB)
      //                std::cout << std::endl << "( " << a << ", " << b << ", "
      //                << an + bn << " ), ";
      //            else if (a == endA && b < endB)
      //                std::cout << std::endl << "(  , " << b << ", " << an +
      //                bn << " ), ";
      //            else if (b == endB && a < endA)
      //                std::cout << std::endl << "( " << a << ",  , " << an +
      //                bn << " ), ";

      if (a < endA && b < endB) {
        clause.push_back((_outputs[a] << 1) ^ direction);
        clause.push_back((_outputs[b] << 1) ^ direction);
        clause.push_back((W[an + bn] << 1) ^ !direction);
      } else if (a == endA && b < endB) {
        clause.push_back((_outputs[b] << 1) ^ direction);
        clause.push_back((W[an + bn] << 1) ^ !direction);
      } else if (b == endB && a < endA) {
        clause.push_back((_outputs[a] << 1) ^ direction);
        clause.push_back((W[an + bn] << 1) ^ !direction);
      }
      _dgpw->AddClause(clause);
      clause.clear();
    }
  }

  for (uint32_t a = beginA; a < endB; ++a) {
    assert((a - beginA) < W.size());
    _outputs[a] = W[a - beginA];
  }
  //    if (beginA == 0 && endA==5 && endB==10)
  //    {
  //        std::cout << std::endl;
  //        exit(0);
  //    }
}

void Sorter::CreateTotalizerEncodeTree() {
  if (_outputTree != nullptr)
    return;

  _outputTree = new TotalizerEncodeTree(static_cast<uint32_t>(_outputs.size()));
  // give boundaries!

  _outputTree->CreateOutputTreeReturnMaxDepth(
      0, static_cast<uint32_t>(_outputs.size()), &_outputs);
}

uint32_t Sorter::TotalizerEncodeOnes(TotalizerEncodeTree *tree,
                                     uint32_t outputIndex, uint32_t outputVar) {
  //    std::cout << __PRETTY_FUNCTION__ << ", " <<  tree->_size << ", " <<
  //    tree->_depth << ", " << outputIndex << std::endl;
  // TOBI: really to be bigger equal and not bigger?
  assert(tree->_size >= 1);

  bool direction(true);
  
  if (outputVar == 0) {
    outputVar = _dgpw->NewVariable();
    
    uint32_t countingLit = (outputVar << 1) ^ 1;
    if (tree->_exponent != UINT32_MAX && tree->_exponent != 0) {
      // case we are in the bottom bucket
      // we need all soft clause relaxation literals
      // we need all tare variables
      std::vector<uint32_t> litsC;
      std::vector<uint64_t> wghtsC;
      tree->GetAllLeavesAndWeights(litsC, wghtsC, tree->_exponent);
      // for (int index = tree->_exponent; index >= 0; index--) {
      //   auto softClauses = *_dgpw->_mainCascade->_structure[(unsigned)index]->_softClauses;
      //   for (auto softclause : softClauses) {
      //     wghtsC.push_back(1 << index);
      //     litsC.push_back(softclause->relaxationLit ^ 1);
      //   }
      //   if (!_dgpw->_mainCascade->_structure[(unsigned)index]->_isLastBucket) {
      //     wghtsC.push_back(1 << index);
      //     auto tares = _dgpw->_mainCascade->_structure[(unsigned)index]->_tares;
      //     litsC.push_back((tares[0] << 1) ^ 1);
      //   }
      // }
      _dgpw->_mainCascade->vPL->write_comment("reification of bottom bucket EncodeOnes Variable: " + std::to_string(outputVar));
      _dgpw->_mainCascade->vPL->reificationLiteralRightImpl(
          countingLit, litsC, wghtsC, (outputIndex + 1) * (1ULL << tree->_exponent), true);
      _dgpw->_mainCascade->vPL->reificationLiteralLeftImpl(
          countingLit, litsC, wghtsC, (outputIndex + 1) * (1ULL << tree->_exponent),
          true);
    } else {
      // case we are in the top bucket
      std::vector<uint32_t> leaves;
      tree->GetAllLeaves(leaves);
      _dgpw->_mainCascade->vPL->write_comment("reification of top bucket EncodeOnes Variable: " + std::to_string(outputVar));
      _dgpw->_mainCascade->vPL->reificationLiteralRightImpl(
          countingLit, leaves, outputIndex + 1, true);
      _dgpw->_mainCascade->vPL->reificationLiteralLeftImpl(
          countingLit, leaves, outputIndex + 1,
          true);
    }
  }

  uint32_t sizeA = tree->_child1->_size / tree->_child1->_everyNthOutput;
  uint32_t sizeB = tree->_child2->_size / tree->_child2->_everyNthOutput;
  ;
  uint32_t beginA = 0, beginB = 0;

  //    sizeA = 5;
  //    sizeB = 5;
  //    outputIndex = 0;

  //    std::cout << "TreeSize: " << tree->_size << "  sizeA: " << sizeA << "
  //    sizeB: " << sizeB << std::endl;

  std::vector<uint32_t> clause;
  // to assure ones at the ending.
  uint32_t beginIndex = (outputIndex < sizeB) ? 0 : outputIndex - sizeB + 1;
  uint32_t endIndex = (outputIndex < sizeA) ? outputIndex + 1 : sizeA;
  uint32_t BIndexHelper = (beginIndex == 0) ? outputIndex : sizeB - 1;
  //    std::cout << std::endl << "  beginIndex: " << beginIndex << "  endIndex:
  //    " << endIndex << "  BIndexHelper: " << BIndexHelper << "   " <<
  //    std::endl; std::cout << "1's           size A: "<< sizeA << "    size B:
  //    " << sizeB << "     outputIndex: " << outputIndex << std::endl;
  for (uint32_t index = beginIndex; index <= endIndex; index++) {
    uint32_t a = beginA + index - 1;
    uint32_t b = BIndexHelper + beginIndex - index;

           if (a != beginA - 1)
               std::cout << "( " << a << ", ";
           else
               std::cout << "(  , ";
           if (b != beginB - 1)
               std::cout << b << ", ";
           else
               std::cout << " , ";
           std::cout << outputIndex << " ), ";
           std::cout << " 1's" << std::endl;

    ////        if (_dgpw->_featureTest)
    clause.push_back((outputVar << 1) ^ direction);


    uint32_t vara = 0;
    uint32_t varb = 0;

    if (a != beginA - 1) {
      vara = tree->_child1->ReturnOutputEncodeIfNecessary(a, this, true);
      clause.push_back( ((vara << 1) ^ !direction) );
    }
    if (b != beginB - 1) {
      varb = tree->_child2->ReturnOutputEncodeIfNecessary(b, this, true);
      clause.push_back( ((varb << 1) ^ !direction) );
    }

    ////        if (!_dgpw->_featureTest)
    ////            clause.push_back((outputVar << 1) ^ direction);

    if (tree->_exponent != UINT32_MAX and tree->_exponent != 0)
      _dgpw->_mainCascade->vPL->write_comment("we are in the bottom EncodeOnes bucket");
    else
      _dgpw->_mainCascade->vPL->write_comment("we are in the top EncodeOnes bucket");

    _dgpw->_mainCascade->vPL->write_comment("clause EncodeOnes for PW Encoding");
    
    write_vPBproof_dgpwclause(outputVar, vara, varb, a, sizeA, b, sizeB, tree, clause, _dgpw->_mainCascade->vPL, true);
    //_dgpw->_mainCascade->vPL->unchecked_assumption(clause);
    _dgpw->AddClause(clause);
    clause.clear();
  }

  return outputVar;
}

uint32_t Sorter::TotalizerEncodeOutput(TotalizerEncodeTree *tree,
                                       uint32_t outputIndex) {
  //    std::cout << __PRETTY_FUNCTION__ << ", " <<  tree->_size << ", " <<
  //    tree->_depth << ", " << outputIndex << std::endl; std::cout <<
  //    "Tree->child1.size; Tree->child2.size: " << tree->_child1->_size << ", "
  //    << tree->_child1->_size << std::endl;
  // TOBI: really to be bigger equal and not bigger?
  assert(tree->_size >= 1);

  bool direction(true);
  uint32_t outputVar = _dgpw->NewVariable();
  uint32_t countingLit = (outputVar << 1) ^ 1;

  if (tree->_exponent != UINT32_MAX and tree->_exponent != 0) {
    // case we are in the bottom bucket
    // we need all soft clause relaxation literals
    // we need all tare variables
    std::vector<uint32_t> litsC;
    std::vector<uint64_t> wghtsC;
    tree->GetAllLeavesAndWeights(litsC, wghtsC, tree->_exponent);
    // std::cout << "Var, lit: " << outputVar << ", " << countingLit << std::endl;
    // for (unsigned int i = 0; i < litsC.size(); ++i) {
    //   std::cout << "Lit, Weight: " << litsC[i] << ", " << wghtsC[i] << std::endl;
    // }
    // for (int index = tree->_exponent; index >= 0; index--) {
    //   auto softClauses = *_dgpw->_mainCascade->_structure[(unsigned)index]->_softClauses;
    //   for (auto softclause : softClauses) {
    //     wghtsC.push_back(1 << index);
    //     litsC.push_back(softclause->relaxationLit ^ 1);
    //   }
    //   if (!_dgpw->_mainCascade->_structure[(unsigned)index]->_isLastBucket) {
    //     wghtsC.push_back(1 << index);
    //     auto tares = _dgpw->_mainCascade->_structure[(unsigned)index]->_tares;
    //     litsC.push_back((tares[0] << 1) ^ 1);
    //   }
    // }
    _dgpw->_mainCascade->vPL->write_comment("reification of bottom bucket EncodeZeros Variable: " + std::to_string(outputVar));
    _dgpw->_mainCascade->vPL->write_comment("outputIndex = " + std::to_string(outputIndex) + " exponent = " + std::to_string(tree->_exponent));
    _dgpw->_mainCascade->vPL->reificationLiteralRightImpl(
        countingLit, litsC, wghtsC, (static_cast<uint64_t>(outputIndex) + 1) * (1ULL << static_cast<uint64_t>(tree->_exponent)), true);
    _dgpw->_mainCascade->vPL->reificationLiteralLeftImpl(
        countingLit, litsC, wghtsC, (static_cast<uint64_t>(outputIndex) + 1) * (1ULL << static_cast<uint64_t>(tree->_exponent)),
        true);
  } else {
    // case we are in the top bucket
    std::vector<uint32_t> leaves;
    tree->GetAllLeaves(leaves);
    _dgpw->_mainCascade->vPL->write_comment("reification of top bucket EncodeZeros Variable: " + std::to_string(outputVar));
    _dgpw->_mainCascade->vPL->reificationLiteralRightImpl(
        countingLit, leaves, outputIndex + 1, true);
    _dgpw->_mainCascade->vPL->reificationLiteralLeftImpl(
        countingLit, leaves, outputIndex + 1,
        true);
  }

  //    std::cout << "tree->_size: " << tree->_size << std::endl;
  //    std::cout << "child1:      " << tree->_child1->_size << std::endl;
  //    std::cout << "child1Nth:   " << tree->_child1->_everyNthOutput <<
  //    std::endl;
  uint32_t sizeA = tree->_child1->_size / tree->_child1->_everyNthOutput;
  //    std::cout << "child2: " << std::endl;
  uint32_t sizeB = tree->_child2->_size / tree->_child2->_everyNthOutput;
  ;
  uint32_t beginA = 0, beginB = 0;

  //    sizeA = 5;
  //    sizeB = 5;
  //    outputIndex = 0;

  //    std::cout << "TreeSize: " << tree->_size << "  sizeA: " << sizeA << "
  //    sizeB: " << sizeB << std::endl;

  std::vector<uint32_t> clause;
  uint32_t beginIndex = (outputIndex < sizeB) ? 0 : outputIndex - sizeB;
  uint32_t endIndex = (outputIndex < sizeA) ? outputIndex : sizeA;
  uint32_t BIndexHelper = (beginIndex == 0) ? outputIndex : sizeB;
  //    std::cout << "outputIndex: " << outputIndex << "  beginIndex: " <<
  //    beginIndex << "  endIndex: " << endIndex << "  BIndexHelper: " <<
  //    BIndexHelper << "   " << std::endl; std::cout << std::endl << "0's size
  //    A: "<< sizeA << "    size B: " << sizeB << "     outputIndex: " <<
  //    outputIndex << std::endl;
  for (uint32_t index = beginIndex; index <= endIndex; index++) {
    uint32_t a = beginA + index;
    uint32_t b = BIndexHelper + beginIndex - index;

    //        std::cout << "a,b,outInd: " << a<< ", " << b << ", " <<
    //        outputIndex << std::endl; std::cout << std::endl << std::endl; if
    //        (a != sizeA)
    //            std::cout << "( " << a << ", ";
    //        else
    //            std::cout << "(  , ";
    //        if (b != sizeB)
    //            std::cout << b << ", ";
    //        else
    //            std::cout << " , ";
    //        std::cout << outputIndex << " ), ";
    //        std::cout << " 0's" << std::endl;

    //        if (_dgpw->_featureTest)
    //        uint32_t outputv = (outputVar << 1) ^ !direction;
    //        std::cout << "NewVar: " <<  outputv;
    clause.push_back((outputVar << 1) ^ !direction);

    assert(tree != NULL);
    assert(tree->_child1 != NULL);
    assert(tree->_child2 != NULL);

    uint32_t vara = 0, varb = 0;

    if (a != sizeA) {
      //            uint32_t child1
      //            =(tree->_child1->ReturnOutputEncodeIfNecessary(a, this) <<
      //            1) ^ direction; std::cout << ", Child1: " << child1;
      vara = tree->_child1->ReturnOutputEncodeIfNecessary(a, this);
      clause.push_back(
          (vara << 1) ^
          direction);
    }

    if (b != sizeB) {
      //            uint32_t child2
      //            =(tree->_child2->ReturnOutputEncodeIfNecessary(b, this) <<
      //            1) ^ direction; std::cout << ", Child2: " << child2;
      varb = tree->_child2->ReturnOutputEncodeIfNecessary(b, this);
      clause.push_back(
          (varb << 1) ^
          direction);
    }
    //        std::cout << std::endl;

    //        std::cout << "clause: (  ";
    //        for (uint32_t i = 0; i < clause.size(); i++)
    //            std::cout << clause[i] << "  ";
    //        std::cout << ")" << std::endl;

    _dgpw->_mainCascade->vPL->write_comment("clause for PW EncodeZeros Encoding");
    // Derivation of the clause in the proof
    write_vPBproof_dgpwclause(outputVar, vara, varb, a, sizeA, b, sizeB, tree, clause, _dgpw->_mainCascade->vPL);

    _dgpw->AddClause(clause);
    clause.clear();
  }

  if (!_setting->encode01)
    return outputVar;

  // to assure ones at the ending.
  beginIndex = (outputIndex < sizeB) ? 0 : outputIndex - sizeB + 1;
  endIndex = (outputIndex < sizeA) ? outputIndex + 1 : sizeA;
  BIndexHelper = (beginIndex == 0) ? outputIndex : sizeB - 1;
  //    std::cout << std::endl << "  beginIndex: " << beginIndex << "  endIndex:
  //    " << endIndex << "  BIndexHelper: " << BIndexHelper << "   " <<
  //    std::endl; std::cout << "1's           size A: "<< sizeA << "    size B:
  //    " << sizeB << "     outputIndex: " << outputIndex << std::endl;
  for (uint32_t index = beginIndex; index <= endIndex; index++) {
    uint32_t a = beginA + index - 1;
    uint32_t b = BIndexHelper + beginIndex - index;

    //        if (a != beginA - 1)
    //            std::cout << "( " << a << ", ";
    //        else
    //            std::cout << "(  , ";
    //        if (b != beginB - 1)
    //            std::cout << b << ", ";
    //        else
    //            std::cout << " , ";
    //        std::cout << outputIndex << " ), ";
    //        std::cout << " 1's" << std::endl;

    ////        if (_dgpw->_featureTest)
    clause.push_back((outputVar << 1) ^ direction);

    if (a != beginA - 1)
      clause.push_back(
          (tree->_child1->ReturnOutputEncodeIfNecessary(a, this) << 1) ^
          !direction);
    if (b != beginB - 1)
      clause.push_back(
          (tree->_child2->ReturnOutputEncodeIfNecessary(b, this) << 1) ^
          !direction);

    ////        if (!_dgpw->_featureTest)
    ////            clause.push_back((outputVar << 1) ^ direction);

    _dgpw->AddClause(clause);
    clause.clear();
  }
  return outputVar;
}

constraintid Sorter::write_vPBproof_dgpwclause(uint32_t outputVar, uint32_t vara, uint32_t varb, uint32_t a, uint32_t sizeA, uint32_t b, uint32_t sizeB, TotalizerEncodeTree* tree, std::vector<uint32_t>& clause, VeriPbProofLogger* vPL, bool encodeOnes){
    
  cuttingplanes_derivation cpder = vPL->CP_constraintid(
                                  encodeOnes ?  vPL->getReifiedConstraintLeftImpl(outputVar) :
                                                vPL->getReifiedConstraintRightImpl(outputVar));

  std::string commentTares = "Tares Child1: ";
  for(int i = 0; i < tree->_child1->_tares.size(); i++) commentTares += " " + vPL->var_name(tree->_child1->_tares[i]);
  commentTares += " Tares Child2: ";
  for(int i = 0; i < tree->_child2->_tares.size(); i++) commentTares += " " + vPL->var_name(tree->_child2->_tares[i]);
  vPL->write_comment(commentTares);
  vPL->write_comment(" isBottomBucket = " + std::to_string(tree->_isBottomBucket) + " exponent = " + std::to_string(tree->_exponent));
  vPL->write_comment("vara = " + vPL->var_name(vara) + " a = " + std::to_string(a) + " sizeA = " + std::to_string(sizeA)  + " child1 exponent = " + std::to_string(tree->_child1->_exponent) + " child1 isbottombucket = " + std::to_string(tree->_child1->_isBottomBucket));
  vPL->write_comment("varb = " + vPL->var_name(varb) + " b = " + std::to_string(b) + " sizeB = " + std::to_string(sizeB)  + " child2 exponent = " + std::to_string(tree->_child2->_exponent) + " child2 isbottombucket = " + std::to_string(tree->_child2->_isBottomBucket));
  
  write_vPBproof_for_child_dgpwclause(cpder, vara, a, sizeA, tree->_isBottomBucket, tree->_exponent, tree->_child1, vPL, encodeOnes);
  write_vPBproof_for_child_dgpwclause(cpder, varb, b, sizeB, tree->_isBottomBucket, tree->_exponent, tree->_child2, vPL, encodeOnes);
  
  if(!encodeOnes && (tree->_child1->_isBottomBucket && a == sizeA || tree->_child2->_isBottomBucket && b == sizeB))
    cpder = vPL->CP_division(cpder, wght_max);
  else
    cpder = vPL->CP_saturation(cpder);
  
  constraintid c =  vPL->write_CP_derivation(cpder);

  vPL->check_last_constraint(clause);
  
  return c;
}

void Sorter::write_vPBproof_for_child_dgpwclause(cuttingplanes_derivation& cpder, uint32_t var, uint32_t index, uint32_t  size, bool bottombucket, uint32_t exp, TotalizerEncodeTree* child, VeriPbProofLogger* vPL, bool encodeOnes){
  if(index == (encodeOnes ? -1ULL : size)){
    uint64_t mult = 1; 
    
    for(int i = 0; i < child->_leaves.size(); i++){
      uint32_t leaf = neg(child->_leaves[i]);

      if(bottombucket){
        // Assumption: root node top bucket for 2^0 has _isBottomBucket true.
        if(child->_isBottomBucket){
          mult = child->_leavesWeights[i];
        }
        else{
          mult = 1ULL << exp;
        }
      }

      // Assumption: leaf is positive relaxation literal. 
      cpder = vPL->CP_weakening(cpder, leaf, mult);
    }

    for(int i = 0; i < child->_tares.size(); i++){
      uint32_t tarelit = create_literal(child->_tares[i], false);

      if(bottombucket){
        // Assumption: root node top bucket for 2^0 has _isBottomBucket true.
        if(child->_isBottomBucket){
          mult = 1ULL << ((child->_tares.size()-1) - i);
        }
        else{
          mult = 1ULL << exp;
        }
      }

      cpder = vPL->CP_weakening(cpder, tarelit, mult);
    }
  }
  else if(child->_encodedOutputs.size() > 1){
    constraintid cxntoadd = encodeOnes ? vPL->getReifiedConstraintRightImpl(var) : vPL->getReifiedConstraintLeftImpl(var);
    if(bottombucket && !child->_isBottomBucket){
      // Assumption: root node top bucket for 2^0 has _isBottomBucket true.
      uint64_t mult = 1ULL << exp;
      cpder = vPL->CP_addition(cpder, 
                            vPL->CP_multiplication(vPL->CP_constraintid(cxntoadd), mult));
    }
    else{
      cpder = vPL->CP_addition(cpder, 
                            vPL->CP_constraintid(cxntoadd));
    }
  }
}


void Sorter::AddClausesToIndex(bool direction, uint32_t outputInd,
                               uint32_t sizeA, uint32_t beginA, uint32_t sizeB,
                               uint32_t endA, uint32_t endB, uint32_t outputVar,
                               uint32_t beginB) {
  std::vector<uint32_t> clause;
  uint32_t beginIndex = (outputInd < sizeB) ? 0 : outputInd - sizeB;
  uint32_t endIndex = (outputInd < sizeA) ? outputInd : sizeA;
  uint32_t BIndexHelper = (beginIndex == 0) ? outputInd + endA : endB;
  // std::cout << "outputInd: " << outputInd << "  beginIndex: " << beginIndex
  // << "  endIndex: " << endIndex << "  BIndexHelper: " << BIndexHelper << " "
  // << std::endl;
  for (uint32_t index = beginIndex; index <= endIndex; index++) {
    uint32_t a = beginA + index;
    uint32_t b = BIndexHelper + beginIndex - index;

    //        if (a != endA)
    //            std::cout << "( " << a << ", ";
    //        else
    //            std::cout << "(  , ";
    //        if (b != endB)
    //            std::cout << b << ", ";
    //        else
    //            std::cout << " , ";
    //        std::cout << outputInd << " ), ";
    //        std::cout << std::endl;

    clause.push_back((outputVar << 1) ^ !direction);

    if (a != endA)
      clause.push_back((_outputs[a] << 1) ^ direction);
    if (b != endB)
      clause.push_back((_outputs[b] << 1) ^ direction);

    // std::csshout << "clause.size(): " << clause.size() << std::endl;
    _dgpw->AddClause(clause);
    clause.clear();
  }

  if (!_setting->encode01)
    return;

  // to assure ones at the ending.
  beginIndex = (outputInd < sizeB) ? 0 : outputInd - sizeB + 1;
  endIndex = (outputInd < sizeA) ? outputInd + 1 : sizeA;
  BIndexHelper = (beginIndex == 0) ? outputInd + endA : endB - 1;
  // std::cout << "outputInd: " << outputInd << "  beginIndex: " << beginIndex
  // << "  endIndex: " << endIndex << "  BIndexHelper: " << BIndexHelper01 << "
  // " << std::endl;
  //    std::cout << std::endl << "1's           size A: "<< sizeA << "    size
  //    B: " << sizeB << std::endl;
  for (uint32_t index = beginIndex; index <= endIndex; index++) {
    uint32_t a = beginA + index - 1;
    uint32_t b = BIndexHelper + beginIndex - index;

    //        if (a != beginA - 1)
    //            std::cout << "( " << a << ", ";
    //        else
    //            std::cout << "(  , ";
    //        if (b != beginB - 1)
    //            std::cout << b << ", ";
    //        else
    //            std::cout << " , ";
    //        std::cout << outputInd << " ), ";
    //        std::cout << std::endl;

    clause.push_back((outputVar << 1) ^ direction);

    if (a != beginA - 1)
      clause.push_back((_outputs[a] << 1) ^ !direction);
    if (b != beginB - 1)
      clause.push_back((_outputs[b] << 1) ^ !direction);

    _dgpw->AddClause(clause);
    clause.clear();
  }
}

void Sorter::EncodingUVDiagonal(uint32_t beginA, uint32_t endA, uint32_t endB,
                                std::vector<uint32_t> *outputList) {
  // Variables----------------------------------------------------
  // if direction is true 0's at the beginning, 1's at the end
  // if direction is false the other way around.
  bool direction(true);
  // std::cout << "beginA: " << beginA << "   endA: " << endA << "   endB: " <<
  // endB << "   outputs.size(): " << _outputs.size() << std::endl;
  assert(beginA <= endA);
  assert(endA < endB);
  assert(endB <= _outputs.size());
  // -------------------------------------------------------------

  //    std::vector<uint32_t> outputList(endB - beginA);
  //    std::int32_t n(0);
  //    std::generate(std::begin(outputList), std::end(outputList), [&]{ return
  //    n++; });

  uint32_t sizeA = endA - beginA;
  uint32_t sizeB = endB - endA;
  uint32_t beginB = endA;

  std::vector<uint32_t> W;

  // for every output
  //    std::cout << std::endl << std::endl;
  //    std::cout << "_outputs.size(): " << _outputs.size() << std::endl;
  //    std::cout << "sizeA: " << sizeA << "  beginA: " << beginA << "  endA: "
  //    << endA << std::endl; std::cout << "sizeB: " << sizeB << "  endB: " <<
  //    endB << std::endl;
  if (outputList == nullptr) {
    for (uint32_t outputInd = 0; outputInd < endB - beginA; outputInd++) {
      W.push_back(_dgpw->NewVariable());
      // to assure zeros at the beginning.
      AddClausesToIndex(direction, outputInd, sizeA, beginA, sizeB, endA, endB,
                        W[outputInd], beginB);
    }
  } else {
    for (auto outputInd : *outputList) {
      W.push_back(_dgpw->NewVariable());
      // to assure zeros at the beginning.
      AddClausesToIndex(direction, outputInd, sizeA, beginA, sizeB, endA, endB,
                        W.back(), beginB);
    }
  }

  //    for (uint32_t a = beginA; a < endB; ++a)
  //    {
  //        assert( (a-beginA) < W.size() );
  //        _outputs[a] = W[a-beginA];
  //    }

  if (outputList == NULL) {
    for (uint32_t outputInd = 0; outputInd < endB - beginA; outputInd++)
    // for (uint32_t a = beginA; a < endB; ++a)
    {
      // assert( (a-beginA) < W.size() );
      _outputs[outputInd + beginA] = W[outputInd];
    }
  } else {
    for (auto outputInd : *outputList) {
      _outputs[outputInd + beginA] = W[outputInd];
    }
  }
}
} // namespace DGPW
} // namespace Pacose
