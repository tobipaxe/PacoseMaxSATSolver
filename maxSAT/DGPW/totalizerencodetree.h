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

#ifndef TOTALIZERENCODETREE_H
#define TOTALIZERENCODETREE_H

#include "sorter.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

namespace Pacose {
namespace DGPW {

struct TotalizerEncodeTree {
  TotalizerEncodeTree(uint32_t size)
      : _encodedOutputs(size, 0), _leaves(), _leavesWeights(), _tares(),
        _size(size), _depth(0), _howOftenUsed(0), _maxPos(0),
        _allOutputsEncoded(false), _hasBeenBucketBefore(false),
        _onesEncoded(false), _isBottomBucket(true), _everyNthOutput(1),
        _exponent(UINT32_MAX), _verbosity(0), _child1(nullptr),
        _child2(nullptr) {}

  ~TotalizerEncodeTree() {
    delete _child1;
    delete _child2;
  }

  // the inputs are the _encodedOutputs of the children.
  std::vector<uint32_t> _encodedOutputs;
  // PL: the inputs are all leaves
  std::vector<uint32_t> _leaves;
  std::vector<uint64_t> _leavesWeights;
  std::vector<uint32_t> _tares;

  uint32_t _size;
  uint32_t _depth;
  uint32_t _howOftenUsed;
  uint32_t _maxPos;
  uint32_t _exponent;
  uint32_t _verbosity;

  bool _allOutputsEncoded;
  bool _hasBeenBucketBefore;
  bool _onesEncoded;
  bool _isBottomBucket;

  // as standard every output counts.
  uint32_t _everyNthOutput;

  TotalizerEncodeTree *_child1;
  TotalizerEncodeTree *_child2;

  void CutVectorAbove(uint32_t maxPos) {
    assert(_size == _encodedOutputs.size());
    if (_verbosity > 1)
      std::cout << std::endl
                << "Size before Vector Cut ABOVE EncodeTree: "
                << _encodedOutputs.size() << "  SIZE: " << _size << std::endl;
    _maxPos = maxPos;
    while (_size > maxPos) {
      _encodedOutputs.pop_back();
      _size--;
    }
    if (_verbosity > 1)
      std::cout << "Size after Vector Cut ABOVE EncodeTree: "
                << _encodedOutputs.size() << "  SIZE: " << _size << std::endl
                << std::endl;
  }

  void CutVectorBelow(uint32_t minPos) {
    assert(_size == _encodedOutputs.size());

    if (_verbosity > 1)
      std::cout << std::endl
                << "Size before Vector Cut BELOW EncodeTree: "
                << _encodedOutputs.size() << "  SIZE: " << _size << std::endl;
    if (_size > minPos)
      _encodedOutputs.erase(_encodedOutputs.begin(),
                            _encodedOutputs.begin() + minPos);

    _size = _encodedOutputs.size();

    if (_verbosity > 1)
      std::cout << "Size after Vector Cut BELOW EncodeTree: "
                << _encodedOutputs.size() << "  SIZE: " << _size << std::endl
                << std::endl;
  }

  uint32_t ReturnOutputEncodeIfNecessary(uint32_t index, Sorter *sorter,
                                         bool encodeOnlyOnes = false) {
    //        std::cout << __PRETTY_FUNCTION__ << std::endl;
    //        std::cout << std::endl << "NodeSize: " << _size;
    //        if (_everyNthOutput > 1)
    //            std::cout << "/" << _everyNthOutput;
    //        else if (_hasBeenBucketBefore)
    //            std::cout << "*" << _howOftenUsed;
    //        std::cout << "  Encoded Outputs:  (  ";
    //        for (uint32_t i = 0; i < _encodedOutputs.size(); i++)
    //            std::cout << _encodedOutputs[i] << "  ";
    //        std::cout << ")" << std::endl;

    //        std::cout << "requested Index: " << index;
    //        std::cout << "   _everyNthOutput; " << _everyNthOutput << " _size:
    //        " << _size
    //                  << "   index: " << index << std::endl;

    index = (_everyNthOutput * (index + 1)) - 1;

    // TOBI: WHY??
    //  std::cout << "size: " << _size << "  index: " << index
    //            << "  encodedOutputs.size(): " << _encodedOutputs.size()
    //            << std::endl;
    assert(_size > index);
    assert(_size == _encodedOutputs.size());

    if (encodeOnlyOnes && !_onesEncoded && _size != 1) {
      _encodedOutputs[index] =
          sorter->TotalizerEncodeOnes(this, index, _encodedOutputs[index]);
      _onesEncoded = true;
    } else if (_encodedOutputs[index] != 0) {
      return _encodedOutputs[index];
    } else {
      _encodedOutputs[index] = sorter->TotalizerEncodeOutput(this, index);
    }

    //        std::cout << "                                 One Index more
    //        Encoded!! NodeLabel:            " << _size; if(_everyNthOutput >
    //        1)
    //            std::cout << "/" << _everyNthOutput;
    //        else if (_hasBeenBucketBefore)
    //            std::cout << "*" << _howOftenUsed;
    //        std::cout << "  Encoded Outputs:  (  ";
    //        for (uint32_t i = 0; i < _encodedOutputs.size(); i++)
    //            std::cout << _encodedOutputs[i] << "  ";
    //        std::cout << ")" << std::endl;

    return _encodedOutputs[index];
  }

  /**
   * @brief ActualizeValues
   *          After connecting two subbuckets, the values have to been
   * actualized.
   */
  void ActualizeValues() {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    // std::cout << "size() " << _encodedOutputs.size() << std::endl;
    // std::cout << "_child1->_size " << _child1->_size << std::endl;
    // std::cout << "_child2->_size " << _child2->_size << std::endl;
    _size = _child1->_size / _child1->_everyNthOutput;
    _size += _child2->_size / _child2->_everyNthOutput;
    // std::cout << "size() " << _encodedOutputs.size() << std::endl;

    //        std::cout << "child1: " << _child1->_size /
    //        _child1->_everyNthOutput << std::endl; std::cout << "child2: " <<
    //        _child2->_size / _child2->_everyNthOutput << std::endl; std::cout
    //        << "ACTUALIZE VALUES SIZE: " << _size << std::endl;

    _depth = (_child1->_depth > _child2->_depth) ? _child1->_depth + 1
                                                 : _child2->_depth + 1;
    _encodedOutputs.resize(_size, 0);
  }

  void combineLeaves() {
    if (_leaves.size() == 0) {
      if (_verbosity > 10) {
        std::cout << "Only leaves before resizing: ";
        for (uint32_t i = 0; i < _leaves.size(); i++) {
          std::cout << _leaves[i] << ", ";
        }
        std::cout << std::endl;
      }
      // _leaves.resize(_child1->_leaves.size() + _child2->_leaves.size()); //
      // Resize _leaves vector before inserting elements std::cout << "resized
      // to: " << _leaves.size() << std::endl; std::cout << "Only leaves without
      // content: "; for (uint32_t i = 0; i < _leaves.size(); i++) {
      //   std::cout << _leaves[i] << ", ";
      // }
      // std::cout << std::endl;

      if (_verbosity > 10) {
        std::cout << "Only leaves from child1: ";
        for (uint32_t i = 0; i < _child1->_leaves.size(); i++) {
          std::cout << _child1->_leaves[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "Only leaves from child2: ";
        for (uint32_t i = 0; i < _child2->_leaves.size(); i++) {
          std::cout << _child2->_leaves[i] << ", ";
        }
        std::cout << std::endl;

        std::cout << "child1/2 enO: " << _child1->_everyNthOutput << ", "
                  << _child2->_everyNthOutput << std::endl;
      }
      // assert(_child1->_everyNthOutput != _child2->_everyNthOutput);
      for (uint32_t i = 0; i < _child1->_leaves.size(); i++) {
        _leaves.push_back(_child1->_leaves[i]);
        if (_child1->_everyNthOutput > 1)
          _leavesWeights.push_back(2);
        else
          _leavesWeights.push_back(1);
      }
      for (uint32_t i = 0; i < _child2->_leaves.size(); i++) {
        _leaves.push_back(_child2->_leaves[i]);
        if (_child2->_everyNthOutput > 1)
          _leavesWeights.push_back(2);
        else
          _leavesWeights.push_back(1);
      }
      if (_verbosity > 10) {
        std::cout << "_child1->_leaves.size() = " << _child1->_leaves.size()
                  << std::endl;
        std::cout << "_child2->_leaves.size() = " << _child2->_leaves.size()
                  << std::endl;
        std::cout << "_leaves.size() = " << _leaves.size() << std::endl;
        std::cout << "_leavesWeights.size() = " << _leavesWeights.size()
                  << std::endl;

        // assert(_leaves.size() == _child1->_leaves.size() +
        // _child2->_leaves.size()); assert(_leaves.size() ==
        // _leavesWeights.size());

        std::cout << "Combined leaves and weights: ";
        for (uint32_t i = 0; i < _leaves.size(); i++) {
          std::cout << _leaves[i] << "(" << _leavesWeights[i] << "), ";
        }
        std::cout << std::endl;
      }
    }
  }

  void CombineLeavesForBottomBucket(TotalizerEncodeTree *bb,
                                    TotalizerEncodeTree *tb) {
    if (_verbosity > 10)
      std::cout << __PRETTY_FUNCTION__ << std::endl;
    // Note that bb's leaves are already sorted by the recursive call!
    std::sort(tb->_leaves.begin(), tb->_leaves.end());

    uint32_t i = 0, j = 0;

    _leaves.reserve(bb->_leaves.size() + tb->_leaves.size());
    if (_verbosity > 10)
      std::cout << "Before while loop..." << std::endl;

    while (i < tb->_leaves.size() &&
           j < bb->_leaves
                   .size()) { // Literal in both in bottom bucket and top bucket
      if (tb->_leaves[i] == bb->_leaves[j]) {
        _leaves.push_back(tb->_leaves[i]);
        _leavesWeights.push_back(bb->_leavesWeights[j] + (1ULL << _exponent));
        i++;
        j++;
      } else if (tb->_leaves[i] <
                 bb->_leaves[j]) { // Literal only in top bucket
        _leaves.push_back(tb->_leaves[i]);
        _leavesWeights.push_back(1ULL << _exponent);
        i++;
      } else { // Literal only in bottom bucket
        _leaves.push_back(bb->_leaves[j]);
        _leavesWeights.push_back(bb->_leavesWeights[j]);
        j++;
      }
    }
    if (_verbosity > 10)
      std::cout << "after (almost) combining the two vectors..." << std::endl;

    while (i < tb->_leaves.size()) {
      _leaves.push_back(tb->_leaves[i]);
      _leavesWeights.push_back(1ULL << _exponent);
      i++;
    }
    if (_verbosity > 10) {
      std::cout << "after (almost) combining the two vectors... i, j: " << i
                << ", " << j << std::endl;
      std::cout << "BB leaves size :" << bb->_leaves.size()
                << " TB leaves size : " << tb->_leaves.size()
                << " BB leavesWeights size: " << bb->_leavesWeights.size()
                << " TB leavesWeights size: " << tb->_leavesWeights.size()
                << " current exponent: " << _exponent
                << " current expValue: " << (1ULL << _exponent)
                << std::endl;
      std::cout << "BB Leaves:";

      for (int i = 0; i < bb->_leaves.size(); i++) {
        std::cout << " " << bb->_leavesWeights[i] << " lit(" << bb->_leaves[i]
                  << ")";
      }
      std::cout << std::endl;
      std::cout << "TB Leaves:";
      for (int i = 0; i < tb->_leaves.size(); i++) {
        std::cout << " "
                  << (1ULL << _exponent)
                  << " lit(" << tb->_leaves[i] << ")";
      }
      std::cout << std::endl;
    }

    while (j < bb->_leaves.size()) {
      _leaves.push_back(bb->_leaves[j]);
      _leavesWeights.push_back(bb->_leavesWeights[j]);
      j++;
    }
    if (_verbosity > 10) {
      std::cout << "after combining the two vectors...2" << std::endl;

      std::cout << "Leaves:";
      for (int i = 0; i < _leaves.size(); i++) {
        std::cout << " " << _leavesWeights[i] << " lit(" << _leaves[i] << ")";
      }
      std::cout << std::endl;
    }
  }

  void AddBookkeepingForPL(bool isBottomBucket, uint32_t lowTare) {
    _isBottomBucket = isBottomBucket;
    if (_verbosity > 10) {
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cout << "Low tare: " << lowTare << std::endl;
      std::cout << "Current size: " << _size << std::endl;
      if (_size == 1) {
        std::cout << "Leave: " << _encodedOutputs[0] << std::endl;
      }
    }

    if (_isBottomBucket && _child1 && _child1->_everyNthOutput == 2) {

      if (_child1->_leaves.size() == 0)
        _child1->AddBookkeepingForPL(true, lowTare);
      if (_child2->_leaves.size() == 0)
        _child2->AddBookkeepingForPL(false, lowTare);

      if (_verbosity > 10) {
        std::cout << "Visiting node " << this << std::endl;
        std::cout << "Child 1 is bottom bucket. " << std::endl;
        std::cout << "_isBottomBucket = " << _isBottomBucket
                  << " _child1->_isBottomBucket = " << _child1->_isBottomBucket
                  << "_child2->_isBottomBucket = " << _child2->_isBottomBucket
                  << std::endl;
      }

      if (_exponent == UINT32_MAX)
        _exponent = _child1->_exponent + 1;

      CombineLeavesForBottomBucket(_child1, _child2);

    } else if (_isBottomBucket && _child2 && _child2->_everyNthOutput == 2) {

      if (_child1->_leaves.size() == 0)
        _child1->AddBookkeepingForPL(false, lowTare);
      if (_child2->_leaves.size() == 0)
        _child2->AddBookkeepingForPL(true, lowTare);

      if (_verbosity > 10) {
        std::cout << "Visiting node " << this << std::endl;
        std::cout << "Child 2 is bottom bucket. " << std::endl;
        std::cout << "_isBottomBucket = " << _isBottomBucket
                  << " _child1->_isBottomBucket = " << _child1->_isBottomBucket
                  << "_child2->_isBottomBucket = " << _child2->_isBottomBucket
                  << std::endl;
      }
      if (_exponent == UINT32_MAX)
        _exponent = _child2->_exponent + 1;

      CombineLeavesForBottomBucket(_child2, _child1);

    } else {

      if (_child1 && _child1->_leaves.size() == 0)
        _child1->AddBookkeepingForPL(false, lowTare);
      if (_child2 && _child2->_leaves.size() == 0)
        _child2->AddBookkeepingForPL(false, lowTare);

      if (_verbosity > 10) {
        std::cout << "Visiting node " << this << std::endl;
        std::cout << "We are in a top bucket" << std::endl;
      }

      if (_child1) {
        if (_verbosity > 10) {
          std::cout << "child1: " << _child1 << std::endl;
          std::cout << "Leaves child1: ";
        }
        for (auto leaf : _child1->_leaves) {
          if (_verbosity > 10)
            std::cout << leaf << ", ";
          _leaves.push_back(leaf);
        }
        if (_verbosity > 10)
          std::cout << std::endl;
      }
      if (_child2) {
        if (_verbosity > 10) {
          std::cout << "child2: " << _child2 << std::endl;
          std::cout << "Leaves child2: ";
        }
        for (auto leaf : _child2->_leaves) {
          if (_verbosity > 10)
            std::cout << leaf << ", ";
          _leaves.push_back(leaf);
        }
        if (_verbosity > 10)
          std::cout << std::endl;
      }

      if (!_child1 && !_child2) { // We are in a leaf!
        if (_verbosity > 10) {
          std::cout << "Indeed, we are in a leaf. " << std::endl;
          assert(_size == 1);
          std::cout << "_leaves.size() = " << _leaves.size() << std::endl;
        }

        if (_encodedOutputs[0] < lowTare)
          _leaves.push_back(_encodedOutputs[0] << 1 ^ 1);
        else
          _tares.push_back(_encodedOutputs[0]);
        if (_verbosity > 10) {
          std::cout << "Leaves: ";
          for (auto leaf : _leaves) {
            std::cout << leaf << ", ";
          }
          std::cout << std::endl;
        }
      }

      if (_isBottomBucket) {
        if (_verbosity > 10)
          std::cout << "We are at the 2^0 Bottom Bucket. " << std::endl;
        // then we are in the top bottom bucket with exponent 0
        std::sort(_leaves.begin(), _leaves.end());
        if (_exponent == UINT32_MAX)
          _exponent = 0;
        for (int i = 0; i < _leaves.size(); i++) {
          _leavesWeights.push_back(1ULL << _exponent);
        }
      }
    }

    if (_encodedOutputs.size() > 1) {
      if (_verbosity > 10)
        std::cout << "Starting to copy the tares. " << std::endl;
      std::copy(_child1->_tares.begin(), _child1->_tares.end(),
                std::back_inserter(_tares));
      std::copy(_child2->_tares.begin(), _child2->_tares.end(),
                std::back_inserter(_tares));
    }
    if (_verbosity > 10)
      std::cout << "End visiting Node." << std::endl;
  }

  // create empty tree for given size()
  /**
   * @brief CreateOutputTreeReturnMaxDepth
   *                    create empty tree; leaves contain input values.
   * @param lo          Points to first element of inputVector to look at.
   *                    A hard lower bound can be given!
   * @param hi          Points to last element + 1 of inputVector to look at.
   *                    A hard upper bound can be given!
   * @param inputVector Pointer to whole given inputVector.
   * @return            Depth of tree, root has greatest depth.
   */
  uint32_t CreateOutputTreeReturnMaxDepth(uint32_t lo, uint32_t hi,
                                          std::vector<uint32_t> *inputVector) {

    if ((hi - lo) > 1) {
      assert(hi > lo);
      uint32_t m = ((hi - lo) >> 1);
      _child1 = new TotalizerEncodeTree(m);
      uint32_t depth1 = _child1->CreateOutputTreeReturnMaxDepth(
          lo, lo + m, inputVector);
      _child2 = new TotalizerEncodeTree(_size - m);
      uint32_t depth2 = _child2->CreateOutputTreeReturnMaxDepth(
          lo + m, hi, inputVector);
      _depth = (depth1 > depth2) ? depth1 : depth2;
    } else {
      assert(hi - lo == 1);
      assert(hi > lo);
      _encodedOutputs[0] = (*inputVector)[lo];
      _allOutputsEncoded = true;
      _depth = 0;
    }
    return _depth + 1;
  }

  /// @brief Calculates the _exponent variable for all everyNthOutput children
  /// recursively
  void CalculateExponents() {
    uint32_t currentExponent = CalculateMaxExponent();
    SetExponentsRecursively(currentExponent);
  }

  //    uint32_t AttachTwoChildrenReturnMaxDepth(TotalizerEncodeTree*
  //    firstChild, TotalizerEncodeTree* secondChild)
  //    {
  //        assert(firstChild->_size == 0 || secondChild->_size == 0);

  //        _child1 = firstChild;
  //        _child2 = secondChild;
  //        _depth = (firstChild->_depth > secondChild->_depth) ?
  //        firstChild->_depth + 1 : secondChild->_depth + 1;

  //        return _depth;
  //    }

  /**
   * @brief DumpOutputTree
   *          The output is an Trivial Graph format output - *.tgf!!
   *          To visualize it, use a program like yed.
   *          At first print out all existing node numbers with their label.
   *          Followed by the # symbol.
   *          Then cout all connections between the nodes.
   */
  void DumpOutputTree(std::string filename, bool labelWithOutputs = false) {
    std::ofstream out(filename);
    std::streambuf *coutbuf = std::cout.rdbuf(); // save old buf
    std::cout.rdbuf(out.rdbuf()); // redirect std::cout to out.txt!

    DumpNodes(1, labelWithOutputs);
    std::cout << "#" << std::endl;
    DumpConnections(1);
    std::cout.rdbuf(coutbuf); // reset to standard output again
    std::cout << "c " << filename << " file written." << std::endl;
  }

  /**
   * @brief DumpNodes
   *          Cout all existing node numbers with their label.
   * @param nodeNumber
   *          Call with next not yet given node number.
   * @return  Highest node subNode number.
   */
  uint32_t DumpNodes(uint32_t nodeNumber, bool labelWithOutputs) {
    if (_size == 1) {
      std::cout << nodeNumber << " " << _encodedOutputs[0] << std::endl;
      return nodeNumber;
    }
    std::cout << nodeNumber << " ";
    if (labelWithOutputs) {
      std::cout << "(";
      for (uint32_t i = 0; i < _encodedOutputs.size(); i++) {
        if (i != _encodedOutputs.size() - 1)
          std::cout << _encodedOutputs[i] << ",";
        else
          std::cout << _encodedOutputs[i] << ")";
      }
    } else {
      std::cout << "s" << _size;
      if (_everyNthOutput > 1)
        std::cout << "/nth" << _everyNthOutput;
      else if (_hasBeenBucketBefore)
        std::cout << "*" << _howOftenUsed;
      std::cout << "|BB" << _isBottomBucket << "|E" << _exponent;
      std::cout << "|L";
      for (uint32_t i = 0; i < _leaves.size(); i++) {
        std::cout << "," << _leaves[i];
      }
      std::cout << "|W";
      for (uint32_t i = 0; i < _leavesWeights.size(); i++) {
        std::cout << "," << _leavesWeights[i];
      }
      std::cout << "|T";
      for (uint32_t i = 0; i < _tares.size(); i++) {
        std::cout << "," << _tares[i];
      }
    }
    std::cout << std::endl;
    uint32_t nodeNumber1 = _child1->DumpNodes(nodeNumber + 1, labelWithOutputs);
    uint32_t nodeNumber2 =
        _child2->DumpNodes(nodeNumber1 + 1, labelWithOutputs);
    return nodeNumber2;
  }

  /**
   * @brief DumpConnections
   *          Cout all possible connections between nodes.
   * @param nodeNumber
   *          Call with next not yet given node number.
   * @return  Highest node subNode number.
   */
  uint32_t DumpConnections(uint32_t nodeNumber) {
    if (_size == 1)
      return nodeNumber;

    std::cout << nodeNumber << " " << nodeNumber + 1 << std::endl;
    uint32_t nodeNumber1 = _child1->DumpConnections(nodeNumber + 1);
    std::cout << nodeNumber << " " << nodeNumber1 + 1 << std::endl;
    uint32_t nodeNumber2 = _child2->DumpConnections(nodeNumber1 + 1);
    return nodeNumber2;
  }

  void GetAllLeaves(std::vector<uint32_t> &leaves) {
    if (_size == 1) {
      leaves.push_back((_encodedOutputs[0] << 1) ^ 1);
    } else {
      _child1->GetAllLeaves(leaves);
      _child2->GetAllLeaves(leaves);
    }
  }

  void GetAllLeavesAndWeights(std::vector<uint32_t> &leaves,
                              std::vector<uint64_t> &weights) {
    GetAllLeavesAndWeights(leaves, weights, _exponent);
  }

  void GetAllLeavesAndWeights(std::vector<uint32_t> &leaves,
                              std::vector<uint64_t> &weights,
                              uint32_t currentExponent) {
    if (_size == 1) {
      leaves.push_back((_encodedOutputs[0] << 1) ^ 1);
      weights.push_back(1ULL << currentExponent);
    } else {
      assert(_child1);
      assert(_child2);
      if (_child1->_everyNthOutput > 1)
        _child1->GetAllLeavesAndWeights(leaves, weights, currentExponent - 1);
      else
        _child1->GetAllLeavesAndWeights(leaves, weights, currentExponent);
      if (_child2->_everyNthOutput > 1)
        _child2->GetAllLeavesAndWeights(leaves, weights, currentExponent - 1);
      else
        _child2->GetAllLeavesAndWeights(leaves, weights, currentExponent);
    }
  }

private:
  // Copy constructor.
  TotalizerEncodeTree(const TotalizerEncodeTree &) = default;

  // Assignment operator.
  TotalizerEncodeTree &operator=(const TotalizerEncodeTree &) = default;

  void SetExponentsRecursively(uint32_t currentExponent) {
    assert(currentExponent != UINT32_MAX);
    if (_exponent == UINT32_MAX)
      _exponent = currentExponent;
    if (_child1 && _child1->_everyNthOutput > 1) {
      _child1->SetExponentsRecursively(_exponent - 1);
      return;
    } else if (_child2 && _child2->_everyNthOutput > 1) {
      _child2->SetExponentsRecursively(_exponent - 1);
      return;
    }
    assert(currentExponent == 0);
  }

  uint32_t CalculateMaxExponent(uint32_t maxExponent = 0) {
    assert(!(_child1 && _child1->_everyNthOutput > 1 && _child2 &&
             _child2->_everyNthOutput > 1));

    if (_child1 && _child1->_everyNthOutput > 1) {
      maxExponent = _child1->CalculateMaxExponent(maxExponent + 1);
    } else if (_child2 && _child2->_everyNthOutput > 1) {
      maxExponent = _child2->CalculateMaxExponent(maxExponent + 1);
    }
    return maxExponent;
  }
};

} // namespace DGPW
} // Namespace Pacose

#endif // TOTALIZERENCODETREE_H
