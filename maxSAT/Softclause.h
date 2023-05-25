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
 
#ifndef SOFTCLAUSE_H
#define SOFTCLAUSE_H

#include <iostream>
#include <vector>

namespace Pacose {
// For MAX-SAT we store the soft clauses separately

struct SoftClause {
  SoftClause(const std::vector<uint32_t> &clause, uint64_t weight = 1)
      : relaxationLit(0),
        positiveForcePossible(false),
        clause(clause),
        lastassignment(0),
        weight(weight),
        originalWeight(weight) {}

  // Constructor
  SoftClause(uint32_t relaxlit, const std::vector<uint32_t> &clause,
             uint64_t weight = 1)
      : relaxationLit(relaxlit),
        positiveForcePossible(false),
        clause(clause),
        lastassignment(0),
        weight(weight),
        originalWeight(weight) {}

  // Destructor
  ~SoftClause(void) {}

  uint32_t relaxationLit;
  
  // is it possible to force the clause to be true with the relaxation literal?
  bool positiveForcePossible;
  
  std::vector<uint32_t> clause;

  // Stores last found assignment for the relaxation lit
  uint32_t lastassignment;

  //  uint32_t contra;
  uint64_t weight;
  uint64_t originalWeight;

  static bool bigger(const SoftClause *s1, const SoftClause *s2) {
    return s1->weight > s2->weight;
  }
  static bool ids(const SoftClause *s1, const SoftClause *s2) {
    return s1->relaxationLit > s2->relaxationLit;
  }
};
} // Namespace Pacose


// struct SoftClauseContraSorter {
//    bool operator()(SoftClause* s1, SoftClause* s2) const;
//};

// inline bool SoftClauseContraSorter::operator()(SoftClause* s1, SoftClause*
// s2) const {
//    return s1->contra > s2->contra;
//}

#endif  // SOFTCLAUSE_H
