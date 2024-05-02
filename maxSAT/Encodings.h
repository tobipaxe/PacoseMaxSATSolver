/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson
Changes are made for QMaxSAT by the QMaxSAT team.
Updated by Tobias Paxian 2020

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
**************************************************************************************************/
#ifndef ENCODINGS_H
#define ENCODINGS_H

//#include "Solver.h"
//#include "Pacose.h"
//#include "SolverTypes.h"
#include "Settings.h"
#include "vector"


class SATSolverProxy;
namespace Pacose {

// class Pacose;

class Encodings {
 public:
  explicit Encodings(Settings *settings) : _settings(settings), _relaxLit(0) {}

  // nothing to delete, settings is destroyed by Pacose
  ~Encodings(){};

  // private:
  Settings *_settings;
  uint32_t _relaxLit;

  // Warners adder encoding
  void genWarnersHalf(uint32_t &a, uint32_t &b, uint32_t &carry, uint32_t &sum,
                      int comp, SATSolverProxy &S, std::vector<uint32_t> &lits);
  void genWarnersFull(uint32_t &a, uint32_t &b, uint32_t &c, uint32_t &carry,
                      uint32_t &sum, int comp, SATSolverProxy &S,
                      std::vector<uint32_t> &lits);
  void genWarners(std::vector<long long> &weights,
                  std::vector<uint32_t> &blockings, long long max, int k,
                  int comp, SATSolverProxy &S, const uint32_t zero,
                  std::vector<uint32_t> &lits,
                  std::vector<uint32_t> &linkingint);
  void genWarners0(std::vector<long long> &weights,
                   std::vector<uint32_t> &blockings, long long max, long long k,
                   int comp, SATSolverProxy &S, std::vector<uint32_t> &lits,
                   std::vector<uint32_t> &linkingint);

  // Bailleux totalizer encoding
  void genBailleux(std::vector<long long> &weights,
                   std::vector<uint32_t> &blockings, long long total,
                   uint32_t zero, uint32_t one, int comp, SATSolverProxy &S,
                   std::vector<uint32_t> &lits,
                   std::vector<uint32_t> &linkingint, long long UB);
  void genBailleux0(std::vector<long long> &weights,
                    std::vector<uint32_t> &blockings, long long max,
                    long long k, int comp, SATSolverProxy &S,
                    std::vector<uint32_t> &lits,
                    std::vector<uint32_t> &linkingint);

  // [Asin et. al 2011] encoding
  void sComparator(uint32_t &a, uint32_t &b, uint32_t &c1, uint32_t &c2,
                   int comp, SATSolverProxy &S, std::vector<uint32_t> &lits);
  void genSMerge(std::vector<uint32_t> &linkA, std::vector<uint32_t> &linkB,
                 uint32_t zero, uint32_t one, int comp, SATSolverProxy &S,
                 std::vector<uint32_t> &lits, std::vector<uint32_t> &linkingint,
                 long long UB);
  void genKCard(std::vector<long long> &weights,
                std::vector<uint32_t> &blockings, long long total,
                uint32_t zero, uint32_t one, int comp, SATSolverProxy &S,
                std::vector<uint32_t> &lits, std::vector<uint32_t> &linkingint,
                long long UB);
  void genAsin(std::vector<long long> &weights,
               std::vector<uint32_t> &blockings, long long max, long long k,
               int comp, SATSolverProxy &S, std::vector<uint32_t> &lits,
               std::vector<uint32_t> &linkingint);

  // Ogawa Modulo Totalizer
  void genOgawa(long long weightX, std::vector<uint32_t> &linkingX,
                long long weightY, std::vector<uint32_t> &linkingY,
                long long &total, long long divisor, uint32_t zero,
                uint32_t one, int comp, SATSolverProxy &S,
                std::vector<uint32_t> &lits, std::vector<uint32_t> &linkingint,
                long long UB);
  void genOgawa(std::vector<long long> &weights,
                std::vector<uint32_t> &blockings, long long &total,
                long long divisor, uint32_t zero, uint32_t one, int comp,
                SATSolverProxy &S, std::vector<uint32_t> &lits,
                std::vector<uint32_t> &linkingint, long long UB);
  void genOgawa0(std::vector<long long> &weights,
                 std::vector<uint32_t> &blockings, long long max, long long k,
                 long long &divisor, int comp, SATSolverProxy &S,
                 std::vector<uint32_t> &lits,
                 std::vector<uint32_t> &linkingint);

  // BailleuxW2
  void genBailleuxW2(std::vector<long long> &weights,
                     std::vector<uint32_t> &blockings, long long total,
                     uint32_t zero, uint32_t one, int comp, SATSolverProxy &S,
                     std::vector<uint32_t> &lits,
                     std::vector<uint32_t> &linkingint,
                     std::vector<long long> &linkingW, long long UB);
  void genBailleuxW20(std::vector<long long> &weights,
                      std::vector<uint32_t> &blockings, long long max,
                      long long k, int comp, SATSolverProxy &S,
                      std::vector<uint32_t> &lits,
                      std::vector<uint32_t> &linkingint,
                      std::vector<long long> &linkingWeight);
  void genCCl(uint32_t a, SATSolverProxy &S, std::vector<uint32_t> &lits,
              int varZero);
  void genCCl(uint32_t a, uint32_t b, SATSolverProxy &S,
              std::vector<uint32_t> &lits, int varZero);
  void genCCl(uint32_t a, uint32_t b, uint32_t c, SATSolverProxy &S,
              std::vector<uint32_t> &lits, int varZero);
  void genCCl1(uint32_t a, uint32_t b, uint32_t c, SATSolverProxy &S,
               std::vector<uint32_t> &lits, int varZero);
  void genCCl(uint32_t a, uint32_t b, uint32_t c, uint32_t d, SATSolverProxy &S,
              std::vector<uint32_t> &lits, int varZero);
  void genCCl(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t e,
              SATSolverProxy &S, std::vector<uint32_t> &lits, int varZero);

  // Weighted MaxSAT Totalizer
  void genKWMTO(std::vector<long long> &weights,
                std::vector<uint32_t> &blockings,
                std::vector<long long> &weightsTable, int from, int to, int div,
                uint32_t zero, std::vector<uint32_t> &lower,
                std::vector<long long> &lowerW, std::vector<uint32_t> &upper,
                std::vector<long long> &upperW, SATSolverProxy &S, long long ub,
                std::vector<uint32_t> &lits, int varZero);
  void genKWMTO0(std::vector<long long> &weights,
                 std::vector<uint32_t> &blockings, long long max, long long k,
                 std::vector<long long> &divisors, SATSolverProxy &S,
                 std::vector<uint32_t> &lits,
                 std::vector<std::vector<uint32_t>> &linkingints,
                 std::vector<std::vector<long long>> &linkingWeights);

  // For all incremental QMaxSAT Encodings calls
  void lessthanMR(std::vector<std::vector<uint32_t>> &linkings,
                  std::vector<std::vector<long long>> &linkingWeights,
                  long long ok, long long k, std::vector<long long> &divisors,
                  std::vector<long long> &, SATSolverProxy &S,
                  std::vector<uint32_t> &lits, EncodingType encoding);

  void lessthan(std::vector<uint32_t> &linking,
                std::vector<long long> &linkingWeight, long long ok,
                long long k, long long divisor, std::vector<long long> &cc,
                SATSolverProxy &S, EncodingType encoding);

  // Mixed Radix Weighted Totalizer 2017 competition version
  void genMRWTO(std::vector<long long> &weights,
                std::vector<uint32_t> &blockings,
                std::vector<long long> &weightsTable, int from, int to,
                std::vector<long long> &divisors, uint32_t zero,
                std::vector<std::vector<uint32_t>> &linkingints,
                std::vector<std::vector<long long>> &linkingWeights,
                SATSolverProxy &S, long long ub, std::vector<uint32_t> &lits,
                int varZero);

  void genMRWTO0(std::vector<long long> &weights,
                 std::vector<uint32_t> &blockings, long long max, long long k,
                 std::vector<long long> &divisors, SATSolverProxy &S,
                 std::vector<uint32_t> &lits,
                 std::vector<std::vector<uint32_t>> &linkingints,
                 std::vector<std::vector<long long>> &linkingWeights,
                 EncodingType encoding);

  // Mixed Radix Weighted Totalizer 2019 competition version
  void genMRWTO19(std::vector<long long> &weights,
                  std::vector<uint32_t> &blockings,
                  std::vector<long long> &weightsTable, int from, int to,
                  std::vector<long long> &divisors, uint32_t zero,
                  std::vector<std::vector<uint32_t>> &linkingints,
                  std::vector<std::vector<long long>> &linkingWeights,
                  SATSolverProxy &S, long long ub, std::vector<uint32_t> &lits,
                  int varZero);
  void genMRWTO19_0(std::vector<long long> &weights,
                    std::vector<uint32_t> &blockings, long long max,
                    long long k, std::vector<long long> &divisors,
                    SATSolverProxy &S, std::vector<uint32_t> &lits,
                    std::vector<std::vector<uint32_t>> &linkingints,
                    std::vector<std::vector<long long>> &linkingWeights,
                    EncodingType encoding);

 private:
  long long med3(long long x, long long y, long long z);
  void quicksort(long long a[], int left, int right);
  long long sumWeight(std::vector<long long> &weights);
};
} // Namespace Pacose

#endif  // ENCODINGS_H
