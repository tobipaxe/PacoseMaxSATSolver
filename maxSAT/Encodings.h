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
  unsigned _relaxLit;

  // Warners adder encoding
  void genWarnersHalf(unsigned &a, unsigned &b, unsigned &carry, unsigned &sum,
                      int comp, SATSolverProxy &S, std::vector<unsigned> &lits);
  void genWarnersFull(unsigned &a, unsigned &b, unsigned &c, unsigned &carry,
                      unsigned &sum, int comp, SATSolverProxy &S,
                      std::vector<unsigned> &lits);
  void genWarners(std::vector<long long> &weights,
                  std::vector<unsigned> &blockings, long long max, int k,
                  int comp, SATSolverProxy &S, const unsigned zero,
                  std::vector<unsigned> &lits,
                  std::vector<unsigned> &linkingint);
  void genWarners0(std::vector<long long> &weights,
                   std::vector<unsigned> &blockings, long long max, long long k,
                   int comp, SATSolverProxy &S, std::vector<unsigned> &lits,
                   std::vector<unsigned> &linkingint);

  // Bailleux totalizer encoding
  void genBailleux(std::vector<long long> &weights,
                   std::vector<unsigned> &blockings, long long total,
                   unsigned zero, unsigned one, int comp, SATSolverProxy &S,
                   std::vector<unsigned> &lits,
                   std::vector<unsigned> &linkingint, long long UB);
  void genBailleux0(std::vector<long long> &weights,
                    std::vector<unsigned> &blockings, long long max,
                    long long k, int comp, SATSolverProxy &S,
                    std::vector<unsigned> &lits,
                    std::vector<unsigned> &linkingint);

  // [Asin et. al 2011] encoding
  void sComparator(unsigned &a, unsigned &b, unsigned &c1, unsigned &c2,
                   int comp, SATSolverProxy &S, std::vector<unsigned> &lits);
  void genSMerge(std::vector<unsigned> &linkA, std::vector<unsigned> &linkB,
                 unsigned zero, unsigned one, int comp, SATSolverProxy &S,
                 std::vector<unsigned> &lits, std::vector<unsigned> &linkingint,
                 long long UB);
  void genKCard(std::vector<long long> &weights,
                std::vector<unsigned> &blockings, long long total,
                unsigned zero, unsigned one, int comp, SATSolverProxy &S,
                std::vector<unsigned> &lits, std::vector<unsigned> &linkingint,
                long long UB);
  void genAsin(std::vector<long long> &weights,
               std::vector<unsigned> &blockings, long long max, long long k,
               int comp, SATSolverProxy &S, std::vector<unsigned> &lits,
               std::vector<unsigned> &linkingint);

  // Ogawa Modulo Totalizer
  void genOgawa(long long weightX, std::vector<unsigned> &linkingX,
                long long weightY, std::vector<unsigned> &linkingY,
                long long &total, long long divisor, unsigned zero,
                unsigned one, int comp, SATSolverProxy &S,
                std::vector<unsigned> &lits, std::vector<unsigned> &linkingint,
                long long UB);
  void genOgawa(std::vector<long long> &weights,
                std::vector<unsigned> &blockings, long long &total,
                long long divisor, unsigned zero, unsigned one, int comp,
                SATSolverProxy &S, std::vector<unsigned> &lits,
                std::vector<unsigned> &linkingint, long long UB);
  void genOgawa0(std::vector<long long> &weights,
                 std::vector<unsigned> &blockings, long long max, long long k,
                 long long &divisor, int comp, SATSolverProxy &S,
                 std::vector<unsigned> &lits,
                 std::vector<unsigned> &linkingint);

  // BailleuxW2
  void genBailleuxW2(std::vector<long long> &weights,
                     std::vector<unsigned> &blockings, long long total,
                     unsigned zero, unsigned one, int comp, SATSolverProxy &S,
                     std::vector<unsigned> &lits,
                     std::vector<unsigned> &linkingint,
                     std::vector<long long> &linkingW, long long UB);
  void genBailleuxW20(std::vector<long long> &weights,
                      std::vector<unsigned> &blockings, long long max,
                      long long k, int comp, SATSolverProxy &S,
                      std::vector<unsigned> &lits,
                      std::vector<unsigned> &linkingint,
                      std::vector<long long> &linkingWeight);
  void genCCl(unsigned a, SATSolverProxy &S, std::vector<unsigned> &lits,
              int varZero);
  void genCCl(unsigned a, unsigned b, SATSolverProxy &S,
              std::vector<unsigned> &lits, int varZero);
  void genCCl(unsigned a, unsigned b, unsigned c, SATSolverProxy &S,
              std::vector<unsigned> &lits, int varZero);
  void genCCl1(unsigned a, unsigned b, unsigned c, SATSolverProxy &S,
               std::vector<unsigned> &lits, int varZero);
  void genCCl(unsigned a, unsigned b, unsigned c, unsigned d, SATSolverProxy &S,
              std::vector<unsigned> &lits, int varZero);
  void genCCl(unsigned a, unsigned b, unsigned c, unsigned d, unsigned e,
              SATSolverProxy &S, std::vector<unsigned> &lits, int varZero);

  // Weighted MaxSAT Totalizer
  void genKWMTO(std::vector<long long> &weights,
                std::vector<unsigned> &blockings,
                std::vector<long long> &weightsTable, int from, int to, int div,
                unsigned zero, std::vector<unsigned> &lower,
                std::vector<long long> &lowerW, std::vector<unsigned> &upper,
                std::vector<long long> &upperW, SATSolverProxy &S, long long ub,
                std::vector<unsigned> &lits, int varZero);
  void genKWMTO0(std::vector<long long> &weights,
                 std::vector<unsigned> &blockings, long long max, long long k,
                 std::vector<long long> &divisors, SATSolverProxy &S,
                 std::vector<unsigned> &lits,
                 std::vector<std::vector<unsigned>> &linkingints,
                 std::vector<std::vector<long long>> &linkingWeights);

  // Mixed Radix Weighted Totalizer
  void genMRWTO(std::vector<long long> &weights,
                std::vector<unsigned> &blockings,
                std::vector<long long> &weightsTable, int from, int to,
                std::vector<long long> &divisors, unsigned zero,
                std::vector<std::vector<unsigned>> &linkingints,
                std::vector<std::vector<long long>> &linkingWeights,
                SATSolverProxy &S, long long ub, std::vector<unsigned> &lits,
                int varZero);

  void genMRWTO0(std::vector<long long> &weights,
                 std::vector<unsigned> &blockings, long long max, long long k,
                 std::vector<long long> &divisors, SATSolverProxy &S,
                 std::vector<unsigned> &lits,
                 std::vector<std::vector<unsigned>> &linkingints,
                 std::vector<std::vector<long long>> &linkingWeights,
                 EncodingType encoding);

  // For all incremental QMaxSAT Encodings calls
  void lessthanMR(std::vector<std::vector<unsigned>> &linkings,
                  std::vector<std::vector<long long>> &linkingWeights,
                  long long ok, long long k, std::vector<long long> &divisors,
                  std::vector<long long> &, SATSolverProxy &S,
                  std::vector<unsigned> &lits, EncodingType encoding);

  void lessthan(std::vector<unsigned> &linking,
                std::vector<long long> &linkingWeight, long long ok,
                long long k, long long divisor, std::vector<long long> &cc,
                SATSolverProxy &S, EncodingType encoding);

  // Mixed Radix Weighted Totalizer 2019 competition version
  void genMRWTO19(std::vector<long long> &weights,
                  std::vector<unsigned> &blockings,
                  std::vector<long long> &weightsTable, int from, int to,
                  std::vector<long long> &divisors, unsigned zero,
                  std::vector<std::vector<unsigned>> &linkingints,
                  std::vector<std::vector<long long>> &linkingWeights,
                  SATSolverProxy &S, long long ub, std::vector<unsigned> &lits,
                  int varZero);
  void genMRWTO19_0(std::vector<long long> &weights,
                    std::vector<unsigned> &blockings, long long max,
                    long long k, std::vector<long long> &divisors,
                    SATSolverProxy &S, std::vector<unsigned> &lits,
                    std::vector<std::vector<unsigned>> &linkingints,
                    std::vector<std::vector<long long>> &linkingWeights,
                    EncodingType encoding);

 private:
  long long med3(long long x, long long y, long long z);
  void quicksort(long long a[], int left, int right);
  long long sumWeight(std::vector<long long> &weights);
};
} // Namespace Pacose

#endif  // ENCODINGS_H
