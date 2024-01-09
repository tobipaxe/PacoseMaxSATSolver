/****************************************************************************************[Dimacs.h]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson
Updated by Tobias Paxian 2020

Permission is hereby granted, free of charge, to any person obtaining a copy of
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
**************************************************************************************************/
/* koshi 20140124
#ifndef Minisat_Dimacs_h
#define Minisat_Dimacs_h
*/
#ifndef MAXSATDIMACS_H
#define MAXSATDIMACS_H

#include "../maxSAT/Pacose.h"
#include "../maxSAT/Softclause.h"
#include "ClauseDB.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h> // timing
#include <sys/time.h>     // timing
#include <zlib.h>

namespace Pacose {

//-------------------------------------------------------------------------------------------------
// A simple buffered character stream class:
// static const int buffer_size = 1048576;

class StreamBuffer {
  gzFile in;
  unsigned char buf[1048576];
  int pos;
  int size;

  void assureLookahead() {
    if (pos >= size) {
      pos = 0;
      size = gzread(in, buf, sizeof(buf));
    }
  }

public:
  explicit StreamBuffer(gzFile i) : in(i), pos(0), size(0) {
    assureLookahead();
  }

  int operator*() const { return (pos >= size) ? EOF : buf[pos]; }
  void operator++() {
    pos++;
    assureLookahead();
  }
  int position() const { return pos; }
};

//-------------------------------------------------------------------------------------------------
// End-of-file detection functions for StreamBuffer and char*:

static inline bool isEof(StreamBuffer &in) { return *in == EOF; }
static inline bool isEof(const char *in) { return *in == '\0'; }

//-------------------------------------------------------------------------------------------------
// Generic parse functions parametrized over the input-stream type.

template <class B> static void skipWhitespace(B &in) {
  while ((*in >= 9 && *in <= 13) || *in == 32)
    ++in;
}

template <class B> static void skipLine(B &in) {
  for (;;) {
    if (isEof(in))
      return;
    if (*in == '\n') {
      ++in;
      return;
    }
    ++in;
  }
}

template <class B> static int parseInt(B &in) {
  int val = 0;
  bool neg = false;
  skipWhitespace(in);

  if (*in == '-')
    neg = true, ++in;
  else if (*in == '+')
    ++in;
  if (*in < '0' || *in > '9') {
    std::cerr << "PARSE ERROR!!! Unexpected char in parseInt: " << *in
              << std::endl;
    exit(3);
  }
  while (*in >= '0' && *in <= '9') {
    val = val * 10 + (*in - '0');
    ++in;
  }
  assert(val != 0);
  return neg ? -val : val;
}

template <class B> static uint64_t parseWeight(B &in) {

  uint64_t val = 0;
  skipWhitespace(in);
  if (*in == 'h') {
    ++in;
    if (*in != ' ') {
      std::cerr << "o PARSE ERROR! Unexpected char in parseWeight: " << *in
                << std::endl;
      exit(3);
    }
    return UINT64_MAX;
  }
  if (*in == '-') {
    std::cerr << "o PARSE ERROR! Unexpected negative weight" << std::endl;
    exit(3);
  } else if (*in == '+')
    ++in;
  if (*in < '0' || *in > '9') {
    std::cerr << "o PARSE ERROR! Unexpected char in parseWeight: " << *in
              << std::endl;
    exit(3);
  }
  while (*in >= '0' && *in <= '9')
    val = val * 10 + (*in - '0'), ++in;
  return val;
}

// String matching: consumes characters eagerly, but does not require random
// access iterator.
template <class B> static bool eagerMatch(B &in, const char *str) {
  for (; *str != '\0'; ++str, ++in)
    if (*str != *in)
      return false;
  return true;
}

// void add_new_clause(const std::vector<int> &literals, uint64_t &weight) {
//   size_t size = literals.size();
//   size_t bytes = sizeof(struct clause) + size * sizeof(int);
//   clause *cl = (clause *)new char[bytes];

//   cl->size = size;

//   for (uint32_t i = 0; i < literals.size(); i++)
//     cl->literals[i] = literals[i];

//   // debug(cl, 3);
//   cl->weight = weight;

//   assert(clauses.size() <= (size_t)INT_MAX);
//   if (weight == UINT64_MAX)
//     clauses.push_back(cl);
//   else
//     sclauses.push_back(cl);
// }

void parseWCNF(const std::string &wcnfFile, ClauseDB &clauseDB) {
  // dout2 << __PRETTY_FUNCTION__ << std::endl;

  gzFile input_stream = gzopen(wcnfFile.c_str(), "rb");
  if (input_stream == nullptr) {
    std::cout << "c ERROR! Could not open file: " << wcnfFile << std::endl;
    return;
  }
  StreamBuffer in(input_stream);
  // ClauseDB::clause clause = new ClauseDB::clause;
  std::vector<int> cl;
  uint64_t weight;

  while (true) {
    if (isEof(in)) {
      break;
    } else if (*in == 'c' || *in == '\n') {
      skipLine(in);
      continue;
    }

    weight = parseWeight(in);
    if (weight != UINT64_MAX) {
      if (weight > clauseDB.maxWeight)
        clauseDB.maxWeight = weight;
      clauseDB.sumOfSoftWeights += weight;
      clauseDB.nbSoftClauses++;
    } else
      clauseDB.nbHardClauses++;
    cl.clear();
    while (true) {
      skipWhitespace(in);
      if (*in == '0') {
        skipLine(in);
        break;
      }
      // cl->push_back(clauseDB.SignedTouint32_tLit(parseInt(in)));
      cl.push_back(parseInt(in));
      if (abs(cl.back()) > clauseDB.nbVars)
        clauseDB.nbVars = abs(cl.back());
    }

    if(cl.size() == 1 && weight != UINT64_MAX)
      clauseDB.nbUnitSoftClauses++;

    clauseDB.clauses.push_back(cl);
    clauseDB.weights.push_back(weight);

    // if (weight == UINT64_MAX) {

    //   // std::cout << clauseDB.clauses.back()->literals.front() << std::endl;
    // } else {
    //   clauseDB.sclauses.push_back(std::make_pair(cl,weight));
    //   // std::cout << clauseDB.sclauses.back()->literals.front() <<
    //   std::endl;
    // }
  }
  gzclose(input_stream);
  uint64_t top = clauseDB.sumOfSoftWeights + 1;
  for (size_t i = 0; i < clauseDB.weights.size(); ++i)
    if (clauseDB.weights[i] == UINT64_MAX)
      clauseDB.weights[i] = top;
  // std::cout << clauseDB.nbVars << std::endl;
  // std::cout << "HardClauses: " << clauseDB.clauses.size() << std::endl;
  // std::cout << "SoftClauses: " << clauseDB.sclauses.size() << std::endl;
}

} // Namespace Pacose
//=================================================================================================

#endif // MAXSATDIMACS_H
