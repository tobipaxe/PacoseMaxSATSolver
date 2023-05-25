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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>  // timing
#include <sys/time.h>      // timing
#include <zlib.h>
#include "Pacose.h"
#include "Softclause.h"

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

template <class B>
static void skipWhitespace(B &in) {
  while ((*in >= 9 && *in <= 13) || *in == 32) ++in;
}

template <class B>
static void skipLine(B &in) {
  for (;;) {
    if (isEof(in)) return;
    if (*in == '\n') {
      ++in;
      return;
    }
    ++in;
  }
}

template <class B>
static double parseDouble(B &in) {  // only in the form X.XXXXXe-XX
  bool neg = false;
  double accu = 0.0;
  double currentExponent = 1;
  int exponent;

  skipWhitespace(in);
  if (*in == EOF) return 0;
  if (*in == '-')
    neg = true, ++in;
  else if (*in == '+')
    ++in;
  if (*in < '1' || *in > '9')
    printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  accu = (double)(*in - '0');
  ++in;
  if (*in != '.') printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  ++in;  // skip dot
  currentExponent = 0.1;
  while (*in >= '0' && *in <= '9') {
    accu = accu + currentExponent * ((double)(*in - '0'));
    currentExponent /= 10;
    ++in;
  }
  if (*in != 'e') printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  ++in;                     // skip dot
  exponent = parseInt(in);  // read exponent
  accu *= pow(10, exponent);
  return neg ? -accu : accu;
}

template <class B>
static long long int parseInt(B &in) {
  long long int val = 0;
  bool neg = false;
  skipWhitespace(in);
  if (*in == '-')
    neg = true, ++in;
  else if (*in == '+')
    ++in;
  if (*in < '0' || *in > '9')
    fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  while (*in >= '0' && *in <= '9') val = val * 10 + (*in - '0'), ++in;
  return neg ? -val : val;
}

// String matching: in case of a match the input iterator will be advanced the
// corresponding number of characters.
template <class B>
static bool match(B &in, const char *str) {
  int i;
  for (i = 0; str[i] != '\0'; i++)
    if (in[i] != str[i]) return false;

  in += i;

  return true;
}

// String matching: consumes characters eagerly, but does not require random
// access iterator.
template <class B>
static bool eagerMatch(B &in, const char *str) {
  for (; *str != '\0'; ++str, ++in)
    if (*str != *in) return false;
  return true;
}

//=================================================================================================
// DIMACS Parser:

// koshi 20140106
template <class B>
static bool readClause(B &in, std::vector<unsigned> &lits,
                       long long int &weight, long long int top,
                       long long int &vars, Pacose *pacose) {
  long long int parsed_lit;  // koshi 20140106

  bool soft = false;
  lits.clear();
  if (top == 0) {
    // unweighted MaxSAT
    // koshi 20140106 (ms: all clauses are 1-weighted soft)
    soft = true;
    weight = 1;
  } else {
    // weighted MaxSAT
    // koshi 20140106
    parsed_lit = parseInt(in);
    if ((1 <= parsed_lit && parsed_lit < top) || top == -1) {
      // soft clause
      // top == -1 indicates all clauses are weighted
      soft = true;
      weight = parsed_lit;
    }
  }

  for (;;) {
    parsed_lit = parseInt(in);
    if (parsed_lit == 0) break;
    unsigned var = std::abs(parsed_lit);
    lits.push_back((var << 1) ^ (parsed_lit <= 0));
    // compensate a too small variable number in header!
    if (vars < var) {
      std::cout << "c ERROR: wrong number or vars in the header - can lead to "
                   "wrong result!"
                << std::endl;
      exit(0);

      //      std::cout << "vars: " << vars << "  var: " << var << std::endl;
      //      vars++;
      //      pacose->_satSolver->NewVariable();
    }
  }

  return soft;
}
// koshi 20140106
template <class B, class Pacose>
static long long int parse_DIMACS_main(B &in, Pacose &pacose,
                                       unsigned verbosity) {
  //    int out_nbvar = 0;
  long long int out_top = 0;
  int out_nbsoft = 0;

  std::vector<unsigned> lits;

  long long int vars = 0;
  long long int clauses = 0;
  long long int cnt = 0;

  // koshi 20140106
  int uscnt = 0;
  long long int sWeight = 0;   // sum of weights of soft clauses
  long long int usWeight = 0;  // sum of weights of unit soft clauses
  //  weights.clear(); blockings.clear();
  long long int top = -2;

  for (;;) {
    skipWhitespace(in);
    if (*in == EOF) {
      break;
    } else if (*in == 'c' || *in == '\n') {
      skipLine(in);
    } else if (*in == 'p') {
      top = 0;
      if (eagerMatch(in, "p cnf")) {  // koshi 20140106 (Unweighted MaxSAT)
        vars = parseInt(in);
        clauses = parseInt(in);
        // SATRACE'06 hack
        // if (clauses > 4000000)
        //     S.eliminate(true);
        top = 0;  // all clauses are 1-weighted soft
      } else if (eagerMatch(in, "wcnf")) {
        vars = parseInt(in);
        clauses = parseInt(in);
        while ((*in == 9 || *in == 32))
          ++in;  // koshi 20140117 skip space and tab
        // skipWhitespace(in);
        if (*in >= '0' && *in <= '9') {
          top = parseInt(in);
          //          printf("c top = %lld\n", top);
        } else {
          top = -1;  // weighted Max-SAT (no hard clause)
                     //          printf("c weighted Max-SAT (wms)\n");
        }
      } else {
        printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
      }
      //            out_nbvar   = vars;
      out_top = top;
      //            while (vars > pacose.nVars()) pacose.newVar();
      pacose->_satSolver->NewVariables(vars);

      if (verbosity > 0) {
        printf(
            "c |  Number of variables:    %-12lld                              "
            "  "
            "     |\n",
            vars);
        printf(
            "c |  Number of clauses:      %-12lld                              "
            "  "
            "     |\n",
            clauses);
        if (top == -1)
          printf(
              "c |  all clauses are weigthed soft                              "
              "  "
              "         "
              "    |\n");
        else if (top == 0)
          printf(
              "c |  all clauses are 1-weigthed soft                            "
              "  "
              "         "
              "    |\n");
        else
          printf(
              "c |  Weight of hard clauses: %-12lld                            "
              "  "
              "       |\n",
              top);
      }
    } else {
      // soft or hard clause
      long long int weight = 0;
      cnt++;
      // this should be always true, then out_top can be replaced by top!
      assert(out_top == top);
      if (readClause(in, lits, weight, out_top, vars,
                     pacose)) {  // soft clause
        out_nbsoft++;
        sWeight += weight;
        pacose->AddSoftClause(lits, weight);
        if (lits.size() == 1) {  // unit soft clause
          uscnt++;
          usWeight += weight;
        }
        continue;
      } else if (!pacose->_hasHardClauses) {
        // is hard clause
        pacose->_hasHardClauses = true;
      }
      if (lits.size() > 0) {
        pacose->_satSolver->ResetClause();
        pacose->_satSolver->NewClause();
        for (unsigned lit : lits) {
          if (lit == 0) {
            std::cout << "c Invalid literal 0!" << std::endl;
          }
          pacose->_satSolver->AddLiteral(&lit);
        }
        pacose->_satSolver->CommitClause();
        //        if (saveHardClauses) {
        //          hardClauses.push
        //        }
        //        pacose->_satSolver->AddClause(lits);
      } else {
        std::cout << "c ERROR: Trying to add a empty clause!!! -- skip!"
                  << std::endl;
        assert(1);
      }
    }
  }

  pacose->SetSumOfSoftWeights(sWeight);

  printf("c #Clauses...............: %-12lld\n", cnt);
  printf("c #SoftClauses...........: %-12d\n", out_nbsoft);
  printf("c Weight of SoftClauses..: %-12lld\n", sWeight);

  //  printf(
  //      "c |  Number of unit soft clauses: %-12d "
  //      "|\n",
  //      uscnt);
  // koshi 2013.05.21

  // koshi 20140107
  //    if (vars != S.nVars())
  //        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of
  //        variables.\n");

  if (cnt != clauses) {
    fprintf(stderr,
            "WARNING! DIMACS header mismatch: wrong number of clauses.\n");
  }
  return vars;
}

// Inserts problem into solver.
//
template <class Pacose>
static long long int parse_DIMACS(std::string *maxCNFFile, Pacose *pacose,
                                  unsigned verbosity) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;

  gzFile input_stream = gzopen(maxCNFFile->c_str(), "rb");
  if (input_stream == nullptr) {
    std::cout << "c ERROR! Could not open file: " << *maxCNFFile << std::endl;
    return 0;
  }

  StreamBuffer in(input_stream);

  long long int noVars = parse_DIMACS_main(in, pacose, verbosity);

  gzclose(input_stream);

  return noVars;
}

} // Namespace Pacose
//=================================================================================================

#endif  // MAXSATDIMACS_H
