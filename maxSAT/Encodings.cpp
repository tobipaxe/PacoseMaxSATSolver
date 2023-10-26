/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson
Changes are made for QMaxSAT by the QMaxSAT team.
Updated by Tobias Paxian 2020

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
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
**************************************************************************************************/
#include <errno.h> 

#include <signal.h>
#include <zlib.h>

//#include "utils/System.h"
//#include "utils/ParseUtils.h" // koshi 20170630
//#include "utils/Options.h"
//#include "MaxSATDimacs.h" // koshi 20170630
//#include "core/Solver.h"
//#include "Pacose.h"
#include <math.h>  // pow
#include "Encodings.h"
#include "../solver-proxy/SATSolverProxy.h"
#include "assert.h"


namespace Pacose {
/*
  koshi 20140106
  based on minisat2-070721/maxsat0.2e
 */
// koshi 20140106
#define TOTALIZER 128

//=================================================================================================

/*
  Cardinality Constraints:
  Joost P. Warners, "A linear-time transformation of linear inequalities
  into conjunctive normal form",
  Information Processing Letters 68 (1998) 63-69
 */
// koshi 2013.04.16
// koshi 13.04.05, 13.06.28, 13.07.01, 13.10.04
void Encodings::lessthan(std::vector<uint32_t> &linking,
                         std::vector<long long int> &linkingWeight,
                         long long int ok, long long int k,
                         long long int divisor,  // koshi 13.10.04
                         std::vector<long long int> &cc, SATSolverProxy &S,
                         EncodingType encoding) {
  assert(k > 0);
  _relaxLit = S.NewVariable() << 1;
  //  std::cout << "new relaxlit: " << _relaxLit << std::endl;
  //  if (_settings->verbosity > 3) {
  //    std::cout << "LKs: " << linking.size() << ": ";
  //    for (auto link : linking) {
  //      std::cout << link << ", ";
  //    }
  //    std::cout << std::endl;
  //    std::cout << "LWs: " << linkingWeight.size() << ": ";
  //    for (auto linkw : linkingWeight) {
  //      std::cout << linkw << ", ";
  //    }
  //    std::cout << std::endl;
  //    std::cout << "OK : " << ok << std::endl;
  //    std::cout << "k  : " << k << std::endl;
  //    std::cout << "DIV: " << divisor << std::endl;
  //    std::cout << "ccs: " << cc.size() << ": ";
  //    for (auto cpart : cc) {
  //      std::cout << cpart << ", ";
  //    }
  //    std::cout << std::endl;
  //  }

  if (linking.size() == 0) {
  } else                           // koshi 20140124 20140129
      if (encoding == BAILLEUX) {  // Bailleux encoding (Totalizer)
    for (long long int i = k; i < linking.size() && i < ok; i++) {
      S.ResetClause();
      S.NewClause();
      S.AddLiteral(linking[i] ^ 1);
      S.AddLiteral(_relaxLit);
      S.CommitClause();
    }
  } else if (encoding == ASIN) {  // Asin encoding
    for (long long int i = k - 1; i < linking.size() && i < ok; i++) {
      S.ResetClause();
      S.NewClause();
      S.AddLiteral(linking[i] ^ 1);
      S.AddLiteral(_relaxLit);
      S.CommitClause();
    }

  } else if (encoding == BAILLEUXW2) {  // Weighted  Bailleux encoding (Weighted
                                        // Totalizer) hayata 2014/12/17

    for (long long int i = 1; i < linking.size(); i++) {
      long long int tmp_w = linkingWeight[i];
      if (tmp_w >= k && tmp_w < ok) {
        S.ResetClause();
        S.NewClause();
        S.AddLiteral(linking[i] ^ 1);
        S.AddLiteral(_relaxLit);
        S.CommitClause();
        // if(S.verbosity > 1)
        //	printf("[%lld]\n" , tmp_w);
      }
      // printf("[-%d]" , var(lits[0])+1);
    }
  } else if (encoding == OGAWA) {  // Ogawa encoding (Modulo Totalizer)
    long long int upper = (k - 1) / divisor;
    long long int lower = k % divisor;
    long long int oupper = ok / divisor;
    //    printf("upper = %lld, oupper = %lld\n", upper,oupper);
    if (upper < oupper)
      for (long long int i = divisor + upper + 1; i < divisor + oupper + 1;
           i++) {
        if (linking.size() <= i)
          break;
        else {
          //  printf("linking i = %lld ",i);
          S.ResetClause();
          S.NewClause();
          S.AddLiteral(linking[i] ^ 1);
          S.AddLiteral(_relaxLit);
          S.CommitClause();
        }
      }
    upper = k / divisor;
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(linking[divisor + upper] ^ 1);
    S.AddLiteral(linking[lower] ^ 1);
    S.AddLiteral(_relaxLit);
    //    printf("divisor+upper = %lld, lower = %lld\n",divisor+upper,lower);
    S.CommitClause();
  } else if (encoding == WARNERS) {  // Warners encoding
    std::vector<long long int> cls;
    cls.clear();

    k--;
    if (k % 2 == 0) cls.push_back(1);
    k = k / 2;
    int cnt = 1;
    long long int pos = 0x0002LL;
    while (k > 0) {
      if (k % 2 == 0) cls.push_back(pos);
      //    else if (cls.size() == 0) cls.push_back(pos);
      else
        for (int i = 0; i < cls.size(); i++) cls[i] = cls[i] | pos;
      pos = pos << 1;
      k = k / 2;
      cnt++;
    }
    for (int i = cnt; i < linking.size(); i++) {
      cls.push_back(pos);
      pos = pos << 1;
    }
    for (int i = 0; i < cls.size(); i++) {
      long long int x = cls[i];
      bool found = false;
      for (int j = 0; j < cc.size(); j++) {
        if (x == cc[j]) {
          found = true;
          break;
        }
      }
      if (!found) {
        cc.push_back(x);  // koshi 2013.10.04
        S.ResetClause();
        S.NewClause();
        int j = 0;
        while (x > 0) {
          if ((x & 0x0001L) == 0x0001L) {
            S.AddLiteral(linking[j] ^ 1);
            //            std::cout << (linking[j] ^ 1) << ", ";
          }
          x = x >> 1;
          j++;
        }
        S.AddLiteral(_relaxLit);
        //        std::cout << "+ " << _relaxLit << std::endl << std::endl;
        S.CommitClause();
      }
    }
  }

  _relaxLit = _relaxLit ^ 1;
  S.AddAssumption(&_relaxLit);
  _relaxLit = _relaxLit ^ 1;
  //  std::cout << "assumption: " << _relaxLit << std::endl;
}

// uemura 20161128
void Encodings::lessthanMR(
    std::vector<std::vector<uint32_t>> &linkings,
    std::vector<std::vector<long long int>> &linkingWeights, long long int ok,
    long long int k, std::vector<long long int> &divisors,
    std::vector<long long int> & /* cc */, SATSolverProxy &S,
    std::vector<uint32_t> &lits, EncodingType encoding) {
  assert(k > 0);
  _relaxLit = S.NewVariable() << 1;

  if (linkings[0].size() == 0) {
  } else {  // uemura 20161112
    int ndigit = linkings.size();
    long long int tmp_k = k;
    long long int tmp_ok = ok;
    // okは前回のk

    std::vector<uint32_t> control;

    //各桁を計算して表示する
    int *sp_k = new int[ndigit];
    int *sp_ok = new int[ndigit];
    for (int i = 0; i < ndigit - 1; i++) {
      sp_k[i] = tmp_k % divisors[i];
      sp_ok[i] = tmp_ok % divisors[i];
      tmp_k = tmp_k / divisors[i];
      tmp_ok = tmp_ok / divisors[i];
    }
    sp_k[ndigit - 1] = tmp_k;
    sp_ok[ndigit - 1] = tmp_ok;

    if (_settings->verbosity > 1) {
      printf("c k = %lld(", k);
      for (int i = ndigit - 1; i >= 0; i--) {
        printf("%d/", sp_k[i]);
      }
      printf("\b)");
      for (int i = ndigit - 2; i >= 0; i--) {
        printf("p%d = %lld ", i, divisors[i]);
      }
      printf("\n");
    }
    /* 自分より下の桁が全て0である場合、自分の否定も問題に加える
     * そうでない場合は、自分より大きなものの否定を問題に加える。
     */
    for (int i = ndigit - 1; i >= 0; i--) {
      int cnr_k = 0;
      int cnr_ok = 0;
      for (int k = 0; k < i; k++) {
        if (sp_k[k] > 0) {
          cnr_k = 1;
        }
        if (sp_ok[k] > 0) {
          cnr_ok = 1;
        }
      }
      if (cnr_k == 1) {
        sp_k[i]++;
      }
      if (cnr_ok == 1) {
        sp_ok[i]++;
      }
    }

    for (int cdigit = ndigit - 1; cdigit >= 0; cdigit--) {
      int checknextdigit = 0;
      int sp_k2 = -1;
      long long int tmp_max = 0;
      //現在の桁より下の桁がすべて0出ないことのチェック
      for (int i = cdigit - 1; i >= 0; i--) {
        if (sp_k[i] > 0) {
          checknextdigit = 1;
        }
      }
      //最上位桁の処理=======================================================================================================
      if (cdigit == ndigit - 1) {
        if (linkings[cdigit].size() == 0) {
          fprintf(stderr, "ERROR0 : link size digit[%d] = 0\t@less thanMR\n",
                  cdigit);
          exit(1);
        } else {  // uemura 20161112

          for (long long int i = 1; i < linkings[cdigit].size(); i++) {
            if (linkingWeights[cdigit][i] >= sp_k[cdigit]) {
              if (linkingWeights[cdigit][i] < sp_ok[cdigit]) {  // tmp_ok=>sp_ok
                S.ResetClause();
                S.NewClause();
                S.AddLiteral(linkings[cdigit][i] ^ 1);
                S.AddLiteral(_relaxLit);
                S.CommitClause();
              }
            } else if (checknextdigit == 1) {
              if (tmp_max < linkingWeights[cdigit][i]) {
                tmp_max = linkingWeights[cdigit][i];
                sp_k2 = i;
              }
            }
          }
          if (sp_k[cdigit] > 1 && sp_k2 > 0) {
            control.push_back(linkings[cdigit][sp_k2] ^ 1);
          }
        }
      }

      else if (cdigit >= 0) {
        //最上位より下の桁の処理================================================================================================-
        if (linkings[cdigit].size() == 0) {
          fprintf(stderr, "ERROR1 : link size digit[%d] = 0\t@less thanMR\n",
                  cdigit);
          exit(1);
        } else if (sp_k[cdigit] > 0) {  // uemura 20161112
          for (long long int i = 1; i < linkings[cdigit].size(); i++) {
            if (linkingWeights[cdigit][i] >= sp_k[cdigit]) {
              S.ResetClause();
              S.NewClause();
              for (int ctr = 0; ctr < control.size(); ++ctr) {
                S.AddLiteral(control[ctr]);
              }
              S.AddLiteral(linkings[cdigit][i] ^ 1);
              S.AddLiteral(_relaxLit);
              S.CommitClause();
            } else if (checknextdigit == 1) {
              if (tmp_max < linkingWeights[cdigit][i]) {
                tmp_max = linkingWeights[cdigit][i];
                sp_k2 = i;
              }
            }
          }
          if (sp_k2 > 0) {
            control.push_back(linkings[cdigit][sp_k2] ^ 1);
          } else {
            control.push_back(linkings[cdigit][0] ^ 1);
          }
        }
      }
    }
    control.clear();
    delete[] sp_k;
    delete[] sp_ok;
  }
  _relaxLit = _relaxLit ^ 1;
  S.AddAssumption(&_relaxLit);
  _relaxLit = _relaxLit ^ 1;
}

void Encodings::genWarnersHalf(uint32_t &a, uint32_t &b, uint32_t &carry,
                               uint32_t &sum, int comp, SATSolverProxy &S,
                               std::vector<uint32_t> &lits) {
  //  std::cout << "GWH " << comp << std::endl;
  // carry
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(b ^ 1);
  S.AddLiteral(carry);
  S.CommitClause();

  //  std::cout << (a ^ 1) << ", " << (b ^ 1) << ", " << carry << std::endl;
  // sum
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a);
  S.AddLiteral(b ^ 1);
  S.AddLiteral(sum);
  S.CommitClause();
  //  std::cout << (a) << ", " << (b ^ 1) << ", " << sum << std::endl;

  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(b);
  S.AddLiteral(sum);
  S.CommitClause();
  //  std::cout << (a ^ 1) << ", " << (b) << ", " << sum << std::endl;
  //
  if (comp == 1 || comp == 2 || comp == 21 || comp == 30 || comp == 99) {
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry);
    S.AddLiteral(sum);
    S.AddLiteral(a ^ 1);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry);
    S.AddLiteral(sum);
    S.AddLiteral(b ^ 1);
    S.CommitClause();
    //    std::cout << carry << ", " << sum << ", " << (a ^ 1) << std::endl;
    //    std::cout << carry << ", " << sum << ", " << (b ^ 1) << std::endl;
  }
  if (comp == 2 || comp == 20 || comp == 21 || comp == 30 || comp == 99) {
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry ^ 1);
    S.AddLiteral(sum ^ 1);
    S.CommitClause();
    S.ResetClause();

    if (comp != 21 || comp != 20) {
      S.ResetClause();
      S.NewClause();
      S.AddLiteral(carry ^ 1);
      S.AddLiteral(sum);
      S.AddLiteral(a);
      S.CommitClause();

      S.ResetClause();
      S.NewClause();
      S.AddLiteral(carry ^ 1);
      S.AddLiteral(sum);
      S.AddLiteral(b);
      S.CommitClause();
    }
    //    std::cout << (carry ^ 1) << ", " << (sum ^ 1) << std::endl;
    //    std::cout << (carry ^ 1) << ", " << (sum) << a << std::endl;
    //    std::cout << (carry ^ 1) << ", " << (sum) << b << std::endl;
  }
  // koshi 2013.05.31
  if (comp == 10 || comp == 11 || comp == 20 || comp == 21 || comp == 30 ||
      comp == 99) {  // [Warners 1996]

    // carry
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(a);
    S.AddLiteral(carry ^ 1);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(b);
    S.AddLiteral(carry ^ 1);
    S.CommitClause();

    if (comp != 21 || comp != 20) {
      // sum
      S.ResetClause();
      S.NewClause();
      S.AddLiteral(a ^ 1);
      S.AddLiteral(b ^ 1);
      S.AddLiteral(sum ^ 1);
      S.CommitClause();

      S.ResetClause();
      S.NewClause();
      S.AddLiteral(a);
      S.AddLiteral(b);
      S.AddLiteral(sum ^ 1);
      S.CommitClause();
    }
  }
}

// koshi 2013.04.16
void Encodings::genWarnersFull(uint32_t &a, uint32_t &b, uint32_t &c,
                               uint32_t &carry, uint32_t &sum, int comp,
                               SATSolverProxy &S, std::vector<uint32_t> &lits) {
  //  std::cout << "GWF " << comp << std::endl;
  // carry
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(b ^ 1);
  S.AddLiteral(carry);

  S.CommitClause();
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(c ^ 1);
  S.AddLiteral(carry);
  S.CommitClause();

  S.ResetClause();
  S.NewClause();
  S.AddLiteral(b ^ 1);
  S.AddLiteral(c ^ 1);
  S.AddLiteral(carry);
  S.CommitClause();
  //  std::cout << (a ^ 1) << ", " << (b ^ 1) << ", " << carry << std::endl;
  //  std::cout << (a ^ 1) << ", " << (c ^ 1) << ", " << carry << std::endl;
  //  std::cout << (b ^ 1) << ", " << (c ^ 1) << ", " << carry << std::endl;

  // sum
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a);
  S.AddLiteral(b);
  S.AddLiteral(c ^ 1);
  S.AddLiteral(sum);
  S.CommitClause();

  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a);
  S.AddLiteral(b ^ 1);
  S.AddLiteral(c);
  S.AddLiteral(sum);
  S.CommitClause();

  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(b);
  S.AddLiteral(c);
  S.AddLiteral(sum);
  S.CommitClause();

  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(b ^ 1);
  S.AddLiteral(c ^ 1);
  S.AddLiteral(sum);
  S.CommitClause();
  //  std::cout << (a) << ", " << (b) << ", " << (c ^ 1) << ", " << sum
  //            << std::endl;
  //  std::cout << (a) << ", " << (b ^ 1) << ", " << (c) << ", " << sum
  //            << std::endl;
  //  std::cout << (a ^ 1) << ", " << (b) << ", " << (c) << ", " << sum
  //            << std::endl;
  //  std::cout << (a ^ 1) << ", " << (b ^ 1) << ", " << (c ^ 1) << ", " << sum
  //            << std::endl;

  if (comp == 1 || comp == 2 || comp == 21 || comp == 30 || comp == 99) {
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry);
    S.AddLiteral(sum);
    S.AddLiteral(a ^ 1);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry);
    S.AddLiteral(sum);
    S.AddLiteral(b ^ 1);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry);
    S.AddLiteral(sum);
    S.AddLiteral(c ^ 1);
    S.CommitClause();
  }

  if (comp == 2 || comp == 30 || comp == 99) {
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry ^ 1);
    S.AddLiteral(sum ^ 1);
    S.AddLiteral(a);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry ^ 1);
    S.AddLiteral(sum ^ 1);
    S.AddLiteral(b);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(carry ^ 1);
    S.AddLiteral(sum ^ 1);
    S.AddLiteral(c);
    S.CommitClause();
  }
  // koshi 2013.05.31
  if (comp == 10 || comp == 11 || comp == 30 || comp == 99) {  // [Warners 1996]
    // carry
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(a);
    S.AddLiteral(b);
    S.AddLiteral(carry ^ 1);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(a);
    S.AddLiteral(c);
    S.AddLiteral(carry ^ 1);
    S.CommitClause();

    S.ResetClause();
    S.NewClause();
    S.AddLiteral(b);
    S.AddLiteral(c);
    S.AddLiteral(carry ^ 1);
    S.CommitClause();

    if (comp != 31) {
      // sum
      S.ResetClause();
      S.NewClause();
      S.AddLiteral(a);
      S.AddLiteral(b);
      S.AddLiteral(c);
      S.AddLiteral(sum ^ 1);
      S.CommitClause();

      S.ResetClause();
      S.NewClause();
      S.AddLiteral(a ^ 1);
      S.AddLiteral(b ^ 1);
      S.AddLiteral(c);
      S.AddLiteral(sum ^ 1);
      S.CommitClause();

      S.ResetClause();
      S.NewClause();
      S.AddLiteral(a ^ 1);
      S.AddLiteral(b);
      S.AddLiteral(c ^ 1);
      S.AddLiteral(sum ^ 1);
      S.CommitClause();

      S.ResetClause();
      S.NewClause();
      S.AddLiteral(a);
      S.AddLiteral(b ^ 1);
      S.AddLiteral(c ^ 1);
      S.AddLiteral(sum ^ 1);
      S.CommitClause();
    }
  }
}

#define wbsplit(half, wL, wR, ws, bs, wsL, bsL, wsR, bsR) \
  wsL.clear();                                            \
  bsL.clear();                                            \
  wsR.clear();                                            \
  bsR.clear();                                            \
  int ii = 0;                                             \
  int wsSizeHalf = ws.size() / 2;                         \
  for (; ii < wsSizeHalf; ii++) {                         \
    wsL.push_back(ws[ii]);                                \
    bsL.push_back(bs[ii]);                                \
    wL += ws[ii];                                         \
  }                                                       \
  for (; ii < ws.size(); ii++) {                          \
    wsR.push_back(ws[ii]);                                \
    bsR.push_back(bs[ii]);                                \
    wR += ws[ii];                                         \
  }

// koshi 2013.03.25
// Parallel counter
// koshi 2013.04.16, 2013.05.23
void Encodings::genWarners(std::vector<long long int> &weights,
                           std::vector<uint32_t> &blockings, long long int max,
                           int k, int comp, SATSolverProxy &S,
                           const uint32_t zero, std::vector<uint32_t> &lits,
                           std::vector<uint32_t> &linkingVar) {
  //  std::cout << weights.size() << ", " << blockings.size() << ", " << max <<
  //  ", "
  //            << k << ", " << zero << ", " << lits.size() << ", "
  //            << linkingVar.size() << std::endl;
  linkingVar.clear();
  bool dvar = (comp == 11) ? false : true;

  if (weights.size() == 1) {
    long long int weight = weights[0];
    std::vector<bool> pn;
    pn.clear();
    while (weight > 0) {
      if (weight % 2 == 0)
        pn.push_back(false);
      else
        pn.push_back(true);
      weight /= 2;
    }
    for (int i = 0; i < pn.size(); i++) {
      if (pn[i])
        linkingVar.push_back(blockings[0]);
      else
        linkingVar.push_back(zero);
    }
    pn.clear();
  } else if (weights.size() > 1) {
    long long int weightL = 0;
    long long int weightR = 0;
    std::vector<long long int> weightsL, weightsR;
    std::vector<uint32_t> blockingsL, blockingsR;

    const long long int half = max / 2;
    (void)half;  // suppress warning about unused variable
                 //    printf("Max: %lld  Half: %lld\n", max, half);
    wbsplit(half, weightL, weightR, weights, blockings, weightsL, blockingsL,
            weightsR, blockingsR);

    std::vector<uint32_t> alpha;
    std::vector<uint32_t> beta;
    uint32_t sum = S.NewVariable() << 1 /*(true,dvar)*/;
    uint32_t carry = S.NewVariable() << 1 /*(true,dvar)*/;
    genWarners(weightsL, blockingsL, weightL, k, comp, S, zero, lits, alpha);
    genWarners(weightsR, blockingsR, weightR, k, comp, S, zero, lits, beta);
    weightsL.clear();
    weightsR.clear();
    blockingsL.clear();
    blockingsR.clear();

    bool lessthan = (alpha.size() < beta.size());
    std::vector<uint32_t> &smalls = lessthan ? alpha : beta;
    std::vector<uint32_t> &larges = lessthan ? beta : alpha;
    assert(smalls.size() <= larges.size());

    genWarnersHalf(smalls[0], larges[0], carry, sum, comp, S, lits);
    linkingVar.push_back(sum);

    int i = 1;
    uint32_t carryN;
    for (; i < smalls.size(); i++) {
      sum = S.NewVariable() << 1 /*(true,dvar)*/;
      carryN = S.NewVariable() << 1 /*(true,dvar)*/;
      genWarnersFull(smalls[i], larges[i], carry, carryN, sum, comp, S, lits);
      linkingVar.push_back(sum);
      carry = carryN;
    }
    for (; i < larges.size(); i++) {
      sum = S.NewVariable() << 1 /*(true,dvar)*/;
      carryN = S.NewVariable() << 1 /*(true,dvar)*/;
      genWarnersHalf(larges[i], carry, carryN, sum, comp, S, lits);
      linkingVar.push_back(sum);
      carry = carryN;
    }
    linkingVar.push_back(carry);
    alpha.clear();
    beta.clear();
  }
  int lsize = linkingVar.size();
  for (int i = k; i < lsize; i++) {  // koshi 2013.05.27
    //    printf("resize: k = %d, lsize = %d\n",k,lsize);
    S.ResetClause();
    S.NewClause();
    lits.clear();
    lits.push_back(linkingVar[i] ^ 1);
    S.AddLiteral(lits[0]);
    S.CommitClause();
    //    std::cout << (linkingVar[i] ^ 1) << std::endl;
  }
  for (int i = k; i < lsize; i++) linkingVar.pop_back();  // koshi 2013.05.27
  //  std::cout << "LV: " << linkingVar.size() << std::endl;
}

// koshi 2013.06.28
void Encodings::genWarners0(std::vector<long long int> &weights,
                            std::vector<uint32_t> &blockings, long long int max,
                            long long int k, int comp, SATSolverProxy &S,
                            std::vector<uint32_t> &lits,
                            std::vector<uint32_t> &linkingVar) {
  // koshi 20140109
  //  printf("c Warners' encoding for Cardinality Constraints\n");

  //  printf("c Warners' encoding for Cardinality Constraints k=%lld\n", k);
  int logk = 1;
  while ((k >>= 1) > 0) logk++;

  //  std::cout << "logk: " << logk << std::endl;

  //  printf("c Warners' encoding for Cardinality Constraints logk=%d\n", logk);
  uint32_t zero = S.NewVariable() << 1;
  S.ResetClause();
  S.NewClause();
  lits.push_back(zero ^ 1);
  S.AddLiteral(lits[0]);
  S.CommitClause();

  genWarners(weights, blockings, max, logk, comp, S, zero, lits, linkingVar);
}

/*
  Cardinaltiy Constraints:
  Olivier Bailleux and Yacine Boufkhad,
  "Efficient CNF Encoding of Boolean Cardinality Constraints",
  CP 2003, LNCS 2833, pp.108-122, 2003
 */
// koshi 10.01.08
// 10.01.15 argument UB is added
void Encodings::genBailleux(std::vector<long long int> &weights,
                            std::vector<uint32_t> &blockings,
                            long long int total, uint32_t zero, uint32_t one,
                            int comp, SATSolverProxy &S,
                            std::vector<uint32_t> &lits,
                            std::vector<uint32_t> &linkingVar,
                            long long int UB) {
  assert(weights.size() == blockings.size());

  linkingVar.clear();
  //  bool dvar = (comp == 11) ? false : true;

  std::vector<uint32_t> linkingAlpha;
  std::vector<uint32_t> linkingBeta;

  if (blockings.size() == 1) {  // koshi 20140121
    long long int weight = weights[0];
    assert(weight < UB);
    linkingVar.push_back(one);
    for (int i = 0; i < weight; i++) linkingVar.push_back(blockings[0]);
    linkingVar.push_back(zero);
  } else if (blockings.size() > 1) {
    long long int weightL = 0;
    long long int weightR = 0;
    std::vector<long long int> weightsL, weightsR;
    std::vector<uint32_t> blockingsL, blockingsR;
    //    const long long int half = total / 2;
    //    (void)half;  // suppress warning about unused variable
    wbsplit(1, weightL, weightR, weights, blockings, weightsL, blockingsL,
            weightsR, blockingsR);
    //    S.NewVariable();

    genBailleux(weightsL, blockingsL, weightL, zero, one, comp, S, lits,
                linkingAlpha, UB);
    genBailleux(weightsR, blockingsR, weightR, zero, one, comp, S, lits,
                linkingBeta, UB);

    weightsL.clear();
    blockingsL.clear();
    weightsR.clear();
    blockingsR.clear();

    linkingVar.clear();
    linkingVar.push_back(one);
    for (int i = 0; i < total && i <= UB; i++)
      linkingVar.push_back(S.NewVariable() << 1 /*(true,dvar)*/);
    linkingVar.push_back(zero);
    if (_settings->verbosity > 1)
      std::cout << "LV.size: " << linkingVar.size() << "  total: " << total
                << " UB: " << UB << std::endl;

    for (long long int sigma = 0; sigma <= total && sigma <= UB; sigma++) {
      for (long long int alpha = 0;
           alpha < linkingAlpha.size() - 1 && alpha <= UB; alpha++) {
        long long int beta = sigma - alpha;
        if (0 <= beta && beta < linkingBeta.size() - 1 && beta <= UB) {
          S.ResetClause();
          S.NewClause();
          S.AddLiteral(linkingAlpha[alpha] ^ 1);
          S.AddLiteral(linkingBeta[beta] ^ 1);
          S.AddLiteral(linkingVar[sigma]);
          S.CommitClause();
          if (comp >= 10) {
            S.ResetClause();
            S.NewClause();
            S.AddLiteral(linkingAlpha[alpha + 1]);
            S.AddLiteral(linkingBeta[beta + 1]);
            S.AddLiteral(linkingVar[sigma + 1] ^ 1);
            S.CommitClause();
          }
        }
      }
    }
  }
  linkingAlpha.clear();
  linkingBeta.clear();
}

void Encodings::genBailleux0(std::vector<long long int> &weights,
                             std::vector<uint32_t> &blockings,
                             long long int max, long long int k, int comp,
                             SATSolverProxy &S, std::vector<uint32_t> &lits,
                             std::vector<uint32_t> &linkingVar) {
  // koshi 20140109
  if (_settings->verbosity > 0)
    printf("c Bailleux's encoding for Cardinailty Constraints k = %lld\n", k);

  uint32_t one = S.NewVariable() << 1;
  uint32_t two = one ^ 1;
  S.ResetClause();
  S.NewClause();
  lits.push_back(one);
  S.AddLiteral(lits[0]);
  S.CommitClause();

  genBailleux(weights, blockings, max, two, one, comp, S, lits, linkingVar, k);
}

/*
  Cardinaltiy Constraints:
  Robert Asin, Robert Nieuwenhuis, Albert Oliveras, Enric Rodriguez-Carbonell
  "Cardinality Networks: a theoretical and empirical study",
  Constraints (2011) 16:195-221
 */
// koshi 2013.07.01
inline void Encodings::sComparator(uint32_t &a, uint32_t &b, uint32_t &c1,
                                   uint32_t &c2, int comp, SATSolverProxy &S,
                                   std::vector<uint32_t> &lits) {
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(b ^ 1);
  S.AddLiteral(c2);
  S.CommitClause();
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(a ^ 1);
  S.AddLiteral(c1);
  S.CommitClause();
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(b ^ 1);
  S.AddLiteral(c1);
  S.CommitClause();
  if (comp >= 10) {
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(a);
    S.AddLiteral(b);
    S.AddLiteral(c1 ^ 1);
    S.CommitClause();
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(a);
    S.AddLiteral(c2 ^ 1);
    S.CommitClause();
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(b);
    S.AddLiteral(c2 ^ 1);
    S.CommitClause();
  }
}

// koshi 2013.07.01
void Encodings::genSMerge(std::vector<uint32_t> &linkA,
                          std::vector<uint32_t> &linkB, uint32_t zero,
                          uint32_t one, int comp, SATSolverProxy &S,
                          std::vector<uint32_t> &lits,
                          std::vector<uint32_t> &linkingVar, long long int UB) {
  /* koshi 2013.12.10
assert(UB > 0); is violated when k <= 1
*/

  bool lessthan = (linkA.size() <= linkB.size());
  std::vector<uint32_t> &tan = lessthan ? linkA : linkB;
  std::vector<uint32_t> &tyou = lessthan ? linkB : linkA;
  assert(tan.size() <= tyou.size());

  linkingVar.clear();
  bool dvar = (comp == 11) ? false : true;

  if (tan.size() == 0)
    for (long long int i = 0; i < tyou.size(); i++)
      linkingVar.push_back(tyou[i]);
  else if (tan.size() == 1 && tyou.size() == 1) {
    uint32_t c1 = S.NewVariable() << 1 /*(true,dvar)*/;
    uint32_t c2 = S.NewVariable() << 1 /*(true,dvar)*/;
    linkingVar.push_back(c1);
    linkingVar.push_back(c2);
    sComparator(tan[0], tyou[0], c1, c2, comp, S, lits);
  } else {
    std::vector<uint32_t> oddA, oddB, evenA, evenB;
    oddA.clear();
    oddB.clear();
    evenA.clear();
    evenB.clear();

    long long int i;
    for (i = 0; i < tan.size(); i++) {
      if (i % 2 == 0) {
        evenA.push_back(tan[i]);
        evenB.push_back(tyou[i]);
      } else {
        oddA.push_back(tan[i]);
        oddB.push_back(tyou[i]);
      }
    }
    for (; i < tyou.size(); i++) {
      if (i % 2 == 0) {
        evenA.push_back(zero);
        evenB.push_back(tyou[i]);
      } else {
        oddA.push_back(zero);
        oddB.push_back(tyou[i]);
      }
    }

    // koshi 2013.07.04
    long long int UBceil = UB / 2 + UB % 2;
    long long int UBfloor = UB / 2;
    assert(UBfloor <= UBceil);
    std::vector<uint32_t> d, e;
    genSMerge(evenA, evenB, zero, one, comp, S, lits, d, UBceil);
    genSMerge(oddA, oddB, zero, one, comp, S, lits, e, UBfloor);
    oddA.clear();
    oddB.clear();
    evenA.clear();
    evenB.clear();

    linkingVar.push_back(d[0]);

    assert(d.size() >= e.size());

    while (d.size() > e.size()) e.push_back(zero);
    for (i = 0; i < e.size() - 1; i++) {
      uint32_t c2i = S.NewVariable() << 1 /*(true,dvar)*/;
      uint32_t c2ip1 = S.NewVariable() << 1 /*(true,dvar)*/;
      linkingVar.push_back(c2i);
      linkingVar.push_back(c2ip1);
      sComparator(d[i + 1], e[i], c2i, c2ip1, comp, S, lits);
    }

    linkingVar.push_back(e[i]);

    for (long long int i = UB + 1; i < linkingVar.size(); i++) {
      S.ResetClause();
      S.NewClause();
      S.AddLiteral(linkingVar[i] ^ 1);
      S.CommitClause();
    }
    long long int ssize = linkingVar.size() - UB - 1;
    for (int j = 0; j < ssize; j++) {
      linkingVar.pop_back();
    }

    d.clear();
    e.clear();
  }
}

// koshi 2013.07.01
void Encodings::genKCard(std::vector<long long int> &weights,
                         std::vector<uint32_t> &blockings, long long int total,
                         uint32_t zero, uint32_t one, int comp,
                         SATSolverProxy &S, std::vector<uint32_t> &lits,
                         std::vector<uint32_t> &linkingVar, long long int UB) {
  linkingVar.clear();

  if (blockings.size() == 1) {
    long long int weight = weights[0];
    assert(weight <= UB);
    // koshi 20140121
    for (int i = 0; i < weight; i++) linkingVar.push_back(blockings[0]);
  } else if (blockings.size() > 1) {
    std::vector<uint32_t> linkingAlpha;
    std::vector<uint32_t> linkingBeta;

    long long int weightL = 0;
    long long int weightR = 0;
    std::vector<long long int> weightsL, weightsR;
    std::vector<uint32_t> blockingsL, blockingsR;
    const long long int half = total / 2;
    (void)half;  // suppress warning about unused variable
    wbsplit(half, weightL, weightR, weights, blockings, weightsL, blockingsL,
            weightsR, blockingsR);

    genKCard(weightsL, blockingsL, weightL, zero, one, comp, S, lits,
             linkingAlpha, UB);
    genKCard(weightsR, blockingsR, weightR, zero, one, comp, S, lits,
             linkingBeta, UB);

    genSMerge(linkingAlpha, linkingBeta, zero, one, comp, S, lits, linkingVar,
              UB);

    linkingAlpha.clear();
    linkingBeta.clear();
  }
}

// koshi 2013.07.01
void Encodings::genAsin(std::vector<long long int> &weights,
                        std::vector<uint32_t> &blockings, long long int max,
                        long long int k, int comp, SATSolverProxy &S,
                        std::vector<uint32_t> &lits,
                        std::vector<uint32_t> &linkingVar) {
  // koshi 20140109
  if (_settings->verbosity > 0)
    printf("c Asin's encoding for Cardinailty Constraints\n");

  uint32_t one = S.NewVariable() << 1;
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(one);
  S.CommitClause();

  genKCard(weights, blockings, max, one ^ 1, one, comp, S, lits, linkingVar, k);
}

/*
  Cardinaltiy Constraints:
  Toru Ogawa, YangYang Liu, Ryuzo Hasegawa, Miyuki Koshimura, Hiroshi Fujita,
  "Modulo Based CNF Encoding of Cardinality Constraints and Its Application to
   MaxSAT Solvers",
  ICTAI 2013.
 */
// koshi 2013.10.03
void Encodings::genOgawa(long long int weightX, std::vector<uint32_t> &linkingX,
                         long long int weightY, std::vector<uint32_t> &linkingY,
                         long long int &total, long long int divisor,
                         uint32_t /* zero */, uint32_t one, int /* comp */,
                         SATSolverProxy &S, std::vector<uint32_t> &lits,
                         std::vector<uint32_t> &linkingVar,
                         long long int /* UB */) {
  total = weightX + weightY;
  if (weightX == 0)
    for (int i = 0; i < linkingY.size(); i++) linkingVar.push_back(linkingY[i]);
  else if (weightY == 0)
    for (int i = 0; i < linkingX.size(); i++) linkingVar.push_back(linkingX[i]);
  else {
    long long int upper = total / divisor;
    long long int divisor1 = divisor - 1;
    /*
printf("weightX = %lld, linkingX.size() = %d ", weightX,linkingX.size());
printf("weightY = %lld, linkingY.size() = %d\n", weightY,linkingY.size());
printf("upper = %lld, divisor1 = %lld\n", upper,divisor1);
*/

    linkingVar.push_back(one);
    for (int i = 0; i < divisor1; i++)
      linkingVar.push_back(S.NewVariable() << 1);
    linkingVar.push_back(one);
    for (int i = 0; i < upper; i++) linkingVar.push_back(S.NewVariable() << 1);
    uint32_t carry = S.NewVariable() << 1;

    // lower part
    for (int i = 0; i < divisor; i++)
      for (int j = 0; j < divisor; j++) {
        int ij = i + j;
        S.ResetClause();
        S.NewClause();
        S.AddLiteral(linkingX[i] ^ 1);
        S.AddLiteral(linkingY[j] ^ 1);
        if (ij < divisor) {
          S.AddLiteral(linkingVar[ij]);
          S.AddLiteral(carry);
        } else if (ij == divisor)
          S.AddLiteral(carry);
        else if (ij > divisor)
          S.AddLiteral(linkingVar[ij % divisor]);
        S.CommitClause();
      }

    // upper part
    for (int i = divisor; i < linkingX.size(); i++)
      for (int j = divisor; j < linkingY.size(); j++) {
        int ij = i + j - divisor;
        S.ResetClause();
        S.NewClause();
        S.AddLiteral(linkingX[i] ^ 1);
        S.AddLiteral(linkingY[j] ^ 1);
        if (ij < linkingVar.size()) S.AddLiteral(linkingVar[ij]);
        S.CommitClause();
        //	printf("ij = %lld, linkingVar.size() =
        //%lld\n",ij,linkingVar.size());
        S.ResetClause();
        S.NewClause();
        S.AddLiteral(carry ^ 1);
        S.AddLiteral(linkingX[i] ^ 1);
        S.AddLiteral(linkingY[j] ^ 1);
        if (ij + 1 < linkingVar.size()) S.AddLiteral(linkingVar[ij + 1]);
        S.CommitClause();
      }
  }
  linkingX.clear();
  linkingY.clear();
}

void Encodings::genOgawa(std::vector<long long int> &weights,
                         std::vector<uint32_t> &blockings, long long int &total,
                         long long int divisor, uint32_t zero, uint32_t one,
                         int comp, SATSolverProxy &S,
                         std::vector<uint32_t> &lits,
                         std::vector<uint32_t> &linkingVar, long long int UB) {
  linkingVar.clear();

  std::vector<uint32_t> linkingAlpha;
  std::vector<uint32_t> linkingBeta;

  if (total < divisor) {
    std::vector<uint32_t> linking;
    genBailleux(weights, blockings, total, zero, one, comp, S, lits, linking,
                UB);
    total = linking.size() - 2;
    for (int i = 0; i < divisor; i++)
      if (i < linking.size())
        linkingVar.push_back(linking[i]);
      else
        linkingVar.push_back(zero);
    linkingVar.push_back(one);
    linking.clear();
    //    printf("total = %lld, linkngVar.size() = %d\n",
    //    total,linkingVar.size());
  } else if (blockings.size() == 1) {
    long long int weight = weights[0];
    if (weight < UB) {
      long long int upper = weight / divisor;
      long long int lower = weight % divisor;
      long long int pad = divisor - lower - 1;
      linkingVar.push_back(one);
      for (int i = 0; i < lower; i++) linkingVar.push_back(blockings[0]);
      for (int i = 0; i < pad; i++) linkingVar.push_back(zero);
      linkingVar.push_back(one);
      for (int i = 0; i < upper; i++) linkingVar.push_back(blockings[0]);
      total = weight;
    } else {
      S.ResetClause();
      S.NewClause();
      S.AddLiteral(blockings[0] ^ 1);
      S.CommitClause();
      total = 0;
    }
  } else if (blockings.size() > 1) {
    long long int weightL = 0;
    long long int weightR = 0;
    std::vector<long long int> weightsL, weightsR;
    std::vector<uint32_t> blockingsL, blockingsR;
    const long long int half = total / 2;
    (void)half;  // suppress warning about unused variable
    wbsplit(half, weightL, weightR, weights, blockings, weightsL, blockingsL,
            weightsR, blockingsR);

    genOgawa(weightsL, blockingsL, weightL, divisor, zero, one, comp, S, lits,
             linkingAlpha, UB);
    genOgawa(weightsR, blockingsR, weightR, divisor, zero, one, comp, S, lits,
             linkingBeta, UB);

    weightsL.clear();
    blockingsL.clear();
    weightsR.clear();
    blockingsR.clear();

    genOgawa(weightL, linkingAlpha, weightR, linkingBeta, total, divisor, zero,
             one, comp, S, lits, linkingVar, UB);
  }
  // koshi 2013.11.12
  long long int upper = (UB - 1) / divisor;
  for (long long int i = divisor + upper + 1; i < linkingVar.size(); i++) {
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(linkingVar[i] ^ 1);
    S.CommitClause();
  }
  while (divisor + upper + 2 < linkingVar.size()) linkingVar.pop_back();
}

void Encodings::genOgawa0(std::vector<long long int> &weights,
                          std::vector<uint32_t> &blockings, long long int max,
                          long long int k, long long int &divisor, int comp,
                          SATSolverProxy &S, std::vector<uint32_t> &lits,
                          std::vector<uint32_t> &linkingVar) {
  //  koshi 20140327 assert(max >= TOTALIZER);

  /* koshi 2013.11.11
long long int max0 = max;
*/
  long long int k0 = k;
  long long int odd = 1;
  divisor = 0;
  /* koshi 2013.11.11
while (max0 > 0) {
  divisor++;
  max0 -= odd;
  odd += 2;
}
*/
  while (k0 > 0) {
    divisor++;
    k0 -= odd;
    odd += 2;
  }
  if (_settings->verbosity > 1)
    printf("c max = %lld, divisor = %lld\n", max, divisor);

  // koshi 2013.12.24
  if (divisor <= 2) {
    printf("c divisor is less than or equal to 2 ");
    printf("so we use Warner's encoding, i.e. -card=warn\n");
    _settings->SetEncoding(WARNERS);
    genWarners0(weights, blockings, max, k, comp, S, lits, linkingVar);

  } else {
    // koshi 20140109
    printf("c Ogawa's encoding for Cardinality Constraints\n");

    uint32_t one = S.NewVariable() << 1;
    S.ResetClause();
    S.NewClause();
    S.AddLiteral(one);
    S.CommitClause();
    genOgawa(weights, blockings, max, divisor, one ^ 1, one, comp, S, lits,
             linkingVar, k);
  }
}

// TODO BailW2 K-WTO
void Encodings::genBailleuxW2(std::vector<long long int> &weights,
                              std::vector<uint32_t> &blockings,
                              long long int total, uint32_t zero, uint32_t one,
                              int comp, SATSolverProxy &S,
                              std::vector<uint32_t> &lits,
                              std::vector<uint32_t> &linkingVar,
                              std::vector<long long int> &linkingW,
                              long long int UB) {
  assert(weights.size() == blockings.size());

  linkingVar.clear();
  linkingW.clear();
  bool dvar = (comp == 11) ? false : true;

  std::vector<uint32_t> linkingAlpha;
  std::vector<uint32_t> linkingBeta;

  std::vector<long long int> linkingWA;
  std::vector<long long int> linkingWB;

  if (blockings.size() == 1) {  // koshi 20140121

    // 1個のとき

    long long int weight = weights[0];

    if (weight >= UB) {
      printf("weight(%lld) is over %lld\n", weight, UB);
      exit(1);
    }
    // assert(weight < UB);

    linkingVar.push_back(one);
    linkingW.push_back(0);

    linkingVar.push_back(blockings[0]);
    linkingW.push_back(weights[0]);

  } else if (blockings.size() > 1) {
    // 2個以上のとき

    long long int weightL = 0;
    long long int weightR = 0;
    std::vector<long long int> weightsL, weightsR;
    std::vector<uint32_t> blockingsL, blockingsR;
    const long long int half = total / 2;
    (void)half;  // suppress warning about unused variable 'half'
    // weightsとblockingsを半分に分ける
    wbsplit(half, weightL, weightR, weights, blockings, weightsL, blockingsL,
            weightsR, blockingsR);

    // LEFT
    genBailleuxW2(weightsL, blockingsL, weightL, zero, one, comp, S, lits,
                  linkingAlpha, linkingWA, UB);

    // RIGHT
    genBailleuxW2(weightsR, blockingsR, weightR, zero, one, comp, S, lits,
                  linkingBeta, linkingWB, UB);

    weightsL.clear();
    blockingsL.clear();
    weightsR.clear();
    blockingsR.clear();

    long long int top = ((UB < total) ? UB : total + 1);
    int *table = new int[top];

    table[0] = 1;
    for (int i = 1; i < top; i++) {
      table[i] = 0;
    }

    int a_size = linkingWA.size();
    int b_size = linkingWB.size();

    linkingW.clear();
    linkingVar.clear();

    linkingVar.push_back(one);
    linkingW.push_back(0);

    for (int b = 1; b < b_size; ++b) {
      // 2015 02 07
      if (linkingWB[b] < top) {
        linkingVar.push_back(S.NewVariable() << 1 /*(true,dvar)*/);  //変数生成
        linkingW.push_back(linkingWB[b]);

        //新しく節を生成して追加
        S.ResetClause();
        S.NewClause();
        S.AddLiteral(linkingBeta[b] ^ 1);
        S.AddLiteral(linkingVar[linkingVar.size() - 1]);
        S.CommitClause();

        // printf("[ %d ]" , var(linkingVar[linkingVar.size()-1]));

        table[linkingWB[b]] =
            linkingVar
                .size();  // 1になっていたのをlinkingVar.size()に修正　2015
                          // 01 24

      } else {
        S.ResetClause();
        S.NewClause();
        S.AddLiteral(linkingBeta[b] ^ 1);
        S.CommitClause();
      }
    }

    for (int a = 1; a < a_size; ++a) {
      long long int wa = linkingWA[a];

      if (wa >= top) {
        S.ResetClause();
        S.NewClause();
        S.AddLiteral(linkingAlpha[a] ^ 1);
        S.CommitClause();
        continue;
      }

      for (long long int b = 0; b < b_size; ++b) {
        long long int wb = linkingWB[b];

        if (wa + wb < top) {
          if (table[wa + wb] == 0) {  //新しい重みの和
            linkingVar.push_back(S.NewVariable()
                                 << 1 /*(true,dvar)*/);  ////変数生成
            linkingW.push_back(wa + wb);
            table[wa + wb] =
                linkingVar
                    .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
            // printf("[ %d ]" , var(linkingVar[linkingVar.size()-1]));
          }

          //新しく節を生成して追加
          S.ResetClause();
          S.NewClause();
          S.AddLiteral(linkingAlpha[a] ^ 1);
          S.AddLiteral(linkingBeta[b] ^ 1);
          S.AddLiteral(linkingVar[table[wa + wb] - 1]);
          S.CommitClause();

        } else {
          S.ResetClause();
          S.NewClause();
          S.AddLiteral(linkingAlpha[a] ^ 1);
          S.AddLiteral(linkingBeta[b] ^ 1);
          S.CommitClause();
        }
      }
    }

    delete[] table;
  }

  linkingAlpha.clear();
  linkingBeta.clear();
  linkingWA.clear();
  linkingWB.clear();
}

void Encodings::genBailleuxW20(std::vector<long long int> &weights,
                               std::vector<uint32_t> &blockings,
                               long long int max, long long int k, int comp,
                               SATSolverProxy &S, std::vector<uint32_t> &lits,
                               std::vector<uint32_t> &linkingVar,
                               std::vector<long long int> &linkingWeight) {
  // hayata 2014/12/17
  // printf("\nTOを構築
  // =====================================================\n")

  // printf("\n[bailW]\n");
  if (_settings->verbosity > 0)
    printf("c WTO encoding for Cardinailty Constraints\n");

  uint32_t one = S.NewVariable() << 1;
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(one);
  S.CommitClause();

  // printf("one = %d\n" , var(one)+1);

  genBailleuxW2(weights, blockings, max, one ^ 1, one, comp, S, lits,
                linkingVar, linkingWeight, k);
}

void Encodings::genCCl(uint32_t a, SATSolverProxy &S,
                       std::vector<uint32_t> &lits,
                       int varZero) {  // ogawa 2013/04/02 uemura 20161129
  S.ResetClause();
  S.NewClause();  // lits and varZero defined as global vars
  if (a >> 1 == varZero) {
    if ((a & 1) == 0) return;
  } else
    S.AddLiteral(a);
  S.CommitClause();
}

void Encodings::genCCl(uint32_t a, uint32_t b, SATSolverProxy &S,
                       std::vector<uint32_t> &lits,
                       int varZero) {  // ogawa 2013/04/02 uemura 20161129
  S.ResetClause();
  S.NewClause();  // lits and varZero defined as global vars
  if (a >> 1 == varZero) {
    if ((a & 1) == 0) return;
  } else
    S.AddLiteral(a);
  if (b >> 1 == varZero) {
    if ((b & 1) == 0) return;
  } else
    S.AddLiteral(b);
  S.CommitClause();
}

void Encodings::genCCl(uint32_t a, uint32_t b, uint32_t c, SATSolverProxy &S,
                       std::vector<uint32_t> &lits,
                       int varZero) {  // ogawa 2013/04/02 uemura 20161129
  S.ResetClause();
  S.NewClause();  // lits and varZero defined as global vars
  if (a >> 1 == varZero) {
    if ((a & 1) == 0) return;
  } else
    S.AddLiteral(a);
  if (b >> 1 == varZero) {
    if ((b & 1) == 0) return;
  } else
    S.AddLiteral(b);
  if (c >> 1 == varZero) {
    if ((c & 1) == 0) return;
  } else
    S.AddLiteral(c);
  S.CommitClause();
}

void Encodings::genCCl1(uint32_t a, uint32_t b, uint32_t c, SATSolverProxy &S,
                        std::vector<uint32_t> &lits,
                        int varZero) {  // ogawa 2013/04/02 uemura 20161129
  S.ResetClause();
  S.NewClause();  // lits and varZero defined as global vars
  printf("fe");
  if (a >> 1 == varZero) {
    if ((a & 1) == 0) return;
  } else
    S.AddLiteral(a);
  if (b >> 1 == varZero) {
    if ((b & 1) == 0) return;
  } else
    S.AddLiteral(b);
  if (c >> 1 == varZero) {
    if ((c & 1) == 0) return;
  } else
    S.AddLiteral(c);
  S.CommitClause();
}

void Encodings::genCCl(uint32_t a, uint32_t b, uint32_t c, uint32_t d,
                       SATSolverProxy &S, std::vector<uint32_t> &lits,
                       int varZero) {  // ogawa 2013/04/02 uemura 20161129
  S.ResetClause();
  S.NewClause();  // lits and varZero defined as global vars
  if (a >> 1 == varZero) {
    if ((a & 1) == 0) return;
  } else
    S.AddLiteral(a);
  if (b >> 1 == varZero) {
    if ((b & 1) == 0) return;
  } else
    S.AddLiteral(b);
  if (c >> 1 == varZero) {
    if ((c & 1) == 0) return;
  } else
    S.AddLiteral(c);
  if (d >> 1 == varZero) {
    if ((d & 1) == 0) return;
  } else
    S.AddLiteral(d);

  S.CommitClause();
}

void Encodings::genCCl(uint32_t a, uint32_t b, uint32_t c, uint32_t d,
                       uint32_t e, SATSolverProxy &S,
                       std::vector<uint32_t> &lits,
                       int varZero) {  // ogawa 2013/04/02 uemura 20161129
  S.ResetClause();
  S.NewClause();  // lits and varZero defined as global vars
  if (a >> 1 == varZero) {
    if ((a & 1) == 0) return;
  } else
    S.AddLiteral(a);
  if (b >> 1 == varZero) {
    if ((b & 1) == 0) return;
  } else
    S.AddLiteral(b);
  if (c >> 1 == varZero) {
    if ((c & 1) == 0) return;
  } else
    S.AddLiteral(c);
  if (d >> 1 == varZero) {
    if ((d & 1) == 0) return;
  } else
    S.AddLiteral(d);
  if (e >> 1 == varZero) {
    if ((e & 1) == 0) return;
  } else
    S.AddLiteral(e);
  S.CommitClause();
}

// uemura 20161129
void Encodings::genKWMTO(
    std::vector<long long int> &weights, std::vector<uint32_t> &blockings,
    std::vector<long long int> &weightsTable, int from, int to, int div,
    uint32_t zero, std::vector<uint32_t> &lower,
    std::vector<long long int> &lowerW, std::vector<uint32_t> &upper,
    std::vector<long long int> &upperW, SATSolverProxy &S, long long int ub,
    std::vector<uint32_t> &lits, int varZero) {
  int inputsize = to - from + 1;
  lower.clear();
  lowerW.clear();
  upper.clear();
  upperW.clear();

  if (inputsize == 1) {
    // 1個のとき

    long long int weight = weights[from];

    int low = weight % div;
    int up = weight / div;

    lower.push_back(zero);
    lowerW.push_back(0);

    if (low > 0) {
      lower.push_back(blockings[from]);
      lowerW.push_back(low);
    }

    upper.push_back(zero);
    upperW.push_back(0);

    if (up > 0) {
      upper.push_back(blockings[from]);
      upperW.push_back(up);
    }

  } else {
    int middle = inputsize / 2;
    std::vector<uint32_t> alphaLow;
    std::vector<uint32_t> betaLow;
    std::vector<long long int> WalphaLow;
    std::vector<long long int> WbetaLow;

    std::vector<uint32_t> alphaUp;
    std::vector<uint32_t> betaUp;
    std::vector<long long int> WalphaUp;
    std::vector<long long int> WbetaUp;

    genKWMTO(weights, blockings, weightsTable, from, from + middle - 1, div,
             zero, alphaLow, WalphaLow, alphaUp, WalphaUp, S, ub, lits,
             varZero);

    genKWMTO(weights, blockings, weightsTable, from + middle, to, div, zero,
             betaLow, WbetaLow, betaUp, WbetaUp, S, ub, lits, varZero);

    long long int total = weightsTable[to] - weightsTable[from] + weights[from];

    // LOWERの処理=====================================================================================

    int *tableLOW = new int[div];

    tableLOW[0] = 1;
    for (int i = 1; i < div; i++) {
      tableLOW[i] = 0;
    }

    int a_size = WalphaLow.size();
    int b_size = WbetaLow.size();

    lowerW.clear();
    lower.clear();

    lower.push_back(zero);
    lowerW.push_back(0);

    uint32_t C = S.NewVariable() << 1;

    for (int a = 0; a < a_size; ++a) {
      long long int wa = WalphaLow[a];
      // printf("wa = %d\n",wa);

      for (long long int b = 0; b < b_size; ++b) {
        long long int wb = WbetaLow[b];
        // printf("wb = %d\n",wb);

        long long int wab = (wa + wb) % div;

        if (wa + wb < div) {
          if (tableLOW[wab] == 0) {  //新しい重みの和
            lower.push_back(S.NewVariable() << 1);
            lowerW.push_back(wab);
            tableLOW[wab] =
                lower
                    .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                              // printf("lower.size = %d\n",lower.size());
          }

          genCCl(alphaLow[a] ^ 1, betaLow[b] ^ 1, lower[tableLOW[wab] - 1], C,
                 S, lits, varZero);
          // printf("ClauseLOW[-%d(a=%d) -%d(b=%d) %d c]\n"
          // ,alphaLow[a],a,betaLow[b]^1, b,lower[tableLOW[wab]-1]);//arimura

          /*for(int i = 0 ; i < uint32_ts.size() ; ++i){
                printf("%s%d %s " , sign(uint32_ts[i]) == 1 ? "-" : "" ,
            var(uint32_ts[i]) , i == uint32_ts.size()-1 ? "\n" : "v");
            }*/
        } else if (wab == 0) {
          if (a != 0 || b != 0) {
            genCCl(alphaLow[a] ^ 1, betaLow[b] ^ 1, C, S, lits, varZero);
            // printf("LOwerwab==0\n");
          }
          /*for(int i = 0 ; i < uint32_ts.size() ; ++i){
                printf("%s%d %s " , sign(uint32_ts[i]) == 1 ? "-" : "" ,
            var(uint32_ts[i]) , i == uint32_ts.size()-1 ? "\n" : "v");
            }*/
        } else {  // wa + wb > div

          if (tableLOW[wab] == 0) {  //新しい重みの和

            lower.push_back(S.NewVariable() << 1);
            lowerW.push_back(wab);
            tableLOW[wab] =
                lower
                    .size();  //重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
            // printf("lower.size = %d\n",lower.size());
          }

          genCCl(alphaLow[a] ^ 1, betaLow[b] ^ 1, lower[tableLOW[wab] - 1], S,
                 lits, varZero);
          genCCl(alphaLow[a] ^ 1, betaLow[b] ^ 1, C, S, lits, varZero);
          // printf("ClauseLOW[-%d(a=%d) -%d(b=%d) %d c]\n"
          // ,alphaLow[a],a,betaLow[b]^1, b,lower[tableLOW[wab]-1]);//arimura
          // printf("ClauseLOW[-%d(a=%d) -%d(b=%d) c]\n"
          // ,alphaLow[a],a,betaLow[b]^1, b);//arimura
          /*for(int i = 0 ; i < uint32_ts.size() ; ++i){
                printf("%s%d %s " , sign(uint32_ts[i]) == 1 ? "-" : "" ,
            var(uint32_ts[i]) , i == uint32_ts.size()-1 ? "\n" : "v");
            }*/
        }
      }
    }

    delete[] tableLOW;

    WalphaLow.clear();
    WbetaLow.clear();
    alphaLow.clear();
    betaLow.clear();

    // UPPERの処理=====================================================================================

    // long long int UBU = _min(ub , total)/div;//upperの上限値
    // long long int UBU = _min(total , total)/div + 1;//upperの上限値
    long long int UBU = total / div + 1;  // upperの上限値 uemura 20161129

    int *tableUP = new int[UBU + 1];

    tableUP[0] = 1;
    for (int i = 1; i <= UBU; i++) {
      tableUP[i] = 0;
    }

    a_size = WalphaUp.size();
    b_size = WbetaUp.size();

    upperW.clear();
    upper.clear();

    upper.push_back(zero);
    upperW.push_back(0);

    for (int a = 0; a < a_size; ++a) {
      long long int wa = WalphaUp[a];

      for (long long int b = 0; b < b_size; ++b) {
        long long int wb = WbetaUp[b];

        long long int wab = wa + wb;  //キャリーなしの場合

        if (UBU < wab) {  //超えてる場合

          //新しく節を生成して追加
          genCCl(alphaUp[a] ^ 1, betaUp[b] ^ 1, S, lits, varZero);
          /* for(int i = 0 ; i < uint32_ts.size() ; ++i){
                printf("%s%d %s " , sign(uint32_ts[i]) == 1 ? "-" : "" ,
            var(uint32_ts[i]) , i == uint32_ts.size()-1 ? "\n" : "v");
            }*/
        } else {
          if (wab > 0) {
            if (tableUP[wab] == 0) {  //新しい重みの和
              upper.push_back(S.NewVariable() << 1);
              upperW.push_back(wab);
              tableUP[wab] =
                  upper
                      .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
              // printf("[ %d ]" , var(linkingVar[linkingVar.size()-1]));
            }

            //新しく節を生成して追加
            genCCl(alphaUp[a] ^ 1, betaUp[b] ^ 1, upper[tableUP[wab] - 1], S,
                   lits, varZero);
            /*for(int i = 0 ; i < uint32_ts.size() ; ++i){
                  printf("%s%d %s " , sign(uint32_ts[i]) == 1 ? "-" : "" ,
              var(uint32_ts[i]) , i == uint32_ts.size()-1 ? "\n" : "v");
              }*/
          }
        }

        wab = wa + wb + 1;  //キャリーつきの場合

        if (UBU < wab) {  //超えてる場合

          genCCl(alphaUp[a] ^ 1, betaUp[b] ^ 1, C ^ 1, S, lits, varZero);
          /*for(int i = 0 ; i < uint32_ts.size() ; ++i){
                    printf("%s%d %s " , sign(uint32_ts[i]) == 1 ? "-" : "" ,
             var(uint32_ts[i]) , i == uint32_ts.size()-1 ? "\n" : "v");
                }*/
        } else {
          if (tableUP[wab] == 0) {  //新しい重みの和
            upper.push_back(S.NewVariable() << 1);
            upperW.push_back(wab);
            tableUP[wab] =
                upper
                    .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
          }

          genCCl(alphaUp[a] ^ 1, betaUp[b] ^ 1, C ^ 1, upper[tableUP[wab] - 1],
                 S, lits, varZero);
          if (wab == UBU) {  // test 2015 03 19
            genCCl(upper[tableUP[wab] - 1] ^ 1, C ^ 1, S, lits, varZero);
            // printf("OVER CARRY\n");
          }
          /*for(int i = 0 ; i < uint32_ts.size() ; ++i){
                printf("%s%d %s " , sign(uint32_ts[i]) == 1 ? "-" : "" ,
            var(uint32_ts[i]) , i == uint32_ts.size()-1 ? "\n" : "v");
            }*/
        }
      }
    }

    /*fprintf(stderr , "LOW(%d) ",lowerW.size());
    for(int i = 1 ; i < lowerW.size(); ++i){
        fprintf(stderr , "%lld " , lowerW[i]);

    }
    fprintf(stderr , "\nUP(%d) " , upperW.size());
    for(int i = 1 ; i < upperW.size(); ++i){
        fprintf(stderr , "%lld " , upperW[i]);
    }
    fprintf(stderr , "\n");*/
    // if(carry){

    //}
    // printf("linkingVarUP.size = %d\n" , linkingVarUP.size());

    delete[] tableUP;

    WalphaUp.clear();
    WbetaUp.clear();
    alphaUp.clear();
    betaUp.clear();

    // printf("C = %d " , c>>1);
  }

  /*printf("\nU = { ");
  for(int i = 0 ; i < upper.size() ; ++i)
      printf("%lld " , upperW[i]);

  printf("} L = { ");
  for(int i = 0 ; i < lower.size() ; ++i)
      printf("%lld " , lowerW[i]);

  printf("}");*/

  /*printf("U = { ");
  for(int i = 0 ; i < upper.size() ; ++i)
      printf("%d(%lld) " , var(upper[i]) , upperW[i]);

  printf("} L = { ");
  for(int i = 0 ; i < lower.size() ; ++i)
      printf("%d(%lld) " , var(lower[i]) , lowerW[i]);

  printf("}\n");*/
}

void Encodings::genKWMTO0(
    std::vector<long long int> &weights, std::vector<uint32_t> &blockings,
    long long int max, long long int k, std::vector<long long int> &divisors,
    SATSolverProxy &S, std::vector<uint32_t> &lits,
    std::vector<std::vector<uint32_t>> &linkingVars,
    std::vector<std::vector<long long int>> &linkingWeights) {
  if (_settings->verbosity > 0)
    printf("c WMTO encoding for Cardinailty Constraints\n");

  for (int i = 0; i < 2; i++) {
    linkingVars.resize(linkingVars.size() + 1);
    linkingWeights.resize(linkingWeights.size() + 1);
    //        linkingVars.resize();
    //        linkingWeights.resize();
  }

  divisors.push_back(static_cast<long long>(pow(max, 1.0 / 2.0)));
  if (_settings->verbosity > 1) printf("c p = %lld\n", divisors[0]);

  std::vector<long long int> weightsTable;
  long long int tmp = 0;
  int size = blockings.size();

  for (int i = 0; i < size; ++i) {
    tmp += weights[i];
    weightsTable.push_back(tmp);
  }

  uint32_t zero = S.NewVariable() << 1;
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(zero);
  S.CommitClause();

  int varZero = zero >> 1;

  genKWMTO(weights, blockings, weightsTable, 0, size - 1,
           static_cast<int>(divisors[0]), zero, linkingVars[0],
           linkingWeights[0], linkingVars[1], linkingWeights[1], S, k, lits,
           varZero);
}

// MRWTO UEMURA 20161112
void Encodings::genMRWTO(
    std::vector<long long int> &weights, std::vector<uint32_t> &blockings,
    std::vector<long long int> &weightsTable, int from, int to,
    std::vector<long long int> &divisors, uint32_t zero,
    std::vector<std::vector<uint32_t>> &linkingVars,
    std::vector<std::vector<long long int>> &linkingWeights, SATSolverProxy &S,
    long long int ub, std::vector<uint32_t> &lits, int varZero) {
  int ndigit = linkingVars.size();

  int inputsize = to - from + 1;

  for (int i = 0; i < ndigit; i++) {
    linkingVars[i].clear();
    linkingWeights[i].clear();
  }
  if (inputsize == 1) {
    // 1個のとき

    long long int tmpw = weights[from];
    int *digit = new int[ndigit];

    for (int i = 0; i < ndigit - 1; i++) {
      digit[i] = tmpw % divisors[i];
      tmpw = tmpw / divisors[i];
    }
    digit[ndigit - 1] = tmpw;

    for (int i = 0; i < ndigit; i++) {
      linkingVars[i].push_back(zero);
      linkingWeights[i].push_back(0);
      if (digit[i] > 0) {
        linkingVars[i].push_back(blockings[from]);
        linkingWeights[i].push_back(digit[i]);
      }
    }

    delete[] digit;
  } else {
    int middle = inputsize / 2;

    std::vector<std::vector<uint32_t>> alphalinkingVars;
    std::vector<std::vector<uint32_t>> betalinkingVars;
    std::vector<std::vector<long long int>> alphalinkingWeights;
    std::vector<std::vector<long long int>> betalinkingWeights;
    for (int i = 0; i < ndigit; i++) {
      alphalinkingVars.resize(alphalinkingVars.size() + 1);
      betalinkingVars.resize(betalinkingVars.size() + 1);
      alphalinkingWeights.resize(alphalinkingWeights.size() + 1);
      betalinkingWeights.resize(betalinkingWeights.size() + 1);
    }

    genMRWTO(weights, blockings, weightsTable, from, from + middle - 1,
             divisors, zero, alphalinkingVars, alphalinkingWeights, S, ub, lits,
             varZero);
    genMRWTO(weights, blockings, weightsTable, from + middle, to, divisors,
             zero, betalinkingVars, betalinkingWeights, S, ub, lits, varZero);

    long long int total = weightsTable[to] - weightsTable[from] + weights[from];
    // uint32_t *C = new uint32_t[ndigit];
    std::vector<uint32_t> C;

    for (int cdigit = 0; cdigit < ndigit; cdigit++) {
      //最下位桁の処理=====================================================================================
      if (cdigit == 0) {
        C.resize(C.size() + 1);
        int div = divisors[cdigit];
        int *tableLOW = new int[div];
        tableLOW[0] = 1;
        for (int i = 1; i < div; i++) {
          tableLOW[i] = 0;
        }
        int a_size = alphalinkingWeights[cdigit].size();
        int b_size = betalinkingWeights[cdigit].size();
        linkingVars[cdigit].clear();
        linkingWeights[cdigit].clear();
        linkingVars[cdigit].push_back(zero);
        linkingWeights[cdigit].push_back(0);
        C[cdigit] = S.NewVariable() << 1;

        for (int a = 0; a < a_size; ++a) {
          long long int wa = alphalinkingWeights[cdigit][a];

          for (long long int b = 0; b < b_size; ++b) {
            long long int wb = betalinkingWeights[cdigit][b];
            long long int wab = (wa + wb) % div;
            if (wa + wb < div) {
              if (tableLOW[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableLOW[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVars[cdigit]の何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableLOW[wab] - 1], C[cdigit], S, lits,
                     varZero);

            } else if (wab == 0) {
              if (a != 0 || b != 0) {
                genCCl(alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                       varZero);
              }

            } else {                     // wa + wb > div
              if (tableLOW[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableLOW[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableLOW[wab] - 1], S, lits, varZero);
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                     varZero);
            }
          }
        }
        delete[] tableLOW;
        // alphalinkingWeights[cdigit].clear();
        // betalinkingWeights[cdigit].clear();
        // alphalinkingVars[cdigit].clear();
        // betalinkingVars[cdigit].clear();
      }

      //中間桁の処理====================================================================================================
      else if (cdigit > 0 && cdigit < ndigit - 1) {
        C.resize(C.size() + 1);
        int div = divisors[cdigit];
        int *tableMIDDLE = new int[div];
        tableMIDDLE[0] = 1;
        for (int i = 1; i < div; i++) {
          tableMIDDLE[i] = 0;
        }

        int a_size = alphalinkingWeights[cdigit].size();
        int b_size = betalinkingWeights[cdigit].size();
        linkingVars[cdigit].clear();
        linkingWeights[cdigit].clear();
        linkingVars[cdigit].push_back(zero);
        linkingWeights[cdigit].push_back(0);
        C[cdigit] = S.NewVariable() << 1;
        for (int a = 0; a < a_size; ++a) {
          long long int wa = alphalinkingWeights[cdigit][a];

          for (long long int b = 0; b < b_size; ++b) {
            long long int wb = betalinkingWeights[cdigit][b];
            long long int wab = (wa + wb) % div;

            if (wa + wb < div) {
              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], C[cdigit], S,
                     lits, varZero);

            } else if (wab == 0) {
              if (a != 0 || b != 0) {
                genCCl(alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                       varZero);
              }
            } else {  // wa + wb >= div

              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], S, lits,
                     varZero);
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                     varZero);
            }

            wab = (wa + wb + 1) % div;  // carryがある場合
            if (wa + wb + 1 < div) {
              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
              }
              genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], C[cdigit], S,
                     lits, varZero);

            } else if (wab == 0) {
              if (a != 0 || b != 0) {
                genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                       varZero);
              }

            } else {  // wa + wb + 1 > div

              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
              }
              genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], S, lits,
                     varZero);
              genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                     varZero);
            }
          }
        }
        delete[] tableMIDDLE;
        // alphalinkingWeights[cdigit].clear();
        // betalinkingWeights[cdigit].clear();
        // alphalinkingVars[cdigit].clear();
        // betalinkingVars[cdigit].clear();
      }

      //最上位桁の処理=====================================================================================
      else if (cdigit == ndigit - 1) {
        // long long int UBU = _min(ub , total)/div;
        // //linkingVars[cdigit]の上限値
        long long int UBU = total;  // linkingVars[cdigit]の上限値
        for (int i = 0; i < cdigit; i++) {
          UBU /= divisors[i];
        }
        UBU++;

        int *tableUP = new int[UBU + 1];
        tableUP[0] = 1;
        for (int i = 1; i <= UBU; i++) {
          tableUP[i] = 0;
        }

        int a_size = alphalinkingWeights[cdigit].size();
        int b_size = betalinkingWeights[cdigit].size();
        linkingWeights[cdigit].clear();
        linkingVars[cdigit].clear();
        linkingVars[cdigit].push_back(zero);
        linkingWeights[cdigit].push_back(0);

        for (int a = 0; a < a_size; ++a) {
          long long int wa = alphalinkingWeights[cdigit][a];

          for (long long int b = 0; b < b_size; ++b) {
            long long int wb = betalinkingWeights[cdigit][b];

            long long int wab = wa + wb;  //キャリーなしの場合

            if (UBU < wab) {  //超えてる場合
                              //新しく節を生成して追加
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, S, lits, varZero);
            } else {
              if (wab > 0) {
                if (tableUP[wab] == 0) {  //新しい重みの和
                  linkingVars[cdigit].push_back(S.NewVariable() << 1);
                  linkingWeights[cdigit].push_back(wab);
                  tableUP[wab] =
                      linkingVars[cdigit]
                          .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                }

                //新しく節を生成して追加
                genCCl(alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1,
                       linkingVars[cdigit][tableUP[wab] - 1], S, lits, varZero);
              }
            }

            wab = wa + wb + 1;  //キャリーつきの場合
            if (UBU < wab) {    //超えてる場合
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit - 1] ^ 1, S, lits,
                     varZero);
            } else {
              if (tableUP[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableUP[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit - 1] ^ 1,
                     linkingVars[cdigit][tableUP[wab] - 1], S, lits, varZero);
              if (wab == UBU) {  // test 2015 03 19
                genCCl(linkingVars[cdigit][tableUP[wab] - 1] ^ 1,
                       C[cdigit - 1] ^ 1, S, lits, varZero);
              }
            }
          }
        }
        delete[] tableUP;
      }
    }
    C.clear();
    alphalinkingWeights.clear();
    betalinkingWeights.clear();
    alphalinkingVars.clear();
    betalinkingVars.clear();
  }
}

void Encodings::genMRWTO0(
    std::vector<long long int> &weights, std::vector<uint32_t> &blockings,
    long long int max, long long int k, std::vector<long long int> &divisors,
    SATSolverProxy &S, std::vector<uint32_t> &lits,
    std::vector<std::vector<uint32_t>> &linkingVars,
    std::vector<std::vector<long long int>> &linkingWeights,
    EncodingType encoding) {
  if (_settings->verbosity > 0)
    printf("c MRWTO encoding for Cardinailty Constraints\n");

  int ndigit = 2;  // mrwtoに用いる変数 基数の数 uemura 2016.11.08
  int mrdiv = 2;
  if (encoding == MRWTO) {
    for (long long int i = 512; i < max; i *= 16) {
      ndigit++;
    }
    if (_settings->verbosity > 1) printf("c number of digit = %d\n", ndigit);
  } else if (encoding == MRWTO2) {
    ndigit = 1;
    long long int wtmp = max;
    while (wtmp > mrdiv - 1) {
      wtmp /= mrdiv;
      ndigit++;
    }
    if (_settings->verbosity > 1) printf("c number of digit = %d\n", ndigit);
  }

  // linkingVars,linkingWeightSのサイズの設定
  //  for (int i = 0; i < ndigit; i++) {
  //    linkingVars.resize(linkingVars.size() + 1);
  //    linkingWeights.resize(linkingWeights.size() + 1);
  //  }

  linkingVars.resize(linkingVars.size() + ndigit);
  linkingWeights.resize(linkingWeights.size() + ndigit);

  // uemura MRWTO用のdivisorsの計算
  if (_settings->verbosity > 1) printf("c ");
  if (encoding == MRWTO) {
    for (int k = 0; k < ndigit - 1; k++) {
      divisors.push_back(static_cast<long long int>(pow(max, (1.0 / ndigit))));
      if (_settings->verbosity > 1) printf("p%d= %lld\t", k, divisors[k]);
    }
    if (_settings->verbosity > 1) printf("\n");
  } else if (encoding == MRWTO2) {
    if (_settings->verbosity > 1) printf("c ");
    for (int i = 0; i < ndigit - 1; i++) {
      divisors.push_back(mrdiv);
      if (_settings->verbosity > 1) printf("p%d=\t%lld\t", i, divisors[i]);
    }
    if (_settings->verbosity > 1) printf("\n");
  }

  std::vector<long long int> weightsTable;
  long long int tmp = 0;
  int size = blockings.size();

  for (int i = 0; i < size; ++i) {
    tmp += weights[i];
    weightsTable.push_back(tmp);
  }

  const uint32_t zero = S.NewVariable() << 1;
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(zero);
  S.CommitClause();

  const int varZero = zero >> 1;

  genMRWTO(weights, blockings, weightsTable, 0, size - 1, divisors, zero,
           linkingVars, linkingWeights, S, k, lits, varZero);
}

// MRWTO UEMURA 20161112
void Encodings::genMRWTO19(
    std::vector<long long int> &weights, std::vector<uint32_t> &blockings,
    std::vector<long long int> &weightsTable, int from, int to,
    std::vector<long long int> &divisors, uint32_t zero,
    std::vector<std::vector<uint32_t>> &linkingVars,
    std::vector<std::vector<long long int>> &linkingWeights, SATSolverProxy &S,
    long long int ub, std::vector<uint32_t> &lits, int varZero) {
  int ndigit = linkingVars.size();

  int inputsize = to - from + 1;

  for (int i = 0; i < ndigit; i++) {
    linkingVars[i].clear();
    linkingWeights[i].clear();
  }
  if (inputsize == 1) {
    // 1個のとき

    long long int tmpw = weights[from];
    int *digit = new int[ndigit];

    for (int i = 0; i < ndigit - 1; i++) {
      digit[i] = tmpw % divisors[i];
      tmpw = tmpw / divisors[i];
    }
    digit[ndigit - 1] = tmpw;

    for (int i = 0; i < ndigit; i++) {
      linkingVars[i].push_back(zero);
      linkingWeights[i].push_back(0);
      if (digit[i] > 0) {
        linkingVars[i].push_back(blockings[from]);
        linkingWeights[i].push_back(digit[i]);
      }
    }

    delete[] digit;
  } else {
    int middle = inputsize / 2;

    std::vector<std::vector<uint32_t>> alphalinkingVars;
    std::vector<std::vector<uint32_t>> betalinkingVars;
    std::vector<std::vector<long long int>> alphalinkingWeights;
    std::vector<std::vector<long long int>> betalinkingWeights;

    alphalinkingVars.resize(alphalinkingVars.size() + ndigit);
    betalinkingVars.resize(betalinkingVars.size() + ndigit);
    alphalinkingWeights.resize(alphalinkingWeights.size() + ndigit);
    betalinkingWeights.resize(betalinkingWeights.size() + ndigit);
    //    for (int i = 0; i < ndigit; i++) {
    //      alphalinkingVars.push();
    //      betalinkingVars.push();
    //      alphalinkingWeights.push();
    //      betalinkingWeights.push();
    //    }

    genMRWTO19(weights, blockings, weightsTable, from, from + middle - 1,
               divisors, zero, alphalinkingVars, alphalinkingWeights, S, ub,
               lits, varZero);
    genMRWTO19(weights, blockings, weightsTable, from + middle, to, divisors,
               zero, betalinkingVars, betalinkingWeights, S, ub, lits, varZero);

    long long int total = weightsTable[to] - weightsTable[from] + weights[from];
    // Lit *C = new Lit[ndigit];
    std::vector<uint32_t> C;

    long long int prodiv;
    for (int cdigit = 0; cdigit < ndigit; cdigit++) {
      // CyouRyuu
      if (cdigit != 0) {
        prodiv = 1;
        for (int CR1 = 0; CR1 < cdigit; ++CR1) {
          prodiv = prodiv * divisors[CR1];
        }
      }
      //最下位桁の処理=====================================================================================
      if (cdigit == 0) {
        //        C.push();
        C.resize(C.size() + 1);
        int div = divisors[cdigit];
        int *tableLOW = new int[div];
        tableLOW[0] = 1;
        for (int i = 1; i < div; i++) {
          tableLOW[i] = 0;
        }
        int a_size = alphalinkingWeights[cdigit].size();
        int b_size = betalinkingWeights[cdigit].size();
        linkingVars[cdigit].clear();
        linkingWeights[cdigit].clear();
        linkingVars[cdigit].push_back(zero);
        linkingWeights[cdigit].push_back(0);
        C[cdigit] = S.NewVariable() << 1;

        for (int a = 0; a < a_size; ++a) {
          long long int wa = alphalinkingWeights[cdigit][a];

          for (long long int b = 0; b < b_size; ++b) {
            long long int wb = betalinkingWeights[cdigit][b];
            long long int wab = (wa + wb) % div;
            if (wa + wb < div && wa + wb <= ub) {
              if (tableLOW[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableLOW[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVars[cdigit]の何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableLOW[wab] - 1], C[cdigit], S, lits,
                     varZero);

            } else if (wab == 0 && wa + wb <= ub) {
              if (a != 0 || b != 0) {
                genCCl(alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                       varZero);
              }

            } else if (wa + wb <= ub) {  // wa + wb > div
              if (tableLOW[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableLOW[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableLOW[wab] - 1], S, lits, varZero);
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                     varZero);
            } else {
              genCCl(alphalinkingVars[cdigit][a],
                     betalinkingVars[cdigit][b] ^ 1, S, lits, varZero);
            }
          }
        }
        delete[] tableLOW;
        // alphalinkingWeights[cdigit].clear();
        // betalinkingWeights[cdigit].clear();
        // alphalinkingVars[cdigit].clear();
        // betalinkingVars[cdigit].clear();
      }

      //中間桁の処理====================================================================================================
      else if (cdigit > 0 && cdigit < ndigit - 1) {
        C.resize(C.size() + 1);
        int div = divisors[cdigit];
        int *tableMIDDLE = new int[div];
        tableMIDDLE[0] = 1;
        for (int i = 1; i < div; i++) {
          tableMIDDLE[i] = 0;
        }

        int a_size = alphalinkingWeights[cdigit].size();
        int b_size = betalinkingWeights[cdigit].size();
        linkingVars[cdigit].clear();
        linkingWeights[cdigit].clear();
        linkingVars[cdigit].push_back(zero);
        linkingWeights[cdigit].push_back(0);
        C[cdigit] = S.NewVariable() << 1;
        for (int a = 0; a < a_size; ++a) {
          long long int wa = alphalinkingWeights[cdigit][a];

          for (long long int b = 0; b < b_size; ++b) {
            long long int wb = betalinkingWeights[cdigit][b];
            long long int wab = (wa + wb) % div;

            if (wa + wb < div && (wa + wb) * prodiv <= ub) {
              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], C[cdigit], S,
                     lits, varZero);

            } else if (wab == 0 && (wa + wb) * prodiv <= ub) {
              if (a != 0 || b != 0) {
                genCCl(alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                       varZero);
              }
            } else if ((wa + wb) * prodiv <= ub) {  // wa + wb >= div

              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], S, lits,
                     varZero);
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                     varZero);
            } else {
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, S, lits, varZero);
            }

            wab = (wa + wb + 1) % div;  // carryがある場合
            if (wa + wb + 1 < div && (wa + wb + 1) * prodiv <= ub) {
              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
              }
              genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], C[cdigit], S,
                     lits, varZero);

            } else if (wab == 0 && (wa + wb + 1) * prodiv <= ub) {
              if (a != 0 || b != 0) {
                genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                       varZero);
              }

            } else if ((wa + wb + 1) * prodiv <= ub) {  // wa + wb + 1 > div

              if (tableMIDDLE[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableMIDDLE[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)%divがlinkingVarの何番目に対応するかを記録
              }
              genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1,
                     linkingVars[cdigit][tableMIDDLE[wab] - 1], S, lits,
                     varZero);
              genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit], S, lits,
                     varZero);
            } else {
              genCCl(C[cdigit - 1] ^ 1, alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, S, lits, varZero);
            }
          }
        }
        delete[] tableMIDDLE;
        // alphalinkingWeights[cdigit].clear();
        // betalinkingWeights[cdigit].clear();
        // alphalinkingVars[cdigit].clear();
        // betalinkingVars[cdigit].clear();
      }

      //最上位桁の処理=====================================================================================
      else if (cdigit == ndigit - 1) {
        // long long int UBU = (total<ub?total:ub);
        // //linkingVars[cdigit]の上限値
        long long int UBU = total;
        for (int i = 0; i < cdigit; i++) {
          UBU /= divisors[i];
        }
        UBU++;

        int *tableUP = new int[UBU + 1];
        tableUP[0] = 1;
        for (int i = 1; i <= UBU; i++) {
          tableUP[i] = 0;
        }

        int a_size = alphalinkingWeights[cdigit].size();
        int b_size = betalinkingWeights[cdigit].size();
        linkingWeights[cdigit].clear();
        linkingVars[cdigit].clear();
        linkingVars[cdigit].push_back(zero);
        linkingWeights[cdigit].push_back(0);

        for (int a = 0; a < a_size; ++a) {
          long long int wa = alphalinkingWeights[cdigit][a];

          for (long long int b = 0; b < b_size; ++b) {
            long long int wb = betalinkingWeights[cdigit][b];

            long long int wab = wa + wb;  //キャリーなしの場合

            if (UBU < wab || (wa + wb) * prodiv > ub) {  //超えてる場合
              //新しく節を生成して追加
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, S, lits, varZero);
            } else {
              if (wab > 0) {
                if (tableUP[wab] == 0) {  //新しい重みの和
                  linkingVars[cdigit].push_back(S.NewVariable() << 1);
                  linkingWeights[cdigit].push_back(wab);
                  tableUP[wab] =
                      linkingVars[cdigit]
                          .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
                }

                //新しく節を生成して追加
                genCCl(alphalinkingVars[cdigit][a] ^ 1,
                       betalinkingVars[cdigit][b] ^ 1,
                       linkingVars[cdigit][tableUP[wab] - 1], S, lits, varZero);
              }
            }

            wab = wa + wb + 1;  //キャリーつきの場合
            if (UBU < wab || (wa + wb + 1) * prodiv > ub) {  //超えてる場合
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit - 1] ^ 1, S, lits,
                     varZero);
            } else {
              if (tableUP[wab] == 0) {  //新しい重みの和
                linkingVars[cdigit].push_back(S.NewVariable() << 1);
                linkingWeights[cdigit].push_back(wab);
                tableUP[wab] =
                    linkingVars[cdigit]
                        .size();  //重み(wa+wb)がlinkingVarの何番目に対応するかを記録
              }
              genCCl(alphalinkingVars[cdigit][a] ^ 1,
                     betalinkingVars[cdigit][b] ^ 1, C[cdigit - 1] ^ 1,
                     linkingVars[cdigit][tableUP[wab] - 1], S, lits, varZero);
              if (wab == UBU) {  // test 2015 03 19
                genCCl(linkingVars[cdigit][tableUP[wab] - 1] ^ 1,
                       C[cdigit - 1] ^ 1, S, lits, varZero);
              }
            }
          }
        }
        delete[] tableUP;
      }
    }
    C.clear();
    alphalinkingWeights.clear();
    betalinkingWeights.clear();
    alphalinkingVars.clear();
    betalinkingVars.clear();
  }
}

void Encodings::genMRWTO19_0(
    std::vector<long long int> &weights, std::vector<uint32_t> &blockings,
    long long int max, long long int k, std::vector<long long int> &divisors,
    SATSolverProxy &S, std::vector<uint32_t> &lits,
    std::vector<std::vector<uint32_t>> &linkingVars,
    std::vector<std::vector<long long int>> &linkingWeights,
    EncodingType encoding) {
  if (_settings->verbosity > 0)
    printf("c MRWTO 2019 encoding for Cardinailty Constraints\n");

  // CyouRyuu
  int CR1, CR2;
  long long int *weightsSort = new long long int[blockings.size()];
  for (CR1 = 0; CR1 < blockings.size(); ++CR1) {
    weightsSort[CR1] = weights[CR1];
  }
  quicksort(weightsSort, 0, blockings.size() - 1);
  // upper bound (size limit)
  long long int sizeMax;
  if (k > pow(2, 12)) {
    sizeMax = pow(2, 12);
  } else {
    sizeMax = k;
  }
  if (_settings->verbosity > 0) printf("c sizeMax=%lld\n", sizeMax);
  /*
  //eratosthenes
  //sieve[n]=1 indicates that n is not prime; sieve[n]=0 indicates that n is
  prime. int *sieve = new int[sizeMax+1]; for (CR1 = 0; CR1 < sizeMax+1; ++CR1)
  { sieve[CR1] = 0;
  }
  sieve[0] = sieve[1] = 1;
  int sqrt_MAX = sqrt(sizeMax);
  for (CR1 = 2; CR1 < sqrt_MAX; ++CR1) {
          if (!sieve[CR1]) {
                  for (CR2 = CR1 * CR1; CR2 < sizeMax+1; CR2 += CR1) {
                          sieve[CR2] = 1;
                  }
          }
  }
  */
  // natural number (>2)
  int *natural = new int[sizeMax + 1];
  for (CR1 = 0; CR1 < sizeMax + 1; ++CR1) {
    natural[CR1] = 0;
  }
  natural[0] = natural[1] = 1;
  // prime list completed

  bool pow1IsLastPrime = false;
  int nofDigits = 2;
  int nofDigitsPow1 = 2;
  int primePow1 = 0;
  long long int tmpCR = 0;
  int maxTime = 0;
  int prime = 0;
  int *topCount = new int[sizeMax + 1];
  int *count = new int[sizeMax + 1];
  while (weightsSort[blockings.size() - 1] > 0) {
    // done and rest
    long long int tmProduct = 1;
    for (CR2 = 0; CR2 < divisors.size(); ++CR2) {
      tmProduct = tmProduct * divisors[CR2];
    }
    double restPart = (1.0 * k / tmProduct) + 1;
    // count 0 weights number
    int nofZero = 0;
    while (weightsSort[nofZero] == 0 && nofZero < blockings.size()) {
      nofZero++;
    }
    for (CR1 = 0; CR1 < sizeMax + 1; ++CR1) {
      topCount[CR1] = 0;
    }
    for (CR1 = 0; CR1 < blockings.size(); ++CR1) {
      tmpCR = weightsSort[CR1];
      if (((divisors.size() == 0 && tmpCR <= pow(2, 12)) ||
           (divisors.size() < 3 && tmpCR <= pow(2, 11)) ||
           tmpCR <= pow(2, 10)) &&
          tmpCR != 0) {
        if (CR1 == 0 || (CR1 > 0 && tmpCR != weightsSort[CR1 - 1])) {
          for (CR2 = 0; CR2 < sizeMax + 1; ++CR2) {
            count[CR2] = 0;
          }
          while (tmpCR != 1) {
            for (CR2 = sizeMax; CR2 > 1; --CR2) {
              if (natural[CR2] == 0 &&
                  tmpCR % CR2 == 0) {  // natural[] <-> sieve[]
                tmpCR = tmpCR / CR2;
                count[CR2]++;
                break;
              }
            }
          }
        }
        for (CR2 = 2; CR2 < sizeMax + 1; ++CR2) {
          // If you want to figure out <..> of (if(&& <..>)) below, please see
          // naturNumbDivisor_v2.0 before or other versions if (count[CR2] > 0
          // && divisors.size()+1 >=
          // (log(k)/log(CR2))*(1-1/CR2)-(log(tmProduct)/log(CR2))) {
          if (count[CR2] > 0 &&
              divisors.size() + 1 >= (log(k) / log(CR2)) / CR2) {
            // if (count[CR2] > 0) {
            topCount[CR2]++;
          }
        }
      }
    }
    bool nOkprmie = true;
    while (true) {
      if (!nOkprmie) {
        break;
      }
      maxTime = 0;
      for (CR1 = 2; CR1 < sizeMax + 1; ++CR1) {
        if (topCount[CR1] > maxTime) {
          maxTime = topCount[CR1];
        }
      }
      if (maxTime == 0) {
        break;
      }
      for (CR1 = sizeMax; CR1 > 1; --CR1) {
        if (CR1 <= restPart && topCount[CR1] == maxTime &&
            (log(restPart) / log(CR1)) * maxTime >
                (blockings.size() - nofZero)) {
          // if (topCount[CR1] == maxTime && CR1 <= restPart && 2*maxTime >=
          // (blockings.size()-nofZero)) {
          int tmpDivis, tdM, tdP;
          bool inM, inP;
          nofDigits = 2;
          for (long long int i = 256; i < restPart; i = i * 16) {
            nofDigits++;
          }
          tmpDivis = ((pow(restPart, (1.0 / nofDigits)) + 1) <
                              pow(max, (1.0 / nofDigits))
                          ? (pow(restPart, (1.0 / nofDigits)) + 1)
                          : pow(max, (1.0 / nofDigits)));
          if (CR1 > tmpDivis) {
            tdM = tmpDivis;
            tdP = tmpDivis;
            inM = false;
            inP = false;
            while (CR1 % tdM != 0 && tdM > 1) {
              tdM--;
            }
            if (tdM != 1) inM = true;
            while (CR1 % tdP != 0 && tdP < CR1) {
              tdP++;
            }
            if (tdP != CR1) inP = true;
            int dummyPrime = 1;
            if (inM && inP) {  // both are true
              dummyPrime = ((tmpDivis - tdM) <= (tdP - tmpDivis) ? tdM : tdP);
            } else if (inM && !inP) {  // only tdM is true
              dummyPrime = ((tmpDivis - tdM) <= (CR1 - tmpDivis) ? tdM : CR1);
            } else if (!inM && inP) {  // only tdP is true
              dummyPrime = tdP;
            } else {  // both are false
              dummyPrime = CR1;
            }
            prime = dummyPrime;
          } else {
            prime = CR1;
          }
          if (_settings->verbosity > 0)
            std::cout << "c No." << divisors.size() << " Divisor=" << prime
                      << "\tmaxTime=" << maxTime << "\tnofZero=" << nofZero
                      << "\tdRate="
                      << 100.0 * maxTime / (blockings.size() - nofZero)
                      << "%\tEvl=" << log(k) / log(CR1) / CR1 << std::endl;
          //          printf(
          //              "c No.%d "
          //              "Divisor=%d\tmaxTime=%d\tnofZero=%d\tdRate=%.2f\%\tEvl=%.2f\n",
          //              divisors.size(), prime, maxTime, nofZero,
          //              100.0 * maxTime / (blockings.size() - nofZero),
          //              log(k) / log(CR1) / CR1);
          nOkprmie = false;
          pow1IsLastPrime = false;
          break;
        } else {
          topCount[CR1] = 0;
        }
      }
    }
    if (nOkprmie) {
      if (divisors.size() != 0 && restPart < divisors[0] * divisors[0]) {
        prime = divisors[0];
      } else {
        nofDigits = 2;
        for (long long int i = 256; i < restPart; i = i * 16) {
          nofDigits++;
        }
        prime = ((pow(restPart, (1.0 / nofDigits)) + 1) <
                         pow(max, (1.0 / nofDigits))
                     ? (pow(restPart, (1.0 / nofDigits)) + 1)
                     : pow(max, (1.0 / nofDigits)));
      }
      if (prime <= 1) {
        if (_settings->verbosity > 0)
          printf("c Divisor can not be less than 2 (pow1)\n");
      } else {
        if (_settings->verbosity > 0)
          printf("c No.%zu Divisor=%d\tgenerated by pow1\n", divisors.size(),
                 prime);
        nofDigitsPow1 = nofDigits;
        primePow1 = prime;
        pow1IsLastPrime = true;
      }
    }
    if (prime <= 1) {
      break;
    }
    divisors.push_back(prime);
    for (CR1 = 0; CR1 < blockings.size(); ++CR1) {
      weightsSort[CR1] = weightsSort[CR1] / prime;
    }
    quicksort(weightsSort, 0, blockings.size() - 1);
    prime = 0;
  }

  long long int productPrime = 1;
  if (encoding == MRWTO19) {
    for (int CRk = 0; CRk < divisors.size(); ++CRk) {
      long long int in2Divisors = divisors[CRk];
      productPrime = productPrime * in2Divisors;
    }
  }
  if (productPrime < k) {
    if (pow1IsLastPrime) {
      for (int CRk = 0; CRk < nofDigitsPow1 - 2; CRk++) {
        if (_settings->verbosity > 0)
          printf("c No.%zu Divisor=%d\tgenerated by pow1\n", divisors.size(),
                 primePow1);
        divisors.push_back(primePow1);
      }
    } else {
      double restPart = (1.0 * k / productPrime) + 1;
      nofDigits = 2;
      for (long long int i = 256; i < restPart; i = i * 16) {
        nofDigits++;
      }
      prime =
          ((pow(restPart, (1.0 / nofDigits)) + 1) < pow(max, (1.0 / nofDigits))
               ? (pow(restPart, (1.0 / nofDigits)) + 1)
               : pow(max, (1.0 / nofDigits)));
      for (int CRk = 0; CRk < nofDigits - 1; CRk++) {
        if (_settings->verbosity > 0)
          printf("c No.%zu Divisor=%d\tgenerated by pow2\n", divisors.size(),
                 prime);
        divisors.push_back(prime);
      }
    }
  }

  if (_settings->verbosity > 0) {
    printf("c ");
    for (int CRk = 0; CRk < divisors.size(); ++CRk) {
      printf("p%d= %lld\t", CRk, divisors[CRk]);
    }
    printf("\n");
    if (encoding == MRWTO19) {
      printf("c number of digits = %lu\n", divisors.size() + 1);
    }
  }

  // linkingVars,linkingWeightSのサイズの設定
  //  for (int i = 0; i < divisors.size() + 1; i++) {
  //    linkingVars.push_back();
  //    linkingWeights.push_back();
  //  }

  linkingVars.resize(linkingVars.size() + divisors.size() + 1);
  linkingWeights.resize(linkingWeights.size() + divisors.size() + 1);

  std::vector<long long int> weightsTable;
  long long int tmp = 0;
  int size = blockings.size();

  for (int i = 0; i < size; ++i) {
    tmp += weights[i];
    weightsTable.push_back(tmp);
  }

  const uint32_t zero = S.NewVariable() << 1;
  S.ResetClause();
  S.NewClause();
  S.AddLiteral(zero);
  S.CommitClause();

  const int varZero = zero >> 1;

  genMRWTO19(weights, blockings, weightsTable, 0, size - 1, divisors, zero,
             linkingVars, linkingWeights, S, k, lits, varZero);
}

// koshi 20140106 based on minisat2-070721/maxsat0.2e
// Lit -> mkLit, addClause -> addClause_
// koshi 2013.05.23
long long int Encodings::sumWeight(std::vector<long long int> &weights) {
  long long int sum = 0;
  for (size_t i = 0; i < weights.size(); i++) sum += weights[i];
  return sum;
}

// CyouRyuu (quicksort)
long long int Encodings::med3(long long int x, long long int y,
                              long long int z) {
  if (x < y) {
    if (y < z)
      return y;
    else if (z < x)
      return x;
    else
      return z;
  } else {
    if (z < y)
      return y;
    else if (x < z)
      return x;
    else
      return z;
  }
}

void Encodings::quicksort(long long int a[], int left, int right) {
  if (left < right) {
    int i = left, j = right;
    long long int tmp, pivot = med3(a[i], a[i + (j - i) / 2], a[j]);
    while (1) {
      while (a[i] < pivot) i++;
      while (pivot < a[j]) j--;
      if (i >= j) break;
      tmp = a[i];
      a[i] = a[j];
      a[j] = tmp;
      i++;
      j--;
    }
    quicksort(a, left, i - 1);
    quicksort(a, j + 1, right);
  }
}
} // Namespace Pacose