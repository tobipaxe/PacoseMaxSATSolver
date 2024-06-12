/********************************************************************************************
CadicalSolverProxy.h -- Copyright (c) 2020, Tobias Paxian

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
********************************************************************************************/

#ifndef CADICALSOLVERPROXY_H
#define CADICALSOLVERPROXY_H

using namespace std;

#include "SATSolverProxy.h"
#include "cadical/src/cadical.hpp"

class CadicalSolverProxy : public SATSolverProxy {
 public:
  CadicalSolverProxy(void);
  ~CadicalSolverProxy(void);

  SATSolverType GetSATSolverType(void);
  uint32_t GetModel(int var);
  std::vector<int>* GetWholeModel() {return &_model;};
  int NewVariable();
  void NewVariables(uint32_t number);

  void NewClause();
  void AddLiteral(int *lit);
  void AddLiteral(uint32_t *lit);
  bool CommitClause();
  void ResetClause();

  void AddAssumption(uint32_t *lit);
  void ClearAssumption();

  uint32_t Solve();

  void Reset(void);

  void SetFrozen(int variable);
  void MeltFrozen(int variable);

  uint32_t GetNumberOfVariables();
  uint32_t GetNumberOfClauses();

  void SaveWholeModel();
  void SetPropagationBudget(int64_t propagationBudget);
  uint32_t SolveLimited();
  void Phase(uint32_t *lit);

 protected:
  CaDiCaL::Solver *_cadical;

  uint32_t _vars;
  uint32_t _noClauses;
  bool _hasVars;
  //    int _assumption;
  std::vector<int> _assumptions;
  std::vector<int> _model;
};

#endif /* CADICALSOLVERPROXY_H */
