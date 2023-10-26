/********************************************************************************************
Glucose421SolverProxy.h -- Copyright (c) 2020, Tobias Paxian

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

#ifndef GLUCOSE421SOLVERPROXY_H
#define GLUCOSE421SOLVERPROXY_H

//#include "glucose421/SimpSolver.h"
#include <fstream>  // std::ofstream
#include "SATSolverProxy.h"
#include "glucose-4.2.1/core/Solver.h"

class Glucose421SolverProxy : public SATSolverProxy {
 public:
  Glucose421SolverProxy(void);
  ~Glucose421SolverProxy();

  // pure virtual functions - have to be implemented in the different
  // solverProxy's.
  SATSolverType GetSATSolverType(void);

  uint32_t GetModel(int var);
  int NewVariable();

  void NewClause();
  void AddLiteral(uint32_t *lit);
  bool CommitClause();
  void ResetClause();

  bool AddClause(std::vector<uint32_t> &clause);

  void AddAssumption(uint32_t *lit);
  void ClearAssumption();

  // hard assumptions are sticky - even if clear assumptions is called, hard
  // assumptions are still valid!!!
  void AddHardAssumption(uint32_t *lit);
  void ClearHardAssumption();

  uint32_t Solve();

  void Reset(void);

  // the following functions can be implemented
  // if not it's either translated to the pure virtual functions
  // or it's just not necessary (some extended functionality of SAT solvers)
  uint32_t GetNumberOfClauses();
  uint32_t GetNumberOfVariables();

  //  void AddLiteral(int* lit);
  //  void AddLiteral(int* lit, bool sign);
  //    void AddVariablePrio(uint32_t variable, uint32_t prio);

  //    virtual void SetFrozen(int variable);

  //  uint32_t GetClauseCounter();
  //  uint32_t GetVariableCounter();
  //  void ResetCounter();

  uint32_t Simplify();

  void SetPropagationBudget(int64_t propagationBudget);
  uint32_t SolveLimited();

 protected:
  Glucose421::Solver* _glucose421;
  //    Glucose::Solver* _glucose421;
  Glucose421::vec<Glucose421::Lit> _currentClause;
  Glucose421::vec<Glucose421::Lit> _currentAssumptions;
  Glucose421::vec<Glucose421::Lit> _currentHardAssumptions;
  //  int64_t _propagationBudget;

  //    Minisat::Lit ToLit(uint32_t lit);
};

#endif  // GLUCOSE421SOLVERPROXY_H
