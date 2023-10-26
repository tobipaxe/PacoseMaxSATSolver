/********************************************************************************************
Glucose421SolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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

#include "Glucose421SolverProxy.h"
#include "SATSolverProxy.h"
//#include "glucose-4.2.1/core/Solver.h"
#include "iostream"

Glucose421SolverProxy::Glucose421SolverProxy()
    : _glucose421(new Glucose421::Solver) {
  Reset();
}

Glucose421SolverProxy::~Glucose421SolverProxy() { delete _glucose421; }

SATSolverType Glucose421SolverProxy::GetSATSolverType() {
  return SATSolverType::GLUCOSE421;
}

uint32_t Glucose421SolverProxy::GetModel(int var) {
  int varLbool = toInt(_glucose421->modelValue(var));
  uint32_t rv;
  //    std::cout << "varLbool " << varLbool << std::endl;
  if (varLbool == 0) {
    rv = static_cast<uint32_t>(var << 1);
  } else if (varLbool == 1) {
    rv = static_cast<uint32_t>(var << 1 ^ 1);
  } else if (varLbool == 2) {
    std::cout << "Strange Value undef for variable: " << var << std::endl;
    //         undefined means, both values are possible!
    rv = static_cast<uint32_t>(var);
  } else {
    std::cout << "Strange value vor varLbool: " << varLbool << std::endl;
    std::cout << "Value var: " << var << std::endl;
    assert(false);
  }
  return rv;
}

int Glucose421SolverProxy::NewVariable() {
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
  int v = _glucose421->newVar();
  // std::cout << "Variable: " << v << std::endl;
  return v;
}

void Glucose421SolverProxy::NewClause() { _currentClause.clear(); }

void Glucose421SolverProxy::AddLiteral(uint32_t *lit) {
  // TODO: Ask why toLit is called when SATSolverProxy already converts them
  //  std::cout << "AddLit: " << *lit << "  CC.size(): " <<
  //  _currentClause.size()
  //            << std::endl;
  _currentClause.push(Glucose421::toLit(static_cast<int>(*lit)));
}

bool Glucose421SolverProxy::CommitClause() {
  //  std::cout << "Commit Clause: ";
  //  for (int i = 0; i < _currentClause.size(); i++) {
  //    std::cout << Glucose421::toInt(_currentClause[i]) << ", ";
  //  }
  //  std::cout << std::endl;
  //  if (saveCNF) {
  //    std::vector<uint32_t> clause;
  //    for (int i = 0; i < _currentClause.size(); i++) {
  //      clause.push_back(Glucose421::toInt(_currentClause[i]));
  //    }
  //    CNF.push_back(clause);
  //  }

  // ATTENTION - just for now, to skip current problem!
  if (_currentClause.size() == 0) {
    std::cout << "ERROR: EMPTY CLAUSE - shouldn't occur! - For now it is just "
                 "continued as if nothing had happend! - should be fixed!!!"
              << std::endl;
    return true;
  }

  return _glucose421->addClause(_currentClause);
}

void Glucose421SolverProxy::ResetClause() { _currentClause.clear(); }

bool Glucose421SolverProxy::AddClause(std::vector<uint32_t> &clause) {
  Glucose421::vec<Glucose421::Lit> tmpClause;
  for (uint32_t literal : clause) {
    tmpClause.push(Glucose421::toLit(static_cast<int>(literal)));
    // std::cout << Glucose421::toInt(literal) << " ";
  }
  // std::cout << std::endl;
  return _glucose421->addClause(tmpClause);
}

void Glucose421SolverProxy::AddAssumption(uint32_t *lit) {
  //  std::cout << Glucose421::toInt(*lit) << ";; ";
  _currentAssumptions.push(Glucose421::toLit(static_cast<int>(*lit)));
}

void Glucose421SolverProxy::ClearAssumption() { _currentAssumptions.clear(); }

void Glucose421SolverProxy::AddHardAssumption(uint32_t *lit) {
  //  std::cout << Glucose421::toInt(*lit) << ";; ";
  _currentHardAssumptions.push(Glucose421::toLit(static_cast<int>(*lit)));
}

void Glucose421SolverProxy::ClearHardAssumption() {
  // std::cout << __PRETTY_FUNCTION__ << "  CURRENT HARD ASSUMPTIONS.SIZE: "
  //           << _currentHardAssumptions.size() << std::endl;
  _currentHardAssumptions.clear();
}

uint32_t Glucose421SolverProxy::Solve() {
  _noSolverCalls++;
  EnableTimeLimit();
  EnableMemoryLimit();

  //  for (int i = 0; i < _currentAssumptions.size(); i++) {
  //    std::cout << GLUCOSE421::toInt(_currentAssumptions[i]) << "; ";
  //  }
  //  std::cout << std::endl;
  if (_currentHardAssumptions.size() == 0) {
    // return _glucose421->solve(_currentAssumptions) ? 10 : 20;
    // uint32_t rv = _glucose421->solve(_currentAssumptions) ? 10 : 20;
    // std::cout << "SAT SOLVER SOLVE RESULT: " << rv << std::endl;
    return _glucose421->solve(_currentAssumptions) ? 10 : 20;
  } else {
    Glucose421::vec<Glucose421::Lit> newAssumptions;
    for (int i = 0; i < _currentHardAssumptions.size(); i++)
      newAssumptions.push(_currentHardAssumptions[i]);
    for (int i = 0; i < _currentAssumptions.size(); i++)
      newAssumptions.push(_currentAssumptions[i]);

    // uint32_t rv = _glucose421->solve(newAssumptions) ? 10 : 20;
    // std::cout << "SAT SOLVER WITH HARD ASSUMPTIONS SOLVE RESULT: " << rv <<
    // std::endl;
    return _glucose421->solve(newAssumptions) ? 10 : 20;
  }
}

uint32_t Glucose421SolverProxy::Simplify() {
  return _glucose421->simplify() ? 10 : 20;
}

uint32_t Glucose421SolverProxy::SolveLimited() {
  _noSolverCalls++;
  EnableTimeLimit();
  EnableMemoryLimit();
  int cr = 99;
  if (_currentHardAssumptions.size() == 0) {
    cr = toInt(_glucose421->solveLimited(_currentAssumptions));
  } else {
    Glucose421::vec<Glucose421::Lit> newAssumptions;
    for (int i = 0; i < _currentHardAssumptions.size(); i++)
      newAssumptions.push(_currentHardAssumptions[i]);
    for (int i = 0; i < _currentAssumptions.size(); i++)
      newAssumptions.push(_currentAssumptions[i]);

    cr = toInt(_glucose421->solveLimited(newAssumptions));
  }

  //  _glucose421->printIncrementalStats();
  if (cr == 0)
    return 10;
  else if (cr == 1)
    return 20;
  else if (cr == 2)
    return 0;
  else
    return 99;
}

void Glucose421SolverProxy::SetPropagationBudget(int64_t propagationBudget) {
  _glucose421->setPropBudget(propagationBudget);
}

void Glucose421SolverProxy::Reset() {
  delete _glucose421;
  _glucose421 = new Glucose421::Solver;
  //    _glucose421 = new Glucose421::Solver;
  NewVariable();
  ClearAssumption();
  ResetClause();
}

// void Glucose421SolverProxy::ResetCounter() {
//  _variableCounter = 0;
//  _clauseCounter = 0;
//}

// uint32_t Glucose421SolverProxy::GetVariableCounter() {
//  return _variableCounter;
//}
// uint32_t Glucose421SolverProxy::GetClauseCounter() { return _clauseCounter;
// }

// void Glucose421SolverProxy::SetFrozen(int variable)
//{
//    _glucose421->setFrozen(variable, true);
//}

uint32_t Glucose421SolverProxy::GetNumberOfVariables() {
  return static_cast<uint32_t>(_glucose421->nVars());
}

uint32_t Glucose421SolverProxy::GetNumberOfClauses() {
  return static_cast<uint32_t>(_glucose421->nClauses());
}
