/********************************************************************************************
SATSolverProxy.h -- Copyright (c) 2020, Tobias Paxian

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

#ifndef SATSOLVERPROXY_H
#define SATSOLVERPROXY_H

#include <cstdint>
#include <string>
#include <vector>

/**
 * @brief The SATSolverType enum    which solver is in use
 */
enum class SATSolverType {
  ANTOM = 0,
  GLUCOSE3 = 1,
  GLUCOSE421 = 2,
  CRYPTOMINISAT = 3,
  CADICAL = 4,
  LINGELING = 5,
  MAPLEGLUCOSE = 6,
};


/**
 * @brief SATSolverProxy::InitSATSolver should be a simple SAT solver proxy
 *        only really necessary functions are forced to be implemented in the
 * inherited class
 * @param solverType
 * @return
 */
class SATSolverProxy {
 public:
  SATSolverProxy(uint32_t noOfThreads = 1)
      : _noSolverCalls(0),
        _noOfThreads(noOfThreads),
        _memoryLimit(-1),
        _cpuLimit(-1) {}

  virtual ~SATSolverProxy(){};

  /**
   * @brief InitSATSolver Initializes new proxy
   * @param solverType    which solver to use
   * @return              returns instance of solver
   */
  static SATSolverProxy *InitSATSolver(SATSolverType solverType,
                                       uint32_t noOfThreads = 1);

  // pure virtual functions - have to be implemented in the different
  // solverProxy's.
  /**
   * @brief GetSATSolverType
   * @return                  returns enum which SAT solver is used
   */
  virtual SATSolverType GetSATSolverType() = 0;

  /**
   * @brief GetModel  gets a model for this variable
   * @param var
   * @return          model for this variable, 2x+{0,1}
   */
  virtual uint32_t GetModel(int var) = 0;
  /**
   * @brief NewVariable   adds a new var to the solver
   * @return              if the new variable is successfully added to the
   * solver
   */
  virtual int NewVariable() = 0;

  /**
   * @brief NewClause     test if clause is empty
   */
  virtual void NewClause() = 0;

  /**
   * @brief AddLiteral    adds a new literal to the clause
   * @param lit           if even than it is positive literal
   *                      if odd than it is negative
   *                      the variable is \lfloor lit/2 \rfloor
   *                      sign is lit%2
   */
  virtual void AddLiteral(uint32_t *lit) = 0;
  void AddLiteral(uint32_t lit);

  /**
   * @brief CommitClause  commits the given clause to the solver
   */
  virtual bool CommitClause() = 0;

  /**
   * @brief ResetClause   empties the clause vector such that a new clause can
   * be opened
   */
  virtual void ResetClause() = 0;

  /**
   * @brief AddAssumption     adds another assumption
   * @param lit
   */
  virtual void AddAssumption(uint32_t *lit) = 0;

  /**
   * @brief ClearAssumption   clears assumption vector
   */
  virtual void ClearAssumption() = 0;

  /**
   * @brief Solve         solves the given cnf instance
   * @param assumptions
   * @return
   */
  virtual uint32_t Solve() = 0;

  /**
   * @brief Simplify      removes already SAT clauses
   *                      CURRENTLY ONLY GLUCOSE4
   * @param assumptions
   * @return
   */
  virtual uint32_t Simplify();

  /**
   * @brief Reset     reset the whole solver and start over again from 0
   *                  can be more effective if reset is done by only resetting
   * all variables
   */
  virtual void Reset(void) = 0;

  /**
   * @brief AddClause     adds directly a whole clause
   * @param clause - as in Glucose / Antom (2*var)^sign
   * @return successful added
   */
  virtual bool AddClause(std::vector<uint32_t> &clause);

  /**
   * @brief AddClause     adds directly a whole clause
   * @param clause - as in DIMACS format +-var
   * @return successful added?
   */
  virtual bool AddClause(std::vector<int> &clause);

  // the following functions can be implemented
  // if not it's either translated to the pure virtual functions
  // or it's just not necessary (some extended functionality of SAT solvers)
  /**
   * @brief GetNumberOfClauses
   * @return number of clauses actually in the clause DB
   */
  virtual uint32_t GetNumberOfClauses();

  /**
   * @brief GetNumberOfVariables
   * @return number of variables actually in the clause DB
   */
  virtual uint32_t GetNumberOfVariables();

  /**
   * @brief AddLiteral    negative numbers are negated literals
   * @param lit
   */
  virtual void AddLiteral(int *lit);

  /**
   * @brief AddLiteral    the sign indicates wether a number is negated
   * @param lit
   * @param sign
   */
  virtual void AddLiteral(int *lit, bool sign);

  /**
   * @brief AddAssumptions
   * @param lit
   */
  virtual void AddAssumptions(std::vector<uint32_t> &clause);

  /**
   * @brief AddVariablePriority   antom specific function to add a decision
   * priority to a variable
   * @param variable
   * @param prio
   */
  virtual void AddVariablePrio(uint32_t variable, uint32_t prio);

  virtual void SetNumberOfThreads(uint32_t n);

  virtual void SetFrozen(int variable);

  /**
   * @brief NewVariables  adds number new variables to the solver
   *                      some solver like cms have this directly implemented
   * @param number
   * @return
   */
  virtual void NewVariables(uint32_t number);

  /**
   * @brief SetTimeLimit  set the time limit variable
   *
   * @param limit
   * @return
   */
  virtual void SetTimeLimit(double limit);

  /**
   * @brief SetMemoryLimit    set the memory limit variable
   *
   * @param limit
   * @return
   */
  virtual void SetMemoryLimit(uint32_t limit);

  /**
   * @brief SetTimeLimit  enable the time limit by reading the time limit
   * variable
   * @param limit
   * @return
   */
  virtual void EnableTimeLimit();

  /**
   * @brief SetMemoryLimit    enable the memory limit by reading the memory
   * limit variable
   * @param limit
   * @return
   */
  virtual void EnableMemoryLimit();

  //  // at the moment only for glucose 421 implemented - only for testing
  //  reasons virtual uint32_t GetVariableCounter(); virtual uint32_t
  //  GetClauseCounter(); virtual void ResetCounter();

  SATSolverType ReturnSolverType(const std::string &solver);
  //    void SaveCNF() { saveCNF = true; };
  std::vector<std::vector<uint32_t>> CNF;

  // at the moment only CMS
  virtual void SetReconf(uint32_t reconf);

  virtual void SetPropagationBudget(int64_t propagationBudget);
  virtual uint32_t SolveLimited();

  // at the moment only for glucose 421
  virtual void AddHardAssumption(uint32_t *lit);
  virtual void ClearHardAssumption();
  virtual void Phase(uint32_t *lit);

  uint32_t _noSolverCalls;

 protected:
  //  bool saveCNF;

  uint32_t _noOfThreads;
  int32_t _memoryLimit;
  double _cpuLimit;
};

#endif  // SATSOLVERPROXY_H
