/********************************************************************************************
SATSolverProxy.cpp -- Copyright (c) 2020, Tobias Paxian

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

#include "SATSolverProxy.h"

// list of all solverProxys
//#include "AntomSolverProxy.h"
//#include "MapleGlucoseSolverProxy.h"
#include <csignal>
#include <iostream>
#include <sys/resource.h>
//#include "CadicalSolverProxy.h"
// #include "CryptoMiniSat5SolverProxy.h"
//#include "Glucose3SolverProxy.h"
#include "Glucose421SolverProxy.h"

SATSolverProxy *SATSolverProxy::InitSATSolver(SATSolverType solverType,
                                              unsigned int noOfThreads) {
  //  std::cout << "Create new SAT solver." << std::endl;
  SATSolverProxy *newProxy = nullptr;
  switch (solverType) {
    //        case SATSolverType::ANTOM :
    //        {
    //            SATSolverProxy* newProxy = new AntomSolverProxy();
    //            std::cout << "New AntomSolverProxy: " << newProxy <<
    //            std::endl;

    //            return newProxy;
    //        }
  case SATSolverType::GLUCOSE421: {
    newProxy = new Glucose421SolverProxy();
    //      std::cout << "New Glucose421SolverProxy: " << newProxy <<
    //      std::endl;

    break;
  }
    //        case SATSolverType::MAPLEGLUCOSE :
    //        {
    //            SATSolverProxy* newProxy = new MapleGlucoseSolverProxy();
    //            std::cout << "New MapleGlucoseSolverProxy: " << newProxy <<
    //            std::endl;

    //            return newProxy;
    //        }
    // case SATSolverType::CRYPTOMINISAT: {
    // newProxy = new CryptoMiniSATSolverProxy();
    //      std::cout << "New CryptoMiniSAT ";
    //      newProxy->GetSATSolverType();
    //      std::cout << " SolverProxy: " << newProxy << std::endl;

    //   break;
    // }
    //    case SATSolverType::GLUCOSE3: {
    //      newProxy = new Glucose3SolverProxy();
    //      //      std::cout << "New Glucose3SolverProxy: " << newProxy <<
    //      std::endl;

    //      break;
    //    }
  case SATSolverType::CADICAL: {
    //      newProxy = new CadicalSolverProxy();
    //      std::cout << "New CadicalSolverProxy: " << newProxy << std::endl;
    break;
  }
  default: {
    std::cout << "Solver type" << static_cast<int>(solverType)
              << " currently not available!" << std::endl;
    exit(1);
  }
  }
  return newProxy;
}

SATSolverType SATSolverProxy::ReturnSolverType(const std::string &solver) {
  if (solver == "antom")
    return SATSolverType::ANTOM;
  if (solver == "glucose3")
    return SATSolverType::GLUCOSE3;
  if (solver == "glucose421")
    return SATSolverType::GLUCOSE421;
  if (solver == "cryptominisat")
    return SATSolverType::CRYPTOMINISAT;
  if (solver == "cms")
    return SATSolverType::CRYPTOMINISAT;
  if (solver == "cadical")
    return SATSolverType::CADICAL;
  if (solver == "mapleglucose")
    return SATSolverType::MAPLEGLUCOSE;
  return SATSolverType::CRYPTOMINISAT;
}

unsigned int SATSolverProxy::GetNumberOfClauses() {
  //  std::cout << __PRETTY_FUNCTION__ << std::endl;
  //  std::cout << "This function is not yet implemented for the chosen solver!"
  //            << std::endl;
  return 0;
}

unsigned int SATSolverProxy::GetNumberOfVariables() {
  //  std::cout << __PRETTY_FUNCTION__ << std::endl;
  //  std::cout << "This function is not yet implemented for the chosen solver!"
  //            << std::endl;
  return 0;
}

void SATSolverProxy::AddLiteral(int *lit) {
  unsigned int uLit;
  if (*lit < 0) {
    uLit = static_cast<unsigned int>(((-*lit) << 1) ^ 1);
  } else {
    uLit = static_cast<unsigned int>(*lit << 1);
  }
  AddLiteral(&uLit);
}

void SATSolverProxy::AddLiteral(int *lit, bool sign) {
  unsigned int uLit;
  uLit = static_cast<unsigned int>((*lit << 1) ^ !sign);

  AddLiteral(&uLit);
}

bool SATSolverProxy::AddClause(std::vector<int> &clause) {
  ResetClause();
  NewClause();
  for (int literal : clause) {
    AddLiteral(&literal);
  }
  bool rV = CommitClause();
  ResetClause();

  return rV;
}

void SATSolverProxy::AddLiteral(unsigned int lit) { AddLiteral(&lit); }

unsigned int SATSolverProxy::Simplify() {
  std::cout << "c Simplification is not yet possible for this SAT solver! "
               "Simplification currently only for Glucose4!"
            << std::endl;
  return 10;
}

bool SATSolverProxy::AddClause(std::vector<unsigned int> &clause) {
  //    std::cout << "CB: " << clause.back() << std::endl;
  //  std::cout << "uint" << std::endl;
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  ResetClause();
  NewClause();
  for (unsigned literal : clause) {
    AddLiteral(&literal);
  }
  bool rV = CommitClause();
  ResetClause();
  // std::cout << "return: " << rV << std::endl;
  return rV;
}

void SATSolverProxy::AddAssumptions(std::vector<unsigned int> &assumptions) {
  for (auto literal : assumptions) {
    AddAssumption(&literal);
  }
}

void SATSolverProxy::AddVariablePrio(unsigned int /*variable*/,
                                     unsigned int /*prio*/) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  //    std::cout << "This function is not yet implemented for the chosen
  //    solver!" << std::endl;
}

void SATSolverProxy::SetNumberOfThreads(unsigned int /*n*/) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout
      << "This is not possible for this SAT solver! Only one core is used!"
      << std::endl;
}

void SATSolverProxy::SetFrozen(int /*variable*/) {
  //    std::cout << __PRETTY_FUNCTION__ << std::endl;
  //    std::cout << "This function is not yet implemented for the chosen
  //    solver!" << std::endl;
}

void SATSolverProxy::NewVariables(unsigned number) {
  for (unsigned i = 0; i < number; i++) {
    NewVariable();
  }
}

void SIG_Exit(int signum) {
  switch (signum) {
  case SIGXCPU:
    std::cout << "CPU time limit exceeded." << std::endl;
    exit(signum);
  case SIGINT:
    std::cout << "Memory limit exceeded" << std::endl;
    exit(signum);
  }
}

void SATSolverProxy::SetTimeLimit(double limit) { _cpuLimit = limit; }

void SATSolverProxy::SetMemoryLimit(uint32_t limit) {
  _memoryLimit = static_cast<int>(limit);
}

void SATSolverProxy::EnableTimeLimit() {
  if (_cpuLimit < 0) {
    return;
  }
  signal(SIGXCPU, SIG_Exit);
  if (_cpuLimit != 0 && _cpuLimit != INT32_MAX) {
    struct rlimit rl;

    // First get the time limit on CPU
    getrlimit(RLIMIT_CPU, &rl);

    if (rl.rlim_max == RLIM_INFINITY || (rlim_t)_cpuLimit < rl.rlim_max) {
      // Change the time limit
      rl.rlim_cur = (rlim_t)_cpuLimit;

      // Now call setrlimit() to set the
      // changed value.
      if (setrlimit(RLIMIT_CPU, &rl) == -1) {
        std::cout << "c WARNING! Could not set resource limit: CPU-time."
                  << std::endl;
      }
    }
  }
}

void SATSolverProxy::EnableMemoryLimit() {
  if (_memoryLimit < 0) {
    return;
  }
  signal(SIGINT, SIG_Exit);
  if (_memoryLimit != 0 && _memoryLimit != INT32_MAX) {
    rlim_t new_limit = (rlim_t)_memoryLimit * 1024 * 1024;
    struct rlimit rl;
    getrlimit(RLIMIT_AS, &rl);

    if (rl.rlim_max == RLIM_INFINITY || new_limit < rl.rlim_max) {
      rl.rlim_cur = new_limit;
      if (setrlimit(RLIMIT_AS, &rl) == -1) {
        std::cout << "c WARNING! Could not set resource limit: Virtual memory."
                  << std::endl;
      }
    }
  }
}

void SATSolverProxy::SetReconf(unsigned reconf) {
  std::cout << "c " << __FUNCTION__ << " -- only possible for CMS."
            << std::endl;
}

void SATSolverProxy::SetPropagationBudget(int64_t propagationBudget) {
  std::cout << "c " << __FUNCTION__ << " -- only possible for glucose4."
            << std::endl;
}

unsigned int SATSolverProxy::SolveLimited() {
  std::cout << "c " << __FUNCTION__
            << " -- only possible for glucose4, call normal solve instead."
            << std::endl;
  return Solve();
}

void SATSolverProxy::AddHardAssumption(unsigned int *lit) {
  std::cout
      << "c " << __FUNCTION__
      << " incremental mode (hard assumptions) is not possible for this solver."
      << std::endl;
  exit(1);
}

void SATSolverProxy::ClearHardAssumption() {
  std::cout
      << "c " << __FUNCTION__
      << " incremental mode (hard assumptions) is not possible for this solver."
      << std::endl;
  exit(1);
}
