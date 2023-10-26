/********************************************************************************************
Main.cpp -- Copyright (c) 2020, Tobias Paxian, Nhat Minh Hoang

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

#include <cassert>
#include <cstring>
#include <iostream>
#include "CLI11.hpp"  // Parser
#include "SATSolverProxy.h"
#include "Timer.h"

double timeParseCnfFile;

void skipWhiteSpace(char **p) {
  while (**p == ' ') {
    ++(*p);
  }
}

bool isEndOfLine(const char *p) { return *p == '\0'; }

long long int fast_atoi(char **p) {
  long long int x = 0;
  bool neg = false;
  if (**p == '-') {
    neg = true;
    ++(*p);
  }

  while (**p >= '0' && **p <= '9') {
    x = (x * 10) + (**p - '0');
    ++(*p);
  }

  return neg ? -x : x;
}

long long int getInt(char **p) {
  skipWhiteSpace(p);
  return fast_atoi(p);
}

bool match(char **p, const char *str) {
  for (; *str != '\0'; ++str, ++(*p)) {
    if (*str != **p) return false;
  }

  return true;
}

bool parseCnfFile(SATSolverProxy *_satSolver, const std::string &cnfFile) {
  FunctionTimer timer(&timeParseCnfFile, TimeUnit::MILLISECONDS);
  uint32_t nbVars = 0;
  uint32_t nbClauses;

  std::ifstream source;
  source.open(cnfFile.c_str(), std::ifstream::in);
  assert(source.good());

  std::string line;
  while (source.good() && getline(source, line)) {
    if (line.length() > 0 && line[0] == 'c') continue;

    // Copy immutable c-string to mutable array
    char line_array[line.length() + 1];
    std::strncpy(line_array, line.c_str(), line.length());
    line_array[line.length()] = '\0';
    char *cLine = line_array;

    // Read header
    if (*cLine == 'p') {
      if (!match(&cLine, "p")) {
        std::cout << "c Header must start with p" << std::endl;
        exit(1);
      }
      skipWhiteSpace(&cLine);

      if (!match(&cLine, "cnf")) {
        std::cout << "c File ist not a SAT instance" << std::endl;
        exit(1);
      }
      skipWhiteSpace(&cLine);

      nbVars = static_cast<uint32_t>(getInt(&cLine));
      nbClauses = static_cast<uint32_t>(getInt(&cLine));

      _satSolver->NewVariables(nbVars + 1);
      continue;
    }

    _satSolver->NewClause();

    // Start parsing the clauses
    int literal;

    while (!(isEndOfLine(cLine))) {
      literal = getInt(&cLine);
      if (literal == 0) break;

      assert(static_cast<uint32_t>(abs(literal)) <= nbVars);
      _satSolver->AddLiteral(&literal);
    }

    if (!_satSolver->CommitClause()) return false;

    _satSolver->ResetClause();
  }
  std::cout << "c Successfully read in CNF-File!" << std::endl;
  std::cout << "c With " << nbClauses << " clauses and " << nbVars << " variables." << std::endl;
  return true;
}

SATSolverType returnSolverType(const std::string &solver) {
  if (solver == "antom") return SATSolverType::ANTOM;
  if (solver == "glucose3") return SATSolverType::GLUCOSE3;
  if (solver == "glucose421") return SATSolverType::GLUCOSE421;
  if (solver == "cryptominisat") return SATSolverType::CRYPTOMINISAT;
  if (solver == "cadical") return SATSolverType::CADICAL;
  if (solver == "lingeling") return SATSolverType::LINGELING;

  return SATSolverType::MAPLEGLUCOSE;
}

void printSAT(const int satEnum) {
  switch (satEnum) {
    case 10:
      std::cout << "c SAT" << std::endl;
      break;

    case 20:
      std::cout << "c UNSAT" << std::endl;
      break;

    default:
      std::cout << "c UNKNOWN" << std::endl;
      break;
  }
}

TimeUnit getTimeUnit(const std::string &tunit) {
  if (tunit == "milliseconds") return TimeUnit::MILLISECONDS;
  if (tunit == "microseconds") return TimeUnit::MICROSECONDS;

  return TimeUnit::SECONDS;
}

//=================================================================================================
int main(int argc, char **argv) {
  // Option Parsing with CLI11
  CLI::App app{"\nSAT Solver Proxy\n"};

  // Define options
  std::string maxCnfFile;
  std::string solver = "cadical";
  uint32_t nThreads = 1;
  double cpuTimeLimit = 0;
  uint32_t memoryLimit = 0;
  std::string tunit = "seconds";
  int64_t propBudget = -1;

  app.add_option("-f, --file, 1", maxCnfFile, "The used max CNF file")
      ->required()
      ->check(CLI::ExistingFile);
  app.add_set("-s, --solver, 2", solver,
              {"antom", "glucose3", "glucose421", "cryptominisat", "cadical",
               "lingeling", "mapleglucose"},
              "The chosen solver.")
      ->required();
  app.add_set("-u, --unit, 3", tunit,
              {"seconds", "microseconds", "milliseconds"},
              "Set the time unit.");
  app.add_option("-t, --threads", nThreads, "The number of desired threads.",
                 true);
  app.add_option("-c, --cpuTimeLimit", cpuTimeLimit,
                 "Set the time limit for the cpu.");
  app.add_option("-m, --memoryLimit", memoryLimit,
                 "Set the maximal amount of memory the solver should use.");
  app.add_option(
      "-p,--propagationBudget", propBudget,
      "Set the maximal amount of propagations the solver should use.");

  CLI11_PARSE(app, argc, argv)

  std::cout << "c This is SATSolverProxy 2019" << std::endl;

  SATSolverProxy *satSolverProxy =
      SATSolverProxy::InitSATSolver(returnSolverType(solver), nThreads);
  satSolverProxy->SetTimeLimit(cpuTimeLimit);
  satSolverProxy->SetMemoryLimit(memoryLimit);

  int returnValue = 0;
  if (parseCnfFile(satSolverProxy, maxCnfFile)) {
    std::cout << "c Time spent on reading the cnf file: " << timeParseCnfFile
              << " "
              << "milliseconds." << std::endl;

    std::cout << "c Start solving!" << std::endl;
    CpuTimer cpuClock;
    cpuClock.SetTimeReference();
    if (propBudget != -1) {
      satSolverProxy->SetPropagationBudget(propBudget);
      returnValue = static_cast<int>(satSolverProxy->SolveLimited());
    } else {
      returnValue = static_cast<int>(satSolverProxy->Solve());
    }
    cpuClock.Stop(getTimeUnit(tunit));
    std::cout << "c Time spent on solving: " << cpuClock.TotalTime() << " "
              << tunit << "." << std::endl;

    printSAT(returnValue);
  } else {
    std::cout << "c Error in Parsing CNF File! This should generally mean:"
              << std::endl;
    std::cout << "c UNSAT" << std::endl;
  }

  return returnValue;
}
