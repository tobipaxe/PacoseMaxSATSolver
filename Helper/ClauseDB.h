#ifndef CLAUSEDB_H
#define CLAUSEDB_H

#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>

namespace Pacose {

struct ClauseDB {

std::vector<std::vector<int>> clauses;
std::vector<uint64_t> weights;

uint64_t sumOfSoftWeights = 0;
uint64_t maxWeight = 0;
int nbSoftClauses = 0;
int nbVarsInHard = 0;
int nbHardClauses = 0;
int nbVars = 0;
int nbUnitSoftClauses = 0;

inline uint32_t SignedTouint32_tLit(int literal) {
  return (static_cast<uint32_t>(abs(literal)) << 1) ^ (literal < 0);
}

};

} // Namespace Pacose

#endif // CLAUSEDB_H