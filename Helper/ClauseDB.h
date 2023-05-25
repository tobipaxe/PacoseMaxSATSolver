#ifndef CLAUSEDB_H
#define CLAUSEDB_H

#include <iostream>
#include <stdint.h>
#include <vector>
#include <utility>

namespace Pacose {

struct ClauseDB {

// struct clause {
//   // size_t size;
//   uint64_t weight;

//   std::vector<unsigned> literals;
// };

// std::vector<std::pair<std::vector<unsigned>*,uint64_t>> sclauses;
// std::vector<std::vector<unsigned>*> clauses;
std::vector<std::vector<int>> clauses;
std::vector<uint64_t> weights;


uint64_t sumOfSoftWeights = 0;
uint64_t maxWeight = 0;
int nbSoftClauses = 0;
int nbVarsInHard = 0;
int nbClauses = 0;
int nbVars = 0;

inline unsigned SignedToUnsignedLit(int literal) {
  return (static_cast<unsigned>(abs(literal)) << 1) ^ (literal < 0);
}

};

} // Namespace Pacose

#endif // CLAUSEDB_H