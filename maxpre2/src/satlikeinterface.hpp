#ifndef SATLIKEINTERFACE_HPP
#define SATLIKEINTERFACE_HPP

#include "probleminstance.hpp"
#include "global.hpp"
#include "log.hpp"
#include <vector>
using namespace std;

	
namespace SatlikeInterface {
	// do satlike search, return true if solution is found and if solution is found, the model is in best_solution
	bool do_search(const maxPreprocessor::ProblemInstance& pi, const vector<bool>& initial_solution, vector<bool>& best_solution, double timeLimit);
};

#endif
