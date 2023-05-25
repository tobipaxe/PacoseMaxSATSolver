#include "satlikeinterface.hpp"
#include "satlike/Alg_SATLike.h"
#include "global.hpp"


bool SatlikeInterface::do_search(const maxPreprocessor::ProblemInstance& pi, const vector<bool>& initial_solution, vector<bool>& best_solution, double timeLimit) {
	maxPreprocessor::Timer timer;
	timer.start();
	Satlike satlike;
	
	// init satlike
	lit** clauses;
	int* lit_count;
	long long* weights;
	clauses = new lit*[pi.clauses.size()];
	lit_count = new int[pi.clauses.size()];
	weights = new long long[pi.clauses.size()];
	unsigned clauses_size=0;
	for (unsigned ii=0, i=0; ii<pi.clauses.size(); ++ii, ++i) {
		if (pi.isClauseRemoved(ii)) {--i; continue;}
		unsigned n = pi.clauses[ii].lit.size();
		if (n==0) {--i; continue;} // empty clauses can cause satlike to crash
		
		clauses[i] = new lit[n+1];
		lit_count[i] = n;
		weights[i] = pi.clauses[ii].weight;
		for (unsigned j=0; j<n; ++j) {
			clauses[i][j].clause_num = i;
			clauses[i][j].var_num    = maxPreprocessor::litVariable(pi.clauses[ii].lit[j])+1;
			clauses[i][j].sense      = maxPreprocessor::litValue(pi.clauses[ii].lit[j]);
		}
		clauses[i][n].var_num = 0;
		clauses[i][n].clause_num = -1;
		++clauses_size;
	}
	satlike.build_instance(pi.vars, clauses_size, maxPreprocessor::HARDWEIGHT, clauses, lit_count, weights);
	
	


	vector<int> init_solution;	
	if (initial_solution.size()) {
		init_solution.resize(pi.vars+1);
	}
	for (unsigned i=0; i<initial_solution.size();++i) init_solution[i+1]=initial_solution[i];
	
	
	// do search
	satlike.settings();
	satlike.init(init_solution);
	bool found_solution = 0;
	for (int step = 1; step < satlike.max_flips; ++step) {
		if (satlike.hard_unsat_nb == 0 && (satlike.soft_unsat_weight < satlike.opt_unsat_weight || satlike.best_soln_feasible == 0)) {
			satlike.max_flips = step + satlike.max_non_improve_flip;
			if (satlike.soft_unsat_weight < satlike.opt_unsat_weight) {
				satlike.best_soln_feasible = 1;
				satlike.opt_unsat_weight = satlike.soft_unsat_weight;
				best_solution.resize(satlike.num_vars);
				for (int v = 1; v <= satlike.num_vars; ++v) best_solution[v-1] = satlike.cur_soln[v];
				found_solution = 1;
// 				cerr << "o " << satlike.opt_unsat_weight << endl;
			}
			if (satlike.opt_unsat_weight == 0)
				break;
		}
		int flipvar = satlike.pick_var();
		satlike.flip(flipvar);
		satlike.time_stamp[flipvar] = step;
		
		if (timer.getTime().count()>timeLimit) {
			break;
		}
	}
// 	cout << "c satlike search done!" << endl;
	//satlike.local_search(init_solution);
	
	satlike.free_memory();
	
	return found_solution;
}
