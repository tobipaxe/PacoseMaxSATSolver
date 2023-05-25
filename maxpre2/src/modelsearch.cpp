pair<uint64_t, uint64_t> Preprocessor::modelCostCheck(vector<bool>& model) {
	// calculates the cost of a model
	// updates the upper bound for cost if model is cheaper than current best model
	// returns a pair: first element is the sum of the unsat labels (used in hardening)
	// second element is the actual total cost of the model
	uint64_t cost1=0;
	uint64_t cost2=trace.removedWeight;;
	for (unsigned c=0; c<pi.clauses.size(); ++c) {
		if (pi.isClauseRemoved(c) || pi.clauses[c].isHard()) continue;
		bool sat=0;
		for (int lit : pi.clauses[c].lit) {
			if (model[litVariable(lit)] == litValue(lit)) {
				sat=1;
				break;
			}
		}
		if (!sat) {
			if (pi.clauses[c].lit.size())   cost1+=pi.clauses[c].weight;
			else                            cost2+=pi.clauses[c].weight;
		}
	}
	cost2+=cost1;
	if (cost2 < bestCost) {
		// best overall found model
		bestCost = cost2;
		bestModel = model;
	}
	return {cost1, cost2};
}

int Preprocessor::updateBestModel(vector<bool>& best_model, uint64_t& best_cost, vector<bool>& model) {
	pair<uint64_t, uint64_t> cost = modelCostCheck(model);
	if (cost.F<best_cost) {
		best_cost = cost.F;
		swap(best_model, model);
		return 1;
	}
	return 0;
}


void Preprocessor::findGoodModelSub(vector<pair<uint64_t, int> >& labels, int n, vector<int>& idx, vector<int> base_assumps, uint64_t& best_cost, vector<bool>& best_model, int iters_left, double leaveTime, bool binarySearch) {
	if (n==0) return;
	
	vector<bool> model;
	
// 	log("findGoodModel1 call, n=", n, " idx.size()=", idx.size());
	// binary search for a point where heaviest labels are not satisfiable
	// then extract cores until instance is satisfiable, then continue iterating the process
	int a;

	if (binarySearch) {
		a=0;
		int b=idx.size()-1;
		int nba=base_assumps.size();
		while (a<b) {
			if (!rLog.requestTime(rLog.activeTechnique)) return;
			double timeLimit = rLog.allocatedTimeLeft(rLog.activeTechnique)+leaveTime;
			if (timeLimit<0 && best_cost<(~0ull)) return;
			
			int m=(a+b+1)/2;
			for (int i=idx[m]; i<n; ++i) {
				base_assumps.push_back(labels[i].S);
			}
			stats["asearch_sat_solver_calls"]+=1;
			int rv = satSolver->solveLimited(base_assumps, timeLimit);
			
			base_assumps.resize(nba);
			
			if (rv==-1) continue; // TODO: is this bad?
			if (rv==1) {
				satSolver->getModel(model);
				updateBestModel(best_model, best_cost, model);
				b=m-1;
				if (b==0) return;
			} else {
				a=m;
			}
		}
		
		vector<int> cl={0};
		for (int i=idx[a]; i<n; ++i) {
			if (opt.hardenInModelSearch && a!=(int)idx.size()-1 && i>=idx[a+1]) {
				cl[0]=labels[i].S;
				satSolver->addClause(cl);
				// TODO: sat solver becomes not reusable after this!
				// currently doesn't matter since sat solver is created every time again before new sat solver technique
				// but this is not good anyways since satSolver is actually planned to be reusable at least in some cases...
			} else base_assumps.push_back(labels[i].S);
		}
	} else {
		// linear search
		int sat=1;
		vector<int> cl={0};
		a=idx.size();
		while (sat && --a>=0) {
			if (!rLog.requestTime(rLog.activeTechnique)) return;
			double timeLimit = rLog.allocatedTimeLeft(rLog.activeTechnique)+leaveTime;
			if (timeLimit<0 && best_cost<(~0ull)) return;
			
			for (int i=idx[a]; i<n; ++i) {
				base_assumps.push_back(labels[i].S);
			}
			stats["asearch_sat_solver_calls"]+=1;
			sat = satSolver->solveLimited(base_assumps, timeLimit);
			if (sat == -1) continue;
			if (sat) {
				satSolver->getModel(model);
				updateBestModel(best_model, best_cost, model);
				if (opt.hardenInModelSearch) {
					for (int i=idx[a]; i<n; ++i) {
						cl[0]=labels[i].S;
						satSolver->addClause(cl);
					}
				}
			}
			n=idx[a];
		}
		if (a<0) return;
	}
	int sat_partitions=idx.size()-a;
	
	if (sat_partitions>20 || ((1<<sat_partitions) > (int)idx.size())) {
		stats["sat>log"]+=1;
		if (opt.modelSearchAsearchType == -1) binarySearch = 1;
	} else {
		stats["sat<log"]+=1;		
		if (opt.modelSearchAsearchType == -1) binarySearch=0;
	}

	if (!rLog.requestTime(rLog.activeTechnique)) return;
	double timeLimit = rLog.allocatedTimeLeft(rLog.activeTechnique) - leaveTime;
	if (timeLimit<0 && best_cost<(~0ull)) return;

	// remove cores until satisfiable
	int rv=0;
	vector<int> core;
	stats["asearch_sat_solver_calls"]+=1;
	while ((rv = satSolver->solveLimited(base_assumps, timeLimit)) == 0) {
		satSolver->getCore(core);
		sort(core.begin(), core.end());
		for (unsigned i=0;i<base_assumps.size();++i) {
			if (binary_search(core.begin(), core.end(), litNegation(base_assumps[i]))) {
				base_assumps[i--]=base_assumps.back();
				base_assumps.pop_back();
			}
		}
		if (!rLog.requestTime(rLog.activeTechnique)) return;
		timeLimit = rLog.allocatedTimeLeft(rLog.activeTechnique)-leaveTime;
		if (timeLimit<0 && best_cost<(~0ull)) return;
		stats["asearch_sat_solver_calls"]+=1;
	}
	
	if (rv!=-1) {
		satSolver->getModel(model);
		updateBestModel(best_model, best_cost, model);
		if (iters_left!=1) {
			int nn=idx[a];
			idx.resize(a);
			findGoodModelSub(labels, nn, idx, base_assumps, best_cost, best_model, iters_left-1, leaveTime, binarySearch);
		}
	}
}

uint64_t Preprocessor::findGoodModelA(vector<pair<uint64_t, int> >& labels, vector<bool>& best_model, int iterLimit, double leaveTime) {
	// assumes labels are sorted
	// first try to find a good model
	int p = sqrt(labels.size())*3;
	if (!p) p=1;
	
	vector<int> idx; idx.push_back(0);
	uint64_t sum=0;
	int sz=0;
	for (unsigned i=0; i<labels.size(); ++i) {
		++sz;
		if ((i>1 && labels[i].F>sum) || (sz>=p && labels[i].F != labels[i-1].F)) {
			idx.push_back(i);
			i+=p;
			sz=0;
		}
		sum+=labels[i].F;
	}
	vector<int> base_assumps;
	uint64_t best_cost = ~0ull;
	if (iterLimit) findGoodModelSub(labels, labels.size(), idx, base_assumps, best_cost, best_model, iterLimit, leaveTime, opt.modelSearchAsearchType);
	return best_cost;
}


uint64_t Preprocessor::findGoodModel(vector<bool>& best_model, int assumpsSearchIterLimit, double improveTimeLimit, int satLikeTries, double timeSatLike) {
	prepareSatSolver();
	
	vector<pair<uint64_t, int> > labels;
	for (int lit=0; lit<2*pi.vars; ++lit) {
		if (pi.isVarRemoved(litVariable(lit))) continue;
		if (pi.isLitLabel(lit)) {
			labels.push_back({pi.labelWeight(litVariable(lit)), lit});
		}
	}
	sort(labels.begin(), labels.end());
	
	
	// try to find good model with sat solver+assumptions technique
	Timer timer;
	timer.start();
	double t=timer.getTime().count();
	
	int best_solution_from=-1;
	vector<bool> model;	
	uint64_t best_cost= ~(0ull);
 	if (assumpsSearchIterLimit) {
		double leaveTime = min(rLog.allocatedTimeLeft(rLog.activeTechnique)/2.0, improveTimeLimit);
		stats["asearch_sat_solver_calls"]+=0;
		stats["sat<log"]+=0;
		stats["sat>log"]+=0;
		best_cost=findGoodModelA(labels, best_model, assumpsSearchIterLimit, leaveTime);
		stats["best_model_asearch_time"]+=(timer.getTime().count()-t);
		t=timer.getTime().count();
		if (best_model.size()) {
			best_solution_from=0;
			if (rLog.requestTime(rLog.activeTechnique) && improveTimeLimit>0) { // try to improve the solution
				if (SatlikeInterface::do_search(pi, best_model, model, min(rLog.allocatedTimeLeft(rLog.activeTechnique), improveTimeLimit))) {
					if (updateBestModel(best_model, best_cost, model)) {
						best_solution_from=1;
					}
				}
				stats["best_model_improve_time"]+=(timer.getTime().count()-t);
				t=timer.getTime().count();
			}
		}
	}
	// try random satlike models
	// sometimes they give better results
	// okay, not almost ever they give better results..., 
	for (int i=0;i<satLikeTries;++i) {
		if (!rLog.requestTime(rLog.activeTechnique)) break;
		vector<bool> tmp;
		if (SatlikeInterface::do_search(pi, tmp, model, min(rLog.allocatedTimeLeft(rLog.activeTechnique), timeSatLike))) {
			uint64_t cost = 0;
			for (auto& l : labels) {
				if (model[litVariable(l.S)] != litValue(l.S)) {
					cost += l.F;
				}
			}
			if (cost<best_cost) {
				best_cost=cost;
				swap(best_model, model);
				best_solution_from=2;
			}
		}
		stats["best_model_satlike_time"]+=(timer.getTime().count()-t);
		t=timer.getTime().count();
	}
	
	stats["best_model_from_asearch"]+=0;
	stats["best_model_from_improve"]+=0;
	stats["best_model_from_satlike"]+=0;
	stats["no_best_model"]+=0;
	if (best_solution_from==0) stats["best_model_from_asearch"]+=1;
	else if (best_solution_from==1) stats["best_model_from_improve"]+=1;
	else if (best_solution_from==2) stats["best_model_from_satlike"]+=1;
	else stats["no_best_model"]+=1;
		
	return best_cost;
}

