// TrimMaxSAT Technique
// 


int getM(int n, double z) { // assumes 0 <= z < 1
	if (n<=2) return 1; 
	double p=sqrt(n);
	int l = ((int)(2*p)-1); // size of candidate list
	int i = ((int)(z*l))+1; // index (1-indexed) in the candidate list
	
	if (i<p) return n/i;
	else     return ((int)(2*p))-i;
}

int Preprocessor::tryTMS(vector<int>& neverSat) {
	if (!neverSat.size()) return 0;
	for (unsigned i=0; i<neverSat.size(); ++i) {
		if (canSatLits.count(neverSat[i])) {
			neverSat[i--]=neverSat.back();
			neverSat.pop_back();
		}
	}
	
	vector<bool> model;
	vector<int> assumptions;
	vector<vector<int> > sets;

	//statistics...
	int SATcalls=0;
	int UNSATcalls=0;
	
	double lb=0, z=0.5;
	while (neverSat.size()) {
		if (!rLog.requestTime(rLog.activeTechnique)) {
			log("TMS stopped due to timelimit, ", SATcalls+UNSATcalls, " SAT-solver calls done by TMS, SAT: ", SATcalls, ", UNSAT: ", UNSATcalls);
			stats["TMS_interrupted"]+=1;
			return -1;
		}
		double timeLimit = rLog.allocatedTimeLeft(rLog.activeTechnique);
		
		// create set partition
		int m = getM(neverSat.size(), z);
		
		sets.resize(m);
		for (int i=0; i<m; ++i) sets[i].clear();
		for (int i=0, j=0; i<(int)neverSat.size(); ++i) {
			sets[j].push_back(neverSat[i]);
			if (++j == m) j=0;
		}
		
		// prepare sat solver
		assumptions.resize(m);
		for (int i=0; i<m; ++i) {
			int nv = pi.addVar();
			sets[i].push_back(posLit(nv));
			satSolver->addClause(sets[i]);
			assumptions[i] = negLit(nv);
		}
		int SAT = satSolver->solveLimited(assumptions, timeLimit);
        if (SAT==-1) return 0;
		
		if (SAT) {
			++SATcalls;
			// get model, remove satisfiable clauses from neverSat
			satSolver->getModel(model);
			modelCostCheck(model); // updates UB for cost

			for (unsigned i=0; i<model.size(); ++i) {
				if (model[i]) canSatLits.insert(posLit(i));
				else          canSatLits.insert(negLit(i));
			}
			for (unsigned i=0; i<neverSat.size(); ++i) {
				if (model[litVariable(neverSat[i])] == litValue(neverSat[i])) {
					neverSat[i--]=neverSat.back();
					neverSat.pop_back();
				}
			}
			z -= (z-lb)/2;
		} else {
			++UNSATcalls;
			if (m>1) { // decrease m
				z += (1-z)/2;
			} else {   // done, all unsat softs can be hardened
				// remove added clauses...
				vector<int> cl={0};
				for (int i=0; i<m; ++i) {
					cl[0]=sets[i].back();
					satSolver->addClause(cl);
				}
				break;
			}
		}
		
		// remove added clauses...
		vector<int> cl={0};
		for (int i=0; i<m; ++i) {
			cl[0]=sets[i].back();
			satSolver->addClause(cl);
		}
	}
	
	uint64_t w=0;
	for (unsigned i=0; i<neverSat.size(); ++i){
		if (pi.isLabel[litVariable(neverSat[i])]) {
			w += pi.labelWeight(litVariable(neverSat[i]));
			rLog.removeLabel(1);
		} else {
			rLog.removeVariable(1);
			stats["TMS_Vars_removed"]+=1;
		}
		setVariable(litVariable(neverSat[i]), !litValue(neverSat[i]));
	}
	
	log(neverSat.size(), " unsatisfiable softs detected and deleted by TMS, total weight: ", w);
	log(SATcalls+UNSATcalls, " SAT-solver calls done by TMS, SAT: ", SATcalls, ", UNSAT: ", UNSATcalls);
	return neverSat.size();
}




// Use TrimMaxSAT algorithm for general backbone computing
int Preprocessor::doBBTMS(int maxVars) {
	rLog.startTechnique(Log::Technique::BBTMS);
	if (!rLog.requestTime(Log::Technique::BBTMS)) {
		rLog.stopTechnique(Log::Technique::BBTMS);
		return 0;
	}
	prepareSatSolver();
	stats["doBBTMS"]+=1;
	stats["TMS_Vars_removed"]+=0;
		
	vector<pair<int, int> > vars;
	for (int i=0; i<pi.vars; ++i) {
		if (pi.isVarRemoved(i)) continue;
		if (canSatLits.count(posLit(i)) && canSatLits.count(negLit(i))) continue;
		vars.push_back({-((int)pi.litClauses[posLit(i)].size()+pi.litClauses[negLit(i)].size()), i});
	}
	sort(vars.begin(), vars.end());
	
	vector<int> litsPos;
	vector<int> litsNeg;
	for (unsigned i=0; i<vars.size(); ++i) {
		if (maxVars>=0 && ((int)litsPos.size()>=maxVars && (int)litsNeg.size()>=maxVars)) break;
		if ((maxVars<0 || (int)litsPos.size()<maxVars) && !canSatLits.count(posLit(vars[i].S))) {
			litsPos.push_back(posLit(vars[i].S));
		}
		if ((maxVars<0 || (int)litsNeg.size()<maxVars) && !canSatLits.count(negLit(vars[i].S))) {
			litsNeg.push_back(negLit(vars[i].S));
		}
	}
	
	int removedVars=0;
	if (rLog.requestTime(Log::Technique::BBTMS)) removedVars+=tryTMS(litsPos);
	if (rLog.requestTime(Log::Technique::BBTMS)) removedVars+=tryTMS(litsNeg);
	
	rLog.stopTechnique(Log::Technique::BBTMS);
	return removedVars;
}


int Preprocessor::doTMS() {
	rLog.startTechnique(Log::Technique::TMS);
	if (!rLog.requestTime(Log::Technique::TMS)) {
		rLog.stopTechnique(Log::Technique::TMS);
		return 0;
	}
	
	stats["doTMS"]+=1;
	stats["TMS_interrupted"]+=0;
	prepareSatSolver();
	
	vector<int> bvars;
	for (int i=0; i<pi.vars; ++i) {
		if (!pi.isLabel[i] || pi.isVarRemoved(i)) continue;
		if (pi.isLabel[i] == VAR_TRUE) {
			if (!canSatLits.count(posLit(i))) bvars.push_back(posLit(i));
		} else if (pi.isLabel[i] == VAR_FALSE) {
			if (!canSatLits.count(negLit(i))) bvars.push_back(negLit(i));
		}
	}
	
	int removedSofts=0;
	if (rLog.requestTime(Log::Technique::TMS)) removedSofts=tryTMS(bvars);
	
	rLog.stopTechnique(Log::Technique::TMS);
	return removedSofts;
}
