// Failed and redundant literal detection and elimination
// equivalence detection and redundant equivalence detection
//  (x->y and neg(x)->neg(y))
// forced value detection and redundant forced value detection
//  (x->y and neg(x)->y)



inline bool commonLitExists(vector<int>& a, vector<int>& b, int ignore = -1) {
	// assumes a and b are sorted
	if (a.size() < b.size()) {
		for (int l : a) {
			if (l==ignore) continue;
			if (binary_search(b.begin(), b.end(), l)) return 1;
		}
	} else {
		for (int l : b) {
			if (l==ignore) continue;
			if (binary_search(a.begin(), a.end(), l)) return 1;
		}
	}
	return 0;
}

int Preprocessor::tryFLE(int lit, vector<int>& up, bool doRLE) {
	vector<int> a={lit};
	if (!satSolver->propagate(a, up, 2)) {
		stats["FLE_FLE_success"]+=1;
		if (pi.isLabel[litVariable(lit)]) {
			rLog.removeLabel(1);
		} else{
			rLog.removeVariable(1);
		}
		a[0]=litNegation(lit);
		satSolver->addClause(a);
		int removed = setVariable(litVariable(lit), !litValue(lit));
		rLog.removeClause(removed);
		return removed;
	}
	bool flit=0;
	for (unsigned i=0; i<up.size(); ++i) {
		if (pi.isVarRemoved(litVariable(up[i]))) {
			up[i--]=up.back();
			up.pop_back();
		}
		if (up[i]==lit) flit=1;
	}
	if (!flit) up.push_back(lit);
	
	sort(up.begin(), up.end());
	
	if (doRLE && !pi.isLitLabel(lit)) {
		bool allsat = 1;
		for (int c : pi.litClauses[lit]) {
			if (!commonLitExists(up, pi.clauses[c].lit, lit)) {
				allsat=0;
				break;
			}
		}
		if (allsat) {
			stats["FLE_RLE_success"]+=1;
			if (pi.isLabel[litVariable(lit)]) {
				rLog.removeLabel(1);
			} else{
				rLog.removeVariable(1);
			}
			a[0]=litNegation(lit);
			satSolver->addClause(a);
			int removed = setVariable(litVariable(lit), !litValue(lit));
			rLog.removeClause(removed);
			return removed;
		}
	}
	return 0;
}

void Preprocessor::replaceLit(int l, int lit) {
	// replace l with lit on every clause
	// handles also labels
	vector<int> clauses = pi.litClauses[l];
	for (int c : clauses) {
		if (pi.clauses[c].weight != HARDWEIGHT) {
			// handling labels
			
			if (pi.isLitLabel(lit)) {
				// simple add weight
				pi.clauses[pi.litClauses[lit][0]].weight += pi.clauses[c].weight;
				pi.removeClause(c);
			} else if (pi.isLitLabel(litNegation(lit))) {
				// remove weight, three cases
				if (pi.clauses[c].weight < pi.labelWeight(litVariable(lit))) {
					trace.removeWeight(pi.clauses[c].weight);
					
					pi.clauses[pi.litClauses[litNegation(lit)][0]].weight -= pi.clauses[c].weight;
					pi.removeClause(c);
				} else if (pi.labelWeight(litVariable(lit)) < pi.clauses[c].weight) {
					trace.removeWeight(pi.labelWeight(litVariable(lit)));
					
					int cl = pi.litClauses[litNegation(lit)][0];
					pi.replaceLiteralInClause(litNegation(lit), lit, cl);
					swap(pi.litClauses[lit][0], pi.litClauses[lit].back());
					pi.clauses[cl].weight = pi.clauses[c].weight - pi.clauses[cl].weight;
					pi.isLabel[litVariable(lit)] = litPolarity(lit);
					pi.removeClause(c);
				} else {
					trace.removeWeight(pi.labelWeight(litVariable(lit)));
					pi.removeClause(c);
					pi.removeClause(pi.litClauses[litNegation(lit)][0]);
					pi.isLabel[litVariable(lit)] = VAR_UNDEFINED;
				}
			} else {
				// make lit a label
				pi.replaceLiteralInClause(l, lit, c);
				swap(pi.litClauses[lit][0], pi.litClauses[lit].back());
				pi.isLabel[litVariable(lit)] = litPolarity(lit);
			}
			pi.isLabel[litVariable(l)] = VAR_UNDEFINED;
			continue;
		}
		
		// handling nonlabels
		
		if (binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), litNegation(lit))) {
			pi.removeClause(c);
		} else if (binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), lit)) {
			pi.removeLiteralFromClause(l, c);
		} else {
			pi.replaceLiteralInClause(l, lit, c);
		}
	}
}

void Preprocessor::handleEqLits(vector<int>& lits) {
	// select from lits the variable that has most clauses
	int liti=0;
	int nofcl=pi.litClauses[lits[0]].size()+pi.litClauses[litNegation(lits[0])].size();
	for (unsigned i=1; i<lits.size(); ++i) {
		if (pi.litClauses[lits[i]].size() + pi.litClauses[litNegation(lits[i])].size() > (unsigned)nofcl) {
			nofcl=pi.litClauses[lits[i]].size() + pi.litClauses[litNegation(lits[i])].size();
			liti=i;
		}
	}
	
	int lit=lits[liti];
	lits[liti]=lits.back(); lits.pop_back();
	
	// replace all other variables with the selected one
	vector<int> clauses;
	for (int l : lits) {
		if (pi.isLabel[litVariable(l)]) {
			rLog.removeLabel(1);
		} else {
			rLog.removeVariable(1);
		}
		
		replaceLit(l, lit);
		replaceLit(litNegation(l), litNegation(lit));
		// on sat solver, only add l->lit and neg(l)->neg(lit)
		vector<int> cl={l, litNegation(lit)};
		satSolver->addClause(cl);
		cl[0]=litNegation(l);
		cl[1]=lit;
		satSolver->addClause(cl);

		trace.setEqual(lit, l);
		assert(pi.isVarRemoved(litVariable(l)));
	}
}


int Preprocessor::tryFLE(int var, bool doRLE, bool findEqs, bool findRedEqs, bool findForced, bool findRedForced) {
	assert(findEqs || !findRedEqs);
	assert(findForced || !findRedForced);
	
	vector<int> up1;
	vector<int> up2;
	if (tryFLE(posLit(var), up1, doRLE)) return 1;
	if (tryFLE(negLit(var), up2, doRLE)) return 1;
	vector<int> eqs;
	set<int> uvars;
	
	int removed=0;
	if (up1.size() < up2.size() && (findEqs || findForced)) {
		for (int l : up1) {
			if (findEqs && binary_search(up2.begin(), up2.end(), litNegation(l))) {
				eqs.push_back(l);
				uvars.insert(litVariable(l));
				if (l!=posLit(var)) stats["FLE_Eqs"]+=1;
			}else if (findForced && binary_search(up2.begin(), up2.end(), l)) {
				vector<int> a={l};
				satSolver->addClause(a);
				if (pi.isLabel[litVariable(l)]) rLog.removeLabel(1);
				else                            rLog.removeVariable(1);
				removed+=setVariable(litVariable(l), litValue(l));
				uvars.insert(litVariable(l));
				stats["FLE_ForcedValues"]+=1;
			}
		}
	} else if (findEqs || findForced) {
		for (int l : up2) {
			if (findEqs && binary_search(up1.begin(), up1.end(), litNegation(l))) {
				eqs.push_back(litNegation(l));
				uvars.insert(litVariable(l));
				if (l!=negLit(var)) stats["FLE_Eqs"]+=1;
			}
			if (findForced && binary_search(up1.begin(), up1.end(), l)) {
				vector<int> a={l};
				satSolver->addClause(a);
				if (pi.isLabel[litVariable(l)]) rLog.removeLabel(1);
				else                            rLog.removeVariable(1);
				removed+=setVariable(litVariable(l), litValue(l));
				uvars.insert(litVariable(l));
				stats["FLE_ForcedValues"]+=1;
			}
		}
	}
	if ((findRedEqs || findRedForced) && !pi.isLabel[var]) {
		for (int l : up1) {
			if (l==posLit(var)) continue;
			if (uvars.count(litVariable(l))) continue;
			if (pi.isLabel[litVariable(l)]) continue;
			if (findRedForced) {
				// lit -> l is known
				// check if neg(lit) -> l is redundant
				bool allsat=1;
				for (int c : pi.litClauses[litNegation(l)]) {
					if (!commonLitExists(up2, pi.clauses[c].lit, negLit(var))) {
						allsat=0;
						break;
					}
				}
				if (allsat) {
					vector<int> a={l};
					satSolver->addClause(a);
					uvars.insert(litVariable(l));
					rLog.removeVariable(1);
					removed+=setVariable(litVariable(l), litValue(l));
					stats["FLE_RedForcedValues"]+=1;
					continue;
				}
			}
			if (findRedEqs) {
				// lit -> l is known
				// check if neg(lit) -> neg(l) is redundant
				bool allsat=1;
				for (int c : pi.litClauses[l]) {
					if (!commonLitExists(up2, pi.clauses[c].lit, negLit(var))) {
						allsat=0;
						break;
					}
				}
				if (allsat) {
					eqs.push_back(l);
					uvars.insert(litVariable(l));
					stats["FLE_RedEqs"]+=1;
				}
			}
		}
		for (int l : up2) {
			if (l==negLit(var)) continue;
			if (uvars.count(litVariable(l))) continue;
			if (pi.isLabel[litVariable(l)]) continue;
			if (findRedForced) {
				// neg(lit) -> l is known
				// check if lit -> l is redundant
				bool allsat=1;
				for (int c : pi.litClauses[litNegation(l)]) {
					if (!commonLitExists(up1, pi.clauses[c].lit, posLit(var))) {
						allsat=0;
						break;
					}
				}
				if (allsat) {
					uvars.insert(litVariable(l));
					vector<int> a={l};
					satSolver->addClause(a);
					rLog.removeVariable(1);
					removed+=setVariable(litVariable(l), litValue(l));
					stats["FLE_RedForcedValues"]+=1;
					continue;
				}
			}
			if (findRedEqs) {
				// neg(lit) -> l is known
				// check if lit -> neg(l) is redundant
				bool allsat=1;
				for (int c : pi.litClauses[l]) {
					if (!commonLitExists(up1, pi.clauses[c].lit, posLit(var))) {
						allsat=0;
						break;
					}
				}
				if (allsat) {
					eqs.push_back(litNegation(l));
					uvars.insert(litVariable(l));
					stats["FLE_RedEqs"]+=1;
				}
			}
		}
	}
	rLog.removeClause(removed);
	if (findEqs) {
		bool f=0;
		for (int l : eqs) if (l==posLit(var)) f=1;
		assert(f);
		// Similar technique proposed in Probing-Based Preprocessing Techniquesfor Propositional Satisfiability, Inˆes Lynce and Jo ̃ao Marques-Silva, 2003
		if (eqs.size()>1) {
			handleEqLits(eqs);
			return 1;
		}
	}
	return removed;
}





int Preprocessor::doFLE(bool doRLE, bool findEqs, bool findRedEqs, bool findForced, bool findRedForced) {
	prepareSatSolver();
	stats["FLE_FLE_success"]+=0;	
	stats["FLE_RLE_success"]+=0;	
	stats["FLE_Eqs"]+=0;
	stats["FLE_RedEqs"]+=0;
	stats["FLE_ForcedValues"]+=0;
	stats["FLE_RedForcedValues"]+=0;
	

	
	rLog.startTechnique(Log::Technique::FLE);	
	int removed=0;
	stats["doFLE"]+=1;
	vector<int> tvars =  pi.tl.getTouchedVariables("FLE");
	sort(tvars.begin(), tvars.end());
	vector<int> svars;
	for (int vvar=0; vvar<pi.vars; ++vvar) {
		int var=vvar+flePos;
		if (var>=pi.vars) var-=pi.vars;
		if (pi.isVarRemoved(var)) continue;
		if (!binary_search(tvars.begin(), tvars.end(), var)) {
			svars.push_back(var);
			continue;
		}
		
		if (!rLog.requestTime(Log::Technique::FLE)) {
			flePos=var;
			stats["FLE_stop_position"]+=float(vvar-svars.size())/pi.vars;
			rLog.stopTechnique(Log::Technique::FLE);
			return removed;
		}
		removed += tryFLE(var, doRLE, findEqs, findRedEqs, findForced, findRedForced);
	}
	
	for (unsigned i=0; i<svars.size(); ++i) {
		int var=svars[i];
		if (pi.isVarRemoved(var)) continue;
		if (!rLog.requestTime(Log::Technique::FLE)) {
			flePos=var;
			stats["FLE_stop_position"]+=float(pi.vars-svars.size()-i)/pi.vars;
			rLog.stopTechnique(Log::Technique::FLE);
			return removed;
		}
		removed += tryFLE(var, doRLE, findEqs, findRedEqs, findForced, findRedForced);
	}
	
	stats["FLE_stop_position"]+=1;
	rLog.stopTechnique(Log::Technique::FLE);
	return removed;
}
