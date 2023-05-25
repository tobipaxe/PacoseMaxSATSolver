// Detecting redundancies



bool Preprocessor::checkPositiveReduct(const vector<int>& clause, const vector<int>& labels) { // assumes literals are sorted...
	// returns true if clause doesn't affect the satisfiability of the formula
	// and any satisfying assignment of the formula \ clause can be turned to satisfy formula + clause by flipping the literals in the clause
	// which means that if the clause doesn't contain any labels, it can be either added or removed from the formula preserving the cost
	
	// vector labels contains ''extra'' labels that could be added to clause but which are allowed only to be in one polarity
	// in some cases this may reduce the number of clauses added to the solver possibly making the formula easier and even change result from UNSAT to SAT
	
	
	// TODO: don't always initialize a new sat solver...
	SATSolver* solver = (SATSolver*)new Glucose3();
	
	// find clauses satisfied by negations of literals in clause
	vector<int> cls;
	for (int l : clause) {
		for (int c : pi.litClauses[litNegation(l)]) {
			if (pi.isClauseRemoved(c)) continue;
			if (!pi.clauses[c].isHard()) continue;
			cls.push_back(c);
		}
	}
	sort(cls.begin(), cls.end());
	cls.erase(unique(cls.begin(), cls.end()), cls.end());
	
	// check if some clauses can be removed
	vector<int> assumptions;
	for (int l : labels) {
		if (binary_search(clause.begin(), clause.end(), l)) continue;
		bool ok=1;
		for (int c : pi.litClauses[litNegation(l)]) {
			if (pi.isClauseRemoved(c)) continue;
			if (!pi.clauses[c].isHard()) continue;
			if (!binary_search(cls.begin(), cls.end(), c) && pi.isLitLabel(litNegation(l))) {
				ok=0;
				break;
			}
		}
		if (ok) assumptions.push_back(pi.isLitLabel(l) ? l : litNegation(l));
	}
	sort(assumptions.begin(), assumptions.end());
	
	// add touched versions of clauses
	for (int c : cls) {
		vector<int> tcl;
		bool skip = 0;
		for (int l : pi.clauses[c].lit) {
			if (binary_search(clause.begin(), clause.end(), litNegation(l)) || binary_search(clause.begin(), clause.end(), l)) {
				tcl.push_back(l);
			}
			if (binary_search(assumptions.begin(), assumptions.end(), l)) {
				skip=1;
				break;
			}
		}
		if (!skip) solver->addClause(tcl);
	}
	
	solver->addClause(clause);
	
	int sat = solver->solveLimited(rLog.allocatedTimeLeft(rLog.activeTechnique));
	delete solver;
	++redSatSolverCalls;
	stats["check_positive_reduct_size"]+=clause.size();
	stats["check_positive_reduct_calls"]+=1;
	if (sat==-1) {
		stats["check_filtered_positive_reduct_timeout"]+=1;
		stats["check_filtered_positive_reduct_timeout_size"]+=clause.size();
		return false;
	} else if (sat) {
		stats["check_positive_reduct_sat"]+=1;
		stats["check_positive_reduct_sat_size"]+=clause.size();
	} else {
		stats["check_positive_reduct_unsat"]+=1;
		stats["check_positive_reduct_unsat_size"]+=clause.size();
	}
	return sat;
}


bool Preprocessor::checkExtendedPositiveReduct(const vector<int>& clause_, const vector<int>& labels) { // assumes literals are sorted...
	// returns true if clause doesn't affect the satisfiability of the formula
	// and any satisfying assignment of the formula \ clause can be turned to satisfy formula + clause by flipping the literals in the clause
	// which means that if the clause doesn't contain any labels, it can be either added or removed from the formula preserving the cost
	
	// vector labels contains ''extra'' labels that could be added to clause but which are allowed only to be in one polarity
	// in some cases this may reduce the number of clauses added to the solver possibly making the formula easier and even change result from UNSAT to SAT
	
	// Extend clause by UP
	vector<int> alpha;
	alpha.resize(clause_.size());
	for (unsigned i=0; i<clause_.size();++i) {
		alpha[i]=litNegation(clause_[i]);
	}
	vector<int> up;
	bool s = satSolver->propagate(alpha, up, 2);
	if (!s) {
		return true;
	}
	vector<int> clause;
	for (int l : clause_) clause.push_back(l);
	for (unsigned i=0; i<up.size(); ++i) {
		if (!pi.isLabel[litVariable(up[i])] && !binary_search(clause_.begin(), clause_.end(), up[i])) {
			clause.push_back(litNegation(up[i]));
		}
	}
	sort(clause.begin(), clause.end());
	return checkFilteredPositiveReduct(clause, labels);
}


bool Preprocessor::checkFilteredPositiveReduct(vector<int>& clause, const vector<int>& labels, bool f2) { // assumes literals are sorted...
	// returns true if clause doesn't affect the satisfiability of the formula
	// and any satisfying assignment of the formula \ clause can be turned to satisfy formula + clause by flipping the literals in the clause
	// which means that if the clause doesn't contain any labels, it can be either added or removed from the formula preserving the cost
	
	// vector labels contains ''extra'' labels that could be added to clause but which are allowed only to be in one polarity
	// in some cases this may reduce the number of clauses added to the solver possibly making the formula easier and even change result from UNSAT to SAT
	
	// TODO: don't always initialize a new sat solver...
	SATSolver* solver = (SATSolver*)new Glucose3();
	
	// find clauses satisfied by negations of literals in clause
	vector<int> cls;
	for (int l : clause) {
		for (int c : pi.litClauses[litNegation(l)]) {
			if (pi.isClauseRemoved(c)) continue;
			if (!pi.clauses[c].isHard()) continue;
			cls.push_back(c);
		}
	}
	sort(cls.begin(), cls.end());
	cls.erase(unique(cls.begin(), cls.end()), cls.end());
	
	vector<int> alpha;
	alpha.resize(clause.size());
	for (unsigned i=0; i<clause.size();++i) {
		alpha[i]=litNegation(clause[i]);
	}
	vector<int> up;
	satSolver->propagate(alpha, up, 2);
	sort(up.begin(), up.end());
	
	if (f2) {
		for (int l : up) {
			if (binary_search(clause.begin(), clause.end(), litNegation(l))) continue;
			bool bad = false;
			for (int c : pi.litClauses[l]) {
				if (!binary_search(cls.begin(), cls.end(), c)) {
					bad = true;
					break;
				}
			}
			if (!bad) {
				clause.push_back(l);
			}
		}
	}
	
	// check if some clauses can be removed
	vector<int> assumptions;
	for (int l : labels) {
		if (binary_search(clause.begin(), clause.end(), l)) continue;
		bool ok=1;
		for (int c : pi.litClauses[litNegation(l)]) {
			if (pi.isClauseRemoved(c)) continue;
			if (!pi.clauses[c].isHard()) continue;
			if (!binary_search(cls.begin(), cls.end(), c) && pi.isLitLabel(litNegation(l))) {
				ok=0;
				break;
			}
		}
		if (ok) assumptions.push_back(pi.isLitLabel(l) ? l : litNegation(l));
	}
	sort(assumptions.begin(), assumptions.end());
	
	// add touched versions of clauses
	for (int c : cls) {
		vector<int> tcl;
		bool skip = 0;
		for (int l : pi.clauses[c].lit) {
			if (binary_search(clause.begin(), clause.end(), litNegation(l)) || binary_search(clause.begin(), clause.end(), l)) {
				tcl.push_back(l);
			} else if (binary_search(up.begin(), up.end(), l)) {
				skip=1;
				break;
			}
			if (binary_search(assumptions.begin(), assumptions.end(), l)) {
				skip=1;
				break;
			}
		}
		if (!skip) solver->addClause(tcl);
	}
	
	solver->addClause(clause);
	
	int sat = solver->solveLimited(rLog.allocatedTimeLeft(rLog.activeTechnique));
	delete solver;
	++redSatSolverCalls;
	stats["check_filtered_positive_reduct_size"]+=clause.size();
	stats["check_filtered_positive_reduct_calls"]+=1;
	if (sat==-1) {
		stats["check_filtered_positive_reduct_timeout"]+=1;
		stats["check_filtered_positive_reduct_timeout_size"]+=clause.size();
		return false;
	} else if (sat) {
		stats["check_filtered_positive_reduct_sat"]+=1;
		stats["check_filtered_positive_reduct_sat_size"]+=clause.size();
	} else {
		stats["check_filtered_positive_reduct_unsat"]+=1;
		stats["check_filtered_positive_reduct_unsat_size"]+=clause.size();
	}
	return sat;
}


bool Preprocessor::checkTrimmedFilteredPositiveReduct(const vector<int>& clause_, const vector<int>& labels, bool f2) { // assumes literals are sorted...
	vector<int> clause(clause_);
	int removed = trimReductClause(clause);
	rLog.removeClause(removed);
	return checkFilteredPositiveReduct(clause, labels, f2);
}


int Preprocessor::trimReductClause(vector<int>& clause) {
	// check if negations of some literals in clause unit propagate some other negations, in that case the propagated can be removed
	// returns number of clauses removed accidentally, if UP call returns UNSAT
	set<int> toRemove;
	sort(clause.begin(), clause.end());
	for (int l : clause) {
		if (toRemove.count(l)) continue;
		vector<int> a={litNegation(l)};
		vector<int> c;
		bool s = satSolver->propagate(a, c, 2);
		if (s) {
			for (int lit : c) {
				if (litNegation(lit)!=l && binary_search(clause.begin(), clause.end(), litNegation(lit))) {
					toRemove.insert(litNegation(lit));
				}
			}
		} else{
			clause.clear();
			vector<int> cl={l};
			satSolver->addClause(cl);
			return setVariable(litVariable(l), litValue(l));
		}
	}
	
	int jj=0;
	for (unsigned j=0;j<clause.size();++j) {
		if (!toRemove.count(clause[j])) {
			clause[jj++]=clause[j];
		}
	}
	clause.resize(jj);
	
	return 0;
}


unsigned hasMoreInCommon(vector<int>& c1, vector<int>& c2, unsigned lim) {
	// check if clause 1 and clause 2 have more or less lits in common than lim
	// if d=0, check if more, if d=1, check if less.
	if (min(c1.size(), c2.size()) <= lim) {
		return 0;
	}
	unsigned cm=0;
	unsigned i=0, j=0;
	while (i<c1.size() && j<c2.size()) {
		if (c1[i]==c2[j]) {
			++i; ++j; ++cm;
		} else {
			if (c1[i]<c2[j]) {
				++i;
			} else {
				++j;
			}
			if (cm + min(c1.size()-i, c2.size()-j)<= lim) {
				return 0;
			}
		}
	}
	return cm;
}

int Preprocessor::findREDPartitionForLit(int lit, unsigned n, vector<pair<vector<pair<uint64_t, int > >, vector<int> > >& clauses, uint64_t mcost) {
	// tries to create n clauses that subsume all clauses where negLit(lit) is
	// uses greedy algorithm that doesn't necessarily find a solution even if one exists...
	// if succeeded, on clauses there will be n pairs: <clause, labels>, clause containt the literals, assumptions contains labels that may be useful when checking redundancy
	// clauses may contain labels, that have at most cost mcost (should be labelWeight(litVariable[lit]), if lit is label, 0 otherwise).
	if (pi.litClauses[litNegation(lit)].size() < n) return 0;
	
	
	vector<int> cl(pi.litClauses[litNegation(lit)]);
	
	
	clauses.resize(n);
	unsigned theoretically_best = 0;
	uint64_t cost=0;
	for (unsigned i=0; i<n; ++i) {
		if (i>0 && n<pi.litClauses[litNegation(lit)].size()) { // choose base clause heuristically: select the clause that has least maximum literals in common with any earlier base clause
			int j=i;
			unsigned bstd=1e9;
			for (unsigned jj=i; jj<cl.size(); ++jj) {
				unsigned maxd=0;
				for (unsigned ii=0;ii<i;++ii) {
					maxd=max(maxd, hasMoreInCommon(pi.clauses[cl[ii]].lit, pi.clauses[cl[jj]].lit, maxd));
					if (maxd>=bstd) break;
				}
				if (maxd<bstd) {
					bstd=maxd;
					j=jj;
					if (bstd == theoretically_best) break;
				}
			}
			theoretically_best = bstd;
			swap(cl[i], cl[j]);
		}
		
		int c=cl[i];
		
		vector<pair<uint64_t, int> >& clause = clauses[i].F;
		vector<int>& labels = clauses[i].S;
		uint64_t a=HARDWEIGHT, b=HARDWEIGHT;
		for (int l : pi.clauses[c].lit) {
			if (l==litNegation(lit)) continue;
			if (pi.isLabel[litVariable(l)]) {
				if (pi.labelWeight(litVariable(l)) <= mcost) {
					uint64_t w=pi.labelWeight(litVariable(l));
					clause.push_back({w, l});
					if (w<b) {
						b=w;
						if (a>b) swap(a, b);
					}
				}
				labels.push_back(l);
				continue;
			}
			clause.push_back({0, l});
			b=a;
			a=0;
		}
		cost += a+b;
		if (cost>mcost) return 0;
	}
	
	for (unsigned ii=n; ii<cl.size(); ++ii) {
		int i=0; 
		int bstd=0;
		for (unsigned jj=0;jj<n;++jj) {
			int d=hasMoreInCommon(pi.clauses[cl[ii]].lit, pi.clauses[cl[jj]].lit, bstd);
			if (d>bstd) {
				bstd=d;
				i=jj;
			}
		}
		vector<pair<uint64_t, int> >& clause = clauses[i].F;
		vector<int>& labels = clauses[i].S;
		
		int c = cl[ii];
		
		int jj=0;
		uint64_t ao=HARDWEIGHT, bo=HARDWEIGHT;
		uint64_t a=HARDWEIGHT, b=HARDWEIGHT;
		for (unsigned j=0;j<clause.size();++j) {
			if (binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), clause[j].S)) {
				clause[jj++]=clause[j];
				if (clause[j].F<b) {
					b=clause[j].F;
					if (a>b) swap(a, b);
					if (clause[j].F<bo) {
						bo=clause[j].F;
						if (ao>bo) swap(ao, bo);
					}
				}
			} else {
				if (clause[j].F<bo) {
					bo=clause[j].F;
					if (ao>bo) swap(ao, bo);
				}
			}
		}
		if (jj < 2) return 0;
		clause.resize(jj);
		
		cost += a+b-ao-bo;
		if (cost>mcost) return 0;
		
		jj=0;
		for (unsigned j=0;j<labels.size();++j) {
			if (binary_search(pi.clauses[c].lit.begin(), pi.clauses[c].lit.end(), labels[j])) {
				labels[jj++]=labels[j];
			}
		}
		labels.resize(jj);
	}
	
	// sort clauses in ascending order to the size to reduce the number of sat solver calls, short clauses first, they are most likely to not be redundant
	auto cmp = [&](const pair<vector<pair<uint64_t, int> >, vector<int> >& a, const pair<vector<pair<uint64_t, int> >, vector<int> >& b) {
		return a.F.size() < b.F.size();
	};
	
	sort(clauses.begin(), clauses.end(), cmp);
	
	return cost>0?1:2;
}

int Preprocessor::tryREDOnLit(int lit) {
	// find literals that are common in all clauses where neg(lb) is
	// and check if that clause is redundant, meaning it can be added and since it subsumes all clauses where neg(lb) is,
	// lb can then be hardened
	// case 1: no labels in redundant clause: harden lb, add reduct, subsume clauses
	// case 2: labels in redundant clause that have total weight less than the weight of label, harden lb
	
	uint64_t mcost = 0;
	if (pi.isLitLabel(lit)) {
		mcost = pi.labelWeight(litVariable(lit));
	} else if (pi.isLitLabel(litNegation(lit))) {
		return 0;
	}
	
	if (!pi.litClauses[litNegation(lit)].size()) return 0;	
	
	int m=pi.litClauses[litNegation(lit)].size();
	// try to find n clauses that would subsume clauses where litNegation(lb) is
	// if found, we can harden
	int UNSATChecksLeft = opt.LRED_maxUNSATReductChecksPerLabel;
	if (UNSATChecksLeft==0) return 0;
	
	for (int n=min(m, opt.LRED_minPartitions); n<=m; n = ((n < m && (n<<1) >= m) ? m : (n<<1) ) ) {
		// log("try RED partition with ", n, " parts");
		// try without allowing labels...
		vector<pair<vector<pair<uint64_t, int> >, vector<int> > > clauses;
		int v=findREDPartitionForLit(lit, n, clauses, mcost);
		// log("returned ", v);
		if (v) {
			bool ok= (v>1);
			uint64_t cost=0;
			vector<int> cost_labels;
			if (ok) {
				// try without labels...
				vector<int> satisfied;
				for (int i=0; i<n; ++i) {
					if (!rLog.requestTime(rLog.activeTechnique)) return 0;
					
					vector<int> clause;
					uint64_t min_cost = 0;
					int min_cost_l = 0;
					for (auto& l : clauses[i].F) {
						if (!l.F) clause.push_back(l.S);
						else if (!min_cost) {
							min_cost=l.F;
							min_cost_l=l.S;
						}
					}
					if (checkExtendedPositiveReduct(clause, clauses[i].S)) {
						
						int sz = clause.size();
						if (int removed = trimReductClause(clause)) {
							rLog.removeClause(removed);
						}
						stats["redundant_clauses"]+=1;
						if (clause.size()) {
							stats["added_reduct_clauses"]+=1;
							if ((int)clause.size() != sz) {
								stats["trimmed_reduct_clauses"]+=1;
								stats["reduct_clause_size_before_trimming"]+=sz;
								stats["reduct_clause_size_after_trimming"]+=clause.size();
							}
							pi.addClause(clause, HARDWEIGHT);
							satSolver->addClause(clause);
							satisfied.push_back(i);
						} else {
							stats["trim_reduct_clause_conflicts"]+=1;
						}
					} else {
						if (--UNSATChecksLeft==0) return 0;
						
						cost += min_cost;
						cost_labels.push_back(min_cost_l);
						if (cost>mcost) break;
					}
				}
				satisfied.push_back(clauses.size());
				for (unsigned i=0, k=0, j=0; k<satisfied.size(); ++k) {
					while ((int)i<satisfied[k]) {
						swap(clauses[j++], clauses[i++]);
					}
					++i;
				}
				clauses.resize(clauses.size()-satisfied.size()+1);
			}
			if (!clauses.size()) { // can do subsuming
				// remove original clauses
				vector<int> toRemove(pi.litClauses[litNegation(lit)]);
				for (int c : toRemove) {
					pi.removeClause(c);
				}
				// harden lb
				vector<int> cl={lit};
				satSolver->addClause(cl);
				int removed = setVariable(litVariable(lit), litValue(lit));

				rLog.removeClause(removed);
				rLog.removeLabel(1);
				return n;
			}
			if (cost > mcost) continue;
			ok=1;
			sort(cost_labels.begin(), cost_labels.end());
			for (unsigned i=0; i<clauses.size(); ++i) {
				bool found = 0;
				vector<int> clause;
				for (auto& l : clauses[i].F) {
					if (!rLog.requestTime(rLog.activeTechnique)) return 0;
					if (!l.F) clause.push_back(l.S);
					else {
						if (!binary_search(cost_labels.begin(), cost_labels.end(), l.S)) {
							cost += l.F;
							if (cost>mcost) break;
						}
						clause.push_back(l.S);
						
						if (checkExtendedPositiveReduct(clause, clauses[i].S)) {
							found = 1;
							break;
						} else {
							if (--UNSATChecksLeft==0) return 0;
						}
					}
				}
				if (!found) {
					ok=0;
					break;
				}
			}
			if (ok) { // can harden label
				vector<int> cl={lit};
				satSolver->addClause(cl);
				int removed=setVariable(litVariable(lit), litValue(lit));
				
				rLog.removeClause(removed);
				rLog.removeLabel(1);
				return -n;
			}
		}
	}
	return 0;
}

int Preprocessor::tryREDOnClause(int c) {
	// check if clause c is redundant and if it is, remove it
	// TODO: reconstruction
	// TODO: doesn't work properly with UP checks... 
	vector<int> cl;
	for (int l : pi.clauses[c].lit) {
		if (pi.isLabel[litVariable(l)]) continue;
		cl.push_back(l);
	}
	vector<int> tmp;
	if (checkFilteredPositiveReduct(cl, tmp)) {
		pi.removeClause(c);
		rLog.removeClause(1);
		return 1;
	}
	return 0;
}

int Preprocessor::modelCuttingNewLit(vector<int>& clause, vector<int>& bclause, vector<int>& up, int litsInClause, int lit) {
	if (lit!=-1) {
		if (opt.MRED_trimBefore && binary_search(up.begin(), up.end(), litNegation(lit))) return 0;
		clause.push_back(lit);
	}
	if ((int)clause.size()==litsInClause || lit==-1) {
		vector<int> tmp;
		if (checkExtendedPositiveReduct(clause, tmp)) {
			stats["modelcutting_clauses"]+=1;
			stats["modelcutting_clause_sizes"]+=clause.size();
			
			int sz = clause.size();
			int removed=0;
			if (!opt.MRED_trimBefore && (removed = trimReductClause(clause))) {
				rLog.removeClause(removed);
			}
			
			stats["redundant_clauses"]+=1;
			if (clause.size()) {
				stats["added_reduct_clauses"]+=1;
				if ((int)clause.size() != sz) {
					stats["trimmed_reduct_clauses"]+=1;
					stats["before_trimming"]+=sz;
					stats["after_trimming"]+=clause.size();
				}
				rLog.removeClause(-1);
				pi.addClause(clause, HARDWEIGHT);
				satSolver->addClause(clause);
			} else {
				stats["trim_reduct_clause_conflicts"]+=1;
			}
			clause.clear();
			bclause.clear();
			up.clear();
			return 1;
		}
		clause.clear();
		bclause.clear();
		up.clear();
	} else if (opt.MRED_trimBefore) {
		bclause.push_back(litNegation(lit));
		up.clear();
		satSolver->propagate(bclause, up, 2);
		sort(up.begin(), up.end());
	}
	return 0;
}

int Preprocessor::tryModelCuttingRED() {
	vector<bool> model;
	findGoodModel(model, opt.MRED_asearchIterLimit, opt.MRED_improveTimeLimit, opt.MRED_satLikeTries, opt.MRED_satLikeTimeLimit);
	if (model.empty()) return 0;
	
	stats["modelcutting_models"]+=0;
	stats["modelcutting_clauses"]+=0;
	stats["modelcutting_clause_sizes"]+=0;
	stats["tryModelCuttingRED_modelFound"]+=1;

	int nofReducts=0;
	for (int litsInClause = opt.MRED_minLitsInClause; nofReducts==0 && litsInClause < opt.MRED_maxLitsInClause; litsInClause*=2) {
		log("try to cut a model, litsInClause=", litsInClause);
		vector<int> clause;
		vector<int> bclause;
		vector<int> tmp;
		vector<int> up;
		if (!opt.MRED_randomizedTries) {
			for (int v=0; v<pi.vars; ++v) {
				if (!rLog.requestTime(Log::Technique::MRED)) return nofReducts;
				if (pi.isVarRemoved(v) || pi.isLabel[v]) continue;
				int lit = model[v] ? negLit(v) : posLit(v);
				bool nf = nofReducts==0;
				nofReducts += modelCuttingNewLit(clause, bclause, up, litsInClause, lit);
				if (nf && nofReducts) stats["modelcutting_models"]+=1;
			}
			if (clause.size()) {
				bool nf = nofReducts==0;
				nofReducts += modelCuttingNewLit(clause, bclause, up, litsInClause, -1);
				if (nf && nofReducts) stats["modelcutting_models"]+=1;
			}
		} else {
			vector<int> lits;
			for (int v=0; v<pi.vars;++v) {
				if (pi.isVarRemoved(v) || pi.isLabel[v]) continue;
				lits.push_back(model[v]?negLit(v):posLit(v));
			}
			for (int tries=0; nofReducts==0 && tries<opt.MRED_randomizedTries; ++tries) {
				if (!rLog.requestTime(Log::Technique::MRED)) return nofReducts;
				random_shuffle(lits.begin(), lits.end());
				clause.clear();
				for (int lit : lits) {
					if (!rLog.requestTime(Log::Technique::MRED)) return nofReducts;
					bool nf = nofReducts==0;
					nofReducts += modelCuttingNewLit(clause, bclause, up, litsInClause, lit);
					if (nf && nofReducts) stats["modelcutting_models"]+=1;
				}
				if (clause.size()) {
					bool nf = nofReducts==0;
					nofReducts += modelCuttingNewLit(clause, bclause, up, litsInClause, -1);
					if (nf && nofReducts) stats["modelcutting_models"]+=1;
				}
			}
		}
	}
	return nofReducts;
}


int Preprocessor::tryUPLitRED(int lit) {
	if (pi.isLitLabel(litNegation(lit)) || pi.isVarRemoved(litVariable(lit))) return 0;
	
	vector<int> a={lit};
	vector<int> l;
	if (checkExtendedPositiveReduct(a,l)) {
		if (pi.isLabel[litVariable(lit)]) {
			rLog.removeLabel(1);
		} else {
			rLog.removeVariable(1);
		}
		satSolver->addClause(a);
		rLog.removeClause(setVariable(litVariable(lit), litValue(lit)));
		return 1;
	}
	return 0;
}

int Preprocessor::doLabelRED() {
	rLog.startTechnique(Log::Technique::LRED);
	prepareSatSolver();
	
	int satb=redSatSolverCalls;
	int redLabels1=0;
	int redLabels2=0;
	int redLits=0;
	stats["doLabelRED"]+=1;
	stats["red_literals"]+=0;
	stats["red_labels"]+=0;
	stats["red_labels_a"]+=0;
	stats["red_labels_b"]+=0;
	stats["added_reduct_clauses"]+=0;
	stats["reduct_clause_size_before_trimming"]+=0;
	stats["reduct_clause_size_after_trimming"]+=0;
	stats["trimmed_reduct_clauses"]+=0;
	stats["trim_reduct_clause_conflicts"]+=0;
	stats["check_filtered_positive_reduct_size"]+=0;
	stats["check_filtered_positive_reduct_timeout"]+=0;
	stats["check_filtered_positive_reduct_timeout_size"]+=0;
	stats["check_filtered_positive_reduct_sat"]+=0;
	stats["check_filtered_positive_reduct_sat_size"]+=0;
	stats["check_filtered_positive_reduct_unsat"]+=0;
	stats["check_filtered_positive_reduct_unsat_size"]+=0;
	stats["redundant_clauses"]+=0;
	stats["extended_good"]+=0;
	stats["extended_bad"]+=0;
	stats["trimmed_filtered_good"]+=0;
	stats["trimmed_filtered_bad"]+=0;
	stats["trimmed_filtered_f2_good"]+=0;
	stats["trimmed_filtered_f2_bad"]+=0;
	stats["filtered_good"]+=0;
	stats["filtered_f2_good"]+=0;
	for (int l=0; l<2*pi.vars; ++l) {
		if (!rLog.requestTime(Log::Technique::LRED)) {
			rLog.stopTechnique(Log::Technique::LRED);
			stats["LabelRed_timeouts_pos"]+=float(l)/2/pi.vars;
			log(redLabels1+redLabels2, " labels hardened by LRED: ", redLabels1, ", ", redLabels2, ", ", redLits, ", sat solver calls: ", redSatSolverCalls-satb);
			return redLabels1+redLabels2;
		}
	
		if (pi.isLitLabel(l)) {
			int v=tryREDOnLit(l);
			if (v>0) {
				log("Label hardened by LRED, technique a, n=", v);
				stats["red_labels"]+=1;
				stats["red_labels_a"]+=1;
				++redLabels1;
			} else if (v<0) {
				stats["red_labels"]+=1;
				stats["red_labels_b"]+=1;
				++redLabels2;
			}
		} else if (!pi.isLitLabel(litNegation(l))) {
			continue;
			int v=tryREDOnLit(l);
			if (v) {
				stats["red_literals"]+=1;
// 				log("Literal hardened by LRED, n=", v);
				++redLits;
			}
		}
	}
	stats["LabelRed_timeouts_pos"]+=1;
	rLog.stopTechnique(Log::Technique::LRED);
	log(redLabels1+redLabels2, " labels hardened by LRED: ", redLabels1, ", ", redLabels2, ", ", redLits, " literals hardened by LRED, sat solver calls: ", redSatSolverCalls-satb);
	return redLabels1 + redLabels2 + redLits;
}

int Preprocessor::doClauseRED() {
	rLog.startTechnique(Log::Technique::CRED);
	int redClauses=0;
	int satb=redSatSolverCalls;
	for (unsigned c=0; c<pi.clauses.size(); ++c) {
		if (!rLog.requestTime(Log::Technique::CRED)) {
			rLog.stopTechnique(Log::Technique::CRED);
			log(redClauses, " clauses removed by CRED, sat solver calls: ", redSatSolverCalls-satb);
			return redClauses;
		}
		if (pi.isClauseRemoved(c)) continue;
		if (pi.clauses[c].lit.size()<20) continue; // magic constant, TODO: tune
		redClauses+=tryREDOnClause(c);
	}
	rLog.stopTechnique(Log::Technique::CRED);
	log(redClauses, " clauses removed by RED, sat solver calls: ", redSatSolverCalls-satb);
	return redClauses;
}


int Preprocessor::doModelCuttingRED() {
	rLog.startTechnique(Log::Technique::MRED);
	if (!rLog.requestTime(Log::Technique::MRED)) {
		rLog.stopTechnique(Log::Technique::MRED);
		return 0;
	}
	prepareSatSolver();
	int rv = tryModelCuttingRED();
	rLog.stopTechnique(Log::Technique::MRED);
	return rv;

}

int Preprocessor::doUPLitRED() {
	rLog.startTechnique(Log::Technique::URED);
	if (!rLog.requestTime(Log::Technique::URED)) {
		rLog.stopTechnique(Log::Technique::URED);
		return 0;
	}
	
	prepareSatSolver();
	
	
	vector<pair<int, int> > vars;
	for (int var=0; var<pi.vars; ++var) {
		if (pi.isVarRemoved(var)) continue;
		vars.push_back({-((int)pi.litClauses[posLit(var)].size()+pi.litClauses[negLit(var)].size()), var});
	}
	sort(vars.begin(), vars.end());
	
	int removedVars = 0;
	for (unsigned i=0; i<vars.size(); ++i) {
		if (rLog.requestTime(Log::Technique::URED)) removedVars += tryUPLitRED(posLit(vars[i].S));
		if (rLog.requestTime(Log::Technique::URED)) removedVars += tryUPLitRED(negLit(vars[i].S));
	}
	
	rLog.stopTechnique(Log::Technique::URED);
	return removedVars;
}
