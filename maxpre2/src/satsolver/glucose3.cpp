#include <iostream>


#include "satsolverinterface.hpp"
#include "solvers/glucose3/core/Solver.h"

#ifndef GLUCOSE3_CPP
#define GLUCOSE3_CPP


class Glucose3: SATSolver {
private:
	GlucosePre::Solver solver;
	maxPreprocessor::Timer timer;
	
	int topv;
	
	void intsToLits(const vector<int>& lits, GlucosePre::vec<GlucosePre::Lit>& rLits) {
		rLits.capacity(lits.size());
		for (unsigned i=0;i<lits.size();++i) {
			while ((lits[i]>>1)>topv) {solver.newVar();++topv;}
			rLits.push_(GlucosePre::toLit(lits[i]));
		}
	}
	void litsToInts(GlucosePre::vec<GlucosePre::Lit>& lits, vector<int>& rLits) {
		rLits.resize(lits.size());
		for (int i=0; i<lits.size(); ++i) rLits[i]=GlucosePre::toInt(lits[i]);
	}
	
	void lboolsToBools(GlucosePre::vec<GlucosePre::lbool>& lbools, vector<bool>& bools) {
		bools.resize(lbools.size());
		for (int i=0;i<lbools.size(); ++i) bools[i] = (lbools[i] == l_TrueP ? 1 : 0);
	}
	int SATcalls;
	int UNSATcalls;
	int outOfBudgetCalls;
public:
	void addClause(const vector<int>& clause) {
		for (unsigned i=0;i<clause.size();++i) while ((clause[i]>>1)>topv) {solver.newVar();++topv;}
		
		if (clause.size()==1) {
			solver.addClause(GlucosePre::toLit(clause[0]));
		} else if (clause.size() == 2) {
			solver.addClause(GlucosePre::toLit(clause[0]), GlucosePre::toLit(clause[1]));
		} else if (clause.size() == 3) {
			solver.addClause(GlucosePre::toLit(clause[0]), GlucosePre::toLit(clause[1]), GlucosePre::toLit(clause[2]));
		} else {
			GlucosePre::vec<GlucosePre::Lit> clause_;
			intsToLits(clause, clause_);
			solver.addClause_(clause_);
		}
	}
	

	bool solve(vector<int>& assumptions) {
		GlucosePre::vec<GlucosePre::Lit> assumps;
		intsToLits(assumptions, assumps);
		
		timer.start();
		bool rv = solver.solve(assumps);
		timer.stop();
		
		if (rv) ++SATcalls;
		else    ++UNSATcalls;
		return rv;
	}
	
	int solveLimited(vector<int>& assumptions, int64_t maxPropagations, int64_t maxConflicts) {
		GlucosePre::vec<GlucosePre::Lit> assumps;
		intsToLits(assumptions, assumps);
		int rv = -1;
		solver.setConfBudget(maxConflicts);
		solver.setPropBudget(maxPropagations);
		
		timer.start();
		//Glucose::lbool srv = solver.solve(assumps)?l_True:l_False;
		GlucosePre::lbool srv = solver.solveLimited(assumps);
		timer.stop();
		
		if (srv == l_TrueP)     rv=1;
		else if (srv==l_FalseP) rv=0;
		else if (srv==l_UndefP) rv=-1;
		
		if (rv==1) ++SATcalls;
		else if (rv==0)   ++UNSATcalls;
		else if (rv==-1) ++outOfBudgetCalls;
		return rv;
		
	}
	
	int solveLimited(vector<int>& assumptions, double time_) { // estimated time seconds maximum time
		double time=time_;
		if (time/2 + 10 < time) time=time/2+10; // magic magic
		
		double atime=timer.getTime().count();
		
		double ttime=0.01+atime;
		int64_t propagations = solver.propagations+10000 ;
		int64_t conflicts = solver.conflicts+100;
		return solveLimited(assumptions, time/ttime * propagations, time/ttime * conflicts - 1);
	}
	
	bool solve() {
		
		timer.start();
		bool rv = solver.solve();
		timer.stop();
		
		return rv;
	}


	void getCore(vector<int>& retCore) {
		litsToInts(solver.conflict, retCore);
	}
	
	void getModel(vector<bool>& retModel) {
		lboolsToBools(solver.model, retModel);
	}
	
	bool propagate(vector<int>& assumptions, vector<int>& result, int phase_saving) {
		GlucosePre::vec<GlucosePre::Lit> assumps;
		intsToLits(assumptions, assumps);
		GlucosePre::vec<GlucosePre::Lit> res;
		bool rv = solver.prop_check(assumps, res, phase_saving);
		litsToInts(res, result);
		return rv;
	}
	
	bool testUPConflict(vector<int>& assumptions, int phase_saving) {
		GlucosePre::vec<GlucosePre::Lit> assumps;
		intsToLits(assumptions, assumps);
		return solver.prop_confl_check(assumps, phase_saving);
	}
	
	int numberOfSATcalls() {
		return SATcalls;
	}
	int numberOfUNSATcalls() {
		return UNSATcalls;
	}
	
	Glucose3() : topv(-1), SATcalls(0), UNSATcalls(0), outOfBudgetCalls(0) {
		
	}
	
	~Glucose3() {
		
	}
};


#endif
