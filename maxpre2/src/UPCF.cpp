// Failed literal detection


int Preprocessor::tryFLD(int lit) {
	vector<int> a={lit};
	vector<int> c;
	if (!satSolver->propagate(a, c, 2)) {
		stats["FLD_harden"]+=1;
		if (pi.isLabel[litVariable(lit)] != VAR_UNDEFINED) {
			rLog.removeLabel(1);
		} else{
			rLog.removeVariable(1);
		}
		int removed = setVariable(litVariable(lit), !litValue(lit));
		rLog.removeClause(removed);
		return removed;
	}
	return 0;
}

int Preprocessor::doFLD() {
	prepareSatSolver();
	stats["FLD_harden"]+=0;	
	
	rLog.startTechnique(Log::Technique::FLD);	
	int redClauses=0;
	int removed=0;
	stats["doFLD_calls"]+=1;
	for (int l=0; l<2*pi.vars; ++l) {
		if (!rLog.requestTime(Log::Technique::FLD)) {
			stats["FLD_stop_position"]+=float(l)/2/pi.vars;
			rLog.stopTechnique(Log::Technique::FLD);
			return removed;
		}
		if (pi.isVarRemoved(litVariable(l))) continue;
		removed+=tryFLD(l);
	}
	stats["FLD_stop_position"]+=1;
	rLog.stopTechnique(Log::Technique::FLD);
	return removed;
}
