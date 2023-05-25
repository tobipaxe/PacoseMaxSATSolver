// Label hardening, uses SAT-solver-calls to find models...



int Preprocessor::tryHARD() {	
	vector<pair<uint64_t, int> > labels;
	for (int lit=0; lit<2*pi.vars; ++lit) {
		if (pi.isVarRemoved(litVariable(lit))) continue;
		if (pi.isLitLabel(lit)) {
			labels.push_back({pi.labelWeight(litVariable(lit)), lit});
		}
	}
	sort(labels.begin(), labels.end());
	
	vector<bool> best_model;

	
	uint64_t best_cost = findGoodModel(best_model, opt.HARD_asearchIterLimit, opt.HARD_improveTimeLimit, opt.HARD_satLikeTries, opt.HARD_satLikeTimeLimit);
	
	int removed=0;
	int removedLabels=0;
	for (int i=(int)labels.size()-1;i>=0 && labels[i].F>best_cost;--i) {
		removed+=setVariable(litVariable(labels[i].S), litValue(labels[i].S));
		++removedLabels;
	}
	rLog.removeLabel(removedLabels);
	log("Hardening done, removed ", removedLabels, " labels, ", removed, " clauses.");
	return removed;
}




int Preprocessor::doHARD() {
	rLog.startTechnique(Log::Technique::HARD);
	if (!rLog.requestTime(Log::Technique::HARD)) {
		rLog.stopTechnique(Log::Technique::HARD);
		return 0;
	}
	
	stats["doHARD"]+=1;
	
	
	int removedClauses=0;
	if (rLog.requestTime(Log::Technique::HARD)) removedClauses=tryHARD();
	
	rLog.stopTechnique(Log::Technique::HARD);
	return removedClauses;
}
