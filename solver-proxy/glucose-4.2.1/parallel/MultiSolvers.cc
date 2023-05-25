/***************************************************************************************[MultiSolvers.cc]
 Glucose -- Copyright (c) 2009-2014, Gilles Audemard, Laurent Simon
                                CRIL - Univ. Artois, France
                                LRI  - Univ. Paris Sud, France (2009-2013)
                                Labri - Univ. Bordeaux, France

 Syrup (Glucose Parallel) -- Copyright (c) 2013-2014, Gilles Audemard, Laurent Simon
                                CRIL - Univ. Artois, France
                                Labri - Univ. Bordeaux, France

Glucose sources are based on MiniSat (see below MiniSat copyrights). Permissions and copyrights of
Glucose (sources until 2013, Glucose 3.0, single core) are exactly the same as Minisat on which it
is based on. (see below).

Glucose-Syrup sources are based on another copyright. Permissions and copyrights for the parallel
version of Glucose-Syrup (the "Software") are granted, free of charge, to deal with the Software
without restriction, including the rights to use, copy, modify, merge, publish, distribute,
sublicence, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

- The above and below copyrights notices and this permission notice shall be included in all
copies or substantial portions of the Software;
- The parallel version of Glucose (all files modified since Glucose 3.0 releases, 2013) cannot
be used in any competitive event (sat competitions/evaluations) without the express permission of
the authors (Gilles Audemard / Laurent Simon). This is also the case for any competitive event
using Glucose Parallel as an embedded SAT engine (single core or not).


--------------- Original Minisat Copyrights

Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#include <pthread.h>
#include "../parallel/MultiSolvers.h"
#include "../mtl/Sort.h"
#include "../utils/System.h"
#include "../simp/SimpSolver.h"
#include <errno.h>
#include <string.h>
#include "../parallel/SolverConfiguration.h"

using namespace Glucose421;

extern const char *_parallel;
extern const char *_cunstable;

// Shared with ClausesBuffer.cc
BoolOption opt_whenFullRemoveOlder(_parallel, "removeolder", "When the FIFO for exchanging clauses between threads is full, remove older clauses", false);
IntOption opt_fifoSizeByCore(_parallel, "fifosize", "Size of the FIFO structure for exchanging clauses between threads, by threads", 100000);

// Shared options with Solver.cc
BoolOption opt_dontExportDirectReusedClauses(_cunstable, "reusedClauses", "Don't export directly reused clauses", false);
BoolOption opt_plingeling(_cunstable, "plingeling", "plingeling strategy for sharing clauses (exploratory feature)", false);

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


static inline double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double) ru.ru_utime.tv_sec + (double) ru.ru_utime.tv_usec / 1000000;
}


void MultiSolvers::informEnd(lbool res) {
    result = res;
    pthread_cond_broadcast(&cfinished);
}

MultiSolvers::MultiSolvers(size_t threads):
        use_simplification(true), ok(true), nbsolvers(threads),
        allClonesAreBuilt(0), showModel(false), winner(-1), var_decay(1 / 0.95), clause_decay(1 / 0.999), cla_inc(1), var_inc(1), random_var_freq(0.02), restart_first(100),
        restart_inc(1.5), learntsize_factor((double) 1 / (double) 3), learntsize_inc(1.1), expensive_ccmin(true), polarity_mode(polarity_false), maxmemory(0),
        verb(0), verbEveryConflicts(10000), numvar(0), numclauses(0) {
    result = l_Undef;
    this->sharedcomp = new SharedCompanion(threads);

    for (size_t i = 0u; i < nbsolvers; ++i)
	{
		ParallelSolver* s = new ParallelSolver(i);
		s->verbosity = 0; // No reportf in solvers... All is done in MultiSolver
		s->sharedcomp = this->sharedcomp;
		this->sharedcomp->addSolver(s);
		solvers.push(s);
	}

    SolverConfiguration::configure(this, nbsolvers);

    pthread_mutex_init(&m, NULL);  //PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&mfinished, NULL); //PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_init(&cfinished, NULL);
}

MultiSolvers::~MultiSolvers()
{
    for (size_t i = 0u; i < nbsolvers; ++i)
	{
		delete solvers[i];
	}

	delete this->sharedcomp;

	pthread_cond_destroy(&cfinished);
	pthread_mutex_destroy(&m);
	pthread_mutex_destroy(&mfinished);
}

/**
 * Generate All solvers
 */

void MultiSolvers::generateAllSolvers() {
    for(int i = solvers.size(); i < nbsolvers; i++) {
        ParallelSolver *s = (ParallelSolver *) solvers[0]->clone();
        solvers.push(s);
        s->verbosity = 0; // No reportf in solvers... All is done in MultiSolver
        s->setThreadNumber(i);
        s->sharedcomp = this->sharedcomp;
        this->sharedcomp->addSolver(s);
        assert(solvers[i]->threadNumber() == i);
    }

    SolverConfiguration::configure(this, nbsolvers);
}

/**
 * Choose solver for threads i (if no given in command line see above)
 */
ParallelSolver *MultiSolvers::retrieveSolver(int i) {
    return new ParallelSolver(i);
}

Var MultiSolvers::newVar(bool sign, bool dvar) {
    assert(solvers[0] != NULL);
    numvar++;
    int v;
    sharedcomp->newVar(sign);
	for(int i = 0; i < nbsolvers; i++) {
		v = solvers[i]->newVar(sign, dvar);
	}
    return numvar;
}

bool MultiSolvers::addClause_(vec<Lit> &ps) {
    assert(solvers[0] != NULL); // There is at least one solver.
    // Check if clause is satisfied and remove false/duplicate literals:
    if(!okay()) return false;

    sort(ps);
    Lit p;
    int i, j;
    for(i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if(solvers[0]->value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if(solvers[0]->value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);

    if(ps.size() == 0) {
        return ok = false;
    }
    else if(ps.size() == 1) {
        assert(solvers[0]->value(ps[0]) == l_Undef);

        for(int i = 0; i < nbsolvers; i++)
            solvers[i]->uncheckedEnqueue(ps[0]);

        return ok = ((solvers[0]->propagate()) == CRef_Undef); // checks only main solver here for propagation constradiction
    } else {
        for(int i = 0; i < nbsolvers; i++) {
            solvers[i]->addClause(ps);
        }
        numclauses++;
    }
    return true;
}


bool MultiSolvers::simplify() {
    assert(solvers.size() == 1);

    if(!okay()) return false;
    return ok = solvers[0]->simplify();
}

bool MultiSolvers::eliminate()
{
    assert(solvers.size() == 1);

    SimpSolver *s = (SimpSolver *) getPrimarySolver();
    s->use_simplification = use_simplification;
    if(!use_simplification) return true;

    return s->eliminate(true);
}

// TODO: Use a template here
void *localLaunch(void *arg) {
    ParallelSolver *s = (ParallelSolver *) arg;

	printf("Starting solver thread %d\n", s->threadNumber());
    (void) s->solve_(s->_do_simp, s->_turn_off_simp);
	printf("Finished solver thread %d\n", s->threadNumber());

    pthread_exit(NULL);
}

void MultiSolvers::printStats() {
    static int nbprinted = 1;
    double cpu_time = cpuTime();
    printf("c\n");

    printf("c |-------------------------------------------------------------------------------------------------------|\n");
    printf("c | id | starts | decisions  |  confls    |  Init T  |  learnts | exported | imported | promoted |    %%   | \n");
    printf("c |-------------------------------------------------------------------------------------------------------|\n");

    //printf("%.0fs | ",cpu_time);
    for(int i = 0; i < solvers.size(); i++) {
        solvers[i]->reportProgress();
        //printf(" %2d: %12ld confl. |", i,  (long int) solvers[i]->conflicts);
    }
    long long int totalconf = 0;
    long long int totalprop = 0;
    for(int i = 0; i < solvers.size(); i++) {
        totalconf += (long int) solvers[i]->conflicts;
        totalprop += solvers[i]->propagations;
    }
    printf("c \n");

    printf("c synthesis %11lld conflicts %11lld propagations %8.0f conflicts/sec %8.0f propagations/sec\n",
           totalconf, totalprop, (double) totalconf / cpu_time, (double) totalprop / cpu_time);


    nbprinted++;
}

// Still a ugly function... To be rewritten with some statistics class some day
void MultiSolvers::printFinalStats() {
    sharedcomp->printStats();
    printf("c\nc\n");
    printf("c\n");
    printf("c |---------------------------------------- FINAL STATS --------------------------------------------------|\n");
    printf("c\n");

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Threads       |      Total      ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10d ", i);
    }
    printf("|\n");

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");


//--
    printf("c | Conflicts     ");
    long long int totalconf = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalconf += solvers[i]->conflicts;
    printf("| %15lld ", totalconf);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Decisions     ");
    long long int totaldecs = 0;
    for(int i = 0; i < solvers.size(); i++)
        totaldecs += solvers[i]->decisions;
    printf("| %15lld ", totaldecs);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->decisions);
    printf("|\n");

    //--
    printf("c | Propagations  ");
    long long int totalprops = 0;
    for(int i = 0; i < solvers.size(); i++)
        totalprops += solvers[i]->propagations;
    printf("| %15lld ", totalprops);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->propagations);
    printf("|\n");


    printf("c | Avg_Trail     ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumTrail] / solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Avg_DL        ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumDecisionLevels] / solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Avg_Res       ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumRes] / solvers[i]->conflicts);
    printf("|\n");

    //--
    printf("c | Avg_Res_Seen  ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->conflicts==0 ? 0 : solvers[i]->stats[sumResSeen] / solvers[i]->conflicts);
    printf("|\n");

    //--

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Exported      ");
    uint64_t exported = 0;
    for(int i = 0; i < solvers.size(); i++)
        exported += solvers[i]->stats[nbexported];
    printf("| %15" PRIu64" ", exported);

    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbexported]);
    printf("|\n");
//--
    printf("c | Imported      ");
    uint64_t imported = 0;
    for(int i = 0; i < solvers.size(); i++)
        imported += solvers[i]->stats[nbimported];
    printf("| %15" PRIu64" ", imported);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbimported]);
    printf("|\n");
//--

    printf("c | Good          ");
    uint64_t importedGood = 0;
    for(int i = 0; i < solvers.size(); i++)
        importedGood += solvers[i]->stats[nbImportedGoodClauses];
    printf("| %15" PRIu64" ", importedGood);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbImportedGoodClauses]);
    printf("|\n");
//--

    printf("c | Purge         ");
    uint64_t importedPurg = 0;
    for(int i = 0; i < solvers.size(); i++)
        importedPurg += solvers[i]->stats[nbimportedInPurgatory];
    printf("| %15" PRIu64" ", importedPurg);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbimportedInPurgatory]);
    printf("|\n");
//--

    printf("c | Promoted      ");
    uint64_t promoted = 0;
    for(int i = 0; i < solvers.size(); i++)
        promoted += solvers[i]->stats[nbPromoted];
    printf("| %15" PRIu64" ", promoted);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbPromoted]);
    printf("|\n");
//--

    printf("c | Remove_Imp    ");
    uint64_t removedimported = 0;
    for(int i = 0; i < solvers.size(); i++)
        removedimported += solvers[i]->stats[nbRemovedUnaryWatchedClauses];
    printf("| %15" PRIu64" ", removedimported);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->stats[nbRemovedUnaryWatchedClauses]);
    printf("|\n");
//--

    printf("c | Blocked_Reuse ");
    uint64_t blockedreused = 0;
    for(int i = 0; i < solvers.size(); i++)
        blockedreused += solvers[i]->nbNotExportedBecauseDirectlyReused;
    printf("| %15" PRIu64" ", blockedreused);
    for(int i = 0; i < solvers.size(); i++)
        printf("| %10" PRIu64" ", solvers[i]->nbNotExportedBecauseDirectlyReused);
    printf("|\n");
//--
    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Unaries       ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[nbUn]);
    }
    printf("|\n");
//--

    printf("c | Binaries      ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[nbBin]);
    }
    printf("|\n");
//--


    printf("c | Glues         ");
    printf("|                 ");
    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[nbDL2]);
    }
    printf("|\n");
//--

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");

    printf("c | Orig_Seen     ");
    uint64_t origseen = 0;

    for(int i = 0; i < solvers.size(); i++) {
        origseen += solvers[i]->stats[originalClausesSeen];
    }
    printf("| %13" PRIu64" %% ", (nClauses() > 0 && solvers.size() > 0) ?origseen * 100 / nClauses() / solvers.size() : 0);

    for(int i = 0; i < solvers.size(); i++) {
        printf("| %10" PRIu64" ", solvers[i]->stats[originalClausesSeen]);
    }

    printf("|\n");


    int winner = -1;
    for(int i = 0; i < solvers.size(); i++) {
        if(sharedcomp->winner() == solvers[i])
            winner = i;
    }

//--
    if(winner != -1) {
        printf("c | Diff Orig seen");
        printf("|                 ");

        for(int i = 0; i < solvers.size(); i++) {
            if(i == winner) {
                printf("|      X     ");
                continue;
            }
            if(solvers[i]->stats[originalClausesSeen] > solvers[winner]->stats[originalClausesSeen])
                printf("| %10" PRIu64" ", solvers[i]->stats[originalClausesSeen] - solvers[winner]->stats[originalClausesSeen]);
            else
                printf("| -%9" PRIu64" ", solvers[winner]->stats[originalClausesSeen] - solvers[i]->stats[originalClausesSeen]);

        }

        printf("|\n");
    }


//--

    if(winner != -1) {
        int sum = 0;
        printf("c | Hamming       ");
        for(int i = 0; i < solvers.size(); i++) {
            if(i == winner)
                continue;
            int nb = 0;
            for(int j = 0; j < nVars(); j++) {
                if(solvers[i]->valuePhase(j) != solvers[winner]->valuePhase(j)) nb++;
            }
            sum += nb;

        }
        sum = sum / (solvers.size() > 1 ? solvers.size() - 1 : 1);

        printf("| %13d %% ", sum * 100 / nVars());

        for(int i = 0; i < solvers.size(); i++) {
            if(i == winner) {
                printf("|      X     ");
                continue;
            }
            int nb = 0;
            for(int j = 0; j < nVars(); j++) {
                if(solvers[i]->valuePhase(j) != solvers[winner]->valuePhase(j)) nb++;
            }
            printf("| %10d ", nb);
            sum += nb;

        }
        printf("|\n");
    }

    printf("c |---------------|-----------------");
    for(int i = 0; i < solvers.size(); i++)
        printf("|------------");
    printf("|\n");
}

lbool MultiSolvers::solve_(bool do_simp, bool turn_off_simp) {
    pthread_attr_t thAttr;
    int i;

	for(i = 0; i < nbsolvers; i++)
	{
		solvers[i]->_do_simp = do_simp;
		solvers[i]->_turn_off_simp = turn_off_simp;
	}

	generateAllSolvers();

    model.clear();

    /* Initialize and set thread detached attribute */
    pthread_attr_init(&thAttr);
    pthread_attr_setdetachstate(&thAttr, PTHREAD_CREATE_JOINABLE);

    // Launching all solvers
    for(i = 0; i < nbsolvers; i++) {
        pthread_t *pt = (pthread_t *) malloc(sizeof(pthread_t));
        threads.push(pt);
        solvers[i]->pmfinished = &mfinished;
        solvers[i]->pcfinished = &cfinished;
        pthread_create(threads[i], &thAttr, &localLaunch, (void *) solvers[i]);
    }

    bool done = false;
    bool adjustedlimitonce = false;
	uint64_t spinCounter = 0u;

    (void) pthread_mutex_lock(&m);
    while(!done) {
        struct timespec timeout;
		clock_gettime(CLOCK_REALTIME, &timeout);
        timeout.tv_sec += 0;
        timeout.tv_nsec += 1000000ul; // 1 ms
		while (timeout.tv_nsec > 1000000000ul)
		{
			timeout.tv_nsec -= 1000000000ul;
			timeout.tv_sec++;
		}

		int waitResult = pthread_cond_timedwait(&cfinished, &mfinished, &timeout);
		switch (waitResult)
		{
			case 0:
				if (!sharedcomp->jobFinished())
				{
					continue;
				}
				done = true;
				break;
			case ETIMEDOUT:
				if (!sharedcomp->jobFinished())
				{
					if ((++spinCounter % 5000ul) == 0u) // Wait 5 seconds
					{
						spinCounter = 0u;
						printStats();
					}
					continue;
				}
				done = true;
				break;
			case EINVAL:
				printf("Invalid arguments");
				break;
			case EPERM:
				printf("Non-owning thread");
				break;
		}

		float mem = memUsed();
        if(verb >= 1) printf("c Total Memory so far : %.2fMb\n", mem);
        if((maxmemory > 0) && (mem > maxmemory) && !sharedcomp->panicMode)
            printf("c ** reduceDB switching to Panic Mode due to memory limitations !\n");
			sharedcomp->panicMode = true;

        if(!done && !adjustedlimitonce) {
            uint64_t sumconf = 0;
            uint64_t sumimported = 0;
            for(int i = 0; i < nbsolvers; i++) {
                sumconf += solvers[i]->conflicts;
                sumimported += solvers[i]->stats[nbimported];
            }
            if(sumconf > 10000000 && sumimported > 4 * sumconf) { // too many many imported clauses (after a while)
                for(int i = 0; i < nbsolvers; i++) { // we have like 32 threads, so we need to export just very good clauses
                    solvers[i]->goodlimitlbd -= 2;
                    solvers[i]->goodlimitsize -= 4;
                }
                adjustedlimitonce = true;
                printf("c adjusting (once) the limits to send fewer clauses.\n");
            }
        }
    }
    (void) pthread_mutex_unlock(&m);

    for(i = 0; i < nbsolvers; i++) { // Wait for all threads to finish
        pthread_join(*threads[i], NULL);
    }

	for(i = 0; i < nbsolvers; i++) {
		delete threads[i];
	}

    assert(sharedcomp != NULL);
    result = sharedcomp->jobStatus;
    if(result == l_True) {
        sharedcomp->jobFinishedBy->extendModel();
        int n = sharedcomp->jobFinishedBy->nVars();
        model.growTo(n);
        for(int i = 0; i < n; i++) {
            model[i] = sharedcomp->jobFinishedBy->model[i];
            assert(model[i] != l_Undef);
        }
    }

    return result;
}

lbool MultiSolvers::solve(                                bool do_simp, bool turn_off_simp)
{
	for (size_t i = 0u; i < nbsolvers; ++i)
	{
		solvers[i]->budgetOff();
		solvers[i]->assumptions.clear();
	}
	return solve_(do_simp, turn_off_simp);
}

lbool MultiSolvers::solve(Lit p,                          bool do_simp, bool turn_off_simp)
{
	for (size_t i = 0u; i < nbsolvers; ++i)
	{
		solvers[i]->budgetOff();
		solvers[i]->assumptions.clear();
		solvers[i]->assumptions.push(p);
	}
	return solve_(do_simp, turn_off_simp);
}

lbool MultiSolvers::solve(Lit p, Lit q,                   bool do_simp, bool turn_off_simp)
{
	for (size_t i = 0u; i < nbsolvers; ++i)
	{
		solvers[i]->budgetOff();
		solvers[i]->assumptions.clear();
		solvers[i]->assumptions.push(p);
		solvers[i]->assumptions.push(q);
	}
	return solve_(do_simp, turn_off_simp);
}

lbool MultiSolvers::solve(Lit p, Lit q, Lit r,            bool do_simp, bool turn_off_simp)
{
	for (size_t i = 0u; i < nbsolvers; ++i)
	{
		solvers[i]->budgetOff();
		solvers[i]->assumptions.clear();
		solvers[i]->assumptions.push(p);
		solvers[i]->assumptions.push(q);
		solvers[i]->assumptions.push(r);
	}
	return solve_(do_simp, turn_off_simp);
}

lbool MultiSolvers::solve(const vec<Lit>& assumps,        bool do_simp, bool turn_off_simp)
{
	for (size_t i = 0u; i < nbsolvers; ++i)
	{
		solvers[i]->budgetOff();
		solvers[i]->assumptions.clear();
		assumps.copyTo(solvers[i]->assumptions);
	}
	return solve_(do_simp, turn_off_simp);
}

lbool MultiSolvers::solveLimited(const vec<Lit>& assumps, bool do_simp, bool turn_off_simp)
{
	for (size_t i = 0u; i < nbsolvers; ++i)
	{
		assumps.copyTo(solvers[i]->assumptions);
	}
	return solve_(do_simp, turn_off_simp);
}

ParallelSolver* MultiSolvers::getWinningSolver(void)
{
	return sharedcomp->jobFinishedBy;
}

ParallelSolver *MultiSolvers::getSolver(int thread)
{
	return solvers[thread];
}

void MultiSolvers::setSolverCount(int solvers)
{
	nbsolvers = solvers;
	this->sharedcomp->setNbThreads(solvers);
	generateAllSolvers();
}

int MultiSolvers::getSolverCount()
{
	return nbsolvers;
}

void MultiSolvers::interrupt()
{
	for (size_t i = 0u; i < nbsolvers; ++i)
	{
		solvers[i]->interrupt();
	}
}
