#ifndef Alg_SATLIke_h
#define Alg_SATLIke_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h> //these two h files are for timing in linux
#include <unistd.h>

using namespace std;

#define mypop(stack) stack[--stack##_fill_pointer]
#define mypush(item, stack) stack[stack##_fill_pointer++] = item

const float MY_RAND_MAX_FLOAT = 10000000.0;
const int MY_RAND_MAX_INT = 10000000;
const float BASIC_SCALE = 0.0000001; //1.0f/MY_RAND_MAX_FLOAT;

// Define a data structure for a literal.
struct lit
{
    int clause_num; //clause num, begin with 0
    int var_num;    //variable num, begin with 1
    bool sense;     //is 1 for true literals, 0 for false literals.
};

static struct tms start_time;
static double get_runtime()
{
    struct tms stop;
    times(&stop);
    return (double)(stop.tms_utime - start_time.tms_utime + stop.tms_stime - start_time.tms_stime) / sysconf(_SC_CLK_TCK);
}
static void start_timing()
{
    times(&start_time);
}

class Satlike
{
  public:
    /***********non-algorithmic information ****************/
    int problem_weighted;
    int partial; //1 if the instance has hard clauses, and 0 otherwise.
    int pure_sat;

    int max_clause_length;
    int min_clause_length;

    //size of the instance
    int num_vars;    //var index from 1 to num_vars
    int num_clauses; //clause index from 0 to num_clauses-1
    int num_hclauses;
    int num_sclauses;

    //steps and time
    int tries;
    int max_tries;
    int max_flips;
    int max_non_improve_flip;
    int step;

    int print_time;
    int cutoff_time;
    int prioup_time;
    double opt_time;

    /**********end non-algorithmic information*****************/
    /* literal arrays */
    lit **var_lit;         //var_lit[i][j] means the j'th literal of var i.
    int *var_lit_count;    //amount of literals of each var
    lit **clause_lit;      //clause_lit[i][j] means the j'th literal of clause i.
    int *clause_lit_count; // amount of literals in each clause

    /* Information about the variables. */
    long long *score;
    long long *time_stamp;
    int **var_neighbor;
    int *var_neighbor_count;
    int *neighbor_flag;
    int *temp_neighbor;
    bool if_using_neighbor;

    /* Information about the clauses */
    unsigned long long top_clause_weight;
    long long *org_clause_weight;
    long long total_soft_weight;
    long long *clause_weight;
    int *sat_count;
    int *sat_var;
    long long *clause_selected_count;
    int *best_soft_clause;

    //original unit clause stack
    lit *unit_clause;
    int unit_clause_count;

    //unsat clauses stack
    int *hardunsat_stack;          //store the unsat clause number
    int *index_in_hardunsat_stack; //which position is a clause in the unsat_stack
    int hardunsat_stack_fill_pointer;

    int *softunsat_stack;          //store the unsat clause number
    int *index_in_softunsat_stack; //which position is a clause in the unsat_stack
    int softunsat_stack_fill_pointer;

    //variables in unsat clauses
    int *unsatvar_stack;
    int unsatvar_stack_fill_pointer;
    int *index_in_unsatvar_stack;
    int *unsat_app_count; //a varible appears in how many unsat clauses

    //good decreasing variables (dscore>0 and confchange=1)
    int *goodvar_stack;
    int goodvar_stack_fill_pointer;
    int *already_in_goodvar_stack;

    /* Information about solution */
    int *cur_soln; //the current solution, with 1's for True variables, and 0's for False variables
    int *best_soln;
    int *local_opt_soln;
    int best_soln_feasible; //when find a feasible solution, this is marked as 1.
    int local_soln_feasible;
    int hard_unsat_nb;
    unsigned long long soft_unsat_weight;
    unsigned long long opt_unsat_weight;
    unsigned long long local_opt_unsat_weight;

    //clause weighting
    int *large_weight_clauses;
    int large_weight_clauses_count;
    int large_clause_count_threshold;

    int *soft_large_weight_clauses;
    int *already_in_soft_large_weight_stack;
    int soft_large_weight_clauses_count;
    int soft_large_clause_count_threshold;

    //tem data structure used in algorithm
    int *best_array;
    int best_count;
    int *temp_lit;

    //parameters used in algorithm
    float rwprob;
    float rdprob;
    float smooth_probability;
    int hd_count_threshold;
    int h_inc;
    int softclause_weight_threshold;

    //function used in algorithm
    void build_neighbor_relation();
    void allocate_memory();
    bool verify_sol();
    void increase_weights();
    void smooth_weights();
    void update_clause_weights();
    void unsat(int clause);
    void sat(int clause);
    void init(vector<int> &init_solution);
    void flip(int flipvar);
    void update_goodvarstack1(int flipvar);
    void update_goodvarstack2(int flipvar);
    int pick_var();

    Satlike();
    void settings();
//     void build_instance(char *filename);
    void build_instance(int numVars, int numClauses, unsigned long long topClauseweight, lit **satlike_clause, int *satlike_clause_lit_count, long long *satlike_clause_weight);
    void local_search(vector<int> &init_solution);
    void local_search_with_decimation(char *inputfile);
    void simple_print();
    void print_best_solution();
    void free_memory();
};

Satlike::Satlike() {}

void Satlike::settings()
{
    cutoff_time = 10;
    if (problem_weighted == 1)
    {
        cout << "c problem weighted = 1" << endl;
        max_tries = 100000000;
        max_flips = 10000000;
        max_non_improve_flip = max_flips;

        large_clause_count_threshold = 0;
        soft_large_clause_count_threshold = 0;

        rdprob = 0.01;
        hd_count_threshold = 15;
        rwprob = 0.1;
        smooth_probability = 0.01;

        if ((top_clause_weight / num_sclauses) > 10000)
        {
            h_inc = 300;
            softclause_weight_threshold = 500;
        }
        else
        {
            h_inc = 3;
            softclause_weight_threshold = 0;
        }

        if (num_vars > 2000)
        {
            rdprob = 0.01;
            hd_count_threshold = 15;
            rwprob = 0.1;
            smooth_probability = 0.0001;
        }
    }
    else
    {
        cout << "c problem weighted = 0" << endl;
        max_tries = 100000000;
        max_flips = 200000000;
        max_non_improve_flip = 10000000;

        large_clause_count_threshold = 0;
        soft_large_clause_count_threshold = 0;

        rdprob = 0.01;
        hd_count_threshold = 42;
        rwprob = 0.091;
        smooth_probability = 0.000003;

        h_inc = 1;
        softclause_weight_threshold = 400;

        if (num_vars < 1100) //for somall instances
        {
            h_inc = 1;
            softclause_weight_threshold = 0;
            rdprob = 0.01;
            hd_count_threshold = 15;
            rwprob = 0;
            smooth_probability = 0.01;
            return;
        }
    }
}

void Satlike::allocate_memory()
{
    int malloc_var_length = num_vars + 10;
    int malloc_clause_length = num_clauses + 10;

    unit_clause = new lit[malloc_clause_length];

    var_lit = new lit *[malloc_var_length];
    var_lit_count = new int[malloc_var_length];
    //clause_lit = new lit *[malloc_clause_length];
    //clause_lit_count = new int[malloc_clause_length];

    score = new long long[malloc_var_length];
    var_neighbor = new int *[malloc_var_length];
    var_neighbor_count = new int[malloc_var_length];
    time_stamp = new long long[malloc_var_length];
    neighbor_flag = new int[malloc_var_length];
    temp_neighbor = new int[malloc_var_length];

    //org_clause_weight = new long long[malloc_clause_length];
    clause_weight = new long long[malloc_clause_length];
    sat_count = new int[malloc_clause_length];
    sat_var = new int[malloc_clause_length];
    clause_selected_count = new long long[malloc_clause_length];
    best_soft_clause = new int[malloc_clause_length];

    hardunsat_stack = new int[malloc_clause_length];
    index_in_hardunsat_stack = new int[malloc_clause_length];
    softunsat_stack = new int[malloc_clause_length];
    index_in_softunsat_stack = new int[malloc_clause_length];

    unsatvar_stack = new int[malloc_var_length];
    index_in_unsatvar_stack = new int[malloc_var_length];
    unsat_app_count = new int[malloc_var_length];

    goodvar_stack = new int[malloc_var_length];
    already_in_goodvar_stack = new int[malloc_var_length];

    cur_soln = new int[malloc_var_length];
    best_soln = new int[malloc_var_length];
    local_opt_soln = new int[malloc_var_length];

    large_weight_clauses = new int[malloc_clause_length];
    soft_large_weight_clauses = new int[malloc_clause_length];
    already_in_soft_large_weight_stack = new int[malloc_clause_length];

    best_array = new int[malloc_var_length];
    temp_lit = new int[malloc_var_length];
}

void Satlike::free_memory()
{
    int i;
    for (i = 0; i < num_clauses; i++)
        delete[] clause_lit[i];

    for (i = 1; i <= num_vars; ++i)
    {
        delete[] var_lit[i];
        delete[] var_neighbor[i];
    }

    delete[] var_lit;
    delete[] var_lit_count;
    delete[] clause_lit;
    delete[] clause_lit_count;

    delete[] score;
    delete[] var_neighbor;
    delete[] var_neighbor_count;
    delete[] time_stamp;
    delete[] neighbor_flag;
    delete[] temp_neighbor;

    delete[] org_clause_weight;
    delete[] clause_weight;
    delete[] sat_count;
    delete[] sat_var;
    delete[] clause_selected_count;
    delete[] best_soft_clause;

    delete[] hardunsat_stack;
    delete[] index_in_hardunsat_stack;
    delete[] softunsat_stack;
    delete[] index_in_softunsat_stack;

    delete[] unsatvar_stack;
    delete[] index_in_unsatvar_stack;
    delete[] unsat_app_count;

    delete[] goodvar_stack;
    delete[] already_in_goodvar_stack;

    //delete [] fix;
    delete[] cur_soln;
    delete[] best_soln;
    delete[] local_opt_soln;

    delete[] large_weight_clauses;
    delete[] soft_large_weight_clauses;
    delete[] already_in_soft_large_weight_stack;

    delete[] best_array;
    delete[] temp_lit;
}

void Satlike::build_neighbor_relation()
{
    int i, j, count;
    int v, c, n;
    int temp_neighbor_count;

    for (v = 1; v <= num_vars; ++v) neighbor_flag[v]=0;
		
    for (v = 1; v <= num_vars; ++v)
    {
        neighbor_flag[v] = 1;
        temp_neighbor_count = 0;

        for (i = 0; i < var_lit_count[v]; ++i)
        {
            c = var_lit[v][i].clause_num;
            for (j = 0; j < clause_lit_count[c]; ++j)
            {
                n = clause_lit[c][j].var_num;
                if (neighbor_flag[n] != 1)
                {
                    neighbor_flag[n] = 1;
                    temp_neighbor[temp_neighbor_count++] = n;
                }
            }
        }

        neighbor_flag[v] = 0;

        var_neighbor[v] = new int[temp_neighbor_count];
        var_neighbor_count[v] = temp_neighbor_count;

        count = 0;
        for (i = 0; i < temp_neighbor_count; i++)
        {
            var_neighbor[v][count++] = temp_neighbor[i];
            neighbor_flag[temp_neighbor[i]] = 0;
        }
    }
}

/*
void Satlike::build_instance(char *filename)
{
    cout << "org build instance function" << endl;
}
*/

void Satlike::build_instance(int numVars, int numClauses, unsigned long long topClauseweight, lit **satlike_clause, int *satlike_clause_lit_count, long long *satlike_clause_weight)
{
    istringstream iss;
    string line;

    /*** build problem data structures of the instance ***/
    start_timing();

    num_vars = numVars;
    num_clauses = numClauses;
    top_clause_weight = topClauseweight;
    clause_lit = satlike_clause;
    clause_lit_count = satlike_clause_lit_count;
    org_clause_weight = satlike_clause_weight;

    allocate_memory();
    int v;

    for (v = 1; v <= num_vars; ++v)
    {
        var_lit_count[v] = 0;
        var_lit[v] = NULL;
        var_neighbor[v] = NULL;
    }

    problem_weighted = 0;
    partial = 0;
    num_hclauses = num_sclauses = 0;
    max_clause_length = 0;
    min_clause_length = 100000000;
    unit_clause_count = 0;
	total_soft_weight = 0;
    //int *redunt_test = new int[num_vars + 1];
    //memset(redunt_test, 0, sizeof(int) * num_vars + 1);
    //Now, read the clauses, one at a time.
    for (int i = 0; i < num_clauses; ++i)
    {
        for (int j = 0; j < clause_lit_count[i]; ++j)
        {
            //int temv = clause_lit[i][j].var_num;
            //int temsense = clause_lit[i][j].sense;
            var_lit_count[clause_lit[i][j].var_num]++;
        }

        if (org_clause_weight[i] != (long long)top_clause_weight)
        {
            if (org_clause_weight[i] != 1)
                problem_weighted = 1;
            total_soft_weight += org_clause_weight[i];
            num_sclauses++;
        }
        else
        {
            num_hclauses++;
            partial = 1;
        }
        if (clause_lit_count[i] == 1)
            unit_clause[unit_clause_count++] = clause_lit[i][0];

        if (clause_lit_count[i] > max_clause_length)
            max_clause_length = clause_lit_count[i];

        if (clause_lit_count[i] < min_clause_length)
            min_clause_length = clause_lit_count[i];
    }
    cout << "c befor creat var literal arrays " << endl;
    //creat var literal arrays
    for (v = 1; v <= num_vars; ++v)
    {
        var_lit[v] = new lit[var_lit_count[v] + 1];
        var_lit_count[v] = 0; //reset to 0, for build up the array
    }
    //scan all clauses to build up var literal arrays
    for (int i = 0; i < num_clauses; ++i)
    {
        for (int j = 0; j < clause_lit_count[i]; ++j)
        {
            v = clause_lit[i][j].var_num;
            var_lit[v][var_lit_count[v]] = clause_lit[i][j];
            ++var_lit_count[v];
        }
    }

    cout << "c before build neighbor" << endl;

    for (v = 1; v <= num_vars; ++v)
        var_lit[v][var_lit_count[v]].clause_num = -1;

    cout << "c build instime is " << get_runtime() << endl;

    if (get_runtime() > 1.0 || num_clauses > 10000000)
        if_using_neighbor = false;
    else
    {
        cout << "c using neighbor " << endl;
        if_using_neighbor = true;
        build_neighbor_relation();
    }

    best_soln_feasible = 0;
    opt_unsat_weight = total_soft_weight + 1;
}

void Satlike::init(vector<int> &init_solution)
{
    soft_large_weight_clauses_count = 0;
    //Initialize clause information
    for (int c = 0; c < num_clauses; c++)
    {
        already_in_soft_large_weight_stack[c] = 0;
        clause_selected_count[c] = 0;

        if (org_clause_weight[c] == (long long)top_clause_weight)
            clause_weight[c] = 1;
        else
        {
            clause_weight[c] = org_clause_weight[c];
            if (clause_weight[c] > 1 && already_in_soft_large_weight_stack[c] == 0)
            {
                already_in_soft_large_weight_stack[c] = 1;
                soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
            }
        }
    }

    //init solution
    if (init_solution.size() == 0)
    {
        for (int v = 1; v <= num_vars; v++)
        {
            cur_soln[v] = rand() % 2;
            time_stamp[v] = 0;
            unsat_app_count[v] = 0;
        }
    }
    else
    {
        for (int v = 1; v <= num_vars; v++)
        {
            cur_soln[v] = init_solution[v];
            if (cur_soln[v] != 0 && cur_soln[v] != 1)
                cur_soln[v] = rand() % 2;
            time_stamp[v] = 0;
            unsat_app_count[v] = 0;
        }
    }

    //init stacks
    hard_unsat_nb = 0;
    soft_unsat_weight = 0;
    hardunsat_stack_fill_pointer = 0;
    softunsat_stack_fill_pointer = 0;
    unsatvar_stack_fill_pointer = 0;
    large_weight_clauses_count = 0;

    /* figure out sat_count, sat_var and init unsat_stack */
    for (int c = 0; c < num_clauses; ++c)
    {
        sat_count[c] = 0;
        for (int j = 0; j < clause_lit_count[c]; ++j)
        {
            if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
            {
                sat_count[c]++;
                sat_var[c] = clause_lit[c][j].var_num;
            }
        }
        if (sat_count[c] == 0)
        {
            unsat(c);
        }
    }

    /*figure out score*/
    for (int v = 1; v <= num_vars; v++)
    {
        score[v] = 0;
        for (int i = 0; i < var_lit_count[v]; ++i)
        {
            int c = var_lit[v][i].clause_num;
            if (sat_count[c] == 0)
                score[v] += clause_weight[c];
            else if (sat_count[c] == 1 && var_lit[v][i].sense == cur_soln[v])
                score[v] -= clause_weight[c];
        }
    }

    //init goodvars stack
    goodvar_stack_fill_pointer = 0;
    for (int v = 1; v <= num_vars; v++)
    {
        if (score[v] > 0)
        {
            already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
            mypush(v, goodvar_stack);
        }
        else
            already_in_goodvar_stack[v] = -1;
    }
}

void Satlike::increase_weights()
{
    int i, c, v;
    for (i = 0; i < hardunsat_stack_fill_pointer; ++i)
    {
        c = hardunsat_stack[i];
        clause_weight[c] += h_inc;

        if (clause_weight[c] == (h_inc + 1))
            large_weight_clauses[large_weight_clauses_count++] = c;

        for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
        {
            score[v] += h_inc;
            if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
            {
                already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v, goodvar_stack);
            }
        }
    }
    for (i = 0; i < softunsat_stack_fill_pointer; ++i)
    {
        c = softunsat_stack[i];
        if (clause_weight[c] > softclause_weight_threshold)
            continue;
        else
            clause_weight[c]++;

        if (clause_weight[c] > 1 && already_in_soft_large_weight_stack[c] == 0)
        {
            already_in_soft_large_weight_stack[c] = 1;
            soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
        }
        for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
        {
            score[v]++;
            if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
            {
                already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v, goodvar_stack);
            }
        }
    }
}

void Satlike::smooth_weights()
{
    int i, clause, v;

    for (i = 0; i < large_weight_clauses_count; i++)
    {
        clause = large_weight_clauses[i];
        if (sat_count[clause] > 0)
        {
            clause_weight[clause] -= h_inc;

            if (clause_weight[clause] == 1)
            {
                large_weight_clauses[i] = large_weight_clauses[--large_weight_clauses_count];
                i--;
            }
            if (sat_count[clause] == 1)
            {
                v = sat_var[clause];
                score[v] += h_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }

    for (i = 0; i < soft_large_weight_clauses_count; i++)
    {
        clause = soft_large_weight_clauses[i];
        if (sat_count[clause] > 0)
        {
            clause_weight[clause]--;
            if (clause_weight[clause] == 1 && already_in_soft_large_weight_stack[clause] == 1)
            {
                already_in_soft_large_weight_stack[clause] = 0;
                soft_large_weight_clauses[i] = soft_large_weight_clauses[--soft_large_weight_clauses_count];
                i--;
            }
            if (sat_count[clause] == 1)
            {
                v = sat_var[clause];
                score[v]++;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
}

void Satlike::update_clause_weights()
{
    if (((rand() % MY_RAND_MAX_INT) * BASIC_SCALE) < smooth_probability && large_weight_clauses_count > large_clause_count_threshold)
    {
        smooth_weights();
    }
    else
    {
        increase_weights();
    }
}

int Satlike::pick_var()
{
    int i, v;
    int best_var;

    if (goodvar_stack_fill_pointer > 0)
    {
        if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rdprob)
            return goodvar_stack[rand() % goodvar_stack_fill_pointer];

        if (goodvar_stack_fill_pointer < hd_count_threshold)
        {
            best_var = goodvar_stack[0];
            for (i = 1; i < goodvar_stack_fill_pointer; ++i)
            {
                v = goodvar_stack[i];
                if (score[v] > score[best_var])
                    best_var = v;
                else if (score[v] == score[best_var])
                {
                    if (time_stamp[v] < time_stamp[best_var])
                        best_var = v;
                }
            }
            return best_var;
        }
        else
        {
            best_var = goodvar_stack[rand() % goodvar_stack_fill_pointer];
            for (i = 1; i < hd_count_threshold; ++i)
            {
                v = goodvar_stack[rand() % goodvar_stack_fill_pointer];
                if (score[v] > score[best_var])
                    best_var = v;
                else if (score[v] == score[best_var])
                {
                    if (time_stamp[v] < time_stamp[best_var])
                        best_var = v;
                }
            }
            return best_var;
        }
    }

    update_clause_weights();

    int sel_c;
    lit *p;

    if (hardunsat_stack_fill_pointer > 0)
    {
        sel_c = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
    }
    else
    {
        sel_c = softunsat_stack[rand() % softunsat_stack_fill_pointer];
    }
    if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rwprob)
        return clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;

    best_var = clause_lit[sel_c][0].var_num;
    p = clause_lit[sel_c];
    for (p++; (v = p->var_num) != 0; p++)
    {
        if (score[v] > score[best_var])
            best_var = v;
        else if (score[v] == score[best_var])
        {
            if (time_stamp[v] < time_stamp[best_var])
                best_var = v;
        }
    }

    return best_var;
}

void Satlike::local_search(vector<int> &init_solution)
{
    cout << "c changing to SATLike solver! " << endl;
    settings();
    for (tries = 1; tries < max_tries; ++tries)
    {
        init(init_solution);
        for (step = 1; step < max_flips; ++step)
        {
			int counter=0;
			if (hard_unsat_nb == 0) ++counter;
			if (soft_unsat_weight == 0) ++counter;
			if (opt_unsat_weight == 0) ++counter;
			if (best_soln_feasible == 0) ++counter;
            if (hard_unsat_nb == 0 && (soft_unsat_weight < opt_unsat_weight || best_soln_feasible == 0))
            {
                max_flips = step + max_non_improve_flip;
                if (soft_unsat_weight < opt_unsat_weight)
                {
                    best_soln_feasible = 1;
                    opt_unsat_weight = soft_unsat_weight;
                    opt_time = get_runtime();
                    for (int v = 1; v <= num_vars; ++v)
                        best_soln[v] = cur_soln[v];
                }
                if (opt_unsat_weight == 0)
                    return;
            }
            int flipvar = pick_var();
            flip(flipvar);
            time_stamp[flipvar] = step;
        }
    }
}

void Satlike::update_goodvarstack1(int flipvar)
{
    int v;
    //remove the vars no longer goodvar in goodvar stack
    for (int index = goodvar_stack_fill_pointer - 1; index >= 0; index--)
    {
        v = goodvar_stack[index];
        if (score[v] <= 0)
        {
            goodvar_stack[index] = mypop(goodvar_stack);
            already_in_goodvar_stack[goodvar_stack[index]] = index;
            already_in_goodvar_stack[v] = -1;
        }
    }
    if (if_using_neighbor)
    {
        //add goodvar
        for (int i = 0; i < var_neighbor_count[flipvar]; ++i)
        {
            v = var_neighbor[flipvar][i];
            if (score[v] > 0)
            {
                if (already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
    else
    {
        lit *clause_c;
        int c;
        for (int i = 0; i < var_lit_count[flipvar]; ++i)
        {
            c = var_lit[flipvar][i].clause_num;
            clause_c = clause_lit[c];
            for (lit *p = clause_c; (v = p->var_num) != 0; p++)
            {
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }

        if (already_in_goodvar_stack[flipvar] != 0 && score[flipvar] > 0)
        {
            if (goodvar_stack[already_in_goodvar_stack[flipvar]] == flipvar)
            {
                int tem_v = mypop(goodvar_stack);
                goodvar_stack[already_in_goodvar_stack[flipvar]] = tem_v;
                already_in_goodvar_stack[tem_v] = already_in_goodvar_stack[flipvar];
                already_in_goodvar_stack[flipvar] = -1;
            }
        }
    }
}
void Satlike::update_goodvarstack2(int flipvar)
{
    if (score[flipvar] > 0 && already_in_goodvar_stack[flipvar] == -1)
    {
        already_in_goodvar_stack[flipvar] = goodvar_stack_fill_pointer;
        mypush(flipvar, goodvar_stack);
    }
    else if (score[flipvar] <= 0 && already_in_goodvar_stack[flipvar] != -1)
    {
        int index = already_in_goodvar_stack[flipvar];
        int last_v = mypop(goodvar_stack);
        goodvar_stack[index] = last_v;
        already_in_goodvar_stack[last_v] = index;
        already_in_goodvar_stack[flipvar] = -1;
    }
    int i, v;
    for (i = 0; i < var_neighbor_count[flipvar]; ++i)
    {
        v = var_neighbor[flipvar][i];
        if (score[v] > 0)
        {
            if (already_in_goodvar_stack[v] == -1)
            {
                already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v, goodvar_stack);
            }
        }
        else if (already_in_goodvar_stack[v] != -1)
        {
            int index = already_in_goodvar_stack[v];
            int last_v = mypop(goodvar_stack);
            goodvar_stack[index] = last_v;
            already_in_goodvar_stack[last_v] = index;
            already_in_goodvar_stack[v] = -1;
        }
    }
}

void Satlike::flip(int flipvar)
{
    int i, v, c;
    lit *clause_c;

    int org_flipvar_score = score[flipvar];
    cur_soln[flipvar] = 1 - cur_soln[flipvar];
    for (i = 0; i < var_lit_count[flipvar]; ++i)
    {
		
        c = var_lit[flipvar][i].clause_num;
        clause_c = clause_lit[c];
        if (cur_soln[flipvar] == var_lit[flipvar][i].sense)
        {
            ++sat_count[c];
            if (sat_count[c] == 2) //sat_count from 1 to 2
            {
                score[sat_var[c]] += clause_weight[c];
            }
            else if (sat_count[c] == 1) // sat_count from 0 to 1
            {
                sat_var[c] = flipvar; //record the only true lit's var
                for (lit *p = clause_c; (v = p->var_num) != 0; p++)
                {
                    score[v] -= clause_weight[c];
                }
                sat(c);
            }
        }
        else // cur_soln[flipvar] != cur_lit.sense
        {
            --sat_count[c];
            if (sat_count[c] == 1) //sat_count from 2 to 1
            {
                for (lit *p = clause_c; (v = p->var_num) != 0; p++)
                {
                    if (p->sense == cur_soln[v])
                    {
                        score[v] -= clause_weight[c];
                        sat_var[c] = v;
                        break;
                    }
                }
            }
            else if (sat_count[c] == 0) //sat_count from 1 to 0
            {
                for (lit *p = clause_c; (v = p->var_num) != 0; p++)
                {
                    score[v] += clause_weight[c];
                }
                unsat(c);
            } //end else if
        }     //end else
    }

    //update information of flipvar
    score[flipvar] = -org_flipvar_score;
    update_goodvarstack1(flipvar);
}

void Satlike::print_best_solution()
{
    if (best_soln_feasible == 0)
        return;

    printf("v");
    for (int i = 1; i <= num_vars; i++)
    {
        printf(" ");
        if (best_soln[i] == 0)
            printf("-");
        printf("%d", i);
    }
    printf("\n");
}

bool Satlike::verify_sol()
{
    int c, j, flag;
    long long verify_unsat_weight = 0;

    for (c = 0; c < num_clauses; ++c)
    {
        flag = 0;
        for (j = 0; j < clause_lit_count[c]; ++j)
        {
            if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
            {
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            if (org_clause_weight[c] == (long long)top_clause_weight) //verify hard clauses
            {
                //output the clause unsatisfied by the solution
                cout << "c Error: hard clause " << c << " is not satisfied" << endl;

                cout << "c ";
                for (j = 0; j < clause_lit_count[c]; ++j)
                {
                    if (clause_lit[c][j].sense == 0)
                        cout << "-";
                    cout << clause_lit[c][j].var_num << " ";
                }
                cout << endl;
                cout << "c ";
                for (j = 0; j < clause_lit_count[c]; ++j)
                    cout << cur_soln[clause_lit[c][j].var_num] << " ";
                cout << endl;
                return 0;
            }
            else
            {
                verify_unsat_weight += org_clause_weight[c];
            }
        }
    }

    if (verify_unsat_weight == (long long)opt_unsat_weight)
    {
        cout << "c yes " << verify_unsat_weight << endl;
    }
    else
    {
        cout << "c Error: find opt=" << opt_unsat_weight << ", but verified opt=" << verify_unsat_weight << endl;
    }
    return 0;
}

void Satlike::simple_print()
{
    if (best_soln_feasible == 1)
    {
        if (verify_sol() == 1)
            cout << opt_unsat_weight << '\t' << opt_time << endl;
        else
            cout << "solution is wrong " << endl;
    }
    else
        cout << -1 << '\t' << -1 << endl;
}

inline void Satlike::unsat(int clause)
{
    if (org_clause_weight[clause] == (long long)top_clause_weight)
    {
        index_in_hardunsat_stack[clause] = hardunsat_stack_fill_pointer;
        mypush(clause, hardunsat_stack);
        hard_unsat_nb++;
    }
    else
    {
        index_in_softunsat_stack[clause] = softunsat_stack_fill_pointer;
        mypush(clause, softunsat_stack);
        soft_unsat_weight += org_clause_weight[clause];
    }
}

inline void Satlike::sat(int clause)
{
    int index, last_unsat_clause;

    if (org_clause_weight[clause] == (long long)top_clause_weight)
    {
        last_unsat_clause = mypop(hardunsat_stack);
        index = index_in_hardunsat_stack[clause];
        hardunsat_stack[index] = last_unsat_clause;
        index_in_hardunsat_stack[last_unsat_clause] = index;

        hard_unsat_nb--;
    }
    else
    {
        last_unsat_clause = mypop(softunsat_stack);
        index = index_in_softunsat_stack[clause];
        softunsat_stack[index] = last_unsat_clause;
        index_in_softunsat_stack[last_unsat_clause] = index;

        soft_unsat_weight -= org_clause_weight[clause];
    }
}

#endif
