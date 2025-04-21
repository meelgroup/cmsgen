/******************************************
Copyright (c) 2016, Mate Soos

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#pragma once

#include <atomic>
#include <vector>
#include <string>
#include "solvertypesmini.h"

namespace CMSGen {
    struct CMSatPrivateData;
    #ifdef _WIN32
    class __declspec(dllexport) SATSolver
    #else
    class SATSolver
    #endif
    {
    public:
        SATSolver(void* config = NULL
        , std::atomic<bool>* interrupt_asap = NULL
        , uint32_t* seed = NULL
        );
        ~SATSolver();

        ////////////////////////////
        // Adding variables and clauses
        ////////////////////////////
        void new_var(); //add a new variable to the solver
        void new_vars(const size_t n); //and many new variables to the solver -- much faster
        unsigned nVars() const; //get number of variables inside the solver
        bool add_clause(const std::vector<Lit>& lits);
        bool add_xor_clause(const std::vector<unsigned>& vars, bool rhs);
        void set_var_weight(Lit lit, double weight);

        ////////////////////////////
        // Solving and simplifying
        // You can call solve() multiple times: incremental mode is supported!
        ////////////////////////////

        lbool solve(const std::vector<Lit>* assumptions = 0, bool only_indep_solution = false); //solve the problem, optionally with assumptions. If only_indep_solution is set, only the independent variables set with set_independent_vars() are returned in the solution
        const std::vector<lbool>& get_model() const; //get model that satisfies the problem. Only makes sense if previous solve()/simplify() call was l_True
        const std::vector<Lit>& get_conflict() const; //get conflict in terms of the assumptions given in case the previous call to solve() was l_False
        bool okay() const; //the problem is still solveable, i.e. the empty clause hasn't been derived

        ////////////////////////////
        // Debug all calls for later replay with --debuglit FILENAME
        ////////////////////////////
        void log_to_file(std::string filename);

        ////////////////////////////
        // Configuration
        // -- Note that nothing else can be changed, only these.
        // -- The main.cpp has access to the internal config, but it changes
        // -- all the time and hence exposing it to the outside world would
        // -- be very brittle.
        ////////////////////////////
        void set_allow_otf_gauss(); //allow on-the-fly gaussian elimination
        void set_max_time(double max_time); //max time to run to on next solve() call
        void set_max_confl(int64_t max_confl); //max conflict to run to on next solve() call
        void set_verbosity(unsigned verbosity = 0); //default is 0, silent
        void set_no_simplify(); //never simplify
        void set_sampling_vars(std::vector<uint32_t>* sampl_vars);
        void set_timeout_all_calls(double secs); //max timeout on all subsequent solve() or simplify

        ////////////////////////////
        // Get generic info
        ////////////////////////////
        static const char* get_version(); //get solver version in string format
        static const char* get_version_sha1(); //get SHA1 version string of the solver
        static const char* get_compilation_env(); //get compilation environment string
        std::string get_text_version_info();  //get printable version and copyright text


        ////////////////////////////
        // Get info about only the last solve() OR simplify() call
        // summed for all threads
        ////////////////////////////
        uint64_t get_last_conflicts(); //get total number of conflicts of last solve() or simplify() call of all threads
        uint64_t get_last_propagations();  //get total number of propagations of last solve() or simplify() call made by all threads
        uint64_t get_last_decisions(); //get total number of decisions of last solve() or simplify() call made by all threads

        ////////////////////////////
        //Get info about total sum of all time of all threads
        ////////////////////////////

        uint64_t get_sum_conflicts(); //get total number of conflicts of all time of all threads
        uint64_t get_sum_propagations();  //get total number of propagations of all time made by all threads
        uint64_t get_sum_decisions(); //get total number of decisions of all time made by all threads

        void print_stats() const; //print solving stats. Call after solve()/simplify()
        void interrupt_asap(); //call this asynchronously, and the solver will try to cleanly abort asap
        void add_in_partial_solving_stats(); //used only by Ctrl+C handler. Ignore.

    private:

        ////////////////////////////
        // Do not bother with this, it's private
        ////////////////////////////
        CMSatPrivateData *data;
    };
}
