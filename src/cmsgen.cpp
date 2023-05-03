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

#include "constants.h"
#include "cmsgen.h"
#include "solver.h"
#include "drat.h"
#include "shareddata.h"
#include <fstream>

#include <thread>
#include <mutex>
#include <atomic>
using std::thread;

#define CACHE_SIZE 10ULL*1000ULL*1000UL
#ifndef LIMITMEM
#define MAX_VARS (1ULL<<28)
#else
#define MAX_VARS 3000
#endif

using namespace CMSGen;

static bool print_thread_start_and_finish = false;

namespace CMSGen {
    struct CMSatPrivateData {
        explicit CMSatPrivateData(std::atomic<bool>* _must_interrupt)
        {
            must_interrupt = _must_interrupt;
            if (must_interrupt == NULL) {
                must_interrupt = new std::atomic<bool>(false);
                must_interrupt_needs_delete = true;
            }
        }
        ~CMSatPrivateData()
        {
            for(Solver* this_s: solvers) {
                delete this_s;
            }
            if (must_interrupt_needs_delete) {
                delete must_interrupt;
            }

            delete shared_data;
        }
        CMSatPrivateData(const CMSatPrivateData&) = delete;
        CMSatPrivateData& operator=(const CMSatPrivateData&) = delete;

        vector<Solver*> solvers;
        SharedData *shared_data = NULL;
        int which_solved = 0;
        std::atomic<bool>* must_interrupt;
        bool must_interrupt_needs_delete = false;
        bool okay = true;
        int sql = 0;
        double timeout = std::numeric_limits<double>::max();
        bool interrupted = false;

        //variables and clauses added/to add
        unsigned cls = 0;
        unsigned vars_to_add = 0;
        vector<Lit> cls_lits;

        //For single call setup
        uint32_t num_solve_simplify_calls = 0;

        //stats
        uint64_t previous_sum_conflicts = 0;
        uint64_t previous_sum_propagations = 0;
        uint64_t previous_sum_decisions = 0;
        vector<double> cpu_times;
    };
}

struct DataForThread
{
    explicit DataForThread(CMSatPrivateData* data, const vector<Lit>* _assumptions = NULL) :
        solvers(data->solvers)
        , cpu_times(data->cpu_times)
        , lits_to_add(&(data->cls_lits))
        , vars_to_add(data->vars_to_add)
        , assumptions(_assumptions)
        , update_mutex(new std::mutex)
        , which_solved(&(data->which_solved))
        , ret(new lbool(l_Undef))
    {
    }

    ~DataForThread()
    {
        delete update_mutex;
        delete ret;
    }
    vector<Solver*>& solvers;
    vector<double>& cpu_times;
    vector<Lit> *lits_to_add;
    uint32_t vars_to_add;
    const vector<Lit> *assumptions;
    std::mutex* update_mutex;
    int *which_solved;
    lbool* ret;
};

DLL_PUBLIC SATSolver::SATSolver(
    void* config
    , std::atomic<bool>* interrupt_asap
    )
{
    data = new CMSatPrivateData(interrupt_asap);

    if (config && ((SolverConf*) config)->verbosity) {
        //NOT SAFE
        //yes -- this system will use a lock, but the solver itself won't(!)
        //so things will get mangled and printed wrongly
        //print_thread_start_and_finish = true;
    }

    data->solvers.push_back(new Solver((SolverConf*) config, data->must_interrupt));
    data->cpu_times.push_back(0.0);
}

DLL_PUBLIC SATSolver::~SATSolver()
{
    delete data;
}

struct OneThreadAddCls
{
    OneThreadAddCls(DataForThread& _data_for_thread, size_t _tid) :
        data_for_thread(_data_for_thread)
        , tid(_tid)
    {
    }

    void operator()()
    {
        Solver& solver = *data_for_thread.solvers[tid];
        solver.new_external_vars(data_for_thread.vars_to_add);

        vector<Lit> lits;
        vector<uint32_t> vars;
        bool ret = true;
        size_t at = 0;
        const vector<Lit>& orig_lits = (*data_for_thread.lits_to_add);
        const size_t size = orig_lits.size();
        while(at < size && ret) {
            if (orig_lits[at] == lit_Undef) {
                lits.clear();
                at++;
                for(; at < size
                    && orig_lits[at] != lit_Undef
                    && orig_lits[at] != lit_Error
                    ; at++
                ) {
                    lits.push_back(orig_lits[at]);
                }
                ret = solver.add_clause_outer(lits);
            } else {
                vars.clear();
                at++;
                bool rhs = orig_lits[at].sign();
                at++;
                for(; at < size
                    && orig_lits[at] != lit_Undef
                    && orig_lits[at] != lit_Error
                    ; at++
                ) {
                    vars.push_back(orig_lits[at].var());
                }
                ret = solver.add_xor_clause_outer(vars, rhs);
            }
        }

        if (!ret) {
            data_for_thread.update_mutex->lock();
            *data_for_thread.ret = l_False;
            data_for_thread.update_mutex->unlock();
        }
    }

    DataForThread& data_for_thread;
    const size_t tid;
};

static bool actually_add_clauses_to_threads(CMSatPrivateData* data)
{
    DataForThread data_for_thread(data);
    std::vector<std::thread> thds;
    for(size_t i = 0; i < data->solvers.size(); i++) {
        thds.push_back(thread(OneThreadAddCls(data_for_thread, i)));
    }
    for(std::thread& thread : thds){
        thread.join();
    }
    bool ret = (*data_for_thread.ret == l_True);

    //clear what has been added
    data->cls_lits.clear();
    data->vars_to_add = 0;

    return ret;
}

DLL_PUBLIC void SATSolver::set_max_time(double max_time)
{
  for (size_t i = 0; i < data->solvers.size(); ++i) {
    Solver& s = *data->solvers[i];
    if (max_time >= 0) {
      // the main loop in solver.cpp is checks `maxTime`
      // against `cpuTime`, so we specify `s.conf.maxTime`
      // as an offset from `cpuTime`.
      s.conf.maxTime = cpuTime() + max_time;

      //don't allow for overflow
      if (s.conf.maxTime < max_time)
          s.conf.maxTime = max_time;
    }
  }
}

DLL_PUBLIC void SATSolver::set_max_confl(int64_t max_confl)
{
  for (size_t i = 0; i < data->solvers.size(); ++i) {
    Solver& s = *data->solvers[i];
    if (max_confl >= 0) {
      s.conf.max_confl = s.get_stats().conflStats.numConflicts + max_confl;

      //don't allow for overflow
      if (s.conf.max_confl < max_confl)
          s.conf.max_confl = max_confl;
    }
  }
}

DLL_PUBLIC void SATSolver::set_no_simplify()
{
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        s.conf.doRenumberVars = false;
        s.conf.simplify_at_startup = false;
        s.conf.simplify_at_every_startup = false;
        s.conf.full_simplify_at_startup = false;
        s.conf.perform_occur_based_simp = false;
        s.conf.do_simplify_problem = false;
    }
}

DLL_PUBLIC void SATSolver::set_allow_otf_gauss()
{
    #ifndef USE_GAUSS
    std::cerr << "ERROR: CryptoMiniSat was not compiled with GAUSS" << endl;
    exit(-1);
    #else
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        //s.conf.reconfigure_at = 0;
        //s.conf.reconfigure_val = 15;
        s.conf.gaussconf.max_num_matrices = 10;
        s.conf.gaussconf.autodisable = false;
        s.conf.allow_elim_xor_vars = false;
    }
    #endif
}

DLL_PUBLIC void SATSolver::set_no_simplify_at_startup()
{
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        s.conf.simplify_at_startup = false;
    }
}

DLL_PUBLIC void SATSolver::set_no_equivalent_lit_replacement()
{
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        s.conf.doFindAndReplaceEqLits = false;
    }
}

DLL_PUBLIC void SATSolver::set_sampling_vars(vector<uint32_t>* sampl_vars)
{
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        s.conf.sampling_vars = sampl_vars;
    }
}

DLL_PUBLIC void SATSolver::set_verbosity(unsigned verbosity)
{
    if (data->solvers.empty())
        return;

    Solver& s = *data->solvers[0];
    s.conf.verbosity = verbosity;
}

DLL_PUBLIC void SATSolver::set_timeout_all_calls(double timeout)
{
    data->timeout = timeout;
}

DLL_PUBLIC bool SATSolver::add_clause(const vector< Lit >& lits)
{
    bool ret = true;
    if (data->solvers.size() > 1) {
        if (data->cls_lits.size() + lits.size() + 1 > CACHE_SIZE) {
            ret = actually_add_clauses_to_threads(data);
        }

        data->cls_lits.push_back(lit_Undef);
        for(Lit lit: lits) {
            data->cls_lits.push_back(lit);
        }
    } else {
        data->solvers[0]->new_vars(data->vars_to_add);
        data->vars_to_add = 0;

        ret = data->solvers[0]->add_clause_outer(lits);
        data->cls++;
    }

    return ret;
}

DLL_PUBLIC bool SATSolver::add_xor_clause(const std::vector<unsigned>& vars, bool rhs)
{
    bool ret = true;
    if (data->solvers.size() > 1) {
        if (data->cls_lits.size() + vars.size() + 1 > CACHE_SIZE) {
            ret = actually_add_clauses_to_threads(data);
        }

        data->cls_lits.push_back(lit_Error);
        data->cls_lits.push_back(Lit(0, rhs));
        for(uint32_t var: vars) {
            data->cls_lits.push_back(Lit(var, false));
        }
    } else {
        data->solvers[0]->new_vars(data->vars_to_add);
        data->vars_to_add = 0;

        ret = data->solvers[0]->add_xor_clause_outer(vars, rhs);
        data->cls++;
    }

    return ret;
}

struct OneThreadCalc
{
    OneThreadCalc(
        DataForThread& _data_for_thread,
        size_t _tid,
        bool _solve,
        bool _only_sampling_solution
    ) :
        data_for_thread(_data_for_thread)
        , tid(_tid)
        , solve(_solve)
        , only_sampling_solution(_only_sampling_solution)
    {}

    void operator()()
    {
        if (print_thread_start_and_finish) {
            start_time = cpuTime();
            //data_for_thread.update_mutex->lock();
            //cout << "c Starting thread " << tid << endl;
            //data_for_thread.update_mutex->unlock();
        }

        OneThreadAddCls cls_adder(data_for_thread, tid);
        cls_adder();
        lbool ret;
        if (solve) {
            ret = data_for_thread.solvers[tid]->solve_with_assumptions(data_for_thread.assumptions, only_sampling_solution);
        } else {
            ret = data_for_thread.solvers[tid]->simplify_with_assumptions(data_for_thread.assumptions);
        }

        data_for_thread.cpu_times[tid] = cpuTime();
        if (print_thread_start_and_finish) {
            data_for_thread.update_mutex->lock();
            ios::fmtflags f(cout.flags());
            cout << "c Finished thread " << tid << " with result: " << ret
            << " T-diff: " << std::fixed << std::setprecision(2)
            << (data_for_thread.cpu_times[tid]-start_time)
            << endl;
            cout.flags(f);
            data_for_thread.update_mutex->unlock();
        }


        if (ret != l_Undef) {
            data_for_thread.update_mutex->lock();
            *data_for_thread.which_solved = tid;
            *data_for_thread.ret = ret;
            //will interrupt all of them
            data_for_thread.solvers[0]->set_must_interrupt_asap();
            data_for_thread.update_mutex->unlock();
        }
    }

    DataForThread& data_for_thread;
    const size_t tid;
    double start_time;
    bool solve;
    bool only_sampling_solution;
};

lbool calc(
    const vector< Lit >* assumptions,
    bool solve, CMSatPrivateData *data,
    bool only_sampling_solution = false
) {
    //Reset the interrupt signal if it was set
    data->must_interrupt->store(false, std::memory_order_relaxed);

    //Set timeout information
    if (data->timeout != std::numeric_limits<double>::max()) {
        for (size_t i = 0; i < data->solvers.size(); ++i) {
            Solver& s = *data->solvers[i];
            s.conf.maxTime = cpuTime() + data->timeout;
        }
    }

    if (data->solvers.size() > 1 && data->sql > 0) {
        std::cerr
        << "Multithreaded solving and SQL cannot be specified at the same time"
        << endl;
        exit(-1);
    }

    if (data->solvers.size() == 1) {
        data->solvers[0]->new_vars(data->vars_to_add);
        data->vars_to_add = 0;

        lbool ret ;
        if (solve) {
            ret = data->solvers[0]->solve_with_assumptions(assumptions, only_sampling_solution);
        } else {
            ret = data->solvers[0]->simplify_with_assumptions(assumptions);
        }
        data->okay = data->solvers[0]->okay();
        data->cpu_times[0] = cpuTime();
        return ret;
    }

    //Multi-thread from now on.
    DataForThread data_for_thread(data, assumptions);
    std::vector<std::thread> thds;
    for(size_t i = 0
        ; i < data->solvers.size()
        ; i++
    ) {
        thds.push_back(thread(OneThreadCalc(data_for_thread, i, solve, only_sampling_solution)));
    }
    for(std::thread& thread : thds){
        thread.join();
    }
    lbool real_ret = *data_for_thread.ret;

    //This does it for all of them, there is only one must-interrupt
    data_for_thread.solvers[0]->unset_must_interrupt_asap();

    //clear what has been added
    data->cls_lits.clear();
    data->vars_to_add = 0;
    data->okay = data->solvers[*data_for_thread.which_solved]->okay();
    return real_ret;
}

DLL_PUBLIC lbool SATSolver::solve(const vector< Lit >* assumptions, bool only_sampling_solution)
{
    data->num_solve_simplify_calls++;

    //set information data (props, confl, dec)
    data->previous_sum_conflicts = get_sum_conflicts();
    data->previous_sum_propagations = get_sum_propagations();
    data->previous_sum_decisions = get_sum_decisions();

    return calc(assumptions, true, data, only_sampling_solution);
}

DLL_PUBLIC lbool SATSolver::simplify(const vector< Lit >* assumptions)
{
    data->num_solve_simplify_calls++;

    //set information data (props, confl, dec)
    data->previous_sum_conflicts = get_sum_conflicts();
    data->previous_sum_propagations = get_sum_propagations();
    data->previous_sum_decisions = get_sum_decisions();

    return calc(assumptions, false, data);
}

DLL_PUBLIC const vector< lbool >& SATSolver::get_model() const
{
    return data->solvers[data->which_solved]->get_model();
}

DLL_PUBLIC const std::vector<Lit>& SATSolver::get_conflict() const
{

    return data->solvers[data->which_solved]->get_final_conflict();
}

DLL_PUBLIC uint32_t SATSolver::nVars() const
{
    return data->solvers[0]->nVarsOutside() + data->vars_to_add;
}

DLL_PUBLIC void SATSolver::new_var()
{
    new_vars(1);
}

DLL_PUBLIC void SATSolver::new_vars(const size_t n)
{
    if (n >= MAX_VARS
        || (data->vars_to_add + n) >= MAX_VARS
    ) {
        throw CMSGen::TooManyVarsError();
    }

    data->vars_to_add += n;
}

DLL_PUBLIC const char* SATSolver::get_version_sha1()
{
    return Solver::get_version_sha1();
}

DLL_PUBLIC const char* SATSolver::get_version()
{
    return Solver::get_version_tag();
}

DLL_PUBLIC const char* SATSolver::get_compilation_env()
{
    return Solver::get_compilation_env();
}

std::string SATSolver::get_text_version_info()
{
    std::stringstream ss;
    ss << "c CMSGen Copyright Mate Soos (soos.mate@gmail.com)" << endl;
    ss << "c CMSGen SHA revision " << get_version_sha1() << endl;
    ss << "c CMSGen is MIT licensed" << endl;
    ss << "c CMSGen compilation env " << get_compilation_env() << endl;
    #ifdef __GNUC__
    ss << "c CMSGen compiled with gcc version " << __VERSION__ << endl;
    #else
    ss << "c CMSGen compiled with non-gcc compiler" << endl;
    #endif

    return ss.str();
}

DLL_PUBLIC void SATSolver::print_stats() const
{
    double cpu_time_total = cpuTimeTotal();

    double cpu_time;
    if (data->interrupted) {
        //cannot know, we have in fact no idea how much time passed...
        //we have to guess. Shitty guess comes here... :S
        cpu_time = cpuTimeTotal()/(double)data->solvers.size();
    } else {
        cpu_time = data->cpu_times[data->which_solved];
    }

    //If only one thread, then don't confuse the user. The difference
    //is minimal.
    if (data->solvers.size() == 1) {
        cpu_time = cpu_time_total;
    }

    data->solvers[data->which_solved]->print_stats(cpu_time, cpu_time_total);
}

DLL_PUBLIC void SATSolver::interrupt_asap()
{
    data->must_interrupt->store(true, std::memory_order_relaxed);
}

void DLL_PUBLIC SATSolver::add_in_partial_solving_stats()
{
    data->solvers[data->which_solved]->add_in_partial_solving_stats();
    data->interrupted = true;
}

DLL_PUBLIC std::vector<Lit> SATSolver::get_zero_assigned_lits() const
{
    return data->solvers[data->which_solved]->get_zero_assigned_lits();
}

DLL_PUBLIC bool SATSolver::okay() const
{
    return data->okay;
}

DLL_PUBLIC std::vector<std::pair<Lit, Lit> > SATSolver::get_all_binary_xors() const
{
    return data->solvers[0]->get_all_binary_xors();
}

DLL_PUBLIC vector<std::pair<vector<uint32_t>, bool> >
SATSolver::get_recovered_xors(bool elongate) const
{
    vector<std::pair<vector<uint32_t>, bool> > ret;
    Solver& s = *data->solvers[0];

    std::pair<vector<uint32_t>, bool> tmp;
    vector<Xor> xors = s.get_recovered_xors(elongate);
    for(const auto& x: xors) {
        tmp.first = x.get_vars();
        tmp.second = x.rhs;
        ret.push_back(tmp);
    }
    return ret;
}

DLL_PUBLIC uint64_t SATSolver::get_sum_conflicts()
{
    uint64_t conlf = 0;
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        conlf += s.sumConflicts;
    }
    return conlf;
}

DLL_PUBLIC uint64_t SATSolver::get_sum_propagations()
{
    uint64_t props = 0;
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        props += s.sumPropStats.propagations;
    }
    return props;
}

DLL_PUBLIC uint64_t SATSolver::get_sum_decisions()
{
    uint64_t dec = 0;
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        dec += s.sumSearchStats.decisions;
    }
    return dec;
}

DLL_PUBLIC uint64_t SATSolver::get_last_conflicts()
{
    return get_sum_conflicts() - data->previous_sum_conflicts;
}

DLL_PUBLIC uint64_t SATSolver::get_last_propagations()
{
    return get_sum_propagations() - data->previous_sum_propagations;
}

DLL_PUBLIC uint64_t SATSolver::get_last_decisions()
{
    return get_sum_decisions() - data->previous_sum_decisions;
}

DLL_PUBLIC const std::vector<Lit>& SATSolver::get_decisions_reaching_model() const
{
    if (!get_decision_reaching_valid()) {
        cout << "ERROR: you called get_decisions_reaching_model() but it's not a valid decision set!" << endl;
        exit(-1);
    }
    return data->solvers[data->which_solved]->get_decisions_reaching_model();
}

DLL_PUBLIC void SATSolver::set_need_decisions_reaching()
{
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        s.conf.need_decisions_reaching = true;
    }
}

DLL_PUBLIC bool SATSolver::get_decision_reaching_valid() const
{
    return data->solvers[data->which_solved]->get_decision_reaching_valid();
}

DLL_PUBLIC void SATSolver::set_var_weight(Lit lit, double weight)
{
    actually_add_clauses_to_threads(data);
    for (size_t i = 0; i < data->solvers.size(); ++i) {
        Solver& s = *data->solvers[i];
        s.set_var_weight(lit, weight);
    }
}
