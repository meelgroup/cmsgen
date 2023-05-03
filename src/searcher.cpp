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

#include "searcher.h"
#include "occsimplifier.h"
#include "time_mem.h"
#include "solver.h"
#include <iomanip>
#include "varreplacer.h"
#include "clausecleaner.h"
#include "propbyforgraph.h"
#include <algorithm>
#include <sstream>
#include <cstddef>
#include <cmath>
#include <ratio>
#include "reducedb.h"
#include "watchalgos.h"
#include "hasher.h"
#include "solverconf.h"
#include "distillerlong.h"
#include "xorfinder.h"
#include "matrixfinder.h"
#ifdef USE_GAUSS
#include "gaussian.h"
#endif

//#define DEBUG_RESOLV
//#define VERBOSE_DEBUG

using namespace CMSGen;
using std::cout;
using std::endl;

//#define VERBOSE_DEBUG_GEN_CONFL_DOT

#ifdef VERBOSE_DEBUG
#define VERBOSE_DEBUG_GEN_CONFL_DOT
#endif

/**
@brief Sets a sane default config and allocates handler classes
*/
Searcher::Searcher(const SolverConf *_conf, Solver* _solver, std::atomic<bool>* _must_interrupt_inter) :
        HyperEngine(
            _conf
            , _solver
            , _must_interrupt_inter
        )
        , solver(_solver)
        , cla_inc(1)
{
    var_decay_vsids = conf.var_decay_vsids_start;
    step_size = solver->conf.orig_step_size;

    var_inc_vsids = conf.var_inc_vsids_start;
    more_red_minim_limit_binary_actual = conf.more_red_minim_limit_binary;
    more_red_minim_limit_cache_actual = conf.more_red_minim_limit_cache;
    mtrand.seed(conf.origSeed);
    hist.setSize(conf.shortTermHistorySize, conf.blocking_restart_trail_hist_length);
    cur_max_temp_red_lev2_cls = conf.max_temp_lev2_learnt_clauses;
}

Searcher::~Searcher()
{
    #ifdef USE_GAUSS
    clear_gauss_matrices();
    #endif
}

void Searcher::new_var(const bool bva, const uint32_t orig_outer)
{
    PropEngine::new_var(bva, orig_outer);

    var_act_vsids.push_back(0);
    var_act_maple.push_back(0);
    insert_var_order_all((int)nVars()-1);
}

void Searcher::new_vars(size_t n)
{
    PropEngine::new_vars(n);

    var_act_vsids.insert(var_act_vsids.end(), n, 0);
    var_act_maple.insert(var_act_maple.end(), n, 0);
    for(int i = n-1; i >= 0; i--) {
        insert_var_order_all((int)nVars()-i-1);
    }
}

void Searcher::save_on_var_memory()
{
    PropEngine::save_on_var_memory();

    var_act_vsids.resize(nVars());
    var_act_maple.resize(nVars());

    var_act_vsids.shrink_to_fit();
    var_act_maple.shrink_to_fit();

}

void Searcher::updateVars(
    const vector<uint32_t>& /*outerToInter*/
    , const vector<uint32_t>& interToOuter
) {
    updateArray(var_act_vsids, interToOuter);
    updateArray(var_act_maple, interToOuter);
}

template<bool update_bogoprops>
inline void Searcher::add_lit_to_learnt(
    const Lit lit
) {
    const uint32_t var = lit.var();
    assert(varData[var].removed == Removed::none);

    //If var is at level 0, don't do anything with it, just skip
    if (seen[var] || varData[var].level == 0) {
        return;
    }
    seen[var] = 1;

    if (!update_bogoprops) {
        bump_vsids_var_act<update_bogoprops>(var, 0.5);
        implied_by_learnts.push_back(var);

        if (conf.doOTFSubsume) {
            tmp_learnt_clause_size++;
            seen2[lit.toInt()] = 1;
            tmp_learnt_clause_abst |= abst_var(lit.var());
        }
    }

    if (varData[var].level >= decisionLevel()) {
        pathC++;
    } else {
        learnt_clause.push_back(lit);
    }
}

inline void Searcher::recursiveConfClauseMin()
{
    uint32_t abstract_level = 0;
    for (size_t i = 1; i < learnt_clause.size(); i++) {
        //(maintain an abstraction of levels involved in conflict)
        abstract_level |= abstractLevel(learnt_clause[i].var());
    }

    size_t i, j;
    for (i = j = 1; i < learnt_clause.size(); i++) {
        if (varData[learnt_clause[i].var()].reason.isNULL()
            || !litRedundant(learnt_clause[i], abstract_level)
        ) {
            learnt_clause[j++] = learnt_clause[i];
        }
    }
    learnt_clause.resize(j);
}

void Searcher::create_otf_subsuming_implicit_clause(const Clause& cl)
{
    OTFClause newCl;
    newCl.size = 0;
    for(const Lit
        *it = cl.begin(), *end = cl.end()
        ; it != end
        ; ++it
    ) {
        if (seen2[it->toInt()]) {
            assert(newCl.size < 3);
            newCl.lits[newCl.size] = *it;
            newCl.size++;
        }
    }
    otf_subsuming_short_cls.push_back(newCl);
    if (conf.verbosity >= 6) {
        cout << "New implicit clause that subsumes a long clause:";
        for(unsigned  i = 0; i < newCl.size; i++) {
            cout
            << newCl.lits[i] << " ";
        }
        cout  << endl;
    }

    if (drat->enabled() || solver->conf.simulate_drat) {
        *drat << add;
        for(unsigned  i = 0; i < newCl.size; i++) {
            *drat << newCl.lits[i];
        }
        *drat << fin;
    }

    stats.otfSubsumed++;
    stats.otfSubsumedImplicit++;
    stats.otfSubsumedRed += cl.red();
    stats.otfSubsumedLitsGained += cl.size() - newCl.size;
}

void Searcher::create_otf_subsuming_long_clause(
    Clause& cl
    , const ClOffset offset
) {
    (*solver->drat) << deldelay << cl << fin;
    solver->detachClause(cl, false);
    stats.otfSubsumed++;
    stats.otfSubsumedLong++;
    stats.otfSubsumedRed += cl.red();
    stats.otfSubsumedLitsGained += cl.size() - tmp_learnt_clause_size;

    size_t i = 0;
    size_t i2 = 0;
    for (; i < cl.size(); i++) {
        if (seen2[cl[i].toInt()]) {
            cl[i2++] = cl[i];
        }
    }
    cl.shrink(i-i2);
    assert(cl.size() == tmp_learnt_clause_size);
    if (conf.verbosity >= 6) {
        cout
        << "New smaller clause OTF:" << cl << endl;
    }
    *drat << add << cl
    << fin << findelay;
    otf_subsuming_long_cls.push_back(offset);
}

void Searcher::check_otf_subsume(const ClOffset offset, Clause& cl)
{
    size_t num_lits_from_cl = 0;
    for (const Lit lit: cl) {
        if (seen2[lit.toInt()]) {
            num_lits_from_cl++;
        }
    }
    if (num_lits_from_cl != tmp_learnt_clause_size)
        return;

    if (num_lits_from_cl <= 2) {
        create_otf_subsuming_implicit_clause(cl);
    } else {
        create_otf_subsuming_long_clause(cl, offset);
    }
}

void Searcher::normalClMinim()
{
    size_t i,j;
    for (i = j = 1; i < learnt_clause.size(); i++) {
        const PropBy& reason = varData[learnt_clause[i].var()].reason;
        size_t size;
        Clause* cl = NULL;
        PropByType type = reason.getType();
        if (type == null_clause_t) {
            learnt_clause[j++] = learnt_clause[i];
            continue;
        }

        switch (type) {
            case clause_t:
                cl = cl_alloc.ptr(reason.get_offset());
                size = cl->size()-1;
                break;

            case binary_t:
                size = 1;
                break;

            default:
                release_assert(false);
                std::exit(-1);
        }

        for (size_t k = 0; k < size; k++) {
            Lit p;
            switch (type) {
                case clause_t:
                    p = (*cl)[k+1];
                    break;

                case binary_t:
                    p = reason.lit2();
                    break;

                default:
                    release_assert(false);
                    std::exit(-1);
            }

            if (!seen[p.var()] && varData[p.var()].level > 0) {
                learnt_clause[j++] = learnt_clause[i];
                break;
            }
        }
    }
    learnt_clause.resize(j);
}

void Searcher::debug_print_resolving_clause(const PropBy confl) const
{
#ifndef DEBUG_RESOLV
    //Avoid unused parameter warning
    (void) confl;
#else
    switch(confl.getType()) {
        case binary_t: {
            cout << "resolv bin: " << confl.lit2() << endl;
            break;
        }

        case clause_t: {
            Clause* cl = cl_alloc.ptr(confl.get_offset());
            cout << "resolv (long): " << *cl << endl;
            break;
        }

        case xor_t: {
            //in the future, we'll have XOR clauses. Not yet.
            assert(false);
            exit(-1);
            break;
        }

        case null_clause_t: {
            assert(false);
            break;
        }
    }
#endif
}

void Searcher::update_clause_glue_from_analysis(Clause* cl)
{
    assert(cl->red());
    const unsigned new_glue = calc_glue(*cl);

    if (new_glue < cl->stats.glue) {
        if (cl->stats.glue <= conf.protect_cl_if_improved_glue_below_this_glue_for_one_turn) {
            cl->stats.ttl = 1;
        }
        cl->stats.glue = new_glue;

        if (cl->stats.locked_for_data_gen) {
            assert(cl->stats.which_red_array == 0);
        } else if (new_glue <= conf.glue_put_lev0_if_below_or_eq
            && cl->stats.which_red_array >= 1
        ) {
            //move to lev0 if very low glue
            cl->stats.which_red_array = 0;
        } else {
            //move to lev1 if low glue
            if (new_glue <= conf.glue_put_lev1_if_below_or_eq
                && solver->conf.glue_put_lev1_if_below_or_eq != 0
            ) {
                cl->stats.which_red_array = 1;
            }
        }
     }
}

template<bool update_bogoprops>
Clause* Searcher::add_literals_from_confl_to_learnt(
    const PropBy confl
    , const Lit p
) {
    #ifdef VERBOSE_DEBUG
    debug_print_resolving_clause(confl);
    #endif
    sumAntecedents++;

    Clause* cl = NULL;
    switch (confl.getType()) {
        case binary_t : {
            sumAntecedentsLits += 2;
            if (confl.isRedStep()) {
                stats.resolvs.binRed++;
            } else {
                stats.resolvs.binIrred++;
            }
            break;
        }

        case clause_t : {
            cl = cl_alloc.ptr(confl.get_offset());
            sumAntecedentsLits += cl->size();
            if (cl->red()) {
                stats.resolvs.longRed++;
            } else {
                stats.resolvs.longIrred++;
            }

            if (!update_bogoprops
                && cl->red()
                && cl->stats.which_red_array != 0
            ) {
                if (conf.update_glues_on_analyze) {
                    update_clause_glue_from_analysis(cl);
                }
                cl->stats.last_touched = sumConflicts;

                //If stats or predictor, bump all because during final
                //we will need this data and during dump when stats is on
                //we also need this data.
                if (cl->stats.which_red_array == 2) {
                    bump_cl_act<update_bogoprops>(cl);
                }
            }

            break;
        }

        case null_clause_t:
        default:
            assert(false && "Error in conflict analysis (otherwise should be UIP)");
    }

    size_t i = 0;
    bool cont = true;
    Lit x = lit_Undef;
    while(cont) {
        switch (confl.getType()) {
            case binary_t:
                if (i == 0) {
                    x = failBinLit;
                } else {
                    x = confl.lit2();
                    cont = false;
                }
                break;

            case clause_t:
                assert(!cl->getRemoved());
                x = (*cl)[i];
                if (i == cl->size()-1) {
                    cont = false;
                }
                break;
            case null_clause_t:
                assert(false);
        }
        if (p == lit_Undef || i > 0) {
            add_lit_to_learnt<update_bogoprops>(x);
        }
        i++;
    }
    return cl;
}

template<bool update_bogoprops>
inline void Searcher::minimize_learnt_clause()
{
    const size_t origSize = learnt_clause.size();

    toClear = learnt_clause;
    if (conf.doRecursiveMinim) {
        recursiveConfClauseMin();
    } else {
        normalClMinim();
    }
    for (const Lit lit: toClear) {
        if (!update_bogoprops
            && conf.doOTFSubsume
        ) {
            seen2[lit.toInt()] = 0;
        }
        seen[lit.var()] = 0;
    }
    toClear.clear();

    stats.recMinCl += ((origSize - learnt_clause.size()) > 0);
    stats.recMinLitRem += origSize - learnt_clause.size();
}

inline void Searcher::minimize_using_permdiff()
{
    if (conf.doMinimRedMore
        && learnt_clause.size() > 1
    ) {
        stats.permDiff_attempt++;
        stats.moreMinimLitsStart += learnt_clause.size();
        watch_based_learnt_minim();

        stats.moreMinimLitsEnd += learnt_clause.size();
    }
}

inline void Searcher::watch_based_learnt_minim()
{
    MYFLAG++;
    const auto& ws  = watches[~learnt_clause[0]];
    uint32_t nb = 0;
    for (const Watched& w: ws) {
        if (w.isBin()) {
            Lit imp = w.lit2();
            if (permDiff[imp.var()] == MYFLAG && value(imp) == l_True) {
                nb++;
                permDiff[imp.var()] = MYFLAG - 1;
            }
        } else {
            break;
        }
    }
    uint32_t l = learnt_clause.size() - 1;
    if (nb > 0) {
        for (uint32_t i = 1; i < learnt_clause.size() - nb; i++) {
            if (permDiff[learnt_clause[i].var()] != MYFLAG) {
                Lit p = learnt_clause[l];
                learnt_clause[l] = learnt_clause[i];
                learnt_clause[i] = p;
                l--;
                i--;
            }
        }
        learnt_clause.resize(learnt_clause.size()-nb);
        stats.permDiff_success++;
        stats.permDiff_rem_lits+=nb;
    }
}

void Searcher::print_fully_minimized_learnt_clause() const
{
    if (conf.verbosity >= 6) {
        cout << "Final clause: " << learnt_clause << endl;
        for (uint32_t i = 0; i < learnt_clause.size(); i++) {
            cout << "lev learnt_clause[" << i << "]:" << varData[learnt_clause[i].var()].level << endl;
        }
    }
}

size_t Searcher::find_backtrack_level_of_learnt()
{
    if (learnt_clause.size() <= 1)
        return 0;
    else {
        uint32_t max_i = 1;
        for (uint32_t i = 2; i < learnt_clause.size(); i++) {
            if (varData[learnt_clause[i].var()].level > varData[learnt_clause[max_i].var()].level)
                max_i = i;
        }
        std::swap(learnt_clause[max_i], learnt_clause[1]);
        return varData[learnt_clause[1].var()].level;
    }
}

template<bool update_bogoprops>
inline Clause* Searcher::create_learnt_clause(PropBy confl)
{
    pathC = 0;
    int index = trail.size() - 1;
    Lit p = lit_Undef;
    Clause* last_resolved_cl = NULL;

    learnt_clause.push_back(lit_Undef); //make space for ~p
    do {
        #ifdef DEBUG_RESOLV
        cout << "p is: " << p << endl;
        #endif

        //This is for OTF subsumption ("OTF clause improvement" by Han&Somezi)
        //~p is essentially popped from the temporary learnt clause
        if (p != lit_Undef) {
            if (!update_bogoprops && conf.doOTFSubsume) {
                tmp_learnt_clause_size--;
                assert(seen2[(~p).toInt()] == 1);
                seen2[(~p).toInt()] = 0;
            }

            //We MUST under-estimate
            tmp_learnt_clause_abst &= ~(abst_var((~p).var()));
        }

        last_resolved_cl = add_literals_from_confl_to_learnt<update_bogoprops>(confl, p);

        // Select next implication to look at
        while (!seen[trail[index--].var()]);

        p = trail[index+1];
        assert(p != lit_Undef);

        if (!update_bogoprops
            && pathC > 1
            && conf.doOTFSubsume
            //A long clause
            && last_resolved_cl != NULL
            //Good enough clause to try to minimize
            && (!last_resolved_cl->red() || last_resolved_cl->stats.glue <= conf.doOTFSubsumeOnlyAtOrBelowGlue)
            //Must subsume, so must be smaller
            && last_resolved_cl->size() > tmp_learnt_clause_size
            //Must not be a temporary clause
            && !last_resolved_cl->gauss_temp_cl()
            && !last_resolved_cl->used_in_xor()
        ) {
            last_resolved_cl->recalc_abst_if_needed();
            //Everything in learnt_cl_2 seems to be also in cl
            if ((last_resolved_cl->abst & tmp_learnt_clause_abst) ==  tmp_learnt_clause_abst
            ) {
                check_otf_subsume(confl.get_offset(), *last_resolved_cl);
            }
        }

        confl = varData[p.var()].reason;
        assert(varData[p.var()].level > 0);

        //This clears out vars that haven't been added to learnt_clause,
        //but their 'seen' has been set
        seen[p.var()] = 0;

        //Okay, one more path done
        pathC--;
    } while (pathC > 0);
    assert(pathC == 0);
    learnt_clause[0] = ~p;

    if (conf.doOTFSubsume
        && !update_bogoprops
    ) {
        for(const Lit lit: learnt_clause) {
            seen2[lit.toInt()] = 0;
        }
    }

    return last_resolved_cl;
}

void Searcher::simple_create_learnt_clause(
    PropBy confl,
    vector<Lit>& out_learnt,
    bool True_confl
) {
    int until = -1;
    int mypathC = 0;
    Lit p = lit_Undef;
    int index = trail.size() - 1;
    assert(decisionLevel() == 1);

    do {
        if (!confl.isNULL()) {
            if (confl.getType() == binary_t) {
                if (p == lit_Undef && True_confl == false) {
                    Lit q = failBinLit;
                    if (!seen[q.var()]) {
                        seen[q.var()] = 1;
                        mypathC++;
                    }
                }
                Lit q = confl.lit2();
                if (!seen[q.var()]) {
                    seen[q.var()] = 1;
                    mypathC++;
                }
            } else {
                const Clause& c = *solver->cl_alloc.ptr(confl.get_offset());

                // if True_confl==true, then choose p begin with the 1st index of c
                for (uint32_t j = (p == lit_Undef && True_confl == false) ? 0 : 1
                    ; j < c.size()
                    ; j++
                ) {
                    Lit q = c[j];
                    assert(q.var() < seen.size());
                    if (!seen[q.var()]) {
                        seen[q.var()] = 1;
                        mypathC++;
                    }
                }
            }
        } else {
            assert(confl.isNULL());
            out_learnt.push_back(~p);
        }
        // if not break, while() will come to the index of trail blow 0, and fatal error occur;
        if (mypathC == 0) {
            break;
        }
        // Select next clause to look at:
        while (!seen[trail[index--].var()]);
        // if the reason cr from the 0-level assigned var, we must break avoid move forth further;
        // but attention that maybe seen[x]=1 and never be clear. However makes no matter;
        if ((int)trail_lim[0] > index + 1
            && until == -1
        ) {
            until = out_learnt.size();
        }
        p = trail[index + 1];
        confl = varData[p.var()].reason;

        //under normal circumstances this does not happen, but here, it can
        //reason is undefined for level 0
        if (varData[p.var()].level == 0) {
            confl = PropBy();
        }
        seen[p.var()] = 0;
        mypathC--;
    } while (mypathC >= 0);

    if (until != -1)
        out_learnt.resize(until);
}

Clause* Searcher::otf_subsume_last_resolved_clause(Clause* last_resolved_cl)
{
    //We can only on-the-fly subsume with clauses that are not 2- or 3-long
    //furthermore, we cannot subsume a clause that is marked for deletion
    //due to its high glue value
    if (!conf.doOTFSubsume
        //Last was a lont clause
        || last_resolved_cl == NULL
        //Final clause will not be implicit
        || learnt_clause.size() <= 2
        //Larger or equivalent clauses cannot subsume the clause
        || learnt_clause.size() >= last_resolved_cl->size()
    ) {
        return NULL;
    }

    //Does it subsume?
    if (!subset(learnt_clause, *last_resolved_cl))
        return NULL;

    //on-the-fly subsumed the original clause
    stats.otfSubsumed++;
    stats.otfSubsumedLong++;
    stats.otfSubsumedRed += last_resolved_cl->red();
    stats.otfSubsumedLitsGained += last_resolved_cl->size() - learnt_clause.size();
    return last_resolved_cl;
}

void Searcher::print_debug_resolution_data(const PropBy confl)
{
#ifndef DEBUG_RESOLV
    //Avoid unused parameter warning
    (void) confl;
#else
    cout << "Before resolution, trail is: " << endl;
    print_trail();
    cout << "Conflicting clause: " << confl << endl;
    cout << "Fail bin lit: " << failBinLit << endl;
#endif
}

template<bool update_bogoprops>
Clause* Searcher::analyze_conflict(
    const PropBy confl
    , uint32_t& out_btlevel
    , uint32_t& glue
) {
    //Set up environment
    learnt_clause.clear();
    assert(toClear.empty());
    implied_by_learnts.clear();
    otf_subsuming_short_cls.clear();
    otf_subsuming_long_cls.clear();
    tmp_learnt_clause_size = 0;
    tmp_learnt_clause_abst = 0;
    assert(decisionLevel() > 0);

    print_debug_resolution_data(confl);
    Clause* last_resolved_cl = create_learnt_clause<update_bogoprops>(confl);
    stats.litsRedNonMin += learnt_clause.size();
    minimize_learnt_clause<update_bogoprops>();
    stats.litsRedFinal += learnt_clause.size();

    //further minimisation 1 -- short, small glue clauses
    glue = std::numeric_limits<uint32_t>::max();
    if (learnt_clause.size() <= conf.max_size_more_minim) {
        glue = calc_glue(learnt_clause);
        if (glue <= conf.max_glue_more_minim) {
            minimize_using_permdiff();
        }
    }
    if (glue == std::numeric_limits<uint32_t>::max()) {
        glue = calc_glue(learnt_clause);
    }
    print_fully_minimized_learnt_clause();

    if (learnt_clause.size() > conf.max_size_more_minim
        && glue <= (conf.glue_put_lev0_if_below_or_eq+2)
        && conf.doMinimRedMoreMore
    ) {
        minimise_redundant_more_more(learnt_clause);
    }

    out_btlevel = find_backtrack_level_of_learnt();
    if (!update_bogoprops) {
        bump_var_activities_based_on_implied_by_learnts<update_bogoprops>(out_btlevel);
    }
    sumConflictClauseLits += learnt_clause.size();

    return otf_subsume_last_resolved_clause(last_resolved_cl);

}

bool Searcher::litRedundant(const Lit p, uint32_t abstract_levels)
{
    #ifdef DEBUG_LITREDUNDANT
    cout << "c " << __func__ << " called" << endl;
    #endif

    analyze_stack.clear();
    analyze_stack.push(p);

    size_t top = toClear.size();
    while (!analyze_stack.empty()) {
        #ifdef DEBUG_LITREDUNDANT
        cout << "At point in litRedundant: " << analyze_stack.top() << endl;
        #endif

        const PropBy reason = varData[analyze_stack.top().var()].reason;
        PropByType type = reason.getType();
        analyze_stack.pop();

        //Must have a reason
        assert(!reason.isNULL());

        size_t size;
        Clause* cl = NULL;
        switch (type) {
            case clause_t:
                cl = cl_alloc.ptr(reason.get_offset());
                size = cl->size()-1;
                break;

            case binary_t:
                size = 1;
                break;

            case null_clause_t:
            default:
                release_assert(false);
        }

        for (size_t i = 0
            ; i < size
            ; i++
        ) {
            Lit p2;
            switch (type) {
                case clause_t:
                    p2 = (*cl)[i+1];
                    break;

                case binary_t:
                    p2 = reason.lit2();
                    break;

                case null_clause_t:
                default:
                    release_assert(false);
            }
            stats.recMinimCost++;

            if (!seen[p2.var()] && varData[p2.var()].level > 0) {
                if (!varData[p2.var()].reason.isNULL()
                    && (abstractLevel(p2.var()) & abstract_levels) != 0
                ) {
                    seen[p2.var()] = 1;
                    analyze_stack.push(p2);
                    toClear.push_back(p2);
                } else {
                    //Return to where we started before function executed
                    for (size_t j = top; j < toClear.size(); j++) {
                        seen[toClear[j].var()] = 0;
                    }
                    toClear.resize(top);

                    return false;
                }
            }
        }
    }

    return true;
}
template Clause* Searcher::analyze_conflict<true>(const PropBy confl
    , uint32_t& out_btlevel
    , uint32_t& glue
);
template Clause* Searcher::analyze_conflict<false>(const PropBy confl
    , uint32_t& out_btlevel
    , uint32_t& glue
);

bool Searcher::subset(const vector<Lit>& A, const Clause& B)
{
    //Set seen
    for (uint32_t i = 0; i != B.size(); i++)
        seen[B[i].toInt()] = 1;

    bool ret = true;
    for (uint32_t i = 0; i != A.size(); i++) {
        if (!seen[A[i].toInt()]) {
            ret = false;
            break;
        }
    }

    //Clear seen
    for (uint32_t i = 0; i != B.size(); i++)
        seen[B[i].toInt()] = 0;

    return ret;
}

void Searcher::analyze_final_confl_with_assumptions(const Lit p, vector<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push_back(p);

    if (decisionLevel() == 0) {
        return;
    }

    //It's been set at level 0. The seen[] may not be large enough to do
    //seen[p.var()] -- we might have mem-saved that
    if (varData[p.var()].level == 0) {
        return;
    }

    seen[p.var()] = 1;

    assert(!trail_lim.empty());
    for (int64_t i = (int64_t)trail.size() - 1; i >= (int64_t)trail_lim[0]; i--) {
        const uint32_t x = trail[i].var();
        if (seen[x]) {
            const PropBy reason = varData[x].reason;
            if (reason.isNULL()) {
                assert(varData[x].level > 0);
                out_conflict.push_back(~trail[i]);
            } else {
                switch(reason.getType()) {
                    case PropByType::clause_t : {
                        const Clause& cl = *cl_alloc.ptr(reason.get_offset());
                        assert(value(cl[0]) == l_True);
                        for(const Lit lit: cl) {
                            if (varData[lit.var()].level > 0) {
                                seen[lit.var()] = 1;
                            }
                        }
                        break;
                    }
                    case PropByType::binary_t: {
                        const Lit lit = reason.lit2();
                        if (varData[lit.var()].level > 0) {
                            seen[lit.var()] = 1;
                        }
                        break;
                    }

                    default:
                        assert(false);
                        break;
                }
            }
            seen[x] = 0;
        }
    }
    seen[p.var()] = 0;
}

void Searcher::update_assump_conflict_to_orig_outside(vector<Lit>& out_conflict)
{
    if (assumptions.empty()) {
        return;
    }

    vector<AssumptionPair> inter_assumptions;
    for(const auto& ass: assumptions) {
        inter_assumptions.push_back(
            AssumptionPair(map_outer_to_inter(ass.lit_outer), ass.lit_orig_outside));
    }

    std::sort(inter_assumptions.begin(), inter_assumptions.end());
    std::sort(out_conflict.begin(), out_conflict.end());
    assert(out_conflict.size() <= assumptions.size());
    //They now are in the order where we can go through them linearly

    /*cout << "out_conflict: " << out_conflict << endl;
    cout << "assumptions: ";
    for(AssumptionPair p: assumptions) {
        cout << "inter: " << p.lit_inter << " , outer: " << p.lit_orig_outside << " , ";
    }
    cout << endl;*/

    uint32_t at_assump = 0;
    uint32_t j = 0;
    for(size_t i = 0; i < out_conflict.size(); i++) {
        Lit lit = out_conflict[i];

        //lit_outer is actually INTER here, because we updated above
        while(lit != ~inter_assumptions[at_assump].lit_outer) {
            at_assump++;
            assert(at_assump < inter_assumptions.size() && "final conflict contains literals that are not from the assumptions!");
        }
        assert(lit == ~inter_assumptions[at_assump].lit_outer);

        //in case of symmetry breaking, we can be in trouble
        //then, the orig_outside is actually lit_Undef
        //in these cases, the symmetry breaking literal needs to be taken out
        if (inter_assumptions[at_assump].lit_orig_outside != lit_Undef) {
            //Update to correct outside lit
            out_conflict[j++] = ~inter_assumptions[at_assump].lit_orig_outside;
        }
    }
    out_conflict.resize(j);
}

void Searcher::check_blocking_restart()
{
    if (conf.do_blocking_restart
        && sumConflicts > conf.lower_bound_for_blocking_restart
        && hist.glueHist.isvalid()
        && hist.trailDepthHistLonger.isvalid()
        && decisionLevel() > 0
        && trail_lim.size() > 0
        && trail.size() > hist.trailDepthHistLonger.avg()*conf.blocking_restart_multip
    ) {
        hist.glueHist.clear();
        if (!blocked_restart) {
            stats.blocked_restart_same++;
        }
        blocked_restart = true;
        stats.blocked_restart++;
    }
}

lbool Searcher::search()
{
    assert(ok);
    #ifdef SLOW_DEBUG
    check_no_duplicate_lits_anywhere();
    check_order_heap_sanity();
    #endif

    //Stats reset & update
    stats.numRestarts++;
    stats.clauseID_at_start_inclusive = clauseID;
    hist.clear();
    hist.reset_glue_hist_size(conf.shortTermHistorySize);

    assert(solver->prop_at_head());

    //Loop until restart or finish (SAT/UNSAT)
    blocked_restart = false;
    PropBy confl;
    lbool dec_ret = l_Undef;

    while (!params.needToStopSearch
        || !confl.isNULL() //always finish the last conflict
    ) {
        #ifdef USE_GAUSS
        gqhead = qhead;
        #endif
        confl = propagate_any_order_fast();

        if (!confl.isNULL()) {
            //manipulate startup parameters
            if (
                ((stats.conflStats.numConflicts & 0xfff) == 0xfff) &&
                var_decay_vsids < conf.var_decay_vsids_max
            ) {
                var_decay_vsids += 0.01;
            }
            print_restart_stat();
            hist.trailDepthHistLonger.push(trail.size());
            if (!handle_conflict<false>(confl)) {
                return l_False;
            }
            check_need_restart();
        } else {
            assert(ok);
            #ifdef USE_GAUSS
            gauss_ret ret = gauss_jordan_elim();
            //cout << "ret: " << ret << " -- " << endl;
            if (ret == gauss_ret::g_cont) {
                //cout << "g_cont" << endl;
                check_need_restart();
                continue;
            //TODO conflict should be goto-d to "confl" label
            }

            if (ret == gauss_ret::g_false) {
                return l_False;
            }

            assert(ret == gauss_ret::g_nothing);
            //cout << "g_nothing" << endl;
            #endif //USE_GAUSS

            if (decisionLevel() == 0
                && !clean_clauses_if_needed()
            ) {
                return l_False;
            };
            reduce_db_if_needed();
            dec_ret = new_decision<false>();
            if (dec_ret != l_Undef) {
                return dec_ret;
            }
        }
    }
    max_confl_this_phase -= (int64_t)params.conflictsDoneThisRestart;

    cancelUntil<true, false>(0);
    confl = propagate<false>();
    if (!confl.isNULL()) {
        ok = false;
        return l_False;
    }
    assert(solver->prop_at_head());
    return l_Undef;
}

/**
@brief Picks a new decision variable to branch on

@returns l_Undef if it should restart instead. l_False if it reached UNSAT
         (through simplification)
*/
template<bool update_bogoprops>
lbool Searcher::new_decision()
{
    Lit next = lit_Undef;
    while (decisionLevel() < assumptions.size()) {
        // Perform user provided assumption:
        Lit p = map_outer_to_inter(assumptions[decisionLevel()].lit_outer);
        assert(varData[p.var()].removed == Removed::none);

        if (value(p) == l_True) {
            // Dummy decision level:
            new_decision_level();
        } else if (value(p) == l_False) {
            analyze_final_confl_with_assumptions(~p, conflict);
            return l_False;
        } else {
            assert(p.var() < nVars());
            stats.decisionsAssump++;
            next = p;
            break;
        }
    }

    if (next == lit_Undef) {
        // New variable decision:
        next = pickBranchLit();

        //No decision taken, because it's SAT
        if (next == lit_Undef)
            return l_True;

        //Update stats
        stats.decisions++;
        sumDecisions++;
    }

    // Increase decision level and enqueue 'next'
    assert(value(next) == l_Undef);
    new_decision_level();
    enqueue<update_bogoprops>(next);

    return l_Undef;
}

double Searcher::luby(double y, int x)
{
    int size = 1;
    int seq;
    for (seq = 0
        ; size < x + 1
        ; seq++
    ) {
        size = 2 * size + 1;
    }

    while (size - 1 != x) {
        size = (size - 1) >> 1;
        seq--;
        x = x % size;
    }

    return std::pow(y, seq);
}

void Searcher::check_need_restart()
{
    if ((stats.conflStats.numConflicts & 0xff) == 0xff) {
        //It's expensive to check the time the time
        if (cpuTime() > conf.maxTime) {
            params.needToStopSearch = true;
        }

        if (must_interrupt_asap())  {
            if (conf.verbosity >= 3)
                cout << "c must_interrupt_asap() is set, restartig as soon as possible!" << endl;
            params.needToStopSearch = true;
        }
    }

    if ((params.rest_type == Restart::fixed)
        && (int64_t)params.conflictsDoneThisRestart > max_confl_this_phase
    ) {
        params.needToStopSearch = true;
    }

    //Conflict limit reached?
    if (params.conflictsDoneThisRestart > params.max_confl_to_do) {
        if (conf.verbosity >= 3) {
            cout
            << "c Over limit of conflicts for this restart"
            << " -- restarting as soon as possible!" << endl;
        }
        params.needToStopSearch = true;
    }
}

template<bool update_bogoprops>
void Searcher::add_otf_subsume_long_clauses()
{
    //Hande long OTF subsumption
    for(size_t i = 0; i < otf_subsuming_long_cls.size(); i++) {
        const ClOffset offset = otf_subsuming_long_cls[i];
        Clause& cl = *solver->cl_alloc.ptr(offset);

        //Find the l_Undef
        size_t at = std::numeric_limits<size_t>::max();
        for(size_t i2 = 0; i2 < cl.size(); i2++) {
            if (value(cl[i2]) == l_Undef) {
                at = i2;
                break;
            }
        }
        assert(at != std::numeric_limits<size_t>::max());
        std::swap(cl[at], cl[0]);
        assert(value(cl[0]) == l_Undef);

        //Find another l_Undef or an l_True
        at = 0;
        for(size_t i2 = 1; i2 < cl.size(); i2++) {
            if (value(cl[i2]) == l_Undef || value(cl[i2]) == l_True) {
                at = i2;
                break;
            }
        }
        assert(cl.size() > 2);

        if (at == 0) {
            //If none found, we have a propagating clause_t
            enqueue<update_bogoprops>(cl[0], decisionLevel() == 0 ? PropBy() : PropBy(offset));

            //Drat
            if (decisionLevel() == 0) {
                *drat << add << cl[0]
                << fin;
            }
        } else {
            //We have a non-propagating clause

            std::swap(cl[at], cl[1]);
            assert(value(cl[1]) == l_Undef || value(cl[1]) == l_True);
        }
        solver->attachClause(cl, false);
        cl.setStrenghtened();
    }
    otf_subsuming_long_cls.clear();
}

template<bool update_bogoprops>
void Searcher::add_otf_subsume_implicit_clause()
{
    //Handle implicit OTF subsumption
    for(vector<OTFClause>::iterator
        it = otf_subsuming_short_cls.begin(), end = otf_subsuming_short_cls.end()
        ; it != end
        ; ++it
    ) {
        assert(it->size > 1);
        //Find the l_Undef
        size_t at = std::numeric_limits<size_t>::max();
        for(size_t i2 = 0; i2 < it->size; i2++) {
            if (value(it->lits[i2]) == l_Undef) {
                at = i2;
                break;
            }
        }
        assert(at != std::numeric_limits<size_t>::max());
        std::swap(it->lits[at], it->lits[0]);
        assert(value(it->lits[0]) == l_Undef);

        //Find another l_Undef or an l_True
        at = 0;
        for(size_t i2 = 1; i2 < it->size; i2++) {
            if (value(it->lits[i2]) == l_Undef
                || value(it->lits[i2]) == l_True
            ) {
                at = i2;
                break;
            }
        }

        if (at == 0) {
            //If none found, we have a propagation
            //Calculate reason
            PropBy by = PropBy();

            //if decision level is non-zero, we have to be more careful
            if (decisionLevel() != 0) {
                assert(it->size == 2);
                by = PropBy(it->lits[1], true);
            }

            //Enqueue this literal, finally
            enqueue<update_bogoprops>(
                it->lits[0]
                , by
            );

            //Drat
            if (decisionLevel() == 0) {
                *drat << add << it->lits[0]
                << fin;
            }
        } else {
            //We have a non-propagating clause
            std::swap(it->lits[at], it->lits[1]);
            assert(value(it->lits[1]) == l_Undef
                || value(it->lits[1]) == l_True
            );

            //Attach new binary/tertiary clause
            if (it->size == 2) {
                solver->attach_bin_clause(it->lits[0], it->lits[1], true);
            }
        }
    }
    otf_subsuming_short_cls.clear();
}

void Searcher::update_history_stats(size_t backtrack_level, uint32_t glue)
{
    assert(decisionLevel() > 0);

    //short-term averages
    hist.branchDepthHist.push(decisionLevel());
    hist.branchDepthDeltaHist.push(decisionLevel() - backtrack_level);
    hist.conflSizeHist.push(learnt_clause.size());
    hist.trailDepthDeltaHist.push(trail.size() - trail_lim[backtrack_level]);

    //long-term averages
    hist.backtrackLevelHistLT.push(backtrack_level);
    hist.conflSizeHistLT.push(learnt_clause.size());
    hist.trailDepthHistLT.push(trail.size());
    hist.glueHistLT.push(glue);
    hist.glueHist.push(glue);

    //Global stats from cnf.h
    sumClLBD += glue;
    sumClSize += learnt_clause.size();
}

template<bool update_bogoprops>
void Searcher::attach_and_enqueue_learnt_clause(Clause* cl, bool enq)
{
    switch (learnt_clause.size()) {
        case 0:
            assert(false);
        case 1:
            //Unitary learnt
            stats.learntUnits++;
            if (enq) enqueue(learnt_clause[0]);
            assert(decisionLevel() == 0);
            break;
        case 2:
            //Binary learnt
            stats.learntBins++;
            solver->attach_bin_clause(learnt_clause[0], learnt_clause[1], true, enq);
            if (enq) enqueue(learnt_clause[0], PropBy(learnt_clause[1], true));
            break;

        default:
            //Long learnt
            stats.learntLongs++;
            solver->attachClause(*cl, enq);
            if (enq) enqueue(learnt_clause[0], PropBy(cl_alloc.get_offset(cl)));
            for(uint32_t i = 0; i < solver->conf.bump_new_learnt_cls; i++) {
                bump_cl_act<update_bogoprops>(cl);
            }

            break;
    }
}

inline void Searcher::print_learning_debug_info() const
{
    #ifndef VERBOSE_DEBUG
    return;
    #else
    cout
    << "Learning:" << learnt_clause
    << endl
    << "reverting var " << learnt_clause[0].var()+1
    << " to " << !learnt_clause[0].sign()
    << endl;
    #endif
}

void Searcher::print_learnt_clause() const
{
    if (conf.verbosity >= 6) {
        cout
        << "c learnt clause: "
        << learnt_clause
        << endl;
    }
}

Clause* Searcher::handle_last_confl_otf_subsumption(
    Clause* cl
    , const uint32_t glue
    , bool decision_cl
) {
    if (learnt_clause.size() <= 2 ||
        cl == NULL ||
        cl->gauss_temp_cl() ||
        !conf.doOTFSubsume
    ) {
        //Cannot make a non-implicit into an implicit
        if (learnt_clause.size() <= 2) {
            *drat << add << learnt_clause
            << fin;
            cl = NULL;
        } else {
            cl = cl_alloc.Clause_new(learnt_clause
            , sumConflicts
            );
            cl->makeRed(glue);
            ClOffset offset = cl_alloc.get_offset(cl);
            unsigned which_arr = 2;

            if (cl->stats.locked_for_data_gen) {
                which_arr = 0;
            } else if (glue <= conf.glue_put_lev0_if_below_or_eq) {
                which_arr = 0;
            } else if (
                glue <= conf.glue_put_lev1_if_below_or_eq
                && conf.glue_put_lev1_if_below_or_eq != 0
            ) {
                which_arr = 1;
            } else {
                which_arr = 2;
            }

            if (which_arr == 0) {
                stats.red_cl_in_which0++;
            }

            cl->stats.which_red_array = which_arr;
            cl->stats.is_decision_cl = decision_cl;
            solver->longRedCls[cl->stats.which_red_array].push_back(offset);

            *drat << add << *cl
            << fin;
        }
    } else {
        //On-the-fly subsumption
        assert(cl->size() > 2);
        *drat << deldelay << *cl << fin;
        solver->detachClause(*cl, false);

        //Shrink clause
        assert(cl->size() > learnt_clause.size());
        for (uint32_t i = 0; i < learnt_clause.size(); i++) {
            (*cl)[i] = learnt_clause[i];
        }
        cl->resize(learnt_clause.size());
        assert(cl->size() == learnt_clause.size());

        //Update stats
        if (cl->red() && cl->stats.glue > glue) {
            cl->stats.glue = glue;
        }

        *drat << add << *cl
         << fin << findelay;
    }

    return cl;
}

template<bool update_bogoprops>
bool Searcher::handle_conflict(const PropBy confl)
{
    if (!update_bogoprops) {
        stats.conflStats.numConflicts++;
        sumConflicts++;

        if (sumConflicts == 100000 && //TODO magic constant
            longRedCls[0].size() < 100 &&
            //so that in case of some "standard-minisat behavriour" config
            //we don't override it
            conf.glue_put_lev0_if_below_or_eq != 0
        ) {
            conf.glue_put_lev0_if_below_or_eq += 2; //TODO magic constant
        }
        params.conflictsDoneThisRestart++;
    }

    if (decisionLevel() == 0)
        return false;

    uint32_t backtrack_level;
    uint32_t glue;
    Clause* subsumed_cl = analyze_conflict<update_bogoprops>(
        confl
        , backtrack_level  //return backtrack level here
        , glue             //return glue here
    );
    print_learnt_clause();

    //Add decision-based clause in case it's short
    decision_clause.clear();
    if (!update_bogoprops
        && conf.do_decision_based_cl
        && learnt_clause.size() > conf.decision_based_cl_min_learned_size
        && decisionLevel() <= conf.decision_based_cl_max_levels
        && decisionLevel() >= 2
    ) {
        for(int i = (int)trail_lim.size()-1; i >= 0; i--) {
            Lit l = ~trail[trail_lim[i]];
            if (!seen[l.toInt()]) {
                decision_clause.push_back(l);
                seen[l.toInt()] = 1;
            }
        }
        for(Lit l: decision_clause) {
            seen[l.toInt()] = 0;
        }
    }

    if (!update_bogoprops) {
        update_history_stats(backtrack_level, glue);
    }
    uint32_t plus = learnt_clause.size() > 2;
    plus += decision_clause.size() > 2;
    cancelUntil<true, update_bogoprops>(backtrack_level, plus);

    add_otf_subsume_long_clauses<update_bogoprops>();
    add_otf_subsume_implicit_clause<update_bogoprops>();
    print_learning_debug_info();
    assert(value(learnt_clause[0]) == l_Undef);
    glue = std::min<uint32_t>(glue, std::numeric_limits<uint32_t>::max());
    Clause* cl = handle_last_confl_otf_subsumption(
        subsumed_cl
        , glue
        , false //decision clause?
    );
    assert(learnt_clause.size() <= 2 || cl != NULL);
    attach_and_enqueue_learnt_clause<update_bogoprops>(cl);

    //Add decision-based clause
    if (!update_bogoprops
        && decision_clause.size() > 0
    ) {
        int i = decision_clause.size();
        while(--i >= 0) {
            if (value(decision_clause[i]) == l_True
                || value(decision_clause[i]) == l_Undef
            ) {
                break;
            }
        }
        std::swap(decision_clause[0], decision_clause[i]);
        learnt_clause = decision_clause;
        sumDecisionBasedCl++;
        cl = handle_last_confl_otf_subsumption(
            NULL                   //orig cl to minimise
            , learnt_clause.size() //glue
            , true //decision clause?
        );
        attach_and_enqueue_learnt_clause<update_bogoprops>(cl, false);
    }

    if (!update_bogoprops) {
        varDecayActivity();
        decayClauseAct<update_bogoprops>();
    }

    return true;
}
template bool Searcher::handle_conflict<true>(const PropBy confl);
template bool Searcher::handle_conflict<false>(const PropBy confl);

void Searcher::resetStats()
{
    startTime = cpuTime();

    //Rest solving stats
    stats.clear();
    propStats.clear();
    lastCleanZeroDepthAssigns = trail.size();
}

void Searcher::print_restart_header()
{
    //Print restart output header
    if ((lastRestartPrintHeader == 0 || (lastRestartPrintHeader + 1600000) < sumConflicts)
        && conf.verbosity
    ) {
        cout
        << "c"
        << " " << std::setw(6) << "type"
        << " " << std::setw(5) << "VSIDS"
        << " " << std::setw(5) << "rest"
        << " " << std::setw(5) << "conf"
        << " " << std::setw(5) << "freevar"
        << " " << std::setw(5) << "IrrL"
        << " " << std::setw(5) << "IrrB"
        << " " << std::setw(7) << "l/longC"
        << " " << std::setw(7) << "l/allC";

        for(size_t i = 0; i < longRedCls.size(); i++) {
            cout << " " << std::setw(4) << "RedL" << i;
        }

        cout
        << " " << std::setw(5) << "RedB"
        << " " << std::setw(7) << "l/longC"
        << " " << std::setw(7) << "l/allC"
        << endl;
        lastRestartPrintHeader = sumConflicts;
    }
}

void Searcher::print_restart_stat_line() const
{
    print_restart_stats_base();
    if (conf.print_full_restart_stat) {
        solver->print_clause_stats();
        hist.print();
    } else {
        solver->print_clause_stats();
    }

    cout << endl;
}

void Searcher::print_restart_stats_base() const
{
    cout << "c"
         << " " << std::setw(6) << restart_type_to_short_string(params.rest_type);
    cout << " " << std::setw(5) << sumRestarts();

    if (sumConflicts >  20000) {
        cout << " " << std::setw(4) << sumConflicts/1000 << "K";
    } else {
        cout << " " << std::setw(5) << sumConflicts;
    }

    cout << " " << std::setw(7) << solver->get_num_free_vars();
}

struct MyInvSorter {
    bool operator()(size_t num, size_t num2)
    {
        return num > num2;
    }
};

struct MyPolarData
{
    MyPolarData (size_t _pos, size_t _neg, size_t _flipped) :
        pos(_pos)
        , neg(_neg)
        , flipped(_flipped)
    {}

    size_t pos;
    size_t neg;
    size_t flipped;

    bool operator<(const MyPolarData& other) const
    {
        return (pos + neg) > (other.pos + other.neg);
    }
};

void Searcher::print_restart_stat()
{
    //Print restart stat
    if (conf.verbosity
        && !conf.print_all_restarts
        && ((lastRestartPrint + conf.print_restart_line_every_n_confl)
          < sumConflicts)
    ) {
        print_restart_stat_line();
        lastRestartPrint = sumConflicts;
    }
}

void Searcher::reset_temp_cl_num()
{
    cur_max_temp_red_lev2_cls = conf.max_temp_lev2_learnt_clauses;
}

void Searcher::reduce_db_if_needed()
{
    if (conf.every_lev1_reduce != 0
        && sumConflicts >= next_lev1_reduce
    ) {
        solver->reduceDB->handle_lev1();
        next_lev1_reduce = sumConflicts + conf.every_lev1_reduce;
    }

    if (conf.every_lev2_reduce != 0) {
        if (sumConflicts >= next_lev2_reduce) {
            solver->reduceDB->handle_lev2();
            cl_alloc.consolidate(solver);
            next_lev2_reduce = sumConflicts + conf.every_lev2_reduce;
        }
    } else {
        if (longRedCls[2].size() > cur_max_temp_red_lev2_cls) {
            solver->reduceDB->handle_lev2();
            cur_max_temp_red_lev2_cls *= conf.inc_max_temp_lev2_red_cls;
            cl_alloc.consolidate(solver);
        }
    }
}

bool Searcher::clean_clauses_if_needed()
{
    assert(decisionLevel() == 0);

    if (!ok || !propagate_any_order_fast().isNULL()) {
        return ok = false;
    }

    const size_t newZeroDepthAss = trail.size() - lastCleanZeroDepthAssigns;
    if (newZeroDepthAss > 0
        && simpDB_props < 0
        && newZeroDepthAss > ((double)nVars()*0.05)
    ) {
        if (conf.verbosity >= 2) {
            cout << "c newZeroDepthAss : " << newZeroDepthAss
            << " -- "
            << (double)newZeroDepthAss/(double)nVars()*100.0
            << " % of active vars"
            << endl;
        }
        lastCleanZeroDepthAssigns = trail.size();
        solver->clauseCleaner->remove_and_clean_all();

        cl_alloc.consolidate(solver);
        rebuildOrderHeap(); //TODO only filter is needed!
        simpDB_props = (litStats.redLits + litStats.irredLits)<<5;
    }

    return true;
}

void Searcher::rebuildOrderHeap()
{
    vector<uint32_t> vs;
    for (uint32_t v = 0; v < nVars(); v++) {
        if (varData[v].removed != Removed::none
            //NOTE: the level==0 check is needed because SLS calls this
            //when there is a solution already, but we should only skip
            //level 0 assignements
            || (value(v) != l_Undef && varData[v].level == 0)
        ) {
            continue;
        } else {
            vs.push_back(v);
        }
    }
    order_heap_vsids.build(vs);
}

bool Searcher::must_abort(const lbool status) {
    if (status != l_Undef) {
        if (conf.verbosity >= 6) {
            cout
            << "c Returned status of search() is non-l_Undef at confl:"
            << sumConflicts
            << endl;
        }
        return true;
    }

    if (stats.conflStats.numConflicts >= max_confl_per_search_solve_call) {
        if (conf.verbosity >= 3) {
            cout
            << "c search over max conflicts"
            << endl;
        }
        return true;
    }

    if (cpuTime() >= conf.maxTime) {
        if (conf.verbosity >= 3) {
            cout
            << "c search over max time"
            << endl;
        }
        return true;
    }

    if (solver->must_interrupt_asap()) {
        if (conf.verbosity >= 3) {
            cout
            << "c search interrupting as requested"
            << endl;
        }
        return true;
    }

    return false;
}

lbool Searcher::solve(
    const uint64_t _max_confls
) {
    assert(ok);
    assert(qhead == trail.size());
    max_confl_per_search_solve_call = _max_confls;
    num_search_called++;
    #ifdef SLOW_DEBUG
    //When asking for a lot of simple soluitons, search() gets called a lot
    check_no_removed_or_freed_cl_in_watch();
    #endif

    if (solver->conf.verbosity >= 6) {
        cout
        << "c Searcher::solve() called"
        << endl;
    }

    resetStats();
    lbool status = l_Undef;
    if (conf.restartType == Restart::fixed) {
        params.rest_type = Restart::fixed;
        max_confl_this_phase = conf.fixed_restart_num_confl;
    }

    #ifdef USE_GAUSS
    if (!solver->xor_clauses_updated) {
        if (conf.verbosity >= 3) {
            cout << "c [find&init matx] XORs not updated, and either (XORs are not detached OR assumps does not contain clash variable) -> or not performing matrix init. Matrices: " << gmatrices.size() << endl;
        }
    } else {
        if (conf.verbosity >= 1) cout << "c [find&init matx] performing matrix init" << endl;
        clear_gauss_matrices();
        {
            MatrixFinder finder(solver);
            ok = finder.findMatrixes();
            if (!ok) {
                status = l_False;
                goto end;
            }
        }
        if (!solver->init_all_matrices()) {
            return l_False;
        }
    }

    #ifdef SLOW_DEBUG
    for(size_t i = 0; i< solver->gmatrixes.size(); i++) {
        if (solver->gmatrixes[i]) {
            solver->gmatrixes[i]->check_watchlist_sanity();
            assert(solver->gmatrixes[i]->get_matrix_no() == i);
        }
    }
    #endif

    #endif //USE_GAUSS

    assert(solver->check_order_heap_sanity());
    while(stats.conflStats.numConflicts < max_confl_per_search_solve_call
        && status == l_Undef
    ) {
        #ifdef SLOW_DEBUG
        assert(solver->check_order_heap_sanity());
        #endif

        assert(watches.get_smudged_list().empty());

        lastRestartConfl = sumConflicts;
        params.clear();
        params.max_confl_to_do = max_confl_per_search_solve_call-stats.conflStats.numConflicts;
        status = search();
        if (status == l_Undef) {
            adjust_phases_restarts();
        }

        if (must_abort(status)) {
            goto end;
        }

        if (status == l_Undef &&
            solver->conf.do_distill_clauses &&
            sumConflicts > next_distill
        ) {
            if (!solver->distill_long_cls->distill(true, false)) {
                status = l_False;
                goto end;
            }
            next_distill = std::min<double>(sumConflicts * 0.2 + sumConflicts + 3000,
                                    sumConflicts + 50000);
        }
    }

    end:
    finish_up_solve(status);

    return status;
}

void Searcher::adjust_phases_restarts()
{
    //Haven't finished the phase. Keep rolling.
    if (max_confl_this_phase > 0)
        return;

    switch(conf.restartType) {
    case Restart::fixed:
        assert(params.rest_type == Restart::fixed);
        max_confl_this_phase = conf.fixed_restart_num_confl;
        break;
    }
}

void Searcher::print_solution_varreplace_status() const
{
    for(size_t var = 0; var < nVarsOuter(); var++) {
        if (varData[var].removed == Removed::replaced
            || varData[var].removed == Removed::elimed
        ) {
            assert(value(var) == l_Undef || varData[var].level == 0);
        }

        if (conf.verbosity >= 6
            && varData[var].removed == Removed::replaced
            && value(var) != l_Undef
        ) {
            cout
            << "var: " << var
            << " value: " << value(var)
            << " level:" << varData[var].level
            << " type: " << removed_type_to_string(varData[var].removed)
            << endl;
        }
    }
}

void Searcher::print_solution_type(const lbool status) const
{
    if (conf.verbosity >= 6) {
        if (status == l_True) {
            cout << "Solution from Searcher is SAT" << endl;
        } else if (status == l_False) {
            cout << "Solution from Searcher is UNSAT" << endl;
            cout << "OK is: " << okay() << endl;
        } else {
            cout << "Solutions from Searcher is UNKNOWN" << endl;
        }
    }
}

void Searcher::finish_up_solve(const lbool status)
{
    print_solution_type(status);

    if (status == l_True) {
        #ifdef SLOW_DEBUG
        check_order_heap_sanity();
        #endif
        model = assigns;

        if (conf.need_decisions_reaching) {
            for(size_t i = 0; i < trail_lim.size(); i++) {
                size_t at = trail_lim[i];

                //we need this due to dummy decision levels
                //that could create a situation where new_decision_level()
                //has been called, but then no variable needs to be decided
                //for SAT.
                if (at < trail.size()) {
                    decisions_reaching_model.push_back(trail[at]);
                }
            }
        }

        cancelUntil(0);
        print_solution_varreplace_status();
    } else if (status == l_False) {
        if (conflict.size() == 0) {
            ok = false;
        }
        cancelUntil(0);
    }

    stats.cpu_time = cpuTime() - startTime;
    if (conf.verbosity >= 4) {
        cout << "c Searcher::solve() finished"
        << " status: " << status
        << " numConflicts : " << stats.conflStats.numConflicts
        << " SumConfl: " << sumConflicts
        << " max_confl_per_search_solve_call:" << max_confl_per_search_solve_call
        << endl;
    }

    print_iteration_solving_stats();
}

void Searcher::print_iteration_solving_stats()
{
    if (conf.verbosity >= 3) {
        cout << "c ------ THIS ITERATION SOLVING STATS -------" << endl;
        stats.print(propStats.propagations, conf.do_print_times);
        propStats.print(stats.cpu_time);
        print_stats_line("c props/decision"
            , float_div(propStats.propagations, stats.decisions)
        );
        print_stats_line("c props/conflict"
            , float_div(propStats.propagations, stats.conflStats.numConflicts)
        );
        cout << "c ------ THIS ITERATION SOLVING STATS -------" << endl;
    }
}

Lit Searcher::pickBranchLit()
{
    #ifdef VERBOSE_DEBUG
    cout << "picking decision variable, dec. level: " << decisionLevel()
    #endif

    Lit next = lit_Undef;

    // Random decision:
    Heap<VarOrderLt> &order_heap = order_heap_vsids;
    vector<double>& var_act = var_act_vsids;

    if (conf.random_var_freq > 0) {
        double rand = mtrand.randDblExc();
        double frq = conf.random_var_freq;
        if (rand < frq && !order_heap.empty()) {
            uint32_t next_var = var_Undef;
            while(!order_heap.empty()
                && next_var == var_Undef)
            {
                next_var = order_heap.random_element(mtrand);
                if (value(next_var) == l_Undef
                    && solver->varData[next_var].removed == Removed::none
                ) {
                    stats.decisionsRand++;
                    next = Lit(next_var, !pick_polarity(next_var));
                } else {
                    //make this var the top, and remove it
                    assert(var_act.size() > next_var);
                    assert(order_heap.inHeap(next_var));
                    double backup = var_act[order_heap.inspectTop()];
                    var_act[next_var] = var_act[order_heap.inspectTop()]*2+10e2;
                    order_heap.update(next_var);
                    uint32_t removed_var = (uint32_t)order_heap.removeMin();
                    assert(removed_var == next_var);
                    var_act[next_var] = backup;

                    next_var = var_Undef;
                }
            }
        }
    }

    if (next == lit_Undef) {
        uint32_t v = var_Undef;
        while (v == var_Undef || value(v) != l_Undef) {
            //There is no more to branch on. Satisfying assignment found.
            if (order_heap.empty()) {
                return lit_Undef;
            }
            v = order_heap.removeMin();
        }
        next = Lit(v, !pick_polarity(v));
    }

    //No vars in heap: solution found
    #ifdef SLOW_DEBUG
    if (next != lit_Undef) {
        assert(solver->varData[next.var()].removed == Removed::none);
    }
    #endif
    return next;
}

void Searcher::cache_based_morem_minim(vector<Lit>& cl)
{
    int64_t limit = more_red_minim_limit_cache_actual;
    const size_t first_n_lits_of_cl =
        std::min<size_t>(conf.max_num_lits_more_more_red_min, cl.size());
    for (size_t at_lit = 0; at_lit < first_n_lits_of_cl; at_lit++) {
        Lit lit = cl[at_lit];

        //Timeout
        if (limit < 0)
            break;

        //Already removed this literal
        if (seen[lit.toInt()] == 0)
            continue;

        assert(solver->implCache.size() > lit.toInt());
        const TransCache& cache1 = solver->implCache[lit];
        limit -= (int64_t)cache1.lits.size()/2;
        for (const LitExtra litExtra: cache1.lits) {
            assert(seen.size() > litExtra.getLit().toInt());
            if (seen[(~(litExtra.getLit())).toInt()]) {
                stats.cacheShrinkedClause++;
                seen[(~(litExtra.getLit())).toInt()] = 0;
            }
        }
    }
}

void Searcher::binary_based_morem_minim(vector<Lit>& cl)
{
    int64_t limit  = more_red_minim_limit_binary_actual;
    const size_t first_n_lits_of_cl =
        std::min<size_t>(conf.max_num_lits_more_more_red_min, cl.size());
    for (size_t at_lit = 0; at_lit < first_n_lits_of_cl; at_lit++) {
        Lit lit = cl[at_lit];
        //Already removed this literal
        if (seen[lit.toInt()] == 0)
            continue;

        //Watchlist-based minimisation
        watch_subarray_const ws = watches[lit];
        for (const Watched* i = ws.begin() , *end = ws.end()
            ; i != end && limit > 0
            ; i++
        ) {
            limit--;
            if (i->isBin()) {
                if (seen[(~i->lit2()).toInt()]) {
                    stats.binTriShrinkedClause++;
                    seen[(~i->lit2()).toInt()] = 0;
                }
                continue;
            }
            break;
        }
    }
}

void Searcher::minimise_redundant_more_more(vector<Lit>& cl)
{
    /*if (conf.doStamp&& conf.more_more_with_stamp) {
        stamp_based_morem_minim(learnt_clause);
    }*/

    stats.furtherShrinkAttempt++;
    for (const Lit lit: cl) {
        seen[lit.toInt()] = 1;
    }

    if (conf.doCache && conf.more_more_with_cache) {
        cache_based_morem_minim(cl);
    }

    binary_based_morem_minim(cl);

    //Finally, remove the literals that have seen[literal] = 0
    //Here, we can count do stats, etc.
    bool changedClause  = false;
    vector<Lit>::iterator i = cl.begin();
    vector<Lit>::iterator j= i;

    //never remove the 0th literal -- TODO this is a bad thing
    //we should be able to remove this, but I can't figure out how to
    //reorder the clause then
    seen[cl[0].toInt()] = 1;
    for (vector<Lit>::iterator end = cl.end(); i != end; i++) {
        if (seen[i->toInt()]) {
            *j++ = *i;
        } else {
            changedClause = true;
        }
        seen[i->toInt()] = 0;
    }
    stats.furtherShrinkedSuccess += changedClause;
    cl.resize(cl.size() - (i-j));
}

void Searcher::stamp_based_morem_minim(vector<Lit>& cl)
{
    //Stamp-based minimization
    stats.stampShrinkAttempt++;
    const size_t origSize = cl.size();

    Lit firstLit = cl[0];
    std::pair<size_t, size_t> tmp;
    tmp = stamp.stampBasedLitRem(cl, STAMP_RED);
    if (tmp.first || tmp.second) {
        //cout << "Rem RED: " << tmp.first + tmp.second << endl;
    }
    tmp = stamp.stampBasedLitRem(cl, STAMP_IRRED);
    if (tmp.first || tmp.second) {
        //cout << "Rem IRRED: " << tmp.first + tmp.second << endl;
    }

    //Handle removal or moving of the first literal
    size_t at = std::numeric_limits<size_t>::max();
    for(size_t i = 0; i < cl.size(); i++) {
        if (cl[i] == firstLit) {
            at = i;
            break;
        }
    }
    if (at != std::numeric_limits<size_t>::max()) {
        //Make original first lit first in the final clause, too
        std::swap(cl[0], cl[at]);
    } else {
        //Re-add first lit
        cl.push_back(lit_Undef);
        std::swap(cl[0], cl[cl.size()-1]);
        cl[0] = firstLit;
    }

    stats.stampShrinkCl += ((origSize - cl.size()) > 0);
    stats.stampShrinkLit += origSize - cl.size();
}

bool Searcher::VarFilter::operator()(uint32_t var) const
{
    return (cc->value(var) == l_Undef && solver->varData[var].removed == Removed::none);
}

uint64_t Searcher::sumRestarts() const
{
    return stats.numRestarts + solver->get_stats().numRestarts;
}

size_t Searcher::hyper_bin_res_all(const bool check_for_set_values)
{
    size_t added = 0;

    for(std::set<BinaryClause>::const_iterator
        it = solver->needToAddBinClause.begin()
        , end = solver->needToAddBinClause.end()
        ; it != end
        ; ++it
    ) {
        lbool val1 = value(it->getLit1());
        lbool val2 = value(it->getLit2());

        if (conf.verbosity >= 6) {
            cout
            << "c Attached hyper-bin: "
            << it->getLit1() << "(val: " << val1 << " )"
            << ", " << it->getLit2() << "(val: " << val2 << " )"
            << endl;
        }

        //If binary is satisfied, skip
        if (check_for_set_values
            && (val1 == l_True || val2 == l_True)
        ) {
            continue;
        }

        if (check_for_set_values) {
            assert(val1 == l_Undef && val2 == l_Undef);
        }

        solver->attach_bin_clause(it->getLit1(), it->getLit2(), true, false);
        added++;
    }
    solver->needToAddBinClause.clear();

    return added;
}

#ifdef USE_GAUSS
Searcher::gauss_ret Searcher::gauss_jordan_elim()
{
    #ifdef VERBOSE_DEBUG
    cout << "Gauss searcher::Gauss_elimination called, declevel: " << decisionLevel() << endl;
    #endif
    if (gqueuedata.empty() || !solver->conf.gaussconf.enabled) {
        return gauss_ret::g_nothing;
    }

    for(uint32_t i = 0; i < gqueuedata.size(); i++) {
        auto& gqd = gqueuedata[i];
        gqd.reset();

        if (gqd.engaus_disable) {
            continue;
        }

        if (solver->conf.gaussconf.autodisable &&
            (gqd.num_entered_mtx & 0xff) == 0xff && //only check once in a while
            gqd.num_entered_mtx > 1000
        ) {
            uint32_t limit = (double)gqd.num_entered_mtx*0.01;
            uint32_t useful = 2*gqd.num_conflicts+gqd.num_props;
            if (useful < limit) {
                const double perc = stats_line_percent(gqd.num_conflicts*2+gqd.num_props, gqd.num_entered_mtx);
                if (solver->conf.verbosity) {
                    cout << "c [gauss] <" <<  i <<  "> Disabling GJ-elim in this round. "
                    " Usefulness was: "
                    << std::setprecision(2) << std::fixed << perc
                    <<  "%" << endl;
                }
                gqd.engaus_disable = true;
            }
        }
    }
    assert(qhead == trail.size());
    assert(gqhead <= qhead);

    bool confl_in_gauss = false;
    while (gqhead <  qhead
        && !confl_in_gauss
    ) {
        const Lit p = trail[gqhead++];
        assert(gwatches.size() > p.var());
        vec<GaussWatched>& ws = gwatches[p.var()];
        GaussWatched* i = ws.begin();
        GaussWatched* j = i;
        const GaussWatched* end = ws.end();
        #ifdef VERBOSE_DEBUG
        cout << "New GQHEAD: " << p << endl;
        #endif

        for (; i != end; i++) {
            if (gqueuedata[i->matrix_num].engaus_disable) {
                //remove watch and continue
                continue;
            }

            gqueuedata[i->matrix_num].enter_matrix = true;
            if (gmatrices[i->matrix_num]->find_truths2(
                i, j, p.var(), i->row_id, gqueuedata[i->matrix_num])
            ) {
                continue;
            } else {
                confl_in_gauss = true;
                i++;
                break;
            }
        }

        for (; i != end; i++) {
            *j++ = *i;
        }
        ws.shrink(i-j);

        for (size_t g = 0; g < gqueuedata.size(); g++) {
            if (gqueuedata[g].engaus_disable)
                continue;

            if (gqueuedata[g].do_eliminate) {
                gmatrices[g]->eliminate_col2(p.var(), gqueuedata[g]);
                confl_in_gauss |= (
                    gqueuedata[g].ret == gauss_res::long_confl ||
                    gqueuedata[g].ret == gauss_res::bin_confl);
            }
        }
    }

    gauss_ret finret = gauss_ret::g_nothing;
    for (GaussQData& gqd: gqueuedata) {
        if (gqd.engaus_disable)
            continue;

        if (gqd.enter_matrix) {
            gqueuedata[0].num_entered_mtx++;
            sum_gauss_entered_mtx++;
        }

        //There was a conflict but this is not that matrix.
        //Just skip.
        if (confl_in_gauss
            && gqd.ret != gauss_res::long_confl
            && gqd.ret != gauss_res::bin_confl
        ) {
            continue;
        }

        switch (gqd.ret) {
            case gauss_res::bin_confl :{
                //assert(confl.getType() == PropByType::binary_t && "this should hold, right?");
                bool ret = handle_conflict<false>(gqd.confl);
                #ifdef VERBOSE_DEBUG
                cout << "Handled binary GJ conflict"
                << " conf level:" <<  varData[gqd.confl.lit2().var()].level
                << " conf value: " << value(gqd.confl.lit2())
                << " failbin level: " << varData[solver->failBinLit.var()].level
                << " failbin value: " << value(solver->failBinLit)
                << endl;
                #endif

                gqd.num_conflicts++;
                sum_gauss_confl++;

                if (!ret) return gauss_ret::g_false;
                return gauss_ret::g_cont;
            }

            case gauss_res::long_confl :{
                gqd.num_conflicts++;
                sum_gauss_confl++;

                #ifdef DEBUG_GAUSS
                for(uint32_t i = 0; i < gqd.conflict_clause_gauss.size(); i++) {
                    Lit& l = gqd.conflict_clause_gauss[i];
                    if (value(l) != l_False) {
                        cout << "about to fail, size: " << gqd.conflict_clause_gauss.size() << " i = " << i << " val: " << value(l) << endl;
                    }
                    assert(value(l) == l_False);
                }
                #endif

                Clause* conflPtr = solver->cl_alloc.Clause_new(
                    gqd.conflict_clause_gauss,
                    sumConflicts
                );

                conflPtr->set_gauss_temp_cl();
                gqd.confl = PropBy(solver->cl_alloc.get_offset(conflPtr));
                gqhead = qhead = trail.size();

                bool ret = handle_conflict<false>(gqd.confl);
                #ifdef VERBOSE_DEBUG
                cout << "Handled long GJ conflict" << endl;
                #endif

                solver->free_cl(gqd.confl.get_offset());
                if (!ret) return gauss_ret::g_false;
                return gauss_ret::g_cont;
            }

            case gauss_res::prop:
                gqd.num_props++;
                sum_gauss_prop++;
                finret = gauss_ret::g_cont;

            case gauss_res::none:
                //nothing
                break;

            default:
                assert(false);
                return gauss_ret::g_nothing;
        }
    }
    #ifdef VERBOSE_DEBUG
    cout << "Exiting GJ" << endl;
    #endif
    return finret;
}
#endif //USE_GAUSS

std::pair<size_t, size_t> Searcher::remove_useless_bins(bool except_marked)
{
    size_t removedIrred = 0;
    size_t removedRed = 0;

    if (conf.doTransRed) {
        for(std::set<BinaryClause>::iterator
            it = uselessBin.begin()
            , end = uselessBin.end()
            ; it != end
            ; ++it
        ) {
            propStats.otfHyperTime += 2;
            if (solver->conf.verbosity >= 10) {
                cout << "Removing binary clause: " << *it << endl;
            }
            propStats.otfHyperTime += solver->watches[it->getLit1()].size()/2;
            propStats.otfHyperTime += solver->watches[it->getLit2()].size()/2;
            bool removed;
            if (except_marked) {
                bool rem1 = removeWBin_except_marked(solver->watches, it->getLit1(), it->getLit2(), it->isRed());
                bool rem2 = removeWBin_except_marked(solver->watches, it->getLit2(), it->getLit1(), it->isRed());
                assert(rem1 == rem2);
                removed = rem1;
            } else {
                removeWBin(solver->watches, it->getLit1(), it->getLit2(), it->isRed());
                removeWBin(solver->watches, it->getLit2(), it->getLit1(), it->isRed());
                removed = true;
            }

            if (!removed) {
                continue;
            }

            //Update stats
            if (it->isRed()) {
                solver->binTri.redBins--;
                removedRed++;
            } else {
                solver->binTri.irredBins--;
                removedIrred++;
            }
            *drat << del << it->getLit1() << it->getLit2() << fin;

            #ifdef VERBOSE_DEBUG_FULLPROP
            cout << "Removed bin: "
            << it->getLit1() << " , " << it->getLit2()
            << " , red: " << it->isRed() << endl;
            #endif
        }
    }
    uselessBin.clear();

    return std::make_pair(removedIrred, removedRed);
}

template<bool update_bogoprops>
PropBy Searcher::propagate() {
    const size_t origTrailSize = trail.size();

    PropBy ret;
    ret = propagate_any_order<update_bogoprops>();

    //Drat -- If declevel 0 propagation, we have to add the unitaries
    if (decisionLevel() == 0 &&
        (drat->enabled() || solver->conf.simulate_drat)
    ) {
        for(size_t i = origTrailSize; i < trail.size(); i++) {
            #ifdef DEBUG_DRAT
            if (conf.verbosity >= 6) {
                cout
                << "c 0-level enqueue:"
                << trail[i]
                << endl;
            }
            #endif
            *drat << add << trail[i]
            << fin;
        }
        if (!ret.isNULL()) {
            *drat << add
            << fin;
        }
    }

    return ret;
}
template PropBy Searcher::propagate<true>();
template PropBy Searcher::propagate<false>();

size_t Searcher::mem_used() const
{
    size_t mem = HyperEngine::mem_used();
    mem += otf_subsuming_short_cls.capacity()*sizeof(OTFClause);
    mem += otf_subsuming_long_cls.capacity()*sizeof(ClOffset);
    mem += var_act_vsids.capacity()*sizeof(uint32_t);
    mem += var_act_maple.capacity()*sizeof(uint32_t);
    mem += order_heap_vsids.mem_used();
    mem += learnt_clause.capacity()*sizeof(Lit);
    mem += hist.mem_used();
    mem += conflict.capacity()*sizeof(Lit);
    mem += model.capacity()*sizeof(lbool);
    mem += analyze_stack.mem_used();
    mem += assumptions.capacity()*sizeof(Lit);

    if (conf.verbosity >= 3) {
        cout
        << "c otfMustAttach bytes: "
        << otf_subsuming_short_cls.capacity()*sizeof(OTFClause)
        << endl;

        cout
        << "c toAttachLater bytes: "
        << otf_subsuming_long_cls.capacity()*sizeof(ClOffset)
        << endl;

        cout
        << "c toclear bytes: "
        << toClear.capacity()*sizeof(Lit)
        << endl;

        cout
        << "c trail bytes: "
        << trail.capacity()*sizeof(Lit)
        << endl;

        cout
        << "c trail_lim bytes: "
        << trail_lim.capacity()*sizeof(Lit)
        << endl;

        cout
        << "c order_heap_vsids bytes: "
        << order_heap_vsids.mem_used()
        << endl;

        cout
        << "c learnt clause bytes: "
        << learnt_clause.capacity()*sizeof(Lit)
        << endl;

        cout
        << "c hist bytes: "
        << hist.mem_used()
        << endl;

        cout
        << "c conflict bytes: "
        << conflict.capacity()*sizeof(Lit)
        << endl;

        cout
        << "c Stack bytes: "
        << analyze_stack.capacity()*sizeof(Lit)
        << endl;
    }

    return mem;
}

void Searcher::fill_assumptions_set()
{
    #ifdef SLOW_DEBUG
    for(auto x: varData) {
        assert(x.assumption == l_Undef);
    }
    #endif

    for(const AssumptionPair lit_pair: assumptions) {
        const Lit lit = map_outer_to_inter(lit_pair.lit_outer);
        varData[lit.var()].assumption = lit.sign() ? l_False : l_True;
    }
}

void Searcher::unfill_assumptions_set()
{
    for(const AssumptionPair lit_pair: assumptions) {
        const Lit lit = map_outer_to_inter(lit_pair.lit_outer);
        varData[lit.var()].assumption = l_Undef;
    }

    #ifdef SLOW_DEBUG
    for(auto x: varData) {
        assert(x.assumption == l_Undef);
    }
    #endif
}

inline void Searcher::varDecayActivity()
{
    var_inc_vsids *= (1.0 / var_decay_vsids);
}

void Searcher::update_var_decay_vsids()
{
    if (var_decay_vsids >= conf.var_decay_vsids_max) {
        var_decay_vsids = conf.var_decay_vsids_max;
    }
}

void Searcher::consolidate_watches(const bool full)
{
    double t = cpuTime();
    if (full) {
        watches.full_consolidate();
    } else {
        watches.consolidate();
    }
    double time_used = cpuTime() - t;

    if (conf.verbosity) {
        cout
        << "c [consolidate] "
        << (full ? "full" : "mini")
        << conf.print_times(time_used)
        << endl;
    }

    std::stringstream ss;
    ss << "consolidate " << (full ? "full" : "mini") << " watches";
}

//Normal running
template
void Searcher::cancelUntil<true, false>(uint32_t level, uint32_t clid_plus);

//During inprocessing, dont update anyting really (probing, distilling)
template
void Searcher::cancelUntil<false, true>(uint32_t level, uint32_t clid_plus);

template<bool do_insert_var_order, bool update_bogoprops>
void Searcher::cancelUntil(uint32_t level
    , uint32_t
) {
    #ifdef VERBOSE_DEBUG
    cout << "Canceling until level " << level;
    if (level > 0) cout << " sublevel: " << trail_lim[level];
    cout << endl;
    #endif

    if (decisionLevel() > level) {
        #ifdef USE_GAUSS
        for (EGaussian* gauss: gmatrices)
            if (gauss) {
                gauss->canceling(trail_lim[level]);
            }
        #endif //USE_GAUSS

        //Go through in reverse order, unassign & insert then
        //back to the vars to be branched upon
        for (int sublevel = trail.size()-1
            ; sublevel >= (int)trail_lim[level]
            ; sublevel--
        ) {
            #ifdef VERBOSE_DEBUG
            cout
            << "Canceling lit " << trail[sublevel]
            << " sublevel: " << sublevel
            << endl;
            #endif

            const uint32_t var = trail[sublevel].var();
            assert(value(var) != l_Undef);

            assigns[var] = l_Undef;
            if (do_insert_var_order) {
                insert_var_order(var);
            }
        }
        qhead = trail_lim[level];
        trail.resize(trail_lim[level]);
        trail_lim.resize(level);
    }

    #ifdef VERBOSE_DEBUG
    cout
    << "Canceling finished. Now at level: " << decisionLevel()
    << " sublevel: " << trail.size()-1
    << endl;
    #endif
}

inline bool Searcher::check_order_heap_sanity() const
{
    if (conf.sampling_vars) {
        for(uint32_t outside_var: *conf.sampling_vars) {
            uint32_t outer_var = map_to_with_bva(outside_var);
            outer_var = solver->varReplacer->get_var_replaced_with_outer(outer_var);
            uint32_t int_var = map_outer_to_inter(outer_var);
            assert(varData[int_var].removed == Removed::none);

            if (int_var < nVars() &&
                varData[int_var].removed == Removed::none &&
                value(int_var) == l_Undef
            ) {
                assert(order_heap_vsids.inHeap(int_var));
            }
        }
    }

    for(size_t i = 0; i < nVars(); i++)
    {
        if (varData[i].removed == Removed::none
            && value(i) == l_Undef)
        {
            if (!order_heap_vsids.inHeap(i)) {
                cout << "ERROR var " << i+1 << " not in VSIDS heap."
                << " value: " << value(i)
                << " removed: " << removed_type_to_string(varData[i].removed)
                << endl;
                return false;
            }
        }
    }
    assert(order_heap_vsids.heap_property());

    return true;
}

#ifdef USE_GAUSS
void Searcher::clear_gauss_matrices()
{
    for(uint32_t i = 0; i < gqueuedata.size(); i++) {
        auto gqd = gqueuedata[i];
        if (solver->conf.verbosity && gqd.num_entered_mtx > 0) {
            cout << "c [gauss] < " << i << " > "
            << "entered mtx    : " << std::left << print_value_kilo_mega(gqd.num_entered_mtx, false)
            << endl;

            cout << "c [gauss] < " << i << " > "
            << "confl triggered: "
            << stats_line_percent(gqd.num_conflicts, gqd.num_entered_mtx)
            << " %" <<endl;

            cout << "c [gauss] < " << i << " > "
            << "prop  triggered: "
            << stats_line_percent(gqd.num_props, gqd.num_entered_mtx)
            << " %" <<endl;
        }

        if (solver->conf.verbosity >= 3 && gqd.num_entered_mtx > 0) {
            cout
            << "c [gauss] num_props       : "<< print_value_kilo_mega(gqd.num_props) << endl
            << "c [gauss] num_conflicts   : "<< print_value_kilo_mega(gqd.num_conflicts)  << endl;
        }
        gqd.reset_stats();
    }

    if (solver->conf.verbosity >= 3 && sum_gauss_entered_mtx > 0) {
        cout
        << "c [gauss] sum_gauss_prop: " << print_value_kilo_mega(sum_gauss_prop) << endl
        << "c [gauss] sum_gauss_confl : " << print_value_kilo_mega(sum_gauss_confl) << endl
        << "c [gauss] sum_gauss_entered_mtx    : " << print_value_kilo_mega(sum_gauss_entered_mtx) << endl;
    }

    //cout << "Clearing matrices" << endl;
    for(EGaussian* g: gmatrices) {
        delete g;
    }
    for(auto& w: gwatches) {
        w.clear();
    }
    gmatrices.clear();
    gqueuedata.clear();
    xor_clauses_updated = true;
}
#endif

void Searcher::check_assumptions_sanity()
{
    for(AssumptionPair& lit_pair: assumptions) {
        Lit inter_lit = map_outer_to_inter(lit_pair.lit_outer);
        assert(inter_lit.var() < varData.size());
        assert(varData[inter_lit.var()].removed == Removed::none);
        if (varData[inter_lit.var()].assumption == l_Undef) {
            cout << "Assump " << inter_lit << " has .assumption : "
            << varData[inter_lit.var()].assumption << endl;
        }
        assert(varData[inter_lit.var()].assumption != l_Undef);
    }
}
