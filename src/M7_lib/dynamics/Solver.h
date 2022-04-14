//
// Created by rja on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/util/Timer.h>
#include <M7_lib/observables/RefExcits.h>
#include <M7_lib/observables/UniformTwf.h>
#include <M7_lib/observables/RefExcits.h>
#include <M7_lib/bilinear/Bilinears.h>
#include <M7_lib/io/FciqmcStats.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/io/TimingStats.h>
#include <M7_lib/io/FciqmcStats.h>
#include <M7_lib/io/ParallelStats.h>
#include <M7_lib/mae/Maes.h>
#include <M7_lib/wavefunction/Reference.h>
#include <M7_lib/wavefunction/DeterministicSubspace.h>

#include "Propagator.h"
#include "Annihilator.h"

/**
 * This central class brings together wavefunctions, propagator, expectation values, and statistics.
 *
 * The main loop of the projector method is carried out in the execute method, and each iteration of the loop is
 * referred to as a cycle.
 *
 * Many subordinate objects have a state which is cyclic in nature, and so confusion can arise in recognising at any
 * given point in the algorithm which quantities correspond to the "last", "current", and "next" cycles.
 *
 * The basic structure of a Solver::execute cycle is as follows
 * call begin_cycle
 * do WF(i+1) = Prop(i) x WF(i)
 * call end_cycle
 * write stats
 *
 * begin_cycle and end_cycle defined here call in turn the begin_cycle and end_cycle methods of subordinate objects. One
 * might look to put all resets of cyclical variables in a single method, but it is actually necessary to split these
 * procedures up. E.g. the number of walkers must be zeroed in begin_cycle, before being accumulated in the loop over
 * occupied MBFs, and then being sum-reduced over all MPI ranks. This reduced value is required by other modules
 * including notably the stats writer.
 *
 * begin_cycle should be the more logically intricate method, handling:
 *  - updates to variables derived from cyclical (reducable) variables
 *  - cyclical variable rezeroing
 *  - interactive variable changes
 *
 * where as end_cycle should be relatively lightweight, performing only the reductions themselves, allowing access to
 * the reduced quantities in the update method call in the next cycle before rezeroing.
 */
class Solver {

    size_t m_icycle = 0ul;
    Propagator &m_prop;
    const fciqmc_config::Document &m_opts;
    Wavefunction &m_wf;
    References m_refs;

    std::unique_ptr<FciqmcStats> m_stats = nullptr;
    std::unique_ptr<TimingStats> m_timing_stats = nullptr;
    std::unique_ptr<ParallelStats> m_parallel_stats = nullptr;

    /*
     * Timers for the main parts of the solver
     */
    // whole cycle
    Timer m_cycle_timer;
    // whole loop over occupied rows
    Timer m_propagate_timer;
    // individual iterations over an occupied row
    Timer m_spawning_timer;
    // time waited at MPI_Barrier
    Timer m_synchronization_timer;
    // time taken to complete communication of spawned walkers
    Timer m_communicate_timer;
    // time taken to complete whole annihilation loop
    Timer m_annihilate_timer;

    /*
     * Sanity checking variables
     */
    defs::wf_t m_chk_nwalker_local = 0.0;
    size_t m_chk_ninitiator_local = 0ul;

    InteractiveVariable<bool> m_exit;
public:
    Maes m_maes;

private:
    Annihilator m_annihilator;
    Archive m_archive;
    DeterministicSubspaces m_detsubs;

public:

    Solver(const fciqmc_config::Document& opts, Propagator &prop, Wavefunction &wf, std::vector<TableBase::Loc> ref_locs);

    Solver(const fciqmc_config::Document& opts, Propagator &prop, Wavefunction &wf, TableBase::Loc ref_loc):
        Solver(opts, prop, wf, std::vector<TableBase::Loc>(wf.npart(), ref_loc)){}

    /**
     * Perform ncycle iterations of the solver.
     *
     * Each iteration of the main loop advances the m_icycle counter. This method delegates to the main subalgorithms in
     * the FCIQMC scheme. Since iterative updates to data structures occur at various points in the body of this main
     * loop, it is of paramount importance that the state of these data structures is clearly defined with respect to
     * the action of the different subroutines invoked.
     *
     * Stats are numeric quantities associated with both the wavefunction and its propagation. Those stats pertaining to
     * the instantaneous state of the wavefunction are termed "state" stats, whereas those pertaining to propagation, or
     * changes in the wavefunction under the action of projectors, are called "delta" stats
     *
     * m_wf represents the initial wavefunction Psi_0
     * for each icycle in range [0, ncycle)
     *  // stats are zeroed.
     *  // wavefunction represents Psi_icycle
     *  loop over occupied MBFs
     *  // wavefunction represents intermediate state between Psi_icycle and Psi_(icycle+1)
     *  // since the death/cloning step has been applied but not the annihilation
     *  // local "state" stats reflect the wavefunction Psi_icycle
     *  // local "delta" stats partially reflect the change from Psi_icycle to Psi_(icycle+1)
     *  communicate spawned walkers
     *  loop over spawned walkers
     *  // wavefunction represents Psi_(icycle+1)
     *  // local "delta" stats partially reflect the change from Psi_icycle to Psi_(icycle+1)
     *  // all local stats are reduced to global stats
     * @param ncycle
     */
    void execute(size_t ncycle = 1);

    /**
     * reset variables and those of member objects for a new solver iteration
     */
    void begin_cycle();

    /**
     * Perform diagonal and off-diagonal propagation via m_prop for a valid row of m_wf.m_store. The row being
     * propagated from is the one currently pointed to by m_wf.m_store.m_row
     * @param ipart
     *  flat index of m_wf.m_format being propagated
     */
    void propagate_row(const size_t &ipart);

    /**
     * Loop over all rows in m_wf.m_store which have a non-zero MBF field
     */
    void loop_over_occupied_mbfs();

    /**
     * Loop over all rows in m_wf.m_store which have a non-zero MBF field but perform no propagation, just add any
     * required weight-averaged contributions to the MEVs
     */
    void finalizing_loop_over_occupied_mbfs(size_t icycle);

    void loop_over_spawned();

    void end_cycle();

    void output_stats();

};

#endif //M7_SOLVER_H
