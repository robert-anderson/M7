//
// Created by Robert J. Anderson on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/util/Timer.h>
#include <M7_lib/observables/UniformTwf.h>
#include <M7_lib/observables/HfExcits.h>
#include <M7_lib/bilinear/Bilinears.h>
#include <M7_lib/io/FciqmcStats.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/io/TimingStats.h>
#include <M7_lib/io/FciqmcStats.h>
#include <M7_lib/io/ParallelStats.h>
#include <M7_lib/mae/Maes.h>
#include <M7_lib/wavefunction/Reference.h>
#include <M7_lib/wavefunction/Wavefunction.h>
#include <M7_lib/wavefunction/DeterministicSubspace.h>

#include "M7_lib/propagator/Propagator.h"
#include "M7_lib/observables/InstEsts.h"
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
    /**
     * current cycle number
     */
    uint_t m_icycle = 0ul;
    /**
     * configuration document read from the YAML file provided as the command line argument
     */
    const conf::Document &m_opts;
    /**
     * generates the diagonal and off-diagonal updates to the solution vector due to application of the integrator
     */
    Propagator &m_prop;
    /**
     * solution vector storing multiple eigenvectors of H in distributed memory
     */
    wf::Fci &m_wf;
    /**
     * reference many-body basis functions (MBFs)
     */
    wf::Refs m_refs;
    /**
     * Hartree-Fock basis function, allocated only when H satisfies the Brillouin theorem wrt the initial reference
     */
    const std::unique_ptr<shared_rows::Walker> m_hf = nullptr;
    /**
     * statistics relating to the propagation of the walker population
     */
    std::unique_ptr<FciqmcStats> m_stats = nullptr;
    /**
     * statistics relating to the walltime of major operations
     */
    std::unique_ptr<TimingStats> m_timing_stats = nullptr;
    /**
     * statistics relating to parallelization and load balancing
     */
    std::unique_ptr<ParallelStats> m_parallel_stats = nullptr;

    /*
     * Timers for the main parts of the solver
     */
    /**
     * timer for a whole FCIQMC cycle
     */
    Timer m_cycle_timer;
    /**
     * timer for a whole loop over occupied rows
     */
    Timer m_propagate_timer;
    /**
     * timer for individual iterations over an occupied row
     */
    Timer m_spawning_timer;
    /**
     * time waited at MPI_Barrier
     */
    Timer m_synchronization_timer;
    /**
     * time taken to complete communication of spawned walkers
     */
    Timer m_communicate_timer;
    /**
     * time taken to complete whole annihilation loop
     */
    Timer m_annihilate_timer;

    /**
     * Sanity checking variables
     */
    wf_t m_chk_nwalker_local = 0.0;

    /**
     * listens for a file requesting a "soft exit"
     */
    InteractiveVariable<bool> m_exit;
public:
    /**
     * Multidimensional averaging estimators i.e. reference connections, RDMs, spectral moments
     */
    Maes m_maes;

    /**
     * Instantaneous estimators. Optional stats that are linear in the wavefunction (e.g. S^2 operator)
     */
    InstEsts m_inst_ests;

private:
    /**
     * instance of helper class to handle annihilation of received spawns with each other and the existing walkers
     */
    Annihilator m_annihilator;
    /**
     * selections of MBF in which semi-stochastic propagation is performed
     */
    DeterministicSubspaces m_detsubs;

    std::unique_ptr<shared_rows::Walker> make_hf() const;

public:

    Solver(const conf::Document& opts, Propagator &prop, wf::Fci &wf, v_t<TableBase::Loc> ref_locs);

    Solver(const conf::Document& opts, Propagator &prop, wf::Fci &wf, TableBase::Loc ref_loc):
        Solver(opts, prop, wf, v_t<TableBase::Loc>(wf.npart(), ref_loc)){}

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
    void execute(uint_t ncycle = 1);

    /**
     * reset variables and those of member objects for a new solver iteration
     */
    void begin_cycle();

    /**
     * Perform diagonal and off-diagonal propagation via m_prop for a valid row of m_wf.m_store. The row being
     * propagated from is the one currently pointed to by m_wf.m_store.m_row
     * @param walker
     *
     * @param ipart
     *  flat index of m_wf.m_format being propagated
     */
    void propagate_row(Walker& walker, uint_t ipart);

    /**
     * Loop over all rows in m_wf.m_store which have a non-zero MBF field
     */
    void loop_over_occupied_mbfs();

    /**
     * Loop over all rows in m_wf.m_store which have a non-zero MBF field but perform no propagation, just add any
     * required weight-averaged contributions to the MEVs
     */
    void finalizing_loop_over_occupied_mbfs(uint_t icycle);

    void loop_over_spawned();

    void end_cycle();

    void output_stats();

};

#endif //M7_SOLVER_H
