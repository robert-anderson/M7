//
// Created by rja on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/util/Timer.h>
#include <src/core/observables/RefExcits.h>
#include <src/core/observables/UniformTwf.h>
#include <src/core/observables/WeightedTwf.h>
#include <src/core/observables/RefExcits.h>
#include <src/core/bilinear/Bilinears.h>
#include <src/core/io/FciqmcStats.h>
#include <src/core/io/Archivable.h>
#include <src/core/io/TimingStats.h>
#include <src/core/mae/Maes.h>
#include "src/core/wavefunction/Reference.h"
#include "src/core/io/FciqmcStats.h"
#include "src/core/io/ParallelStats.h"
#include "Propagator.h"
#include "src/core/wavefunction/DeterministicSubspace.h"
#include "Annihilator.h"

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
    std::unique_ptr<UniformTwf> m_uniform_twf;
    std::unique_ptr<WeightedTwf> m_weighted_twf;
    Maes m_maes;
    Annihilator m_annihilator;
    Archive m_archive;

    std::unique_ptr<DeterministicSubspace> m_detsub = nullptr;

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
    void finalizing_loop_over_occupied_mbfs();

    void annihilate_row(const size_t &dst_ipart, const field::Mbf &dst_mbf, const defs::wf_t &delta_weight,
                        bool allow_initiation, bool src_deterministic, const size_t &irow_store);

    void annihilate_row(const size_t &dst_ipart, const field::Mbf &dst_mbf, const defs::wf_t &delta_weight,
                        bool allow_initiation, bool src_deterministic);

    void make_mev_contribs_from_unique_src_mbfs(SpawnTableRow &row_current, SpawnTableRow &row_block_start,
                                                const size_t &irow_block_end, const size_t &irow_store);

    void loop_over_spawned();

    void end_cycle();

    void output_mevs(size_t icycle);

    void output_mevs();

    void output_stats();

};

#endif //M7_SOLVER_H
