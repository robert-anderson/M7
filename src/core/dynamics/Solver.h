//
// Created by rja on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/util/Timer.h>
#include <src/core/observables/AverageCoefficients.h>
#include <src/core/observables/UniformTwf.h>
#include <src/core/observables/WeightedTwf.h>
#include <src/core/observables/MevGroup.h>
#include <src/core/io/FciqmcStats.h>
#include "Reference.h"
#include "src/core/table/Communicator.h"
#include "src/core/io/FciqmcStats.h"
#include "src/core/io/ParallelStats.h"
#include "Propagator.h"
#include "DeterministicSubspace.h"

class Solver {

    size_t m_icycle = 0ul;
    Propagator &m_prop;
    const Options &m_opts;
    Wavefunction &m_wf;
    References m_refs;

    std::unique_ptr<FciqmcStats> m_stats = nullptr;
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
    MevGroup m_mevs;

    DeterministicSubspace2 m_detsub;

public:

    Solver(Propagator &prop, Wavefunction &wf, std::vector<TableBase::Loc> ref_locs);

    Solver(Propagator &prop, Wavefunction &wf, TableBase::Loc ref_loc): Solver(prop, wf, std::vector<TableBase::Loc>(wf.npart(), ref_loc)){}

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
     *  loop over occupied ONVs
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
     * Loop over all rows in m_wf.m_store which have a non-zero ONV field
     */
    void loop_over_occupied_onvs();

    /**
     * Loop over all rows in m_wf.m_store which have a non-zero ONV field but perform no propagation, just add any
     * required weight-averaged contributions to the MEVs
     */
    void finalizing_loop_over_occupied_onvs();

    void annihilate_row(const size_t &dst_ipart, const fields::Onv<> &dst_onv, const defs::wf_t &delta_weight,
                        bool allow_initiation, bool src_deterministic, const size_t &irow_store);

    void annihilate_row(const size_t &dst_ipart, const fields::Onv<> &dst_onv, const defs::wf_t &delta_weight,
                        bool allow_initiation, bool src_deterministic);

    /**
     * Make all contributions to MEVs from the current occupied ONV row. These contributions always include the
     * diagonals, where the bra and ket ONVs are the same. Explicit contributions from connections to the Hartree-Fock
     * ONV are also optionally included here - in that case it is taken to be true that the single excitations are never
     * generated due to the Brillouin theorem
     *
     * Care is needed here to avoid off-by-one-like errors. Such considerations are required in a few different methods,
     * but all relevant details are summarised here.
     *
     * Five different types of event must be considered in the proper accumulation of averaged contributions:
     *  1. a row is created in the walker table
     *  2. the mev accumulation epoch begins
     *  3. a row is about to be removed from the walker table
     *  4. the boundary between two block averaging periods is reached
     *  5. the end of the calculation is reached (equivalent to 3. for every row in the table)
     *
     * these are each handled in the following manner:
     *  1. if row creation is performed in the annhilation step of cycle i, at least one component of its weight array
     *     will be set to a non-zero value immediately afterwards in the same cycle. the next cycle (i+1) will be
     *     treated as the first cycle of the lifetime of this new row. thus, on creation of the row, the m_icycle_occ
     *     member of the row will be set to i, and the m_average_weight will be zeroed. then, if contributions were
     *     "averaged" every cycle (i.e. m_opts.ncycle_mev_period=1), the newly added instantaneous row.m_weight would be
     *     summed into row.m_average_weight and a call to row.occupied_ncycle(m_icycle) would return 1, since m_icycle
     *     would be incremented to i+1. therefore, the added contribution would be correct. assuming that an arbitrary
     *     element of the weight then remains unchanged, on cycle i+x row.m_average_weight would be
     *     x * row.m_weight, and row.occupied_ncycle(m_icycle) would give x, the correct normalization
     *
     *  2. when the accumulation epoch begins while the wavefunction already contains an occupied set, each row must be
     *     treated as though it were created on the previous iteration. In the loop over occupied ONVs, MC cycle i, if
     *     the accumulation epoch has begun at the beginning of cycle i, then row.m_icycle_occ must be set to i-1, and
     *     the average zeroed. Then, immediately afterwards the cyclic summation of row.m_weight into
     *     row.m_average_weight would occur, and the same sanity checks as described in 1. would pass, since the number
     *     of occupied cycles for m_opts.ncycle_mev_period=1 would evaluate to 1, and that's exactly the number of
     *     past row.m_weights that have been summed into row.m_average_weights
     *
     *  3. if the conditions have been met for a row to be removed, namely all elements of the weight have become zero,
     *     then it is necessary that immediately prior to its deletion from m_wf.m_store, any contributions it owes to
     *     MEV elements which take contributions from products of averaged weights are made. if a row were added to
     *     m_wf.m_store in the annihilation loop of cycle i, and then removed in the loop over occupied ONVs of cycle
     *     i+1, then the number of occupied cycles would evaluate to 1, but the average weight of the row would not
     *     have yet received any contributions, and so the contributions to MEVs would be zero, correctly. If however
     *     the row survived for one more cycle before meeting the criteria for deletion, the initially added weight
     *     would have contributed once, but not the weight on the iteration of deletion. However, this weight is
     *     necessarily zero, and so the contribution would be correct for row.occupied_ncycle(m_icycle)=2.
     *
     *  4. here, the row is treated as though it becomes unoccupied on cycle i, with its average value zeroed. thence,
     *     in the loop over occupied ONVs of cycle i+1, the instantaneous weight is summed in and the normalization is
     *     correct.
     *
     *  5. if the calculation ends and the number of iterations contributing to the unnormalized coefficient averages
     *     is not an integral multiple of the block averaging period, the contributions owed to the MEV estimates must
     *     be added in a special "finalizing" loop over occupied ONVs. crucially, this is done *before* the
     *     instantaneous weight is summed into the average, since this was already done in the previous iteration.
     */
    void make_average_weight_mev_contribs(const size_t& icycle);

    /**
     * Make all MEV contributions due to products of instantaneous walker weights.
     * @param src_onv
     *  source ONV - from which the spawned walker originated
     * @param src_weight
     *  source weight - the number of walkers on the source ONV
     * @param dst_ipart
     *  the WF part index for which this contribution is bound
     */
    void make_instant_mev_contribs(const fields::Onv<> &src_onv, const defs::wf_t &src_weight, const size_t &dst_ipart);

    void make_mev_contribs_from_unique_src_onvs(SpawnTableRow &row_current, SpawnTableRow &row_block_start,
                                                const size_t &irow_block_end, const size_t &irow_store);

    void loop_over_spawned();

    void end_cycle();

    void output_mevs();

    void output_stats();

    const MevGroup& mevs() const;
};

#endif //M7_SOLVER_H
