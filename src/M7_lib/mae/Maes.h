//
// Created by rja on 15/08/2021.
//

#ifndef M7_MAES_H
#define M7_MAES_H

#include <M7_lib/io/MaeStats.h>
#include <M7_lib/bilinear/Bilinears.h>
#include <M7_lib/observables/RefExcits.h>

/**
 * A group of all Multidimensional Averaging Estimators needed for a calculation
 */
struct Maes {
    /**
     * changes state on the cycle on which accumulation (making contributions) begins
     */
    Epoch m_accum_epoch;
    /**
     * group of MAEs that are bilinear in the wavefunction (e.g. RDMs, MPRT2 intermediates, spectral moments)
     */
    Bilinears m_bilinears;
    /**
     * averaged amplitudes of MBFs that are excitations of the reference MBF
     */
    RefExcits m_ref_excits;
    /**
     * number of cycles between consecutive averaging operations and stats output
     */
    const size_t m_period;
    /**
     * cycle on which the current period started
     */
    size_t m_icycle_period_start = ~0ul;
    /**
     * stats output object
     */
    std::unique_ptr<MaeStats> m_stats = nullptr;

    Maes(const conf::AvEsts &opts, sys::Size extents, size_t nelec, size_t nroot);

    operator bool() const;

    bool all_stores_empty() const;

    bool is_period_cycle(size_t icycle);

    void end_cycle();

    /**
     * Make all contributions to MAEs from the current occupied MBF row.
     *
     * Currently MAEs consist of the average reference excitations and the bilinears, of which the latter entail the
     * most careful handling.
     *
     * Averaged contributions to the bilinear MAEs always include the diagonals, where the bra and ket MBFs are the
     * same. Explicit contributions from connections to the Hartree-Fock MBF are also optionally included here - in that
     * case it is taken to be true that the 1100 excitations are never generated due to the Brillouin theorem
     *
     * Care is needed here to avoid off-by-one-like errors. Such considerations are required in a few different methods,
     * but all relevant details are summarised here.
     *
     * Five different types of event must be considered in the proper accumulation of averaged contributions:
     *  1. a row is created in the walker table
     *  2. the bilinear MAE accumulation epoch begins
     *  3. a row is about to be removed from the walker table
     *  4. the boundary between two block averaging periods is reached
     *  5. the end of the calculation is reached (equivalent to #3 for every row in the table)
     *
     * these are each handled in the following manner:
     *  1. if row creation is performed in the annhilation step of cycle i, at least one component of its weight array
     *     will be set to a non-zero value immediately afterwards in the same cycle. the next cycle (i+1) will be
     *     treated as the first cycle of the lifetime of this new row. thus, on creation of the row, the m_icycle_occ
     *     member of the row will be set to i, and the m_average_weight will be zeroed. then, if contributions were
     *     "averaged" every cycle the newly added instantaneous row.m_weight would be summed into row.m_average_weight
     *     and a call to row.occupied_ncycle(m_icycle) would return 1, since m_icycle would be incremented to i+1.
     *     therefore, the added contribution would be correct. assuming that an arbitrary element of the weight then
     *     remains unchanged, on cycle i+x row.m_average_weight would be x * row.m_weight, and
     *     row.occupied_ncycle(m_icycle) would give x, the correct normalization
     *
     *  2. when the accumulation epoch begins while the wavefunction already contains an occupied set, each row must be
     *     treated as though it were created on the previous iteration. In the loop over occupied MBFs, MC cycle i, if
     *     the accumulation epoch has begun at the beginning of cycle i, then row.m_icycle_occ must be set to i-1, and
     *     the average zeroed. Then, immediately afterwards the cyclic summation of row.m_weight into
     *     row.m_average_weight would occur, and the same sanity checks as described in 1. would pass, since the number
     *     of occupied cycles for m_opts.ncycle_mev_period=1 would evaluate to 1, and that's exactly the number of
     *     past row.m_weights that have been summed into row.m_average_weights
     *
     *  3. if the conditions have been met for a row to be removed, namely all elements of the weight have become zero,
     *     then it is necessary that immediately prior to its deletion from m_wf.m_store, any contributions it owes to
     *     MEV elements which take contributions from products of averaged weights are made. if a row were added to
     *     m_wf.m_store in the annihilation loop of cycle i, and then removed in the loop over occupied MBFs of cycle
     *     i+1, then the number of occupied cycles would evaluate to 1, but the average weight of the row would not
     *     have yet received any contributions, and so the contributions to MEVs would be zero, correctly. If however
     *     the row survived for one more cycle before meeting the criteria for deletion, the initially added weight
     *     would have contributed once, but not the weight on the iteration of deletion. However, this weight is
     *     necessarily zero, and so the contribution would be correct for row.occupied_ncycle(m_icycle)=2.
     *
     *  4. here, the row is treated as though it becomes unoccupied on cycle i, with its average value zeroed. thence,
     *     in the loop over occupied MBFs of cycle i+1, the instantaneous weight is summed in and the normalization is
     *     correct.
     *
     *  5. if the calculation ends and the number of iterations contributing to the unnormalized coefficient averages
     *     is not an integral multiple of the block averaging period, the contributions owed to the MEV estimates must
     *     be added in a special "finalizing" loop over occupied MBFs. crucially, this is done *before* the
     *     instantaneous weight is summed into the average, since this was already done in the previous iteration.
     */
    void make_average_contribs(WalkerTableRow &row, const References &refs, const size_t &icycle);

    void output(size_t icycle, const Hamiltonian& ham, bool final=false);
};


#endif //M7_MAES_H
