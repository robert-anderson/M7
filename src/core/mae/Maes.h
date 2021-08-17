//
// Created by rja on 15/08/2021.
//

#ifndef M7_MAES_H
#define M7_MAES_H

#include "src/core/bilinear/Bilinears.h"
#include "src/core/observables/RefExcits.h"

struct Maes {
    Epoch m_accum_epoch;
    Bilinears m_bilinears;
    RefExcits m_ref_excits;
    const size_t m_period;
    size_t m_icycle_period_start = ~0ul;
    /**
     * work space
     */
    conn::Mbf m_conn_work;
    FrmOps m_com_work;

    Maes(const fciqmc_config::AvEsts &opts, const Hamiltonian &ham) :
            m_accum_epoch("MAE accumulation"), m_bilinears(opts, ham, m_accum_epoch),
            m_ref_excits(opts.m_ref_excits, ham.nsite()), m_period(opts.m_stats_period),
            m_conn_work(ham.nsite()), m_com_work(ham.nsite()) {}

    operator bool() const {
        return m_bilinears || m_ref_excits;
    }

    bool all_stores_empty() const {
        return m_bilinears.all_stores_empty() && m_ref_excits.all_stores_empty();
    }

    bool is_period_cycle(size_t icycle) {
        if (!m_accum_epoch) return false;
        if (!m_period) return false;
        if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
            m_icycle_period_start = icycle;
            return false;
        }
        return !((icycle - m_icycle_period_start) % m_period);
    }

    void end_cycle() {
        m_bilinears.end_cycle();
    }

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
    void make_average_contribs(WalkerTableRow &row, const References &refs, const size_t &icycle) {
        if (!m_accum_epoch) return;
        // the current cycle should be included in the denominator
        if (!row.occupied_ncycle(icycle)) {
            DEBUG_ASSERT_TRUE(row.m_average_weight.is_zero(), "average value should have been rezeroed");
            return;
        }
        defs::wf_comp_t ncycle_occ = row.occupied_ncycle(icycle);

        for (size_t ipart = 0ul; ipart < row.m_wf_format.m_nelement; ++ipart) {
            auto &ref = refs[ipart];
            auto &ref_mbf = ref.get_mbf();
            auto ipart_replica = row.ipart_replica(ipart);
            /*
             * if contributions are coming from two replicas, we should take the mean
             */
            double dupl_fac = (ipart_replica == ipart) ? 1.0 : 0.5;
            /*
             * the "average" weights actually refer to the unnormalized average. The averages are obtained by dividing
             * each by the number of cycles for which the row is occupied.
             */
            const auto av_weight = row.m_average_weight[ipart] / ncycle_occ;

            m_conn_work.connect(ref_mbf, row.m_mbf, m_com_work);
            auto exsig_from_ref = m_conn_work.exsig();
            auto is_ref_conn = !exsig_from_ref && m_bilinears.m_rdms.takes_contribs_from(exsig_from_ref);

            /*
             * accumulate contributions to reference excitations if required
             */
            m_ref_excits.make_contribs(m_conn_work, dupl_fac * ncycle_occ * av_weight, ipart);

            if (m_bilinears.m_rdms) {
                auto av_weight_rep = row.m_average_weight[ipart_replica] / ncycle_occ;
                /*
                 * scale up the product by a factor of the number of instantaneous contributions being accounted for in this
                 * single averaged contribution (ncycle_occ)
                 */
                m_bilinears.make_contribs(row.m_mbf, dupl_fac * ncycle_occ * av_weight * av_weight_rep);
                if (m_bilinears.m_rdms.m_explicit_ref_conns && is_ref_conn) {
                    const auto av_weight_ref = ref.norm_average_weight(icycle, ipart);
                    const auto av_weight_ref_rep = ref.norm_average_weight(icycle, ipart_replica);
                    m_bilinears.m_rdms.make_contribs(ref_mbf, row.m_mbf,
                                                     dupl_fac * ncycle_occ * av_weight_ref * av_weight_rep);
                    m_bilinears.m_rdms.make_contribs(row.m_mbf, ref_mbf,
                                                     dupl_fac * ncycle_occ * av_weight * av_weight_ref_rep);
                }
                row.m_average_weight = 0;
                row.m_icycle_occ = icycle;
            }
        }
    }
};


#endif //M7_MAES_H
