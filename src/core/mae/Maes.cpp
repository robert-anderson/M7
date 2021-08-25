//
// Created by rja on 15/08/2021.
//

#include "Maes.h"

Maes::Maes(const fciqmc_config::AvEsts &opts, size_t nsite, size_t nelec) :
        m_accum_epoch("MAE accumulation"), m_bilinears(opts, nsite, nelec, m_accum_epoch),
        m_ref_excits(opts.m_ref_excits, nsite), m_period(opts.m_stats_period),
        m_conn_work(nsite), m_com_work(nsite) {
    if (*this) {
        m_stats = mem_utils::make_unique<MaeStats>(
                "M7.maes.stats",
                "FCIQMC Multidimensional Averaged Estimators",
                MaeStatsRow(m_bilinears.m_rdms, m_bilinears.m_spec_moms, false));
    }
}

Maes::operator bool() const {
    return m_bilinears || m_ref_excits;
}

bool Maes::all_stores_empty() const {
    return m_bilinears.all_stores_empty() && m_ref_excits.all_stores_empty();
}

bool Maes::is_period_cycle(size_t icycle) {
    if (!m_accum_epoch) return false;
    if (!m_period) return false;
    if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
        m_icycle_period_start = icycle;
        return false;
    }
    return !((icycle - m_icycle_period_start) % m_period);
}

void Maes::end_cycle() {
    m_bilinears.end_cycle();
}

void Maes::make_average_contribs(WalkerTableRow &row, const References &refs, const size_t &icycle) {
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

            auto exsig_from_ref = ref.exsig(row.m_mbf);
            auto is_ref_conn = exsig_from_ref && m_bilinears.m_rdms.takes_contribs_from(exsig_from_ref);
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

void Maes::output(size_t icycle, const Hamiltonian &ham) {
    if (!*this) return;
    if (!is_period_cycle(icycle)) return;
    std::cout << m_bilinears.m_rdms.get_energy(ham.m_frm) << std::endl;
}