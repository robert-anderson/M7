//
// Created by Robert J. Anderson on 15/08/2021.
//

#include "Maes.h"

Maes::Maes(const conf::Mae &opts, sys::Sector sector, uint_t nroot) :
        m_accum_epoch("MAE accumulation"),
        m_rdms(opts.m_rdm, sector, m_accum_epoch),
        m_spec_moms(opts.m_spec_mom, sector, m_accum_epoch),
        m_hf_excits(opts.m_hf_excits, sector.size(), nroot), m_period(opts.m_stats_period) {
    if (*this) {
        m_stats = ptr::smart::make_unique<MaeStats>(
                opts.m_stats_path, "FCIQMC Multidimensional Averaged Estimators",
                MaeStatsRow(m_rdms), 1ul);
    }
}

Maes::operator bool() const {
    return m_rdms || m_spec_moms || m_hf_excits;
}

bool Maes::all_stores_empty() const {
    return m_rdms.all_stores_empty() && m_hf_excits.all_stores_empty();
}

bool Maes::is_period_cycle(uint_t icycle) {
    if (!m_accum_epoch) return false;
    if (!m_period) return false;
    if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
        m_icycle_period_start = icycle;
        return false;
    }
    return !((icycle - m_icycle_period_start) % m_period);
}

void Maes::end_cycle() {
    m_rdms.end_cycle();
}

void Maes::make_average_contribs(Walker &row, const shared_rows::Walker* hf, uint_t icycle) {
    if (!m_accum_epoch) return;
    // the current cycle should be included in the denominator
    if (!row.occupied_ncycle(icycle)) {
        DEBUG_ASSERT_TRUE(row.m_average_weight.is_zero(), "average value should have been rezeroed");
        return;
    }
    wf_comp_t ncycle_occ = row.occupied_ncycle(icycle);

    for (uint_t ipart = 0ul; ipart < row.m_wf_format.m_nelement; ++ipart) {
        auto ipart_replica = row.ipart_replica(ipart);
        const auto iroot = ipart / row.nreplica();
        /*
         * the "average" weights actually refer to the unnormalized average. The averages are obtained by dividing
         * each by the number of cycles for which the row is occupied.
         */
        const auto av_weight = row.m_average_weight[ipart] / ncycle_occ;

        /*
         * accumulate contributions to reference excitations if required
         */
        if (hf) m_hf_excits.make_contribs(row.m_mbf, hf->mbf(), ncycle_occ * av_weight, iroot);

        if (m_rdms) {
            auto av_weight_rep = row.m_average_weight[ipart_replica] / ncycle_occ;
            /*
             * scale up the product by a factor of the number of instantaneous contributions being accounted for in this
             * single averaged contribution (ncycle_occ)
             */
            m_rdms.make_contribs(row.m_mbf, row.m_mbf, ncycle_occ * av_weight * av_weight_rep);

            if (hf) {
                auto exsig_from_hf = mbf::exsig(hf->mbf(), row.m_mbf);
                if ((exsig_from_hf != opsig::c_zero) && m_rdms.takes_contribs_from(exsig_from_hf)) {
                    const auto av_weight_hf = hf->norm_average_weight(icycle, ipart);
                    const auto av_weight_hf_rep = hf->norm_average_weight(icycle, ipart_replica);
                    m_rdms.make_contribs(hf->mbf(), row.m_mbf,
                                                     ncycle_occ * av_weight_hf * av_weight_rep);
                    m_rdms.make_contribs(row.m_mbf, hf->mbf(),
                                                     ncycle_occ * av_weight * av_weight_hf_rep);
                }
            }
        }
    }
    row.m_average_weight = 0;
    row.m_icycle_occ = icycle + 1;
}

void Maes::output(uint_t icycle, const Hamiltonian &ham, bool final) {
    if (!*this) return;
    if (!is_period_cycle(icycle) && !final) return;
    auto& stats_row = m_stats->m_row;

    ham_comp_t rdm_energy = 0.0;
    if (m_rdms.is_energy_sufficient(ham)) rdm_energy = m_rdms.get_energy(ham);

    if (mpi::i_am_root()) {
        stats_row.m_icycle = icycle;
        if (m_rdms) {
            stats_row.m_total_norm = m_rdms.m_total_norm.m_reduced;
            stats_row.m_rdm_energy = rdm_energy;
        }
        m_stats->commit();
    }
}
