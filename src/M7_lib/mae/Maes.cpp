//
// Created by Robert J. Anderson on 15/08/2021.
//

#include "Maes.h"
#include "NotfMaeFiller.h"

Maes::Maes(const conf::Mae &opts, const wf::Vectors& wf) :
        m_accum_epoch("MAE accumulation"),
        m_rdms(opts.m_rdm, wf, m_accum_epoch),
        m_spec_moms(opts.m_spec_mom, wf, m_accum_epoch), m_opts(opts),
        m_on_the_fly(m_opts.m_filling_algorithm.m_value == "on_the_fly"){
    if (*this) {
        m_stats = ptr::smart::make_unique<MaeStats>(
                opts.m_stats_path, "FCIQMC Multidimensional Averaged Estimators",
                MaeStatsRow(m_rdms), 1ul);
    }
}

Maes::operator bool() const {
    return m_rdms || m_spec_moms;
}

bool Maes::all_stores_empty() const {
    return m_rdms.all_stores_empty();
}

bool Maes::is_period_cycle(uint_t icycle) {
    if (!m_on_the_fly) return false;
    if (!m_accum_epoch) return false;
    if (!m_opts.m_stats_period) return false;
    if (m_icycle_period_start == ~0ul || m_icycle_period_start == icycle) {
        m_icycle_period_start = icycle;
        return false;
    }
    return !((icycle - m_icycle_period_start) % m_opts.m_stats_period);
}

void Maes::end_cycle() {
    m_rdms.end_cycle();
}

void Maes::make_otf_average_contribs(Walker &row, const shared_rows::Walker* hf, uint_t icycle) {
    if (!m_on_the_fly) return;
    if (!m_accum_epoch) return;
    // the current cycle should be included in the denominator
    if (!row.occupied_ncycle(icycle)) {
        DEBUG_ASSERT_TRUE(row.m_average_weight.is_zero(), "average value should have been rezeroed");
        return;
    }
    wf_comp_t ncycle_occ = row.occupied_ncycle(icycle);

    for (uint_t ipart = 0ul; ipart < row.m_wf_format.m_nelement; ++ipart) {
        auto ipart_replica = row.ipart_replica(ipart);
        // const auto iroot = ipart / row.nreplica();
        /*
         * the "average" weights actually refer to the unnormalized average. The averages are obtained by dividing
         * each by the number of cycles for which the row is occupied.
         */
        const auto av_weight = row.m_average_weight[ipart] / ncycle_occ;

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
                    m_rdms.make_contribs(hf->mbf(), row.m_mbf, ncycle_occ * av_weight_hf * av_weight_rep);
                    m_rdms.make_contribs(row.m_mbf, hf->mbf(), ncycle_occ * av_weight * av_weight_hf_rep);
                }
            }
        }
    }
    row.m_average_weight = 0;
    row.m_icycle_occ = icycle + 1;
}

void Maes::fill_from_wf_hist(const Table<MbfWeightRow>& hist) {
    REQUIRE_TRUE_ALL(all_stores_empty(), "stores should be empty if no on-the-fly contributions have been made");

    logging::info("Filling MAEs using histogrammed partial CI vector composed of {} MBFs", hist.nrow_in_use());

    if (m_opts.m_filling_algorithm.m_value == "outer_product") {
        const auto displ = mpi::evenly_shared_displ(hist.nrow_in_use());
        const auto count = mpi::evenly_shared_count(hist.nrow_in_use());

        auto bra = hist.m_row;
        auto ket = bra;
        for (bra.restart(displ); bra.in_range(displ + count); ++bra) {
            for (ket.restart(); ket; ++ket) {
                const auto exsig = mbf::exsig(bra.m_mbf, ket.m_mbf);
                if (!m_rdms.takes_contribs_from(exsig)) continue;
                const auto contrib = bra.m_weight[0] * ket.m_weight[0];
                m_rdms.make_contribs(bra.m_mbf, ket.m_mbf, contrib);
            }
        }
    }
    else if (m_opts.m_filling_algorithm.m_value == "bitset_isect_hashmap_ri")
        NotfMaeFiller::fill(hist, NotfMaeFiller::Hashmap, &m_rdms);
    else if (m_opts.m_filling_algorithm.m_value == "bitset_isect_pair_loop_ri")
        NotfMaeFiller::fill(hist, NotfMaeFiller::PairLoop, &m_rdms);

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