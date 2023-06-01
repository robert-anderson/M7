//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"

std::unique_ptr<guide::Wavefunction> Propagator::make_imp_samp_guide(const conf::GuideWavefunction& opts) const {
    const str_t fmt = "Guiding propagation via importance sampling with {} wavefunction";
    if (opts.m_gutzwiller_like.m_enabled) {
        logging::info(fmt, "Gutzwiller-like");
        return ptr::smart::make_poly_unique<guide::Wavefunction, guide::GutzwillerLike>(
                m_ham, opts.m_gutzwiller_like.m_fac.m_value);
    }
    else if (opts.m_suppress_multi_occ.m_enabled) {
        logging::info(fmt, "multiple occupation suppressing");
        return ptr::smart::make_poly_unique<guide::Wavefunction, guide::SuppressMultiOcc>(
                opts.m_suppress_multi_occ.m_fac.m_value);
    }
    return nullptr;
}


void Propagator::diagonal(wf::Vectors &wf, Walker &walker, uint_t ipart) {
    const ham_comp_t& hdiag = walker.m_hdiag;
    DEBUG_ASSERT_NEAR_EQ(hdiag, m_ham.get_energy(walker.m_mbf), "incorrect diagonal H element cached");
    auto death_rate = (hdiag - m_shifts[walker]->m_values[ipart]) * tau();
    if (death_rate == 0.0) return;
    wf.scale_weight(walker, ipart, 1.0 - death_rate);
    /*
     * enhancement = exp (-sum_k tau_k * (S'_k - S_k))
     * where k extends over cycle indices for which the walker lives
     * S'_k is the shift to which the walker was exposed on cycle k
     * S_k is the shift to which the walker should properly have been exposed on cycle k
     */
    walker.m_log_enhancement_fac -= tau() * (m_shifts[walker]->m_values[ipart] - m_shifts.m_spaces[0]->m_values[ipart]);
}

void Propagator::update(uint_t icycle, const wf::Vectors& wf) {
    m_shifts.update(wf, icycle, tau());
}

void Propagator::imp_samp_delta(wf_t& delta, ham_t src_ovlp, const Mbf& dst_mbf) const {
    if (m_imp_samp_guide) {
        delta/=src_ovlp;
        delta*=m_imp_samp_guide->overlap(dst_mbf);
    }
}

void Propagator::imp_samp_delta(wf_t& delta, const Mbf& src_mbf, const Mbf& dst_mbf) const {
    if (m_imp_samp_guide) imp_samp_delta(delta, m_imp_samp_guide->overlap(src_mbf), dst_mbf);
}

hash::digest_t Propagator::checksum() const {
    v_t<hash::digest_t> all;
    const auto local = checksum_();
    mpi::all_gather(local, all);
    hash::digest_t tot = 0;
    for (auto i: all) tot+=i;
    return tot;
}