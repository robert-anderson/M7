//
// Created by Robert J. Anderson on 27/02/2020.
//

#include "ExactLinear.h"

ExactLinear::ExactLinear(
        const Hamiltonian& ham, const conf::Document& opts,
        const Wavefunction& wf, bool only_nonzero_h_spawns) :
        Propagator(opts, ham, wf), m_only_nonzero_h_spawns(only_nonzero_h_spawns),
        m_conn_iters(ham),
        m_mag_log(opts.m_propagator.m_max_bloom, 0, 1, opts.m_propagator.m_static_tau, true,
                  opts.m_propagator.m_tau_min, opts.m_propagator.m_tau_max, 0.0, opts.m_propagator.m_period) {}

void ExactLinear::off_diagonal(Wavefunction& wf, const Walker& walker, const uint_t& ipart) {
    auto& src_mbf = walker.m_mbf;
    const wf_t& weight = walker.m_weight[ipart];
    bool src_initiator = walker.is_initiator(ipart, m_nadd_initiator);
    bool src_deterministic = walker.m_deterministic.get(wf.iroot_part(ipart));
    src_mbf.m_decoded.clear();

    DEBUG_ASSERT_NE(weight, 0.0, "shouldn't be trying to propagate off-diagonal from zero weight");
    conn::Mbf conn(src_mbf);
    auto body = [&]() {
        DEBUG_ASSERT_NE(conn.exsig(), 0ul, "diagonal connection generated");
        auto helement = m_ham.get_element(src_mbf, conn);
        if (m_only_nonzero_h_spawns && ! ham::is_significant(helement)) return;
        auto& dst_mbf = m_dst[src_mbf];
        conn.apply(src_mbf, dst_mbf);
        m_mag_log.log(0, helement, 1.0);
        auto delta = -weight * tau() * helement;
        imp_samp_delta(delta, src_mbf, dst_mbf);
        wf.add_spawn(dst_mbf, delta, src_initiator, src_deterministic, ipart, src_mbf, weight);
    };
    m_conn_iters.loop(conn, src_mbf, body);
}

void ExactLinear::diagonal(Wavefunction& wf, Walker& walker, const uint_t& ipart) {
    const ham_comp_t& hdiag = walker.m_hdiag;
    DEBUG_ASSERT_NEAR_EQ(hdiag, m_ham.get_energy(walker.m_mbf), "incorrect diagonal H element cached");
    wf.scale_weight(walker, ipart, 1 - (hdiag - m_shift[ipart]) * tau());
}

void ExactLinear::update(const uint_t& icycle, const Wavefunction &wf) {
    Propagator::update(icycle, wf);
    m_mag_log.update(icycle, m_tau);
}
