//
// Created by rja on 27/02/2020.
//

#include <M7_lib/enumerator/ContainerCombinationEnumerator.h>

#include "ExactPropagator.h"

ExactPropagator::ExactPropagator(
        const Hamiltonian &ham, const fciqmc_config::Document &opts,
        const NdFormat<defs::ndim_wf> &wf_fmt, bool only_nonzero_h_spawns) :
        Propagator(opts, ham, wf_fmt), m_only_nonzero_h_spawns(only_nonzero_h_spawns),
        m_conn_iters(ham),
        m_mag_log(opts.m_propagator.m_max_bloom, 0, 1, opts.m_propagator.m_static_tau, true,
                  opts.m_propagator.m_tau_min, opts.m_propagator.m_tau_max, 0.0, opts.m_propagator.m_period) {}

void ExactPropagator::off_diagonal(Wavefunction &wf, const size_t &ipart) {
    const auto &row = wf.m_store.m_row;
    auto &src_mbf = row.m_mbf;
    const defs::wf_t &weight = row.m_weight[ipart];
    bool src_initiator = row.m_initiator.get(ipart);
    bool src_deterministic = row.m_deterministic.get(wf.iroot_part(ipart));
    src_mbf.m_decoded.clear();

    DEBUG_ASSERT_TRUE(weight,"shouldn't be trying to propagate off-diagonal from zero weight");

    auto body = [&](const conn::Mbf &conn) {
        DEBUG_ASSERT_NE(conn.exsig(), 0ul, "diagonal connection generated");
        auto helement = m_ham.get_element(src_mbf, conn);
        if (m_only_nonzero_h_spawns && consts::nearly_zero(helement, helem_tol)) return;
        auto& dst_mbf = m_dst[src_mbf];
        conn.apply(src_mbf, dst_mbf);
        m_mag_log.log(0, helement, 1.0);
        auto delta = -weight * tau() * helement;
        imp_samp_delta(delta, src_mbf, dst_mbf, row.m_hdiag);
        //std::cout << 123 << std::endl; exit(0);
        wf.add_spawn(dst_mbf, delta, src_initiator, src_deterministic, ipart, src_mbf, weight);
    };
    m_conn_iters.loop(src_mbf, body);
}

void ExactPropagator::diagonal(Wavefunction &wf, const size_t &ipart) {
    auto &row = wf.m_store.m_row;
    const defs::ham_comp_t &hdiag = row.m_hdiag;
    DEBUG_ASSERT_NEARLY_EQ(hdiag, m_ham.get_energy(row.m_mbf), consts::eps(hdiag), "incorrect diagonal H element cached");
    wf.scale_weight(ipart, 1 - (hdiag - m_shift[ipart]) * tau());
}

void ExactPropagator::update(const size_t &icycle, const Wavefunction &wf) {
    Propagator::update(icycle, wf);
    m_mag_log.update(icycle, m_tau);
}
