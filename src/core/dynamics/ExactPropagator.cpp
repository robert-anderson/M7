//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "ExactPropagator.h"
#include "FciqmcCalculation.h"

ExactPropagator::ExactPropagator(
        const Hamiltonian &ham, const fciqmc_config::Document &opts,
        const NdFormat<defs::ndim_wf> &wf_fmt, bool only_nonzero_h_spawns) :
        Propagator(opts, ham, wf_fmt), m_only_nonzero_h_spawns(only_nonzero_h_spawns),
        m_connections(foreach_conn::make(ham)),
        m_mag_log(opts.m_propagator.m_max_bloom, 0, 1, opts.m_propagator.m_static_tau, true,
                  opts.m_propagator.m_tau_min, opts.m_propagator.m_tau_max, 0.0, opts.m_propagator.m_period) {}

void ExactPropagator::off_diagonal(Wavefunction &wf, const size_t &ipart) {
    const auto &row = wf.m_store.m_row;
    auto &src_mbf = row.m_mbf;
    const defs::wf_t &weight = row.m_weight[ipart];
    bool src_initiator = row.m_initiator.get(ipart);
    bool src_deterministic = row.m_deterministic.get(ipart);
    OccupiedOrbitals occs(src_mbf);
    ASSERT(occs.size() > 0);
    VacantOrbitals vacs(src_mbf);
    ASSERT(vacs.size() > 0);

    ASSERT(!consts::float_is_zero(weight));

    auto body = [&](const conn::Mbf &conn, const fields::Mbf &dst_onv, defs::ham_t helement) {
        m_mag_log.log(0, helement, 1.0);
        const auto delta = -weight * tau() * helement;
        wf.add_spawn(dst_onv, delta, src_initiator, src_deterministic, ipart, src_mbf, weight);
    };
    m_connections->foreach(src_mbf, body, m_only_nonzero_h_spawns);
}

void ExactPropagator::diagonal(Wavefunction &wf, const size_t &ipart) {
    auto &row = wf.m_store.m_row;
    const defs::ham_comp_t &hdiag = row.m_hdiag;
    ASSERT(hdiag == m_ham.get_energy(row.m_mbf));
    wf.scale_weight(ipart, 1 - (hdiag - m_shift[ipart]) * tau());
}

void ExactPropagator::update(const size_t &icycle, const Wavefunction &wf) {
    Propagator::update(icycle, wf);
    m_mag_log.update(icycle, m_tau);
}
