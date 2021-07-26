//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "ExactPropagator.h"
#include "FciqmcCalculation.h"

ExactPropagator::ExactPropagator(const Hamiltonian<> &ham, const fciqmc_config::Document &opts,
                                 const NdFormat<defs::ndim_wf>& wf_fmt, bool only_nonzero_h_spawns) :
                                 Propagator(opts, ham, wf_fmt), m_only_nonzero_h_spawns(only_nonzero_h_spawns) {}

void ExactPropagator::off_diagonal(Wavefunction &wf, const size_t &ipart) {
//    const auto &row = wf.m_store.m_row;
//    auto& src_onv = row.m_onv;
//    const defs::wf_t &weight = row.m_weight[ipart];
//    bool src_initiator = row.m_initiator.get(ipart);
//    bool src_deterministic = row.m_deterministic.get(ipart);
//    OccupiedOrbitals occs(src_onv);
//    ASSERT(occs.size() > 0);
//    VacantOrbitals vacs(src_onv);
//    ASSERT(vacs.size() > 0);
//
//    ASSERT(!consts::float_is_zero(weight));
//
//    auto body = [&](const suite::Conns &conn, const Wavefunction::mbf_t& dst_onv, const defs::ham_t &helement){
//        const auto delta = -weight * tau() * helement;
//        wf.add_spawn(dst_onv, delta, src_initiator, src_deterministic, ipart, src_onv, weight);
//    };
    //m_ham.foreach_connection(src_onv, body, true, m_only_nonzero_h_spawns, false);
}

void ExactPropagator::diagonal(Wavefunction &wf, const size_t &ipart) {
    auto &row = wf.m_store.m_row;
    const defs::ham_comp_t &hdiag = row.m_hdiag;
    ASSERT(hdiag == m_ham.get_energy(row.m_mbf));
    wf.scale_weight(ipart, 1 - (hdiag - m_shift[ipart]) * tau());
}

bool ExactPropagator::is_exact() const {
    return true;
}
