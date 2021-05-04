//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "ExactPropagator.h"
#include "FciqmcCalculation.h"

void ExactPropagator::off_diagonal(Wavefunction &wf, const size_t &ipart) {
    const auto &row = wf.m_store.m_row;
    const auto &src_onv = row.m_onv;
    const defs::wf_t &weight = row.m_weight[ipart];
    bool src_initiator = row.m_initiator.get(ipart);
    OccupiedOrbitals occs(src_onv);
    ASSERT(occs.size() > 0);
    VacantOrbitals vacs(src_onv);
    ASSERT(vacs.size() > 0);

    ASSERT(!consts::float_is_zero(weight));

    auto body = [&](const conn::Antisym<0> &conn, const fields::FermionOnv& dst_onv, const defs::ham_t &helement){
        auto delta = -weight * tau() * helement;
        if (consts::float_is_zero(delta)) return;
        wf.add_spawn(dst_onv, delta, src_initiator, false, ipart, src_onv, weight);
    };

    m_ham.foreach_connection(src_onv, body, true, true);
}

void ExactPropagator::diagonal(Wavefunction &wf, const size_t &ipart) {
    auto &row = wf.m_store.m_row;
    const defs::ham_comp_t &hdiag = row.m_hdiag;
    ASSERT(hdiag == m_ham.get_energy(row.m_onv));
    wf.scale_weight(ipart, 1 - (hdiag - m_shift[ipart]) * tau());
}

bool ExactPropagator::is_exact() const {
    return true;
}

ExactPropagator::ExactPropagator(const Hamiltonian<> &ham, const Options &opts,
                                 const NdFormat<defs::ndim_wf> wf_fmt) : Propagator(opts, ham, wf_fmt) {}
