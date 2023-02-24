//
// Created by rja on 03/12/22.
//

#include "CommutingObservable.h"

commuting_obs::Estimator::Estimator(const Hamiltonian* op, const wf::References* refs) :
        m_op(op), m_refs(refs), m_proj_num(refs->m_proj_energy_nums.m_format){
}

void commuting_obs::Estimator::make_numerator_contribs(const Walker& walker) {
    const auto& format = m_proj_num.m_local.m_format;
    const auto nroot = format.m_shape[0ul];
    const auto nreplica = format.m_shape[1ul];
    for (uint_t iroot=0ul; iroot < nroot; ++iroot) {
        const auto ipart_0 = format.flatten(iroot, 0ul);
        const auto op_elem = m_op->get_element(walker.m_mbf, (*m_refs)[ipart_0].mbf());
        if (!ham::is_significant(op_elem)) continue;
        for (uint_t ireplica=0ul; ireplica < nreplica; ++ireplica){
            const auto ipart = format.flatten(iroot, ireplica);
            m_proj_num.m_local[{iroot, ireplica}] += op_elem*walker.m_weight[ipart];
        }
    }
}

void commuting_obs::Estimator::begin_cycle(uint_t /*icycle*/) {
    m_proj_num.m_local.zero();
}

void commuting_obs::Estimator::end_cycle(uint_t /*icycle*/) {
    m_proj_num.all_sum();
}
