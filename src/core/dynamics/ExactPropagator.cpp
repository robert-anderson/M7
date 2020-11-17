//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "ExactPropagator.h"
#include "FciqmcCalculation.h"


void ExactPropagator::off_diagonal(Wavefunction &m_wf, const size_t &irow) {
    auto src_onv = m_wf.m_walkers.m_onv(irow);
    const auto weight = m_wf.m_walkers.m_weight(irow, 0, 0);
    bool src_initiator = m_wf.m_walkers.m_flags.m_initiator(irow, 0, 0);
    OccupiedOrbitals occs(src_onv);
    ASSERT(occs.m_nind > 0);
    VacantOrbitals vacs(src_onv);
    ASSERT(vacs.m_nind > 0);

    ASSERT(!consts::float_is_zero(weight));

    for (size_t iocc = 0ul; iocc < occs.m_nind; ++iocc) {
        auto occ = occs.m_inds[iocc];

        for (size_t ivac = 0ul; ivac < vacs.m_nind; ++ivac) {
            auto vac = vacs.m_inds[ivac];

            m_aconn.zero();
            m_aconn.add(occ, vac);
            m_dst_onv.clear();
            m_aconn.apply(src_onv, m_dst_onv);
            auto helement = m_ham.get_element_1(m_aconn);
            if (consts::float_is_zero(helement)) continue;

            auto delta = -weight * tau() * helement;
            if (consts::float_is_zero(delta)) continue;
            m_wf.add_spawn(m_dst_onv, delta, src_initiator, false);
        }
    }

    ContainerCombinationEnumerator<defs::det_work> occ_enumerator(occs.m_inds, occs.m_nind, 2);
    defs::inds occ_inds(2);

    while (occ_enumerator.next(occ_inds)) {
        {
            ContainerCombinationEnumerator<defs::det_work> vac_enumerator(vacs.m_inds, vacs.m_nind, 2);
            defs::inds vac_inds(2);
            while (vac_enumerator.next(vac_inds)) {

                m_aconn.zero();
                m_aconn.add(occ_inds[0], occ_inds[1], vac_inds[0], vac_inds[1]);
                m_dst_onv.clear();
                m_aconn.apply(src_onv, m_dst_onv);
                auto helement = m_ham.get_element_2(m_aconn);
                if (consts::float_is_zero(helement)) continue;

                auto delta = -weight * tau() * helement;
                if (consts::float_is_zero(delta)) continue;
                m_wf.add_spawn(m_dst_onv, delta, src_initiator, false);
            }
        }
    }
}


void ExactPropagator::diagonal(Wavefunction &m_wf, const size_t &irow) {
    auto hdiag = m_wf.m_walkers.m_hdiag(irow);
    m_wf.scale_weight(irow, 1 - (hdiag - m_shift) * tau());
}
