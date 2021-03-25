//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "ExactPropagator.h"
#include "FciqmcCalculation.h"

void ExactPropagator::off_diagonal(Wavefunction &wf) {
    const auto& row = wf.m_store.m_row;
    const auto& ipart = wf.m_ipart;
    const auto& src_onv = row.m_onv;
    const defs::wf_t& weight = row.m_weight(ipart);
    bool src_initiator = row.m_initiator.get(ipart);
    OccupiedOrbitals occs(src_onv);
    ASSERT(occs.size() > 0);
    VacantOrbitals vacs(src_onv);
    ASSERT(vacs.size() > 0);

    ASSERT(!consts::float_is_zero(weight));

    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        const auto occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            auto vac = vacs[ivac];

            m_aconn.zero();
            m_aconn.add(occ, vac);
            m_dst_onv.zero();
            m_aconn.apply(src_onv, m_dst_onv);
            auto helement = m_ham.get_element_1(m_aconn);
            if (consts::float_is_zero(helement)) continue;

            auto delta = -weight * tau() * helement;
            if (consts::float_is_zero(delta)) continue;
            wf.add_spawn(m_dst_onv, delta, src_initiator, false, ipart, src_onv, weight);
        }
        defs::ham_t delta;
        delta = -weight * tau()* off_diagonal_bosons(m_ham, m_aconn, src_onv, m_dst_onv, occ, 1);
        if (!consts::float_is_zero(delta))
            wf.add_spawn(m_dst_onv, delta, src_initiator, false, ipart, src_onv, weight);
        delta = -weight * tau()* off_diagonal_bosons(m_ham, m_aconn, src_onv, m_dst_onv, occ, -1);
        if (!consts::float_is_zero(delta))
            wf.add_spawn(m_dst_onv, delta, src_initiator, false, ipart, src_onv, weight);
    }

    if (m_ham.int_2e_rank()==2) {
        ContainerCombinationEnumerator<defs::inds> occ_enumerator(occs.inds(), occs.size(), 2);
        defs::inds occ_inds(2);

        while (occ_enumerator.next(occ_inds)) {
            {
                ContainerCombinationEnumerator<defs::inds> vac_enumerator(vacs.inds(), vacs.size(), 2);
                defs::inds vac_inds(2);
                while (vac_enumerator.next(vac_inds)) {

                    m_aconn.zero();
                    m_aconn.add(occ_inds[0], occ_inds[1], vac_inds[0], vac_inds[1]);
                    m_dst_onv.zero();
                    m_aconn.apply(src_onv, m_dst_onv);
                    auto helement = m_ham.get_element_2(m_aconn);
                    if (consts::float_is_zero(helement)) continue;

                    auto delta = -weight * tau() * helement;
                    if (consts::float_is_zero(delta)) continue;
                    wf.add_spawn(m_dst_onv, delta, src_initiator, false, ipart, src_onv, weight);
                }
            }
        }
    }
}

void ExactPropagator::diagonal(Wavefunction &wf) {
    auto& row = wf.m_store.m_row;
    const defs::ham_comp_t& hdiag = row.m_hdiag;
    ASSERT(hdiag==m_ham.get_energy(row.m_onv));
    wf.scale_weight(1 - (hdiag - m_shift) * tau());
}
