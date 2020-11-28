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
            m_dst_onv.clear();
            m_aconn.apply(src_onv, m_dst_onv);
            auto helement = m_ham.get_element_1(m_aconn);
            if (consts::float_is_zero(helement)) continue;

            auto delta = -weight * tau() * helement;
            if (consts::float_is_zero(delta)) continue;
            m_wf.add_spawn(m_dst_onv, delta, src_initiator, false);
        }
        defs::ham_t delta;
        delta = -weight * tau()*off_diagonal_bosons(src_onv, occ, 1);
        if (!consts::float_is_zero(delta))
            m_wf.add_spawn(m_dst_onv, delta, src_initiator, false);
        delta = -weight * tau()*off_diagonal_bosons(src_onv, occ, -1);
        if (!consts::float_is_zero(delta))
            m_wf.add_spawn(m_dst_onv, delta, src_initiator, false);
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
}

void ExactPropagator::diagonal(Wavefunction &m_wf, const size_t &irow) {
    auto hdiag = m_wf.m_walkers.m_hdiag(irow);
    ASSERT(hdiag==m_ham.get_energy(m_wf.m_walkers.m_onv(irow)));
    m_wf.scale_weight(irow, 1 - (hdiag - m_shift) * tau());
}

defs::ham_t ExactPropagator::off_diagonal_bosons(const views::FermionOnv &src_onv, const size_t &occ, int change) {
    return 0.0;
}

defs::ham_t ExactPropagator::off_diagonal_bosons(const views::FermiBosOnv &src_onv, const size_t &occ, int change) {
    const size_t imode = occ < m_ham.nsite() ? occ : occ-m_ham.nsite();
    if (src_onv.m_bonv(imode)==0 && (change<0)) return 0.0;
    else if (src_onv.m_bonv(imode)==m_ham.nboson_cutoff() && (change>0)) return 0.0;

//    m_dst_onv = src_onv;
//    m_dst_onv.m_bonv(imode)+=change;

    m_aconn.zero();
    m_aconn.m_bonvconn.add(imode, change);
    m_dst_onv.clear();
    m_aconn.apply(src_onv, m_dst_onv);
    ASSERT(src_onv.m_fonv==m_dst_onv.m_fonv);
    auto com = m_dst_onv.m_bonv(imode);
    if (change<0) com+=change;
    auto helem = m_ham.bc().get_element_1(imode, imode, com);

#ifndef DNDEBUG
    auto chk_helem = m_ham.bc().v(imode, imode, imode)*std::sqrt(com+1);
    ASSERT(consts::floats_equal(helem, chk_helem));
#endif
    return helem;
}
