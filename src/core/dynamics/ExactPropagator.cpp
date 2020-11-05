//
// Created by rja on 27/02/2020.
//

#if 0

#include <src/core/enumerator/ContainerCombinationEnumerator.h>
#include "ExactPropagator.h"
#include "FciqmcCalculation.h"

ExactPropagator::ExactPropagator(FciqmcCalculation *fciqmc) : Propagator(fciqmc) {}

void ExactPropagator::off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                                   SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) {

    OccupiedOrbitals occs(src_det);
    ASSERT(occs.m_nind > 0);
    VacantOrbitals vacs(src_det);
    ASSERT(vacs.m_nind > 0);

    ASSERT(!consts::float_is_zero(*weight));

    FermionOnv dst_det(src_det.nsite());
    for (size_t iocc = 0ul; iocc < occs.m_nind; ++iocc) {
        auto occ = occs.m_inds[iocc];

        for (size_t ivac = 0ul; ivac < vacs.m_nind; ++ivac) {
            auto vac = vacs.m_inds[ivac];

            m_aconn.zero();
            m_aconn.add(occ, vac);
            m_aconn.apply(src_det, dst_det);
            auto helement = m_ham->get_element_1(m_aconn);
            if (consts::float_is_zero(helement)) continue;

            auto delta = -*weight * tau() * helement;
            if (consts::float_is_zero(delta)) continue;
            spawn(spawn_list, dst_det, delta, flag_initiator, flag_deterministic);
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
                m_aconn.apply(src_det, dst_det);
                auto helement = m_ham->get_element_2(m_aconn);
                if (consts::float_is_zero(helement)) continue;

                auto delta = -*weight * tau() * helement;
                if (consts::float_is_zero(delta)) continue;
                spawn(spawn_list, dst_det, delta, flag_initiator, flag_deterministic);
            }
        }
    }
}

void ExactPropagator::diagonal(const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight,
                               bool flag_deterministic,
                               defs::ham_comp_t &delta_square_norm, defs::ham_comp_t &delta_nw) {
    auto tmp = *weight;
    weight *= (1.0 - (*hdiag - m_shift) * tau());
    delta_square_norm += std::pow(std::abs(*weight), 2) - std::pow(std::abs(tmp), 2);
    delta_nw += std::abs(*weight) - std::abs(tmp);
}


#endif