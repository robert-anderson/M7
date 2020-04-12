//
// Created by rja on 27/02/2020.
//

#include <src/core/enumerator/ContainerCombinationEnumerator.h>
#include "ExactPropagator.h"
#include "FciqmcCalculation.h"

ExactPropagator::ExactPropagator(FciqmcCalculation *fciqmc): Propagator(fciqmc) {}

void ExactPropagator::off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                                   SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) {

    auto anticonn = m_fciqmc->m_scratch->anticonn->get(0);
    OccupiedOrbitals occs(src_det);
    assert(occs.m_nind>0);
    VacantOrbitals vacs(src_det);
    assert(vacs.m_nind>0);

    assert(!consts::float_is_zero(*weight));

    Determinant dst_det(src_det.nsite());
    for (size_t iocc = 0ul; iocc < occs.m_nind; ++iocc) {
        auto occ = occs.m_inds[iocc];

        for (size_t ivac = 0ul; ivac < vacs.m_nind; ++ivac) {
            auto vac = vacs.m_inds[ivac];

            anticonn.zero();
            anticonn.add(occ, vac);
            anticonn.apply(src_det, dst_det);
            auto helement = m_ham->get_element_1(anticonn);
            if (consts::float_is_zero(helement)) continue;

            auto delta = -*weight * m_tau * helement;
            auto irank = m_rank_allocator.get_rank(dst_det);
            if (consts::float_is_zero(delta)) continue;
            spawn_list.add(irank, dst_det, delta, flag_initiator);
        }
    }

    ContainerCombinationEnumerator<defs::det_work> occ_enumerator(occs.m_inds, occs.m_nind, 2);
    defs::inds occ_inds(2);

    while (occ_enumerator.next(occ_inds)) {
        {
            ContainerCombinationEnumerator<defs::det_work> vac_enumerator(vacs.m_inds, vacs.m_nind, 2);
            defs::inds vac_inds(2);
            while (vac_enumerator.next(vac_inds)) {

                anticonn.zero();
                anticonn.add(occ_inds[0], occ_inds[1], vac_inds[0], vac_inds[1]);
                anticonn.apply(src_det, dst_det);
                auto helement = m_ham->get_element_2(anticonn);
                if (consts::float_is_zero(helement)) continue;

                auto delta = -*weight * m_tau * helement;
                auto irank = m_rank_allocator.get_rank(dst_det);
                if (consts::float_is_zero(delta)) continue;
                spawn_list.add(irank, dst_det, delta, flag_initiator);
            }
        }
    }
}
