//
// Created by rja on 27/02/2020.
//

#include <src/core/enumerator/VectorCombinationEnumerator.h>
#include "ExactPropagator.h"

ExactPropagator::ExactPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                                 const RankAllocator<DeterminantElement> &rank_allocator)
    : Propagator(input, ham, rank_allocator) {}


void ExactPropagator::off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                                   SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) {

    AntisymConnection connection(src_det);
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

            connection.zero();
            connection.add(occ, vac);
            connection.apply(src_det, dst_det);
            auto helement = m_ham->get_element_1(connection);
            if (consts::float_is_zero(helement)) continue;

            auto delta = -*weight * m_tau * helement;
            auto irank = m_rank_allocator.get_rank(dst_det);
            if (consts::float_is_zero(delta)) continue;
            spawn_list.add(irank, dst_det, delta, flag_initiator);
        }
    }

    VectorCombinationEnumerator occ_enumerator(occs.m_inds, occs.m_nind, 2);
    defs::inds occ_inds(2);

    while (occ_enumerator.next(occ_inds)) {
        {
            VectorCombinationEnumerator vac_enumerator(vacs.m_inds, vacs.m_nind, 2);
            defs::inds vac_inds(2);
            while (vac_enumerator.next(vac_inds)) {

                connection.zero();
                connection.add(occ_inds[0], occ_inds[1], vac_inds[0], vac_inds[1]);
                connection.apply(src_det, dst_det);
                auto helement = m_ham->get_element_2(connection);
                if (consts::float_is_zero(helement)) continue;

                auto delta = -*weight * m_tau * helement;
                auto irank = m_rank_allocator.get_rank(dst_det);
                if (consts::float_is_zero(delta)) continue;
                spawn_list.add(irank, dst_det, delta, flag_initiator);
            }
        }
    }
}
