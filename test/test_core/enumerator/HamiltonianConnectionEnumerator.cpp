//
// Created by rja on 06/07/2020.
//

#include <M7_lib/enumerator/ContainerCombinationEnumerator.h>
#include "gtest/gtest.h"

#if 0
class ConnectionList : public MappedList<DeterminantElement> {
public:
    DeterminantField m_determinant;
    NumericField<defs::ham_t> m_helement;

    ConnectionList(const FermionHamiltonian &h, const FermionOnv &ref, size_t nbucket, const defs::ham_comp_t eps) :
            MappedList("test connection list", m_determinant, nbucket),
            m_determinant(this, 1, ref.nsite()), m_helement(this) {
        OccupiedOrbitals occs(ref);
        VacantOrbitals vacs(ref);
        AntisymFermionOnvConnection connection(ref);

        FermionOnv excited(ref.nsite());
        for (size_t iocc = 0ul; iocc < occs.m_nind; ++iocc) {
            const auto &occ = occs.m_inds[iocc];
            for (size_t ivac = 0ul; ivac < vacs.m_nind; ++ivac) {
                const auto &vac = vacs.m_inds[ivac];
                connection.zero();
                connection.add(occ, vac);
                connection.apply(ref, excited);
                auto helement = h.get_element(connection);
                if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                    size_t irow = expand_push(excited);
                    m_helement(irow) = helement;
                }
            }
        }

        ContainerCombinationEnumerator<defs::det_work> occ_enumerator(occs.m_inds, occs.m_nind, 2);
        defs::inds occ_inds(2);
        while (occ_enumerator.next(occ_inds)) {
            {
                ContainerCombinationEnumerator<defs::det_work> vac_enumerator(vacs.m_inds, vacs.m_nind, 2);
                defs::inds vac_inds(2);
                while (vac_enumerator.next(vac_inds)) {
                    connection.zero();
                    connection.add(occ_inds[0], occ_inds[1], vac_inds[0], vac_inds[1]);
                    connection.apply(ref, excited);
                    auto helement = h.get_element(connection);
                    if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                        size_t irow = expand_push(excited);
                        m_helement(irow) = helement;
                        ASSERT(lookup_irow(m_determinant(irow)) == irow);
                    }
                }
            }
        }
    }
};

TEST(HamiltonianConnectionEnumerator, TestHfDet) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());

    auto det = ham.guess_reference(0);
    HamiltonianConnectionEnumerator enumerator(ham, det);

    auto excited = det;
    MatrixElement<defs::ham_t> matrix_element(det);

    ConnectionList conn_list(ham, det, 1000, 1e-12);
    size_t count = 0;
    while (enumerator.next(matrix_element)) {
        ASSERT_NE(matrix_element.aconn.nexcit(), 1); // Brillouin Theorem
        matrix_element.aconn.apply(det, excited);
        ASSERT_EQ(excited.nsetbit(), det.nsetbit());
        ASSERT_EQ(excited, conn_list.m_determinant(count));
        count++;
    }
    ASSERT_EQ(count, 47);
}


TEST(HamiltonianConnectionEnumerator, ExcitedDet) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());

    FermionOnv det(ham.nsite());
    det.set(defs::inds{0, 1, 2, 6, 7, 9});
    HamiltonianConnectionEnumerator enumerator(ham, det);

    auto excited = det;
    MatrixElement<defs::ham_t> matrix_element(det);

    ConnectionList conn_list(ham, det, 1000, 1e-12);
    size_t count = 0;
    while (enumerator.next(matrix_element)) {
        matrix_element.aconn.apply(det, excited);
        ASSERT_EQ(excited.nsetbit(), det.nsetbit());
        ASSERT_EQ(excited, conn_list.m_determinant(count));
        count++;
    }
    ASSERT_EQ(count, 31);
}
#endif