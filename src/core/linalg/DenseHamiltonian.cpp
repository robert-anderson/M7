//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/enumerator/CombinationEnumerator.h>
#include <src/core/enumerator/DeterminantEnumerator.h>
#include "DenseHamiltonian.h"

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &source) :
        Matrix<defs::ham_t>(source.nci()) {
    elements::Determinant bra(source.nsite());
    elements::Determinant ket(source.nsite());

    size_t ibra = ~0ul;
    DeterminantEnumerator bra_enum(source.nsite(), source.nelec());

    while (bra_enum.next(bra, ibra)) {
        size_t iket = ~0ul;
        DeterminantEnumerator ket_enum(source.nsite(), source.nelec());
        while (ket_enum.next(ket, iket)) {
            auto h_elem = source.get_element(bra, ket);
            if (!consts::float_is_zero(h_elem)) {
                (*this)(ibra, iket) = h_elem;
            } else ASSERT(consts::floats_nearly_equal(h_elem, (*this)(ibra, iket)));
        }
    }
}

//DenseHamiltonian::DenseHamiltonian(const Hamiltonian &source, DeterminantList &detlist):
//    Matrix<defs::ham_t>(detlist.high_water_mark(0)) {
//    for (size_t ibra=0ul; ibra<m_nrow; ++ibra) {
//        auto bra = detlist.m_determinant(ibra);
//        for (size_t iket = 0ul; iket < m_nrow; ++iket) {
//            auto ket = detlist.m_determinant(iket);
//            auto h_elem = source.get_element(bra, ket);
//            if (!consts::float_is_zero(h_elem)) (*this)(ibra, iket) = h_elem;
//            else ASSERT(consts::floats_nearly_equal(h_elem, (*this)(ibra, iket)));
//        }
//    }
//}

//DenseHamiltonian::DenseHamiltonian(const Hamiltonian &source, const BosonCouplings &bc):
//        Matrix<defs::ham_t>(source.nci()){
//    Determinant dbra(source.nsite());
//    Determinant dket(source.nsite());
//    Permanent pbra(bc.nmode(), bc.nocc_cutoff());
//    Permanent pket(bc.nmode(), bc.nocc_cutoff());
//    // TODO James: generate all matrix elements
//}
//
