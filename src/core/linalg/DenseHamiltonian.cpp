//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/enumerator/Enumerators.h>
#include "DenseHamiltonian.h"

DenseHamiltonian::DenseHamiltonian(const FermionHamiltonian &source) :
        Matrix<defs::ham_t>(source.nci()) {
    elements::FermionOnv bra(source.nsite());
    elements::FermionOnv ket(source.nsite());

    size_t ibra = ~0ul;
    enums::FermionOnv bra_enum(source.nsite(), source.nelec());

    while (bra_enum.next(bra, ibra)) {
        size_t iket = ~0ul;
        enums::FermionOnv ket_enum(source.nsite(), source.nelec());
        while (ket_enum.next(ket, iket)) {
            auto h_elem = source.get_element(bra, ket);
            if (!consts::float_is_zero(h_elem)) {
                (*this)(ibra, iket) = h_elem;
            } else ASSERT(consts::floats_nearly_equal(h_elem, (*this)(ibra, iket)));
        }
    }
}

DenseHamiltonian::DenseHamiltonian(const FermiBosHamiltonian &source) :
        Matrix<defs::ham_t>(source.nci()) {
    const auto nsite = source.nsite();
    const auto nelec = source.nelec();
    const auto nmode = source.nmode();
    const auto nboson_cutoff = source.nboson_cutoff();
    elements::FermiBosOnv bra(nsite, nmode);
    elements::FermiBosOnv ket(nsite, nmode);

    size_t ibra = ~0ul;
    enums::FermiBosOnv bra_enum(nsite, nelec, nmode, nboson_cutoff);

    while (bra_enum.next(bra, ibra)) {
        size_t iket = ~0ul;
        enums::FermiBosOnv ket_enum(nsite, nelec, nmode, nboson_cutoff);
        while (ket_enum.next(ket, iket)) {
            auto h_elem = source.get_element(bra, ket);
            if (!consts::float_is_zero(h_elem)) {
                (*this)(ibra, iket) = h_elem;
            } else ASSERT(consts::floats_nearly_equal(h_elem, (*this)(ibra, iket)));
        }
    }
}


//DenseHamiltonian::DenseHamiltonian(const FermionHamiltonian &source, DeterminantList &detlist):
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

//DenseHamiltonian::DenseHamiltonian(const FermionHamiltonian &source, const BosonCouplings &bc):
//        Matrix<defs::ham_t>(source.nci()){
//    FermionOnv dbra(source.nsite());
//    FermionOnv dket(source.nsite());
//    Permanent pbra(bc.nmode(), bc.nocc_cutoff());
//    Permanent pket(bc.nmode(), bc.nocc_cutoff());
//    // TODO James: generate all matrix elements
//}
//
