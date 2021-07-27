//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/enumerator/Enumerators.h>
#include "DenseHamiltonian.h"

void DenseHamiltonian::setup_frm(const Hamiltonian &source) {

    buffered::FrmOnv bra(source.nsite());
    buffered::FrmOnv ket(source.nsite());

    size_t ibra = ~0ul;
    enums::FermionOnv bra_enum(source.nsite(), source.nelec());
    conn::FrmOnv conn(source.nsite());

    while (bra_enum.next(bra, ibra)) {
        size_t iket = ~0ul;
        enums::FermionOnv ket_enum(source.nsite(), source.nelec());
        while (ket_enum.next(ket, iket)) {
            conn.connect(bra, ket);
            auto h_elem = source.get_element(bra, conn);
            if (!consts::float_is_zero(h_elem)) {
                (*this)(ibra, iket) = h_elem;
            } else ASSERT(consts::floats_nearly_equal(h_elem, (*this)(ibra, iket)));
        }
    }
}

void DenseHamiltonian::setup_frmbos(const Hamiltonian &source) {
    buffered::FrmBosOnv bra(source.nsite());
    buffered::FrmBosOnv ket(source.nsite());

    size_t ibra = ~0ul;
    enums::FermiBosOnv bra_enum(source.nsite(), source.nelec(), source.m_bos.m_nmode, source.m_bos.m_nboson_max);
    conn::FrmBosOnv conn(source.nsite());

    while (bra_enum.next(bra, ibra)) {
        size_t iket = ~0ul;
        enums::FermiBosOnv ket_enum(source.nsite(), source.nelec(), source.m_bos.m_nmode, source.m_bos.m_nboson_max);
        while (ket_enum.next(ket, iket)) {
            conn.connect(bra, ket);
            auto h_elem = source.get_element(bra, conn);
            if (!consts::float_is_zero(h_elem)) {
                (*this)(ibra, iket) = h_elem;
            } else ASSERT(consts::floats_nearly_equal(h_elem, (*this)(ibra, iket)));
        }
    }
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &source) : Matrix<defs::ham_t>(source.nci()) {
    if (source.m_bos.m_nboson_max) setup_frmbos(source);
    else setup_frm(source);
}
