//
// Created by Robert John Anderson on 2020-01-18.
//

#include "DenseHamiltonian.h"
#include "../enumerators/CombinationEnumerator.h"
#include "../enumerators/BitfieldEnumerator.h"
#include "../fermion/DeterminantConnection.h"

DenseHamiltonian::DenseHamiltonian(const AbInitioHamiltonian &source):
Matrix<defs::ham_t, true>(source.m_nci)
{
    Determinant bra(source.nspatorb());
    Determinant ket(source.nspatorb());

    size_t ibra{~0ul}; defs::inds bra_setinds(source.nelec());
    CombinationEnumerator bra_enum(source.norb(), source.nelec());
    while(bra_enum.next(bra_setinds, ibra)){
        bra.zero();
        bra.set(bra_setinds);
        {
            size_t iket{~0ul}; defs::inds ket_setinds(source.nelec());
            CombinationEnumerator ket_enum(source.norb(), source.nelec());
            while(ket_enum.next(ket_setinds, iket)){
                ket.zero();
                ket.set(ket_setinds);
                auto h_elem = source.get_element(bra, ket);
                if (!consts::float_is_zero(h_elem)) set(ibra, iket, h_elem);
                else assert(consts::floats_nearly_equal(h_elem, get(ibra, iket)));
            }
        }
    }
}

