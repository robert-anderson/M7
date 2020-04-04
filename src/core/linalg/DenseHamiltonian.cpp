//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/enumerator/CombinationEnumerator.h>
#include "DenseHamiltonian.h"

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &source) :
        Matrix<defs::ham_t>(source.nci()) {
    Determinant bra(source.nsite());
    Determinant ket(source.nsite());

    size_t ibra = ~0ul;
    defs::inds bra_setinds(source.nelec());
    CombinationEnumerator bra_enum(source.nsite() * 2, source.nelec());

    while (bra_enum.next(bra_setinds, ibra)) {
        bra.zero();
        bra.set(bra_setinds);
        {
            size_t iket = ~0ul;
            defs::inds ket_setinds(source.nelec());
            CombinationEnumerator ket_enum(source.nsite() * 2, source.nelec());
            while (ket_enum.next(ket_setinds, iket)) {
                ket.zero();
                ket.set(ket_setinds);
                auto h_elem = source.get_element(bra, ket);
                if (!consts::float_is_zero(h_elem))
                    (*this)(ibra, iket) = h_elem;
                else
                    assert(consts::floats_nearly_equal(h_elem, (*this)(ibra, iket)));
            }
        }
    }
}
