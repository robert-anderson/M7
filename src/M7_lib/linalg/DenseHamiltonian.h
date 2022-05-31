//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include <M7_lib/defs.h>
#include <M7_lib/foreach/MbfForeach.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>

#include "Dense.h"


using namespace mbf_foreach;
class DenseHamiltonian : public dense::SquareMatrix<defs::ham_t> {
    typedef std::unique_ptr<PairBase> unique_t;
    /**
     * examine the requirements of the given Hamiltonian, and build an iterator over the row and column spaces
     * of its matrix representation in a physically appropriate, complete space of MBFs
     * @param h
     *  hamiltonian object
     * @param force_general
     *  if true, use the full combinatorial basis with no symmetries assumed or imposed (for testing)
     * @return
     *  a type-agnostic iterator which fills the dense matrix representation on a call to its loop method
     */
    std::unique_ptr<PairBase> make_pair_iterator(const Hamiltonian& h, sys::Particles particles, bool force_general);

    /**
     * fall-through case if no specific iterator is found or if general iterator is requested
     */
    std::unique_ptr<PairBase> make_pair_iterator(const Hamiltonian& h, sys::Particles particles);

    /**
     * @param h
     *  hamiltonian object
     * @param particles
     *  particle number sector in which to construct the CI matrix
     * @param force_general
     *  if true, use the full combinatorial basis with no symmetries assumed or imposed (for testing)
     * @return
     *  number of rows (and columns) in the matrix representation of H
     */
    size_t nrow(const Hamiltonian& h, sys::Particles particles, bool force_general);

    template<typename mbf_t>
    void loop_over_pair_iterator(PairBase* foreach, const Hamiltonian& h){
        auto fn = [this, &h](const mbf_t &bra, size_t ibra, const mbf_t &ket, size_t iket) {
            (*this)(ibra, iket) = h.get_element(bra, ket);
        };
        foreach->loop(fn);
    }

public:
    DenseHamiltonian(const Hamiltonian& h, sys::Particles particles, bool force_general=false);
    DenseHamiltonian(const Hamiltonian& h, bool force_general=false);
};

#endif //M7_DENSEHAMILTONIAN_H
