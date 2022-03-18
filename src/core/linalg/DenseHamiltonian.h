//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include <foreach/MbfForeach.h>
#include "defs.h"
#include "Dense.h"
#include "hamiltonian/Hamiltonian.h"


using namespace mbf_foreach;
class DenseHamiltonian : public dense::SquareMatrix<defs::ham_t> {

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
    std::unique_ptr<PairBase> make_pair_iterator(const Hamiltonian& h, bool force_general);

    /**
     * fall-through case if no specific iterator is found or if general iterator is requested
     */
    std::unique_ptr<PairBase> make_pair_iterator(const Hamiltonian& h);

    size_t nrow(const Hamiltonian& h, bool force_general);

public:
    DenseHamiltonian(const Hamiltonian& h, bool force_general=false);
};

#endif //M7_DENSEHAMILTONIAN_H
