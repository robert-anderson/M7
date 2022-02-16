//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include <src/core/basis/MbfForeach.h>
#include "src/defs.h"
#include "Dense.h"
#include "src/core/hamiltonian/Hamiltonian.h"


using namespace mbf_foreach;
class DenseHamiltonian : public dense::SquareMatrix<defs::ham_t> {
    /**
     * examine the requirements of the given Hamiltonian, and build an iterator over the row and column spaces
     * of its matrix representation in a physically appropriate, complete space of MBFs
     * @param h
     *  hamiltonian object
     * @param general
     *  if true, use the full combinatorial fermion-boson basis with no symmetries assumed or imposed (for testing)
     * @return
     *  a type-agnostic iterator which fills the dense matrix representation on a call to its loop method
     */
    std::unique_ptr<PairBase> make_pair_iterator(const Hamiltonian& h, bool general);

    size_t nrow(const Hamiltonian& h, bool general=false);

public:
    DenseHamiltonian(const Hamiltonian& h, bool general=false);
};

#endif //M7_DENSEHAMILTONIAN_H
