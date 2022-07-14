//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include <M7_lib/defs.h>
#include <M7_lib/foreach/MbfForeach.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>

#include "Dense.h"
#include "Fci.h"

struct DenseHamiltonianBase {
    const Hamiltonian& m_h;
    fci::BasisIters m_basis_iters;
    DenseHamiltonianBase(const Hamiltonian& h, sys::Particles particles, bool force_general);
};

using namespace mbf_foreach;
class DenseHamiltonian : DenseHamiltonianBase, public dense::SquareMatrix<ham_t> {

    template<typename mbf_t>
    void loop_over_pair_iterator(mbf_t& work_bra, mbf_t& work_ket){
        auto fn = [this, &work_bra, &work_ket](uint_t ibra, uint_t iket) {
            (*this)(ibra, iket) = m_h.get_element(work_bra, work_ket);
        };
        m_basis_iters.m_pair->loop(work_bra, work_ket, fn);
    }

public:
    DenseHamiltonian(const Hamiltonian& h, sys::Particles particles, bool force_general=false);
    DenseHamiltonian(const Hamiltonian& h, bool force_general=false);
};

#endif //M7_DENSEHAMILTONIAN_H
