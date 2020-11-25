//
// Created by rja on 27/02/2020.
//

#ifndef M7_EXACTPROPAGATOR_H
#define M7_EXACTPROPAGATOR_H

#include "Propagator.h"

class ExactPropagator : public Propagator {

    defs::ham_t off_diagonal_bosons(const views::FermionOnv &src_onv, const size_t &occ, int change);

    defs::ham_t off_diagonal_bosons(const views::FermiBosOnv &src_onv, const size_t &occ, int change);

public:
    ExactPropagator(const Hamiltonian &ham, const Options &opts) : Propagator(ham, opts) {}

    void diagonal(Wavefunction &m_wf, const size_t &irow) override;

    void off_diagonal(Wavefunction &m_wf, const size_t &irow) override;

};

#endif //M7_EXACTPROPAGATOR_H
