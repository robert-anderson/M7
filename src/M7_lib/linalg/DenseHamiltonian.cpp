//
// Created by Robert John Anderson on 2020-01-18.
//

#include "DenseHamiltonian.h"


DenseHamiltonianBase::DenseHamiltonianBase(const Hamiltonian& h, sys::Particles particles, bool force_general):
    m_h(h), m_basis_iters(fci::BasisIters::make(h, particles, force_general)){}


DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, sys::Particles particles, bool force_general) :
        DenseHamiltonianBase(h, particles, force_general),
        dense::SquareMatrix<ham_t>(m_basis_iters.niter_single()){

    if (h.m_basis.m_frm && h.m_basis.m_bos) {
        buffered::FrmBosOnv bra(h.m_basis);
        auto ket = bra;
        loop_over_pair_iterator(bra, ket);
    }
    else if (h.m_basis.m_frm) {
        buffered::FrmOnv bra(h.m_basis.m_frm);
        auto ket = bra;
        loop_over_pair_iterator(bra, ket);
    }
    else if (h.m_basis.m_bos) {
        buffered::BosOnv bra(h.m_basis.m_bos);
        auto ket = bra;
        loop_over_pair_iterator(bra, ket);
    }
}

DenseHamiltonian::DenseHamiltonian(const Hamiltonian &h, bool force_general) :
        DenseHamiltonian(h, h.default_particles(), force_general){}
