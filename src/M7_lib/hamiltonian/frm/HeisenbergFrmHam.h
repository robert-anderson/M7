//
// Created by anderson on 2/5/22.
//

#ifndef M7_HEISENBERGFRMHAM_H
#define M7_HEISENBERGFRMHAM_H

#include "M7_lib/basis/Lattice.h"

#include "M7_lib/hamiltonian/frm/FrmHam.h"

/**
 * N-dimensional Quantum Heisenberg model with nearest neighbor-only exchange interactions and no external field. J
 * parameter is redundant (only scales the energy) if the model Hamiltonian is the only term, but it is included to
 * handle situations in which its strength relative to other interactions becomes important
 *
 * H = sum<ij>  Jij  Si . Sj
 *
 * where Sk is the vector spin operator of site k. This can be expanded in terms of spin operators that can be directly
 * related to Fock space fermion operators in z-projected spin quantization:
 *
 * H = (1/2) sum<ij>  Jij  (S+i S-j  +  S-i S+j  +  2Szi Szj)
 *
 * so oppositely-spinned fermions occupying distinct, neighboring orbitals exchange sites under the action of the
 * off-diagonal part
 * and for the diagonal part: nearest neighbors of similar spin produce a positive multiple of the interaction strength,
 * whereas those of opposite spin produce a negative factor.
 *
 * with positive J, similarly spinned neighbors incur an energy penalty, and so the model is antiferromagnetic in this
 * regime.
 */
struct HeisenbergFrmHam : SpinModelFrmHam {
    /**
     * interaction strength
     */
    const defs::ham_t m_j;
    Lattice m_lattice;

    HeisenbergFrmHam(defs::ham_t j, Lattice lattice);

    explicit HeisenbergFrmHam(const fciqmc_config::FermionHamiltonian &opts);

    defs::ham_t get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const override;

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

};


#endif //M7_HEISENBERGFRMHAM_H
