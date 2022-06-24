//
// Created by Oskar Weser on 17/3/22.
//

#ifndef M7_SPINSQUAREFRMHAM_H
#define M7_SPINSQUAREFRMHAM_H

#include "M7_lib/hamiltonian/frm/FrmHam.h"

struct SpinSquareFrmHam : FrmHam, ElecSpecTerm {

    // This is Sz * (Sz - 1) which stays constant
    const defs::ham_comp_t m_sz_term;

    /*
     * electron number and Sz sector information is required in order to store the required conserved m_sz_term
     */
    SpinSquareFrmHam(const sys::frm::Basis& basis, const sys::frm::Electrons& elecs);

    defs::ham_t get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const override;

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;
};

#endif //M7_SPINSQUAREFRMHAM_H
