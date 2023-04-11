//
// Created by Oskar Weser on 17/3/22.
//

#ifndef M7_SPINSQUAREFRMHAM_H
#define M7_SPINSQUAREFRMHAM_H

#include "M7_lib/hamiltonian/frm/FrmHam.h"

/**
 * in atomic units:
 *  S^2 = S+S- + Sz * (Sz - 1)
 */
struct SpinSquareFrmHam : FrmHam, ElecSpecTerm {

    // This is Sz * (Sz - 1) which stays constant
    const ham_comp_t m_sz_term;

    /*
     * electron number and Sz sector information is required in order to store the required conserved m_sz_term
     */
    SpinSquareFrmHam(const sys::frm::Sector& sector);

    ham_t get_coeff_2200(uint_t a, uint_t b, uint_t i, uint_t j) const override;

    ham_t get_element_0000(const field::FrmOnv &onv) const override;

    ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;
};

#endif //M7_SPINSQUAREFRMHAM_H
