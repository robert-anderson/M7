//
// Created by Oskar Weser on 17/3/22.
//

#ifndef M7_GENERALFRMHAM_H
#define M7_GENERALFRMHAM_H

#include "FrmHam.h"

struct SpinHam : FrmHam {

    SpinHam(size_t nelec, size_t nsite, int ms2_restrict)
        : FrmHam(nelec, nsite, ms2_restrict){}

    SpinHam(const FrmHam &in_ham)
        : SpinHam(in_ham.m_nelec, in_ham.m_nsite, in_ham.m_ms2_restrict) {}

    explicit SpinHam(const fciqmc_config::FermionHamiltonian &opts);

    defs::ham_t get_coeff_1100(size_t i, size_t j) const override;

    defs::ham_t get_coeff_2200(size_t i, size_t j, size_t k, size_t l) const override;

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;
};


#endif //M7_GENERALFRMHAM_H
