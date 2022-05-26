//
// Created by anderson on 01/04/2022.
//

#ifndef M7_TCFRMHAM_H
#define M7_TCFRMHAM_H

#include "GeneralFrmHam.h"
#include "TcHam.h"


struct TcFrmHam : TcHam, GeneralFrmHam {
    explicit TcFrmHam(const fciqmc_config::FermionHamiltonian &opts)
        : TcHam(), GeneralFrmHam(opts) {}

    TcFrmHam(std::string fcidump_fname, bool spin_major, int charge = 0)
        : TcHam(), GeneralFrmHam(fcidump_fname, spin_major, charge) {};

    defs::ham_t get_coeff_3300(size_t a, size_t b, size_t c, size_t i,
                               size_t j, size_t k) const override;

    defs::ham_t get_element_3300(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

    // defs::ham_t get_coeff_2200(size_t i, size_t j, size_t k, size_t l) const
    // override;

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

};

#endif  // M7_TCFRMHAM_H
