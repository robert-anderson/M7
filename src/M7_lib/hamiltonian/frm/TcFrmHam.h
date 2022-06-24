//
// Created by anderson on 01/04/2022.
//

#ifndef M7_TCFRMHAM_H
#define M7_TCFRMHAM_H

#include "GeneralFrmHam.h"
#include "M7_lib/hamiltonian/TcHam.h"


struct TcFrmHam : TcHam, GeneralFrmHam {

    TcFrmHam(const FcidumpInfo& info, bool spin_major): TcHam(), GeneralFrmHam(info, spin_major) {};

    explicit TcFrmHam(opt_pair_t opts): TcHam(), GeneralFrmHam(opts) {}

    ham_t get_coeff_3300(uint_t a, uint_t b, uint_t c, uint_t i, uint_t j, uint_t k) const override;

    ham_t get_element_3300(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

    ham_t get_element_0000(const field::FrmOnv &onv) const override;

    ham_t get_element_1100(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

    ham_t get_element_2200(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

};

#endif  // M7_TCFRMHAM_H
