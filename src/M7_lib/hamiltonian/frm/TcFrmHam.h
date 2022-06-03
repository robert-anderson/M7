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

    defs::ham_t get_coeff_3300(size_t a, size_t b, size_t c, size_t i, size_t j, size_t k) const override;

    defs::ham_t get_element_3300(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv,
                                 const conn::FrmOnv &conn) const override;

};

#endif  // M7_TCFRMHAM_H
