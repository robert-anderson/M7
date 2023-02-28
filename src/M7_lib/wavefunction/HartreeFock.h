//
// Created by rja on 26/01/23.
//

#ifndef M7_HARTREEFOCK_H
#define M7_HARTREEFOCK_H

#include "M7_lib/communication/SharedRows.h"
#include "HfC2Accumulation.h"

struct HartreeFock : shared_rows::Walker {
    hf_excit_coeffs::HfExcitCoeffs m_accum;
    HartreeFock(const shared_rows::Walker::src_t& src, TableBase::Loc loc):
        shared_rows::Walker("Hartree-Fock ONV", src, loc), m_accum(this){}

    void update() {
        shared_rows::Walker::update();
        m_accum.update();
    }
};


#endif //M7_HARTREEFOCK_H
