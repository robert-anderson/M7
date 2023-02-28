//
// Created by rja on 26/01/23.
//

#ifndef M7_HARTREEFOCK_H
#define M7_HARTREEFOCK_H

#include "M7_lib/communication/SharedRows.h"
#include "HfExcitHists.h"

struct HartreeFock : shared_rows::Walker {
    hf_excit_hist::Accumulators m_excit_accums;
    HartreeFock(const shared_rows::Walker::src_t& src, TableBase::Loc loc, const conf::HfExcits& excit_opts):
            shared_rows::Walker("Hartree-Fock ONV", src, loc), m_excit_accums(excit_opts, this){}

    void update() {
        shared_rows::Walker::update();
    }
};


#endif //M7_HARTREEFOCK_H
