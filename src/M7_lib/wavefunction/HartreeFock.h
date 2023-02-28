//
// Created by rja on 26/01/23.
//

#ifndef M7_HARTREEFOCK_H
#define M7_HARTREEFOCK_H

#include "M7_lib/communication/SharedRows.h"

struct HartreeFock : shared_rows::Walker {
    HartreeFock(const shared_rows::Walker::src_t& src, TableBase::Loc loc):
        shared_rows::Walker("Hartree-Fock ONV", src, loc) {}
};


#endif //M7_HARTREEFOCK_H
