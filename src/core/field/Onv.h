//
// Created by rja on 21/05/2021.
//

#ifndef M7_ONV_H
#define M7_ONV_H

#include "Fields.h"

namespace onv {

    static size_t nsite(const fields::Onv<0>& onv){
        return onv.m_nsite;
    }
    static size_t nsite(const fields::Onv<1>& onv){
        return onv.m_frm.m_nsite;
    }

};


#endif //M7_ONV_H
