//
// Created by anderson on 2/4/22.
//

#ifndef M7_SPINONVFIELD_H
#define M7_SPINONVFIELD_H

#include <src/core/basis/BasisDims.h>
#include "BitsetField.h"

struct SpinOnvField : BitsetField<size_t, 1> {
    typedef BitsetField<size_t, 1> base_t;
    using base_t::get;
    using base_t::set;
    using base_t::clr;
    using base_t::put;
    using base_t::operator=;
    using base_t::inds_t;

    SpinOnvField(Row* row, size_t nsite, std::string name=""): base_t(row, {nsite}, name){}

    SpinOnvField(Row* row, BasisDims bd, std::string name=""): SpinOnvField(row, bd.m_nsite, name){}

    SpinOnvField(const SpinOnvField& other) : base_t(other){}

    SpinOnvField& operator=(const SpinOnvField& other) {
        base_t::operator=(other);
        return *this;
    }

    void exchange(const size_t& isite, const size_t& jsite){
        bool ibit = get(isite), jbit = get(jsite);
        put(isite, jbit);
        put(jsite, ibit);
    }
};

#endif //M7_SPINONVFIELD_H
