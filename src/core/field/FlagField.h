//
// Created by rja on 17/02/2021.
//

#ifndef M7_FLAGFIELD_H
#define M7_FLAGFIELD_H

#include "BitsetField.h"

struct FlagsField : BitsetFieldBase<uint8_t> {
    typedef BitsetFieldBase<uint8_t> base_t;
    FlagsField(Row* row, size_t nflag) : base_t(row, 1, nflag) {}

    bool get(const size_t& iflag) const {
        return base_t::base_get(iflag);
    }

    void set(const size_t& iflag) {
        base_t::base_set(iflag);
    }

    void clr(const size_t& iflag) {
        base_t::base_clr(iflag);
    }

    void put(const size_t& iflag, bool v) {
        base_t::base_put(iflag, v);
    }
};


struct FlagField : BitsetFieldBase<uint8_t> {
    typedef BitsetFieldBase<uint8_t> base_t;
    FlagField(Row* row) : base_t(row, 1, 1) {}

    operator bool() const {
        return get();
    }

    bool get() const {
        return base_t::base_get(0);
    }

    void set() {
        base_t::base_set(0);
    }

    void clr() {
        base_t::base_clr(0);
    }

    void put(bool v) {
        base_t::base_put(0, v);
    }
};


#endif //M7_FLAGFIELD_H
