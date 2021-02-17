//
// Created by rja on 17/02/2021.
//

#ifndef M7_FLAGFIELDZ_H
#define M7_FLAGFIELDZ_H

#include "BitsetFieldZ.h"

/*
 * two specializations for pretty access syntax, need dummy template arg to make
 * use of partial template specialization
 */

template<typename dummy_t, size_t nind_element>
struct FlagFieldZ : BitsetFieldBaseZ<uint8_t, 0ul, nind_element> {
    typedef BitsetFieldBaseZ<uint8_t, 0ul, nind_element> base_t;
    using base_t::m_element_format;
    using base_t::nbit_dword;
    using base_t::begin;
    FlagFieldZ(std::array<size_t, nind_element> shape) : base_t(shape){}

    bool get(const std::array<size_t, nind_element>& einds) const {
        auto const& ibit = m_element_format.flatten(einds);
        return bit_utils::get(((uint8_t*)begin())[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void set(const std::array<size_t, nind_element>& einds) const {
        auto const& ibit = m_element_format.flatten(einds);
        bit_utils::set(((uint8_t*)begin())[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void clr(const std::array<size_t, nind_element>& einds) const {
        auto const& ibit = m_element_format.flatten(einds);
        bit_utils::clr(((uint8_t*)begin())[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void put(const std::array<size_t, nind_element>& einds, bool v) const {
        v ? set(einds) : clr(einds);
    }
};

template<typename dummy_t>
struct FlagFieldZ<dummy_t, 0ul> : BitsetFieldBaseZ<uint8_t, 0ul, 0ul> {
    typedef BitsetFieldBaseZ<uint8_t, 0ul, 0ul> base_t;
    using base_t::m_element_format;
    using base_t::nbit_dword;
    using base_t::begin;
    FlagFieldZ() : base_t({}){}

    bool get() const {
        return bit_utils::get(((uint8_t*)begin())[0], 0);
    }

    void set() const {
        bit_utils::set(((uint8_t*)begin())[0], 0);
    }

    void clr() const {
        bit_utils::clr(((uint8_t*)begin())[0], 0);
    }

    void put(bool v) const {
        v ? set() : clr();
    }
};


#endif //M7_FLAGFIELDZ_H
