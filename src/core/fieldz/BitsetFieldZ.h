//
// Created by rja on 09/02/2021.
//

#ifndef M7_BITSETFIELDZ_H
#define M7_BITSETFIELDZ_H

#include "FieldBaseZ.h"
#include "src/core/util/utils.h"
#include "src/core/parallel/MPIAssert.h"
#include "NdFieldBaseZ.h"

template<typename T, size_t nind_item, size_t nind_element>
struct BitsetFieldBaseZ : FullyFormattedFieldBaseZ<T, nind_item, nind_element> {
    static_assert(std::is_integral<T>::value, "Basis for bitset field must be an integral type");
    static constexpr size_t nbyte_dword() {return sizeof(T);}
    static constexpr size_t nbit_dword() {return sizeof(T) * CHAR_BIT;}

    // total number of bits stored in one element
    const size_t m_nbit;
    // total number of data words of type T required to store each item
    const size_t m_item_dsize;
    // number of bits unused in last dataword of an item
    const size_t m_nbit_spare;

    using FieldBaseZ::m_item_size;
    using FieldBaseZ::m_size;
    using FieldBaseZ::begin;
    using FieldBaseZ::end;
    using FieldBaseZ::zero;
    typedef FullyFormattedFieldBaseZ<T, nind_item, nind_element> base_t;
    using base_t::m_item_format;
    using base_t::m_element_format;
    using base_t::nelement_all;

public:
    BitsetFieldBaseZ(std::array<size_t, nind_element> shape) :
            base_t(integer_utils::divceil(NdFormat<nind_element>(shape).nelement(), nbit_dword()) * sizeof(T), shape),
            m_nbit(m_element_format.nelement()), m_item_dsize(m_item_size/sizeof(T)),
            m_nbit_spare(m_item_dsize * nbit_dword() - m_nbit){
    }

    bool flat_get(const size_t &iitem, const size_t &ibit) const {
        ASSERT(ibit < m_nbit);
        return bit_utils::get(begin(iitem)[iitem * ibit / nbit_dword()], ibit % nbit_dword());
    }

    void flat_set(const size_t &iitem, const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::set(begin(iitem)[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void flat_clr(const size_t &iitem, const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::clr(begin(iitem)[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void flat_set(const size_t &iitem, const defs::inds &setinds) {
        zero();
        for (auto ibit: setinds) flat_set(iitem, ibit);
    }

    T get_dataword(const size_t &idataword) const {
        auto tmp = ((T *) begin())[idataword];
        if ((idataword + 1)%m_item_dsize==0) {
            tmp = bit_utils::truncate(tmp, m_nbit_spare);
        }
        return tmp;
    }

    T get_antidataword(const size_t &idataword) const {
        auto tmp = ~(((T *) begin())[idataword]);
        if ((idataword + 1)%m_item_dsize==0) {
            tmp = bit_utils::truncate(tmp, m_nbit_spare);
        }
        return tmp;
    }

    std::string to_string_element(const size_t& iitem) const override {
        std::string res;
        res.reserve(m_nbit);
        for (size_t i = 0ul; i < m_nbit; ++i)
            res += flat_get(iitem, i) ? "1" : "0";
        return res;
    }
};


#endif //M7_BITSETFIELDZ_H