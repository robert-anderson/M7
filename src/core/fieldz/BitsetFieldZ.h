//
// Created by rja on 09/02/2021.
//

#ifndef M7_BITSETFIELDZ_H
#define M7_BITSETFIELDZ_H

#include "FieldBaseZ.h"
#include "src/core/util/utils.h"
#include "src/core/parallel/MPIAssert.h"
#include "NdFieldBaseZ.h"

template<typename T>
struct BitsetFieldBaseZ : FieldBaseZ {
    static_assert(std::is_integral<T>::value, "Basis for bitset field must be an integral type");
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

public:
    BitsetFieldBaseZ(RowZ* row, size_t nitem, size_t nbit) :
            FieldBaseZ(row, integer_utils::divceil(nbit, nbit_dword())*sizeof(T), nitem, typeid(T)),
            m_nbit(nbit), m_item_dsize(m_item_size/sizeof(T)),
            m_nbit_spare(m_item_dsize * nbit_dword() - m_nbit){
    }

protected:
    bool base_get(const size_t &iitem, const size_t &ibit) const {
        ASSERT(ibit < m_nbit);
        return bit_utils::get(((const T*)begin(iitem))[iitem * ibit / nbit_dword()], ibit % nbit_dword());
    }

    bool base_get(const size_t &ibit) const {
        ASSERT(ibit < m_nbit);
        return bit_utils::get(*(const T*)begin(), ibit % nbit_dword());
    }

    void base_set(const size_t &iitem, const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::set(((T*)begin(iitem))[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void base_set(const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::set(*(T*)begin(), ibit % nbit_dword());
    }

    void base_clr(const size_t &iitem, const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::clr(((T*)begin(iitem))[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void base_clr(const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::clr(*(T*)begin(), ibit % nbit_dword());
    }

    void base_put(const size_t &iitem, const size_t &ibit, bool v) {
        v ? base_set(iitem, ibit) : base_clr(iitem, ibit);
    }

    void base_put(const size_t &ibit, bool v) {
        v ? base_set(ibit) : base_clr(ibit);
    }

    void base_set(const size_t &iitem, const defs::inds &setinds) {
        zero();
        for (auto ibit: setinds) base_set(iitem, ibit);
    }

    void base_set(const defs::inds &setinds) {
        zero();
        for (auto ibit: setinds) base_set(ibit);
    }

public:
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

    size_t nsetbit() const {
        size_t result = 0;
        for (size_t idataword = 0ul; idataword<m_item_dsize; ++idataword){
            result+=bit_utils::nsetbit(get_dataword(idataword));
        }
        return result;
    }

    std::string to_string_element(const size_t& iitem) const override {
        std::string res;
        res.reserve(m_nbit);
        for (size_t i = 0ul; i < m_nbit; ++i)
            res += base_get(iitem, i) ? "1" : "0";
        return res;
    }
};

#endif //M7_BITSETFIELDZ_H