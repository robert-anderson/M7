//
// Created by Robert John Anderson on 2020-02-03.
//

#ifndef M7_BITSETENUMERATOR_H
#define M7_BITSETENUMERATOR_H

#include <x86intrin.h>
#include <src/core/basis/DeterminantField.h>
#include "Enumerator.h"
#include "src/core/table/BitsetField.h"
#include "src/core/util/defs.h"

enum BitfieldOps {
    null_op, not_op, and_op, and_not_op, xor_op
};

template<enum BitfieldOps op>
class BitsetEnumerator : public Enumerator<size_t> {
protected:
    const BitsetElement &m_data1, &m_data2;
    size_t m_offset;
    size_t m_idata = ~0ul;
    // working variables
    defs::data_t m_work = 0ul;

public:
    BitsetEnumerator(
        const BitsetElement &data1,
        const BitsetElement &data2,
        Enumerator *subsequent = nullptr,
        size_t offset = 0) :
        Enumerator<size_t>(subsequent),
        m_data1(data1), m_data2(data2), m_offset(offset) {
        ASSERT(data1.compatible_with(data2));
    }

    inline bool next_element(size_t &result) override {
        if (!m_work) {
            if (m_idata + 1 == m_data1.dsize()) return false;
            m_work = get_work(++m_idata);
            return next_element(result);
        }
        result = bit_utils::next_setbit(m_work);
        result+=m_idata*sizeof(defs::data_t)*CHAR_BIT;
        result+=m_offset;
        if (op == not_op) {
            /*
             * enumeration over cleared bit positions must in all cases be terminated
             * at the length of the represented bitfield.
             */
            if (result >= m_data1.nbit() + m_offset) return false;
        }
        return true;
    }

    inline defs::data_t get_work(const size_t &idata) const {
        switch (op) {
            case null_op:
                return m_data1.get_dataword(idata);
            case not_op:
                return ~m_data1.get_dataword(idata);
            case and_op:
                return m_data1.get_dataword(idata) & m_data2.get_dataword(idata);
            case and_not_op:
                return m_data1.get_dataword(idata) & ~m_data2.get_dataword(idata);
            case xor_op:
                return m_data1.get_dataword(idata) ^ m_data2.get_dataword(idata);
        }
    }
};

class BitsetSetEnumerator : public BitsetEnumerator<null_op> {
public:
    BitsetSetEnumerator(const BitsetElement &data1,
                        Enumerator<size_t> *subsequent = nullptr, size_t offset = 0)
        : BitsetEnumerator<null_op>(data1, data1, subsequent, offset) {}
};

class BitsetClrEnumerator : public BitsetEnumerator<not_op> {
public:
    BitsetClrEnumerator(const BitsetElement &data1,
                        Enumerator<size_t> *subsequent = nullptr, size_t offset = 0)
        : BitsetEnumerator<not_op>(data1, data1, subsequent, offset) {}
};

typedef BitsetEnumerator<and_op> BitsetAndEnumerator;
typedef BitsetEnumerator<and_not_op> BitsetAndNotEnumerator;
typedef BitsetEnumerator<xor_op> BitsetXorEnumerator;

template<enum BitfieldOps op>
class DeterminantEnumerator : public BitsetEnumerator<op> {
public:
    DeterminantEnumerator(
        const DeterminantElement &data1,
        const DeterminantElement &data2,
        Enumerator<size_t> *subsequent = nullptr, size_t offset = 0) :
        BitsetEnumerator<op>(data1, data2, nullptr, offset) {}

    inline bool next_element(size_t &result) override {
        auto tmp = BitsetEnumerator<op>::next_element(result);
        const auto nbit = BitsetEnumerator<op>::m_data1.nbit();
        if (result >= nbit) {
            result += 50;//nbit;
            result -= 64;//CHAR_BIT * (BitsetEnumerator<op>::m_data1.m_spec->element_dsize() / 2) * sizeof(defs::data_t);
        }
        return tmp;
    }
};

class DeterminantSetEnumerator : public DeterminantEnumerator<null_op> {
public:
    DeterminantSetEnumerator(const DeterminantElement &data1, Enumerator<size_t> *subsequent = nullptr,
                             size_t offset = 0) : DeterminantEnumerator<null_op>(data1, data1, subsequent, offset) {}
};

class DeterminantClrEnumerator : public DeterminantEnumerator<not_op> {
public:
    DeterminantClrEnumerator(const DeterminantElement &data1, Enumerator<size_t> *subsequent = nullptr,
                             size_t offset = 0) :
        DeterminantEnumerator<not_op>(data1, data1, subsequent, offset) {}
};

typedef DeterminantEnumerator<and_op> DeterminantAndEnumerator;
typedef DeterminantEnumerator<and_not_op> DeterminantAndNotEnumerator;
typedef DeterminantEnumerator<xor_op> DeterminantXorEnumerator;

#endif //M7_BITSETENUMERATOR_H
