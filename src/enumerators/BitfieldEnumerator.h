//
// Created by Robert John Anderson on 2020-02-03.
//

#ifndef M7_BITFIELDENUMERATOR_H
#define M7_BITFIELDENUMERATOR_H

#include <x86intrin.h>
#include <assert.h>
#include "Enumerator.h"
#include "../data/BitfieldNew.h"
#include "../fermion/Determinant.h"

enum BitfieldOps {
    null_op, not_op, and_op, and_not_op, xor_op
};

template<enum BitfieldOps op>
class BitfieldEnumerator : public Enumerator<size_t> {
protected:
    const BitfieldNew &m_data1, m_data2;
    size_t m_offset;
    size_t m_idata = (size_t) -1;
    // working variables
    defs::data_t m_work = 0ul;

public:
    BitfieldEnumerator(
            const BitfieldNew &data1,
            const BitfieldNew &data2,
            Enumerator *subsequent = nullptr,
            size_t offset = 0) :
            Enumerator<size_t>(subsequent),
            m_data1(data1), m_data2(data2), m_offset(offset) {
        assert(data1.m_ndataword == data2.m_ndataword);
    }

    inline bool next_element(size_t &result) override {
        if (!m_work) {
            if (m_idata + 1 == m_data1.m_ndataword) return false;
            m_work = get_work(++m_idata);
            return next_element(result);
        }
        if (sizeof(defs::data_t) == 8) {
            result = __tzcnt_u64(m_work) + 64 * m_idata;
        } else if (sizeof(defs::data_t) == 4) {
            result = _tzcnt_u32(m_work) + 32 * m_idata;
        }
        BitfieldNew::clr_bit(m_work, result);
        result += m_offset;
        if (op==not_op){
            /*
             * enumeration over cleared bit positions must in all cases be terminated
             * at the length of the represented bitfield.
             */
            if (result>=m_data1.m_nbit + m_offset) return false;
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

class BitfieldSetEnumerator : public BitfieldEnumerator<null_op> {
public:
    BitfieldSetEnumerator(const BitfieldNew &data1,
                          Enumerator<size_t> *subsequent = nullptr, size_t offset = 0)
            : BitfieldEnumerator<null_op>(data1, data1, subsequent, offset) {}
};

class BitfieldClrEnumerator : public BitfieldEnumerator<not_op> {
public:
    BitfieldClrEnumerator(const BitfieldNew &data1,
                          Enumerator<size_t> *subsequent = nullptr, size_t offset = 0)
            : BitfieldEnumerator<not_op>(data1, data1, subsequent, offset) {}
};

typedef BitfieldEnumerator<and_op> BitfieldAndEnumerator;
typedef BitfieldEnumerator<and_not_op> BitfieldAndNotEnumerator;
typedef BitfieldEnumerator<xor_op> BitfieldXorEnumerator;

template<enum BitfieldOps op>
class DeterminantEnumerator : public BitfieldEnumerator<op> {
    BitfieldEnumerator<op> second_bitfield_enumerator;
public:
    DeterminantEnumerator(
            const Determinant &data1,
            const Determinant &data2,
            Enumerator<size_t> *subsequent = nullptr, size_t offset = 0) :
            BitfieldEnumerator<op>(data1.m_bitfields[0], data2.m_bitfields[0], nullptr, offset),
            second_bitfield_enumerator(
                    BitfieldEnumerator<op>(data1.m_bitfields[1],
                                           data2.m_bitfields[1], subsequent,
                                           offset + data1.m_bitfields[0].m_nbit)) {
        this->set_subsequent(&second_bitfield_enumerator);
    }
};

class DeterminantSetEnumerator : public DeterminantEnumerator<null_op> {
public:
    DeterminantSetEnumerator(const Determinant &data1, Enumerator<size_t> *subsequent = nullptr,
                             size_t offset = 0) : DeterminantEnumerator<null_op>(data1, data1, subsequent, offset) {}
};

class DeterminantClrEnumerator : public DeterminantEnumerator<not_op> {
public:
    DeterminantClrEnumerator(const Determinant &data1, Enumerator<size_t> *subsequent = nullptr,
                             size_t offset = 0) :
            DeterminantEnumerator<not_op>(data1, data1, subsequent, offset) {}
};

typedef DeterminantEnumerator<and_op> DeterminantAndEnumerator;
typedef DeterminantEnumerator<and_not_op> DeterminantAndNotEnumerator;
typedef DeterminantEnumerator<xor_op> DeterminantXorEnumerator;

#endif //M7_BITFIELDENUMERATOR_H
