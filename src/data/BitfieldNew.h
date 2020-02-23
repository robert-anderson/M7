//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_BITFIELDNEW_H
#define M7_BITFIELDNEW_H

#include "src/defs.h"
#include <cstring>

class BitfieldNew {
    /*
     * only initialize m_data_internal if we cannot
     * depend on a RawBitfield owned by another object.
     */
    std::vector<defs::data_t> m_data_internal;
    defs::data_t *m_data_external = nullptr;
    defs::data_t *m_data = nullptr;
public:
    const size_t m_nbit, m_ndataword;

    BitfieldNew(const size_t &nbit);

    BitfieldNew(const size_t &nbit, defs::data_t *data_external);

    void set(const size_t &i);

    void set(const defs::inds &i);

    void clr(const size_t &i);

    bool get(const size_t &i) const;

    defs::data_t get_dataword(const size_t &i) const;

    void zero();

    bool is_zero() const;

    bool is_null() const;

    std::string to_string(size_t padding = 0) const;

    void print() const;

    size_t nsetbits() const;

    size_t nsetbits_common(const BitfieldNew &other) const;

    size_t nsetbits_cleared(const BitfieldNew &other) const;

    static size_t ndataword(const size_t &nbit) {
        return nbit == 0 ? 0 : 1 + (nbit - 1) / (8 * sizeof(defs::data_t));
    }

    template<typename T>
    static inline void clr_bit(T &x, size_t i) {
        x &= ~((T) 1ul << i);
    }

    template<typename T>
    static inline void set_bit(T &x, size_t i) {
        x |= ((T) 1ul << i);
    }

    template<typename T>
    static inline bool get_bit(T &x, size_t i) {
        return (x >> i) & T(1ul);
    }

    size_t hash(const size_t modular_divisor = ~0ul) const;

    inline int compare(const BitfieldNew &rhs) const {
        return memcmp((void *) (this->m_data), (void *) (rhs.m_data), m_ndataword * sizeof(defs::data_t));
    }

    BitfieldNew &operator=(const BitfieldNew &rhs) {
        memcpy((void *) (this->m_data), (void *) (rhs.m_data), m_ndataword * sizeof(defs::data_t));
        return *this;
    }

    bool operator==(const BitfieldNew &rhs) const {
        return compare(rhs) == 0;
    }

    bool operator>(const BitfieldNew &rhs) const {
        return compare(rhs) > 0;
    }

    bool operator<(const BitfieldNew &rhs) const {
        return compare(rhs) < 0;
    }

    bool operator>=(const BitfieldNew &rhs) const {
        return !(compare(rhs) < 0);
    }

    bool operator<=(const BitfieldNew &rhs) const {
        return !(compare(rhs) > 0);
    }

    defs::inds setinds() const;

    defs::inds clrinds() const;

};

#endif //M7_BITFIELDNEW_H
