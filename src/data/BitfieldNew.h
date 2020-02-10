//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_BITFIELDNEW_H
#define M7_BITFIELDNEW_H

#include "../defs.h"

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

    std::string to_string() const;

    void print() const;

    size_t nsetbits() const;

    size_t nsetbits_common(const BitfieldNew &other) const;

    size_t nsetbits_cleared(const BitfieldNew &other) const;

    static size_t ndataword(const size_t &nbit) {
        return nbit == 0 ? 0 : 1 + (nbit - 1) / (8 * sizeof(defs::data_t));
    }
    template <typename T>
    static inline void clr_bit(T &x, size_t i){
        x &= ~((T)1ul << i);
    }

    template <typename T>
    static inline void set_bit(T &x, size_t i){
        x |= ((T)1ul << i);
    }

    template <typename T>
    static inline bool get_bit(T &x, size_t i) {
        return (x >> i) & T(1ul);
    }

    size_t hash(const size_t modular_divisor=~0ul) const;
};

#endif //M7_BITFIELDNEW_H
