//
// Created by Robert John Anderson on 2020-01-26.
//

#ifndef M7_DETERMINANTCONNECTION_H
#define M7_DETERMINANTCONNECTION_H


#include "../defs.h"
#include "Determinant.h"

struct DeterminantConnection {
    // inds (zero-based positions) of the bits set in src but not in dst
    defs::inds m_removed_inds;
    size_t m_nremoved{};
    // inds (zero-based positions) of the bits set in dst but not in src
    defs::inds m_inserted_inds;
    size_t m_ninserted{};
    // inds (zero-based positions) of the bits set in both src and dst
    defs::inds m_common_inds;
    size_t m_ncommon{};
    bool m_parity{};

    DeterminantConnection(const size_t &nbit);

    //DeterminantConnection(const Determinant &det);

    /*

    template <bool get_common=true, bool get_parity=true>
    void update(const Bitfield &src, const Bitfield &dst){
        size_t iparity;
        if constexpr (get_parity) iparity = 0;
        if constexpr (get_parity | get_common) m_ncommon = 0;
        m_nremoved = 0;
        m_ninserted = 0;

        assert(src.m_nbit == dst.m_nbit);
        for (auto i{0ul}; i < src.m_nchar; ++i) {
            auto isrc = src.get_char(i);
            auto idst = dst.get_char(i);
            auto nset = chars.nsetbits[isrc & ~idst];
            auto start = chars.setbitinds[isrc & ~idst];
            for (auto j{0ul}; j < nset; ++j) {
                m_removed_inds[m_nremoved] = start[j] + 8 * i;
                m_nremoved++;
                if constexpr (get_parity) {
                    iparity += chars.parities[isrc][idst];
                    // only flip the parity if the number of moved bits and
                    // the number of bits they move through are both ODD.
                    if ((m_ncommon&1ul) & (nset&1ul)) iparity++;
                }
            }
            nset = chars.nsetbits[idst & ~isrc];
            start = chars.setbitinds[idst & ~isrc];
            for (auto j{0ul}; j < nset; ++j) {
                m_inserted_inds[m_ninserted] = start[j] + 8 * i;
                m_ninserted++;
                if constexpr (get_parity) {
                    iparity += chars.parities[idst][isrc];
                    if ((m_ncommon&1ul) & (nset&1ul)) iparity++;
                }
            }
            if constexpr (get_parity | get_common)
                nset = chars.nsetbits[isrc & idst];
            if constexpr (get_common) {
                start = chars.setbitinds[isrc & idst];
                for (auto j{0ul}; j < nset; ++j) {
                    m_common_inds[m_ncommon+j] = start[j] + 8 * i;
                }
            }
            if constexpr (get_parity | get_common)
                m_ncommon+=nset;
        }
        if constexpr (get_parity) m_parity = iparity & 1ul;
    }
     */
};

#endif //M7_DETERMINANTCONNECTION_H
