//
// Created by rja on 05/11/2020.
//

#ifndef M7_BOSONCOUPLINGS_H
#define M7_BOSONCOUPLINGS_H

#include "FermionHamiltonian.h"
#include "src/core/basis/Connections.h"

class BosonCouplings {
    const size_t m_nboson_cutoff, m_nmode;
    defs::ham_t m_v, m_omega;

public:
    BosonCouplings(size_t nmode, size_t nboson_cutoff, defs::ham_t v, defs::ham_t omega) :
            m_nboson_cutoff(nboson_cutoff), m_nmode(nmode), m_v(v), m_omega(omega) {}

    const size_t &nboson_cutoff() const {
        return m_nboson_cutoff;
    }

    const size_t &nmode() const {
        return m_nmode;
    }

    defs::ham_t v(const size_t &p, const size_t &q, const size_t &n) const {
        return (p == q && p == n) ? m_v : 0.0;
    }

    size_t nci() const {
        return ci_utils::boson_dim(m_nmode, m_nboson_cutoff);
    }

    defs::ham_t get_element_0(const conn::AsFermionOnv &aconn, const conn::BosonOnv &bonvconn) const {
        defs::ham_t res = 0;
        if (aconn.nexcit()) return res;
        for (size_t imode = 0ul; imode < m_nmode; ++imode)
            res += m_omega * bonvconn.com(imode);
        return res;
    }

    defs::ham_t get_element_0(const conn::AsFermiBosOnv &conn) const {
        return get_element_0(conn.m_aconn, conn.m_bonvconn);
    }

    defs::ham_t get_element_1(const conn::AsFermionOnv &aconn, const conn::BosonOnv &bonvconn) const {
        const auto imode = bonvconn.changed_mode(0);
        const auto change = bonvconn.changes(0);
        const auto occ_fac = std::sqrt(bonvconn.com(imode) + 1);
        if (std::abs(change) > 1) return 0.0;
        // bosons don't couple to higher fermion excitations (yet?)
        switch (aconn.nexcit()) {
            case 0: {
                defs::ham_t res = 0;
                for (size_t iocc = 0ul; iocc < aconn.ncom(); ++iocc) {
                    auto p = aconn.com()[iocc] % m_nmode;
                    res += v(p, p, imode) * occ_fac;
                }
                return res;
            }
            case 1: {
                auto p = aconn.cre()[0] % m_nmode;
                auto q = aconn.ann()[0] % m_nmode;
                ASSERT(p != q) // spin conservation
                ASSERT(v(p, q, imode) == 0)
                return v(p, q, imode);
            }
            default:
                return 0;
        }
    }

    defs::ham_t get_element_1(const conn::AsFermiBosOnv &conn) const {
        return get_element_1(conn.m_aconn, conn.m_bonvconn);
    }

    defs::ham_t get_element(const conn::AsFermionOnv &aconn,
                            const conn::BosonOnv &bonvconn) const {
        switch (bonvconn.nchanged_mode()) {
            case 0:
                return get_element_0(aconn, bonvconn);
            case 1:
                return get_element_1(aconn, bonvconn);
            default:
                return 0;
        }
    }

};


#endif //M7_BOSONCOUPLINGS_H
