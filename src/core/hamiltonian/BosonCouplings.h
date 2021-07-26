//
// Created by rja on 05/11/2020.
//

#ifndef M7_BOSONCOUPLINGS_H
#define M7_BOSONCOUPLINGS_H

#include "src/core/connection/Connections.h"

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

    defs::ham_t get_element_0(const fields::BosOnv &onv) const {
        defs::ham_t res = 0;
        for (size_t imode = 0ul; imode < m_nmode; ++imode)
            res += m_omega * static_cast<defs::ham_comp_t>(onv[imode]);
        return res;
    }

    defs::ham_comp_t get_energy(const fields::BosOnv &onv) const {
        return consts::real(get_element_0(onv));
    }

    defs::ham_t get_element_0(const conn::FrmOnv &frm_conn, const conn::BosOnv &bos_conn) const {
        defs::ham_t res = 0;
        if (frm_conn.size()) return res;
        for (size_t imode = 0ul; imode < m_nmode; ++imode)
            res += m_omega * static_cast<defs::ham_comp_t>(bos_conn.com(imode));
        return res;
    }

    defs::ham_t get_element_0(const conn::FrmBosOnv &conn) const {
        return get_element_0(conn.m_frm, conn.m_bos);
    }

    defs::ham_t get_element_1(const size_t& p, const size_t& imode, const size_t& com) const {
        const auto occ_fac = std::sqrt(com + 1);
        return v(p, p, imode) * occ_fac;
    }

    defs::ham_t get_element_1(const conn::FrmOnv &frm_conn, const FrmOps& frm_com, const conn::BosOnv &bos_conn) const {
        const auto imode = bos_conn.changed_mode(0);
        const auto change = bos_conn.changes(0);

        // bosons don't couple to higher fermion excitations (yet?)
        if (std::abs(change) > 1) return 0.0;

        switch (bos_conn.nexcit()) {
            case 0: {
                defs::ham_t res = 0;
                for (auto& iocc: frm_com.inds()){
                    auto isite = iocc % m_nmode;
                    res += get_element_1(isite, imode, bos_conn.com(imode));
                }
                return res;
            }
            case 1: {
                auto p = frm_conn.m_cre[0] % m_nmode;
                auto q = frm_conn.m_ann[0] % m_nmode;
                ASSERT(p != q) // spin conservation
                ASSERT(v(p, q, imode) == 0.0)
                return v(p, q, imode);
            }
            default:
                return 0;
        }
    }

    defs::ham_t get_element_1(const conn::FrmBosOnv &conn, const FrmOps& frm_com) const {
        return get_element_1(conn.m_frm, frm_com, conn.m_bos);
    }

    defs::ham_t get_element(const conn::FrmBosOnv &conn, const FrmOps& frm_com) const {
        return get_element(conn.m_frm, frm_com, conn.m_bos);
    }

    defs::ham_t get_element(const conn::FrmOnv &frm_conn, const FrmOps& frm_com,
                            const conn::BosOnv &bos_conn) const {
        switch (bos_conn.nchanged_mode()) {
            case 0:
                return get_element_0(frm_conn, bos_conn);
            case 1:
                return get_element_1(frm_conn, frm_com, bos_conn);
            default:
                return 0;
        }
    }

};

#endif //M7_BOSONCOUPLINGS_H