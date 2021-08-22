//
// Created by rja on 05/11/2020.
//

#ifndef M7_BOSONCOUPLINGS_H
#define M7_BOSONCOUPLINGS_H

#include <src/core/io/EbdumpFileReader.h>
#include "src/core/connection/Connections.h"
#include "src/core/field/Fields.h"
#include "src/core/parallel/SharedArray.h"
#include "HamiltonianData.h"

class BosonCouplings {
    size_t index(const size_t &n, const size_t &p, const size_t &q) const {
        return n * m_nmode2 + p * m_nmode + q;
    }

public:

    const size_t m_nboson_max, m_nmode;
    const size_t m_nmode2;
    SharedArray<defs::ham_t> m_v;

    ham_data::TermContribs m_contribs_1110;
    ham_data::TermContribs m_contribs_1101;

    BosonCouplings(size_t nmode, size_t nboson_max, std::string fname) :
            m_nboson_max(nboson_max), m_nmode(nmode), m_nmode2(nmode * nmode),
            m_v(m_nboson_max ? m_nmode2 * nmode : 0ul),
            m_contribs_1110(conn_utils::encode_exsig(1, 1, 1, 0)),
            m_contribs_1101(conn_utils::encode_exsig(1, 1, 0, 1)) {
        if (!m_nboson_max) return;

        defs::inds inds(3);
        defs::ham_t value;
        EbdumpFileReader file_reader(fname);
        REQUIRE_EQ_ALL(file_reader.m_nspatorb, m_nmode, "expected number of boson modes not found in file");

        log::info("Reading Boson coupling Hamiltonian coefficients from file \"" + file_reader.m_fname + "\"...");
        while (file_reader.next(inds, value)) {
            if (consts::float_is_zero(value)) continue;
            auto ranksig = file_reader.ranksig(inds);
            auto exsig = file_reader.exsig(inds, ranksig);
            m_contribs_1110.is_nonzero(exsig);
            m_contribs_1101.is_nonzero(conn_utils::hermconj(exsig));
            m_v.set(index(inds[0], inds[1], inds[2]), value);
        }
    }

    defs::ham_t v(const size_t &n, const size_t &p, const size_t &q) const {
        return m_v[index(n, p, q)];
    }

    defs::ham_t get_element_1(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
        if (conn.m_bos.size() != 1) return 0.0;
        bool is_ann = conn.m_bos.m_ann.size();

        const auto imode = is_ann ? conn.m_bos.m_ann[0].m_imode : conn.m_bos.m_cre[0].m_imode;
        const auto com = size_t(onv.m_bos[imode]) - (is_ann ? 1 : 0);

        switch (conn.m_frm.size()) {
            case 0: {
                const auto occ_fac = std::sqrt(com + 1);
                // fermion ONVs do not differ, so sum over occupied spin orbitals
                defs::ham_t res = 0.0;
                auto fn = [&](const size_t &ibit) {
                    auto isite = ibit % m_nmode;
                    res += v(imode, isite, isite) * occ_fac;
                };
                onv.m_frm.foreach(fn);
                return res;
            }
            case 2: {
                auto isite = conn.m_frm.m_cre[0] % m_nmode;
                auto jsite = conn.m_frm.m_ann[0] % m_nmode;
                if (isite == jsite) return 0.0; // spin conservation
                return v(imode, isite, jsite);
            }
            default:
                return 0;
        }
    }

    defs::ham_t get_element(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
        if (conn.m_bos.size() != 1ul) return 0.0;
        return get_element_1(onv, conn);
    }

};

#endif //M7_BOSONCOUPLINGS_H