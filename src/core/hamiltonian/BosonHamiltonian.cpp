//
// Created by rja on 26/07/2021.
//

#include "BosonHamiltonian.h"
#include "src/core/io/BosdumpFileReader.h"

BosonHamiltonian::BosonHamiltonian(const BosdumpHeader &header, size_t nboson_max) :
        m_nmode(header.m_nmode), m_nboson(header.m_nboson), m_nboson_max(m_nboson ? m_nboson : nboson_max),
        m_coeffs_1(m_nboson_max ? m_nmode : 0ul),
        m_coeffs_2(m_nboson_max ? m_nmode : 0ul),
        m_contribs_0011(exsig_utils::ex_0011), m_contribs_0022(exsig_utils::ex_0022) {
    if (!m_nboson_max || !m_nmode) return;

    BosdumpFileReader file_reader(header.m_fname);
    defs::inds inds(4);
    defs::ham_t value;

    log::info("Reading Boson Hamiltonian coefficients from file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::float_is_zero(value)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);
        DEBUG_ASSERT_TRUE(exsig_utils::contribs_to(exsig, ranksig),
                          "excitation does not contribute to this operator rank");
        if (ranksig == exsig_utils::ex_0011) {
            m_contribs_0011.set_nonzero(exsig);
            m_coeffs_1.set(inds[0], inds[1], value);
        } else if (ranksig == exsig_utils::ex_0022) {
            m_contribs_0022.set_nonzero(exsig);
            m_coeffs_2.set(inds[0], inds[1], inds[2], inds[3], value);
        }
    }
    log_data();
}

defs::ham_t BosonHamiltonian::get_element(const field::BosOnv &onv) const {
    if (!m_nboson_max) return 0.0;
    defs::ham_t res = 0;
    for (size_t imode = 0ul; imode < m_nmode; ++imode) {
        if (!onv[imode]) continue;
        defs::ham_comp_t occi = onv[imode];
        res += m_coeffs_1.get(imode, imode) * occi;
        for (size_t jmode = 0ul; jmode < imode; ++jmode) {
            if (!onv[jmode]) continue;
            defs::ham_comp_t occj = onv[jmode];
            // imode and jmode are different
            // i, j -> i, j
            res += 2 * m_coeffs_2.get(imode, imode, jmode, jmode) * occi * occj;
        }
        res += 0.5 * m_coeffs_2.get(imode, imode, imode, imode) * occi * (occi - 1);
    }
    return res;
}

defs::ham_comp_t BosonHamiltonian::get_energy(const field::BosOnv &onv) const {
    return consts::real(get_element(onv));
}

defs::ham_t BosonHamiltonian::get_element(const field::BosOnv &src, const conn::BosOnv &conn) const {
    // this Hamiltonian conserves boson number
    if (conn.m_ann.size() != conn.m_cre.size()) return 0.0;
    // single number-conserving boson operators not implemented;
    if(conn.size() == 2) return 0.0;
    if (!conn.size()) return get_element(src);
    if (conn.size() == 4) {
        auto i = conn.m_cre[0].m_imode;
        auto j = conn.m_cre[0].m_nop == 2 ? i : conn.m_cre[1].m_imode;
        auto k = conn.m_ann[0].m_imode;
        auto l = conn.m_ann[0].m_nop == 2 ? k : conn.m_ann[1].m_imode;
        size_t ni = src[i];
        size_t nj = src[j];
        size_t nk = src[k];
        size_t nl = src[l];

        defs::ham_comp_t occ_fac = 1.0;
        if (i == j) {
            if (k == l) {
                DEBUG_ASSERT_NE(i, k, "ii <- ii case implies diagonal element, not double excitation");
                // ii <- kk
                occ_fac = 0.5*std::sqrt((ni + 2) * (ni + 1) * nk * (nk - 1));
            } else {
                // ii <- kl
                occ_fac = std::sqrt((ni + 2) * (ni + 1) * nk * nl);
            }
        } else {
            if (k == l) {
                // ij <- kk
                occ_fac = std::sqrt((ni + 1) * (nj + 1) * nk * (nk - 1));
            } else {
                // ij <- kl
                occ_fac = 2*std::sqrt((ni + 1) * (nj + 1) * nk * nl);
            }
        }
        return m_coeffs_2.phys_element(i, j, k, l) * occ_fac;
    }
    return 0.0;
}

size_t BosonHamiltonian::nci() const {
    return ci_utils::boson_dim(m_nmode, m_nboson_max);
}

void BosonHamiltonian::log_data() const {
    if (!m_contribs_0011.is_nonzero(0ul))
        log::info("1-boson (0011) term has no diagonal (0000) contributions");
    if (!m_contribs_0011.is_nonzero(exsig_utils::ex_0011))
        log::info("1-boson (0011) term has no single-excitation (0011) contributions");
    if (!m_contribs_0022.is_nonzero(0ul))
        log::info("2-boson (0022) term has no diagonal (0000) contributions");
    if (!m_contribs_0022.is_nonzero(exsig_utils::ex_0011))
        log::info("2-boson (0022) term has no single-excitation (0011) contributions");
    if (!m_contribs_0022.is_nonzero(exsig_utils::ex_0022))
        log::info("2-boson (0022) term has no double-excitation (0022) contributions");
}