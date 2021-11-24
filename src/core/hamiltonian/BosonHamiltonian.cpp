//
// Created by rja on 26/07/2021.
//

#include "BosonHamiltonian.h"
#include "src/core/io/BosdumpFileReader.h"

BosonHamiltonian::BosonHamiltonian(const std::string& fname, size_t nboson_max) :
        m_nboson_max(nboson_max), m_nmode(read_nmode(fname)), m_nboson(read_nboson(fname)),
        m_coeffs_1(m_nboson_max ? m_nmode : 0ul),
        m_coeffs_2(m_nboson_max ? m_nmode : 0ul),
        m_contribs_0011(exsig_utils::ex_0011) {
    if (!m_nboson_max || !m_nmode) return;

    BosdumpFileReader file_reader(fname);
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
            res += 2*m_coeffs_2.get(imode, imode, jmode, jmode) * occi * occj;
        }
        res += 0.5*m_coeffs_2.get(imode, imode, imode, imode) * occi * (occi-1);
    }
    return res;
}

defs::ham_comp_t BosonHamiltonian::get_energy(const field::BosOnv &onv) const {
    return consts::real(get_element(onv));
}

defs::ham_t BosonHamiltonian::get_element(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    if (conn.size()) return 0.0;
    return get_element(onv);
}

size_t BosonHamiltonian::nci() const {
    return ci_utils::boson_dim(m_nmode, m_nboson_max);
}

void BosonHamiltonian::log_data() const {
    if (!m_contribs_0011.is_nonzero(0ul))
        log::info("1-boson (0011) term has no diagonal (0000) contributions");
    if (!m_contribs_0011.is_nonzero(exsig_utils::encode(0, 0, 1, 1)))
        log::info("1-boson (0011) term has no single-excitation (0011) contributions");
}