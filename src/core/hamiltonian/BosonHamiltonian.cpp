//
// Created by rja on 26/07/2021.
//

#include "BosonHamiltonian.h"
#include "src/core/io/BosdumpFileReader.h"

BosonHamiltonian::BosonHamiltonian(size_t nmode, size_t nboson_max, std::string fname):
        m_nboson_max(nboson_max), m_nmode(nmode), m_coeffs(m_nboson_max ? m_nmode * m_nmode : 0ul),
        m_contribs_0011(conn_utils::encode_exsig(0, 0, 1, 1)){
    if (!m_nboson_max) return;

    defs::inds inds(2);
    defs::ham_t value;
    BosdumpFileReader file_reader(fname);
    REQUIRE_EQ_ALL(file_reader.m_nspatorb, m_nmode, "expected number of boson modes not found in file");

    log::info("Reading Boson Hamiltonian coefficients from file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::float_is_zero(value)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);
        m_contribs_0011.is_nonzero(exsig);
        m_coeffs.set(index(inds[0], inds[1]), value);
    }
}

defs::ham_t BosonHamiltonian::get_element(const field::BosOnv &onv) const {
    defs::ham_t res = 0;
    for (size_t imode = 0ul; imode < m_nmode; ++imode)
        res += m_coeffs[index(imode, imode)] * static_cast<defs::ham_comp_t>(onv[imode]);
    return res;
}

defs::ham_comp_t BosonHamiltonian::get_energy(const field::BosOnv &onv) const {
    return consts::real(get_element(onv));
}

size_t BosonHamiltonian::nci() const {
    return ci_utils::boson_dim(m_nmode, m_nboson_max);
}
