//
// Created by rja on 26/07/2021.
//

#include "BosonHamiltonian.h"

BosonHamiltonian::BosonHamiltonian(size_t nmode, size_t nboson_max, std::string fname):
        m_nboson_max(nboson_max), m_nmode(nmode), m_omega(0) {}

defs::ham_t BosonHamiltonian::get_element(const field::BosOnv &onv) const {
    defs::ham_t res = 0;
    for (size_t imode = 0ul; imode < m_nmode; ++imode)
        res += m_omega * static_cast<defs::ham_comp_t>(onv[imode]);
    return res;
}

defs::ham_comp_t BosonHamiltonian::get_energy(const field::BosOnv &onv) const {
    return consts::real(get_element(onv));
}

size_t BosonHamiltonian::nci() const {
    return ci_utils::boson_dim(m_nmode, m_nboson_max);
}
