//
// Created by Robert J. Anderson on 12/9/21.
//

#include "HolsteinLadderHam.h"

defs::ham_t HolsteinLadderHam::get_coeff_1110(size_t imode, size_t i, size_t j) const {
    if (imode != m_basis.m_frm.isite(i)) return 0;
    if (imode != m_basis.m_frm.isite(j)) return 0;
    return m_g;
}

defs::ham_t HolsteinLadderHam::get_coeff_1101(size_t imode, size_t i, size_t j) const {
    return get_coeff_1110(imode, j, i);
}

defs::ham_t HolsteinLadderHam::get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return 0;
}

defs::ham_t HolsteinLadderHam::get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return 0;
}
