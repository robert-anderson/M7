//
// Created by Robert J. Anderson on 07/03/2022.
//

#include "InteractingBoseGasBosHam.h"

InteractingBoseGasBosHam::InteractingBoseGasBosHam(uint_t ndim, uint_t nwave, ham_comp_t ek_scale) :
        BosHam(Planewaves::size(ndim, nwave)), m_planewaves(ndim, nwave), m_ek_scale(ek_scale){}

InteractingBoseGasBosHam::InteractingBoseGasBosHam(BosHam::init_opts_t opts) :
        InteractingBoseGasBosHam(
                opts.m_ham.m_interacting_bose_gas.m_ndim,
                opts.m_ham.m_interacting_bose_gas.m_nwave,
                opts.m_ham.m_interacting_bose_gas.m_ek_scale){}

ham_t InteractingBoseGasBosHam::get_coeff_0022(uint_t i, uint_t j, uint_t k, uint_t l) const {
    return m_planewaves.conserving(i, j, k, l) ? 1.0 : 0.0;
}

ham_t InteractingBoseGasBosHam::get_element_0000(const field::BosOnv& onv) const {
    // total linear momentum
    ham_t tot = 0.0;
    for (uint_t imode=0ul; imode<m_basis.m_nmode; ++imode){
        if (!onv[imode]) continue;
        tot+=m_planewaves.kinetic_energy(imode)*onv[imode];
    }
    return tot * m_ek_scale;
}

ham_t InteractingBoseGasBosHam::get_element_0022(const field::BosOnv&, const conn::BosOnv& conn) const {
    const auto i = conn.m_cre.get_imode(0);
    const auto j = conn.m_cre.get_imode(1);
    const auto k = conn.m_ann.get_imode(0);
    const auto l = conn.m_ann.get_imode(1);
    return InteractingBoseGasBosHam::get_coeff_0022(i, j, k, l);
}
