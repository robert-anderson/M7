//
// Created by rja on 04/08/2021.
//

#include "UniformFrmBos.h"

bool UniformFrmBos::draw(const FrmBosOnv &src_onv, const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                         defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
    if(!m_h.m_nboson_max) return false;

    DEBUG_ASSERT_EQ(src_onv.m_bos.nelement(), src_onv.m_frm.m_nsite,
                    "excit gen assumes one boson mode per fermion site");

    auto imode = occs[m_prng.draw_uint(occs.size())] % src_onv.m_frm.m_nsite;
    // if cre, attempt to generate an ONV with an additional boson occupying the mode at imode
    // else, attempt to generate a "de-excited" ONV with one less boson occupying the mode at imode
    bool cre;
    size_t curr_occ = src_onv.m_bos[imode];
    DEBUG_ASSERT_LE(curr_occ, m_h.m_nboson_max, "current occupation of selected mode exceeds cutoff");

    // doubly-occupied sites are twice as likely to be drawn
    prob = src_onv.m_frm.site_nocc(imode);
    prob/=occs.size();

    if (curr_occ == m_h.m_nboson_max){
        cre = false; // always de-excite
    }
    else if (curr_occ == 0){
        cre = true; // always excite
    }
    else{
        // either excite or de-excite with equal probability
        cre = m_prng.draw_uint(2)&1ul; // TODO use one PRNG for both by drawing from [0, 2*occ.m_nind)
        prob *= 0.5;
    }
    DEBUG_ASSERT_LE(size_t(curr_occ+cre), m_h.m_nboson_max, "generated boson occupation exceeds cutoff");

    conn.clear();
    if (cre) conn.m_bos.m_cre.add({imode, 1ul});
    else conn.m_bos.m_ann.add({imode, 1ul});

    helem = m_h.get_element(src_onv, conn);
    return true;
}

std::string UniformFrmBos::description() const {
    return log::format("Boson create {} OR annihilate {} if fermion site occupied", 1ul, 1ul);
}

size_t UniformFrmBos::approx_nconn() const {
    // assume there's one excitation and one de-excitation available per electron
    return m_h.m_frm.nelec()*2;
}
