//
// Created by anderson on 11/29/21.
//

#include "Mbf.h"

void mbf::set_aufbau_mbf(field::FrmOnv &onv) {
    const auto nalpha = onv.m_hs.m_nelec_alpha;
    const auto nbeta = onv.m_hs.m_nelec_beta;
    DEBUG_ASSERT_EQ(nalpha + nbeta, onv.m_hs.m_nelec, "inconsistent na, nb, nelec");
    onv.zero();
    for (size_t i = 0ul; i < nalpha; ++i) onv.set({0, i});
    for (size_t i = 0ul; i < nbeta; ++i) onv.set({1, i});
}

void mbf::set_aufbau_mbf(field::BosOnv &onv) {
    if (!onv.nelement()) return;
    onv.zero();
    onv[0] = onv.m_hs.m_nboson;
}

void mbf::set_neel_mbf(field::FrmOnv &onv) {
    REQUIRE_TRUE(onv.m_hs.ms2_conserved(), "Neel state requires conserved 2*Ms");
    REQUIRE_LE(std::abs(onv.m_hs.m_ms2), 1,
               "Neel state requires overall spin of -1, 0, or 1");
    onv.zero();
    size_t ispin = 0;
    for (size_t isite = 0ul; isite < onv.m_hs.m_nelec; ++isite) {
        onv.set({ispin, isite});
        ispin = !ispin;
    }
}

void mbf::set_from_def_array(field::FrmOnv &mbf, const std::vector<defs::inds> &def, size_t idef) {
    if (def.empty()) return;
    REQUIRE_LT(idef, def.size(), "MBF definition index OOB");
    mbf.zero();
    auto& definds = def[idef];
    if (definds.size() == mbf.m_hs.m_sites.m_nspinorb) {
        // assume the determinant is specified as a bit string
        for (size_t i=0ul; i<definds.size(); ++i) {
            auto flag = definds[i];
            REQUIRE_LE(flag, 1ul, "this MBF definition is being treated as a bitstring but value is not 0 or 1");
            if (flag) mbf.set(i);
        }
    }
    else {
        // assume the determinant is specified as a vector of set spin orbital indices
        for (auto& ind: definds){
            REQUIRE_LT(ind, mbf.m_sites.m_nspinorb, "spin orbital index OOB");
            mbf.set(ind);
        }
    }
}

void mbf::set_from_def_array(field::BosOnv &mbf, const std::vector<defs::inds> &def, size_t idef) {
    if (def.empty()) return;
    REQUIRE_LT(idef, def.size(), "MBF definition index OOB");
    mbf.zero();
    auto& definds = def[idef];
    REQUIRE_EQ(definds.size(), mbf.m_nelement,
               "length of input vector does not match number of boson modes");
    size_t i = 0ul;
    for (auto &occ: definds) mbf[i++] = occ;
}

void mbf::set(field::FrmOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef) {
    if (!def.m_frm.get().empty()) set_from_def_array(mbf, def.m_frm, idef);
    else if (def.m_neel) set_neel_mbf(mbf);
    else set_aufbau_mbf(mbf);
    REQUIRE_EQ(mbf.nsetbit(), mbf.m_hs.m_nelec, "too many electrons in MBF");
    if (mbf.m_hs.ms2_conserved())
        REQUIRE_EQ(mbf.ms2(), mbf.m_hs.m_ms2, "MBF has incorrect total Ms");
}

void mbf::set(field::BosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef) {
    if (!def.m_bos.get().empty()) set_from_def_array(mbf, def.m_bos, idef);
    else set_aufbau_mbf(mbf);
}

void mbf::set(field::FrmBosOnv &mbf, const fciqmc_config::MbfDef &def, size_t idef) {
    set(mbf.m_frm, def, idef);
    set(mbf.m_bos, def, idef);
}