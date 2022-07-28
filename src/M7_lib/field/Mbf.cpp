//
// Created by Robert J. Anderson on 11/29/21.
//

#include "Mbf.h"

void mbf::set_aufbau_mbf(field::FrmOnv &onv, sys::frm::Electrons elecs) {
    onv.zero();
    for (uint_t i = 0ul; i < elecs.m_nalpha; ++i) onv.set({0, i});
    for (uint_t i = 0ul; i < elecs.m_nbeta; ++i) onv.set({1, i});

}

void mbf::set_aufbau_mbf(field::FrmOnv &onv, sys::Particles particles) {
    set_aufbau_mbf(onv, particles.m_frm);
}

void mbf::set_aufbau_mbf(field::BosOnv &onv, sys::bos::Bosons bosons) {
    if (!onv.nelement()) return;
    onv.zero();
    onv[0] = bosons;
}

void mbf::set_aufbau_mbf(field::BosOnv &onv, sys::Particles particles) {
    set_aufbau_mbf(onv, particles.m_bos);
}

void mbf::set_neel_mbf(field::FrmOnv &onv, sys::frm::Electrons elecs) {
    REQUIRE_EQ(size_t(elecs), onv.m_basis.m_nsite, "Neel state requires one electron per site");
    REQUIRE_TRUE(elecs.m_ms2.conserve(), "Neel state requires conserved 2*Ms");
    REQUIRE_LE(std::abs(elecs.m_ms2), 1,
               "Neel state requires overall spin of -1, 0, or 1");
    onv.zero();
    uint_t ispin = elecs.m_ms2 < 0;
    for (uint_t isite = 0ul; isite < elecs; ++isite) {
        onv.set({ispin, isite});
        ispin = !ispin;
    }
}

void mbf::set_from_def_array(field::FrmOnv &mbf, const v_t<uintv_t> &def, uint_t idef) {
    if (def.empty()) return;
    REQUIRE_LT(idef, def.size(), "MBF definition index OOB");
    mbf.zero();
    auto& definds = def[idef];
    if (definds.size() == mbf.m_basis.m_nspinorb) {
        // assume the determinant is specified as a bit string
        for (uint_t i=0ul; i<definds.size(); ++i) {
            auto flag = definds[i];
            REQUIRE_LE(flag, 1ul, "this MBF definition is being treated as a bitstring but value is not 0 or 1");
            if (flag) mbf.set(i);
        }
    }
    else {
        // assume the determinant is specified as a vector of set spin orbital indices
        for (auto& ind: definds){
            REQUIRE_LT(ind, mbf.m_basis.m_nspinorb, "spin orbital index OOB");
            mbf.set(ind);
        }
    }
}

void mbf::set_from_def_array(field::BosOnv &mbf, const v_t<uintv_t> &def, uint_t idef) {
    if (def.empty()) return;
    REQUIRE_LT(idef, def.size(), "MBF definition index OOB");
    mbf.zero();
    auto& definds = def[idef];
    REQUIRE_EQ(definds.size(), mbf.m_nelement,
               "length of input vector does not match number of boson modes");
    uint_t i = 0ul;
    for (auto &occ: definds) mbf[i++] = occ;
}

void mbf::set(field::FrmOnv &mbf, sys::Particles particles, const conf::MbfDef &def, uint_t idef) {
    auto elecs = particles.m_frm;
    if (!def.m_frm.m_value.empty()) set_from_def_array(mbf, def.m_frm, idef);
    else if (def.m_neel) set_neel_mbf(mbf, elecs);
    else set_aufbau_mbf(mbf, elecs);
    REQUIRE_EQ(mbf.nsetbit(), elecs, "too many electrons in MBF");
    if (elecs.m_ms2.conserve())
        REQUIRE_EQ(mbf.ms2(), elecs.m_ms2, "MBF has incorrect total 2*Ms");
}

void mbf::set(field::BosOnv &mbf, sys::Particles particles, const conf::MbfDef &def, uint_t idef) {
    if (!def.m_bos.m_value.empty()) set_from_def_array(mbf, def.m_bos, idef);
    else set_aufbau_mbf(mbf, particles.m_bos);
}

void mbf::set(field::FrmBosOnv &mbf, sys::Particles particles, const conf::MbfDef &def, uint_t idef) {
    set(mbf.m_frm, particles, def, idef);
    set(mbf.m_bos, particles, def, idef);
}
