//
// Created by anderson on 12/9/21.
//

#include "GeneralLadderHam.h"

GeneralLadderHam::GeneralLadderHam(const EbdumpHeader &header, size_t nboson_max) :
        LadderHam({header.m_nsite, header.m_nmode}, nboson_max),
        m_v(m_bd, header.m_uhf), m_v_unc(m_bd.m_nmode, 0.0) {
    if (!m_nboson_max || !(m_bd.m_nsite || m_bd.m_nmode)) return;
    REQUIRE_EQ(m_bd.m_nsite == 0, m_bd.m_nmode == 0,
               "if the number of sites is non-zero, so also must be the number of boson modes. "
               "NMODE definition may be missing from EBDUMP file header");

    defs::inds inds(3);
    defs::ham_t value;
    EbdumpFileReader file_reader(header.m_fname);

    log::info("Reading boson ladder coupled and uncoupled coefficients from file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::nearly_zero(value)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);
        if (ranksig == exsig_utils::ex_0010) {
            m_contribs_0010.set_nonzero(exsig);
            m_contribs_0001.set_nonzero(exsig_utils::hermconj(exsig));
            m_v_unc[inds[0]] = value;
        } else {
            DEBUG_ASSERT_EQ(ranksig, exsig_utils::ex_1110, "ranksig should be either 0010 or 1110");
            m_contribs_1110.set_nonzero(exsig);
            m_contribs_1101.set_nonzero(exsig_utils::hermconj(exsig));
            m_v.set(inds[0], inds[1], inds[2], value);
        }
    }
    log_data();
}

defs::ham_t GeneralLadderHam::get_coeff_0010(size_t imode) const {
    return m_v_unc[imode];
}

defs::ham_t GeneralLadderHam::get_coeff_0001(size_t imode) const {
    return m_v_unc[imode];
}

defs::ham_t GeneralLadderHam::get_coeff_1110(size_t imode, size_t j, size_t i) const {
    return m_v.get(imode, i, j);
}

defs::ham_t GeneralLadderHam::get_coeff_1101(size_t imode, size_t j, size_t i) const {
    return m_v.get(imode, i, j);
}

defs::ham_t GeneralLadderHam::get_element_0010(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    return m_v_unc[conn.m_cre[0].m_imode] * conn.occ_fac(onv);
}

defs::ham_t GeneralLadderHam::get_element_0001(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    return m_v_unc[conn.m_ann[0].m_imode] * conn.occ_fac(onv);
}

defs::ham_t GeneralLadderHam::get_element_pure(const field::FrmBosOnv &onv, size_t imode, bool cre) const {
    const auto occ_fac = std::sqrt(size_t(onv.m_bos[imode]) + cre);
    defs::ham_t res = m_v_unc[imode];
    // fermion ONVs do not differ, so sum over occupied spin orbitals
    auto fn = [&](size_t ibit) {
        auto isite = onv.m_frm.isite(ibit);
        res += m_v.get(imode, isite, isite);
    };
    onv.m_frm.foreach(fn);
    return res * occ_fac;
}

defs::ham_t GeneralLadderHam::get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return get_element_pure(onv, conn.m_bos.m_cre[0].m_imode, true);
}

defs::ham_t GeneralLadderHam::get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return get_element_pure(onv, conn.m_bos.m_ann[0].m_imode, false);
}

defs::ham_t GeneralLadderHam::get_element_coupled(const field::FrmBosOnv &onv,
                                                  const conn::FrmOnv &frm_conn, size_t imode, bool cre) const {
    DEBUG_ASSERT_TRUE(onv.m_frm.get(frm_conn.m_ann[0]), "annihilated op not occupied in ONV")
    DEBUG_ASSERT_FALSE(onv.m_frm.get(frm_conn.m_cre[0]), "created op occupied in ONV")
    const auto occ_fac = std::sqrt(size_t(onv.m_bos[imode]) + cre);
    auto isite = onv.m_frm.isite(frm_conn.m_cre[0]);
    auto jsite = onv.m_frm.isite(frm_conn.m_ann[0]);
    /*
     * respect hermitian conjugation of the fermion-boson operator product: 1110 (boson creation) is the
     * conventionally non-conjugated term
     */
    auto element = cre ? m_v.get(imode, isite, jsite) : m_v.get(imode, jsite, isite);
    element *= occ_fac;
    return frm_conn.phase(onv.m_frm) ? -element : element;
}

defs::ham_t GeneralLadderHam::get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return get_element_coupled(onv, conn.m_frm, conn.m_bos.m_cre[0].m_imode, true);
}

defs::ham_t GeneralLadderHam::get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return get_element_coupled(onv, conn.m_frm, conn.m_bos.m_ann[0].m_imode, false);
}