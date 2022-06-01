//
// Created by Robert J. Anderson on 12/9/21.
//

#include "GeneralLadderHam.h"


GeneralLadderHam::GeneralLadderHam(const EbdumpHeader &header, bool spin_major, size_t bos_occ_cutoff) :
        FrmBosHam({{header.m_nsite, {PointGroup(), header.m_orbsym}, header.m_spin_resolved}, {header.m_nmode, bos_occ_cutoff}}),
        m_v(m_basis.size(), header.m_uhf),
        m_v_unc(m_basis.m_bos.m_nmode, 0.0) {
    if (!m_basis) return;
    REQUIRE_EQ(bool(m_basis.m_frm), bool(m_basis.m_bos),
               "if the number of sites is non-zero, so also must be the number of boson modes. "
               "NMODE definition may be missing from EBDUMP file header");

    defs::inds inds(3);
    defs::ham_t value;
    EbdumpFileReader file_reader(header.m_fname, spin_major);

    log::info("Reading boson ladder coupled and uncoupled coefficients from file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::nearly_zero(value)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);
        if (ranksig == exsig_utils::ex_1110 || ranksig == exsig_utils::ex_1101) {
            DEBUG_ASSERT_EQ(ranksig, exsig_utils::ex_1110, "ranksig should be either 0010 or 1110");
            m_contribs_1110.set_nonzero(exsig);
            m_contribs_1101.set_nonzero(exsig_utils::hermconj(exsig));
            m_v.set(inds[0], inds[1], inds[2], value);
        }
        // else, ignore any 0001 or 0010 coefficients here - they belong in the BosHam
    }
    log_data();
}

defs::ham_t GeneralLadderHam::get_coeff_0010(size_t imode) const {
    return m_v_unc[imode];
}

defs::ham_t GeneralLadderHam::get_coeff_0001(size_t imode) const {
    return m_v_unc[imode];
}

defs::ham_t GeneralLadderHam::get_coeff_1110(size_t imode, size_t i, size_t j) const {
    return m_v.get(imode, i, j);
}

defs::ham_t GeneralLadderHam::get_coeff_1101(size_t imode, size_t i, size_t j) const {
    return m_v.get(imode, j, i);
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
    auto fn = [&](size_t i) {
        res += m_v.get(imode, i, i);
    };
    onv.m_frm.foreach_setbit(fn);
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
    auto i = frm_conn.m_cre[0];
    auto j = frm_conn.m_ann[0];
    /*
     * respect hermitian conjugation of the fermion-boson operator product: 1110 (boson creation) is the
     * conventionally non-conjugated term
     */
    auto element = cre ? m_v.get(imode, i, j) : m_v.get(imode, j, i);
    element *= occ_fac;
    return frm_conn.phase(onv.m_frm) ? -element : element;
}

defs::ham_t GeneralLadderHam::get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return get_element_coupled(onv, conn.m_frm, conn.m_bos.m_cre[0].m_imode, true);
}

defs::ham_t GeneralLadderHam::get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return get_element_coupled(onv, conn.m_frm, conn.m_bos.m_ann[0].m_imode, false);
}