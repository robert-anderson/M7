//
// Created by Robert J. Anderson on 12/9/21.
//

#include "GeneralBosHam.h"


GeneralBosHam::GeneralBosHam(const BosdumpHeader &header, size_t occ_cutoff) :
        BosHam({header.m_nmode, occ_cutoff}),
        m_coeffs_1(m_basis.m_nmode), m_coeffs_2(m_basis.m_nmode) {

    BosdumpFileReader file_reader(header.m_fname);
    defs::inds inds(4);
    defs::ham_t value;

    log::info("Reading Boson Hamiltonian coefficients from file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::nearly_zero(value)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);
        DEBUG_ASSERT_TRUE(utils::exsig::contribs_to(exsig, ranksig),
                          "excitation does not contribute to this operator rank");
        if (ranksig == utils::exsig::ex_0011) {
            m_contribs_0011.set_nonzero(exsig);
            m_coeffs_1.set(inds[0], inds[1], value);
        } else if (ranksig == utils::exsig::ex_0022) {
            m_contribs_0022.set_nonzero(exsig);
            m_coeffs_2.set(inds[0], inds[1], inds[2], inds[3], value);
        }
    }
    log_data();
}

GeneralBosHam::GeneralBosHam(opt_pair_t opts) :
        GeneralBosHam(BosdumpHeader(opts.m_ham.m_bosdump.m_path), opts.m_basis.m_bos_occ_cutoff){}

defs::ham_t GeneralBosHam::get_coeff_0011(size_t i, size_t j) const {
    return m_coeffs_1.get(i, j);
}

defs::ham_t GeneralBosHam::get_coeff_0022(size_t i, size_t j, size_t k, size_t l) const {
    return m_coeffs_2.get(i, j, k, l);
}

defs::ham_t GeneralBosHam::get_element_0000(const field::BosOnv &onv) const {
    defs::ham_t h = 0;
    for (size_t imode = 0ul; imode < m_basis.m_nmode; ++imode) {
        if (!onv[imode]) continue;
        defs::ham_comp_t occi = onv[imode];
        h += m_coeffs_1.get(imode, imode) * occi;
        for (size_t jmode = 0ul; jmode < imode; ++jmode) {
            if (!onv[jmode]) continue;
            defs::ham_comp_t occj = onv[jmode];
            // imode and jmode are different
            // i, j -> i, j
            h += 2 * m_coeffs_2.get(imode, imode, jmode, jmode) * occi * occj;
        }
        h += 0.5 * m_coeffs_2.get(imode, imode, imode, imode) * occi * (occi - 1);
    }
    return h;
}

defs::ham_t GeneralBosHam::get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    DEBUG_ASSERT_NE(conn.m_ann.size(), conn.m_cre.size(), "this Hamiltonian conserves boson number");
    DEBUG_ASSERT_EQ(conn.size(), size_t(2), "incorrectly sized connection passed to get_element_0011");
    // get mode indices
    auto a = conn.m_cre[0].m_imode;
    auto i = conn.m_ann[0].m_imode;
    // get occupation at each index
    auto na = onv[a];
    auto ni = onv[i];

    // I think we cannot have i==a otherwise we are in a diagonal
    DEBUG_ASSERT_NE(i, a, "a <- i case implies a diagonal element, not a single excitation");
    return std::sqrt((na+1)*ni);
}

defs::ham_t GeneralBosHam::get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    DEBUG_ASSERT_NE(conn.m_ann.size(), conn.m_cre.size(), "this Hamiltonian conserves boson number");
    DEBUG_ASSERT_NE(conn.size(), size_t(2), "single number-conserving boson operator passed to get_element_0022");
    DEBUG_ASSERT_NE(conn.size(), size_t(0), "empty connection passed to get_element_0022");
    DEBUG_ASSERT_EQ(conn.size(), size_t(4), "incorrectly sized connection passed to get_element_0022");
    // get mode indices
    auto i = conn.m_cre[0].m_imode;
    auto j = conn.m_cre[0].m_nop == 2 ? i : conn.m_cre[1].m_imode;
    auto k = conn.m_ann[0].m_imode;
    auto l = conn.m_ann[0].m_nop == 2 ? k : conn.m_ann[1].m_imode;
    // get occupation of mode at each index
    size_t ni = onv[i];
    size_t nj = onv[j];
    size_t nk = onv[k];
    size_t nl = onv[l];

    defs::ham_comp_t occ_fac = 1.0;
    if (i == j) {
        if (k == l) {
            DEBUG_ASSERT_NE(i, k, "ii <- ii case implies diagonal element, not double excitation");
            // ii <- kk
            occ_fac = 0.5*std::sqrt((ni + 2) * (ni + 1) * nk * (nk - 1));
        } else {
            // ii <- kl
            occ_fac = std::sqrt((ni + 2) * (ni + 1) * nk * nl);
        }
    } else {
        if (k == l) {
            // ij <- kk
            occ_fac = std::sqrt((ni + 1) * (nj + 1) * nk * (nk - 1));
        } else {
            // ij <- kl
            occ_fac = 2*std::sqrt((ni + 1) * (nj + 1) * nk * nl);
        }
    }
    return m_coeffs_2.phys_element(i, j, k, l) * occ_fac;
}
