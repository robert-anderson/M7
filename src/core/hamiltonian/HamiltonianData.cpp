//
// Created by rja on 21/08/2021.
//

#include "HamiltonianData.h"

size_t ham_data::TermContribs::ind(size_t exsig) const {
    if (!contribs_to(exsig, m_ranksig)) return ~0ul;
    size_t ifrm = decode_nfrm_cre(exsig);
    DEBUG_ASSERT_LE(ifrm, m_nexsig_contrib_frm, "invalid number of like-indexed fermion operators");
    size_t ibos = decode_nbos_cre(exsig);
    DEBUG_ASSERT_LE(ibos, m_nexsig_contrib_bos, "invalid number of like-indexed boson operators");
    return ifrm * m_nexsig_contrib_bos + ibos;
}

ham_data::TermContribs::TermContribs(size_t ranksig) :
        m_ranksig(ranksig), m_basesig(base_exsig(ranksig)),
        m_nexsig_contrib_frm(ncontrib_frm(ranksig)), m_nexsig_contrib_bos(ncontrib_bos(ranksig)),
        m_exsig_nonzero(m_nexsig_contrib_frm * m_nexsig_contrib_bos, false) {}

void ham_data::TermContribs::set_nonzero(size_t exsig) {
    m_exsig_nonzero[ind(exsig)] = true;
}

bool ham_data::TermContribs::is_nonzero(size_t exsig) const {
    return m_exsig_nonzero[ind(exsig)];
}

bool ham_data::FrmModelAttributes::is_hubbard_1d() const {
    return m_on_site_only_doubles && m_nn_only_singles;
}

bool ham_data::FrmModelAttributes::is_hubbard_1d_pbc() const {
    return m_on_site_only_doubles && (m_nnp_only_singles && !m_nn_only_singles);
}

size_t ham_data::FrmModelAttributes::iorb_to_isite(const size_t &iorb, const size_t &nsite) {
    return iorb < nsite ? iorb : iorb + nsite;
}

bool ham_data::FrmModelAttributes::nearest_neighbors(const size_t &nsite, const size_t &iorb, const size_t &jorb,
                                                     bool periodic) {
    auto isite = iorb_to_isite(iorb, nsite);
    auto jsite = iorb_to_isite(jorb, nsite);
    if (isite + 1 == jsite || jsite + 1 == isite) return true;
    if (periodic) {
        return (isite == 0 && jsite == nsite - 1) || (isite == nsite - 1 && jsite == 0);
    }
    return false;
}

bool
ham_data::FrmModelAttributes::on_site(const size_t &nsite, const size_t &iorb, const size_t &jorb, const size_t &korb,
                                      const size_t &lorb) const {
    auto isite = iorb_to_isite(iorb, nsite);
    auto jsite = iorb_to_isite(jorb, nsite);
    auto ksite = iorb_to_isite(korb, nsite);
    auto lsite = iorb_to_isite(lorb, nsite);
    return isite == jsite && jsite == ksite && ksite == lsite;
}

void ham_data::FrmModelAttributes::nonzero(const size_t &nsite, const size_t &i, const size_t &j) {
    if (i == j) {
        // one-electron term is not purely off-diagonal
        m_nn_only_singles = false;
        m_nnp_only_singles = false;
    } else {
        if (!nearest_neighbors(nsite, i, j, false)) m_nn_only_singles = false;
        if (!nearest_neighbors(nsite, i, j, true)) m_nnp_only_singles = false;
    }
}

void ham_data::FrmModelAttributes::nonzero(const size_t &nsite, const size_t &i, const size_t &j, const size_t &k,
                                           const size_t &l) {
    if (!(i == j && j == k && k == l)) m_on_site_only_doubles = false;
}

bool ham_data::KramersAttributes::conserving() const {
    return m_conserving_singles && m_conserving_double;
}
