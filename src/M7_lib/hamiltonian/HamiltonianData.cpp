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
    auto i = ind(exsig);
    REQUIRE_NE(i, ~0ul, "exsig doesn't contribute to this ranksig");
    m_exsig_nonzero[i] = true;
}

bool ham_data::TermContribs::is_nonzero(size_t exsig) const {
    auto i = ind(exsig);
    REQUIRE_NE(i, ~0ul, "exsig doesn't contribute to this ranksig");
    return m_exsig_nonzero[i];
}

bool ham_data::KramersAttributes::conserving() const {
    return m_conserving_singles && m_conserving_double;
}
