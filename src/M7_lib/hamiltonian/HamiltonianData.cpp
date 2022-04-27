//
// Created by Robert J. Anderson on 21/08/2021.
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

ham_data::TermContribs::TermContribs(const ham_data::TermContribs &other) : TermContribs(other.m_ranksig){}

ham_data::TermContribs &ham_data::TermContribs::operator=(const ham_data::TermContribs &other) {
    REQUIRE_EQ(other.m_ranksig, m_ranksig, "incompatible ranks");
    REQUIRE_EQ(other.m_basesig, m_basesig, "incompatible base signatures");
    m_exsig_nonzero = other.m_exsig_nonzero;
    return *this;
}

ham_data::TermContribs::TermContribs(const ham_data::TermContribs &contribs_1, const ham_data::TermContribs &contribs_2)
        : TermContribs(contribs_1.m_ranksig){
    REQUIRE_EQ(contribs_1.m_ranksig, contribs_2.m_ranksig, "incompatible ranks");
    auto base_nfrm_cre = decode_nfrm_cre(m_basesig);
    auto base_nfrm_ann = decode_nfrm_ann(m_basesig);
    auto base_nbos_cre = decode_nbos_cre(m_basesig);
    auto base_nbos_ann = decode_nbos_ann(m_basesig);
    /*
     * loop over all exsigs which may or may not contribute to the terms
     */
    for (size_t nfrm=0ul; nfrm < m_nexsig_contrib_frm; ++nfrm) {
        for (size_t nbos=0ul; nbos < m_nexsig_contrib_frm; ++nbos) {
            auto exsig = encode(base_nfrm_cre + nfrm, base_nfrm_ann + nfrm,
                                base_nbos_cre + nbos, base_nbos_ann + nbos);
            /*
             * if either of the hamiltonians has nonzero contribs from the excitation signature, then the exsig
             * is nonzero in the sum of the two
             */
            if (contribs_1.is_nonzero(exsig) || contribs_2.is_nonzero(exsig)) set_nonzero(exsig);
        }
    }
}

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
    return m_conserving_singles && m_conserving_doubles;
}

ham_data::KramersAttributes::KramersAttributes(const ham_data::KramersAttributes &attrs_1,
                                               const ham_data::KramersAttributes &attrs_2) {
    m_conserving_singles = attrs_1.m_conserving_singles && attrs_2.m_conserving_singles;
    m_conserving_doubles = attrs_1.m_conserving_doubles && attrs_2.m_conserving_doubles;
}
