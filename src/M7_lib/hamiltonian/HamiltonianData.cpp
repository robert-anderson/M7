//
// Created by Robert J. Anderson on 21/08/2021.
//

#include "HamiltonianData.h"

uint_t ham::TermContribs::ind(OpSig exsig) const {
    if (!exsig.contribs_to(m_ranksig)) return ~0ul;
    uint_t ifrm = std::min(exsig.nfrm_cre(), exsig.nfrm_ann());
    DEBUG_ASSERT_LE(ifrm, m_nexsig_contrib_frm, "invalid number of like-indexed fermion operators");
    uint_t ibos = std::min(exsig.nbos_cre(), exsig.nbos_ann());
    DEBUG_ASSERT_LE(ibos, m_nexsig_contrib_bos, "invalid number of like-indexed boson operators");
    return ifrm * m_nexsig_contrib_bos + ibos;
}

ham::TermContribs::TermContribs(OpSig ranksig) :
        m_ranksig(ranksig), m_basesig(ranksig.base()),
        m_nexsig_contrib_frm(ranksig.ncontrib_frm()), m_nexsig_contrib_bos(ranksig.ncontrib_bos()),
        m_exsig_nonzero(m_nexsig_contrib_frm * m_nexsig_contrib_bos, 0) {}

ham::TermContribs &ham::TermContribs::operator=(const ham::TermContribs &other) {
    REQUIRE_EQ(other.m_ranksig, m_ranksig, "incompatible ranks");
    REQUIRE_EQ(other.m_basesig, m_basesig, "incompatible base signatures");
    m_exsig_nonzero = other.m_exsig_nonzero;
    return *this;
}

ham::TermContribs::TermContribs(const ham::TermContribs &other) : TermContribs(other.m_ranksig){
    *this = other;
}

ham::TermContribs::TermContribs(const ham::TermContribs &contribs_1, const ham::TermContribs &contribs_2):
    TermContribs(contribs_1.m_ranksig){
    REQUIRE_EQ(contribs_1.m_ranksig, contribs_2.m_ranksig, "incompatible ranks");
    auto base_nfrm_cre = m_basesig.nfrm_cre();
    auto base_nfrm_ann = m_basesig.nfrm_ann();
    auto base_nbos_cre = m_basesig.nbos_cre();
    auto base_nbos_ann = m_basesig.nbos_ann();
    /*
     * loop over all exsigs which may or may not contribute to the terms
     */
    for (uint_t nfrm=0ul; nfrm < m_nexsig_contrib_frm; ++nfrm) {
        for (uint_t nbos=0ul; nbos < m_nexsig_contrib_bos; ++nbos) {
            OpSig exsig ({base_nfrm_cre + nfrm, base_nfrm_ann + nfrm}, {base_nbos_cre + nbos, base_nbos_ann + nbos});
            /*
             * if either of the hamiltonians has nonzero contribs from the excitation signature, then the exsig
             * is nonzero in the sum of the two
             */
            if (contribs_1.is_nonzero(exsig) || contribs_2.is_nonzero(exsig)) set_nonzero(exsig);
        }
    }
}

void ham::TermContribs::set_nonzero(OpSig exsig) {
    auto i = ind(exsig);
    REQUIRE_NE(i, ~0ul, "exsig doesn't contribute to this ranksig");
    m_exsig_nonzero[i] = 1;
}

bool ham::TermContribs::is_nonzero(OpSig exsig) const {
    const auto i = ind(exsig);
    REQUIRE_NE(i, ~0ul, "exsig doesn't contribute to this ranksig");
    DEBUG_ASSERT_LT(i, m_exsig_nonzero.size(), "contrib index OOB");
    return m_exsig_nonzero[i];
}

bool ham::TermContribs::any_nonzero() const {
    return std::any_of(m_exsig_nonzero.cbegin(), m_exsig_nonzero.cend(), [](bool b){return b;});
}

void ham::TermContribs::bcast(uint_t iroot) {
    mpi::bcast(m_exsig_nonzero, iroot);
}

ham::KramersAttributes::KramersAttributes(const ham::KramersAttributes &attrs_1,
                                          const ham::KramersAttributes &attrs_2) {
    m_conserving_singles = attrs_1.m_conserving_singles && attrs_2.m_conserving_singles;
    m_conserving_doubles = attrs_1.m_conserving_doubles && attrs_2.m_conserving_doubles;
}

bool ham::KramersAttributes::conserving() const {
    return m_conserving_singles && m_conserving_doubles;
}

void ham::KramersAttributes::bcast() {
    mpi::bcast(m_conserving_singles);
    mpi::bcast(m_conserving_doubles);
}
