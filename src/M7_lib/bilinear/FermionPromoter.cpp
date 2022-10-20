//
// Created by Robert J. Anderson on 10/08/2021.
//

#include <M7_lib/foreach/BasicForeach.h>
#include <M7_lib/util/Exsig.h>
#include "FermionPromoter.h"


FermionPromoter::FermionPromoter(uint_t ncom, uint_t exsig, uint_t nop_insert) :
    m_ncom(ncom), m_nexcit(exsig::decode_nfrm_cre(exsig)), m_nop_insert(nop_insert),
    m_ncomb(integer::combinatorial(m_ncom, m_nop_insert)),
    m_all_combs(m_nop_insert * m_ncomb){
    REQUIRE_TRUE(exsig::conserves_nfrm(exsig), "Fermion promotion requires fermion number conservation");
    if (!nop_insert) return;

    basic_foreach::rtnd::Ordered<> foreach_comb(m_ncom, m_nop_insert);
    uint_t icomb = 0ul;
    auto fn = [&](const uintv_t& inds) {
        for (uint_t i = 0ul; i < m_nop_insert; ++i) {
            auto j = icomb * m_nop_insert + i;
            ASSERT(j < m_all_combs.size());
            m_all_combs[j] = inds[i];
        }
        ++icomb;
    };
    foreach_comb.loop(fn);
}

const uint_t *FermionPromoter::begin(uint_t icomb) const {
    DEBUG_ASSERT_LT(icomb, m_ncomb, "combination index OOB");
    return m_all_combs.data() + icomb * m_nop_insert;
}

uint_t FermionPromoter::apply(uint_t icomb, const FrmOps &conn_ops, const FrmOps &com, MaeIndsPartition &mae_inds) const {
    uint_t nexchange = 0ul;
    uint_t icom = 0ul;
    const auto comb = begin(icomb);
    for (uint_t iop = 0ul; iop < m_nexcit; ++iop) {
        auto op = conn_ops[iop];
        // insert all common indices below this excitation index
        for (; icom < m_nop_insert && com[comb[icom]] < op; ++icom) {
            mae_inds[iop + icom] = com[comb[icom]];
            nexchange += iop;
        }
        // insert the excitation index itself
        mae_inds[iop + icom] = op;
    }
    // insert remaining common indices
    for (; icom < m_nop_insert; ++icom) {
        mae_inds[m_nexcit + icom] = com[comb[icom]];
        nexchange += m_nexcit;
    }
    return nexchange;
}

bool FermionPromoter::apply(uint_t icomb, const conn::FrmOnv &conn, const FrmOps &com, MaeIndsPair &frm_inds) const {
    DEBUG_ASSERT_LT(icomb, m_ncomb, "combination index OOB");
    DEBUG_ASSERT_EQ(conn.exsig(), exsig::encode(m_nexcit, m_nexcit, 0, 0),
                    "exsig incompatible with fermion promoter");
    DEBUG_ASSERT_EQ(com.size(), m_ncom, "number of common operators incompatible with fermion promoter");
    const auto nexchange_ann = apply(icomb, conn.m_ann, com, frm_inds.m_ann);
    const auto nexchange_cre = apply(icomb, conn.m_cre, com, frm_inds.m_cre);
    return (nexchange_cre + nexchange_ann) & 1ul;
}
