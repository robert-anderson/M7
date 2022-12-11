//
// Created by rja on 03/11/22.
//

#include "FockRdm4.h"


FockRdm4::FockRdm4(const conf::Rdms &opts, OpSig max_contrib_exsig, sys::Sector sector, uint_t nvalue) :
        ContractedRdm(opts, opsig::c_4400, opsig::c_3300, max_contrib_exsig, sector, nvalue, "4400f"){}


NonDiagFockRdm4::NonDiagFockRdm4(const conf::Rdms &opts, const FockMatrix& fock, sys::Sector sector, uint_t nvalue) :
        FockRdm4(opts, opsig::c_4400, sector, nvalue), m_fock(fock){
    REQUIRE_EQ(fock.nrow(), sector.m_frm.size(), "Incorrectly-sized Fock matrix given");
    if (m_fock.is_diagonal())
        logging::warn("Performing 4RDM contraction of a diagonal Fock matrix without exploiting diagonality");
}

void NonDiagFockRdm4::frm_make_contribs(const FrmOnv &src_onv, const conn::FrmOnv &conn,
                                        const FrmOps &com, wf_t contrib) {
    const auto exlvl = conn.m_cre.size();
    DEBUG_ASSERT_TRUE(conn.m_ann.size() <= m_nfrm_ann && conn.m_cre.size() <= m_nfrm_cre,
                      "this method should not have been delegated given the exsig of the contribution");
    /*
     * number of "inserted" fermion creation/annihilation operator pairs
     */
    const auto rank = m_ranksig.nfrm_cre();
    const auto nins = rank - exlvl;
    /*
     * this determines the precomputed promoter required
     */
    const auto& promoter = m_frm_promoters[nins];
    /*
     * apply each combination of the promoter deterministically
     */
    for (uint_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
        auto phase = promoter.apply(icomb, conn, com, m_full_inds.m_frm);
        /*
         * include the Fermi phase of the excitation
         */
        phase = phase ^ conn.phase(src_onv);
        for (uint_t icre_contract = 0ul; icre_contract < 4ul; ++icre_contract){
            const uint_t icre = m_full_inds.m_frm.m_cre[icre_contract];
            const auto isite_cre = src_onv.m_basis.isite(icre);
            const auto ispin_cre = src_onv.m_basis.ispin(icre);
            for (uint_t iann_contract = 0ul; iann_contract < 4ul; ++iann_contract) {
                const uint_t iann = m_full_inds.m_frm.m_ann[iann_contract];
                const auto isite_ann = src_onv.m_basis.isite(iann);
                const auto ispin_ann = src_onv.m_basis.ispin(iann);

                if (ispin_ann != ispin_cre) continue;

                const auto fock_element = m_fock(isite_cre, isite_ann);
                if (!ham::is_significant(fock_element)) continue;
                /*
                 * fill the uncontracted indices (those identifying the elements of the intermediate)
                 */
                for (uint_t iuncontract = 0ul; iuncontract < 4ul; ++iuncontract) {
                    // the next position in the uncontracted indices to fill
                    uint_t i;
                    if (iuncontract!=icre_contract) {
                        i = iuncontract - (iuncontract > icre_contract ? 1 : 0);
                        m_uncontracted_inds.m_frm.m_cre[i] = m_full_inds.m_frm.m_cre[iuncontract];
                    }
                    if (iuncontract!=iann_contract) {
                        i = iuncontract - (iuncontract > iann_contract ? 1 : 0);
                        m_uncontracted_inds.m_frm.m_ann[i] = m_full_inds.m_frm.m_ann[iuncontract];
                    }
                }

                /*
                 * include the Fermi phase of the contraction and the promotion when computing the overall sign
                 */
                const bool contract_phase = (icre_contract+iann_contract)&1ul;
                add_to_send_table(m_uncontracted_inds,
                                  ((phase ^ contract_phase) ? -contrib : contrib) * fock_element);
            }
        }
    }
}

DiagFockRdm4::DiagFockRdm4(const conf::Rdms& opts, const FockMatrix& fock, sys::Sector sector, uint_t nvalue) :
        FockRdm4(opts, opsig::c_trip, sector, nvalue), m_fock(fock.get_diagonal()){}

void DiagFockRdm4::frm_make_contribs(const FrmOnv& src_onv, const conn::FrmOnv& conn,
                                     const FrmOps& com, wf_t contrib) {
    const auto exlvl = conn.m_cre.size();
    DEBUG_ASSERT_TRUE(conn.m_ann.size() <= m_nfrm_ann && conn.m_cre.size() <= m_nfrm_cre,
                      "this method should not have been delegated given the exsig of the contribution");
    /*
     * quadruples do not contribute
     */
    if (exlvl == 4) return;
    const auto rank = m_ranksig.nfrm_cre();
    /*
     * number of "inserted" fermion creation/annihilation operator pairs
     */
    const auto nins = rank - exlvl;
    /*
     * this determines the precomputed promoter required
     */
    const auto& promoter = m_frm_promoters[nins];
    /*
     * apply each combination of the promoter deterministically
     * starting with the diagonal (Ci * Ci) contribution
     */
    if (!exlvl) {
        for (uint_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
            promoter.apply(icomb, conn, com, m_full_inds.m_frm);
            /*
             * connection and promotion-related Fermi phases are always false
             */
            for (uint_t icontract = 0ul; icontract < 4ul; ++icontract) {
                const uint_t iop = m_full_inds.m_frm.m_cre[icontract];
                REQUIRE_EQ(m_full_inds.m_frm.m_ann[icontract], iop, "not a diagonal contribution");
                const auto isite = src_onv.m_basis.isite(iop);
                const auto fock_element = m_fock[isite];
                if (!ham::is_significant(fock_element)) continue;
                for (uint_t iuncontract = 0ul; iuncontract < 4ul; ++iuncontract) {
                    // the next position in the uncontracted indices to fill
                    uint_t i;
                    if (iuncontract!=icontract) {
                        i = iuncontract - (iuncontract > icontract ? 1 : 0);
                        const uint_t op = m_full_inds.m_frm.m_cre[iuncontract];
                        m_uncontracted_inds.m_frm.m_cre[i] = op;
                        m_uncontracted_inds.m_frm.m_ann[i] = op;
                    }
                }
                add_to_send_table(m_uncontracted_inds, contrib * fock_element);
            }
        }
    }
    else {
        for (uint_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
            auto phase = promoter.apply(icomb, conn, com, m_full_inds.m_frm);
            /*
             * include the Fermi phase of the excitation
             */
            phase = phase ^ conn.phase(src_onv);
            for (uint_t icre_contract = 0ul; icre_contract < 4ul; ++icre_contract){
                const uint_t icre = m_full_inds.m_frm.m_cre[icre_contract];
                const auto isite = src_onv.m_basis.isite(icre);
                for (uint_t iann_contract = 0ul; iann_contract < 4ul; ++iann_contract) {
                    const uint_t iann = m_full_inds.m_frm.m_ann[iann_contract];
                    /*
                     * this is the diagonal Fock matrix optimization
                     */
                    if (iann != icre) continue;

                    const auto fock_element = m_fock[isite];
                    if (!ham::is_significant(fock_element)) continue;
                    /*
                     * fill the uncontracted indices (those identifying the elements of the intermediate)
                     */
                    for (uint_t iuncontract = 0ul; iuncontract < 4ul; ++iuncontract) {
                        // the next position in the uncontracted indices to fill
                        uint_t i;
                        if (iuncontract!=icre_contract) {
                            i = iuncontract - (iuncontract > icre_contract ? 1 : 0);
                            m_uncontracted_inds.m_frm.m_cre[i] = m_full_inds.m_frm.m_cre[iuncontract];
                        }
                        if (iuncontract!=iann_contract) {
                            i = iuncontract - (iuncontract > iann_contract ? 1 : 0);
                            m_uncontracted_inds.m_frm.m_ann[i] = m_full_inds.m_frm.m_ann[iuncontract];
                        }
                    }

                    /*
                     * include the Fermi phase of the contraction and the promotion when computing the overall sign
                     */
                    const bool contract_phase = (icre_contract+iann_contract)&1ul;
                    add_to_send_table(m_uncontracted_inds,
                                      ((phase ^ contract_phase) ? -contrib : contrib) * fock_element);
                }
            }
        }
    }
}
