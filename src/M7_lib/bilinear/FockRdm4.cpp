//
// Created by rja on 03/11/22.
//

#include "FockRdm4.h"


FockRdm4::FockRdm4(const conf::Rdms &opts, const FockMatrix& fock, sys::Sector sector, uint_t nvalue) :
        Rdm(opts, exsig::ex_4400, exsig::ex_3300, sector, nvalue, "4400f"), m_fock(fock){
    REQUIRE_EQ(fock.nrow(), sector.m_frm.size(), "Incorrectly-sized Fock matrix given");
    if (m_fock.is_diagonal())
        logging::warn("Performing 4RDM contraction of a diagonal Fock matrix without exploiting diagonality");
    //logging::info("The given Fock matrix was found to be {}diagonal", (m_diagonal ? "":"non-"));
    //if (m_diagonal) logging::info("Quadruple-fermion (4400) excitations do not contribute to Fock * 4RDM");
    /*
     *
        m_diagonal = m_fock.is_diagonal();
        logging::info("The given Fock matrix was found to be {}diagonal", (m_diagonal ? "":"non-"));
        if (m_diagonal) logging::info("Quadruple-fermion (4400) excitations do not contribute to Fock * 4RDM");
     */
}

void FockRdm4::frm_make_contribs(const FrmOnv &src_onv, const conn::FrmOnv &conn,
                                 const FrmOps &com, const wf_t &contrib) {
    const auto exlvl = conn.m_cre.size();
    DEBUG_ASSERT_TRUE(conn.m_ann.size() <= m_nfrm_ann && conn.m_cre.size() <= m_nfrm_cre,
                      "this method should not have been delegated given the exsig of the contribution");
    /*
     * number of "inserted" fermion creation/annihilation operator pairs
     */
    const auto nins = m_rank - exlvl;
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
        // TODO: diagonal Fock optimisation
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