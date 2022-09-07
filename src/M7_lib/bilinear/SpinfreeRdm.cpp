//
// Created by anderson on 06/09/2022.
//

#include "SpinfreeRdm.h"

uint_t SpinFreeRdm::spinsig(const MaeIndsPartition& inds, const sys::frm::Size& size) {
    uint_t out = 0ul;
    for (uint_t i = 0ul; i<inds.size(); ++i)
        if (size.ispin(inds[i])) bit::set(out, i);
    return out;
}

uint_t SpinFreeRdm::pair_spinsig(const MaeIndsPair& inds, const sys::frm::Size& size) {
    using namespace spin_free_rdm_arrays;
    const auto spinsig_cre = spinsig(inds.m_cre, size);
    DEBUG_ASSERT_LT(spinsig_cre, nspinsig(inds.m_cre.size()), "spinsig OOB");
    const auto spinsig_ann = spinsig(inds.m_ann, size);
    DEBUG_ASSERT_LT(spinsig_ann, nspinsig(inds.m_ann.size()), "spinsig OOB");
    switch (inds.m_cre.size()) {
        case 2: return c_spinsig_pairs_2[spinsig_cre][spinsig_ann];
        case 3: return c_spinsig_pairs_3[spinsig_cre][spinsig_ann];
    }
    return ~0ul;
}

uint_t SpinFreeRdm::spatsig(const MaeIndsPartition& spat_inds) {
    uint_t out = 0ul;
    for (uint_t i=1ul; i<spat_inds.size(); ++i)
        if (spat_inds[i]==spat_inds[i-1]) bit::set(out, i-1);
    return out;
}

void SpinFreeRdm::make_contribs_from_one_row(const MaeRow& row, wf_t norm) {
    const auto elem = row.m_values[0] / norm;
    auto& given_inds = m_uncontracted_inds;
    const auto orbs = m_basis_size.m_frm;
    const bool is_diagonal = (row.m_inds.m_frm.m_cre==row.m_inds.m_frm.m_ann);
    if (m_nfrm_cre_ind==1ul) {
        /*
         * 1-body case is trivial
         */
        auto add = [&](uint_t p0, uint_t q0, wf_t v) {
            m_uncontracted_inds.m_frm.m_cre[0] = p0;
            m_uncontracted_inds.m_frm.m_ann[0] = q0;
            add_to_send_table(m_uncontracted_inds, v);
        };
        const uint_t i0 = row.m_inds.m_frm.m_cre[0];
        const uint_t j0 = row.m_inds.m_frm.m_ann[0];
        const auto p0 = orbs.isite(i0);
        const auto q0 = orbs.isite(j0);
        DEBUG_ASSERT_EQ(orbs.ispin(i0), orbs.ispin(j0), "spin-tracing of Ms non-conserving RDMs is not valid");
        if (is_diagonal) {
            add(p0, q0, elem);
        }
        else {
            // enforce hermiticity symmetry
            add(p0, q0, 0.5*elem);
            add(q0, p0, 0.5*elem);
        }
    }
    else {
        using namespace spin_free_rdm_arrays;
        /*
         * higher-body RDMs are more complicated to spin-trace: make use of precomputed arrays
         */
        DEBUG_ASSERT_TRUE(m_nfrm_cre_ind==2 || m_nfrm_cre_ind==3, "unsupported RDM rank");
        const uint_t rank = m_nfrm_cre_ind;
        /*
         * first, convert to spatial indices
         */
        spinorbs_to_spat(row.m_inds.m_frm, given_inds.m_frm, m_basis_size.m_frm);
        /*
         * from these, the spatial signatures can be obtained
         */
        const auto spatsig_cre = spatsig(given_inds.m_frm.m_cre);
        const auto spatsig_ann = spatsig(given_inds.m_frm.m_ann);
        /*
         * also required is the pair spin signature, this is available from the original indices
         */
        const auto pair_spinsig = SpinFreeRdm::pair_spinsig(row.m_inds.m_frm, m_basis_size.m_frm);
        /*
         * from the information extracted so far, we can obtain the spin-tracing case index
         */
        const auto icase = rank==2 ?
                           c_case_map_2[pair_spinsig][spatsig_ann][spatsig_cre]:
                           c_case_map_3[pair_spinsig][spatsig_ann][spatsig_cre];
        /*
         * each case has an associated slice of permutations it must handle, this slice is specified by the
         * permutation offset array
         */
        const auto iperm_end = rank==2 ? c_perm_end_offsets_2[icase]: c_perm_end_offsets_3[icase];
        const auto iperm_begin = icase ? (rank==2 ? c_perm_end_offsets_2[icase-1]: c_perm_end_offsets_3[icase-1]) : 0ul;
        /*
         * then all that must be done is to loop over all required permutations of the spatial indices
         */
        for (uint_t iperm=iperm_begin; iperm<iperm_end; ++iperm){
            /*
             * get the factor by which the contribution must be multiplied
             */
            const auto factor = rank==2 ? c_factors_2[iperm] : c_factors_3[iperm];
            /*
             * get pointers to the current arrays of permutations
             */
            const uint_t* cre_perm = rank==2 ? c_cre_perms_2[iperm] : c_cre_perms_3[iperm];
            const uint_t* ann_perm = rank==2 ? c_ann_perms_2[iperm] : c_ann_perms_3[iperm];
            /*
             * execute the permutation of given_inds, storing the result in the insertion object
             */
            for (uint_t i=0ul; i<rank; ++i) {
                m_insert_inds.m_frm.m_cre[i] = given_inds.m_frm.m_cre[cre_perm[i]];
                m_insert_inds.m_frm.m_ann[i] = given_inds.m_frm.m_ann[ann_perm[i]];
            }
            if (is_diagonal) {
                // diagonal element
                add_to_send_table(m_insert_inds, elem*factor);
            }
            else {
                // off-diagonal element: enforce hermiticity symmetry by averaging
                add_to_send_table(m_insert_inds, 0.5*elem*factor);
                m_insert_inds.m_frm.conjugate();
                add_to_send_table(m_insert_inds, 0.5*elem*factor);
            }
        }
    }
}

SpinFreeRdm::SpinFreeRdm(const Rdm& src, wf_t norm, uint_t nelem_per_comm) :
        Rdm(src.m_ranksig, src.m_indsig, src.m_basis_size, src.m_dist.nblock(), src.m_nelec,
            src.m_store.m_row.m_values.nelement(),
            {src.m_store.m_hwm, src.m_store.m_bw.get_expansion_factor()},
            {nelem_per_comm, 1.0}, "sf_"+src.name()), m_insert_inds(src.m_indsig){
    REQUIRE_EQ_ALL(m_nfrm_cre, m_nfrm_ann, "spin tracing requires fermion number conservation");
    logging::info("computing the normalized spin-trace of {}", src.name());
    auto row = src.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()) {
        make_contribs_from_one_row(row, norm);
        if (row.index() && !(row.index()%nelem_per_comm)){
            // time to communicate
            Rdm::end_cycle();
        }
    }
    Rdm::end_cycle();
}
