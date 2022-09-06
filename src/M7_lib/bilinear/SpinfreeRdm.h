//
// Created by anderson on 06/09/2022.
//

#ifndef M7_SPINFREERDM_H
#define M7_SPINFREERDM_H

#include "Rdm.h"
#include "SpinsigMap.h"
#include "SpinFreeRdmArrays.h"


struct SpinSigCaseIdentifier {
};


class SpinFreeRdm : public Rdm {

    /**
     * indices (fermion spatial) of the element to be inserted
     */
    buffered::MaeInds m_insert_inds;

    static uint_t spinsig(const MaeIndsPartition& inds, const sys::frm::Size& size) {
        uint_t out = 0ul;
        for (uint_t i = 0ul; i<inds.size(); ++i)
            if (size.ispin(inds[i])) bit::set(out, i);
        return out;
    }

    static uint_t spincase(const MaeIndsPair& inds, const sys::frm::Size& size) {
        return spinsig_map[spinsig(inds.m_cre, size)][spinsig(inds.m_ann, size)];
    }

    static uint_t spinorbs_to_spat(const MaeIndsPartition& in, MaeIndsPartition& out, const sys::frm::Size& size) {
        for (uint_t i = 0ul; i<in.size(); ++i) out[i] = size.isite(in[i]);
    }

    static uint_t spinorbs_to_spat(const MaeIndsPair& in, MaeIndsPair& out, const sys::frm::Size& size) {
        spinorbs_to_spat(in.m_cre, out.m_cre, size);
        spinorbs_to_spat(in.m_ann, out.m_ann, size);
    }

    static uint_t spatsig(const MaeIndsPartition& spat_inds) {
        uint_t out = 0ul;
        for (uint_t i=1ul; i<spat_inds.size(); ++i)
            if (spat_inds[i]==spat_inds[i-1]) bit::set(out, i-1);
        return out;
    }

    void make_contribs_from_one_row(const MaeRow& row) {
        const auto& elem = row.m_values[0];
        auto& given_inds = m_uncontracted_inds;
        const auto orbs = m_basis_size.m_frm;
        if (m_nfrm_cre_ind==1ul) {
            auto add = [&](uint_t p0, uint_t q0, wf_t v) {
                m_uncontracted_inds.m_frm.m_cre[0] = p0;
                m_uncontracted_inds.m_frm.m_ann[0] = q0;
                add_to_send_table(m_uncontracted_inds, v);
            };
            const uint_t i0 = row.m_inds.m_frm.m_cre[0];
            const uint_t j0 = row.m_inds.m_frm.m_ann[0];
            const auto p0 = orbs.isite(i0);
            const auto q0 = orbs.isite(j0);
            const auto s0 = orbs.ispin(i0);
            const auto t0 = orbs.ispin(j0);
            DEBUG_ASSERT_EQ(s0, t0, "spin-tracing of Ms non-conserving RDMs is not valid");
            add(p0, q0, elem);
            // enforce hermiticity symmetry
            add(q0, p0, elem);
        }
        else {
            DEBUG_ASSERT_TRUE(m_nfrm_cre_ind==2 || m_nfrm_cre_ind==3, "unsupported RDM rank");
            const auto ispatsig = spatsig();
            const auto iperm = case_map_3rdm[][][]
            const auto& perm_end_offsets = m_nfrm_cre_ind==2 ? perm_end_offsets_2rdm : perm_end_offsets_3rdm;
            case 2ul: {
                auto add = [&](uint_t p0, uint_t p1, uint_t q0, uint_t q1, wf_t v) {
                    m_uncontracted_inds.m_frm.m_cre[0] = p0;
                    m_uncontracted_inds.m_frm.m_cre[1] = p1;
                    m_uncontracted_inds.m_frm.m_ann[0] = q0;
                    m_uncontracted_inds.m_frm.m_ann[1] = q1;
                    add_to_send_table(m_uncontracted_inds, v);
                };

                const uint_t i0 = row.m_inds.m_frm.m_cre[0];
                const uint_t i1 = row.m_inds.m_frm.m_cre[1];
                const uint_t j0 = row.m_inds.m_frm.m_ann[0];
                const uint_t j1 = row.m_inds.m_frm.m_ann[1];
                const auto p0 = orbs.isite(i0);
                const auto p1 = orbs.isite(i1);
                const auto q0 = orbs.isite(j0);
                const auto q1 = orbs.isite(j1);
                const auto s0 = orbs.ispin(i0);
                const auto s1 = orbs.ispin(i1);
                const auto t0 = orbs.ispin(j0);
                const auto t1 = orbs.ispin(j1);
                DEBUG_ASSERT_EQ(s0 + s1, t0 + t1, "spin-tracing of Ms non-conserving RDMs is not valid");
                if (s0 == t0) {
                    if (s1 == s0) {
                        // aaaa or bbbb
                        add(p0, p1, q0, q1, elem);
                    }
                }
                add_to_send_table(inds, elem);
                break;
            }
            case 3ul:
                break;
            default:
                ABORT("rank is out of range for implemented spin tracers");
        }

    }
public:
    /**
     * @param src
     *  source RDM from which the spin-free version is to be created
     * @param nelem_per_comm
     *  number of elements of src RDM to process before performing an all-to-allv
     */
    SpinFreeRdm(const Rdm& src, uint_t nelem_per_comm=4000ul):
            Rdm(src.m_ranksig, src.m_indsig, src.m_basis_size, src.m_dist.nblock(), src.m_nelec,
                src.m_store.m_row.m_values.nelement(),
                {src.m_store.m_hwm, src.m_store.m_bw.get_expansion_factor()},
                {nelem_per_comm, 1.0}, "sf_"+src.name()), m_insert_inds(src.m_indsig){
        REQUIRE_EQ_ALL(m_nfrm_cre, m_nfrm_ann, "spin tracing requires fermion number conservation");
        auto row = src.m_store.m_row;
        for (row.restart(); row.in_range(); row.step()) {
            make_contribs_from_one_row(row);
        }
    }
};



#endif //M7_SPINFREERDM_H
