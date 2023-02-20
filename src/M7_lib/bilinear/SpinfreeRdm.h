//
// Created by anderson on 06/09/2022.
//

#ifndef M7_SPINFREERDM_H
#define M7_SPINFREERDM_H

#include "Rdm.h"
#include "SpinFreeRdmArrays.h"


class SpinFreeRdm : public Rdm {

    /**
     * indices (fermion spatial) of the element to be inserted
     */
    buffered::RdmInds m_insert_inds;

    static uint_t nspinsig(uint_t rank){
        using namespace spin_free_rdm_arrays;
        return rank==2 ? c_nspinsig_2 : c_nspinsig_3;

    }
    static uint_t spinsig(const RdmIndsPartition& inds, const sys::frm::Size& size);

    static uint_t pair_spinsig(const RdmIndsPair& inds, const sys::frm::Size& size);

    static void spinorbs_to_spat(const RdmIndsPartition& in, RdmIndsPartition& out, const sys::frm::Size& size) {
        for (uint_t i = 0ul; i<in.size(); ++i) out[i] = size.isite(in[i]);
    }

    static void spinorbs_to_spat(const RdmIndsPair& in, RdmIndsPair& out, const sys::frm::Size& size) {
        spinorbs_to_spat(in.m_cre, out.m_cre, size);
        spinorbs_to_spat(in.m_ann, out.m_ann, size);
    }

    /**
     * there are 2 possible distinct spatial signatures for the cre/ann strings of the 2-RDM:
     *  0: spatial indices are distinct
     *  1: spatial indices are the same
     * and for the 3-RDM, there are 3:
     *  0: all spatial indices are distinct
     *  1: first two are the same, last one distinct
     *  2: last two are the same, first one distinct
     * there is no spatial signature 3 for the 3-RDM since there would have to be two electrons in the same spin-orbital
     * @param spat_inds
     *  spatial indices as single cre/ann string
     * @return
     *  the spatial signature
     */
    static uint_t spatsig(const RdmIndsPartition& spat_inds);

    void make_contribs_from_one_row(const RdmRow& row, wf_t norm);
public:
    /**
     * @param src
     *  source RDM from which the spin-free version is to be created
     * @param nelem_per_comm
     *  number of elements of src RDM to process before performing an all-to-allv
     */
    SpinFreeRdm(const Rdm& src, wf_t norm, uint_t nelem_per_comm=4000ul);
};



#endif //M7_SPINFREERDM_H
