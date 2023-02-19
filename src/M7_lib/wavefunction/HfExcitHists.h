//
// Created by rja on 19/02/23.
//

#ifndef M7_HFEXCITHISTS_H
#define M7_HFEXCITHISTS_H


#include "M7_lib/communication/SharedRows.h"

#if 0

struct HfExcitHist {
    const OpSig m_exsig;
    HfExcitHist(uint_t nexcit): m_exsig({nexcit, nexcit}, {0, 0}){}
};
class HfExcitHists {
    const shared_rows::Walker* m_hf;
    const uint_t m_nexcit_max;
    const wf_t m_thresh;

    HfExcitHists(const shared_rows::Walker* hf, uint_t nexcit_max, wf_t thresh):
        m_hf(hf), m_nexcit_max(nexcit_max), m_thresh(thresh){}

    void accumulate(const field::Mbf& mbf, wf_t weight) {
        if (weight < m_thresh)
    }

};


#endif //M7_HFEXCITHISTS_H
#endif //M7_HFEXCITHISTS_H
