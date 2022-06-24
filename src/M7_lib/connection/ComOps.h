//
// Created by Robert J. Anderson on 12/5/21.
//

#ifndef M7_COMOPS_H
#define M7_COMOPS_H

#include "FrmOnvConnection.h"
#include "BosOnvConnection.h"

namespace com_ops {
    using Frm = FrmOps;
    using Bos = BosOps;
    struct FrmBos {
        FrmOps m_frm;
        BosOps m_bos;
        FrmBos(sys::Size size) : m_frm(size.m_frm), m_bos(size.m_bos) {}
    };

    typedef std::tuple<Frm, FrmBos, Bos> mbf_tup_t;

    template<uint_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<c_mbf_type_ind> Mbf;
}


#endif //M7_COMOPS_H
