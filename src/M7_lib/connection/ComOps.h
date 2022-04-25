//
// Created by anderson on 12/5/21.
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
        FrmBos(sys::Size extents) : m_frm(extents.m_sites), m_bos(extents.m_nmode) {}
    };

    typedef std::tuple<Frm, FrmBos, Bos> mbf_tup_t;

    template<size_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<defs::mbf_type_ind> Mbf;
}


#endif //M7_COMOPS_H
