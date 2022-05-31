//
// Created by Robert J. Anderson on 18/08/2021.
//

#include "SpawnTable.h"

SpawnTableRow::SpawnTableRow(const sys::Basis& basis, bool send_parents) :
        m_send_parents(send_parents),
        m_src_mbf(send_parents ? this : nullptr, basis, "source MBF"),
        m_dst_mbf(this, basis, "destination MBF"),
        m_src_weight(send_parents ? this : nullptr, "source weight"),
        m_delta_weight(this, "spawned walker weight"),
        m_src_initiator(this, "source initiator flag"),
        m_src_deterministic(this, "source deterministic flag"),
        m_ipart_dst(this, "WF part index of destination") {}
