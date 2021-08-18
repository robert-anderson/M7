//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_SPAWNTABLE_H
#define M7_SPAWNTABLE_H

#include "src/core/field/Fields.h"
#include "src/core/table/MappedTable.h"

struct SpawnTableRow : public Row {
    const bool m_send_parents;
    field::Mbf m_src_mbf;
    field::Mbf m_dst_mbf;
    field::Number<defs::wf_t> m_src_weight;
    field::Number<defs::wf_t> m_delta_weight;
    field::Flag m_src_initiator;
    field::Flag m_src_deterministic;
    field::Number<uint8_t> m_ipart_dst;

    SpawnTableRow(size_t nsite, bool send_parents) :
            m_send_parents(send_parents),
            m_src_mbf(send_parents ? this : nullptr, nsite, "source MBF"),
            m_dst_mbf(this, nsite, "destination MBF"),
            m_src_weight(send_parents ? this : nullptr, "source weight"),
            m_delta_weight(this, "spawned walker weight"),
            m_src_initiator(this, "source initiator flag"),
            m_src_deterministic(this, "source deterministic flag"),
            m_ipart_dst(this, "WF part index of destination") {}
};

#endif //M7_SPAWNTABLE_H
