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
            m_src_mbf(send_parents ? this : nullptr, nsite),
            m_dst_mbf(this, nsite),
            m_src_weight(send_parents ? this : nullptr),
            m_delta_weight(this),
            m_src_initiator(this),
            m_src_deterministic(this),
            m_ipart_dst(this) {}
};

#endif //M7_SPAWNTABLE_H
