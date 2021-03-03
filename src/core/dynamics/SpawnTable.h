//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_SPAWNTABLE_H
#define M7_SPAWNTABLE_H

#include "src/core/field/Fields.h"
#include "src/core/table/MappedTable.h"

struct SpawnTableRow : public Row {
    fields::Onv<> m_dst_onv;
    fields::Number<defs::wf_t> m_delta_weight;
    fields::Flag m_src_initiator;
    fields::Flag m_src_deterministic;
    fields::Number<uint8_t> m_dst_ipart;

    SpawnTableRow(size_t nsite) :
            m_dst_onv(this, nsite), m_delta_weight(this),
            m_src_initiator(this), m_src_deterministic(this), m_dst_ipart(this){}
};

#endif //M7_SPAWNTABLE_H
