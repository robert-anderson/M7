//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_SPAWNTABLE_H
#define M7_SPAWNTABLE_H

#include <M7_lib/field/Fields.h>
#include <M7_lib/table/MappedTable.h>

struct Spawn : public Row {
    const bool m_send_parents;
    field::Mbf m_src_mbf;
    field::Mbf m_dst_mbf;
    field::Number<wf_t> m_src_weight;
    field::Number<wf_t> m_delta_weight;
    field::Flag m_src_initiator;
    field::Flag m_src_deterministic;
    field::Number<uint8_t> m_ipart_dst;

    Spawn(const sys::Basis& basis, bool send_parents);
};

#endif //M7_SPAWNTABLE_H
