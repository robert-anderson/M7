//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_SPAWNTABLE_H
#define M7_SPAWNTABLE_H

struct SpawnTableRow : public RowZ {
    fieldsz::Onv<> m_dst_onv;
    fieldsz::Number<defs::wf_t> m_delta_weight;
    fieldsz::Flag m_src_initiator;
    fieldsz::Flag m_src_deterministic;
    fieldsz::Number<uint8_t> m_dst_ipart;

    SpawnTableRow(size_t nsite) :
            m_dst_onv(this, nsite), m_delta_weight(this),
            m_src_initiator(this), m_src_deterministic(this), m_dst_ipart(this){}
};

typedef MappedTableZ<SpawnTableRow> SpawnTable;

#endif //M7_SPAWNTABLE_H
