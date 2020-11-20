//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_SPAWNTABLE_H
#define M7_SPAWNTABLE_H

#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"

struct SpawnTable : public TableX {
    fields::Onv m_dst_onv;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_delta_weight;
private:
    struct SpawnTableFlagSet : FlagSet {
        Flags<defs::ndim_wf> m_src_initiator;
        Flag m_src_deterministic;

        SpawnTableFlagSet(fields::Bitset *bitset, size_t nroot, size_t nreplica) :
                FlagSet(bitset),
                m_src_initiator(this, "parent is initiator", nroot, nreplica),
                m_src_deterministic(this, "parent is in deterministic subspace") {}
    };

public:
    fields::Flags<SpawnTableFlagSet> m_flags;

    SpawnTable(size_t nsite, size_t nroot, size_t nreplica) :
            m_dst_onv(this, nsite, "occupation number vectors"),
            m_delta_weight(this, "weights", nroot, nreplica),
            m_flags(this, "flags", nroot, nreplica) {}
};

#endif //M7_SPAWNTABLE_H
