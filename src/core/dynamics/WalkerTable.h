//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"
#include "src/core/io/Options.h"

struct WalkerTable : public Table {
    fields::Onv<> m_onv;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_weight;
    fields::Number<defs::ham_comp_t> m_hdiag;
private:
    struct WalkerTableFlagSet : FlagSet {
        Flags<defs::ndim_wf> m_initiator;
        Flag m_reference_connection;
        Flag m_deterministic;

        WalkerTableFlagSet(size_t nroot, size_t nreplica) :
                m_initiator(this, "is initiator", nroot, nreplica),
                m_reference_connection(this, "is connected to reference ONV"),
                m_deterministic(this, "is in deterministic subspace") {}
    };
public:
    fields::Flags<WalkerTableFlagSet> m_flags;

    WalkerTable(size_t nsite, size_t nroot, size_t nreplica) :
            m_onv(this, nsite, "occupation number vectors"),
            m_weight(this, "weights", nroot, nreplica),
            m_hdiag(this, "hamiltonian diagonal element"),
            m_flags(this, "flags", {nroot, nreplica}){}

    WalkerTable(const Options &opts, size_t nsite):
    WalkerTable(nsite, opts.nroot, opts.nreplica){}
};

struct WalkerMappedTable : public MappedTable<WalkerTable, fields::Onv<>> {
    WalkerMappedTable(size_t nsite, size_t nroot, size_t nreplica, size_t nbucket):
            MappedTable<WalkerTable, fields::Onv<>>(nbucket, m_onv, {nsite, nroot, nreplica})
    {}
};

#endif //M7_WALKERTABLE_H
