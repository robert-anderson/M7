//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include "src/core/io/Options.h"
#include "src/core/fieldz/Row.h"
#include "src/core/fieldz/MappedTable.h"
#include "src/core/fieldz/Fields.h"


struct WalkerTableRow : public Row {
    const size_t m_npart;
    fields::Onv<> m_onv;
    fields::Vector<defs::wf_t> m_weight;
    fields::Number<defs::ham_comp_t> m_hdiag;
    fields::Flags m_initiator;
    fields::Flags m_reference_connection;
    fields::Flags m_deterministic;

    fields::Onv<> &key_field() {
        return m_onv;
    };

    WalkerTableRow(size_t nsite, size_t nroot, size_t nreplica):
    m_npart(nroot*nreplica),
    m_onv(this, nsite),
    m_weight(this, m_npart),
    m_hdiag(this),
    m_initiator(this, m_npart),
    m_reference_connection(this, m_npart),
    m_deterministic(this, m_npart){}

    WalkerTableRow(const Options &opts, size_t nsite):
            WalkerTableRow(nsite, opts.nroot, opts.nreplica){}
};

typedef MappedTable<WalkerTableRow> WalkerTable;

#endif //M7_WALKERTABLE_H
