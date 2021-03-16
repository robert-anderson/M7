//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include "src/core/io/Options.h"
#include "src/core/field/Row.h"
#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"


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

struct UniqueOnvRow : public Row {
    fields::Number<size_t> m_ind;
    fields::Number<size_t> m_nparent;

    fields::Number<size_t> &key_field() {
        return m_ind;
    };

    UniqueOnvRow(): m_ind(this), m_nparent(this){}
};

struct OnvRow : public Row {
    fields::Onv<> m_onv;
    fields::Number<size_t> m_nparent;

    fields::Onv<> &key_field() {
        return m_onv;
    };

    OnvRow(size_t nsite): m_onv(this, nsite), m_nparent(this){}
};



typedef MappedTable<WalkerTableRow> WalkerTable;

#endif //M7_WALKERTABLE_H
