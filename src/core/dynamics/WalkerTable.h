//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include "src/core/io/Options.h"
#include "src/core/fieldz/RowZ.h"
#include "src/core/fieldz/MappedTableZ.h"
#include "src/core/fieldz/FieldsZ.h"


struct WalkerTableRow : public RowZ {
    const size_t m_npart;
    fieldsz::Onv<> m_onv;
    fieldsz::Vector<defs::wf_t> m_weight;
    fieldsz::Number<defs::ham_comp_t> m_hdiag;
    fieldsz::Flags m_initiator;
    fieldsz::Flags m_reference_connection;
    fieldsz::Flags m_deterministic;

    fieldsz::Onv<> &key_field() {
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

typedef MappedTableZ<WalkerTableRow> WalkerTable;

#endif //M7_WALKERTABLE_H
