//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"
#include "src/core/io/Options.h"
#include "src/core/fieldz/RowZ.h"
#include "src/core/fieldz/FieldsZ.h"
#include "src/core/fieldz/MappedTableZ.h"


struct WalkerTableRow : public RowZ {
    fieldsz::Onv<> m_onv;
    fieldsz::Numbers<defs::wf_t, defs::ndim_wf> m_weight;
    fieldsz::Number<defs::ham_comp_t> m_hdiag;
    fieldsz::Flags<defs::ndim_wf> m_initiator;
    fieldsz::Flags<defs::ndim_wf> m_reference_connection;
    fieldsz::Flags<defs::ndim_wf> m_deterministic;
    fieldsz::Onv<> &m_key_field = m_onv;

    WalkerTableRow(size_t nsite, size_t nroot, size_t nreplica):
    m_onv(this, nsite),
    m_weight(this, {nroot, nreplica}),
    m_hdiag(this),
    m_initiator(this, {nroot, nreplica}),
    m_reference_connection(this, {nroot, nreplica}),
    m_deterministic(this, {nroot, nreplica}){}

    WalkerTableRow(const Options &opts, size_t nsite):
            WalkerTableRow(nsite, opts.nroot, opts.nreplica){}
};

typedef MappedTableZ<WalkerTableRow> WalkerTable;

#endif //M7_WALKERTABLE_H
