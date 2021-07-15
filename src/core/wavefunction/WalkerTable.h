//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include <src/core/basis/Connections.h>
#include "src/core/io/Options.h"
#include "src/core/field/Row.h"
#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"


struct WalkerTableRow : public Row {
    const NdFormat<defs::ndim_wf> m_wf_format;
    const NdFormat<defs::ndim_root> m_root_format;
    fields::Onv<> m_onv;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_weight;
    fields::Number<defs::ham_comp_t> m_hdiag;
    fields::Flags<defs::ndim_wf> m_initiator;
    fields::Flags<defs::ndim_root> m_deterministic;
    fields::Flags<defs::ndim_root> m_ref_conn;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_average_weight;
    fields::Number<size_t> m_icycle_occ;

    fields::Onv<> &key_field();;

    WalkerTableRow(size_t nsite, size_t nroot, size_t nreplica, bool average_weights);

    bool is_h5_write_exempt() const override;

    size_t occupied_ncycle(const size_t& icycle_current) const;
};

struct UniqueOnvRow : public Row {
    fields::Number<size_t> m_ind;
    fields::Number<size_t> m_nparent;

    fields::Number<size_t> &key_field() {
        return m_ind;
    };

    UniqueOnvRow() : m_ind(this), m_nparent(this) {}
};

struct OnvRow : public Row {
    fields::Onv<> m_onv;
    fields::Number<size_t> m_nparent;

    fields::Onv<> &key_field() {
        return m_onv;
    };

    OnvRow(size_t nsite) : m_onv(this, nsite), m_nparent(this) {}
};


typedef MappedTable<WalkerTableRow> WalkerTable;

#endif //M7_WALKERTABLE_H
