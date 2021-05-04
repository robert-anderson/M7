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
    const NdFormat<defs::ndim_wf> m_format;
    fields::Onv<> m_onv;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_weight;
    fields::Number<defs::ham_comp_t> m_hdiag;
    fields::Flags<defs::ndim_wf> m_initiator;
    fields::Flags<defs::ndim_wf> m_reference_connection;
    fields::Flags<defs::ndim_wf> m_deterministic;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_average_weight;
    fields::Number<size_t> m_icycle_occ;

    fields::Onv<> &key_field() {
        return m_onv;
    };

    WalkerTableRow(size_t nsite, size_t nroot, size_t nreplica, bool mev_average_weights=defs::enable_mevs) :
            m_format({nroot, nreplica}, {"nroot", "nreplica"}),
            m_onv(this, nsite, "onv"),
            m_weight(this, m_format, "weight"),
            m_hdiag(this),
            m_initiator(this, m_format),
            m_reference_connection(this, m_format),
            m_deterministic(this, m_format),
            m_average_weight(mev_average_weights ? this : nullptr, m_format),
            m_icycle_occ(mev_average_weights ? this : nullptr)
            {}

    WalkerTableRow(const Options &opts, size_t nsite) :
            WalkerTableRow(nsite, opts.nroot, opts.replicate ? 2:1) {}

    bool is_h5_write_exempt() const override {
        return m_onv.is_zero();
    }

    size_t occupied_ncycle(const size_t& icycle_current) const {
        return icycle_current-m_icycle_occ+1;
    }
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
