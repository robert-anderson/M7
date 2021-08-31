//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include <src/core/connection/Connections.h>
#include "src/core/io/Options.h"
#include "src/core/field/Row.h"
#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"


struct WalkerTableRow : public Row {
    const NdFormat<defs::ndim_wf> m_wf_format;
    const NdFormat<defs::ndim_root> m_root_format;
    field::Mbf m_mbf;
    field::Numbers<defs::wf_t, defs::ndim_wf> m_weight;
    field::Number<defs::ham_comp_t> m_hdiag;
    field::Flags<defs::ndim_wf> m_initiator;
    field::Flags<defs::ndim_root> m_deterministic;
    field::Flags<defs::ndim_wf> m_ref_conn;
    field::Numbers<defs::wf_t, defs::ndim_wf> m_average_weight;
    field::Number<size_t> m_icycle_occ;

    field::Mbf &key_field();

    WalkerTableRow(size_t nsite, size_t nroot, size_t nreplica, bool average_weights);

    bool is_h5_write_exempt() const override;

    /**
     * if the current cycle index is the cycle on which the row became occupied, then the occupied ncycle is 1, since
     * the instantaneous weight has been summed-into the average already, hence the "1+"
     * @param icycle_current
     *  currently in-progress MC cycle
     * @return
     *  correct normalization for the average weight
     */
    size_t occupied_ncycle(const size_t& icycle_current) const;

    const size_t& nroot() const;

    const size_t& nreplica() const;

    size_t ipart_replica(const size_t& ipart) const;
};

struct UniqueOnvRow : public Row {
    field::Number<size_t> m_ind;
    field::Number<size_t> m_nparent;

    field::Number<size_t> &key_field() {
        return m_ind;
    };

    UniqueOnvRow() : m_ind(this), m_nparent(this) {}
};

struct OnvRow : public Row {
    field::Mbf m_mbf;
    field::Number<size_t> m_nparent;

    field::Mbf &key_field() {
        return m_mbf;
    };

    OnvRow(size_t nsite) : m_mbf(this, nsite), m_nparent(this) {}
};


typedef MappedTable<WalkerTableRow> WalkerTable;

#endif //M7_WALKERTABLE_H
