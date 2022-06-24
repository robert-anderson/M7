//
// Created by Robert John Anderson on 2020-04-02.
//

#ifndef M7_WALKERTABLE_H
#define M7_WALKERTABLE_H

#include <M7_lib/connection/Connections.h>
#include <M7_lib/io/Options.h>
#include <M7_lib/field/Row.h>
#include <M7_lib/table/MappedTable.h>
#include <M7_lib/field/Fields.h>


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
    field::Number<uint_t> m_icycle_occ;

    field::Mbf &key_field();

    WalkerTableRow(const sys::Basis& basis, uint_t nroot, uint_t nreplica, bool average_weights);

    bool is_h5_write_exempt() const override;

    /**
     * if the current cycle index is the cycle on which the row became occupied, then the occupied ncycle is 1, since
     * the instantaneous weight has been summed-into the average already, hence the "1+"
     * @param icycle_current
     *  currently in-progress MC cycle
     * @return
     *  correct normalization for the average weight
     */
    uint_t occupied_ncycle(const uint_t& icycle_current) const;

    const uint_t& nroot() const;

    const uint_t& nreplica() const;

    uint_t ipart_replica(const uint_t& ipart) const;
};

struct UniqueOnvRow : public Row {
    field::Number<uint_t> m_ind;
    field::Number<uint_t> m_nparent;

    field::Number<uint_t> &key_field() {
        return m_ind;
    };

    UniqueOnvRow() : m_ind(this), m_nparent(this) {}
};

struct OnvRow : public Row {
    field::Mbf m_mbf;
    field::Number<uint_t> m_nparent;

    field::Mbf &key_field() {
        return m_mbf;
    };

    OnvRow(const sys::Sector& sector) : m_mbf(this, sector), m_nparent(this) {}
};


typedef MappedTable<WalkerTableRow> WalkerTable;

#endif //M7_WALKERTABLE_H
