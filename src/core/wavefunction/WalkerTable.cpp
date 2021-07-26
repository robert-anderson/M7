//
// Created by rja on 18/01/2021.
//

#include "WalkerTable.h"

fields::mbf_t &WalkerTableRow::key_field() {
    return m_mbf;
}

WalkerTableRow::WalkerTableRow(size_t nsite, size_t nroot, size_t nreplica, bool average_weights) :
        m_wf_format({nroot, nreplica}, {"nroot", "nreplica"}),
        m_root_format({nroot}, {"nroot"}),
        m_mbf(this, nsite, "onv"),
        m_weight(this, m_wf_format, "weight"),
        m_hdiag(this, "diagonal H element"),
        m_initiator(this, m_wf_format, "initiator status flag"),
        m_deterministic(this, m_root_format, "deterministic subspace flag"),
        m_ref_conn(this, m_root_format, "reference connection flag"),
        m_average_weight(average_weights ? this : nullptr, m_wf_format, "unnormalized average weight"),
        m_icycle_occ(average_weights ? this : nullptr, "cycle index at row creation")
{}

bool WalkerTableRow::is_h5_write_exempt() const {
    return m_mbf.is_zero();
}

size_t WalkerTableRow::occupied_ncycle(const size_t &icycle_current) const {
    return icycle_current-m_icycle_occ;
}
