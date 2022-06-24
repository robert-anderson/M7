//
// Created by Robert J. Anderson on 18/01/2021.
//

#include "WalkerTable.h"

field::Mbf &WalkerTableRow::key_field() {
    return m_mbf;
}

WalkerTableRow::WalkerTableRow(const sys::Basis& basis, uint_t nroot, uint_t nreplica, bool average_weights) :
        m_wf_format({nroot, nreplica}, {"nroot", "nreplica"}),
        m_root_format({nroot}, {"nroot"}),
        m_mbf(this, basis, "many-body basis function"),
        m_weight(this, m_wf_format, "weight"),
        m_hdiag(this, "diagonal H element"),
        m_initiator(this, m_wf_format, "initiator status flag"),
        m_deterministic(this, m_root_format, "deterministic subspace flag"),
        m_ref_conn(this, m_wf_format, "reference connection flag"),
        m_average_weight(average_weights ? this : nullptr, m_wf_format, "unnormalized average weight"),
        m_icycle_occ(average_weights ? this : nullptr, "cycle index at row creation")
{}

bool WalkerTableRow::is_h5_write_exempt() const {
    return m_mbf.is_zero();
}

uint_t WalkerTableRow::occupied_ncycle(const uint_t &icycle_current) const {
    return (1+icycle_current)-m_icycle_occ;
}

const uint_t &WalkerTableRow::nroot() const {
    return m_wf_format.m_shape[0];
}

const uint_t &WalkerTableRow::nreplica() const {
    return m_wf_format.m_shape[1];
}

uint_t WalkerTableRow::ipart_replica(const uint_t &ipart) const {
    return nreplica()==1 ? ipart : (ipart/2)*2+!(ipart&1ul);
}
