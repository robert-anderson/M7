//
// Created by Robert John Anderson on 2020-02-24.
//

#include "WalkerList.h"

WalkerList::WalkerList(size_t nsite, size_t nrow):
PerforableMappedList<Determinant>(
        WalkerListFields(nsite).m_spec, nrow, WalkerListFields(nsite).idet),
        m_fields(nsite){}

Determinant WalkerList::get_determinant(const size_t &irow) {
    return view<Determinant>(irow, m_fields.idet);
}

NumericView<defs::ham_t> WalkerList::get_weight(const size_t &irow) {
    return view<defs::ham_t>(irow, m_fields.iweight);
}

NumericView<defs::ham_comp_t> WalkerList::get_hdiag(const size_t &irow) {
    return view<defs::ham_comp_t>(irow, m_fields.ihdiag);
}

NumericView<bool> WalkerList::get_flag_reference_connection(const size_t &irow) {
    return view<bool>(m_fields.iflag_reference_connection);
}

NumericView<bool> WalkerList::get_flag_initiator(const size_t &irow) {
    return view<bool>(m_fields.iflag_initiator);
}

NumericView<bool> WalkerList::get_flag_deterministic(const size_t &irow) {
    return view<bool>(m_fields.iflag_deterministic);
}
