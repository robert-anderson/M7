//
// Created by Robert John Anderson on 2020-02-24.
//

#include "WalkerList.h"


WalkerListSpecification::WalkerListSpecification(size_t nsite) : Specification(){
    idet = add<Determinant>(nsite);
    iweight = add<defs::ham_t>(1);
    ihdiag = add<defs::ham_comp_t>(1);
    iflag_reference_connection = add<bool>(1);
    iflag_initiator = add<bool>(1);
    iflag_deterministic = add<bool>(1);
}

const WalkerListSpecification &WalkerList::spec() const {
    return static_cast<const WalkerListSpecification &>(m_spec);
}

Determinant WalkerList::get_determinant(const size_t &irow) {
    return view<Determinant>(irow, spec().idet);
}

NumericView<defs::ham_t> WalkerList::get_weight(const size_t &irow) {
    return view<defs::ham_t>(irow, spec().iweight);
}

NumericView<defs::ham_comp_t> WalkerList::get_hdiag(const size_t &irow) {
    return view<defs::ham_comp_t>(irow, spec().ihdiag);
}

NumericView<bool> WalkerList::get_flag_reference_connection(const size_t &irow) {
    return view<bool>(spec().iflag_reference_connection);
}

NumericView<bool> WalkerList::get_flag_initiator(const size_t &irow) {
    return view<bool>(spec().iflag_initiator);
}

NumericView<bool> WalkerList::get_flag_deterministic(const size_t &irow) {
    return view<bool>(spec().iflag_deterministic);
}

WalkerList::WalkerList(const spec_t &spec, size_t nrow) :
PerforableMappedList<Determinant>(spec, nrow, spec.idet) {}

