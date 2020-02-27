//
// Created by rja on 27/02/2020.
//

#include "SpawnList.h"


SpawnList::SpawnList(const spec_t &spec, size_t nrow, defs::data_t *data_external) :
        List(spec, nrow, data_external){}

SpawnList::SpawnList(const spec_t &spec, size_t nrow) : SpawnList(spec, nrow, nullptr){}

const SpawnList::spec_t &SpawnList::spec() const {
    return static_cast<const spec_t &>(m_spec);
}

Determinant SpawnList::get_determinant(const size_t &irow) {
    return view<Determinant>(irow, spec().idet);
}

NumericView<defs::ham_t> SpawnList::get_weight(const size_t &irow) {
    return view<defs::ham_t>(irow, spec().iweight);
}

NumericView<bool> SpawnList::get_flag_parent_initiator(const size_t &irow) {
    return view<bool>(spec().iflag_parent_initiator);
}

