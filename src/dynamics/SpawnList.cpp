//
// Created by rja on 27/02/2020.
//

#include "SpawnList.h"


SpawnList::SpawnList(const spec_T &spec, size_t nrow, defs::data_t *data_external) :
        List(spec, nrow, data_external), m_spec(spec) {}

SpawnList::SpawnList(const spec_T &spec, size_t nrow) :
        List(spec, nrow), m_spec(spec) {}

Determinant SpawnList::get_determinant(const size_t &irow) {
    return view<Determinant>(irow, m_spec.idet);
}

NumericView<defs::ham_t> SpawnList::get_weight(const size_t &irow) {
    return view<defs::ham_t>(irow, m_spec.iweight);
}

NumericView<bool> SpawnList::get_flag_parent_initiator(const size_t &irow) {
    return view<bool>(m_spec.iflag_parent_initiator);
}

