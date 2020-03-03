//
// Created by rja on 27/02/2020.
//

#include "SpawnList.h"


SpawnList::SpawnList(size_t nsite, size_t nrow, defs::data_t *data_external) :
        List(SpawnListFields(nsite).m_spec, nrow, data_external), m_fields(SpawnListFields(nsite)) {}

SpawnList::SpawnList(size_t nsite, size_t nrow) :
        List(SpawnListFields(nsite).m_spec, nrow), m_fields(SpawnListFields(nsite)) {}

Determinant SpawnList::get_determinant(const size_t &irow) {
    return view<Determinant>(irow, m_fields.idet);
}

NumericView<defs::ham_t> SpawnList::get_weight(const size_t &irow) {
    return view<defs::ham_t>(irow, m_fields.iweight);
}

NumericView<bool> SpawnList::get_flag_parent_initiator(const size_t &irow) {
    return view<bool>(m_fields.iflag_parent_initiator);
}
