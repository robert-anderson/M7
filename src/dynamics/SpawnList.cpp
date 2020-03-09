//
// Created by rja on 27/02/2020.
//

#include "SpawnList.h"

SpawnList::SpawnList(size_t nsite, defs::data_t *data_external) :
        List(data_external),
        m_determinant(this, nsite),
        m_weight(this),
        m_parent_initiator(this) {}
