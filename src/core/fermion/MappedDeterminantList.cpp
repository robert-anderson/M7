//
// Created by rja on 07/07/2020.
//

#include "MappedDeterminantList.h"

MappedDeterminantList::MappedDeterminantList(std::string name, size_t nelement, size_t nsite, size_t nbucket) :
        MappedList<DeterminantElement>(name, m_determinant, nbucket), m_determinant(this, nelement, nsite){}
