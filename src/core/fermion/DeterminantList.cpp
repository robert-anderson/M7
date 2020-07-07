//
// Created by rja on 07/07/2020.
//

#include "DeterminantList.h"

DeterminantList::DeterminantList(const size_t &nelement, const size_t &nsite, size_t nsegment) :
        m_determinant(this, nelement, nsite, "determinant"){}
