//
// Created by rja on 07/07/2020.
//

#include "DeterminantList.h"

DeterminantList::DeterminantList(std::string name, const size_t &nelement, const size_t &nsite, size_t nsegment) :
        List(name), m_determinant(this, nelement, nsite, "determinant") {}
