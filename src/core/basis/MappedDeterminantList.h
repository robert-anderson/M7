//
// Created by rja on 07/07/2020.
//

#ifndef M7_MAPPEDDETERMINANTLIST_H
#define M7_MAPPEDDETERMINANTLIST_H


#include <src/core/basis/DeterminantField.h>
#include <src/core/list/MappedList.h>

struct MappedDeterminantList : public MappedList<DeterminantElement>{
    DeterminantField m_determinant;
    MappedDeterminantList(std::string name, size_t nelement, size_t nsite, size_t nbucket);
};


#endif //M7_MAPPEDDETERMINANTLIST_H
