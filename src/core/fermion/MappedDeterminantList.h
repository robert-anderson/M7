//
// Created by rja on 07/07/2020.
//

#ifndef M7_MAPPEDDETERMINANTLIST_H
#define M7_MAPPEDDETERMINANTLIST_H


#include <src/core/table/DeterminantField.h>
#include <src/core/list/MappedList.h>

struct MappedDeterminantList : public MappedList<DeterminantElement>{
    DeterminantField m_determinant;
    MappedDeterminantList(size_t nelement, size_t nsite, size_t nbucket);
};


#endif //M7_MAPPEDDETERMINANTLIST_H
