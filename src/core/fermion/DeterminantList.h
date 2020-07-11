//
// Created by rja on 07/07/2020.
//

#ifndef M7_DETERMINANTLIST_H
#define M7_DETERMINANTLIST_H


#include <src/core/list/List.h>
#include <src/core/table/DeterminantField.h>

struct DeterminantList : public List {
    DeterminantField m_determinant;
    DeterminantList(std::string name, const size_t& nelement, const size_t& nsite, size_t nsegment = 1);
};


#endif //M7_DETERMINANTLIST_H
