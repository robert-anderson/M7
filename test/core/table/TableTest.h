//
// Created by rja on 19/01/2021.
//

#ifndef M7_TABLE_TEST_H
#define M7_TABLE_TEST_H

#include "src/core/table/MappedTable.h"
#include "src/core/fieldz/FieldsZ.h"

namespace table_test {
    struct DetRow : RowZ {
        fieldsz::Onv<0> m_det;
        DetRow(size_t nsite) : m_det(this, nsite){}
    };
}


#endif //M7_TABLE_TEST_H
