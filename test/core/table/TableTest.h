//
// Created by rja on 19/01/2021.
//

#ifndef M7_TABLE_TEST_H
#define M7_TABLE_TEST_H

#include "src/core/field/Fields.h"

namespace table_test {
    struct DetRow : Row {
        field::FrmOnv m_det;
        DetRow(size_t nsite) : m_det(this, nsite){}
    };
}


#endif //M7_TABLE_TEST_H
