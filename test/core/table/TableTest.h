//
// Created by rja on 19/01/2021.
//

#ifndef M7_TABLE_TEST_H
#define M7_TABLE_TEST_H

#include "src/core/table/MappedTable.h"

namespace table_test {

    struct DetTable : Table {
        fields::Det m_config;
        DetTable(size_t nsite) :
                m_config(this, nsite, "configuration") {}
    };

    struct DetMappedTable : MappedTable<DetTable, fields::Det> {
        DetMappedTable(size_t nsite):
                MappedTable<DetTable, fields::Det>(m_config, nsite){}
    };
}


#endif //M7_TABLE_TEST_H
