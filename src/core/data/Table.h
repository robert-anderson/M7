//
// Created by rja on 21/10/2020.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H

#include "src/core/util/defs.h"
#include "Field.h"

struct TableX {
    size_t m_row_size;
    std::vector<FieldX *> m_fields;
    char* m_data;

    TableX(){
    }

    char *begin() {
        return m_data;
    }

    char *begin(const size_t &irow) {
        return begin() + irow * m_row_size;
    }

    size_t add_field(const FieldX& field){
        return 0;// offset
    }
};



#endif //M7_TABLE_H
