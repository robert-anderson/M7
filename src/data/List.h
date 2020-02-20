//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_LIST_H
#define M7_LIST_H

#include "Table.h"

class List : public Table {
    size_t m_high_water_mark{};

public:
    List(Specification spec, size_t nrow) : Table(spec, nrow) {}

    List(Specification spec, size_t nrow, defs::data_t *data_external) : Table(spec, nrow, data_external) {}

    virtual const size_t high_water_mark() const {
        return m_high_water_mark;
    }

    size_t push() {
        size_t tmp;
#pragma omp atomic capture
        tmp = m_high_water_mark++;
        return tmp;
    }

    size_t push(const size_t &irow) {
        size_t tmp;
#pragma omp atomic capture
        tmp = m_high_water_mark += irow;
        return tmp;
    }

};


#endif //M7_LIST_H
