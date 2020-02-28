//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_LIST_H
#define M7_LIST_H

#include "Table.h"

class List : public Table {
    size_t m_high_water_mark{};

public:
    List(const spec_T &spec, size_t nrow);

    List(const spec_T &spec, size_t nrow, defs::data_t *data_external);

    virtual size_t high_water_mark() const;

    void high_water_mark(const size_t& value);

    size_t push();

    size_t push(const size_t &nrow);

    void zero() override;

    void print() const override ;

};


#endif //M7_LIST_H
