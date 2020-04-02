//
// Created by Robert John Anderson on 2020-03-31.
//

#ifndef M7_LIST_H
#define M7_LIST_H

#include "src/core/table/Table.h"

class List : public Table {
    defs::inds m_high_water_mark;
public:
    List(size_t nsegment=1);

    const defs::inds &high_water_mark() const;

    virtual size_t push(const size_t &isegment = 0);

    size_t push(const size_t &isegment, const size_t &nrow);

    void zero() override;
};


#endif //M7_LIST_H
