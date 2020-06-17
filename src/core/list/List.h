//
// Created by Robert John Anderson on 2020-03-31.
//

#ifndef M7_LIST_H
#define M7_LIST_H

#include <src/core/parallel/MPIWrapper.h>
#include "src/core/table/Table.h"

class List : public Table {
    List *m_recv = nullptr;
    defs::inds m_high_water_mark;
public:
    explicit List(size_t nsegment = 1);

    void recv(List *list);

    void expand(size_t delta_nrow) override;

    const defs::inds &high_water_mark() const;

    const size_t &high_water_mark(const size_t isegment) const;

    virtual size_t push(const size_t &isegment = 0);

    size_t push(const size_t &isegment, const size_t &nrow);

    size_t expand_push(const size_t &isegment, const size_t &nrow, double factor = 1.5);

    size_t expand_push(const size_t &isegment=0, double factor = 1.5){
        return expand_push(isegment, 1, factor);
    }

    void zero() override;

    std::string to_string();

    void communicate();

    void all_gather(List &local);
};


#endif //M7_LIST_H
