//
// Created by Robert John Anderson on 2020-03-31.
//

#ifndef M7_LIST_H
#define M7_LIST_H

#include <src/core/parallel/MPIWrapper.h>
#include <src/core/sort/QuickSorter.h>
#include "src/core/table/Table.h"

class List : public Table {
    List *m_recv = nullptr;
    defs::inds m_high_water_mark;
public:
    explicit List(std::string name, size_t nsegment = 1);

    void recv(List *list);

    void expand(size_t delta_nrow) override;

    const defs::inds &high_water_mark() const;

    const size_t &high_water_mark(const size_t isegment) const;

    size_t push(const size_t &isegment=0, const size_t &nrow=1);

    size_t expand_push(const size_t &isegment=0, const size_t &nrow=1, double factor = 1.5);

    void zero() override;

    virtual std::string to_string() const;

    virtual std::string to_string(const defs::inds &nrows) const;

    void communicate();

    void all_gather(List &local);

};

#endif //M7_LIST_H
