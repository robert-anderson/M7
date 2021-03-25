//
// Created by rja on 29/11/2020.
//

#ifndef M7_QUICKSORTER_H
#define M7_QUICKSORTER_H

#include <src/defs.h>
#include <functional>
#include "src/core/table/Table.h"

struct QuickSorter {

    defs::inds m_inds;
    typedef std::function<bool(const size_t &, const size_t &)> comp_t;
    comp_t m_comp_fn;

    QuickSorter(comp_t comp_fn);

    const size_t& operator[](const size_t& i) const;

    void preserve_sort(const size_t &hwm);

    void preserve_sort(const TableBase &table);

    bool is_preserve_sorted(const size_t &hwm);

    bool is_preserve_sorted(const TableBase &table);

    void reorder_sort(TableBase &table);

    bool is_reorder_sorted(const TableBase &table);

private:
    void swap(size_t ii1, size_t ii2);

    size_t partition(size_t iilo, size_t iihi);

    void qs(size_t iilo, size_t iihi);

    void swap(size_t ii1, size_t ii2, TableBase &table);

    size_t partition(size_t iilo, size_t iihi, TableBase &table);

    void qs(size_t iilo, size_t iihi, TableBase &table);

};


#endif //M7_QUICKSORTER_H
