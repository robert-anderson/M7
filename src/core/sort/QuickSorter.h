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

    void sort(const size_t &hwm);

    void sort(const TableBase &table);

    bool is_sorted(const size_t &hwm);

    bool is_sorted(const TableBase &table);

private:
    void swap(size_t ii1, size_t ii2);

    size_t partition(size_t iilo, size_t iihi);

    void qs(size_t iilo, size_t iihi);

};

/*
template <typename row_t, typename field_t>
class TableFieldSorter : public QuickSorter {
    static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from ");
    typedef typename viewable_t::cview_t cview_t;

public:
    TableFieldSorter(std::function<cview_t(const size_t&)> getter_fn, bool max=true, bool abs_val=false):
            QuickSorter(sort_utils::make_compare_fn<viewable_t>(getter_fn, max, abs_val)){}

//    TableFieldSorter(const viewable_t& viewable, bool max=true, bool abs_val=false):
//            TableFieldSorter(
//            std::function<const view_t(const size_t &)>(
//                    [&](const size_t &irow) -> const view_t {
//                        return viewable.get_view(irow, 0);
//                    }
//            ),
//            max, abs_val) {}
};
*/


#endif //M7_QUICKSORTER_H
