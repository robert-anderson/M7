//
// Created by rja on 08/12/2020.
//

#ifndef M7_TABLEEXTREMALVALUES_H
#define M7_TABLEEXTREMALVALUES_H

#include "ExtremalValues.h"

template<typename viewable_t>
class TableExtremalValues : public ExtremalValues {

    static_assert(std::is_base_of<NdFieldGroup<0ul>, viewable_t>::value, "Template arg must be a scalar NdFieldGroup");
    typedef typename viewable_t::view_t view_t;
    typedef typename viewable_t::cview_t cview_t;

public:
    TableExtremalValues(std::function<cview_t(const size_t &)> getter_fn, bool max = true, bool abs_val = false) :
            ExtremalValues(sort_utils::make_compare_fn<viewable_t>(getter_fn, max, abs_val)) {}

    TableExtremalValues(const viewable_t &viewable, bool max = true, bool abs_val = false) :
            TableExtremalValues(
                    [&](const size_t &irow) -> cview_t {
                        return viewable.get_view(irow, 0);
                    }, max, abs_val) {}
};


#endif //M7_TABLEEXTREMALVALUES_H
