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

public:
    TableExtremalValues(std::function<const view_t(const size_t &)> getter_fn, bool max = true, bool abs_val = false) :
            ExtremalValues(sort_utils::make_compare_fn<viewable_t>(getter_fn, max, abs_val)) {}

    TableExtremalValues(const viewable_t &viewable, bool max = true, bool abs_val = false) :
            TableExtremalValues(
                    std::function<const view_t(const size_t &)>(
                            [&](const size_t &irow) -> const view_t {
                                return viewable.get_view(irow, 0);
                            }
                    ),
                    max, abs_val) {}
};


#endif //M7_TABLEEXTREMALVALUES_H
