//
// Created by rja on 29/11/2020.
//

#ifndef M7_EXTREMALVALUES_H
#define M7_EXTREMALVALUES_H


#include "src/defs.h"
#include "src/core/field/TableField.h"
#include <functional>

struct TableX;

class ExtremalValues {
    /*
     * inds will be used as a heap, with popped elements being accumulated in
     * order from the back of the vector
     */
    defs::inds m_inds;
    typedef std::function<bool(const size_t &, const size_t &)> comp_t;

    comp_t m_comp_fn;
    size_t m_nfound;

public:

    ExtremalValues(comp_t comp_fn);

    const size_t& nfound() const;

    const size_t& operator[](const size_t& ifound) const;

    void find(size_t hwm, size_t nfind);

    void find(const TableX &table, size_t nfind);
};

template <typename viewable_t>
class TableExtremalValues : public ExtremalValues {
    static_assert(std::is_base_of<NdFieldGroup<0ul>, viewable_t>::value, "Template arg must be a scalar NdFieldGroup");
    typedef typename viewable_t::const_view_t const_view_t;

public:
    TableExtremalValues(std::function<const_view_t(const size_t&)> getter_fn, bool max=true, bool abs_val=false):
            ExtremalValues(sort_utils::make_compare_fn<viewable_t>(getter_fn, max, abs_val)){}
};

#endif //M7_EXTREMALVALUES_H
