//
// Created by rja on 29/11/2020.
//

#ifndef M7_PARALLELSORTER_H
#define M7_PARALLELSORTER_H

#include <functional>
#include "src/core/field/Elements.h"

/*
 * need a table of "viewable_t" fields, and a method which returns a viewable_t::view_t
 * instance from a possibly multidimensional source field
 *
 * viewable_t::view_t must have all required comparison operators defined (>, <, =>, <=)
 *
 * e.g. sorting WalkerTables by walker weight, but WalkerTables have multidimensional,
 * m_weight fields (defs::ndim_wf) so a ParallelSort would be required for each wf in that array.
 * In that case the only difference would be in the getter function, which would simply
 * target the multidimensional WalkerTables element in question.
 */

//template <typename viewable_t>
//class ParallelSorter {
//    static_assert(std::is_base_of<NdFieldGroup<0ul>, viewable_t>::value, "Template arg must be a scalar NdFieldGroup");
//    BufferedSingleFieldTable<viewable_t> m_table;
//    typedef typename viewable_t::const_view_t const_view_t;
//    std::function<const_view_t(const size_t&)> m_getter_fn;
//
//    ParallelSorter(std::function<const_view_t(const size_t&)> getter_fn): m_getter_fn(getter_fn){}
//};


#endif //M7_PARALLELSORTER_H
