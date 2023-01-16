//
// Created by Robert J. Anderson on 10/07/2021.
//

#ifndef M7_COMPARATORS_H
#define M7_COMPARATORS_H


#include <M7_lib/field/Row.h>
#include <M7_lib/field/Fields.h>

/**
 * type definitions and helper functions for the creation of std::function objects which compare objects in sorting
 * algorithms.
 */
namespace comparators {
    using namespace tag;
    /**
     * the optionality of ordering based on absolute values if not universal,
     * for some types we have the option,
     * for others we have no option but to sort by the unaltered values (unsigned),
     * for others we have no option but to sort by the absolute values (complex)
     * these three categories are encapsulated in these compile time indices
     */
    enum AbsValCategory {
        Any, Invalid, Only
    };
    /**
     * all row index comparators take a pair of indices and return their relative order as boolean
     */
    typedef std::function<bool(uint_t, uint_t)> index_cmp_fn_t;
    /**
     * templated function pointer definition for scalar value comparisons.
     * returns true if the first argument is superior to the second
     */
    template<typename T>
    using value_cmp_fn_t = bool(*)(const T &, const T &);

    /**
     * basic value comparators
     */
    template<typename T> bool abs_gt(const T &v1, const T &v2) { return std::abs(v1) > std::abs(v2); };
    template<typename T> bool abs_lt(const T &v1, const T &v2) { return std::abs(v1) < std::abs(v2); };
    template<typename T> bool gt(const T &v1, const T &v2) { return v1 > v2; };
    template<typename T> bool lt(const T &v1, const T &v2) { return v1 < v2; };

    /**
     * return a function which orders values based on their magnitudes
     * @tparam T
     *  type to sort (signed primitive or complex)
     * @param largest
     *  sense of the sorting inequality operator
     * @return
     *  value comparator
     */
    template<typename T>
    value_cmp_fn_t<T> get_value_cmp_fn(bool largest, Int<true> /*absval*/) {
        return largest ? abs_gt<T> : abs_lt<T>;
    }

    /**
     * return a function which orders values directly
     * @tparam T
     *  type to sort (signed or unsigned primitive)
     * @param largest
     *  sense of the sorting inequality operator
     * @return
     *  value comparator
     */
    template<typename T>
    value_cmp_fn_t<T> get_value_cmp_fn(bool largest, Int<false> /*absval*/) {
        return largest ? gt<T> : lt<T>;
    }

    /**
     * make a value comparison function for each of the three absval categories
     */
    template<typename T>
    value_cmp_fn_t<T> get_value_cmp_fn(bool absval, bool largest, Int<Any> /*category*/) {
        if (absval) return get_value_cmp_fn<T>(largest, Int<true>());
        else return get_value_cmp_fn<T>(largest, Int<false>());
    }

    template<typename T>
    value_cmp_fn_t<T> get_value_cmp_fn(bool /*absval*/, bool largest, Int<Only> /*category*/) {
        return get_value_cmp_fn<T>(largest, Int<true>());
    }

    template<typename T>
    value_cmp_fn_t<T> get_value_cmp_fn(bool /*absval*/, bool largest, Int<Invalid> /*category*/) {
        return get_value_cmp_fn<T>(largest, Int<false>());
    }

    /**
     * make a function which orders values
     * @tparam T
     *  type to sort
     * @param absval
     *  if true, sort the values based on their magnitudes if valid
     * @param largest
     *  sense of the sorting inequality operator
     * @return
     *  value comparator
     */
    template<typename T>
    value_cmp_fn_t<T> get_value_cmp_fn(bool absval, bool largest) {
        auto tmp = gt(T{}, T{}) || lt(T{}, T{}) || abs_gt(T{}, T{}) || abs_lt(T{}, T{});
        (void) tmp;
        // decide category of the given type:
        constexpr uint_t category = dtype::is_complex<T>() ? Only : (std::is_unsigned<T>::value ? Invalid : Any);
        // use this category for tagged dispatch
        return get_value_cmp_fn<T>(absval, largest, Int<category>());
    }

    /**
     * make a function which orders table row indices based on the given comparator
     * @tparam T
     *  numeric type stored in the fields of the compared rows
     * @tparam nind
     *  number of elements in the shape of the compared fields
     * @param field1
     *  the field whose value is superior if value_cmp_fn is true
     * @param field2
     *  the field whose value is superior if value_cmp_fn is false
     * @param value_cmp_fn
     *  the function which determines the relative ordering of values
     * @param inds_to_cmp
     *  vector of the flat element indices to be summed together in the comparison of two fields
     * @return
     *  a comparator which determines the relative ordering of two row indices in a Table
     */
    template<typename T, uint_t nind = 0ul>
    index_cmp_fn_t make_num_field_cmp_fn(
            field::Numbers<T, nind> &field1, field::Numbers<T, nind> &field2,
            value_cmp_fn_t<T> value_cmp_fn, uintv_t inds_to_cmp) {
        auto row1 = field1.m_row;
        auto row2 = field2.m_row;
        REQUIRE_NE(row1, row2, "specified field pair must correspond to different Rows");
        REQUIRE_EQ(row1->m_table, row2->m_table, "specified field pair must correspond to the same Table");
        DEBUG_ASSERT_FALSE(inds_to_cmp.empty(), "need at least one numeric field index for comparison");
        DEBUG_ASSERT_TRUE(std::all_of(inds_to_cmp.cbegin(), inds_to_cmp.cend(),
              [&field1](uint_t i) { return i < field1.nelement(); }), "compared field index OOB");
        return [&field1, &field2, value_cmp_fn, inds_to_cmp] (uint_t irow1, uint_t irow2) -> bool {
            field1.m_row->jump(irow1);
            field2.m_row->jump(irow2);
            return value_cmp_fn(field1.sum_over(inds_to_cmp), field2.sum_over(inds_to_cmp));
        };
    }
};


#endif //M7_COMPARATORS_H
