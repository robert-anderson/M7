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
    typedef std::function<bool(uint_t , uint_t )> index_cmp_fn_t;
    /**
     * templated function definition for scalar value comparisons.
     * returns true if the first argument is "better" than the second
     */
    template<typename T>
    using value_cmp_fn_t = std::function<bool(const T &, const T &)>;

    /**
     * make a function which orders values based on their magnitudes
     * @tparam T
     *  type to sort (signed primitive or complex)
     * @param largest
     *  sense of the sorting inequality operator
     * @return
     *  value comparator
     */
    template<typename T>
    static std::function<bool(const T &, const T &)> make_value_cmp_fn(bool largest, Int<true> /*absval*/) {
        if (largest)
            return [](const T &v, const T &v_cmp) -> bool { return std::abs(v) > std::abs(v_cmp); };
        else
            return [](const T &v, const T &v_cmp) -> bool { return std::abs(v) < std::abs(v_cmp); };
    }

    /**
     * make a function which orders values directly
     * @tparam T
     *  type to sort (signed or unsigned primitive)
     * @param largest
     *  sense of the sorting inequality operator
     * @return
     *  value comparator
     */
    template<typename T>
    static std::function<bool(const T &, const T &)> make_value_cmp_fn(bool largest, Int<false> /*absval*/) {
        if (largest)
            return [](const T &v, const T &v_cmp) -> bool { return v > v_cmp; };
        else
            return [](const T &v, const T &v_cmp) -> bool { return v < v_cmp; };
    }

    /**
     * make a value comparison function for each of the three absval categories
     */
    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool absval, bool largest, Int<Any> /*category*/) {
        if (absval) return make_value_cmp_fn<T>(largest, Int<true>());
        else return make_value_cmp_fn<T>(largest, Int<false>());
    }

    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool /*absval*/, bool largest, Int<Only> /*category*/) {
        return make_value_cmp_fn<T>(largest, Int<true>());
    }

    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool /*absval*/, bool largest, Int<Invalid> /*category*/) {
        return make_value_cmp_fn<T>(largest, Int<false>());
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
    static std::function<bool(const T &, const T &)> make_value_cmp_fn(bool absval, bool largest) {
        // decide category of the given type:
        constexpr uint_t category = dtype::is_complex<T>() ? Only : (std::is_unsigned<T>::value ? Invalid : Any);
        // use this category for tagged dispatch
        return make_value_cmp_fn<T>(absval, largest, Int<category>());
    }

    /**
     * make a function which orders table row indices based on the given comparator
     * @tparam T
     *  numeric type stored in the fields of the compared rows
     * @tparam nind
     *  number of elements in the shape of the compared fields
     * @param field
     *  the field on which the position of row1 relative to row2 is determined
     * @param field_cmp
     *  the field on which the position of row2 relative to row1 is determined
     * @param value_cmp_fn
     *  the function which determines the relative ordering of values
     * @param inds_to_cmp
     *  vector of the flat element indices to be summed together in the comparison of two fields
     * @return
     *  a comparator which determines the relative ordering of two row indices in a Table
     */
    template<typename T, uint_t nind = 0ul>
    static index_cmp_fn_t make_num_field_cmp_fn(
            field::Numbers<T, nind> &field, field::Numbers<T, nind> &field_cmp,
            value_cmp_fn_t<T> value_cmp_fn, uintv_t inds_to_cmp) {
        auto row_ptr = field.m_row;
        auto row_cmp_ptr = field_cmp.m_row;
        REQUIRE_NE(row_ptr, row_cmp_ptr, "specified field pair must correspond to different Rows");
        REQUIRE_EQ(row_ptr->m_table, row_cmp_ptr->m_table, "specified field pair must correspond to the same Table");
        DEBUG_ASSERT_FALSE(inds_to_cmp.empty(), "need at least one numeric field index for comparison");
        DEBUG_ASSERT_TRUE(std::all_of(inds_to_cmp.cbegin(), inds_to_cmp.cend(),
                                      [&field](uint_t i) { return i < field.nelement(); }), "compared field index OOB");
        return [&field, &field_cmp, value_cmp_fn, inds_to_cmp] (uint_t irow, uint_t irow_cmp) -> bool {
            field.m_row->jump(irow);
            field_cmp.m_row->jump(irow_cmp);
            return value_cmp_fn(field.sum_over(inds_to_cmp), field_cmp.sum_over(inds_to_cmp));
        };
    }
};


#endif //M7_COMPARATORS_H
