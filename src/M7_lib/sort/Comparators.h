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
    typedef std::function<bool(const size_t &, const size_t &)> index_cmp_fn_t;
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
    static std::function<bool(const T &, const T &)> make_value_cmp_fn(bool largest, Int<true> absval) {
        if (largest)
            return [](const T &v, const T &v_cmp) { return std::abs(v) > std::abs(v_cmp); };
        else
            return [](const T &v, const T &v_cmp) { return std::abs(v) < std::abs(v_cmp); };
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
    static std::function<bool(const T &, const T &)> make_value_cmp_fn(bool largest, Int<false> absval) {
        if (largest)
            return [](const T &v, const T &v_cmp) { return v > v_cmp; };
        else
            return [](const T &v, const T &v_cmp) { return v < v_cmp; };
    }

    /**
     * make a value comparison function for each of the three absval categories
     */
    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool absval, bool largest, Int<Any> category) {
        if (absval) return make_value_cmp_fn<T>(largest, Int<true>());
        else return make_value_cmp_fn<T>(largest, Int<false>());
    }

    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool absval, bool largest, Int<Only> category) {
        return make_value_cmp_fn<T>(largest, Int<true>());
    }

    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool absval, bool largest, Int<Invalid> category) {
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
        constexpr size_t category = datatype::is_complex<T>() ? Only : (std::is_unsigned<T>::value ? Invalid : Any);
        // use this category for tagged dispatch
        return make_value_cmp_fn<T>(absval, largest, Int<category>());
    }

    /**
     * make a function which orders table row indices based on the given comparator
     * @tparam row_t
     *  Row-derived type
     * @tparam T
     *  numeric type stored in the fields of the compared rows
     * @tparam nind
     *  number of elements in the shape of the compared fields
     * @param row1
     *  a row reference which sorting algorithms can regard as a cursor to move around between data elements to compare
     * @param field1
     *  the field on which the position of row1 relative to row2 is determined
     * @param row2
     *  another row reference to enable comparisons
     * @param field2
     *  the field on which the position of row2 relative to row1 is determined
     * @param value_cmp_fn
     *  the function which determines the relative ordering of values
     * @param inds_to_cmp
     *  vector of the flat element indices to be summed together in the comparison of two fields
     * @return
     *  a comparator which determines the relative ordering of two row indices in a Table
     */
    template<typename row_t, typename T, size_t nind = 0ul>
    static index_cmp_fn_t make_num_field_row_cmp_fn(
            row_t &row1, field::Numbers<T, nind> &field1,
            row_t &row2, field::Numbers<T, nind> &field2,
            value_cmp_fn_t<T> value_cmp_fn, defs::ivec_t inds_to_cmp) {
        REQUIRE_TRUE(static_cast<FieldBase &>(field1).belongs_to_row(row1), "specified row-field pair must correspond");
        REQUIRE_TRUE(static_cast<FieldBase &>(field2).belongs_to_row(row2), "specified row-field pair must correspond");
        DEBUG_ASSERT_FALSE(inds_to_cmp.empty(), "need at least one numeric field index for comparison");
        DEBUG_ASSERT_TRUE(std::all_of(inds_to_cmp.cbegin(), inds_to_cmp.cend(),
                                      [&field1](size_t i) { return i < field1.nelement(); }), "compared field index OOB");
        return [&row1, &field1, &row2, &field2, value_cmp_fn, inds_to_cmp]
                (const size_t &irow, const size_t &irow_cmp) {
            static_cast<const Row &>(row1).jump(irow);
            static_cast<const Row &>(row2).jump(irow_cmp);
            return value_cmp_fn(field1.sum_over(inds_to_cmp), field2.sum_over(inds_to_cmp));
        };
    }
};


#endif //M7_COMPARATORS_H
