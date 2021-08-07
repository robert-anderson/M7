//
// Created by rja on 10/07/2021.
//

#ifndef M7_COMPARATORS_H
#define M7_COMPARATORS_H


#include <src/core/util/utils.h>
#include <src/core/field/Row.h>
#include <src/core/field/Fields.h>

/**
 * type definitions and helper functions for the creation of std::function objects which compare objects in sorting
 * algorithms.
 */
namespace comparators {
    using namespace tags;
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
    static std::function<bool(const T &, const T &)> make_value_cmp_fn(bool largest, Bool<true> absval) {
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
    static std::function<bool(const T &, const T &)> make_value_cmp_fn(bool largest, Bool<false> absval) {
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
    make_value_cmp_fn(bool absval, bool largest, Ind<Any> category) {
        if (absval) return make_value_cmp_fn<T>(largest, Bool<true>());
        else return make_value_cmp_fn<T>(largest, Bool<false>());
    }

    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool absval, bool largest, Ind<Only> category) {
        return make_value_cmp_fn<T>(largest, Bool<true>());
    }

    template<typename T>
    static std::function<bool(const T &, const T &)>
    make_value_cmp_fn(bool absval, bool largest, Ind<Invalid> category) {
        return make_value_cmp_fn<T>(largest, Bool<false>());
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
        constexpr size_t category = consts::is_complex<T>() ? Only : (std::is_unsigned<T>::value ? Invalid : Any);
        // use this category for tagged dispatch
        return make_value_cmp_fn<T>(absval, largest, Ind<category>());
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
     * @param ielement_cmp
     *  the flat element index to be compared between the nind-dimensional fields
     * @return
     *  a comparator which determines the relative ordering of two row indices in a Table
     */
    template<typename row_t, typename T, size_t nind = 0ul>
    static index_cmp_fn_t make_num_field_row_cmp_fn(
            row_t &row1, fields::Numbers<T, nind> &field1,
            row_t &row2, fields::Numbers<T, nind> &field2,
            value_cmp_fn_t<T> value_cmp_fn, size_t ielement_cmp) {
        REQUIRE_TRUE(static_cast<FieldBase &>(field1).belongs_to_row(row1), "specified row-field pair must correspond");
        REQUIRE_TRUE(static_cast<FieldBase &>(field2).belongs_to_row(row2), "specified row-field pair must correspond");
        DEBUG_ASSERT_LE(ielement_cmp, field1.nelement(), "compared field index OOB");
        return [&row1, &field1, &row2, &field2, value_cmp_fn, ielement_cmp]
                (const size_t &irow, const size_t &irow_cmp) {
            static_cast<const Row &>(row1).jump(irow);
            static_cast<const Row &>(row2).jump(irow_cmp);
            return value_cmp_fn(field1[ielement_cmp], field2[ielement_cmp]);
        };
    }
};


#endif //M7_COMPARATORS_H
