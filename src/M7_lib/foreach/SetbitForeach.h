//
// Created by anderson on 25/03/2022.
//

#ifndef M7_SETBITFOREACH_H
#define M7_SETBITFOREACH_H

namespace setbit_foreach {

    /**
     * general foreach iterator over single indices
     * @tparam T
     *  integral host type in which the bit field data is stored
     * @tparam body_fn_t
     *  callable to execute each time a set bit is reached in the work word
     * @tparam get_work_fn_t
     *  callable to retrieve the data words
     * @param dsize
     *  number of datawords to iterate over (number of calls to get_work_fn)
     * @param fn
     *  body_fn_t instance
     * @param get_work_fn
     *  get_work_fn_t instance
     */
    template<typename T, typename body_fn_t, typename get_work_fn_t>
    static void single(size_t dsize, const body_fn_t& fn, const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor_utils::assert_prototype<void(size_t), body_fn_t>();
        functor_utils::assert_prototype<T(size_t), get_work_fn_t>();
        T work;
        for (size_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                size_t ibit = idataword * (CHAR_BIT*sizeof(T)) + bit_utils::next_setbit(work);
                fn(ibit);
            }
        }
    }

    /**
     * general foreach iterator over the inner loop of a pair of set bits: quits when the set bit in the outer loop is
     * reached
     * @tparam T
     *  integral host type in which the bit field data is stored
     * @tparam body_fn_t
     *  callable to execute each time a set bit is reached in the work word
     * @tparam get_work_fn_t
     *  callable to retrieve the data words
     * @param dsize
     *  number of datawords to iterate over (number of calls to get_work_fn)
     * @param ibit
     *  set bit index of the outer loop
     * @param fn
     *  body_fn_t instance
     * @param get_work_fn
     *  get_work_fn_t instance
     */
    template<typename T, typename body_fn_t, typename get_work_fn_t>
    static void pair_inner(size_t dsize, size_t ibit, const body_fn_t &fn, const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor_utils::assert_prototype<void(size_t, size_t), body_fn_t>();
        functor_utils::assert_prototype<T(size_t), get_work_fn_t>();
        T work;
        for (size_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                size_t jbit = idataword * (CHAR_BIT*sizeof(T)) + bit_utils::next_setbit(work);
                if (jbit==ibit) return;
                fn(jbit, ibit);
            }
        }
    }

    /**
     * general foreach iterator over pairs of set bits
     * @tparam T
     *  integral host type in which the bit field data is stored
     * @tparam body_fn_outer_t
     *  callable to execute each time a pair of set bits is reached in a pair of (identical or non-identical) work words
     * @tparam body_fn_inner_t
     *  callable to execute each time a set bit is reached in the work word
     * @tparam get_work_fn_t
     *  callable to retrieve the data words
     * @param dsize
     *  number of datawords to iterate over (number of calls to get_work_fn)
     * @param fn_outer
     *  body_fn_outer_t instance
     * @param fn_inner
     *  body_fn_inner_t instance
     * @param get_work_fn
     *  get_work_fn_t instance
     */
    template<typename T, typename body_fn_outer_t, typename body_fn_inner_t, typename get_work_fn_t>
    static void pair(size_t dsize, const body_fn_outer_t& fn_outer, const body_fn_inner_t& fn_inner,
                     const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor_utils::assert_prototype<void(size_t), body_fn_outer_t>();
        functor_utils::assert_prototype<void(size_t, size_t), body_fn_inner_t>();
        functor_utils::assert_prototype<T(size_t), get_work_fn_t>();
        T work;
        for (size_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                size_t ibit = idataword * (CHAR_BIT*sizeof(T)) + bit_utils::next_setbit(work);
                fn_outer(ibit);
                pair_inner<T>(dsize, ibit, fn_inner, get_work_fn);
            }
        }
    }

    /**
     * convenient wrapper for above definition in the commonly-encountered case where there is no need to call a functor
     * on the single set bits, only the pairs. here the outer functor is set to a null lambda.
     */
    template<typename T, typename body_fn_inner_t, typename get_work_fn_t>
    static void pair(size_t dsize, const body_fn_inner_t& fn_inner, const get_work_fn_t& get_work_fn) {
        auto fn_outer = [](size_t){};
        pair<T>(dsize, fn_outer, fn_inner, get_work_fn);
    }
}

#endif //M7_SETBITFOREACH_H
