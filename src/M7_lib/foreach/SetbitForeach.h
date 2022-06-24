//
// Created by Robert J. Anderson on 25/03/2022.
//

#ifndef M7_SETBITFOREACH_H
#define M7_SETBITFOREACH_H

#include "M7_lib/util/Bit.h"

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
    static void single(uint_t dsize, const body_fn_t& fn, const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor::assert_prototype<void(uint_t), body_fn_t>();
        functor::assert_prototype<T(uint_t), get_work_fn_t>();
        T work;
        for (uint_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                uint_t ibit = idataword * (CHAR_BIT*sizeof(T)) + bit::next_setbit(work);
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
    static void pair_inner(uint_t dsize, uint_t ibit, const body_fn_t &fn, const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor::assert_prototype<void(uint_t, uint_t), body_fn_t>();
        functor::assert_prototype<T(uint_t), get_work_fn_t>();
        T work;
        for (uint_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                uint_t jbit = idataword * (CHAR_BIT*sizeof(T)) + bit::next_setbit(work);
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
    static void pair(uint_t dsize, const body_fn_outer_t& fn_outer, const body_fn_inner_t& fn_inner,
                     const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor::assert_prototype<void(uint_t), body_fn_outer_t>();
        functor::assert_prototype<void(uint_t, uint_t), body_fn_inner_t>();
        functor::assert_prototype<T(uint_t), get_work_fn_t>();
        T work;
        for (uint_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                uint_t ibit = idataword * (CHAR_BIT*sizeof(T)) + bit::next_setbit(work);
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
    static void pair(uint_t dsize, const body_fn_inner_t& fn_inner, const get_work_fn_t& get_work_fn) {
        auto fn_outer = [](uint_t){};
        pair<T>(dsize, fn_outer, fn_inner, get_work_fn);
    }


    template<typename T, typename body_fn_3_t, typename get_work_fn_t>
    static void triple_2(uint_t dsize, uint_t ibit, uint_t jbit, const body_fn_3_t& fn_3,
                         const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor::assert_prototype<void(uint_t, uint_t, uint_t), body_fn_3_t>();
        functor::assert_prototype<T(uint_t), get_work_fn_t>();
        DEBUG_ASSERT_LT(jbit, ibit, "triplet should be strictly ordered");

        T work;
        for (uint_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                uint_t kbit = idataword * (CHAR_BIT*sizeof(T)) + bit::next_setbit(work);
                if (kbit==ibit) return;
                if (kbit==jbit) return;
                fn_3(kbit, jbit, ibit);
            }
        }
    }

    template<typename T, typename body_fn_2_t, typename body_fn_3_t, typename get_work_fn_t>
    static void triple_1(uint_t dsize, uint_t ibit, const body_fn_2_t& fn_2, const body_fn_3_t& fn_3,
                       const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor::assert_prototype<void(uint_t, uint_t), body_fn_2_t>();
        functor::assert_prototype<void(uint_t, uint_t, uint_t), body_fn_3_t>();
        functor::assert_prototype<T(uint_t), get_work_fn_t>();

        T work;
        for (uint_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                uint_t jbit = idataword * (CHAR_BIT*sizeof(T)) + bit::next_setbit(work);
                if (jbit==ibit) return;
                fn_2(jbit, ibit);
                triple_2<T>(dsize, ibit, jbit, fn_3, get_work_fn);
            }
        }
    }


    template<typename T, typename body_fn_1_t, typename body_fn_2_t, typename body_fn_3_t, typename get_work_fn_t>
    static void triple(uint_t dsize, const body_fn_1_t& fn_1, const body_fn_2_t& fn_2, const body_fn_3_t& fn_3,
                     const get_work_fn_t& get_work_fn) {
        static_assert(std::is_integral<T>::value, "buffer type must be integral");
        functor::assert_prototype<void(uint_t), body_fn_1_t>();
        functor::assert_prototype<void(uint_t, uint_t), body_fn_2_t>();
        functor::assert_prototype<void(uint_t, uint_t, uint_t), body_fn_3_t>();
        functor::assert_prototype<T(uint_t), get_work_fn_t>();
        T work;
        for (uint_t idataword = 0; idataword < dsize; ++idataword) {
            work = get_work_fn(idataword);
            while (work) {
                uint_t ibit = idataword * (CHAR_BIT*sizeof(T)) + bit::next_setbit(work);
                fn_1(ibit);
                triple_1<T>(dsize, ibit, fn_2, fn_3, get_work_fn);
            }
        }
    }

    /**
     * convenient wrapper for above definition in the commonly-encountered case where there is no need to call a functor
     * on the single set bits or pairs thereof, only the triples. here the outer and middle functors are set to a null
     * lambda.
     */
    template<typename T, typename body_fn_t, typename get_work_fn_t>
    static void triple(uint_t dsize, const body_fn_t& fn, const get_work_fn_t& get_work_fn) {
        auto fn_1 = [](uint_t){};
        auto fn_2 = [](uint_t, uint_t){};
        triple<T>(dsize, fn_1, fn_2, fn, get_work_fn);
    }

}

#endif //M7_SETBITFOREACH_H
