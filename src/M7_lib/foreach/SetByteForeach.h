//
// Created by rja on 30/07/22.
//

#ifndef M7_SETBYTEFOREACH_H
#define M7_SETBYTEFOREACH_H

#include "M7_lib/util/Bit.h"
#include "M7_lib/util/Functor.h"

namespace setbyte_foreach {

    /**
     * general foreach iterator over single nonzero bytes within a longer word
     * @tparam T
     *  integral host type in which the field data is stored
     * @tparam body_fn_t
     *  callable to execute each time a set byte is reached in the work word
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
                uint_t ibyte = idataword * sizeof(T) + bit::next_setbyte(work);
                fn(ibyte);
            }
        }
    }
}

#endif //M7_SETBYTEFOREACH_H
