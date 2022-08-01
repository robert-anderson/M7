//
// Created by Robert John Anderson on 2020-03-29.
//

#ifndef M7_HASH_H
#define M7_HASH_H

#include <M7_lib/defs.h>
#include "Convert.h"

#include "M7_lib/parallel/MPIAssert.h"
#include <set>

namespace hash {

    typedef uint64_t digest_t;

    template<typename T>
    static T fnv_prime() {
        static_assert(std::is_unsigned<T>::value &&
                      (sizeof(T) == 4 || sizeof(T) == 8), "Invalid hash type");
        switch (sizeof(T)) {
            case 4:
                return 16777619u;
            default:
                return 1099511628211ul;
        }
    }

    template<typename T>
    static T fnv_offset_basis() {
        static_assert(std::is_unsigned<T>::value &&
                      (sizeof(T) == 4 || sizeof(T) == 8), "Invalid hash type");
        switch (sizeof(T)) {
            case 4:
                return 2166136261u;
            default:
                return 14695981039346656037ul;
        }
    }

    static digest_t fnv(const buf_t *begin, const uint_t &size) {
        const auto prime = fnv_prime<digest_t>();
        auto result = fnv_offset_basis<digest_t>();
        for (uint_t ibyte = 0ul; ibyte < size; ++ibyte) {
            result ^= *(begin + ibyte);
            result *= prime;
        }
        return result;
    }

    static digest_t fnv(digest_t v) {
        return fnv(reinterpret_cast<buf_t *>(&v), sizeof(digest_t));
    }

    /**
     * deterministically generate arbitrary testing data: NOT a random number generator
     * @param v
     *  value to be hashed
     * @param lo
     *  inclusive minimum returnable value
     * @param hi
     *  exclusive maximum returnable value
     * @return
     *  hash value
     */
    digest_t in_range(digest_t v, digest_t lo, digest_t hi);

    digest_t in_range(const v_t<digest_t> &v, digest_t lo, digest_t hi);

    /**
     * deterministically generate arbitrary testing data: NOT a random number generator
     * @param v
     *  value to be hashed
     * @param ngen
     *  number of values to generate
     * @param lo
     *  inclusive minimum value included in result
     * @param hi
     *  exclusive maximum value included in result
     * @param sorted
     *  true if the result should be put to ascending order before returning
     * @return
     *  arbitrary integers with repetition allowed in the [lo, hi) range
     */
    v_t<digest_t> in_range(const v_t<digest_t> &v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false);

    v_t<digest_t> in_range(digest_t v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false);

    /**
     * deterministically generate arbitrary testing data: NOT a random number generator
     * @param v
     *  value to be hashed
     * @param ngen
     *  number of values to generate
     * @param lo
     *  inclusive minimum value included in result
     * @param hi
     *  exclusive maximum value included in result
     * @param sorted
     *  true if the result should be put to ascending order before returning
     * @return
     *  unrepeated arbitrary integers in the [lo, hi) range
     */
    v_t<digest_t> unique_digest_in_range(const v_t<digest_t> &v, uint_t ngen,
                                         digest_t lo, digest_t hi, bool sorted = false);


    v_t<digest_t> unique_digest_in_range(digest_t v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false);

    template<typename T=digest_t>
    v_t<T> unique_in_range(const v_t<digest_t> &v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false) {
        static_assert(std::is_integral<T>::value, "vector of unique digests can only be converted to integers");
        return convert::safe_narrow<T>(unique_digest_in_range(v, ngen, lo, hi, sorted));
    }

    template<typename T=digest_t>
    v_t<T> unique_in_range(digest_t v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false) {
        return unique_in_range<T>(v_t<digest_t>{v}, ngen, lo, hi, sorted);
    }
}


#endif //M7_HASH_H
