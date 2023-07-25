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

    typedef uint_t digest_t;

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
    digest_t digest_in_range(digest_t v, digest_t lo, digest_t hi);

    digest_t digest_in_range(const v_t<digest_t> &v, digest_t lo, digest_t hi);

    template<typename T=digest_t>
    T in_range(digest_t v, digest_t lo, digest_t hi) {
        static_assert(std::is_integral<T>::value, "digest can only be converted to integer");
        return convert::safe_narrow<T>(digest_in_range(v, lo, hi));
    }

    template<typename T=digest_t>
    T in_range(const v_t<digest_t> &v, digest_t lo, digest_t hi) {
        static_assert(std::is_integral<T>::value, "digest can only be converted to integer");
        return convert::safe_narrow<T>(digest_in_range(v, lo, hi));
    }

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
    template<typename T=digest_t>
    v_t<T> in_range(const v_t<digest_t> &v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false) {
        v_t<T> out;
        out.reserve(ngen);
        auto vatt = v;
        vatt.push_back(0);
        while (out.size() != ngen) {
            auto r = digest_in_range(vatt, lo, hi);
            out.push_back(r);
            ++vatt.back();
        }
        if (sorted) std::sort(out.begin(), out.end());
        return out;
    }

    template<typename T=digest_t>
    v_t<T> in_range(const digest_t &v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false) {
        v_t<digest_t> vec;
        vec.push_back(v);
        return in_range(vec, ngen, lo, hi, sorted);
    }

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
    template<typename T=digest_t>
    v_t<T> unique_in_range(const v_t<digest_t> &v, uint_t ngen,
                                         digest_t lo, digest_t hi, bool sorted = false) {
        REQUIRE_LE(lo + ngen, hi, "number of unique values can't exceed the range of allowed values");
        v_t<T> out;
        out.reserve(ngen);
        std::set<digest_t> set;
        auto vatt = v;
        vatt.push_back(0);
        while (set.size() != ngen) {
            auto r = in_range(vatt, lo, hi);
            while (set.find(r) != set.end()) if (++r == hi) r = lo;
            out.push_back(r);
            set.insert(r);
            ++vatt.back();
        }
        if (sorted) std::sort(out.begin(), out.end());
        return out;
    }

    template<typename T=digest_t>
    v_t<T> unique_in_range(digest_t v, uint_t ngen, digest_t lo, digest_t hi, bool sorted = false) {
        return unique_in_range<T>(v_t<digest_t>{v}, ngen, lo, hi, sorted);
    }
}


#endif //M7_HASH_H
