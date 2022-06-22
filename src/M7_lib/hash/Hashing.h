//
// Created by Robert John Anderson on 2020-03-29.
//

#ifndef M7_HASHER_H
#define M7_HASHER_H

#include <set>

#include <M7_lib/defs.h>

namespace hashing {

    typedef uint64_t hash_t;
    
    template<typename T>
    static T fnv_prime(){
        static_assert(std::is_unsigned<T>::value &&
                      (sizeof(T) == 4 || sizeof(T) == 8), "Invalid hash type");
        switch (sizeof(T)){
            case 4: return 16777619u;
            default: return 1099511628211ul;
        }
    }

    template<typename T>
    static T fnv_offset_basis(){
        static_assert(std::is_unsigned<T>::value &&
                      (sizeof(T) == 4 || sizeof(T) == 8), "Invalid hash type");
        switch (sizeof(T)){
            case 4: return 2166136261u;
            default: return 14695981039346656037ul;
        }
    }

    static hash_t fnv_hash(const char* begin, const size_t& size){
        const auto prime = fnv_prime<hash_t>();
        auto result = fnv_offset_basis<hash_t>();
        for (size_t ibyte=0ul; ibyte<size; ++ibyte) {
            result ^= *(begin+ibyte);
            result*= prime;
        }
        return result;
    }

    static hash_t fnv_hash(hash_t v){
        return fnv_hash(reinterpret_cast<char*>(&v), sizeof(hash_t));
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
    hash_t in_range(hash_t v, hash_t lo, hash_t hi);

    hash_t in_range(const std::vector<hash_t>& v, hash_t lo, hash_t hi);

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
    std::vector<hash_t> in_range(const std::vector<hash_t>& v, size_t ngen, hash_t lo, hash_t hi, bool sorted=false);

    std::vector<hash_t> in_range(hash_t v, size_t ngen, hash_t lo, hash_t hi, bool sorted=false);

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
    std::vector<hash_t> unique_in_range(const std::vector<hash_t>& v, size_t ngen, hash_t lo, hash_t hi, bool sorted=false);

    std::vector<hash_t> unique_in_range(hash_t v, size_t ngen, hash_t  lo, hash_t hi, bool sorted=false);
};


#endif //M7_HASHER_H
