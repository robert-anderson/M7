//
// Created by Robert John Anderson on 2020-03-29.
//

#ifndef M7_HASHER_H
#define M7_HASHER_H


#include "src/core/util/defs.h"

namespace hashing {

    template<typename hash_T>
    static hash_T fnv_prime(){
        static_assert(std::is_unsigned<hash_T>::value &&
                      (sizeof(hash_T)==4 || sizeof(hash_T)==8), "Invalid hash type");
        switch (sizeof(hash_T)){
            case 4: return 16777619u;
            default: return 1099511628211ul;
        }
    }

    template<typename hash_T>
    static hash_T fnv_offset_basis(){
        static_assert(std::is_unsigned<hash_T>::value &&
                      (sizeof(hash_T)==4 || sizeof(hash_T)==8), "Invalid hash type");
        switch (sizeof(hash_T)){
            case 4: return 2166136261u;
            default: return 14695981039346656037ul;
        }
    }

    static defs::hash_t fnv_hash(char* begin, const size_t& size){
        const auto prime = fnv_prime<defs::hash_t>();
        auto result = fnv_offset_basis<defs::hash_t>();
        for (size_t ibyte=0ul; ibyte<size; ++ibyte) {
            result ^= *(begin+ibyte);
            result*= prime;
        }
        return result;
    }
};


#endif //M7_HASHER_H
