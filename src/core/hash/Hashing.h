//
// Created by Robert John Anderson on 2020-03-29.
//

#ifndef M7_HASHER_H
#define M7_HASHER_H


#include "src/defs.h"

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

    static defs::hash_t fnv_hash(const char* begin, const size_t& size){
        const auto prime = fnv_prime<defs::hash_t>();
        auto result = fnv_offset_basis<defs::hash_t>();
        for (size_t ibyte=0ul; ibyte<size; ++ibyte) {
            result ^= *(begin+ibyte);
            result*= prime;
        }
        return result;
    }

    static defs::hash_t fnv_hash(defs::hash_t v){
        return fnv_hash((char*)&v, sizeof(defs::hash_t));
    }

    static defs::hash_t in_range(defs::hash_t v, const defs::hash_t& lo, const defs::hash_t& hi){
        ASSERT(hi>lo);
        v = fnv_hash(v);
        v%=hi-lo;
        return v+lo;
    }
};


#endif //M7_HASHER_H
