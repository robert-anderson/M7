//
// Created by rja on 07/12/22.
//

#ifndef M7_VECTOR_H
#define M7_VECTOR_H

#include "M7_lib/defs.h"

/**
 * functions acting on std::vectors
 */
namespace vector {

    template<typename T>
    v_t<T> range(T begin, T end, T interval=T(1)){
        v_t<T> tmp;
        tmp.reserve(uint_t((end - begin) / interval));
        for (T v = begin; v < end; v+=interval) tmp.emplace_back(v);
        return tmp;
    }

    template<typename T>
    v_t<T> inserted(const v_t<T>& vec, const v_t<T>& insertion, uint_t pos){
        auto tmp = vec;
        tmp.insert(tmp.begin()+pos, insertion.cbegin(), insertion.cend());
        return tmp;
    }
    template<typename T>
    v_t<T> prepended(const v_t<T>& vec, const v_t<T>& insertion){
        return inserted(vec, insertion, 0ul);
    }
    template<typename T>
    v_t<T> appended(const v_t<T>& vec, const v_t<T>& insertion){
        return inserted(vec, insertion, vec.size());
    }

    template<typename T, typename U>
    v_t<T> inserted(const v_t<T>& vec, const U& insertion, uint_t pos){
        v_t<T> tmp;
        tmp.push_back(T(insertion));
        return inserted(vec, tmp, pos);
    }
    template<typename T, typename U>
    v_t<T> prepended(const v_t<T>& vec, const U& insertion){
        return inserted(vec, insertion, 0ul);
    }
    template<typename T, typename U>
    v_t<T> appended(const v_t<T>& vec, const U& insertion){
        return inserted(vec, insertion, vec.size());
    }
}


#endif //M7_VECTOR_H
