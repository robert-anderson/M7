//
// Created by rja on 15/06/22.
//

#ifndef M7_ARRAY_H
#define M7_ARRAY_H

#include "Convert.h"

namespace array {
    template<typename T, uint_t nind>
    static std::array<T, nind> filled(const T &v) {
        std::array<T, nind> tmp;
        tmp.fill(v);
        return tmp;
    }

    template<typename T, uint_t nind>
    static v_t<T> to_vector(const std::array<T, nind> &array) {
        v_t<T> tmp;
        tmp.assign(array.cbegin(), array.cend());
        return tmp;
    }
}

template<typename T, uint_t nind>
static std::ostream &operator<<(std::ostream &os, const std::array<T, nind> &a) {
    os << convert::to_string(array::to_vector(a));
    return os;
}

#endif //M7_ARRAY_H
