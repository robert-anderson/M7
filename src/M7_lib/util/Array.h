//
// Created by rja on 15/06/22.
//

#ifndef M7_ARRAY_H
#define M7_ARRAY_H

#include "utils.h"
#include "Convert.h"

namespace utils {
    namespace array {
        template<typename T, size_t nind>
        static std::array<T, nind> filled(const T &v) {
            std::array<T, nind> tmp;
            tmp.fill(v);
            return tmp;
        }

        template<typename T, size_t nind>
        static std::vector<T> to_vector(const std::array<T, nind> &array) {
            std::vector<T> tmp;
            tmp.assign(array.cbegin(), array.cend());
            return tmp;
        }
    }
}

template<typename T, size_t nind>
static std::ostream &operator<<(std::ostream &os, const std::array<T, nind> &a) {
    os << utils::convert::to_string(utils::array::to_vector(a));
    return os;
}

#endif //M7_ARRAY_H
