//
// Created by rja on 12/06/22.
//

#ifndef M7_UTIL_ND_H
#define M7_UTIL_ND_H

#include "utils.h"

namespace utils {
    namespace nd {
        template<typename T>
        T nelement(const std::vector<T> &v) {
            T out = 1;
            for (const auto &i: v) out *= i;
            return out;
        };
    }
}


#endif //M7_UTIL_ND_H
