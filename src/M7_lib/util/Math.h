//
// Created by rja on 13/06/22.
//

#ifndef M7_UTIL_MATH_H
#define M7_UTIL_MATH_H

#include "utils.h"

namespace utils {
    namespace math {
        /**
         * exact raising to integer power by recursive squaring
         * @tparam T
         *  type of base and result
         * @param x
         *  number being exponentiated
         * @param exp
         *  integer exponent
         * @return
         *  x to the power y
         */
        template<typename T>
        T pow(T x, size_t exp) {
            if (!exp) return 1;
            if (exp == 1) return x;

            auto tmp = pow(x, exp / 2);
            if (exp & 1ul) return x * tmp * tmp;
            return tmp * tmp;
        }

        /**
         * exponentiation in the case that the exponent is a compile time constant integer
         * @tparam exp
         *  integer exponent
         * @tparam T
         *  type of base and result
         * @param x
         *  number being exponentiated
         * @return
         *  x to the power exp
         *
         */
        template<size_t exp, typename T=void>
        static typename std::enable_if<exp == 0, T>::type pow(T x) {
            return 1ul;
        }

        template<size_t exp, typename T=void>
        static typename std::enable_if<exp != 0, T>::type pow(T x) {
            return x * pow<exp - 1, T>(x);
        }
    }
}

#endif //M7_UTIL_MATH_H