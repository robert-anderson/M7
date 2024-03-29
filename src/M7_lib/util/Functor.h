//
// Created by rja on 15/06/22.
//

#ifndef M7_FUNCTOR_H
#define M7_FUNCTOR_H

#include <functional>
#include "M7_lib/defs.h"


namespace functor {
    template<typename signature_t, typename fn_t>
    void assert_prototype() {
        static_assert(std::is_convertible<fn_t, std::function<signature_t>>::value,
                      "function does not conform to required prototype");
    }

    template<typename signature_t, typename fn_t>
    void assert_prototype(const fn_t &) {
        assert_prototype<signature_t, fn_t>();
    }
}

#endif //M7_FUNCTOR_H
