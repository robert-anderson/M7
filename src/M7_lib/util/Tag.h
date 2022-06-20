//
// Created by rja on 20/06/22.
//

#ifndef M7_TAG_H
#define M7_TAG_H

#include "utils.h"

namespace utils {
    namespace tag {
        template<typename T>
        struct Type {
        };

        template<size_t n>
        struct Int {
            static constexpr size_t value() { return n; }
        };
    }
}


#endif //M7_TAG_H
