//
// Created by rja on 20/06/22.
//

#ifndef M7_TAG_H
#define M7_TAG_H

namespace tag {
    template<typename T>
    struct Type {
    };

    template<uint_t n>
    struct Int {
        static constexpr uint_t value() { return n; }
    };
}


#endif //M7_TAG_H
