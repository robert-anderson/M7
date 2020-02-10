//
// Created by Robert John Anderson on 2020-02-09.
//

#ifndef M7_BITFIELDHASHER_H
#define M7_BITFIELDHASHER_H


#include "BitfieldNew.h"

struct BitfieldHasher {
    size_t operator()(const BitfieldNew &key) const;
};

#endif //M7_BITFIELDHASHER_H
