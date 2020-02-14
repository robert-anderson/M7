//
// Created by Robert John Anderson on 2020-02-13.
//

#ifndef M7_DETERMINANTHASHER_H
#define M7_DETERMINANTHASHER_H

#include "src/fermion/Determinant.h"

struct DeterminantHasher {
    size_t operator()(const Determinant &key) const;
};

#endif //M7_DETERMINANTHASHER_H
