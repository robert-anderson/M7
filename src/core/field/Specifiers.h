//
// Created by rja on 04/11/2020.
//

#ifndef M7_SPECIFIERS_H
#define M7_SPECIFIERS_H

#include "NumericSpecifier.h"
#include "NumericArraySpecifier.h"
#include "BitsetSpecifier.h"
#include "DeterminantSpecifier.h"
#include "BosonOnvSpecifier.h"

namespace specs {
    template<typename T>
    using Number = NumericSpecifier<T>;

    template<typename T, size_t nind>
    using NumberArray = NumericArraySpecifier<T, nind>;

    using Bitset = BitsetSpecifier;
    using Determinant = DeterminantSpecifier;
    using BosonOnv = BosonOnvSpecifier;
}

#endif //M7_SPECIFIERS_H
