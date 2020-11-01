//
// Created by rja on 01/11/2020.
//

#ifndef M7_ELEMENTS_H
#define M7_ELEMENTS_H

#include "BufferedField.h"
#include "BufferedComposite.h"
#include "NumericField.h"
#include "NumericArrayField.h"
#include "Fields.h"

namespace elements {
    template<typename T>
    using Number = BufferedComposite<fields::Number<T>>;

    template<typename T, size_t nind>
    struct NumberArray : BufferedComposite<fields::NumberArray<T, nind>> {
        NumberArray() : BufferedComposite<fields::NumberArray<T, nind>>(
                "Working number array") {}
    };

    struct Bitset : BufferedComposite<fields::Bitset> {
        Bitset(size_t nbit) : BufferedComposite<fields::Bitset>(
                nbit, "Working bitset") {}
    };

    struct Determinant : BufferedComposite<fields::Determinant> {
        Determinant(size_t nsite) : BufferedComposite<fields::Determinant>(
                nsite, "Working determinant") {}
    };

    struct BosonOnv : BufferedComposite<fields::BosonOnv> {
        BosonOnv(size_t nmode) : BufferedComposite<fields::BosonOnv>(
                nmode, "Working Boson ONV") {}
    };

    struct Configuration : BufferedComposite<fields::Configuration> {
        Configuration(size_t nsite, size_t nmode) : BufferedComposite<fields::Configuration>(
                nsite, nmode, "Working configuration") {}
    };
}

#endif //M7_ELEMENTS_H
