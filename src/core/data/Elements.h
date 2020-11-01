//
// Created by rja on 01/11/2020.
//

#ifndef M7_ELEMENTS_H
#define M7_ELEMENTS_H

#include "BufferedField.h"
#include "BufferedComposite.h"
#include "NumericField.h"
#include "NumericArrayField.h"
#include "Field.h"

namespace elements {
    template<typename T>
    using Number = BufferedComposite<fields::Number<T>>;

    template<typename T, size_t nind>
    using NumberArray = BufferedComposite<fields::NumberArray<T, nind>>;

    using Bitset = BufferedComposite<fields::Bitset>;

    using Determinant = BufferedComposite<fields::Determinant>;

    using BosonOnv = BufferedComposite<fields::BosonOnv>;

    using Configuration = BufferedComposite<fields::Configuration>;
}

#endif //M7_ELEMENTS_H
