//
// Created by rja on 02/11/2020.
//

#ifndef M7_VIEWS_H
#define M7_VIEWS_H

#include "FermionOnvSpecifier.h"
#include "NumericSpecifier.h"
#include "NumericArraySpecifier.h"
#include "BosonOnvSpecifier.h"
#include "FermionBosonOnv.h"

namespace views {
    template <typename T>
    using Number = typename NumericSpecifier<T>::view_t;
    template <typename T, size_t nind>
    using NumberArray = typename NumericArraySpecifier<T, nind>::view_t;
    using Bitset = BitsetSpecifier::view_t;
    using BosonOnv = BosonOnvSpecifier::view_t;

    template <bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, fb_onv::View, FermionOnvSpecifier::view_t>::type;

    using Det = Onv<0>;
    using FbOnv = Onv<1>;
}

#endif //M7_VIEWS_H
