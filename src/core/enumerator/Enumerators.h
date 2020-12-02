//
// Created by rja on 05/11/2020.
//

#ifndef M7_ENUMERATORS_H
#define M7_ENUMERATORS_H

#include "FermiBosOnvEnumerator.h"

namespace enums {
    using Combination = CombinationEnumerator;
    using Product = ProductEnumerator;
    using FermionOnv = FermionOnvEnumerator;
    using BosonOnv = BosonOnvEnumerator;
    using FermiBosOnv = FermiBosOnvEnumerator;
    using Onv = std::conditional<defs::enable_bosons, FermiBosOnv, FermionOnv>::type;
}

#endif //M7_ENUMERATORS_H
