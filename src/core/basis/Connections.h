//
// Created by rja on 04/11/2020.
//

#ifndef M7_CONNECTIONS_H
#define M7_CONNECTIONS_H

#include "DeterminantConnection.h"
#include "BosonOnvConnection.h"
#include "FermionBosonConnection.h"

namespace conn {
    using Determinant = DeterminantConnection;
    using Antisym = AntisymConnection;
    using BosonOnv = BosonOnvConnection;
    using FermionBoson = FermionBosonConnection;
    
    template<bool bosons>
    struct ConfigurationSelector {};
    template<> struct ConfigurationSelector<false> {
        typedef Determinant type;
    };
    template<> struct ConfigurationSelector<true> {
        typedef FermionBosonConnection type;
    };

    using Configuration = ConfigurationSelector<defs::bosons>::type;
}

#endif //M7_CONNECTIONS_H
