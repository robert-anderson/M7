//
// Created by rja on 04/11/2020.
//

#ifndef M7_CONNECTIONS_H
#define M7_CONNECTIONS_H

#include "FermionOnvConnection.h"
#include "BosonOnvConnection.h"
#include "FermiBosConnection.h"

namespace conn {

    template <bool enable_bosons=defs::enable_bosons>
    using Basic = typename std::conditional<enable_bosons, AntisymFermiBosConnection, FermionOnvConnection>::type;

    template <bool enable_bosons=defs::enable_bosons>
    using Antisym = typename std::conditional<enable_bosons, AntisymFermiBosConnection, AntisymFermionOnvConnection>::type;

    using Boson = BosonOnvConnection;
//    using FermionOnv = FermionOnvConnection;
//    using AsFermionOnv = AntisymFermionOnvConnection;
//
//    //using FermiBosOnv = FermiBosConnection;
//    using AsFermiBosOnv = AntisymFermiBosConnection;
//    //using Onv = std::conditional<defs::enable_bosons, FermiBosOnv, FermionOnv>::type;
//    using AsOnv = std::conditional<defs::enable_bosons, AsFermiBosOnv, AsFermionOnv>::type;
}

#endif //M7_CONNECTIONS_H
