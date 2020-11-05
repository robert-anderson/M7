//
// Created by rja on 04/11/2020.
//

#ifndef M7_CONNECTIONS_H
#define M7_CONNECTIONS_H

#include "FermionOnvConnection.h"
#include "BosonOnvConnection.h"
#include "FermiBosConnection.h"

namespace conn {
    using FermionOnv = FermionOnvConnection;
    using AsFermionOnv = AntisymFermionOnvConnection;
    using BosonOnv = BosonOnvConnection;
    using FermiBosOnv = FermiBosConnection;
    using AsFermiBosOnv = AntisymFermiBosConnection;
    using Onv = std::conditional<defs::bosons, FermiBosOnv, FermionOnv>::type;
    using AsOnv = std::conditional<defs::bosons, AsFermiBosOnv, AsFermionOnv>::type;
}

#endif //M7_CONNECTIONS_H
