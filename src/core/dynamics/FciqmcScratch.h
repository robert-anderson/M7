//
// Created by Robert John Anderson on 2020-04-12.
//

#ifndef M7_FCIQMCSCRATCH_H
#define M7_FCIQMCSCRATCH_H

#include <src/core/fermion/DecodedDeterminant.h>
#include <src/core/fermion/Connection.h>
#include "src/core/thread/PrivateStore.h"

struct FciqmcScratch {

    /*
     * Thread-private "scratch" storage for objects required in the
     * FCIQMC calculation which have non-trivial constructors
     */
    static const size_t nelement_occ;
    std::unique_ptr<PrivateStore<OccupiedOrbitals>> occ;
    static const size_t nelement_vac;
    std::unique_ptr<PrivateStore<VacantOrbitals>> vac;
    static const size_t nelement_conn;
    std::unique_ptr<PrivateStore<Connection>> conn;
    static const size_t nelement_anticonn;
    std::unique_ptr<PrivateStore<AntisymConnection>> anticonn;

    FciqmcScratch(const DeterminantElement &ref);
};


#endif //M7_FCIQMCSCRATCH_H
