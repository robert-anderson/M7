//
// Created by Robert John Anderson on 2020-04-12.
//

#include "FciqmcScratch.h"

FciqmcScratch::FciqmcScratch(const DeterminantElement &ref) :
    occ(new PrivateStore<OccupiedOrbitals>(ref)),
    vac(new PrivateStore<VacantOrbitals>(ref)),
    conn(new PrivateStore<Connection>(ref)),
    anticonn(new PrivateStore<AntisymConnection>(ref))
{}
