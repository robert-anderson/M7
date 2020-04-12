//
// Created by Robert John Anderson on 2020-04-12.
//

#include "FciqmcScratch.h"


const size_t FciqmcScratch::nelement_occ = 2;
const size_t FciqmcScratch::nelement_vac = 2;
const size_t FciqmcScratch::nelement_conn = 2;
const size_t FciqmcScratch::nelement_anticonn = 2;

FciqmcScratch::FciqmcScratch(const DeterminantElement &ref) :
    occ(std::make_unique<PrivateStore<OccupiedOrbitals>>(nelement_occ, ref)),
    vac(std::make_unique<PrivateStore<VacantOrbitals>>(nelement_occ, ref)),
    conn(std::make_unique<PrivateStore<Connection>>(nelement_occ, ref, ref)),
    anticonn(std::make_unique<PrivateStore<AntisymConnection>>(nelement_occ, ref, ref)) {}
