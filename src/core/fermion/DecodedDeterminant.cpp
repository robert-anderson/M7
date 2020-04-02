//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminant.h"

OccupiedOrbitals::OccupiedOrbitals(const DeterminantElement &det_elem) : DecodedDeterminant(det_elem, false){
    update(det_elem);
}

VacantOrbitals::VacantOrbitals(const DeterminantElement &det_elem) : DecodedDeterminant(det_elem, true){
    update(det_elem);
}
