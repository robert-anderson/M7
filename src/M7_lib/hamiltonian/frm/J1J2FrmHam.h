//
// Created by anderson on 28/07/2022.
//

#ifndef M7_J1J2FRMHAM_H
#define M7_J1J2FRMHAM_H

#include "SumFrmHam.h"
#include "HeisenbergFrmHam.h"

struct J1J2FrmHam : SumFrmHam<HeisenbergFrmHam, HeisenbergFrmHam> {
    J1J2FrmHam(ham_t j2, const std::shared_ptr<lattice::Lattice>& lattice):
        SumFrmHam<HeisenbergFrmHam, HeisenbergFrmHam>(
                {1.0, lattice},{j2, lattice->m_next_nearest},1.0){}
};


#endif //M7_J1J2FRMHAM_H
