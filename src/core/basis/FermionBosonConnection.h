//
// Created by rja on 04/11/2020.
//

#ifndef M7_FERMIONBOSONCONNECTION_H
#define M7_FERMIONBOSONCONNECTION_H

#include "DeterminantConnection.h"
#include "BosonOnvConnection.h"
#include "src/core/field/Views.h"

class FermionBosonConnection {
    DeterminantConnection m_detconn;
    BosonOnvConnection m_permconn;

    FermionBosonConnection(
            const views::FermionBosonConfiguration& ket,
            const views::FermionBosonConfiguration& bra
            ):
            m_detconn(ket.m_det, bra.m_det),
            m_permconn(ket.m_perm, bra.m_perm){}

    void connect(const views::FermionBosonConfiguration &ket, const views::FermionBosonConfiguration &bra){
        m_detconn.connect(ket.m_det, bra.m_det);
        m_permconn.connect(ket.m_perm, bra.m_perm);
    }
    void apply(const views::FermionBosonConfiguration &ket, views::FermionBosonConfiguration &bra){
        m_detconn.apply(ket.m_det, bra.m_det);
        m_permconn.apply(ket.m_perm, bra.m_perm);
    }

};


#endif //M7_FERMIONBOSONCONNECTION_H
