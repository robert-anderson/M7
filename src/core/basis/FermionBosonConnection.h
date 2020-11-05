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
            const views::FermiBosOnv& ket,
            const views::FermiBosOnv& bra
            ):
            m_detconn(ket.m_fonv, bra.m_fonv),
            m_permconn(ket.m_bonv, bra.m_bonv){}

    void connect(const views::FermiBosOnv &ket, const views::FermiBosOnv &bra){
        m_detconn.connect(ket.m_fonv, bra.m_fonv);
        m_permconn.connect(ket.m_bonv, bra.m_bonv);
    }
    void apply(const views::FermiBosOnv &ket, views::FermiBosOnv &bra){
        m_detconn.apply(ket.m_fonv, bra.m_fonv);
        m_permconn.apply(ket.m_bonv, bra.m_bonv);
    }

};


#endif //M7_FERMIONBOSONCONNECTION_H
