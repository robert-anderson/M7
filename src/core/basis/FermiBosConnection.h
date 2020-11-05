//
// Created by rja on 04/11/2020.
//

#ifndef M7_FERMIBOSCONNECTION_H
#define M7_FERMIBOSCONNECTION_H

#include "FermionOnvConnection.h"
#include "BosonOnvConnection.h"
#include "src/core/field/Views.h"

struct FermiBosConnection {
    FermionOnvConnection m_fonvconn;
    BosonOnvConnection m_bonvmconn;

    FermiBosConnection(
            const views::FermiBosOnv &ket,
            const views::FermiBosOnv &bra) :
            m_fonvconn(ket.m_fonv, bra.m_fonv),
            m_bonvmconn(ket.m_bonv, bra.m_bonv) {}

    void connect(const views::FermiBosOnv &ket, const views::FermiBosOnv &bra) {
        m_fonvconn.connect(ket.m_fonv, bra.m_fonv);
        m_bonvmconn.connect(ket.m_bonv, bra.m_bonv);
    }

    void apply(const views::FermiBosOnv &ket, views::FermiBosOnv &bra) {
        m_fonvconn.apply(ket.m_fonv, bra.m_fonv);
        m_bonvmconn.apply(ket.m_bonv, bra.m_bonv);
    }
};


struct AntisymFermiBosConnection {
    AntisymFermionOnvConnection m_aconn;
    BosonOnvConnection m_bonvconn;

    AntisymFermiBosConnection(
            const views::FermiBosOnv &ket,
            const views::FermiBosOnv &bra) :
            m_aconn(ket.m_fonv, bra.m_fonv),
            m_bonvconn(ket.m_bonv, bra.m_bonv) {}

    void connect(const views::FermiBosOnv &ket, const views::FermiBosOnv &bra) {
        m_aconn.connect(ket.m_fonv, bra.m_fonv);
        m_bonvconn.connect(ket.m_bonv, bra.m_bonv);
    }

    void apply(const views::FermiBosOnv &ket, views::FermiBosOnv &bra) {
        m_aconn.apply(ket.m_fonv, bra.m_fonv);
        m_bonvconn.apply(ket.m_bonv, bra.m_bonv);
    }
};

#endif //M7_FERMIBOSCONNECTION_H
