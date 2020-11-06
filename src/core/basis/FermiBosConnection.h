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
            const views::FermiBosOnv &in,
            const views::FermiBosOnv &out) :
            m_fonvconn(in.m_fonv, out.m_fonv),
            m_bonvmconn(in.m_bonv, out.m_bonv) {}

    void connect(const views::FermiBosOnv &in, const views::FermiBosOnv &out) {
        m_fonvconn.connect(in.m_fonv, out.m_fonv);
        m_bonvmconn.connect(in.m_bonv, out.m_bonv);
    }

    void apply(const views::FermiBosOnv &in, views::FermiBosOnv &out) {
        m_fonvconn.apply(in.m_fonv, out.m_fonv);
        m_bonvmconn.apply(in.m_bonv, out.m_bonv);
    }
};


struct AntisymFermiBosConnection {
    AntisymFermionOnvConnection m_aconn;
    BosonOnvConnection m_bonvconn;

    AntisymFermiBosConnection(
            const views::FermiBosOnv &in,
            const views::FermiBosOnv &out) :
            m_aconn(in.m_fonv, out.m_fonv),
            m_bonvconn(in.m_bonv, out.m_bonv) {}

    void connect(const views::FermiBosOnv &in, const views::FermiBosOnv &out) {
        m_aconn.connect(in.m_fonv, out.m_fonv);
        m_bonvconn.connect(in.m_bonv, out.m_bonv);
    }

    void apply(const views::FermiBosOnv &in, views::FermiBosOnv &out) {
        m_aconn.apply(in.m_fonv, out.m_fonv);
        m_bonvconn.apply(in.m_bonv, out.m_bonv);
    }
};

#endif //M7_FERMIBOSCONNECTION_H
