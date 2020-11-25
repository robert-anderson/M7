//
// Created by rja on 04/11/2020.
//

#ifndef M7_FERMIBOSCONNECTION_H
#define M7_FERMIBOSCONNECTION_H

#include "FermionOnvConnection.h"
#include "BosonOnvConnection.h"
#include "src/core/field/Views.h"

//struct FermiBosConnection {
//    FermionOnvConnection m_fonvconn;
//    BosonOnvConnection m_bonvmconn;
//
//    FermiBosConnection(
//            const views::FermiBosOnv &in,
//            const views::FermiBosOnv &out) :
//            m_fonvconn(in.m_fonv, out.m_fonv),
//            m_bonvmconn(in.m_bonv, out.m_bonv) {}
//
//    explicit FermiBosConnection(const views::FermiBosOnv &in) : FermiBosConnection(in, in) {}
//
//    void connect(const views::FermiBosOnv &in, const views::FermiBosOnv &out) {
//        m_fonvconn.connect(in.m_fonv, out.m_fonv);
//        m_bonvmconn.connect(in.m_bonv, out.m_bonv);
//    }
//
//    void apply(const views::FermiBosOnv &in, views::FermiBosOnv &out) {
//        m_fonvconn.apply(in.m_fonv, out.m_fonv);
//        m_bonvmconn.apply(in.m_bonv, out.m_bonv);
//    }
//};


struct AntisymFermiBosConnection : public AntisymFermionOnvConnection {
    BosonOnvConnection m_bonvconn;

    AntisymFermiBosConnection(
            const views::FermiBosOnv &in,
            const views::FermiBosOnv &out);

    explicit AntisymFermiBosConnection(const views::FermiBosOnv &in);

    operator bool() const override;

    void connect(const views::FermiBosOnv &in, const views::FermiBosOnv &out);

    void apply(const views::FermiBosOnv &in, views::FermiBosOnv &out);

    void zero() override;

    bool connected() const;
};

#endif //M7_FERMIBOSCONNECTION_H
