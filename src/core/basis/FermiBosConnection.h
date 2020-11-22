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

    explicit FermiBosConnection(const views::FermiBosOnv &in) : FermiBosConnection(in, in) {}

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

    explicit AntisymFermiBosConnection(const views::FermiBosOnv &in) : AntisymFermiBosConnection(in, in) {}

    void connect(const views::FermiBosOnv &in, const views::FermiBosOnv &out) {
        m_aconn.connect(in.m_fonv, out.m_fonv);
        m_bonvconn.connect(in.m_bonv, out.m_bonv);
    }

    void apply(const views::FermiBosOnv &in, views::FermiBosOnv &out) {
        m_aconn.apply(in.m_fonv, out.m_fonv);
        m_bonvconn.apply(in.m_bonv, out.m_bonv);
    }

    void add(const size_t &ann, const size_t &cre) {
        m_aconn.add(ann, cre);
    }

    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2) {
        m_aconn.add(ann1, ann2, cre1, cre2);
    }

    void zero() {
        m_aconn.zero();
        m_bonvconn.zero();
    }

    bool connected() const {
        return m_aconn.connected() && m_bonvconn.nchanged_mode()<=1;
    }
};

#endif //M7_FERMIBOSCONNECTION_H
