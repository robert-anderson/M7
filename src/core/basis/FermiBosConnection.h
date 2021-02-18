//
// Created by rja on 04/11/2020.
//

#ifndef M7_FERMIBOSCONNECTION_H
#define M7_FERMIBOSCONNECTION_H

#include "FermionOnvConnection.h"
#include "BosonOnvConnection.h"

struct AntisymFermiBosConnection : public AntisymFermionOnvConnection {
    BosonOnvConnection m_bonvconn;

    explicit AntisymFermiBosConnection(size_t nsite);

    AntisymFermiBosConnection(const fieldsz::Onv<1> &in, const fieldsz::Onv<1> &out);

    explicit AntisymFermiBosConnection(const fieldsz::Onv<1> &in);

    operator bool() const override;

    void connect(const fieldsz::Onv<1> &in, const fieldsz::Onv<1> &out);

    void apply(const fieldsz::Onv<1> &in, fieldsz::Onv<1> &out);

    void zero();

    bool connected() const;
};

#endif //M7_FERMIBOSCONNECTION_H
