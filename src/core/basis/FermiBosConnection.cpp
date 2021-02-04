//
// Created by rja on 04/11/2020.
//

#include "FermiBosConnection.h"

AntisymFermiBosConnection::AntisymFermiBosConnection(size_t nsite) :
        AntisymFermionOnvConnection(nsite), m_bonvconn(nsite) {}

AntisymFermiBosConnection::AntisymFermiBosConnection(const views::Onv<1> &in, const views::Onv<1> &out) :
        AntisymFermionOnvConnection(in.m_fonv, out.m_fonv),
        m_bonvconn(in.m_bonv, out.m_bonv) {}

AntisymFermiBosConnection::AntisymFermiBosConnection(const views::Onv<1> &in) : AntisymFermiBosConnection(in, in) {}

AntisymFermiBosConnection::operator bool() const {
    return nexcit() || m_bonvconn;
}

void AntisymFermiBosConnection::connect(const views::Onv<1> &in, const views::Onv<1> &out) {
    AntisymFermionOnvConnection::connect(in.m_fonv, out.m_fonv);
    m_bonvconn.connect(in.m_bonv, out.m_bonv);
}

void AntisymFermiBosConnection::apply(const views::Onv<1> &in, views::Onv<1> &out) {
    AntisymFermionOnvConnection::apply(in.m_fonv, out.m_fonv);
    m_bonvconn.apply(in.m_bonv, out.m_bonv);
}

void AntisymFermiBosConnection::zero() {
    AntisymFermionOnvConnection::zero();
    m_bonvconn.zero();
}

bool AntisymFermiBosConnection::connected() const {
    return (AntisymFermionOnvConnection::connected() and (m_bonvconn.nchanged_mode()==0))
        or (nexcit()==0 and (m_bonvconn.nchanged_mode()<=1));
}
