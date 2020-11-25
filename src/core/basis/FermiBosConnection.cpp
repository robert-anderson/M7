//
// Created by rja on 04/11/2020.
//

#include "FermiBosConnection.h"

AntisymFermiBosConnection::AntisymFermiBosConnection(const views::FermiBosOnv &in, const views::FermiBosOnv &out) :
        AntisymFermionOnvConnection(in.m_fonv, out.m_fonv),
        m_bonvconn(in.m_bonv, out.m_bonv) {}

AntisymFermiBosConnection::AntisymFermiBosConnection(const views::FermiBosOnv &in) : AntisymFermiBosConnection(in, in) {}

AntisymFermiBosConnection::operator bool() const {
    return nexcit() || m_bonvconn;
}

void AntisymFermiBosConnection::connect(const views::FermiBosOnv &in, const views::FermiBosOnv &out) {
    AntisymFermionOnvConnection::connect(in.m_fonv, out.m_fonv);
    m_bonvconn.connect(in.m_bonv, out.m_bonv);
}

void AntisymFermiBosConnection::apply(const views::FermiBosOnv &in, views::FermiBosOnv &out) {
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
