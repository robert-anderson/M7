//
// Created by rja on 04/11/2020.
//

#include "FermiBosConnection.h"

AntisymFermiBosConnection::AntisymFermiBosConnection(size_t nsite) :
        AntisymFermionOnvConnection(nsite), m_bonvconn(nsite) {}

AntisymFermiBosConnection::AntisymFermiBosConnection(const fields::Onv<1> &in, const fields::Onv<1> &out) :
        AntisymFermionOnvConnection(in.m_frm, out.m_frm),
        m_bonvconn(in.m_bos, out.m_bos) {}

AntisymFermiBosConnection::AntisymFermiBosConnection(const fields::Onv<1> &in) : AntisymFermiBosConnection(in, in) {}

AntisymFermiBosConnection::operator bool() const {
    return nexcit() || m_bonvconn;
}

void AntisymFermiBosConnection::connect(const fields::Onv<1> &in, const fields::Onv<1> &out) {
    AntisymFermionOnvConnection::connect(in.m_frm, out.m_frm);
    m_bonvconn.connect(in.m_bos, out.m_bos);
}

void AntisymFermiBosConnection::apply(const fields::Onv<1> &in, fields::Onv<1> &out) {
    AntisymFermionOnvConnection::apply(in.m_frm, out.m_frm);
    m_bonvconn.apply(in.m_bos, out.m_bos);
}

void AntisymFermiBosConnection::apply(const fields::Onv<1> &in) {
    AntisymFermionOnvConnection::apply(in.m_frm);
    m_bonvconn.apply(in.m_bos);
}

void AntisymFermiBosConnection::zero() {
    AntisymFermionOnvConnection::zero();
    m_bonvconn.zero();
}

bool AntisymFermiBosConnection::connected() const {
    return (AntisymFermionOnvConnection::connected() and (m_bonvconn.nchanged_mode()==0))
        or (nexcit()==0 and (m_bonvconn.nchanged_mode()<=1));
}
