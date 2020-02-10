//
// Created by Robert John Anderson on 2020-01-26.
//

#include "DeterminantConnection.h"

DeterminantConnection::DeterminantConnection(const size_t &nbit):
m_removed_inds(defs::inds(nbit)),
m_inserted_inds(defs::inds(nbit)), m_common_inds(defs::inds(nbit)){}

/*
DeterminantConnection::DeterminantConnection(const Determinant &det):
DeterminantConnection::DeterminantConnection(det.m_nbit)
{}*/