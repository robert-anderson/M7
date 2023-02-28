//
// Created by Robert J. Anderson on 11/08/2021.
//

#include "MaeTable.h"

field::MaeInds &MaeRow::key_field() {
    return m_inds;
}

MaeRow::MaeRow(OpSig exsig, uint_t nvalue) :
        m_inds(this, exsig), m_values(this, {nvalue}, "values"){}

field::SpecMomInds &SpecMomsRow::key_field() {
    return m_inds;
}

SpecMomsRow::SpecMomsRow(uint_t nvalue) : m_inds(this), m_values(this, {nvalue}, "values"){}
