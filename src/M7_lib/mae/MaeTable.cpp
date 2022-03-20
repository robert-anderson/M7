//
// Created by rja on 11/08/2021.
//

#include "MaeTable.h"

field::MaeInds &MaeRow::key_field() {
    return m_inds;
}

MaeRow::MaeRow(size_t exsig, size_t nvalue) :
        m_inds(this, exsig), m_values(this, {nvalue}, "values"){}

field::SpecMomInds &SpecMomsRow::key_field() {
    return m_inds;
}

SpecMomsRow::SpecMomsRow(size_t exsig, size_t nvalue) :
        m_inds(this, exsig), m_values(this, {nvalue}, "values"){}
