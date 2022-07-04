//
// Created by Robert J. Anderson on 12/08/2021.
//

#include "SpecMomIndsField.h"

SpecMomIndsField::SpecMomIndsField(Row *row, uint_t exsig, str_t name) :
    base_t(m_left, m_right), m_exsig(exsig),
    m_left(row, exsig, "left " + name), m_right(row, exsig, "right " + name) {}