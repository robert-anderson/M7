//
// Created by Robert J. Anderson on 12/08/2021.
//

#include "SpecMomIndsField.h"

SpecMomIndsField::SpecMomIndsField(Row *row, str_t name) :
    base_t(m_left, m_right), m_left(row, "left " + name), m_right(row, "right " + name) {}