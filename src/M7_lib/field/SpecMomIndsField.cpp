//
// Created by rja on 12/08/2021.
//

#include "SpecMomIndsField.h"

SpecMomIndsField::SpecMomIndsField(Row *row, size_t exsig, std::string name) :
    base_t(m_left, m_right), m_exsig(exsig),
    m_left(row, exsig, "left " + name), m_right(row, exsig, "right " + name) {}