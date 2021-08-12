//
// Created by rja on 12/08/2021.
//

#include "SpecMomIndsField.h"

SpecMomIndsField::SpecMomIndsField(Row *row, size_t exsig, std::string name) :
    base_t(row, {nullptr, exsig, "left " + name}, {nullptr, exsig, "right " + name}),
    m_exsig(exsig), m_name(name), m_left(get<0>()), m_right(get<1>()) {}

SpecMomIndsField::SpecMomIndsField(const SpecMomIndsField &other) :
    SpecMomIndsField(other.row_of_copy(), other.m_exsig, other.m_name) {}
