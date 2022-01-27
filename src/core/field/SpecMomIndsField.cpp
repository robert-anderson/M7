//
// Created by rja on 12/08/2021.
//

#include "SpecMomIndsField.h"

SpecMomIndsField::SpecMomIndsField(Row *row, size_t exsig, std::string name) :
    base_t(row, {nullptr, exsig, "left " + name}, {nullptr, exsig, "right " + name}),
    m_exsig(exsig), m_name(name), m_left(get<0>()), m_right(get<1>()) {}

SpecMomIndsField::SpecMomIndsField(const SpecMomIndsField &other) :
    SpecMomIndsField(other.row_of_copy(), other.m_exsig, other.m_name) {}

#if 0
SpecMomIndsField::SpecMomIndsField(Row *row, size_t exsig, std::string name) :
    base_t(m_left, m_right), m_exsig(exsig),
    m_left(row, exsig, "left " + name), m_right(row, exsig, "right " + name){}

SpecMomIndsField::SpecMomIndsField(const SpecMomIndsField &other) :
        base_t(m_left, m_right), m_exsig(other.m_exsig), m_left(other.m_left), m_right(other.m_right){}

SpecMomIndsField &SpecMomIndsField::operator=(const SpecMomIndsField &other) {
    DEBUG_ASSERT_EQ(m_exsig, other.m_exsig, "excitation signature mismatch");
    m_left = other.m_left;
    m_right = other.m_right;
    return *this;
}

#endif