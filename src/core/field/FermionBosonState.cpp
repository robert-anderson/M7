//
// Created by rja on 04/11/2020.
//

#include "FermionBosonState.h"

fb_state::View::View(DeterminantSpecifier::view_t &&det, NumericArraySpecifier<uint8_t, 1>::view_t &&perm) :
        m_det(std::move(det)), m_perm(std::move(perm)) {}

bool fb_state::View::operator==(const fb_state::View &other) const {
    return m_det==other.m_det && m_perm==other.m_perm;
}

bool fb_state::View::operator!=(const fb_state::View &other) const {
    return !(*this==other);
}

std::string fb_state::View::to_string() {
    return m_det.to_string() + " " + m_perm.to_string();
}

void fb_state::View::print() {
    std::cout << to_string() << std::endl;
}
