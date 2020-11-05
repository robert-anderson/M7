//
// Created by rja on 04/11/2020.
//

#include "FermionBosonOnv.h"

fb_onv::View::View(FermionOnvSpecifier::view_t &&fonv, BosonOnvSpecifier::view_t &&bonv) :
        m_fonv(std::move(fonv)), m_bonv(std::move(bonv)) {}

bool fb_onv::View::operator==(const fb_onv::View &other) const {
    return m_fonv == other.m_fonv && m_bonv == other.m_bonv;
}

bool fb_onv::View::operator!=(const fb_onv::View &other) const {
    return !(*this==other);
}

std::string fb_onv::View::to_string() {
    return m_fonv.to_string() + " " + m_bonv.to_string();
}

void fb_onv::View::print() {
    std::cout << to_string() << std::endl;
}
