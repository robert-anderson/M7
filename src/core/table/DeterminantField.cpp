//
// Created by Robert John Anderson on 2020-03-29.
//

#include "DeterminantField.h"

DeterminantElement::DeterminantElement(DeterminantField *field, char *begin) :
    BitsetElement(field, begin) {
}

std::string DeterminantElement::to_string() const {
    std::string result;
    for (size_t ibit = 0; ibit < m_field->nbit(); ++ibit) {
        result += BitsetElement::get(ibit) ? "1" : "0";
    }
    result += "|";
    for (size_t ibit = 0; ibit < m_field->nbit(); ++ibit) {
        result += BitsetElement::get(nsite() + ibit) ? "1" : "0";
    }
    return result;
}

void DeterminantElement::set(const size_t &ispin, const size_t &iorb) {
    BitsetElement::set(defs::pair{ispin ? nsite() : 0, iorb});
}

void DeterminantElement::set(const defs::inds &ispinorbs) {
    for (const auto &ispinorb: ispinorbs) BitsetElement::set(ispinorb);
}

void DeterminantElement::clr(const size_t &ispin, const size_t &iorb) {
    BitsetElement::set(defs::pair{ispin ? nsite() : 0, iorb});
}

bool DeterminantElement::get(const size_t &ispin, const size_t &iorb) const {
    return BitsetElement::get({ispin ? nsite() : 0, iorb});
}

size_t DeterminantElement::nsite() const {
    return dynamic_cast<DeterminantField *>(m_field)->m_nsite;
}

DeterminantElement DeterminantField::element(const size_t &irow, const size_t &isegment, const size_t &ielement) {
    return DeterminantElement(this, element_begin(irow, isegment, ielement));
}

std::string DeterminantField::to_string(size_t irow, size_t isegment, size_t ibegin, size_t iend) {
    std::string result = "";
    for (size_t ielement = 0ul; ielement < m_nelement; ++ielement) {
        result += element(irow, isegment, ielement).to_string() + " ";
    }
    return result;
}

DeterminantField::DeterminantField(Table *table, size_t nelement, size_t nsite) :
    BitsetField(table, nelement, nsite*2), m_nsite(nsite) {}
