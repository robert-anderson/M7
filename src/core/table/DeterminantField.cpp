//
// Created by Robert John Anderson on 2020-03-29.
//

#include "DeterminantField.h"

DeterminantElement::DeterminantElement(DeterminantField *field, char *begin) :
    BitsetElement(field, begin) {
}

std::string DeterminantElement::to_string() const {
    std::string result;
    for (size_t ibit = 0; ibit < nsite(); ++ibit) {
        result += BitsetElement::get(ibit) ? "1" : "0";
    }
    result += "|";
    for (size_t ibit = 0; ibit < nsite(); ++ibit) {
        result += BitsetElement::get(nsite() + ibit) ? "1" : "0";
    }
    return result;
}

void DeterminantElement::set(const size_t &ispin, const size_t &iorb) {
    BitsetElement::set(defs::pair{0, ispin ? nsite() + iorb : iorb});
}

void DeterminantElement::set(const defs::inds &ispinorbs) {
    for (const auto &ispinorb: ispinorbs) BitsetElement::set(ispinorb);
}

void DeterminantElement::set(const std::string& s){
    zero();
    size_t i=0ul;
    for (auto c: s){
        // divider
        if (c!='|') {
            if (c!='0' && c!='1') throw std::runtime_error("Determinant-defining string must contain only '0', '1', or '|'");
            if (c=='1') set(i);
            ++i;
        }
        else {
            if (i!=nsite()) throw std::runtime_error("Pipe divider is not centralized in Determinant-defining string");
        }
    }
    std::cout << i << " " << 2*nsite() << std::endl;
    if (i<2*nsite()) throw std::runtime_error("Determinant-defining string not long enough");
    if (i>2*nsite()) throw std::runtime_error("Determinant-defining string too long");
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

DeterminantElement DeterminantField::operator()(const size_t &irow, const size_t &isegment, const size_t &ielement) {
    return DeterminantElement(this, element_begin(irow, isegment, ielement));
}

DeterminantField::DeterminantField(Table *table, size_t nelement, size_t nsite, const std::string &description) :
    BitsetField(table, nelement, nsite * 2, description), m_nsite(nsite) {}

std::string DeterminantField::to_string(size_t irow, size_t isegment, size_t ielement) {
    return (*this)(irow, isegment, ielement).to_string();
}

size_t DeterminantElement::AntiDatawordEnumerator::get_dataword(const size_t &idataword) {
    return m_data.get_antidataword(idataword);
}

size_t DeterminantElement::AntiDatawordEnumerator::get_dataword(const size_t &idataword, const size_t &nbit) {
    return m_data.get_antidataword(idataword, nbit);
}
