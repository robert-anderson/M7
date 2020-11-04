//
// Created by rja on 02/11/2020.
//

#include "DeterminantSpecifier.h"

DeterminantSpecifier::DeterminantSpecifier(const size_t &nsite) : BitsetSpecifier(2 * nsite), m_nsite(nsite) {
    m_data.m_details["type"] = "Determinant";
    m_data.m_details["number of sites"] = std::to_string(nsite);
}

std::string DeterminantSpecifier::element_string(char *ptr) const {
    return View(*this, ptr).to_string();
}

DeterminantSpecifier::View DeterminantSpecifier::operator()(char *ptr) const {
    return View(*this, ptr);
}

DeterminantSpecifier::View::View(const DeterminantSpecifier &field, char *ptr) : BitsetSpecifier::View(field, ptr) {}

std::string DeterminantSpecifier::View::to_string() const {
    std::string res;
    res += "(";
    res.reserve(nbit() * 2 + 3);
    size_t i = 0ul;
    for (; i < nbit() / 2; ++i) res += get(i) ? "1" : "0";
    res += ","; // spin channel delimiter
    for (; i < nbit(); ++i) res += get(i) ? "1" : "0";
    res += ")";
    return res;
}

const size_t &DeterminantSpecifier::View::nsite() const {
    return static_cast<const DeterminantSpecifier&>(m_spec).m_nsite;
}

void DeterminantSpecifier::View::set(const size_t &ispin, const size_t &iorb) {
    BitsetSpecifier::View::set(ispin * nsite() + iorb);
}

void DeterminantSpecifier::View::set(const defs::inds &ispinorbs) {
    for (const auto &ispinorb: ispinorbs) set(ispinorb);
}

void DeterminantSpecifier::View::set(const std::string &s) {
    zero();
    size_t i=0ul;
    for (auto c: s){
        // divider
        if (c!=',') {
            if (c!='0' && c!='1') throw std::runtime_error(
                        R"(Determinant-defining string must contain only "0", "1", or ",")");
            if (c=='1') set(i);
            ++i;
        }
        else {
            if (i!=nsite()) throw std::runtime_error("Divider \",\" is not centralized in Determinant-defining string");
        }
    }
    std::cout << i << " " << 2*nsite() << std::endl;
    if (i<2*nsite()) throw std::runtime_error("Determinant-defining string not long enough");
    if (i>2*nsite()) throw std::runtime_error("Determinant-defining string too long");
}

int DeterminantSpecifier::View::spin() const {
    int spin = 0;
    defs::data_t work;
    for (size_t idataword=0ul; idataword<ndataword(); ++idataword){
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword*defs::nbit_data+bit_utils::next_setbit(work);
            if (ibit < nsite()) ++spin;
            else if (ibit >= 2*nsite()) return spin;
            else --spin;
        }
    }
    return spin;
}

int DeterminantSpecifier::View::nalpha() const {
    // number of electrons occupying spinors in the alpha spin channel
    int nalpha = 0;
    defs::data_t work;
    for (size_t idataword=0ul; idataword<ndataword(); ++idataword){
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword*defs::nbit_data+bit_utils::next_setbit(work);
            if (ibit >= nsite()) return nalpha;
            nalpha++;
        }
    }
    return nalpha;
}
