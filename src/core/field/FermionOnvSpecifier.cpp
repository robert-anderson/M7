//
// Created by rja on 02/11/2020.
//

#include "FermionOnvSpecifier.h"

FermionOnvSpecifier::FermionOnvSpecifier(const size_t &nsite) : BitsetSpecifier(2 * nsite), m_nsite(nsite) {
    m_data.m_details["type"] = "FermionOnv";
    m_data.m_details["number of sites"] = std::to_string(nsite);
}

std::string FermionOnvSpecifier::element_string(char *ptr) const {
    return View(*this, ptr).to_string();
}

FermionOnvSpecifier::View FermionOnvSpecifier::operator()(char *ptr) const {
    return View(*this, ptr);
}

FermionOnvSpecifier::View::View(const FermionOnvSpecifier &spec, char *ptr) : BitsetSpecifier::View(spec, ptr) {}

std::string FermionOnvSpecifier::View::to_string() const {
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

const size_t &FermionOnvSpecifier::View::nsite() const {
    return static_cast<const FermionOnvSpecifier&>(m_spec).m_nsite;
}

void FermionOnvSpecifier::View::set(const size_t &ispin, const size_t &iorb) {
    BitsetSpecifier::View::set(ispin * nsite() + iorb);
}

void FermionOnvSpecifier::View::set(const defs::inds &ispinorbs) {
    for (const auto &ispinorb: ispinorbs) set(ispinorb);
}

void FermionOnvSpecifier::View::set(const defs::inds &alpha, const defs::inds &beta) {
    for (const auto &ispinorb: alpha) set(ispinorb);
    for (const auto &ispinorb: beta) set(ispinorb+nsite());
}

void FermionOnvSpecifier::View::set(const std::string &s) {
    zero();
    size_t i=0ul;
    for (auto c: s){
        // divider
        if (c!=',') {
            if (c!='0' && c!='1') throw std::runtime_error(
                        R"(FermionOnv-defining string must contain only "0", "1", or ",")");
            if (c=='1') set(i);
            ++i;
        }
        else {
            if (i!=nsite()) throw std::runtime_error("Divider \",\" is not centralized in FermionOnv-defining string");
        }
    }
    std::cout << i << " " << 2*nsite() << std::endl;
    if (i<2*nsite()) throw std::runtime_error("FermionOnv-defining string not long enough");
    if (i>2*nsite()) throw std::runtime_error("FermionOnv-defining string too long");
}

int FermionOnvSpecifier::View::spin() const {
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

int FermionOnvSpecifier::View::nalpha() const {
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
