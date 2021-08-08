//
// Created by rja on 08/08/2021.
//

#include "FrmOnvField.h"

void FrmOnvField::set_from_string(const std::string &s) {
    zero();
    size_t i = 0ul;
    for (auto c: s) {
        // divider
        if (c != ',') {
            if (c != '0' && c != '1')
                throw std::runtime_error(
                        R"(FermionOnv-defining string must contain only "0", "1", or ",")");
            if (c == '1') set(i);
            ++i;
        } else {
            if (i != m_nsite)
                throw std::runtime_error("Divider \",\" is not centralized in FermionOnv-defining string");
        }
    }
    REQUIRE_GT(i, nbit(), "FermionOnv-defining string not long enough");
    REQUIRE_LT(i, nbit(), "FermionOnv-defining string too long");
}

void FrmOnvField::foreach(const std::function<void(const size_t &)> &body_fn) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            body_fn(ibit);
        }
    }
}

void FrmOnvField::foreach_pair_inner(const size_t &ibit,
                                     const std::function<void(const size_t &, const size_t &)> &body_fn) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t jbit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (jbit==ibit) return;
            body_fn(jbit, ibit);
        }
    }
}

void FrmOnvField::foreach_pair(const std::function<void(const size_t &, const size_t &)> &body_fn) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            foreach_pair_inner(ibit, body_fn);
        }
    }
}

void FrmOnvField::foreach_pair(const std::function<void(const size_t &)> &body_fn_outer,
                               const std::function<void(const size_t &, const size_t &)> &body_fn_inner) const {
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            body_fn_outer(ibit);
            foreach_pair_inner(ibit, body_fn_inner);
        }
    }
}

int FrmOnvField::spin() const {
    int spin = 0;
    size_t work;
    for (size_t idataword = 0; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (ibit < m_nsite) ++spin;
            else if (ibit >= nbit()) return spin;
            else --spin;
        }
    }
    return spin;
}

int FrmOnvField::nalpha() const {
    // number of electrons occupying spinors in the alpha spin channel
    int nalpha = 0;
    size_t work;
    for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
        work = get_dataword(idataword);
        while (work) {
            size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
            if (ibit >= m_nsite) return nalpha;
            nalpha++;
        }
    }
    return nalpha;
}

size_t FrmOnvField::site_nocc(const size_t &isite) const {
    return get({0, isite}) + get({1, isite});
}

std::string FrmOnvField::to_string() const {
    std::string res;
    res += "(";
    res.reserve(nbit() + 3);
    size_t i = 0ul;
    for (; i < m_nsite; ++i)
        res += get(i) ? "1" : "0";
    res += ","; // spin channel delimiter
    for (; i < nbit(); ++i)
        res += get(i) ? "1" : "0";
    res += ")";
    return res;
}

const size_t &FrmOnvField::nsite() const {
    return m_nsite;
}
