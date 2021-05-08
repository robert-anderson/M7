//
// Created by rja on 07/04/2021.
//

#ifndef M7_FERMIONONVFIELD_H
#define M7_FERMIONONVFIELD_H

#include "BitsetField.h"

struct FermionOnvField : BitsetField<defs::data_t, 2> {
    typedef BitsetField<defs::data_t, 2> base_t;
    using base_t::get;
    using base_t::set;
    using base_t::clr;
    using base_t::put;
    using base_t::operator=;
    using base_t::inds_t;

    const size_t m_nsite;

    FermionOnvField(Row* row, size_t nsite, std::string name=""):
        base_t(row, {{2, nsite}, {"spin channel", "site"}}, name), m_nsite(nsite){}


    FermionOnvField(const FermionOnvField& other):
            FermionOnvField(other.row_of_copy(), other.m_format.extent(1), other.m_name){}

    FermionOnvField& operator=(const FermionOnvField& other) {
        FieldBase::operator=(other);
        return *this;
    }

    void set(const defs::inds& setbits_alpha, const defs::inds& setbits_beta) {
        zero();
        for(auto& i: setbits_alpha) set({0, i});
        for(auto& i: setbits_beta) set({1, i});
    }


    void excite(const size_t &i, const size_t &j) {
        auto* dptr = reinterpret_cast<defs::data_t *>(begin());
        clr(dptr, i);
        set(dptr, j);
    }
    void excite(inds_t ann, inds_t cre){
        auto* dptr = reinterpret_cast<defs::data_t *>(begin());
        clr(dptr, ann);
        set(dptr, cre);
    }

    void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        auto* dptr = reinterpret_cast<defs::data_t *>(begin());
        clr(dptr, i);
        clr(dptr, j);
        set(dptr, k);
        set(dptr, l);
    }
    void excite(inds_t ann1, inds_t ann2, inds_t cre1, inds_t cre2) {
        auto* dptr = reinterpret_cast<defs::data_t *>(begin());
        clr(dptr, ann1);
        clr(dptr, ann1);
        set(dptr, cre1);
        set(dptr, cre2);
    }

    void set_from_string(const std::string &s) {
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
        MPI_REQUIRE(i > nbit(), "FermionOnv-defining string not long enough");
        MPI_REQUIRE(i < nbit(), "FermionOnv-defining string too long");
    }

    int spin() const {
        int spin = 0;
        defs::data_t work;
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

    int nalpha() const {
        // number of electrons occupying spinors in the alpha spin channel
        int nalpha = 0;
        defs::data_t work;
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

    std::string to_string() const override {
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

};


#endif //M7_FERMIONONVFIELD_H
