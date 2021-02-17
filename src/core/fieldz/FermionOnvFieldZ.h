//
// Created by rja on 17/02/2021.
//

#ifndef M7_FERMIONONVFIELDZ_H
#define M7_FERMIONONVFIELDZ_H

#include "BitsetFieldZ.h"

struct FermionOnvsFieldZ : BitsetFieldBaseZ<defs::data_t> {
    typedef BitsetFieldBaseZ<defs::data_t> base_t;
    using base_t::m_item_dsize;
    using base_t::get_dataword;
    using base_t::m_nbit;
    using FieldBaseZ::zero;

    const size_t m_nsite;

    FermionOnvsFieldZ(size_t nitem, size_t nsite): base_t(nitem, 2*nsite), m_nsite(nsite){}

    void set(const size_t &iitem, const size_t &ibit) {
        base_t::base_set(iitem, ibit);
    }

    void set(const size_t &iitem, const defs::inds &ibits){
        for (auto i: ibits) base_t::base_set(iitem, i);
    }

    void set(const size_t &iitem, const defs::inds &alpha, const defs::inds &beta) {
        for (auto i: alpha) base_t::base_set(iitem, i);
        for (auto i: beta) base_t::base_set(iitem, i+m_nsite);
    }

    void clr(const size_t &iitem, const size_t &ibit) {
        base_t::base_clr(iitem, ibit);
    }

    void put(const size_t &iitem, const size_t &ibit, bool v) {
        base_t::base_put(iitem, ibit, v);
    }

    bool get(const size_t &iitem, const size_t &ibit) const {
        return base_t::base_get(iitem, ibit);
    }

    void set_from_string(const size_t &iitem, const std::string &s) {
        zero();
        size_t i = 0ul;
        for (auto c: s) {
            // divider
            if (c != ',') {
                if (c != '0' && c != '1')
                    throw std::runtime_error(
                            R"(FermionOnv-defining string must contain only "0", "1", or ",")");
                if (c == '1') set(iitem, i);
                ++i;
            } else {
                if (i != m_nsite)
                    throw std::runtime_error("Divider \",\" is not centralized in FermionOnv-defining string");
            }
        }
        MPI_REQUIRE(i > m_nbit, "FermionOnv-defining string not long enough");
        MPI_REQUIRE(i < m_nbit, "FermionOnv-defining string too long");
    }

    void excite(const size_t &iitem, const size_t &i, const size_t &j) {
        /*
         * single excitation i->j
         */
        base_t::base_clr(iitem, i);
        base_t::base_set(iitem, j);
    }

    void excite(const size_t &iitem, const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        /*
         * double excitation i,j->k,l
         */
        base_t::base_clr(iitem, i);
        base_t::base_clr(iitem, j);
        base_t::base_set(iitem, k);
        base_t::base_set(iitem, l);
    }

    int spin(const size_t &iitem) const {
        int spin = 0;
        defs::data_t work;
        for (size_t idataword = iitem; idataword < m_item_dsize; ++idataword) {
            work = get_dataword(idataword+iitem*m_item_dsize);
            while (work) {
                size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
                if (ibit < m_nsite) ++spin;
                else if (ibit >= m_nbit) return spin;
                else --spin;
            }
        }
        return spin;
    }

    int nalpha(const size_t &iitem) const {
        // number of electrons occupying spinors in the alpha spin channel
        int nalpha = 0;
        defs::data_t work;
        for (size_t idataword = 0ul; idataword < m_item_dsize; ++idataword) {
            work = get_dataword(idataword+iitem*m_item_dsize);
            while (work) {
                size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
                if (ibit >= m_nbit) return nalpha;
                nalpha++;
            }
        }
        return nalpha;
    }

    std::string to_string_element(const size_t& iitem) const override {
        std::string res;
        res += "(";
        res.reserve(m_nbit * 2 + 3);
        size_t i = 0ul;
        for (; i < m_nsite; ++i)
            res += get(iitem, i) ? "1" : "0";
        res += ","; // spin channel delimiter
        for (; i < m_nbit; ++i)
            res += get(iitem, i) ? "1" : "0";
        res += ")";
        return res;
    }
};

struct FermionOnvFieldZ : FermionOnvsFieldZ {

    FermionOnvFieldZ(size_t nsite) : FermionOnvsFieldZ(1, nsite){}

    void set(const size_t &ibit) {
        base_t::base_set(ibit);
    }

    void set(const defs::inds &ibits){
        for (auto i: ibits) base_t::base_set(i);
    }

    void set(const defs::inds &alpha, const defs::inds &beta) {
        for (auto i: alpha) base_t::base_set(i);
        for (auto i: beta) base_t::base_set(i+m_nsite);
    }

    void clr(const size_t &ibit) {
        base_t::base_clr(ibit);
    }

    void put(const size_t &ibit, bool v) {
        base_t::base_put(ibit, v);
    }

    bool get(const size_t &ibit) const {
        return base_t::base_get(ibit);
    }

    void excite(const size_t &i, const size_t &j) {
        /*
         * single excitation i->j
         */
        base_t::base_clr(i);
        base_t::base_set(j);
    }

    void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        /*
         * double excitation i,j->k,l
         */
        base_t::base_clr(i);
        base_t::base_clr(j);
        base_t::base_set(k);
        base_t::base_set(l);
    }

    int spin() const {
        return FermionOnvsFieldZ::spin(0);
    }

    int nalpha() const {
        return FermionOnvsFieldZ::nalpha(0);
    }
};



#endif //M7_FERMIONONVFIELDZ_H