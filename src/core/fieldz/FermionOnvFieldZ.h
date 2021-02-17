//
// Created by rja on 17/02/2021.
//

#ifndef M7_FERMIONONVFIELDZ_H
#define M7_FERMIONONVFIELDZ_H

#include "BitsetFieldZ.h"

template<size_t nind_item>
struct FermionOnvFieldBaseZ : BitsetFieldBaseZ<defs::data_t, nind_item, 2> {
    typedef BitsetFieldBaseZ<defs::data_t, nind_item, 2> base_t;
    using base_t::m_item_format;
    using base_t::m_element_format;
    using base_t::m_item_dsize;
    using base_t::get_dataword;
    using base_t::m_nbit;
    using FieldBaseZ::zero;

    FermionOnvFieldBaseZ(size_t nsite): base_t({2, nsite}){}

    const size_t& nsite() const {
        return m_element_format.extent(1);
    }

    size_t nsetbit() const {
        size_t result = 0;
        for (size_t idataword = 0ul; idataword<m_item_dsize; ++idataword){
            result+=bit_utils::nsetbit(get_dataword(idataword));
        }
        return result;
    }

    void set(const std::array<size_t, nind_item>& inds, const size_t &ispinorb) {
        base_t::flat_set(m_item_format.flatten(inds), ispinorb);
    }

    void set(const std::array<size_t, nind_item>& inds, const defs::inds &alpha, const defs::inds &beta) {
        for (auto i: alpha) base_t::set(inds, {0, i});
        for (auto i: beta) base_t::set(inds, {1, i});
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
                if (i != nsite())
                    throw std::runtime_error("Divider \",\" is not centralized in FermionOnv-defining string");
            }
        }
        MPI_REQUIRE(i > m_nbit, "FermionOnv-defining string not long enough");
        MPI_REQUIRE(i < m_nbit, "FermionOnv-defining string too long");
    }

    void clr(const std::array<size_t, nind_item>& inds, const size_t &ispinorb) {
        base_t::flat_clr(m_item_format.flatten(inds), ispinorb);
    }

    bool get(const std::array<size_t, nind_item>& inds, const size_t &ispinorb) {
        return base_t::flat_get(m_item_format.flatten(inds), ispinorb);
    }


    void excite(const std::array<size_t, nind_item>& inds, const size_t &i, const size_t &j) {
        /*
         * single excitation i->j
         */
        const auto iitem = m_item_format.flatten(inds);
        base_t::flat_clr(iitem, i);
        base_t::flat_set(iitem, j);
    }

    void excite(const std::array<size_t, nind_item>& inds, const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        /*
         * double excitation i,j->k,l
         */
        const auto iitem = m_item_format.flatten(inds);
        base_t::flat_clr(iitem, i);
        base_t::flat_clr(iitem, j);
        base_t::flat_set(iitem, k);
        base_t::flat_set(iitem, l);
    }

    int spin(const std::array<size_t, nind_item>& inds) const {
        const auto iitem = m_item_format.flatten(inds);
        int spin = 0;
        defs::data_t work;
        for (size_t idataword = iitem; idataword < m_item_dsize; ++idataword) {
            work = get_dataword(idataword+iitem*m_item_dsize);
            while (work) {
                size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
                if (ibit < nsite()) ++spin;
                else if (ibit >= 2 * nsite()) return spin;
                else --spin;
            }
        }
        return spin;
    }

    int nalpha(const std::array<size_t, nind_item>& inds) const {
        // number of electrons occupying spinors in the alpha spin channel
        const auto iitem = m_item_format.flatten(inds);
        int nalpha = 0;
        defs::data_t work;
        for (size_t idataword = 0ul; idataword < m_item_dsize; ++idataword) {
            work = get_dataword(idataword+iitem*m_item_dsize);
            while (work) {
                size_t ibit = idataword * base_t::nbit_dword() + bit_utils::next_setbit(work);
                if (ibit >= m_item_dsize) return nalpha;
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
        for (; i < nsite(); ++i)
            res += base_t::flat_get(iitem, i) ? "1" : "0";
        res += ","; // spin channel delimiter
        for (; i < m_nbit; ++i)
            res += base_t::flat_get(iitem, i) ? "1" : "0";
        res += ")";
        return res;
    }
};

template<size_t nind_item>
struct FermionOnvField : FermionOnvFieldBaseZ<nind_item> {
    typedef BitsetFieldBaseZ<defs::data_t, nind_item, 2> base_t;
    using base_t::m_item_format;
    using base_t::m_element_format;
    using base_t::m_item_dsize;
    using base_t::get_dataword;
    using base_t::m_nbit;
    using FieldBaseZ::zero;

    FermionOnvFieldBaseZ(size_t nsite) : base_t({2, nsite}) {}
};



#endif //M7_FERMIONONVFIELDZ_H