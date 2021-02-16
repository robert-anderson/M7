//
// Created by rja on 09/02/2021.
//

#ifndef M7_BITSETFIELDZ_H
#define M7_BITSETFIELDZ_H

#include "FieldBaseZ.h"
#include "src/core/util/utils.h"
#include "src/core/parallel/MPIAssert.h"

template<typename T, size_t nind_item, size_t nind_element>
struct BitsetFieldBaseZ : FullyFormattedFieldBaseZ<T, nind_item, nind_element> {
    static_assert(std::is_integral<T>::value, "Basis for bitset field must be an integral type");
    static constexpr size_t nbyte_word() {return sizeof(T);}
    static constexpr size_t nbit_word() {return sizeof(T)*CHAR_BIT;}

    // total number of bits stored in one element
    const size_t m_nbit;
    // total number of data words of type T required to store each item
    const size_t m_item_dsize;
    // number of bits unused in last dataword of an item
    const size_t m_nbit_spare;

    using FieldBaseZ::m_item_size;
    using FieldBaseZ::m_size;
    using FieldBaseZ::begin;
    using FieldBaseZ::end;
    using FieldBaseZ::zero;
    typedef FullyFormattedFieldBaseZ<T, nind_item, nind_element> base_t;
    using base_t::m_item_format;
    using base_t::m_element_format;
    using base_t::nelement_all;

public:
    BitsetFieldBaseZ(std::array<size_t, nind_element> shape) :
            base_t(integer_utils::divceil(m_nbit, sizeof(T)*CHAR_BIT)*sizeof(T), shape),
            m_nbit(m_element_format.nelement()), m_item_dsize(m_item_size/sizeof(T)),
            m_nbit_spare(m_item_dsize*nbit_word() - m_nbit){
    }

    bool flat_get(const size_t &iitem, const size_t &ibit) const {
        ASSERT(ibit < m_nbit);
        return bit_utils::get(begin(iitem)[iitem * ibit / nbit_word()], ibit % nbit_word());
    }

    void flat_set(const size_t &iitem, const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::set(begin(iitem)[ibit / nbit_word()], ibit % nbit_word());
    }

    void flat_clr(const size_t &iitem, const size_t &ibit) {
        ASSERT(ibit < m_nbit);
        bit_utils::clr(begin(iitem)[ibit / nbit_word()], ibit % nbit_word());
    }

    void flat_set(const size_t &iitem, const defs::inds &setinds) {
        zero();
        for (auto ibit: setinds) flat_set(iitem, ibit);
    }

    T get_dataword(const size_t &idataword) const {
        auto tmp = ((T *) begin())[idataword];
        if ((idataword + 1)%m_item_dsize==0) {
            tmp = bit_utils::truncate(tmp, m_nbit_spare);
        }
        return tmp;
    }

    T get_antidataword(const size_t &idataword) const {
        auto tmp = ~(((T *) begin())[idataword]);
        if ((idataword + 1)%m_item_dsize==0) {
            tmp = bit_utils::truncate(tmp, m_nbit_spare);
        }
        return tmp;
    }

    std::string to_string_element(const size_t& iitem) const override {
        std::string res;
        res.reserve(m_nbit);
        for (size_t i = 0ul; i < m_nbit; ++i)
            res += flat_get(iitem, i) ? "1" : "0";
        return res;
    }
};

template<size_t nind_item>
struct FermionOnvFieldZ : BitsetFieldBaseZ<defs::data_t, nind_item, 2> {
    typedef BitsetFieldBaseZ<defs::data_t, nind_item, 2> base_t;
    using base_t::m_item_format;
    using base_t::m_element_format;
    using base_t::m_item_dsize;
    using base_t::get_dataword;
    using base_t::m_nbit;
    using FieldBaseZ::zero;

    FermionOnvFieldZ(size_t nsite): base_t({2, nsite}){}

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
                size_t ibit = idataword * base_t::nbit_word() + bit_utils::next_setbit(work);
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
                size_t ibit = idataword * base_t::nbit_word() + bit_utils::next_setbit(work);
                if (ibit >= m_item_dsize) return nalpha;
                nalpha++;
            }
        }
        return nalpha;
    }

    std::string to_string_element(const size_t& iitem) const override {
        std::string res;
        res += "(";
        res.reserve(m_nbit() * 2 + 3);
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

/*
 * two specializations for pretty access syntax, need dummy template arg to make
 * use of partial template specialization
 */

template<typename dummy_t, size_t nind_element>
struct FlagFieldZ : BitsetFieldBaseZ<uint8_t, 0ul, nind_element> {
    typedef BitsetFieldBaseZ<uint8_t, 0ul, nind_element> base_t;
    using base_t::m_element_format;
    using base_t::nbit_word;
    using base_t::begin;
    FlagFieldZ(std::array<size_t, nind_element> shape) : base_t(shape){}

    bool get(const std::array<size_t, nind_element>& einds) const {
        auto const& ibit = m_element_format.flatten(einds);
        return bit_utils::get(((uint8_t*)begin())[ibit / nbit_word()], ibit % nbit_word());
    }

    void set(const std::array<size_t, nind_element>& einds) const {
        auto const& ibit = m_element_format.flatten(einds);
        bit_utils::set(((uint8_t*)begin())[ibit / nbit_word()], ibit % nbit_word());
    }

    void clr(const std::array<size_t, nind_element>& einds) const {
        auto const& ibit = m_element_format.flatten(einds);
        bit_utils::clr(((uint8_t*)begin())[ibit / nbit_word()], ibit % nbit_word());
    }
};

template<typename dummy_t>
struct FlagFieldZ<dummy_t, 0ul> : BitsetFieldBaseZ<uint8_t, 0ul, 0ul> {
    typedef BitsetFieldBaseZ<uint8_t, 0ul, 0ul> base_t;
    using base_t::m_element_format;
    using base_t::nbit_word;
    using base_t::begin;
    FlagFieldZ() : base_t({}){}

    bool get() const {
        return bit_utils::get(((uint8_t*)begin())[0], 0);
    }

    void set() const {
        bit_utils::set(((uint8_t*)begin())[0], 0);
    }

    void clr() const {
        bit_utils::clr(((uint8_t*)begin())[0], 0);
    }
};



#endif //M7_BITSETFIELDZ_H
