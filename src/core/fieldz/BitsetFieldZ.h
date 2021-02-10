//
// Created by rja on 09/02/2021.
//

#ifndef M7_BITSETFIELDZ_H
#define M7_BITSETFIELDZ_H

#include "FieldBaseZ.h"
#include "src/core/util/utils.h"
#include "src/core/parallel/MPIAssert.h"

template<typename T, size_t nind>
struct BitsetFieldBaseZ : FieldBaseZ {
    static_assert(std::is_integral<T>::value, "Basis for bitset field must be an integral type");
    static constexpr size_t nbyte_word() {return sizeof(T);}
    static constexpr size_t nbit_word() {return sizeof(T)*CHAR_BIT;}
    NdFormat<nind> m_format;
    const size_t m_nbit;
    const size_t m_element_dsize;

public:
    BitsetFieldBaseZ(const NdFormat<nind> &format, bool is_key=false) :
            FieldBaseZ(integer_utils::divceil(format.nelement(), nbit_word()) * nbyte_word(),
                       typeid(BitsetFieldBaseZ), is_key),
            m_format(format), m_nbit(m_format.nelement()), m_element_dsize(m_element_size * nbyte_word()) {}

    const size_t &nbit() const {
        return m_format.nelement();
    }

protected:
    bool flat_get(const size_t &ibit) const {
        ASSERT(ibit < nbit());
        return bit_utils::get(((defs::data_t *) raw_view())[ibit / nbit_word()], ibit % nbit_word());
    }

    void flat_set(const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::set(((defs::data_t *) raw_view())[ibit / nbit_word()], ibit % nbit_word());
    }

    void flat_clr(const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::clr(((defs::data_t *) raw_view())[ibit / nbit_word()], ibit % nbit_word());
    }

public:
    void flat_set(const defs::inds &setinds) {
        zero();
        for (auto ibit: setinds) flat_set(ibit);
    }

    void flat_clr(const defs::inds &setinds) {
        for (auto ibit: setinds) flat_clr(ibit);
    }

    template<typename ...Args>
    bool get(Args... inds) const {
        return flat_get(m_format.flatten(inds...));
    }

    template<typename ...Args>
    void set(Args... inds) {
        flat_set(m_format.flatten(inds...));
    }

    template<typename ...Args>
    void clr(Args... inds) {
        flat_clr(m_format.flatten(inds...));
    }

    T get_dataword(const size_t &idataword) const {
        auto tmp = ((T *) raw_view())[idataword];
        if (idataword + 1 == m_element_dsize) {
            auto n = nbit() - nbit_word() * idataword;
            tmp = bit_utils::truncate(tmp, n);
        }
        return tmp;
    }

    T get_antidataword(const size_t &idataword) const {
        auto tmp = ~(((T *) raw_view())[idataword]);
        if (idataword + 1 == m_element_dsize) {
            auto n = nbit() - nbit_word() * idataword;
            tmp = bit_utils::truncate(tmp, n);
        }
        return tmp;
    }

    std::string to_string() const override {
        std::string res;
        res.reserve(nbit());
        for (size_t i = 0ul; i < nbit(); ++i) res += flat_get(i) ? "1" : "0";
        return res;
    }
};


struct FermionOnvFieldZ : BitsetFieldBaseZ<defs::data_t, 2> {
    FermionOnvFieldZ(size_t nsite, bool is_key=false):
            BitsetFieldBaseZ<defs::data_t, 2>({2, nsite}, is_key){}

    const size_t& nsite() const {
        return m_format.extent(1);
    }

    size_t nsetbit() const {
        size_t result = 0;
        for (size_t idataword = 0ul; idataword<m_element_dsize; ++idataword){
            result+=bit_utils::nsetbit(get_dataword(idataword));
        }
        return result;
    }

    void set(const size_t &ispinorb) {
        BitsetFieldBaseZ<defs::data_t, 2>::flat_set(ispinorb);
    }

    void set(const defs::inds &alpha, const defs::inds &beta) {
        for (auto i: alpha) BitsetFieldBaseZ<defs::data_t, 2>::set(0, i);
        for (auto i: beta) BitsetFieldBaseZ<defs::data_t, 2>::set(1, i);
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
                if (i != nsite())
                    throw std::runtime_error("Divider \",\" is not centralized in FermionOnv-defining string");
            }
        }
        MPI_REQUIRE(i > m_nbit, "FermionOnv-defining string not long enough");
        MPI_REQUIRE(i < m_nbit, "FermionOnv-defining string too long");
    }

    void clr(const size_t &ispinorb) {
        BitsetFieldBaseZ<defs::data_t, 2>::flat_clr(ispinorb);
    }

    void clr(const defs::inds &alpha, const defs::inds &beta) {
        for (auto i: alpha) BitsetFieldBaseZ<defs::data_t, 2>::clr(0, i);
        for (auto i: beta) BitsetFieldBaseZ<defs::data_t, 2>::clr(1, i);
    }

    bool get(const size_t &ispinorb) const {
        return BitsetFieldBaseZ<defs::data_t, 2>::flat_get(ispinorb);
    }

    void excite(const size_t &i, const size_t &j) {
        /*
         * single excitation i->j
         */
        clr(i);
        set(j);
    }

    void excite(const size_t &i, const size_t &j, const size_t &k, const size_t &l) {
        /*
         * double excitation i,j->k,l
         */
        clr(i);
        clr(j);
        set(k);
        set(l);
    }

    int spin() const {
        int spin = 0;
        defs::data_t work;
        for (size_t idataword = 0ul; idataword < m_element_dsize; ++idataword) {
            work = get_dataword(idataword);
            while (work) {
                size_t ibit = idataword * nbit_word() + bit_utils::next_setbit(work);
                if (ibit < nsite()) ++spin;
                else if (ibit >= 2 * nsite()) return spin;
                else --spin;
            }
        }
        return spin;
    }

    int nalpha() const {
        // number of electrons occupying spinors in the alpha spin channel
        int nalpha = 0;
        defs::data_t work;
        for (size_t idataword = 0ul; idataword < m_element_dsize; ++idataword) {
            work = get_dataword(idataword);
            while (work) {
                size_t ibit = idataword * nbit_word() + bit_utils::next_setbit(work);
                if (ibit >= m_element_dsize) return nalpha;
                nalpha++;
            }
        }
        return nalpha;
    }

    std::string to_string() const override {
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
};

template<size_t nind>
struct FlagFieldZ : BitsetFieldBaseZ<char, nind> {
    FlagFieldZ(const NdFormat<nind> &format) : BitsetFieldBaseZ<char, nind>(format, false){}
};


#endif //M7_BITSETFIELDZ_H
