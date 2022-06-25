//
// Created by Robert J. Anderson on 06/04/2021.
//

#ifndef M7_BITSETFIELD_H
#define M7_BITSETFIELD_H

#include <M7_lib/foreach/SetbitForeach.h>
#include "FieldBase.h"

template<typename T, uint_t nind>
struct BitsetField : FieldBase {
    static_assert(std::is_integral<T>::value, "Basis for bitset field must be an integral type");

    typedef const uinta_t<nind> &inds_t;

    struct BitView {
        BitsetField &m_field;
        const uint_t m_ibit = 0;

        BitView(BitsetField &field, uint_t ibit): m_field(field), m_ibit(ibit){}

        operator bool() const {
            return m_field.get(m_ibit);
        }

        BitView &operator=(bool v) {
            m_field.put(m_ibit, v);
        }
    };


    static constexpr uint_t nbit_dword() { return sizeof(T) * CHAR_BIT; }

    const NdFormat<nind> m_format;
    // total number of data words of type T required
    const uint_t m_dsize;
    // number of bits unused in last dataword
    const uint_t m_nbit_in_last_dword;

    using FieldBase::zero;
    using FieldBase::begin;

    BitsetField(Row *row, NdFormat<nind> format, std::string name="") :
        FieldBase(row, integer::divceil(format.m_nelement, nbit_dword()) * sizeof(T), typeid(T), name),
        m_format(format), m_dsize(m_size / sizeof(T)), m_nbit_in_last_dword(nbit() - (m_dsize - 1) * nbit_dword()) {}

    BitsetField(const BitsetField &other) : FieldBase(other),
        m_format(other.m_format), m_dsize(other.m_dsize), m_nbit_in_last_dword(other.m_nbit_in_last_dword){}

    BitsetField &operator=(const BitsetField &other) {
        DEBUG_ASSERT_EQ(m_format.m_nelement, other.m_format.m_nelement,
                        "cannot copy bitset: length mismatch");
        FieldBase::operator=(other);
        return *this;
    }

    BitsetField &operator=(const uintv_t &setbits) {
        // prezero the element
        zero();
        for (const auto &ind: setbits) set(ind);
        return *this;
    }

    const uint_t &nbit() const {
        return m_format.m_nelement;
    }

    T *dbegin() const {
        return reinterpret_cast<T *>(begin());
    }

    BitView operator[](const uint_t &ibit) {
        return {*this, ibit};
    }

    const BitView operator[](const uint_t &ibit) const {
        return {*this, ibit};
    }

    BitView operator[](inds_t inds) {
        return {*this, m_format.flatten(inds)};
    }

    const BitView operator[](inds_t inds) const {
        return {*this, m_format.flatten(inds)};
    }

    bool get(const T *dptr, const uint_t &ibit) const {
        ASSERT(ibit < nbit());
        return bit::get(dptr[ibit / nbit_dword()], ibit % nbit_dword());
    }

    bool get(const uint_t &ibit) const {
        return get(dbegin(), ibit);
    }

    bool get(const T *dptr, inds_t inds) const {
        return get(dptr, m_format.flatten(inds));
    }

    bool get(inds_t inds) const {
        return get(m_format.flatten(inds));
    }

    void set(T *dptr, const uint_t &ibit) {
        ASSERT(ibit < nbit());
        bit::set(dptr[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void set(const uint_t &ibit) {
        set(reinterpret_cast<T *>(begin()), ibit);
    }

    void set(T *dptr, inds_t inds) {
        set(dptr, m_format.flatten(inds));
    }

    void set(inds_t inds) {
        set(m_format.flatten(inds));
    }

    void put_range(const uint_t &ibegin, const uint_t &iend, bool set) {
        T mask;
        auto iword_begin = ibegin/nbit_dword();
        auto iword_end = iend/nbit_dword();
        auto ibitword_begin = ibegin-iword_begin*nbit_dword();
        auto ibitword_end = iend-iword_end*nbit_dword();
        if (iword_begin==iword_end) {
            // begin and end bits are in the same word, so apply mask and return
            bit::apply_mask(dbegin()[iword_begin], ibitword_begin, ibitword_end, set);
            return;
        }
        // begin part: act from first bit to end of its word
        bit::apply_mask(dbegin()[iword_begin], ibitword_begin, nbit_dword(), set);
        // end part: act to the last bit from the beginning of its word
        bit::apply_mask(dbegin()[iword_end], 0, ibitword_end, set);
        // make the all-bits mask
        bit::make_range_mask(mask, 0, nbit_dword());
        for (uint_t iword=iword_begin+1; iword<iword_end; ++iword){
            // apply it to all words between begin and end
            bit::apply_mask(dbegin()[iword], mask, set);
        }
    }

    void set_range(const uint_t &ibegin, const uint_t &iend) {
        put_range(ibegin, iend, true);
    }

    void set() {
        set_range(0, m_format.m_nelement);
    }

    void clr_range(const uint_t &ibegin, const uint_t &iend) {
        put_range(ibegin, iend, false);
    }

    void clr(T *dptr, const uint_t &ibit) {
        ASSERT(ibit < nbit());
        bit::clr(dptr[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void clr(const uint_t &ibit) {
        clr(reinterpret_cast<T *>(begin()), ibit);
    }

    void clr(T *dptr, inds_t inds) {
        clr(dptr, m_format.flatten(inds));
    }

    void clr(inds_t inds) {
        clr(m_format.flatten(inds));
    }

    void put(T *dptr, const uint_t &ibit, bool v) {
        v ? set(dptr, ibit) : clr(dptr, ibit);
    }

    void put(const uint_t &ibit, bool v) {
        v ? set(ibit) : clr(ibit);
    }

    void put(T *dptr, inds_t inds, bool v) {
        put(dptr, m_format.flatten(inds), v);
    }

    void put(inds_t inds, bool v) {
        put(m_format.flatten(inds), v);
    }

    T get_dataword(const uint_t &idataword) const {
        DEBUG_ASSERT_LT(idataword, m_dsize, "dataword index OOB");
        auto *dptr = reinterpret_cast<T *>(begin());
        auto tmp = dptr[idataword];
        if (idataword + 1 == m_dsize) {
            DEBUG_ASSERT_EQ(tmp, bit::truncate(tmp, m_nbit_in_last_dword),
                       "trailing bits were not clear: possible corruption");
            tmp = bit::truncate(tmp, m_nbit_in_last_dword);
        }
        return tmp;
    }

    T get_antidataword(const uint_t &idataword) const {
        DEBUG_ASSERT_LT(idataword, m_dsize, "dataword index OOB");
        auto *dptr = reinterpret_cast<T *>(begin());
        auto tmp = ~dptr[idataword];
        if ((idataword + 1) == m_dsize) {
            tmp = bit::truncate(tmp, m_nbit_in_last_dword);
        }
        return tmp;
    }

    template<typename fn_t>
    void foreach_setbit(const fn_t& fn) const {
        auto get_work_fn = [this](uint_t idataword){return get_dataword(idataword);};
        setbit_foreach::single<T>(m_dsize, fn, get_work_fn);
    }

    template<typename fn_outer_t, typename fn_inner_t>
    void foreach_setbit_pair(const fn_outer_t& fn_outer, const fn_inner_t& fn_inner) const {
        auto get_work_fn = [this](uint_t idataword){return get_dataword(idataword);};
        setbit_foreach::pair<T>(m_dsize, fn_outer, fn_inner, get_work_fn);
    }

    template<typename fn_inner_t>
    void foreach_setbit_pair(const fn_inner_t& fn_inner) const {
        auto get_work_fn = [this](uint_t idataword){return get_dataword(idataword);};
        setbit_foreach::pair<T>(m_dsize, fn_inner, get_work_fn);
    }

    template<typename fn_1_t, typename fn_2_t, typename fn_3_t>
    void foreach_setbit_triple(const fn_1_t& fn_1, const fn_2_t& fn_2, const fn_3_t& fn_3) const {
        auto get_work_fn = [this](uint_t idataword){return get_dataword(idataword);};
        setbit_foreach::triple<T>(m_dsize, fn_1, fn_2, fn_3, get_work_fn);
    }

    template<typename fn_inner_t>
    void foreach_setbit_triple(const fn_inner_t& fn_inner) const {
        auto get_work_fn = [this](uint_t idataword){return get_dataword(idataword);};
        setbit_foreach::triple<T>(m_dsize, fn_inner, get_work_fn);
    }

    uint_t nsetbit() const {
        uint_t result = 0;
        for (uint_t idataword = 0ul; idataword < m_dsize; ++idataword) {
            result += bit::nsetbit(get_dataword(idataword));
        }
        return result;
    }

    std::string to_string() const override {
        std::string res;
        res.reserve(nbit());
        for (uint_t i = 0ul; i < nbit(); ++i)
            res += get(i) ? "1" : "0";
        return res;
    }

    void h5_write_attrs(const hdf5::NodeWriter& node) const override {
        node.write_attr("bitset shape", m_format.shape_vector());
        node.write_attr("bitset dim names", m_format.dim_names_vector());
    }

    uintv_t h5_shape() const override {
        /*
         * bitsets are stored flat
         */
        return {m_dsize};
    }

    std::vector<std::string> h5_dim_names() const override {
        if (!nind) return {};
        return {};
    }

    hid_t h5_type() const override {
        return hdf5::type<T>();
    }
};

template<typename T>
struct BitField : BitsetField<T, 0> {
    typedef BitsetField<T, 0> base_t;

    BitField(Row *row, std::string name="") : base_t(row, {}, name) {}

    BitField &operator=(bool v) {
        base_t::put(0, v);
        return *this;
    }

    operator bool() const {
        return base_t::get(0);
    }
};

#endif //M7_BITSETFIELD_H
