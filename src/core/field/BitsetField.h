//
// Created by rja on 06/04/2021.
//

#ifndef M7_BITSETFIELD_H
#define M7_BITSETFIELD_H

#include "FieldBase.h"

template<typename T, size_t nind>
struct BitsetField : FieldBase {
    static_assert(std::is_integral<T>::value, "Basis for bitset field must be an integral type");

    typedef const std::array<size_t, nind> &inds_t;

    struct BitView {
        BitsetField &m_field;
        const size_t m_ibit = 0;

        BitView(BitsetField &field, size_t ibit): m_field(field), m_ibit(ibit){}

        operator bool() const {
            return m_field.get(m_ibit);
        }

        BitView &operator=(bool v) {
            m_field.put(m_ibit, v);
        }
    };


    static constexpr size_t nbit_dword() { return sizeof(T) * CHAR_BIT; }

    const NdFormat<nind> m_format;
    // total number of data words of type T required
    const size_t m_dsize;
    // number of bits unused in last dataword
    const size_t m_nbit_in_last_dword;

    using FieldBase::zero;
    using FieldBase::begin;

    BitsetField(Row *row, NdFormat<nind> format, std::string name="") :
        FieldBase(row, integer_utils::divceil(format.m_nelement, nbit_dword()) * sizeof(T), typeid(T), name),
        m_format(format), m_dsize(m_size / sizeof(T)), m_nbit_in_last_dword(nbit() - (m_dsize - 1) * nbit_dword()) {}

    BitsetField(const BitsetField &other) : FieldBase(other),
        m_format(other.m_format), m_dsize(other.m_dsize), m_nbit_in_last_dword(other.m_nbit_in_last_dword){}

    BitsetField &operator=(const BitsetField &other) {
        DEBUG_ASSERT_EQ(m_format.m_nelement, other.m_format.m_nelement,
                        "cannot copy bitset: length mismatch");
        FieldBase::operator=(other);
        return *this;
    }

    BitsetField &operator=(const defs::inds &setbits) {
        // prezero the element
        zero();
        for (const auto &ind: setbits) set(ind);
        return *this;
    }

    BitsetField(BitsetField &&other) : FieldBase(std::move(other)),
       m_format(other.m_format), m_dsize(other.m_dsize), m_nbit_in_last_dword(other.m_nbit_in_last_dword){}

    BitsetField &operator=(BitsetField &&other) {
        *this = other;
        return *this;
    }

    const size_t &nbit() const {
        return m_format.m_nelement;
    }

    T *dbegin() const {
        return reinterpret_cast<T *>(begin());
    }

    BitView operator[](const size_t &ibit) {
        return {*this, ibit};
    }

    const BitView operator[](const size_t &ibit) const {
        return {*this, ibit};
    }

    BitView operator[](inds_t inds) {
        return {*this, m_format.flatten(inds)};
    }

    const BitView operator[](inds_t inds) const {
        return {*this, m_format.flatten(inds)};
    }

    bool get(const T *dptr, const size_t &ibit) const {
        ASSERT(ibit < nbit());
        return bit_utils::get(dptr[ibit / nbit_dword()], ibit % nbit_dword());
    }

    bool get(const size_t &ibit) const {
        return get(dbegin(), ibit);
    }

    bool get(const T *dptr, inds_t inds) const {
        return get(dptr, m_format.flatten(inds));
    }

    bool get(inds_t inds) const {
        return get(m_format.flatten(inds));
    }

    void set(T *dptr, const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::set(dptr[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void set(const size_t &ibit) {
        set(reinterpret_cast<T *>(begin()), ibit);
    }

    void set(T *dptr, inds_t inds) {
        set(dptr, m_format.flatten(inds));
    }

    void set(inds_t inds) {
        set(m_format.flatten(inds));
    }


    void clr(T *dptr, const size_t &ibit) {
        ASSERT(ibit < nbit());
        bit_utils::clr(dptr[ibit / nbit_dword()], ibit % nbit_dword());
    }

    void clr(const size_t &ibit) {
        clr(reinterpret_cast<T *>(begin()), ibit);
    }

    void clr(T *dptr, inds_t inds) {
        clr(dptr, m_format.flatten(inds));
    }

    void clr(inds_t inds) {
        clr(m_format.flatten(inds));
    }

    void put(T *dptr, const size_t &ibit, bool v) {
        v ? set(dptr, ibit) : clr(dptr, ibit);
    }

    void put(const size_t &ibit, bool v) {
        v ? set(ibit) : clr(ibit);
    }

    void put(T *dptr, inds_t inds, bool v) {
        put(dptr, m_format.flatten(inds), v);
    }

    void put(inds_t inds, bool v) {
        put(m_format.flatten(inds), v);
    }

    T get_dataword(const size_t &idataword) const {
        T *dptr = reinterpret_cast<T *>(begin());
        auto tmp = dptr[idataword];
        if (idataword + 1 == m_dsize) {
            DEBUG_ASSERT_EQ(tmp, bit_utils::truncate(tmp, m_nbit_in_last_dword),
                       "trailing bits were not clear: possible corruption");
            tmp = bit_utils::truncate(tmp, m_nbit_in_last_dword);
        }
        return tmp;
    }

    T get_antidataword(const size_t &idataword) const {
        T *dptr = reinterpret_cast<T *>(begin());
        auto tmp = ~dptr[idataword];
        if ((idataword + 1) == m_dsize) {
            tmp = bit_utils::truncate(tmp, m_nbit_in_last_dword);
        }
        return tmp;
    }

    size_t nsetbit() const {
        size_t result = 0;
        for (size_t idataword = 0ul; idataword < m_dsize; ++idataword) {
            result += bit_utils::nsetbit(get_dataword(idataword));
        }
        return result;
    }

    std::string to_string() const override {
        std::string res;
        res.reserve(nbit());
        for (size_t i = 0ul; i < nbit(); ++i)
            res += get(i) ? "1" : "0";
        return res;
    }

    void h5_write_attrs(hid_t parent_handle) override {
        hdf5::AttributeWriterBase::write(parent_handle, "bitset shape", m_format.shape_vector());
        hdf5::AttributeWriterBase::write(parent_handle, "bitset dim names", m_format.dim_names_vector());
    }

    defs::inds h5_shape() const override {
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
