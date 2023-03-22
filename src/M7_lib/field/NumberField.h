//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_NUMBERFIELD_H
#define M7_NUMBERFIELD_H

#include "FieldBase.h"
#include "M7_lib/hdf5/Field.h"


struct NumberFieldBase : FieldBase {
    const uint_t m_element_size, m_nelement;
    const bool m_is_complex;

    NumberFieldBase(Row *row, uint_t element_size, uint_t nelement, bool is_complex,
                    const std::type_info &type_info, str_t name = "", bool force_own_words=false);

    NumberFieldBase(const NumberFieldBase& other);

    NumberFieldBase& operator=(const NumberFieldBase& other) {
        FieldBase::operator=(other);
        return *this;
    }
    
    virtual str_t format_string() const = 0;
};


template<typename T, uint_t nind>
struct NdNumberField : NumberFieldBase {
    typedef const uinta_t<nind> &inds_t;
    const NdFormat<nind> m_format;

    const uint_t &nelement() const {
        return m_format.m_nelement;
    }

    NdNumberField(Row *row, NdFormat<nind> format, str_t name = "", bool force_own_words=false) :
            NumberFieldBase(row, sizeof(T), format.m_nelement, dtype::is_complex<T>(),
                    typeid(T), name, force_own_words), m_format(format) {}

    NdNumberField(const NdNumberField &other) : NumberFieldBase(other), m_format(other.m_format){}

    NdNumberField& operator=(const NdNumberField &other) {
        DEBUG_ASSERT_TRUE(m_format==other.m_format, "cannot copy NdNumberField: format mismatch");
        NumberFieldBase::operator=(other);
        return *this;
    }

    NdNumberField &operator=(const T &v) {
        auto begin = begin_as<T>();
        std::fill(begin, begin+nelement(), v);
        return *this;
    }

    NdNumberField &operator=(const v_t<T> &v) {
        DEBUG_ASSERT_EQ(v.size(), nelement(), "Vector size does not match that of numeric field");
        std::copy(v.data(), v.data()+nelement(), begin_as<T>());
        return *this;
    }

    const T* ctbegin() const {
        return cbegin_as<T>();
    }

    T* tbegin() {
        return begin_as<T>();
    }

    const T* ctend() const {
        return cend_as<T>();
    }

    template<typename U>
    bool operator==(const v_t<U> &v) const {
        DEBUG_ASSERT_EQ(v.size(), nelement(), "Vector size does not match that of numeric field");
        for (uint_t i=0ul; i<m_nelement; ++i) if (v[i]!=(*this)[i]) return false;
        return true;
    }

    bool operator==(const NdNumberField &other) const {
        return FieldBase::operator==(other);
    }

    /**
     * float types have two values for zero: 0 and -0. The latter of which does not correspond to cleared bytes, so it
     * does not suffice to use is_clear for numeric types. IEEE 754 demands that -0 and 0 are equal, even if they differ
     * in the signbit
     * @return
     */
    bool is_zero() const {
        for (auto ptr=ctbegin(); ptr!=ctbegin()+m_nelement; ++ptr) {
            if (*ptr != T(0)) return false;
        }
        return true;
    }

    T sum_over(const uintv_t& inds) const {
        T tot{};
        for (const auto& ind: inds) tot+=(*this)[ind];
        return tot;
    }

    bool is_ordered(uint_t ibegin, uint_t iend, bool strict, bool ascending) const {
        if (!nelement()) return true;
        DEBUG_ASSERT_TRUE(ibegin<m_nelement || ibegin==iend, "first element OOB");
        DEBUG_ASSERT_LE(iend, m_nelement, "last element OOB");
        DEBUG_ASSERT_LE(ibegin, iend, "first element must be before last element");
        if (ibegin==iend) return true;
        T last_value = (*this)[ibegin];
        for (uint_t i = ibegin+1; i < iend; ++i) {
            auto this_value = (*this)[i];
            if (strict) {
                if (ascending) {
                    if (this_value <= last_value) return false;
                } else {
                    if (this_value >= last_value) return false;
                }
            } else {
                if (ascending) {
                    if (this_value < last_value) return false;
                } else {
                    if (this_value > last_value) return false;
                }
            }
        }
        return true;
    }

    bool is_ordered(bool strict, bool ascending) const {
        return is_ordered(0, m_nelement, strict, ascending);
    }

private:
    template<typename U>
    void copy_to(v_t<U>& v) const {
        if (v.size() < m_nelement) v.resize(m_nelement);
        // can't copy since the target type doesn't match: must dereference element-wise and attempt to convert
        for (uint_t i=0ul; i<m_nelement; ++i) v[i] = U((*this)[i]);
    }

public:
    void copy_to(v_t<T>& v) const {
        if (v.size() < m_nelement) v.resize(m_nelement);
        std::memcpy(v.data(), cbegin(), m_size);
    }

    template<typename U=T>
    v_t<U> to_vector() const {
        v_t<U> tmp(m_nelement);
        copy_to(tmp);
        return tmp;
    }

    /*
     * math ops
     */

    void add_to(v_t<T>& v) const {
        DEBUG_ASSERT_EQ(v.size(), m_nelement, "incorrect vector size");
        for (uint_t i = 0ul; i < m_nelement; ++i) v[i] += (*this)[i];
    }

    template<typename U>
    NdNumberField &add_scaled(const U &factor, const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] += other[i] * factor;
        return *this;
    }

    NdNumberField &add_abs(const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] += std::abs(other[i]);
        return *this;
    }

    NdNumberField &operator+=(const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] += other[i];
        return *this;
    }

    NdNumberField &operator+=(const T &v) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] += v;
        return *this;
    }

    NdNumberField &sub_abs(const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= std::abs(other[i]);
        return *this;
    }

    template<typename U>
    NdNumberField &sub_scaled(const U &factor, const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= other[i] * factor;
        return *this;
    }

    NdNumberField &operator-=(const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= other[i];
        return *this;
    }

    NdNumberField &operator-=(const T &v) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= v;
        return *this;
    }

    NdNumberField &operator*=(const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] *= other[i];
        return *this;
    }

    NdNumberField &operator*=(const T &v) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] *= v;
        return *this;
    }

    NdNumberField &operator/=(const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] /= other[i];
        return *this;
    }

    NdNumberField &operator/=(const T &v) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] /= v;
        return *this;
    }

    NdNumberField &operator%=(const NdNumberField &other) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] %= other[i];
        return *this;
    }

    NdNumberField &operator%=(const T &v) {
        for (uint_t i = 0ul; i < m_nelement; ++i) (*this)[i] %= v;
        return *this;
    }

    NdNumberField &to_sqrt() {
        for (uint_t ielement = 0ul; ielement < m_nelement; ++ielement)
            (*this)[ielement] = std::sqrt((*this)[ielement]);
        return *this;
    }

    T sum() const {
        T tmp = 0;
        for (uint_t ielement = 0ul; ielement < m_nelement; ++ielement) tmp += (*this)[ielement];
        return tmp;
    }

    /**
     * to prevent overflow, we may elect to accumulate the sum in a larger type
     */
    template<typename U>
    void sum(U& v) const {
        v = 0;
        for (uint_t ielement = 0ul; ielement < m_nelement; ++ielement) v += U((*this)[ielement]);
    }

    uint_t imax() const {
        return std::distance(cbegin_as<T>(), std::max_element(cbegin_as<T>(), cend_as<T>()));
    }
    T max() const {
        return (*this)[imax()];
    }
    uint_t imin() const {
        return std::distance(cbegin_as<T>(), std::min_element(cbegin_as<T>(), cend_as<T>()));
    }
    T min() const {
        return (*this)[imin()];
    }

    /*
     * access methods
     */
    T &operator[](const uint_t &ielement) {
        DEBUG_ASSERT_LT(ielement, m_nelement, "Numeric field access OOB");
        return begin_as<T>()[ielement];
    }

    const T &operator[](const uint_t &ielement) const {
        DEBUG_ASSERT_LT(ielement, m_nelement, "Numeric field access OOB");
        return cbegin_as<T>()[ielement];
    }

    T &operator[](inds_t inds) {
        return begin_as<T>()[m_format.flatten(inds)];
    }

    const T &operator[](inds_t inds) const {
        return cbegin_as<T>()[m_format.flatten(inds)];
    }

    str_t to_string() const override {
        str_t tmp;
        if (nind > 0) tmp += "[";
        for (uint_t ielement = 0ul; ielement < nelement(); ++ielement)
            tmp += convert::to_string((*this)[ielement]) + " ";
        if (nind > 0) tmp += "]";
        return tmp;
    }

    str_t format_string() const override {
        const str_t complex_dim_string = "real/imag (2)";
        if (!nind && m_is_complex) return complex_dim_string;
        return m_format.to_string() + (m_is_complex ? complex_dim_string : "");
    }

    void save_fn(const hdf5::NodeWriter& nw, const str_t& name, bool this_rank, uint_t max_nitem_per_op) const override {
        auto shape = convert::to_vector(m_format.m_shape.data(), nind);
        auto dim_names = convert::to_vector(m_format.m_dim_names.data(), nind);
        hdf5::field::save<T>(*this, nw, name, shape, dim_names, this_rank, max_nitem_per_op);
    }
};

/**
 * no practical use for this class - included only for testing
 */
struct StringField : NdNumberField<char, 1ul> {
    using NdNumberField<char, 1ul>::operator==;
    typedef NdNumberField<char, 1ul> base_t;
    StringField(Row *row, uint_t length, const str_t& name = "", bool force_own_words=false);

    StringField(const StringField& other);

    StringField& operator=(const StringField& other);

    StringField& operator=(const char* str);

    StringField& operator=(const str_t& str);

    bool operator==(const char* str) const;
    bool operator==(const str_t& str) const;
    bool operator!=(const char* str) const;
    bool operator!=(const str_t& str) const;

    str_t to_string() const override;
};

template<typename T>
struct NumberField : NdNumberField<T, 0ul> {
    typedef NdNumberField<T, 0ul> base_t;
    using base_t::operator=;
    using base_t::operator+=;

    NumberField(Row *row, str_t name = "", bool force_own_words=false) :
        base_t(row, {}, name, force_own_words) {}

    NumberField(const NumberField& other): base_t(other){}

    NumberField& operator=(const NumberField& other){
        base_t::operator=(other);
        return *this;
    }

    bool operator==(const T& v) const{
        return (*this)[0] == v;
    }

    operator T&() {
        return *base_t::tbegin();
    }

    operator const T&() const {
        return *base_t::ctbegin();
    }

    arith::comp_t<T> real() const {
        return arith::real(T(*this));
    }

    arith::comp_t<T> imag() const {
        return arith::imag(T(*this));
    }
};

template<typename T>
T operator+(const T& lhs, const NumberField<T>& rhs){
    return lhs+static_cast<const T&>(rhs);
}
template<typename T>
T& operator+=(T& lhs, const NumberField<T>& rhs){
    lhs+=static_cast<const T&>(rhs);
    return lhs;
}

#endif //M7_NUMBERFIELD_H
