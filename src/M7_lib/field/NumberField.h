//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_NUMBERFIELD_H
#define M7_NUMBERFIELD_H

#include "FieldBase.h"


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

    T* dbegin() {
        return reinterpret_cast<T*>(begin());
    }

    const T* dbegin() const {
        return reinterpret_cast<const T*>(begin());
    }

    T* dend() {
        return reinterpret_cast<T*>(end());
    }

    const T* dend() const {
        return reinterpret_cast<const T*>(end());
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
        std::fill(dbegin(), dend(), v);
        return *this;
    }

    NdNumberField &operator=(const v_t<T> &v) {
        DEBUG_ASSERT_EQ(v.size(), nelement(), "Vector size does not match that of numeric field");
        std::copy(v.data(), v.data()+nelement(), dbegin());
        return *this;
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
        DEBUG_ASSERT_EQ(v.size(), m_nelement, "incorrect vector size");
        // can't copy since the target type doesn't match: must dereference element-wise and attempt to convert
        for (uint_t i=0ul; i<m_nelement; ++i) v[i] = U((*this)[i]);
    }

public:
    void copy_to(v_t<T>& v) const {
        DEBUG_ASSERT_EQ(v.size(), m_nelement, "incorrect vector size");
        std::memcpy(v.data(), begin(), m_size);
    }

    template<typename U=T>
    v_t<U> to_vector() const {
        v_t<U> tmp(m_nelement);
        copy_to(tmp);
        return tmp;
    }

    void to_buffer(v_t<T>& buf, uint_t irow_begin, uint_t nitem_max, std::set<uint_t> irows_empty) const {
        DEBUG_ASSERT_LT(irow_begin, irow_end, "invalid row range");
        buf.resize(nitem_max);
        const auto nitem = FieldBase::to_buffer(static_cast<buf_t*>(buf.data()), irow_begin, nitem_max, irows_empty);
        buf.resize(m_nelement*nitem);
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
        return std::distance(dbegin(), std::max_element(dbegin(), dend()));
    }
    T max() const {
        return (*this)[imax()];
    }
    uint_t imin() const {
        return std::distance(dbegin(), std::min_element(dbegin(), dend()));
    }
    T min() const {
        return (*this)[imin()];
    }

    /*
     * access methods
     */
    T &operator[](const uint_t &ielement) {
        DEBUG_ASSERT_LT(ielement, m_nelement, "Numeric field access OOB");
        return dbegin()[ielement];
    }

    const T &operator[](const uint_t &ielement) const {
        DEBUG_ASSERT_LT(ielement, m_nelement, "Numeric field access OOB");
        return dbegin()[ielement];
    }

    T &operator[](inds_t inds) {
        return dbegin()[m_format.flatten(inds)];
    }

    const T &operator[](inds_t inds) const {
        return dbegin()[m_format.flatten(inds)];
    }

    template<typename U=T>
    operator const typename std::enable_if<!nind, U>::type &() const {
        return (*this)[0];
    }

    template<typename U=T>
    operator typename std::enable_if<!nind, U>::type &() {
        return (*this)[0];
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

    /*
     * HDF5 related
     */
    uintv_t h5_shape() const override {
        return {m_format.m_shape.cbegin(), m_format.m_shape.cend()};
    }

    strv_t h5_dim_names() const override {
        if (!nind) return {};
        return {m_format.m_dim_names.cbegin(), m_format.m_dim_names.cend()};
    }

    hdf5::Type h5_type() const override {
        return {static_cast<T*>(nullptr)};
    }
};

/**
 * no practical use for this class - included only for testing
 */
struct StringField : NdNumberField<char, 1ul> {
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
        return *base_t::dbegin();
    }

    operator const T&() const {
        return *base_t::dbegin();
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
