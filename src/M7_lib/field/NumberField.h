//
// Created by Robert J. Anderson on 09/02/2021.
//

#ifndef M7_NUMBERFIELD_H
#define M7_NUMBERFIELD_H

#include "FieldBase.h"


struct NumberFieldBase : FieldBase {
    const size_t m_element_size, m_nelement;
    const bool m_is_complex;

    NumberFieldBase(Row *row, size_t element_size, size_t nelement, bool is_complex,
                    const std::type_info &type_info, std::string name = "");

    NumberFieldBase(const NumberFieldBase& other);

    NumberFieldBase& operator=(const NumberFieldBase& other) {
        FieldBase::operator=(other);
        return *this;
    }
    
    virtual std::string format_string() const = 0;
};


template<typename T, size_t nind>
struct NdNumberField : NumberFieldBase {
    typedef const std::array<size_t, nind> &inds_t;
    const NdFormat<nind> m_format;

    const size_t &nelement() const {
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
        return reinterpret_cast<const T*>(dend());
    }

    NdNumberField(Row *row, NdFormat<nind> format, std::string name = "") :
            NumberFieldBase(row, sizeof(T), format.m_nelement,
                            datatype::is_complex<T>(), typeid(T), name), m_format(format) {}

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

    NdNumberField &operator=(const std::vector<T> &v) {
        DEBUG_ASSERT_EQ(v.size(), nelement(), "Vector size does not match that of numeric field");
        std::copy(v.data(), v.data()+nelement(), dbegin());
        return *this;
    }

    template<typename U>
    bool operator==(const std::vector<U> &v) const {
        DEBUG_ASSERT_EQ(v.size(), nelement(), "Vector size does not match that of numeric field");
        for (size_t i=0ul; i<m_nelement; ++i) if (v[i]!=(*this)[i]) return false;
        return true;
    }

    bool operator==(const NdNumberField &other) const {
        return FieldBase::operator==(other);
    }

    T sum_over(const defs::uintv_t& inds) const {
        T tot{};
        for (const auto& ind: inds) tot+=(*this)[ind];
        return tot;
    }

    bool is_ordered(size_t ibegin, size_t iend, bool strict, bool ascending) const {
        if (!nelement()) return true;
        DEBUG_ASSERT_TRUE(ibegin<m_nelement || ibegin==iend, "first element OOB");
        DEBUG_ASSERT_LE(iend, m_nelement, "last element OOB");
        DEBUG_ASSERT_LE(ibegin, iend, "first element must be before last element");
        if (ibegin==iend) return true;
        T last_value = (*this)[ibegin];
        for (size_t i = ibegin+1; i < iend; ++i) {
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
    void copy_to(std::vector<U>& v) const {
        DEBUG_ASSERT_EQ(v.size(), m_nelement, "incorrect vector size");
        // can't copy since the target type doesn't match: must dereference element-wise and attempt to convert
        for (size_t i=0ul; i<m_nelement; ++i) v[i] = U((*this)[i]);
    }

public:
    void copy_to(std::vector<T>& v) const {
        DEBUG_ASSERT_EQ(v.size(), m_nelement, "incorrect vector size");
        std::memcpy(v.data(), begin(), m_size);
    }

    template<typename U=T>
    std::vector<U> to_vector() const {
        std::vector<U> tmp(m_nelement);
        copy_to(tmp);
        return tmp;
    }

    /*
     * math ops
     */

    void add_to(std::vector<T>& v) const {
        DEBUG_ASSERT_EQ(v.size(), m_nelement, "incorrect vector size");
        for (size_t i = 0ul; i < m_nelement; ++i) v[i] += (*this)[i];
    }

    template<typename U>
    NdNumberField &add_scaled(const U &factor, const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] += other[i] * factor;
        return *this;
    }

    NdNumberField &add_abs(const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] += std::abs(other[i]);
        return *this;
    }

    NdNumberField &operator+=(const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] += other[i];
        return *this;
    }

    NdNumberField &operator+=(const T &v) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] += v;
        return *this;
    }

    NdNumberField &sub_abs(const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= std::abs(other[i]);
        return *this;
    }

    template<typename U>
    NdNumberField &sub_scaled(const U &factor, const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= other[i] * factor;
        return *this;
    }

    NdNumberField &operator-=(const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= other[i];
        return *this;
    }

    NdNumberField &operator-=(const T &v) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] -= v;
        return *this;
    }

    NdNumberField &operator*=(const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] *= other[i];
        return *this;
    }

    NdNumberField &operator*=(const T &v) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] *= v;
        return *this;
    }

    NdNumberField &operator/=(const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] /= other[i];
        return *this;
    }

    NdNumberField &operator/=(const T &v) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] /= v;
        return *this;
    }

    NdNumberField &operator%=(const NdNumberField &other) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] %= other[i];
        return *this;
    }

    NdNumberField &operator%=(const T &v) {
        for (size_t i = 0ul; i < m_nelement; ++i) (*this)[i] %= v;
        return *this;
    }

    NdNumberField &to_sqrt() {
        for (size_t ielement = 0ul; ielement < m_nelement; ++ielement)
            (*this)[ielement] = std::sqrt((*this)[ielement]);
        return *this;
    }

    T sum() const {
        T tmp = 0;
        for (size_t ielement = 0ul; ielement < m_nelement; ++ielement) tmp += (*this)[ielement];
        return tmp;
    }

    /*
     * access methods
     */
    T &operator[](const size_t &ielement) {
        DEBUG_ASSERT_LT(ielement, m_nelement, "Numeric field access OOB");
        return dbegin()[ielement];
    }

    const T &operator[](const size_t &ielement) const {
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

    std::string to_string() const override {
        std::string tmp;
        if (nind > 0) tmp += "[";
        for (size_t ielement = 0ul; ielement < nelement(); ++ielement)
            tmp += convert::to_string((*this)[ielement]) + " ";
        if (nind > 0) tmp += "]";
        return tmp;
    }

    std::string format_string() const override {
        const std::string complex_dim_string = "real/imag (2)";
        if (!nind && m_is_complex) return complex_dim_string;
        return m_format.to_string() + (m_is_complex ? complex_dim_string : "");
    }

    /*
     * HDF5 related
     */
    defs::uintv_t h5_shape() const override {
        return {m_format.m_shape.cbegin(), m_format.m_shape.cend()};
    }

    std::vector<std::string> h5_dim_names() const override {
        if (!nind) return {};
        return {m_format.m_dim_names.cbegin(), m_format.m_dim_names.cend()};
    }

    hid_t h5_type() const override {
        return hdf5::type<T>();
    }
};

/**
 * no practical use for this class - included only for testing
 */
struct StringField : NdNumberField<char, 1ul> {
    typedef NdNumberField<char, 1ul> base_t;
    StringField(Row *row, size_t length, std::string name = "");

    StringField(const StringField& other);

    StringField& operator=(const StringField& other);

    StringField& operator=(const char* str);

    StringField& operator=(const std::string& str);

    bool operator==(const char* str) const;
    bool operator==(const std::string& str) const;
    bool operator!=(const char* str) const;
    bool operator!=(const std::string& str) const;

    std::string to_string() const override;
};

template<typename T>
struct NumberField : NdNumberField<T, 0ul> {
    typedef NdNumberField<T, 0ul> base_t;
    using base_t::operator=;

    NumberField(Row *row, std::string name = "") : base_t(row, {}, name) {}

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
