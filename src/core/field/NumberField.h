//
// Created by rja on 09/02/2021.
//

#ifndef M7_NUMBERFIELD_H
#define M7_NUMBERFIELD_H

#include "FieldBase.h"


struct NumberFieldBase : FieldBase {
    const size_t m_element_size, m_nelement;
    const bool m_is_complex;

    NumberFieldBase(Row *row, size_t element_size, size_t nelement, bool is_complex,
                    const std::type_info &type_info, std::string name = "") :
            FieldBase(row, element_size * nelement, type_info, name),
            m_element_size(element_size), m_nelement(nelement), m_is_complex(is_complex) {}

    virtual std::string stats_string() const = 0;

    virtual std::string format_string() const = 0;

    template<typename T>
    static std::string stats_string_element(const T &v) {
        return utils::num_to_string(v);
    }

    template<typename T>
    static std::string stats_string_element(const std::complex<T> &v) {
        return utils::num_to_string(v.real()) + " " + utils::num_to_string(v.imag());
    }
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
                            consts::is_complex<T>(), typeid(T), name), m_format(format) {}

    NdNumberField(const NdNumberField &other) :
            NdNumberField(other.row_of_copy(), other.m_format, other.m_name) {}

    NdNumberField &operator=(const T &v) {
        std::fill(dbegin(), dend(), v);
        return *this;
    }

    NdNumberField &operator=(const std::vector<T> &v) {
        DEBUG_ASSERT_EQ(v.size(), nelement(), "Vector size does not match that of numeric field");
        std::copy(dbegin(), dend(), v.data());
        return *this;
    }

    NdNumberField &operator=(const NdNumberField &other) {
        DEBUG_ASSERT_EQ(nelement(), other.nelement(),
                              "Can't assign from incompatible instance");
        static_cast<FieldBase &>(*this) = other;
        return *this;
    }

    bool is_ordered(bool strict, bool ascending) {
        if (!nelement()) return true;
        T last_value = (*this)[0];
        for (size_t i = 1ul; i < nelement(); ++i) {
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


    /*
     * math ops
     */


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
            tmp += utils::num_to_string((*this)[ielement]) + " ";
        if (nind > 0) tmp += "]";
        return tmp;
    }

    std::string stats_string() const override {
        std::string tmp;
        for (size_t ielement = 0ul; ielement < nelement(); ++ielement)
            tmp += stats_string_element((*this)[ielement]) + " ";
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
    defs::inds h5_shape() const override {
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

template<typename T>
struct NumberField : NdNumberField<T, 0ul> {
    typedef NdNumberField<T, 0ul> base_t;
    using base_t::operator=;

    NumberField(Row *row, std::string name = "") : base_t(row, {}, name) {}

    NumberField(const NumberField& other): NumberField(other.row_of_copy(), other.m_name){}

    NumberField& operator=(const NumberField& other){
        base_t::operator=(other);
        return *this;
    }

    operator T&() {
        return *base_t::dbegin();
    }

    operator const T&() const {
        return *base_t::dbegin();
    }
};

#endif //M7_NUMBERFIELD_H
