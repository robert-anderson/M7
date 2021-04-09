//
// Created by rja on 09/02/2021.
//

#ifndef M7_NUMBERFIELD_H
#define M7_NUMBERFIELD_H

#include "FieldBase.h"

template<typename T, size_t nind>
struct NdNumberField : FieldBase {
    typedef const std::array<size_t, nind>& inds_t;
    const NdFormat<nind> m_format;

    const size_t& nelement() const {
        return m_format.nelement();
    }

    NdNumberField(Row* row, NdFormat<nind> format):
            FieldBase(row, sizeof(T) * format.nelement(), typeid(T)), m_format(format){}

    NdNumberField(const NdNumberField& other):
            NdNumberField(other.row_of_copy(), other.m_format){}

    NdNumberField &operator=(const T &v) {
        std::fill((T*)begin(), (T*)end(), v);
        return *this;
    }

    NdNumberField &operator=(const std::vector<T> &v) {
        ASSERT(v.size()==nelement());
        std::memcpy(begin(), v.data(), m_size);
        return *this;
    }

    NdNumberField &operator=(const NdNumberField &v) {
        static_cast<FieldBase&>(*this) = v;
        return *this;
    }

    T& operator[](const size_t& ielement) {
        return ((T *) begin())[ielement];
    }

    const T& operator[](const size_t& ielement) const {
        return ((const T *) begin())[ielement];
    }

    T& operator[](inds_t inds) {
        return ((T *) begin())[m_format.flatten(inds)];
    }

    const T& operator[](inds_t inds) const {
        return ((const T *) begin())[m_format.flatten(inds)];
    }

    std::string to_string() const override {
        std::string tmp;
        if (nind>0) tmp += "[";
        for (size_t ielement = 0ul; ielement<nelement(); ++ielement)
            tmp+=std::to_string((*this)[ielement]) + " ";
        if (nind>0) tmp += "]";
        return tmp;
    }

    defs::inds h5_shape() const override {
        return {m_format.shape().begin(), m_format.shape().end()};
    }

    std::vector<std::string> h5_dim_names() const override {
        if (!nind) return {};
        return {m_format.dim_names().begin(), m_format.dim_names().end()};
    }

    hid_t h5_type() const override {
        return hdf5::type<T>();
    }
};

template<typename T>
struct NumberField : NdNumberField<T, 0ul> {
    typedef NdNumberField<T, 0ul> base_t;
    using base_t::operator=;

    NumberField(Row* row): base_t(row, {}){}

    operator T&() {
        return *(T*)FieldBase::begin();
    }

    operator const T&() const {
        return *(const T*)FieldBase::begin();
    }
};


#endif //M7_NUMBERFIELD_H
