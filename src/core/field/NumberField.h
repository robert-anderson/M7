//
// Created by rja on 09/02/2021.
//

#ifndef M7_NUMBERFIELD_H
#define M7_NUMBERFIELD_H

#include "FieldBase.h"

template<typename T>
struct NumberFieldBase : FieldBase {
    const size_t m_nelement;

    NumberFieldBase(Row* row, size_t nitem, size_t nelement):
            FieldBase(row, sizeof(T) * nelement, nitem, typeid(T), hdf5::type<T>()), m_nelement(nelement){}

    NumberFieldBase(const NumberFieldBase& other):
            NumberFieldBase(other.m_row ? other.m_row->m_child : other.m_row, other.m_nitem, other.m_nelement){}

    NumberFieldBase &operator=(const T &v) {
        std::fill((T*)begin(), (T*)end(), v);
        return *this;
    }

    NumberFieldBase &operator=(const std::vector<T> &v) {
        ASSERT(v.size() == m_nitem*m_nelement);
        std::memcpy(begin(), v.data(), m_size);
        return *this;
    }

    NumberFieldBase &operator=(const NumberFieldBase &v) {
        static_cast<FieldBase&>(*this) = v;
        return *this;
    }

    T& get(const size_t& iitem, const size_t& ielement){
        return ((T *) begin(iitem))[ielement];
    }

    const T& get(const size_t& iitem, const size_t& ielement) const{
        return ((T *) begin(iitem))[ielement];
    }

    std::string to_string_element(const size_t& iitem) const override {
        std::string tmp;
        if (m_nelement>1) tmp += "[";
        for (size_t ielement = 0ul; ielement<m_nelement; ++ielement)
            tmp+=std::to_string(this->get(iitem, ielement)) + " ";
        if (m_nelement>1) tmp += "]";
        return tmp;
    }
};


template<typename T>
struct NumberField : NumberFieldBase<T> {
    using NumberFieldBase<T>::operator=;

    NumberField(Row* row): NumberFieldBase<T>(row, 1, 1){}

    operator T&() {
        return *(T *) FieldBase::begin();
    }

    operator const T&() const {
        return *(const T *) FieldBase::begin();
    }
};


template<typename T>
struct VectorField : NumberFieldBase<T> {
    VectorField(Row* row, size_t nelement): NumberFieldBase<T>(row, 1, nelement){}

    T& operator()(const size_t& ielement) {
        return ((T *) FieldBase::begin())[ielement];
    }

    const T& operator()(const size_t& ielement) const {
        return ((const T *) FieldBase::begin())[ielement];
    }
};


template<typename T>
struct VectorsField : NumberFieldBase<T> {
    VectorsField(Row* row, size_t nitem, size_t nelement): NumberFieldBase<T>(row, nitem, nelement){}

    T& operator()(const size_t& iitem, const size_t& ielement) {
        return ((T *) FieldBase::begin(iitem))[ielement];
    }

    const T& operator()(const size_t& iitem, const size_t& ielement) const {
        return ((const T *) FieldBase::begin(iitem))[ielement];
    }
};

#endif //M7_NUMBERFIELD_H
