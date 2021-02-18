//
// Created by rja on 09/02/2021.
//

#ifndef M7_NUMBERFIELDZ_H
#define M7_NUMBERFIELDZ_H

#include "NdFieldBaseZ.h"

template<typename T>
struct NumberFieldBaseZ : FieldBaseZ {
    const size_t m_nelement;

    NumberFieldBaseZ(RowZ* row, size_t nitem, size_t nelement):
    FieldBaseZ(row, sizeof(T)*nelement, nitem, typeid(T)), m_nelement(nelement){}

    NumberFieldBaseZ(const NumberFieldBaseZ& other):
    NumberFieldBaseZ(other.m_row ? other.m_row->m_child : other.m_row, other.m_nitem, other.m_nelement){}

    NumberFieldBaseZ &operator=(const T &v) {
        std::fill((T*)begin(), (T*)end(), v);
        return *this;
    }

    NumberFieldBaseZ &operator=(const std::vector<T> &v) {
        ASSERT(v.size() == m_nitem*m_nelement);
        std::memcpy(begin(), v.data(), m_size);
        return *this;
    }

    NumberFieldBaseZ &operator=(const NumberFieldBaseZ &v) {
        static_cast<FieldBaseZ&>(*this) = v;
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
struct NumberFieldZ : NumberFieldBaseZ<T> {
    using NumberFieldBaseZ<T>::operator=;

    NumberFieldZ(RowZ* row): NumberFieldBaseZ<T>(row, 1, 1){}

    operator T&() {
        return *(T *) FieldBaseZ::begin();
    }

    operator const T&() const {
        return *(const T *) FieldBaseZ::begin();
    }
};


template<typename T>
struct VectorFieldZ : NumberFieldBaseZ<T> {
    VectorFieldZ(RowZ* row, size_t nelement): NumberFieldBaseZ<T>(row, 1, nelement){}

    T& operator()(const size_t& ielement) {
        return ((T *) FieldBaseZ::begin())[ielement];
    }

    const T& operator()(const size_t& ielement) const {
        return ((const T *) FieldBaseZ::begin())[ielement];
    }
};


template<typename T>
struct VectorsFieldZ : NumberFieldBaseZ<T> {
    VectorsFieldZ(RowZ* row, size_t nitem, size_t nelement): NumberFieldBaseZ<T>(row, nitem, nelement){}

    T& operator()(const size_t& iitem, const size_t& ielement) {
        return ((T *) FieldBaseZ::begin(iitem))[ielement];
    }

    const T& operator()(const size_t& iitem, const size_t& ielement) const {
        return ((const T *) FieldBaseZ::begin(iitem))[ielement];
    }
};

#endif //M7_NUMBERFIELDZ_H
