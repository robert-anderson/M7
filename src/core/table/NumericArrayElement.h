//
// Created by RJA on 11/09/2020.
//

#ifndef M7_NUMERICARRAYELEMENT_H
#define M7_NUMERICARRAYELEMENT_H

#include "NumericField.h"

template <typename T>
class NumericArrayElement {
protected:
    NumericField<T>* m_field;
    const size_t m_array_size;
    const size_t m_irow, m_isegment, m_iarray;

public:
    NumericArrayElement(NumericField<T>* field, const size_t& array_size, const size_t& irow, const size_t& isegment, const size_t& iarray=0):
        m_field(field), m_array_size(array_size), m_irow(irow), m_isegment(isegment), m_iarray(iarray){}

    NumericElement<T> operator()(const size_t &ielement = 0) {
        ASSERT(ielement < m_array_size);
        return (*m_field)(m_irow, m_isegment, m_iarray*m_array_size+ielement);
    }

    NumericElement<T> operator()(const size_t &ielement = 0) const {
        ASSERT(ielement < m_array_size);
        return (*m_field)(m_irow, m_isegment, m_iarray*m_array_size+ielement);
    }

    std::vector<T> to_vector() {
        std::vector<T> tmp;
        tmp.reserve(m_array_size);
        for (size_t i = 0ul; i < m_array_size; ++i) tmp.push_back(*(*this)(i));
        return tmp;
    }


    NumericArrayElement& operator=(const std::vector<T>& v){
        ASSERT(v.size() == m_array_size);
        size_t i = 0;
        for(auto el : v){
            (*this)(i++) = el;
        }
        return *this;
    }

    NumericArrayElement& operator=(const NumericArrayElement<T>& rhs){
        if(&rhs == this) return *this;
        const auto this_begin = (*this)(0).m_begin;
        const auto rhs_begin = rhs(0).m_begin;
        const auto rhs_end = rhs_begin + sizeof(T)*std::min(rhs.m_array_size, m_array_size);
        std::memcpy(this_begin, rhs_begin, rhs_end);
        return *this;
    }

    Field* field() const{
        return m_field;
    }
};



#endif //M7_NUMERICARRAYELEMENT_H
