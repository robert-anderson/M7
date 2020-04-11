//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_NUMERICELEMENT_H
#define SANDBOX2_NUMERICELEMENT_H

#include <src/utils.h>
#include "Element.h"

template<typename T>
class NumericField;

template<typename T>
class NumericElement : public Element {
    T m_internal_value;
public:
    typedef NumericField<T> Field_T;

    NumericElement(NumericField<T> *field, char *begin) :
        Element(field, begin) {}

    NumericElement(const T &value) : NumericElement(nullptr, nullptr) {
        m_internal_value = value;
        m_begin = (char *) &m_internal_value;
    }

    // explicit conversion
    T& operator*() const {return *((T *) m_begin);}

    explicit operator T() const { return *((T *) m_begin); }

    //explicit operator T&() const { return *((T *) m_begin); }

    virtual NumericElement<T> &operator=(const T &v) {
        *((T *) m_begin) = v;
        return *this;
    }

    virtual NumericElement<T> &operator=(const NumericElement<T> &v) {
        std::memcpy(m_begin, v.m_begin, sizeof(T));
        return *this;
    }

    bool operator==(const T &v) const { return *((T *) m_begin) == v; }

    bool operator==(const NumericElement<T> &other) const {
        return T(*this) == T(other);
    }

    bool operator!=(const T &v) const { return !(*this == v); }

    bool operator!=(const NumericElement<T> &other) const {
        return T(*this) != T(other);
    }

    NumericElement<T> &operator*=(const T& v){
        **this*=v;
        return *this;
    }

    NumericElement<T> &operator/=(const T& v){
        **this/=v;
        return *this;
    }

    NumericElement<T> &operator+=(const T& v){
        **this+=v;
        return *this;
    }

    NumericElement<T> &operator-=(const T& v){
        **this*=v;
        return *this;
    }

    bool is_zero() const override {
        return *this == 0;
    }

    size_t size() const override {
        return sizeof(T);
    }

    virtual std::string to_string() {
        return utils::num_to_string(**this);
    }
};

#endif //SANDBOX2_NUMERICELEMENT_H
