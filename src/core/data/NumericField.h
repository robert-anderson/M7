//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICFIELD_H
#define M7_NUMERICFIELD_H

#include "NdField.h"

template<typename T>
struct NumericFieldX : FieldBaseX {
    NumericFieldX() : FieldBaseX(sizeof(T), typeid(NumericFieldX<T>)) {
        m_details["type"] = "Numeric";
        m_details["encoded type"] = consts::type_name<T>();
    }

    typedef T& view_t;
    typedef const T& const_view_t;

    const_raw_view_t convert_to_raw(const T& v) const {
        return {(char*)&v, m_element_size};
    }

    raw_view_t convert_to_raw(const T& v) {
        return {(const char*)&v, m_element_size};
    }

    view_t operator()(char *ptr) const {
        return *(T *) ptr;
    }

    std::string element_string(char *ptr) const override {
        return utils::num_to_string((*this)(ptr));
    }
};


#endif //M7_NUMERICFIELD_H
