//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICFIELD_H
#define M7_NUMERICFIELD_H

#include "NdField.h"

template<typename T>
struct NumericField : FieldBaseX {
    NumericField() : FieldBaseX(sizeof(T), typeid(NumericField<T>)) {}

    typedef T& view_t;
    typedef const T& const_view_t;

    view_t operator()(char *ptr) const {
        return *(T *) ptr;
    }

    std::string element_string(char *ptr) const override {
        return utils::num_to_string((*this)(ptr));
    }

    std::map<std::string, std::string> details() const override {
        auto map = FieldBaseX::details();
        map["field type"] = "Numeric";
        map["encoded type"] = consts::type_name<T>();
        return map;
    }
};


#endif //M7_NUMERICFIELD_H
