//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICSPECIFIER_H
#define M7_NUMERICSPECIFIER_H

#include "ColumnSpecifier.h"
#include "src/core/util/utils.h"

template<typename T>
struct NumericSpecifier : ColumnSpecifier {
    NumericSpecifier() : ColumnSpecifier(sizeof(T), typeid(NumericSpecifier<T>)) {
        m_data.m_details["type"] = "Numeric";
        m_data.m_details["encoded type"] = consts::type_name<T>();
    }

    typedef T& view_t;
    typedef const T& cview_t;

    static defs::hash_t hash(const View& view) {
        return hashing::fnv_hash((char*)&view, sizeof(T));
    }

    view_t operator()(char *ptr) const {
        return *(T *) ptr;
    }

    std::string element_string(char *ptr) const override {
        return utils::num_to_string((*this)(ptr));
    }
};


#endif //M7_NUMERICSPECIFIER_H
