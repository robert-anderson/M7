//
// Created by rja on 02/10/2020.
//

#ifndef M7_NUMERICFIELD_H
#define M7_NUMERICFIELD_H

#include <climits>
#include <src/core/util/utils.h>
#include "Field.h"
#include "Table_NEW.h"

template<typename T, size_t nind>
struct NumericField : public Field_NEW<nind> {
    using FieldBase::m_nelement;
    std::string element_to_string(size_t irow, size_t ielement) const override {
        return utils::num_to_string(flat_get(irow, ielement));
    }

    std::map<std::string, std::string> details() const override {
        auto map = Field_NEW<nind>::details();
        map["field type"] = "Numeric";
        map["encoded type"] = consts::type_name<T>();
        return map;
    }

    template<typename ...Args>
    NumericField(Table_NEW* table, std::string description, Args&& ...shape) :
            Field_NEW<nind>(table, sizeof(T), typeid(T), description, shape...){
        FieldBase::set_offsets();
    }

private:
    /*
     * Flat accessors
     */
    inline T& flat_ptr(const size_t& irow, size_t ielement){
        return ((T*)FieldBase::begin(irow))[ielement];
    }
    inline const T& flat_ptr(const size_t& irow, size_t ielement) const{
        return ((T*)FieldBase::begin(irow))[ielement];
    }
    inline T& flat_get(const size_t& irow, size_t ielement){
        return ((T*)FieldBase::begin(irow))[ielement];
    }
    inline const T& flat_get(const size_t& irow, size_t ielement) const{
        return ((T*)FieldBase::begin(irow))[ielement];
    }

    using Field_NEW<nind>::m_format;

public:
    /*
     * Formatted accessors
     */
    template<typename ...Args>
    T& operator()(const size_t& irow, Args... inds){
        return flat_get(irow, m_format.flat(inds...));
    }
    template<typename ...Args>
    const T& operator()(const size_t& irow, Args... inds) const {
        return flat_get(irow, m_format.flat(inds...));
    }
    template<typename ...Args>
    char* byte_view(const size_t& irow, Args... inds){
        return flat_ptr(irow, m_format.flat(inds...));
    }
    template<typename ...Args>
    const char* byte_view(const size_t& irow, Args... inds) const{
        return flat_ptr(irow, m_format.flat(inds...));
    }
};

#endif //M7_NUMERICFIELD_H
