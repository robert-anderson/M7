//
// Created by rja on 02/10/2020.
//

#ifndef M7_NUMERICFIELD_H
#define M7_NUMERICFIELD_H

#include <climits>
#include <src/core/util/utils.h>
#include "Field.h"
#include "Table.h"

template<typename T, size_t nind>
struct NumericField : public Field<nind> {
    using FieldBase::m_nelement;
    std::string to_string(size_t irow) const override {
        std::string res;
        for (size_t i=0ul; i<m_nelement; ++i) res+=utils::num_to_string(flat_get(irow, i))+" ";
        return res;
    }

    template<typename ...Args>
    NumericField(Table* table, Args&& ...shape) :
    Field<nind>(table, sizeof(T), typeid(T), shape...){
        FieldBase::set_offsets();
    }

private:
    /*
     * Flat accessors
     */
    inline T& flat_get(const size_t& irow, size_t ielement){
        return ((T*)FieldBase::begin(irow))[ielement];
    }
    inline const T& flat_get(const size_t& irow, size_t ielement) const{
        return ((T*)FieldBase::begin(irow))[ielement];
    }

    using Field<nind>::m_format;

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
        return ((T*)FieldBase::begin(irow))[m_format.flat(inds...)];
    }
};

#endif //M7_NUMERICFIELD_H
