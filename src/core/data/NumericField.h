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
    struct View : Field<nind>::View {
        using Field<nind>::View::m_ptr;
        View(const NumericField<T, nind>& field, const size_t& irow, const size_t& ielement):
                Field<nind>::View(field.begin(irow)+ielement*sizeof(T)){
        }
        operator T&() {
            return *((T*)m_ptr);
        }
        operator const T&() const {
            return *((T*)m_ptr);
        }
        View& operator =(const T& v){
            *((T*)m_ptr) = v;
            return *this;
        }
        std::string to_string() const {
            return utils::num_to_string((T)*this);
        }
    };

    using FieldBase::m_nelement;
    std::string to_string(size_t irow) const override {
        std::string res;
        for (size_t i=0ul; i<m_nelement; ++i) res+=View(*this, irow, i).to_string()+" ";
        return res;
    }

    template<typename ...Args>
    NumericField(Table* table, Args&& ...shape) :
    Field<nind>(table, sizeof(T), typeid(T), shape...){
        FieldBase::set_offsets();
    }

    using Field<nind>::m_format;
    template<typename ...Args>
    View operator()(const size_t& irow, Args... inds){
        return View(*this, irow, m_format.flat(inds...));
    }
};

#endif //M7_NUMERICFIELD_H
