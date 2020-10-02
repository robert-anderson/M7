//
// Created by rja on 02/10/2020.
//

#ifndef M7_NUMERICVECTORFIELD_H
#define M7_NUMERICVECTORFIELD_H

#include <src/core/util/utils.h>
#include "Field.h"
#include "NumericField.h"

template<typename T, size_t nind>
struct NumericVectorField : public Field<nind> {
    // number of items in each element i.e. length of the vectors stored
    const size_t m_nitem;

    struct View : Field<nind>::View {
        using Field<nind>::View::m_ptr;
        const size_t& m_nitem;
        View(const NumericVectorField<T, nind>& field, const size_t& irow, const size_t& ielement):
                Field<nind>::View(field.begin(irow)+ielement*sizeof(T)), m_nitem(field.m_nitem){
        }

        T& operator[](const size_t& ientry){
            return ((T*)m_ptr)[ientry];
        }

        const T& operator[](const size_t& ientry) const {
            return ((T*)m_ptr)[ientry];
        }

        View& operator =(const std::vector<T>& v){
            std::copy((char*)v.data(), ((char*)v.data())+std::min(m_nitem, v.size())*sizeof(T), m_ptr);
            return *this;
        }

        std::vector<T> to_vector() const {
            std::vector<T> res(m_nitem);
            std::copy(m_ptr, m_ptr+m_nitem*sizeof(T), (char*)res.data());
            return res;
        }

        std::string to_string() const {
            return utils::to_string(to_vector());
        }
    };

    using FieldBase::m_nelement;
    std::string to_string(size_t irow) const override {
        std::string res;
        for (size_t i=0ul; i<m_nelement; ++i) res+=View(*this, irow, i).to_string()+" ";
        return res;
    }

    template<typename ...Args>
    NumericVectorField(Table* table, size_t nitem, Args&& ...shape) :
            Field<nind>(table, nitem*sizeof(T), typeid(T), shape...), m_nitem(nitem){
        FieldBase::set_offsets();
    }

    using Field<nind>::m_format;
    template<typename ...Args>
    View operator()(const size_t& irow, Args... inds){
        return View(*this, irow, m_format.flat(inds...));
    }
};


#endif //M7_NUMERICVECTORFIELD_H
