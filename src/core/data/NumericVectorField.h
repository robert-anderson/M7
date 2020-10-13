//
// Created by rja on 02/10/2020.
//

#ifndef M7_NUMERICVECTORFIELD_H
#define M7_NUMERICVECTORFIELD_H

#include <src/core/util/utils.h>
#include "Field.h"

template<typename T, size_t nind>
struct NumericVectorField : public Field_NEW<nind> {
    // number of items in each element i.e. length of the vectors stored
    const size_t m_nitem;

    class View : FieldBase::View {
        using Field_NEW<nind>::View::m_ptr;
    public:
        View(const NumericVectorField<T, nind>* field, const size_t& irow, const size_t& ielement):
                Field_NEW<nind>::View(field, irow, ielement){
        }
        inline const size_t& nitem() const{
            return static_cast<const NumericVectorField<T, nind>*>(m_field)->m_nitem;
        }
        T& operator[](const size_t& iitem){
            return ((T*)m_ptr)[iitem];
        }

        const T& operator[](const size_t& iitem) const {
            return ((T*)m_ptr)[iitem];
        }

        View& operator =(const std::vector<T>& v){
            std::copy((char*)v.data(), ((char*)v.data())+std::min(m_nitem, v.size())*sizeof(T), m_ptr);
            return *this;
        }

        std::vector<T> to_vector() const {
            std::vector<T> res(nitem());
            std::copy(m_ptr, m_ptr+nitem()*sizeof(T), (char*)res.data());
            return res;
        }

        std::string to_string() const {
            return utils::to_string(to_vector());
        }
    };

    std::string element_to_string(size_t irow, size_t ielement) const override {
        return View(this, irow, ielement).to_string()+" ";
    }

    std::map<std::string, std::string> details() const override {
        auto map = Field_NEW<nind>::details();
        map["field type"] = "Numeric Vector";
        map["encoded type"] = consts::type_name<T>();
        map["vector length"] = std::to_string(m_nitem);
        return map;
    }

    template<typename ...Args>
    NumericVectorField(Table_NEW* table, size_t nitem, std::string description, Args&& ...shape) :
            Field_NEW<nind>(table, nitem*sizeof(T), typeid(T), description, shape...), m_nitem(nitem){
        FieldBase::set_offsets();
    }

    using Field_NEW<nind>::m_format;
    template<typename ...Args>
    View operator()(const size_t& irow, Args... inds){
        return View(this, irow, m_format.flat(inds...));
    }
};


#endif //M7_NUMERICVECTORFIELD_H
