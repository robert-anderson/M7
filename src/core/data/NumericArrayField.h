//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICARRAYFIELD_H
#define M7_NUMERICARRAYFIELD_H


#include "NdArrayField.h"

template<typename T, size_t nind, size_t nind_view>
struct NumericArrayField : NdArrayField<nind, nind_view> {
    NumericArrayField(TableX *table, std::array<size_t, nind> shape, std::array<size_t, nind_view> view_shape, std::string description) :
            NdArrayField<nind, nind_view>(table, shape, view_shape, sizeof(T), description) {}

    struct View : FieldX::View {
        View(const NumericArrayField& field, char* ptr): FieldX::View(field, ptr){}

        template<typename ...Args>
        T& operator()(Args... inds){
            auto i = static_cast<NumericArrayField&>(m_field).m_view_format.flat(inds...);
            return ((T*)m_ptr)[i];
        }

        View& operator=(const View& other){
            FieldX::View::operator=(other);
            return *this;
        }
    };

    template<typename ...Args>
    View operator()(const size_t &irow, Args... inds) {
        return View(*this, irow, NdFieldX<nind>::m_format.flat(inds...));
    }


};


#endif //M7_NUMERICARRAYFIELD_H
