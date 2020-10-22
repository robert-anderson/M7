//
// Created by rja on 22/10/2020.
//

#ifndef M7_BOSONONVFIELD_H
#define M7_BOSONONVFIELD_H

#include "NumericArrayField.h"
#include <numeric>

template <size_t nind>
struct BosonOnvField : NumericArrayField<uint8_t, nind, 1>{
    BosonOnvField(TableX *table, std::array<size_t, nind> shape, size_t nmode, std::string description) :
            NumericArrayField<uint8_t, nind, 1>(table, shape, {nmode}, sizeof(uint8_t), description) {}

    struct View : NumericArrayField<uint8_t, nind, 1>::View {
        View(const BosonOnvField& field, char* ptr): NumericArrayField<uint8_t, nind, 1>::View(field, ptr){}
        using FieldX::View::m_ptr;
        size_t nmode() const {
            return static_cast<const BosonOnvField<nind>&>(FieldX::View::m_field).m_format.extent(0);
        }
        size_t nboson() const {
            auto ptr = (uint8_t*)m_ptr;
            return std::accumulate(ptr, ptr+nmode(), 0ul);
        }
    };

    template<typename ...Args>
    View operator()(const size_t &irow, Args... inds) {
        return View(*this, irow, NdFieldX<nind>::m_format.flat(inds...));
    }
};

//template<typename T, size_t nind, size_t nind_view>
//struct NumericArrayField : NdArrayField<nind, nind_view> {
//    NumericArrayField(TableX *table, std::array<size_t, nind> shape, std::array<size_t, nind_view> view_shape, std::string description) :
//            NdArrayField<nind, nind_view>(table, shape, view_shape, sizeof(T), description) {}
//
//    struct View : FieldX::View {
//        View(const NumericArrayField& field, char* ptr): FieldX::View(field, ptr){}
//
//        template<typename ...Args>
//        T& operator()(Args... inds){
//            auto i = static_cast<NumericArrayField&>(m_field).m_view_format.flat(inds...);
//            return ((T*)m_ptr)[i];
//        }
//
//        View& operator=(const View& other){
//            FieldX::View::operator=(other);
//            return *this;
//        }
//    };
//
//    template<typename ...Args>
//    View operator()(const size_t &irow, Args... inds) {
//        return NdView(*this, irow, NdFieldX<nind>::m_format.flat(inds...));
//    }
//};


#endif //M7_BOSONONVFIELD_H
