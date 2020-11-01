//
// Created by rja on 22/10/2020.
//

#ifndef M7_BOSONONVFIELD_H
#define M7_BOSONONVFIELD_H

#include "NumericArrayField.h"
#include <numeric>


struct BosonOnvField : NumericArrayField<uint8_t, 1> {
    template<typename ...Args>
    BosonOnvField(size_t nmode):NumericArrayField(nmode) {
        m_details["type"] = "Boson ONV";
    }

    const size_t& nmode() const {
        return m_format.extent(0);
    }

    struct View : NumericArrayField<uint8_t, 1>::View {
        View(const BosonOnvField &field, char *ptr) : NumericArrayField<uint8_t, 1>::View(field, ptr) {}

        size_t nboson() const {
            return std::accumulate(
                    reinterpret_cast<uint8_t *>(m_ptr),
                    reinterpret_cast<uint8_t *>(m_ptr) + nelement(), 0ul);
        }
    };
};

#if 0
template<size_t nind>
struct BosonOnvField : NumericArrayField<uint8_t, nind, 1> {
    typedef NumericArrayField<uint8_t, nind, 1> base_t;
    BosonOnvField(TableX *table, std::array<size_t, nind> shape, size_t nmode, std::string description) :
            base_t(table, shape, {nmode}, description) {}

    size_t nmode() const {
        return NumericArrayField<uint8_t, nind, 1>::m_view_format.extent(0);
    }

    struct View : base_t::View {
        template<typename ...Args>
        View(const BosonOnvField &field, const size_t &irow, const size_t& iflat):
                base_t::View(field, irow, iflat) {}

        using FieldX::View::m_ptr;
        using FieldX::View::m_field;

        size_t nmode() const {
            return static_cast<const BosonOnvField<nind> &>(m_field).nmode();
        }

        size_t nboson() const {
            auto ptr = (uint8_t *) m_ptr;
            return std::accumulate(ptr, ptr + nmode(), 0x0);
        }
    };

    template<typename ...Args>
    View operator()(const size_t &irow, Args... inds) {
        return View(*this, irow, base_t::m_format.flatten(inds...));
    }

    std::map<std::string, std::string> details() const override {
        auto map = NumericArrayField<uint8_t, nind, 1>::details();
        map["field type"] = "Boson Occupation Number Vector";
        map["number of modes"] = std::to_string(nmode());
        return map;
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
//            auto i = static_cast<NumericArrayField&>(m_field).m_view_format.flatten(inds...);
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
//        return NdView(*this, irow, NdFieldX<nind>::m_format.flatten(inds...));
//    }
//};


#endif //M7_BOSONONVFIELD_H
#endif //M7_BOSONONVFIELD_H
