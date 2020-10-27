//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICARRAYFIELD_H
#define M7_NUMERICARRAYFIELD_H


#include "NdArrayField.h"

#if 0
template<typename T>
struct NumericArrayFieldKind{};

template<typename T, size_t nind, size_t nind_view>
struct NumericArrayField : NdArrayField<nind, nind_view> {
    NumericArrayField(TableX *table, std::array<size_t, nind> shape, std::array<size_t, nind_view> view_shape,
                      std::string description) :
            NdArrayField<nind, nind_view>(table, shape, view_shape, sizeof(T),
                                          typeid(NumericArrayFieldKind<T>), description) {}

    struct View : NdArrayField<nind, nind_view>::View {
        template<typename ...Args>
        View(const NumericArrayField &field, const size_t &irow, const size_t& iflat):
                NdArrayField<nind, nind_view>::View(field, irow, iflat) {}

        using FieldX::View::m_field;
        using FieldX::View::m_ptr;
        size_t nelement() const {
            return static_cast<const NumericArrayField &>(m_field).m_view_format.nelement();
        }

        T &operator[](const size_t &iflat) {
            return ((T *) m_ptr)[iflat];
        }

        const T &operator[](const size_t &iflat) const {
            return ((T *) m_ptr)[iflat];
        }

        template<typename ...Args>
        T &operator()(Args... inds) {
            auto i = static_cast<NumericArrayField &>(m_field).m_view_format.flat(inds...);
            return (*this)[i];
        }

        template<typename ...Args>
        const T &operator()(Args... inds) const {
            auto i = static_cast<NumericArrayField &>(m_field).m_view_format.flat(inds...);
            return (*this)[i];
        }

        View(const View &other) : NdArrayField<nind, nind_view>::View(other) {}

        View &operator=(const View &other) {
            FieldX::View::operator=(other);
            return *this;
        }

        std::string to_string() const {
            std::string res = "[";
            for (size_t i = 0ul; i < nelement(); ++i) res += utils::num_to_string((*this)[i])+" ";
            res += "]";
            return res;
        }
    };

    template<typename ...Args>
    View operator()(const size_t &irow, Args... inds) {
        return View(*this, irow, NdArrayField<nind, nind_view>::m_format.flatten(inds...));
    }

    std::string element_string(size_t irow, size_t ielement) const override {
        return View(*this, irow, ielement).to_string();
    }

    std::map<std::string, std::string> details() const override {
        auto map = NdArrayField<nind, nind_view>::details();
        map["field type"] = "Numeric Array";
        map["encoded type"] = consts::type_name<T>();
        return map;
    }
};


#endif //M7_NUMERICARRAYFIELD_H
#endif //M7_NUMERICARRAYFIELD_H
