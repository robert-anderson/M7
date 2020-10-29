//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICARRAYFIELD_H
#define M7_NUMERICARRAYFIELD_H

#include <src/core/util/utils.h>
#include "Field.h"
#include "src/core/nd/NdFormat.h"

template<typename T, size_t nind>
struct NumericArrayField : FieldBaseX {
    NdFormat<nind> m_format;
    const size_t m_nelement;
    template<typename ...Args>
    NumericArrayField(Args... shape) :
    FieldBaseX(NdFormat<nind>(shape...).nelement()*sizeof(T), typeid(NumericArrayField<T, nind>)),
    m_format(shape...), m_nelement(m_format.nelement()){
        m_details["type"] = "Numeric Array";
        m_details["encoded type"] = consts::type_name<T>();
        m_details["element dimensionality"] = std::to_string(nind);
        m_details["element shape"] = utils::to_string(m_format.shape());
    }

    struct View : FieldBaseX::View {
        View(const NumericArrayField &field, char *ptr) : FieldBaseX::View(field, ptr) {}

        const size_t& nelement() const {
            return static_cast<const NumericArrayField&>(m_field).m_nelement;
        }

        T& operator[](const size_t& ind) {
            ASSERT(ind<nelement());
            return ((T*)m_ptr)[ind];
        }
        const T& operator[](const size_t& ind) const {
            ASSERT(ind<nelement());
            return ((T*)m_ptr)[ind];
        }
        template<typename ...Args>
        T& operator()(Args... inds){
            return ((T*)m_ptr)[static_cast<const NumericArrayField&>(m_field).m_format.flatten(inds...)];
        }
        template<typename ...Args>
        const T& operator()(Args... inds) const {
            return ((T*)m_ptr)[static_cast<const NumericArrayField&>(m_field).m_format.flatten(inds...)];
        }

        std::string to_string() const override {
            std::string res;
            res+="[";
            for (size_t i=0ul; i<nelement(); ++i) res+=utils::num_to_string((*this)[i])+" ";
            res+="]";
            return res;
        }
    };

    std::string element_string(char *ptr) const override {
        return View(*this, ptr).to_string();
    }

    typedef View view_t;
    typedef const View const_view_t;

    View operator()(char *ptr) const {
        return View(*this, ptr);
    }
};

#endif //M7_NUMERICARRAYFIELD_H
