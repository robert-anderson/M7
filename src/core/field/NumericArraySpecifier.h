//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICARRAYSPECIFIER_H
#define M7_NUMERICARRAYSPECIFIER_H

#include <src/core/util/utils.h>
#include "FieldSpecifier.h"
#include "src/core/nd/NdFormat.h"

template<typename T, size_t nind>
struct NumericArraySpecifier : FieldSpecifier {
    NdFormat<nind> m_format;
    const size_t m_nelement;
    template<typename ...Args>
    NumericArraySpecifier(Args... shape) :
    FieldSpecifier(NdFormat<nind>(shape...).nelement() * sizeof(T), typeid(NumericArraySpecifier<T, nind>)),
    m_format(shape...), m_nelement(m_format.nelement()){
        m_data.m_details["type"] = "Numeric Array";
        m_data.m_details["encoded type"] = consts::type_name<T>();
        m_data.m_details["element dimensionality"] = std::to_string(nind);
        m_data.m_details["element shape"] = utils::to_string(m_format.shape());
    }

    struct View : FieldSpecifier::View {
        View(const NumericArraySpecifier &spec, char *ptr) : FieldSpecifier::View(spec, ptr) {}

        const size_t& nelement() const {
            return static_cast<const NumericArraySpecifier&>(m_spec).m_nelement;
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
            return ((T*)m_ptr)[static_cast<const NumericArraySpecifier&>(m_spec).m_format.flatten(inds...)];
        }
        template<typename ...Args>
        const T& operator()(Args... inds) const {
            return ((T*)m_ptr)[static_cast<const NumericArraySpecifier&>(m_spec).m_format.flatten(inds...)];
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

#endif //M7_NUMERICARRAYSPECIFIER_H
