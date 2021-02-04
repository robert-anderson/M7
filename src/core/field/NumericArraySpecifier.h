//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICARRAYSPECIFIER_H
#define M7_NUMERICARRAYSPECIFIER_H

#include <src/core/util/utils.h>
#include "ColumnSpecifier.h"
#include "src/core/nd/NdFormat.h"
#include "src/core/nd/NdAccessor.h"

template<typename T, size_t nind>
struct NumericArraySpecifier : ColumnSpecifier {
    NdFormat<nind> m_format;
    const size_t m_nelement;
    template<typename ...Args>
    NumericArraySpecifier(Args... shape) :
    ColumnSpecifier(NdFormat<nind>(shape...).nelement() * sizeof(T), typeid(NumericArraySpecifier<T, nind>)),
    m_format(shape...), m_nelement(m_format.nelement()){
        m_data.m_details["type"] = "Numeric Array";
        m_data.m_details["encoded type"] = consts::type_name<T>();
        m_data.m_details["element dimensionality"] = std::to_string(nind);
        m_data.m_details["element shape"] = utils::to_string(m_format.shape());
    }

    struct View : ColumnSpecifier::View, NdAccessor<T, nind> {
        View(const NumericArraySpecifier &spec, char *ptr) :
                ColumnSpecifier::View(spec, ptr),
                NdAccessor<T, nind>((T*)ptr, spec.m_format){}

        View(const View& other):ColumnSpecifier::View(other), NdAccessor<T, nind>(other){}

        using NdAccessor<T, nind>::nelement;
        using NdAccessor<T, nind>::operator=;
        using ColumnSpecifier::View::operator=;

        std::string to_string() const override {
            std::string res;
            res+="[";
            for (size_t i=0ul; i<nelement(); ++i) res+=utils::num_to_string((*this)[i])+" ";
            res+="]";
            return res;
        }

        View& operator=(const View& other){
            ColumnSpecifier::View::operator=(other);
            return *this;
        }

        template<typename U>
        View& operator=(const std::vector<U> &v){
            NdAccessor<T, nind>::operator=(v);
            return *this;
        }

    };

    std::string element_string(char *ptr) const override {
        return View(*this, ptr).to_string();
    }

    typedef View view_t;

    View operator()(char *ptr) const {
        return View(*this, ptr);
    }
};

#endif //M7_NUMERICARRAYSPECIFIER_H
