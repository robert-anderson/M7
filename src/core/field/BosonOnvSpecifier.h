//
// Created by rja on 22/10/2020.
//

#ifndef M7_BOSONONVSPECIFIER_H
#define M7_BOSONONVSPECIFIER_H

#include "NumericArraySpecifier.h"
#include <numeric>


struct BosonOnvSpecifier : NumericArraySpecifier<uint8_t, 1> {
    BosonOnvSpecifier(size_t nmode) : NumericArraySpecifier(nmode) {
        m_data.m_details["type"] = "Boson ONV";
    }

    const size_t &nmode() const;

    struct View : NumericArraySpecifier<uint8_t, 1>::View {
        View(const BosonOnvSpecifier &spec, char *ptr);

        View(const View& other): NumericArraySpecifier<uint8_t, 1>::View(other){}

        View& operator=(const View& other){
            NumericArraySpecifier<uint8_t, 1>::View::operator=(other);
            return *this;
        }

        template<typename U>
        View& operator=(const std::vector<U> &v){
            NdAccessor<uint8_t, 1>::operator=(v);
            return *this;
        }

        const BosonOnvSpecifier &spec() const {
            return static_cast<const BosonOnvSpecifier &>(m_spec);
        }

        size_t nboson() const;

        size_t nmode() const;

    };

    typedef View view_t;
    typedef const View cview_t;

    View operator()(char *ptr) const;
};

#endif //M7_BOSONONVSPECIFIER_H
