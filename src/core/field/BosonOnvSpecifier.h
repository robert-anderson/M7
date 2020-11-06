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

    struct params_t {
        size_t nmode;
    };

    BosonOnvSpecifier(params_t p) : BosonOnvSpecifier(p.nmode){}

    const size_t &nmode() const;

    struct View : NumericArraySpecifier<uint8_t, 1>::View {
        View(const BosonOnvSpecifier &field, char *ptr);

        const BosonOnvSpecifier &spec() const {
            return static_cast<const BosonOnvSpecifier &>(m_spec);
        }

        size_t nboson() const;

        using NumericArraySpecifier<uint8_t, 1>::View::operator=;
    };

    typedef View view_t;
    typedef const View const_view_t;

    View operator()(char *ptr) const {
        return View(*this, ptr);
    }
};

#endif //M7_BOSONONVSPECIFIER_H
