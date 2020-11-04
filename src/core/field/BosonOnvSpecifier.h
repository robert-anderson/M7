//
// Created by rja on 22/10/2020.
//

#ifndef M7_BOSONONVSPECIFIER_H
#define M7_BOSONONVSPECIFIER_H

#include "NumericArraySpecifier.h"
#include <numeric>


struct BosonOnvSpecifier : NumericArraySpecifier<uint8_t, 1> {
    template<typename ...Args>
    BosonOnvSpecifier(size_t nmode):NumericArraySpecifier(nmode) {
        m_data.m_details["type"] = "Boson ONV";
    }

    const size_t& nmode() const;

    struct View : NumericArraySpecifier<uint8_t, 1>::View {
        View(const BosonOnvSpecifier &field, char *ptr);

        size_t nboson() const;
    };
};

#endif //M7_BOSONONVSPECIFIER_H
