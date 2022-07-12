//
// Created by rja on 10/07/22.
//

#ifndef M7_EXCITGENCONF_H
#define M7_EXCITGENCONF_H

#include "ConfComponents.h"

namespace conf {

    using namespace conf_components;


    struct HubbardPreferDoubleOcc : Section {
        Param<double> m_doub_occ_u_fac;
        HubbardPreferDoubleOcc(Group* parent):
            Section(parent, "hubbard_prefer_double_occ",
                    "prefers to select electrons from doubly-occupied sites",Explicit),
            m_doub_occ_u_fac(this, "doub_occ_u_fac", 1.0,
                             "factor adjusting probability of trying to draw an electron from a doubly occupied site"){}
    };

    struct ExcitGen : Section {
        HubbardPreferDoubleOcc m_hubbard_prefer_double_occ;
        ExcitGen(Group* parent): Section(parent, "excit_gen", "options related to excitation generation", Implicit),
                                 m_hubbard_prefer_double_occ(this){}
    };

}

#endif //M7_EXCITGENCONF_H
