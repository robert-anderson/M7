//
// Created by rja on 04/08/2021.
//

#include "LadderHolstein.h"

size_t LadderHolstein::approx_nconn() const {
    // assume there's one excitation or de-excitation available per electron
    return m_h.m_frm.m_nelec;
}
