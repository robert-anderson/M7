//
// Created by rja on 04/11/2020.
//

#include "FermionOnvEnumerator.h"

FermionOnvEnumerator::FermionOnvEnumerator(size_t nsite, size_t nelec) :
        m_internal_enum(new SpinNonConFonvEnumerator(nsite, nelec)) {}
