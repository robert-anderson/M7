//
// Created by rja on 04/11/2020.
//

#ifndef M7_FERMIONONVENUMERATOR_H
#define M7_FERMIONONVENUMERATOR_H

#include "SpinConFonvEnumerator.h"
#include "SpinNonConFonvEnumerator.h"

class FermionOnvEnumerator : public Enumerator<field::FrmOnv> {
    std::unique_ptr<Enumerator<field::FrmOnv>> m_internal_enum;

public:
    FermionOnvEnumerator(size_t nsite, size_t nelec);

    FermionOnvEnumerator(size_t nsite, size_t nelec, int spin) :
            m_internal_enum(new SpinConFonvEnumerator(nsite, nelec, spin)) {}

    bool next_element(field::FrmOnv &result) override {
        return m_internal_enum->next(result);
    }
};


#endif //M7_FERMIONONVENUMERATOR_H